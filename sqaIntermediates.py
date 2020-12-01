import string
import sys, os, time
import subprocess
import collections
import math
import numpy as np

from sqaTensor import tensor, creOp, desOp, kroneckerDelta
from sqaTensor2 import creDesTensor
from sqaTerm import term, combineTerms
from sqaIndex import index
from sqaOptions import options
from sqaSymmetry import symmetry
from sqaCommutator import commutator


def genIntermediates(input_terms, ind_str = None, factor_depth = 1, transRDM = False):

    if factor_depth < 0 or not isinstance(factor_depth, int):
        raise ValueError('Invalid factor depth provided -- provide integer value that is >= 0') 

    # Make list of integers to append to 'INT' string below
    interm_name_list = np.arange(1, 10000)

    # Create new list of terms that will modify input terms and expand them in terms of intermediates:
    mod_term_list = []

    # Create list for storing intermediates
    intermediates = []

    # Iterate through every term in list of terms
    for in_term in input_terms:

        # Create lists for all tensors
        tensorlist          = []
        tensor_indices_list = []
        credes_list         = []
        pre_factor          = in_term.numConstant

        # Reformat the names of the tensors from SQA
        for t in in_term.tensors:

            # Separate cre/des operators to make into RDMs
            if (isinstance(t, creOp) or isinstance(t, desOp)):
                credes_list.append(t)            

            else:
                tensorlist.append(t)
                tensor_indices_list.append([ind for ind in t.indices])

        # Turn cre/des operators into RDM, if they exist
        if (len(credes_list) > 0):
            rdm_tensor = creDesTensor(credes_list, transRDM)
            tensorlist.append(rdm_tensor)
            tensor_indices_list.append([ind for ind in rdm_tensor.indices])

        ## Create einsum string
        lhs_str = []
        for ten in tensor_indices_list:
            lhs_str.append(''.join([i.name for i in ten]))

        lhs_str = ','.join(lhs_str)            
        rhs_str = '->' + ind_str
        einsum_string = lhs_str + rhs_str

        # Create dictionary w/ indices as keys and value corresponding to dimension of that index
        unique_inds = set(einsum_string) - {',', '-', '>'}
        index_size = [2] * len(unique_inds)
        sizes_dict = dict(zip(unique_inds, index_size))

        # Construct dummy tensors for term in order to assess contraction path
        dummy_tens = []
        
        for t in tensorlist:
            dims = []
            
            for index in t.indices:
                dims.append(sizes_dict[str(index.name)])

            dummy_tens.append(np.empty(tuple(dims)))

        # Compute most efficient contraction path
        path_info = np.einsum_path(einsum_string, *dummy_tens)        

        # TODO: WHICH EINSUM? OPT OR NUMPY?
        # Isolate intermediate contractions from opt_einsum info
        naive          = int(path_info[1].split('\n')[1].split()[-1])
        opt            = int(path_info[1].split('\n')[2].split()[-1])

        # Determine if opt_einsum will improve scaling of contraction
        if naive > opt:

            # Save tuples that indicate optimized order of contracting tensors
            contract_order = [contract for contract in path_info[0][1:-factor_depth]]

            # Determine contraction path and indices of intermediates
            split_path     = path_info[1].split('\n')[10:-factor_depth]
            int_indices    = [str(line).split()[1].split('->')[1] for line in split_path]

            # Define scale outside of loop
            scale_factor_total = 1.0

            # Create intermediate
            for num, contract in enumerate(contract_order):
                
                # Make intermediate name
                tensor_name = 'INT' + str(interm_name_list[0])

                # Determine which tensors from tensorlist are being contracted
                tens_contract = [tensorlist[i] for i in contract]
               
                # Use the external string to modify indexType of the indices in tensors wrt the intermediate term
                new_tensors, def_indices, loop_indices = get_int_indices(tens_contract, int_indices[num])

                # Construct intermediate term w/ updated tensor
                int_term = term(1.0, [], new_tensors)
               
                # Canonicalize term and tensor representation of intermediate and update scale factor
                int_term, scale_factor = make_canonical(int_term, transRDM)
                scale_factor_total *= scale_factor                

                # Update indices after canonicalizing term
                new_tensors, def_indices, loop_indices = get_int_indices(int_term.tensors, int_indices[num])

                # Define intermediate tensor wrt to definition
                int_tensor  = tensor(tensor_name, def_indices, [])

                # Append the first intermeidate automatically
                if not intermediates:
                    intermediates.append([int_term, int_tensor])
                    interm_name_list = interm_name_list[1:]
               
                # Once intermediates list is not empty, check all other intermediates for redundancy against the list
                else:
                    int_term, int_tensor, isRedundant = check_intermediates(intermediates, int_term, int_tensor)

                    # Only append unique intermediate terms to the list of intermediates
                    if not isRedundant:
                        intermediates.append([int_term, int_tensor])
                        interm_name_list = interm_name_list[1:]
             
                # Modify 'tensorlist' for einsum's contract_path function
                tensorlist = [tens for tens in tensorlist if tens not in tens_contract]
                
                # Append representation of INT tensor w/ dummy/external indices defined wrt full contraction
                loop_tensor = tensor(int_tensor.name, loop_indices, [])
                tensorlist.append(loop_tensor)

            pre_factor *= scale_factor_total
            mod_term_list.append(term(pre_factor, [], tensorlist))

        # Scaling of contraction cannot be optimized
        else:
            mod_term_list.append(in_term)
    
    return mod_term_list, intermediates


def make_canonical(int_term, transRDM):

    # Additional canonicalization for RDM tensors
    for tens in [ten for ten in int_term.tensors]:

        if isinstance(tens, creDesTensor):
            int_term = canonicalize_RDM(int_term, transRDM)
            break

    # Create ranking of indices based on the contraction path
    path_rank  = assign_path_rank(int_term)
    space_rank = assign_space_rank(int_term) 

    # Store list of tensors that are in canonical order to create a canonicalized term
    canonTensorList = []
    scale_factor    = int_term.numConstant

    ### RE-ARRANGE THE INDICES WITHIN EACH TENSOR IN THE TERM
    # Continue canonicalizing after modifying RDM tensors    
    for t_ind, t in enumerate(int_term.tensors):

        # Keep track of path rank for each tensor
        loop_path_rank  = path_rank[:len(t.indices)]
        loop_space_rank = space_rank[:len(t.indices)]

        # Modify path rank for RDM des operators
        if isinstance(t, creDesTensor):
           for i, op in enumerate(t.ops):
               if isinstance(op, desOp):
                   loop_path_rank[i] += 100

        # Create final ranking for tensor
        final_rank = []

        for i in range(len(t.indices)):
            final_rank.append(loop_path_rank[i] + loop_space_rank[i])
        
        # Cut used elements of rank lists
        path_rank  = path_rank[len(t.indices):]
        space_rank = space_rank[len(t.indices):]

        # Generate permutation of indices to canonical order
        canon_order = np.argsort(final_rank)         
        
        # Get symmetry information from sqa_tensor
        allowed_sym  = t.symPermutes()

        # Only modify tensors with appropriate symmetry
        if list(canon_order) in allowed_sym[0]:

            # Get index of symmetry for convenience
            ind_sym = allowed_sym[0].index(list(canon_order))

            # Keep track of pre_factor
            scale_factor  *= float(allowed_sym[1][ind_sym])                
            sorted_indices = [t.indices[i] for i in canon_order]
            sorted_tensor  = tensor(t.name, sorted_indices, t.symmetries)
 
            # Append tensor w/ sorted indices
            canonTensorList.append(sorted_tensor)
    
        # Append tensor w/ unpermuted indices
        else:
            canonTensorList.append(t)

    # Form term with tensors w/ canonicalized indices
    canon_index_term = term(1.0, [], canonTensorList)

    ### RE-ARRANGE TENSORS THAT MAKE UP THE TERM
    ## SORT BY RANK
    ranked_term, rank_list = rank_sort_term(canon_index_term)
   
    ## SORT BY NAME
    canon_term = name_sort_term(ranked_term, rank_list)

    return canon_term, scale_factor


def canonicalize_RDM(sqa_term, transRDM):

    # Find path rank of indices in term
    path_rank  = assign_path_rank(sqa_term)

    # Check all tensors in term to find RDM
    for t_ind, t in enumerate(sqa_term.tensors):

        # Keep track of the path rank for each tensor
        loop_path_rank = path_rank[:len(t.indices)]

        ## IF TENSOR IS AN RDM ##
        if isinstance(t, creDesTensor):

            # Extract cre/des operators from creDesTensor
            rdm_ops = [op for op in t.ops]

            # Initialize count variables
            cre_ext = 0
            des_ext = 0
    
            # Count how many cre/des operators have external indices
            for op, rank in zip(rdm_ops, loop_path_rank):
                if isinstance(op, creOp) and rank >= 50:
                    cre_ext += 1
                
                elif isinstance(op, desOp) and rank >= 50:
                    des_ext += 1
    
            # Reverse order of indices if there are more dummy des operators
            if (cre_ext < des_ext) and (not transRDM):
                rdm_ops.reverse()
                reversed_rdm_ops = []
    
                for op in rdm_ops:
                   if isinstance(op, desOp):
                       reversed_rdm_ops.append(creOp(op.indices))
    
                   elif isinstance(op, creOp):
                       reversed_rdm_ops.append(desOp(op.indices))
    
                sqa_term.tensors.pop(t_ind)
                sqa_term.tensors.append(creDesTensor(reversed_rdm_ops, transRDM))
    
        # Remove used path ranks elements
        path_rank = path_rank[len(t.indices):]

    return sqa_term


def assign_path_rank(sqa_term):

    # Store whether indices of tensors in term are being contracted over or not
    tensor_indices = []
    all_ind_list   = []

    for t in sqa_term.tensors:
        tensor_indices.append(''.join([i.name for i in t.indices]))
        all_ind_list.extend(''.join([i.name for i in t.indices]))

    index_dict = {}
    repeat_ind = range(50)
    unique_ind = range(50,100)

    for inds in tensor_indices:
        for char in inds:
            if (''.join(tensor_indices).count(char) == 1) and (len(inds) > 1):
                index_dict[char] = unique_ind[0]
                unique_ind.pop(0)

            # Define unique values for repeating indices
            else:
                if not index_dict.has_key(char):
                    index_dict[char] = repeat_ind[0]
                    repeat_ind.pop(0)

    path_rank_list = [index_dict[ind] for ind in all_ind_list]

    return path_rank_list


def assign_space_rank(sqa_term):

    space_rank_list = []

    # Store subspace information for all tensors in term
    space_indices = []

    for t in sqa_term.tensors:
         space_indices.extend(''.join([i.indType[0][0][0] for i in t.indices]))
 
    for i in space_indices:

        if i == 'c':
            value = 0
    
        elif i == 'a':
            value = 1
    
        elif i == 'v':
            value = 2

        else:
            raise Exception('Index does not have a valid orbital subspace') 

        space_rank_list.append(value * 1000)

    return space_rank_list


def check_intermediates(interm_list, int_term, int_tensor):

    isRedundant = False

    # Check every intermediate in the existing list
    for i, (list_term, list_tensor) in enumerate(interm_list):

        # Compare the amount of tensors in each term
        if len(list_term.tensors) == len(int_term.tensors):

           # Compare the names of the tensors that make up each term
           if [t.name for t in list_term.tensors] == [t.name for t in int_term.tensors]:

               # Last check is to compare the indices of the tensors
               list_p_rank = assign_path_rank(list_term)
               int_p_rank  = assign_path_rank(int_term)

               list_space_rank = assign_space_rank(list_term)
               int_space_rank  = assign_space_rank(int_term)

               # Check subspace of indices of intermediate tensor
               list_tensor_indices = [i.indType[0][0][0] for i in list_tensor.indices]
               int_tensor_indices  = [i.indType[0][0][0] for i in int_tensor.indices]

               if (list_p_rank == int_p_rank) and (list_space_rank == int_space_rank) and (list_tensor_indices == int_tensor_indices):
                
                   # Modify input term and return existing stored intermediate
                   isRedundant     = True
                   int_term        = list_term.copy()
                   int_tensor.name = list_tensor.name
                   break

    return int_term, int_tensor, isRedundant


def get_int_indices(sqa_tensor_list, ext_string):


    # Make copy of tensor list to modify
    new_tensor_list = sqa_tensor_list[:]
    
    # Create list for index objects that are used in INT tensor definition
    ext_ind_list  = []
    loop_ind_list = []

    # Iterate through tensors
    for tens_ind, t in enumerate(new_tensor_list):

        # Iterate through indices
        for ind_ind, i in enumerate(t.indices):

            # Find index objects in tensors that are external wrt intermediate definition
            if i.name in ext_string:

                # Create index list to define a tensor object wrt overall contraction
                if not loop_ind_list:
                    loop_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])

                # If there are indices in the list, ensure there are no duplicates
                else:
                    exist_ind = [ind.name for ind in loop_ind_list]

                    if i.name not in exist_ind:
                        loop_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])

                # Define index as an external index
                new_tensor_list[tens_ind].indices[ind_ind].isSummed = False

                # Create index list to define a tensor object wrt intermediate definition
                if not ext_ind_list:
                    ext_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])

                # If there are indices in the list, ensure there are no duplicates
                else:
                    exist_ind = [ind.name for ind in ext_ind_list]

                    if i.name not in exist_ind:
                        ext_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])

            # If the index is not external, ensure it is a dummy index
            else:
                new_tensor_list[tens_ind].indices[ind_ind].isSummed = True

    return new_tensor_list, ext_ind_list, loop_ind_list


def sort_terms(sorting_list, in_names, in_indices, in_symmetries):

    zipped_info  = zip(sorting_list, in_names, in_indices, in_symmetries)
    sorted_terms = sorted(zipped_info)

    # Create sorted lists
    out_names      = [name for _, name, _, _ in sorted_terms]
    out_indices    = [indices for _, _, indices, _ in sorted_terms]
    out_symmetries = [symmetries for _, _, _, symmetries in sorted_terms]

    return out_names, out_indices, out_symmetries    


def rank_sort_term(sqa_term):

    # Determine rank of tensors in term
    rank_list    = [len(t.indices) for t in sqa_term.tensors]

    if rank_list.count(rank_list[0]) == len(rank_list):
        rank_sorted_term = sqa_term.copy()

    else:    
        names       = [t.name for t in sqa_term.tensors]
        indices     = [t.indices for t in sqa_term.tensors]
        symmetries  = [t.symmetries for t in sqa_term.tensors]

        # Make unique value for each name
        rank_order = np.argsort(rank_list)[::-1]
   
        # Sort tensors in descending order by rank
        rank_names, rank_indices, rank_symmetries = sort_terms(rank_order, names, indices, symmetries)
    
        # Make list of sqa tensors to create new term sorted by rank
        rank_sorted_tensors = []
    
        for name, index, symmetry in zip(rank_names, rank_indices, rank_symmetries):
            rank_sorted_tensors.append(tensor(name, index, symmetry))
    
        rank_sorted_term = term(1.0, [], rank_sorted_tensors)

    return rank_sorted_term, rank_list


def name_sort_term(sqa_term, rank_list):
    
    names       = [t.name for t in sqa_term.tensors]
    indices     = [t.indices for t in sqa_term.tensors]
    symmetries  = [t.symmetries for t in sqa_term.tensors]

    # Make list of sqa tensors to create new term sorted by name within each rank
    sorted_tensors = []
    unique_ranks   = set(rank_list)

    # Preserve descending rank
    while unique_ranks:
        rank_names   = []
        rank_indices = []
        rank_syms    = []
        max_rank     = max(unique_ranks)

        for n, i, s in zip(names, indices, symmetries):
            if int(len(i)) == int(max_rank):
                rank_names.append(n)
                rank_indices.append(i)
                rank_syms.append(s)

        unicode_names = make_names(rank_names)

        # Sort tensors using unicode representation of their names
        uni_names, uni_indices, uni_symmetries = sort_terms(unicode_names, rank_names, rank_indices, rank_syms)
 
        for n, i, s in zip(uni_names, uni_indices, uni_symmetries):
            sorted_tensors.append(tensor(n, i, s))

        unique_ranks.remove(max_rank)

    name_sorted_term = term(1.0, [], sorted_tensors)

    return name_sorted_term


def make_names(name_list):

    # Make unique value for each name
    unicode_names = []
    for n in name_list: 
        p = 0

        for l in range(len(n)):
            p *= 255
            p += ord(n[l])
        
        unicode_names.append(p)

    return unicode_names
