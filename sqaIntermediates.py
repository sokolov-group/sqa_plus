# Copyright 2018-2022 SecondQuantizationAlgebra Developers. All Rights Reserved.
#
# Licensed under the GNU General Public License v3.0;
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Ilia Mazin <ilia.mazin@gmail.com>
#         Donna Odhiambo <donna.odhiambo@proton.me>
#

import sys, time
import itertools

import numpy as np

from sqaTensor import tensor, creOp, desOp, kroneckerDelta, creDesTensor
from sqaTerm import term
from sqaIndex import is_core_index_type, is_active_index_type, is_virtual_index_type, is_cvs_index_type

from sqaMatrixBlock import dummyLabel, reorder_tensor_indices
from sqaOptions import options

from sqaSpinAdapted import convert_credes_to_rdm 

def genIntermediates(input_terms, ind_str = None, custom_path = None):
    "Generate Intermediate Terms for Tensor Rank Reduction."

    # Import options from sqaOptions class
    trans_rdm = options.genIntermediates.trans_rdm
    factor_depth = options.genIntermediates.factor_depth
    opt_einsum = options.genIntermediates.opt_einsum

    if opt_einsum:
        import opt_einsum as oe

    if factor_depth < 0 or not isinstance(factor_depth, int):
        raise ValueError('Invalid factor depth provided -- provide integer value that is >= 0')

    startTime = time.time()
    options.print_header('Generating Intermediate Tensors | Factor Depth = {:}'.format(factor_depth))
    sys.stdout.flush()

    # Make list of integers to append to 'INT' string below
    int_name_list = np.arange(1, 10000)

    # Initialize list for modified input terms that will use intermediate tensors
    modified_term_list = []

    # Initialize list for intermediate tensors
    intermediates = []

    # Convert Cre/Des Objects to RDM Objects
    convert_credes_to_rdm(input_terms, trans_rdm) 
    
    # Iterate through every term in list of terms
    for _term in input_terms:

        # Create lists for tensors
        tensorList = list(_term.tensors)
        tensorIndicesList = [list(t.indices) for t in _term.tensors]
        prefactor = _term.numConstant

        # Create einsum string and dictionary of index sizes
        lhs_str = []
        sizes_dict = {}

        # Loop through lists of indices for each tensor in term
        for ind_list in tensorIndicesList:

            # Make string of indices for left-hand side expression
            lhs_str.append(''.join(i.name for i in ind_list))

            # Define lengths of unique indices
            for ind in ind_list:

                # Only add new indices to dictionary of index sizes
                if ind.name not in sizes_dict:
                    # Weigh size of index by subspace
                    if is_active_index_type(ind):
                        sizes_dict[ind.name] = 2
                    elif is_core_index_type(ind): 
                        sizes_dict[ind.name] = 4
                    elif is_virtual_index_type(ind):
                        sizes_dict[ind.name] = 6
                    else:
                        raise Exception('Index does not belong to a valid orbital subspace')

        # Make einsum string
        einsum_string = str(','.join(lhs_str) + '->' + ind_str)

        # Construct dummy tensors for term in order to assess contraction path
        #dummy_tens = [np.empty(tuple(sizes_dict[idx.name] for idx in t.indices)) for t in tensorList]
        dummy_tens = [np.random.rand(*(sizes_dict[idx.name] for idx in t.indices)) for t in tensorList]

        # Compute most efficient contraction path

        # TODO: WHICH EINSUM? OPT OR NUMPY?
        # Determine if opt_einsum will improve scaling of contraction
        # Isolate intermediate contractions from opt_einsum info
        if opt_einsum:
            # using FLOP count to compare
            path_info = oe.contract_path(einsum_string, *dummy_tens, optimize="greedy")
            naive     = path_info[1].naive_cost
            opt       = path_info[1].opt_cost
        else:
            # using scaling to compare
            path_info = np.einsum_path(einsum_string, *dummy_tens)
            naive     = int(path_info[1].split('\n')[1].split()[-1])
            opt       = int(path_info[1].split('\n')[2].split()[-1])

        if naive > opt:

            # If an order of contracting tensors is specified
            if custom_path:

                # Make copy of user-defined contraction order
                contract_order = custom_path[:]

                ################
                # Check that requested contractions are in-range
                valid_order = []

                for contract in contract_order:

                    # Check that contraction path can be performed
                    i_ind, j_ind = contract

                    if j_ind >= len(lhs_str):
                        options.print_header("WARNING")
                        print("Not enough tensors for contraction: %s. It will be ignored..." % str(contract))
                        options.print_divider()
                    else:
                        valid_order.append(contract)

                contract_order = valid_order
                ################

                # Make list of indices for all intermediates
                int_indices = []

                # Get all indices involved in contraction
                for contract in contract_order:

                    tens_inds       = [lhs_str[i] for i in contract]
                    #contracted_inds = ''.join(lhs_str[i] for i in contract)
                    contracted_inds = ''.join(tens_inds)

                    # Construct string out of indices that appear once (not contracted over)
#                    int_ind = ''
#
#                    for i in contracted_inds:
#                        if contracted_inds.count(i) == 1:
#                            int_ind += i
                    int_ind = ''.join(i for i in contracted_inds if contracted_inds.count(i) == 1)

                    # Append to list of intermediate indices
                    int_indices.append(int_ind)

                    # Modify lhs_string to include intermediate indices and removed contracted ones
#                    lhs_str = [inds for inds in lhs_str if inds not in tens_inds]
                    for i in sorted(contract, reverse=True):
                        lhs_str.pop(i)
                    lhs_str.append(int_ind)

            # Standard procedure for generating contraction path
            else:

                # Save tuples that indicate optimized order of contracting tensors
                if opt_einsum:
                    contract_order = [contract for contract in path_info[0][0:0 + factor_depth]]
                else:
                    contract_order = [contract for contract in path_info[0][1:1 + factor_depth]]

                # Determine contraction path and indices of intermediates
                if opt_einsum:
                    all_int_indices = [inds[2].split('->')[1] for inds in path_info[1].contraction_list]
                    int_indices = all_int_indices[0:0 + factor_depth]
                else:
                    split_path     = path_info[1].split('\n')[10:10 + factor_depth]
                    int_indices    = [str(line).split()[1].split('->')[1] for line in split_path]

            # Define scale outside of loop
            scale_factor_total = 1.0

            # Create intermediate
            for num, contract in enumerate(contract_order):

                # Make intermediate name
                tensor_name = 'INT' + str(int_name_list[0])

                # Determine which tensors from tensorList are being contracted
                tens_contract = [tensorList[i] for i in contract]

                # Modify indexType of tensors to external/dummy based on the intermediate term
                new_tensors, def_indices, loop_indices = get_int_indices(tens_contract, int_indices[num])

                # Construct intermediate term with updated tensor
                int_term = term(1.0, [], new_tensors)

                # Canonicalize term and tensor representation of intermediate and update scale factor
                if options.verbose:
                    print(tensor_name)
                    print(int_term)
                    print('CANONICALIZING...')

                int_term, scale_factor = make_canonical(int_term, trans_rdm)
                scale_factor_total *= scale_factor

                # Update indices after canonicalizing term
                new_tensors, def_indices, loop_indices = get_int_indices(int_term.tensors, int_indices[num])

                # Define intermediate tensor wrt to definition
                int_tensor = tensor(tensor_name, def_indices, [])

                # Append the first intermediate automatically
                if not intermediates:
                    if options.verbose:
                        print(int_tensor.name)
                        print(int_term)
                        print('')
                    intermediates.append([int_term, int_tensor])
                    int_name_list = int_name_list[1:]

                # Once intermediates list is not empty, check all other intermediates for redundancy against the list
                else:
                    int_term, int_tensor, isRedundant = check_intermediates(intermediates, int_term, int_tensor)

                    # Only append unique intermediate terms to the list of intermediates
                    if not isRedundant:
                        if options.verbose:
                            print('FOUND UNIQUE INTERMEDIATE')
                            print(int_tensor.name)
                            print(int_term)
                            print('')

                        intermediates.append([int_term, int_tensor])
                        int_name_list = int_name_list[1:]

                # Modify 'tensorList' for einsum's contract_path function
                tensorList = [tens for tens in tensorList if tens not in tens_contract]

                # Append representation of INT tensor w/ dummy/external indices defined wrt full contraction
                loop_tensor = tensor(int_tensor.name, loop_indices, [])
                tensorList.append(loop_tensor)

            prefactor *= scale_factor_total
            modified_term_list.append(term(prefactor, [], tensorList))

        # Scaling of contraction cannot be optimized
        else:
            modified_term_list.append(_term)

#######
    ## MAKE INDICES EXTERNAL IN MODIFIED TERM LIST
    options.genEinsum.keep_user_defined_dummy_names = True
    for tensor_term in modified_term_list:
        for _tensor in tensor_term.tensors:
            for _index in _tensor.indices:

                index_name = _index.name
                mymap = map(lambda idx: index_name in idx, int_indices)

                if mymap[0]:
                    _index.isSummed = False 
                    _index.userDefined = True
#######

    print("\nTotal intermediates generated: {:}".format(len(intermediates)))
    print("Intermediate generation time :  {:.3f} seconds".format(time.time() - startTime))
    options.print_divider()
    sys.stdout.flush()
    return modified_term_list, intermediates

def make_canonical(int_term, trans_rdm):

    # Additional canonicalization for RDM tensors
    for tens in [ten for ten in int_term.tensors]:

        if isinstance(tens, creDesTensor):
            int_term = canonicalize_rdm(int_term, trans_rdm)
            break

    # Create ranking of indices based on the contraction path
    path_rank  = assign_path_rank(int_term)
    space_rank = assign_space_rank(int_term)

    # Store list of tensors that are in canonical order to create a canonicalized term
    canon_tensor_list = []
    scale_factor      = int_term.numConstant

    ### RE-ARRANGE THE INDICES WITHIN EACH TENSOR IN THE TERM
    # Continue canonicalizing after modifying RDM tensors
#    for t_ind, t in enumerate(int_term.tensors):
    for t in int_term.tensors:
        n_inds = len(t.indices)

        # Slice rank info for this tensor
        loop_path_rank  = path_rank[:n_inds]
        loop_space_rank = space_rank[:n_inds]

##      TODO: figure out what this block does?
##        # Modify path rank for RDM to prioritize destruction operators
##        if isinstance(t, creDesTensor):
##           for i, op in enumerate(t.ops):
##               if isinstance(op, desOp):
##                   loop_path_rank[i] += 100

        # Create final ranking for tensor by combining path & space rank
        final_rank = [path + space for path, space in zip(loop_path_rank, loop_space_rank)]

        # Update ranks
        path_rank  = path_rank[n_inds:]
        space_rank = space_rank[n_inds:]

        # Canonical order by rank
        canon_order = np.argsort(final_rank).tolist()

        # Get symmetry information from sqa_tensor
        symPerms, symFactors  = t.symPermutes()

        # Only modify tensors with appropriate symmetry
        if canon_order in symPerms:

            # Get index of symmetry for convenience
            symInd = symPerms.index(canon_order)

            # Keep track of prefactor
            scale_factor  *= float(symFactors[symInd])

            # If modifying RDM tensor, create the sorted tensor correctly
            if isinstance(t, creDesTensor):
                sorted_ops    = [t.ops[i] for i in canon_order]
                sorted_tensor = creDesTensor(sorted_ops, trans_rdm)
            elif isinstance(t, kroneckerDelta):
                sorted_indices = [t.indices[i] for i in canon_order]
                sorted_tensor  = kroneckerDelta(sorted_indices)
            else:
                sorted_indices = [t.indices[i] for i in canon_order]
                sorted_tensor  = tensor(t.name, sorted_indices, t.symmetries)

            # Append tensor w/ sorted indices
            canon_tensor_list.append(sorted_tensor)

        # Append tensor w/ unpermuted indices
        else:
            canon_tensor_list.append(t)

    # Rebuild term with tensors w/ canonicalized indices
    canon_index_term = term(1.0, [], canon_tensor_list)

    if options.verbose:
        print('----- WRT INDICES IN TENSORS OF TERM -----')
        print(canon_index_term)

    ### RE-ARRANGE TENSORS THAT MAKE UP THE TERM
    ## SORT BY RANK
    ranked_term, rank_list = rank_sort_term(canon_index_term)
    if options.verbose:
        print('----- WRT TENSORS IN TERM BY RANK -----')
        print(ranked_term)

    ## SORT BY NAME
    canon_term = name_sort_term(ranked_term, rank_list)
    if options.verbose:
        print('----- WRT TENSORS IN TERM BY NAME -----')
        print(canon_term)

    return canon_term, scale_factor


def canonicalize_rdm(sqa_term, trans_rdm):

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
            if (cre_ext < des_ext) and (not trans_rdm):
                rdm_ops.reverse()
                reversed_rdm_ops = []

                for op in rdm_ops:
                   if isinstance(op, desOp):
                       reversed_rdm_ops.append(creOp(op.indices))

                   elif isinstance(op, creOp):
                       reversed_rdm_ops.append(desOp(op.indices))

                sqa_term.tensors.pop(t_ind)
                sqa_term.tensors.append(creDesTensor(reversed_rdm_ops, trans_rdm))

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

    # Store subspace information for all tensors in term
    space_rank_list = []

    for t in sqa_term.tensors:
        for i in t.indices:
            if is_core_index_type(i):
                space_rank_list.append(0)
            elif is_active_index_type(i):
                space_rank_list.append(1000)
            elif is_virtual_index_type(i):
                space_rank_list.append(2000)
            else:
                raise Exception('Index does not have a valid orbital subspace')

    return space_rank_list


def check_intermediates(interm_list, int_term, int_tensor):

    isRedundant = False

    if options.verbose:
        print('CHECKING ' + str(int_tensor.name) + '...')

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
                list_tensor_indices = []
                for i_list in list_tensor.indices:
                    if is_core_index_type(i_list):
                        list_tensor_indices.append('c')
                    elif is_active_index_type(i_list):
                        list_tensor_indices.append('a')
                    elif is_virtual_index_type(i_list):
                        list_tensor_indices.append('v')

                int_tensor_indices = []
                for i_int in int_tensor.indices:
                    if is_core_index_type(i_int):
                        int_tensor_indices.append('c')
                    elif is_active_index_type(i_int):
                        int_tensor_indices.append('a')
                    elif is_virtual_index_type(i_int):
                        int_tensor_indices.append('v')

                if (list_p_rank == int_p_rank) and (list_space_rank == int_space_rank) and (list_tensor_indices == int_tensor_indices):

                    # Modify input term and return existing stored intermediate
                    isRedundant     = True
                    if options.verbose:
                        print('------------------------------')
                        print(str(int_tensor.name) + ' IS REDUNDANT. NEXT INT CHECKED WILL HAVE THE SAME NAME')
                        print('STORED INTERMEDIATE TERM')
                        print(list_tensor.name)
                        print(list_term)
                        print('------------------------------')
                        print('')

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

    # Create sets to avoid dublicates
    ext_names = set()
    loop_names = set()

    # Iterate through tensors
    for t in new_tensor_list:
        #print("{:}".format(t))

        # Iterate through indices
        for ind in t.indices:

            # Define index as an external index
            if ind.name in ext_string:
                ind.isSummed = False

                # ext_ind_list: holds tensor object wrt intermediate definition
                if ind.name not in ext_names:
                    ext_ind_list.append(ind)
                    ext_names.add(ind.name)

                # loop_ind_list: holds tensor object wrt overall contraction
                if ind.name not in loop_names:
                    loop_ind_list.append(ind)
                    loop_names.add(ind.name)

            # Define index as a dummy index
            else:
                ind.isSummed = True

#    # Iterate through tensors
#    for tens_ind, t in enumerate(new_tensor_list):
#
#        # Iterate through indices
#        for ind_ind, i in enumerate(t.indices):
#
#            # Find index objects in tensors that are external wrt intermediate definition
#            if i.name in ext_string:
#
#                # Create index list to define a tensor object wrt overall contraction
#                if not loop_ind_list:
#                    loop_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])
#
#                # If there are indices in the list, ensure there are no duplicates
#                else:
#                    exist_ind = [ind.name for ind in loop_ind_list]
#
#                    if i.name not in exist_ind:
#                        loop_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])
#
#                # Define index as an external index
#                new_tensor_list[tens_ind].indices[ind_ind].isSummed = False
#
#                # Create index list to define a tensor object wrt intermediate definition
#                if not ext_ind_list:
#                    ext_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])
#
#                # If there are indices in the list, ensure there are no duplicates
#                else:
#                    exist_ind = [ind.name for ind in ext_ind_list]
#
#                    if i.name not in exist_ind:
#                        ext_ind_list.append(new_tensor_list[tens_ind].indices[ind_ind])
#
#            # If the index is not external, ensure it is a dummy index
#            else:
#                new_tensor_list[tens_ind].indices[ind_ind].isSummed = True

    return new_tensor_list, ext_ind_list, loop_ind_list


def rank_sort_term(sqa_term):

    # Determine rank of tensors in term
    rank_list    = [len(t.indices) for t in sqa_term.tensors]

    if rank_list.count(rank_list[0]) == len(rank_list):
        rank_sorted_term = sqa_term.copy()

    else:
        # Make unique value for each name
        rank_order = np.argsort(rank_list)[::-1]
        rank_sorted_term = term(1.0, [], [sqa_term.tensors[i] for i in rank_order])

    return rank_sorted_term, rank_list


def name_sort_term(sqa_term, rank_list):

    # Make list of sqa tensors to create new term sorted by name within each rank
    sorted_tensors = []
    unique_ranks   = set(rank_list)

    # Preserve descending rank
    while unique_ranks:

        # Make list for tensors with current maximum rank
        max_rank_tensors = []

        # Iterate through tensors
        for tens in sqa_term.tensors:

            # Append tensors with the largest rank to a list
            if len(tens.indices) == max(unique_ranks):
                max_rank_tensors.append(tens)

        # Sort the tensors of max rank by name
        names_to_sort = make_names([tens.name for tens in max_rank_tensors])
        name_sort     = np.argsort(names_to_sort)[::-1]

        max_rank_tensors = [max_rank_tensors[i] for i in name_sort]
        sorted_tensors.extend(max_rank_tensors)

        # Remove maximum rank in this loop
        unique_ranks.remove(max(unique_ranks))

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
