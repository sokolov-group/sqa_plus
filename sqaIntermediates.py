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
#

import numpy as np

from sqaTensor import tensor, creOp, desOp, kroneckerDelta, creDesTensor
from sqaTerm import term
from sqaIndex import is_core_index_type, is_active_index_type, is_virtual_index_type

from sqaOptions import options

def genIntermediates(input_terms, ind_str = None, custom_path = None):

    # Import options from sqaOptions class
    trans_rdm = options.genIntermediates.trans_rdm
    factor_depth = options.genIntermediates.factor_depth

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
            rdm_tensor = creDesTensor(credes_list, trans_rdm)
            tensorlist.append(rdm_tensor)
            tensor_indices_list.append([ind for ind in rdm_tensor.indices])

        # Create einsum string and dictionary of index sizes
        lhs_str = []
        sizes_dict = dict()

        # Loop through lists of indices for each tensor in term
        for ind_list in tensor_indices_list:

            # Make LHS string
            lhs_str.append(''.join([i.name for i in ind_list]))

            # Define lengths of unique indices
            for i in ind_list:

                # Only add new indices to dictionary of index sizes
                if i.name not in sizes_dict.keys():

                    size = None

                    # Weigh size of index by subspace
                    if is_active_index_type(i):
                        size = 2
                    elif is_core_index_type(i):
                        size = 4
                    else:
                        size = 6

                    # Add key/value pair
                    sizes_dict[i.name] = size

        # Make einsum string
        einsum_string = str(','.join(lhs_str) + '->' + ind_str)

        # Construct dummy tensors for term in order to assess contraction path
        dummy_tens = []

        for t in tensorlist:
            dims = []

            for index in t.indices:
                dims.append(sizes_dict[index.name])

            dummy_tens.append(np.empty(tuple(dims)))

        # Compute most efficient contraction path
        path_info = np.einsum_path(einsum_string, *dummy_tens)

        # TODO: WHICH EINSUM? OPT OR NUMPY?
        # Isolate intermediate contractions from opt_einsum info
        naive          = int(path_info[1].split('\n')[1].split()[-1])
        opt            = int(path_info[1].split('\n')[2].split()[-1])

        # Determine if opt_einsum will improve scaling of contraction
        if naive > opt:

            # If an order of contracting tensors is specified
            if custom_path:

                # Make copy of user-defined contraction order
                contract_order = custom_path[:]

                ################
                # Check that requested contractions are in-range
                for con, contract in enumerate(contract_order):

                    # Check that contraction path can be performed
                    i_ind, j_ind = contract

                    if j_ind >= len(lhs_str):
                        options.print_header("WARNING")
                        print("Not enough tensors for contraction: %s. Will be ignored..." % str(contract))
                        options.print_divider()
                        contract_order.pop(con)
                ################

                # Make list of indices for all intermediates
                int_indices = []

                # Get all indices involved in contraction
                for contract in contract_order:

                    tens_inds       = [lhs_str[i] for i in contract]
                    contracted_inds = ''.join([lhs_str[i] for i in contract])

                    # Construct string out of indices not contracted over
                    int_ind = ''

                    for i in contracted_inds:
                        if contracted_inds.count(i) == 1:
                            int_ind += i

                    # Append to list of intermediate indices
                    int_indices.append(int_ind)

                    # Modify lhs_string to include intermediate indices
                    lhs_str = [inds for inds in lhs_str if inds not in tens_inds]
                    lhs_str.append(int_ind)

            # Standard procedure for generating contraction path
            else:

                # Save tuples that indicate optimized order of contracting tensors
                contract_order = [contract for contract in path_info[0][1:1 + factor_depth]]

                # Determine contraction path and indices of intermediates
                split_path     = path_info[1].split('\n')[10:10 + factor_depth]
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
                    interm_name_list = interm_name_list[1:]

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

    # Form term with tensors w/ canonicalized indices
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
