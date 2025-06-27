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

import numpy as np

from sqaTensor import tensor, creOp, desOp, kroneckerDelta, creDesTensor
from sqaTerm import term
from sqaIndex import is_core_index_type, is_active_index_type, is_virtual_index_type, get_spatial_index_type

from sqaOptions import options

from sqaSpinAdapted import convert_credes_to_rdm 

def genIntermediates(input_terms, ind_str = None, custom_path = None):
    """
    Generate Intermediate Terms for Tensor Rank Reduction.

    This function analyzes a list of input tensor contraction terms and generates intermediate tensors
    to reduce computational cost. Using either a user-specified contraction order or NumPy's einsum_path, intermediates
    are formed up to the specified factorization depth.

    Args:
        input_terms (list): List of term objects representing tensor contractions.
        ind_str (str, optional): String specifying the output indices for the contraction.
        custom_path (list, optional): User-defined contraction order as a list of tensor index pairs.

    Returns:
        tuple:
            - modified_term_list (list): List of terms with intermediates inserted.
            - intermediates (list): List of unique intermediate terms and their tensor representations.
    """

    # Import options from sqaOptions class
    trans_rdm = options.genIntermediates.trans_rdm
    factor_depth = options.genIntermediates.factor_depth
    opt_einsum = options.genIntermediates.opt_einsum

    if len(input_terms) == 0:
        raise Exception('List of input terms is empty.')

    if opt_einsum:
        import opt_einsum as oe

    if factor_depth < 0 or not isinstance(factor_depth, int):
        raise ValueError('Invalid factor depth provided -- provide non-negative integer value.')

    startTime = time.time()
    options.print_header('Generating Intermediate Tensors | Factor Depth = {:}'.format(factor_depth))
    options.print_divider()
    sys.stdout.flush()

    # Create list of integers to form unique names of 'INT'
    int_name_list = list(np.arange(1, 10000))

    # Initialize lists
    intermediates = []      # intermediate tensors
    all_int_indices = []    # intermediate index string
    modified_term_list = [] # modified input terms that will use intermediate tensors
    term_data = []          # data for each input term to be processed

    # Convert Cre/Des Objects to RDM Objects
    convert_credes_to_rdm(input_terms, trans_rdm)
    
    # Iterate through input terms to prepare data for contraction
    for _term in input_terms:

        # Create lists for tensors
        tensorList = list(_term.tensors)
        tensorIndicesList = [list(t.indices) for t in _term.tensors]
        prefactor = _term.numConstant

        # Build einsum string and sizes dictionary
        lhs_str, sizes_dict = build_einsum_string(tensorIndicesList)

        term_data.append({
            'term': _term,
            'tensorList': tensorList,
            'prefactor': prefactor,
            'lhs_str': lhs_str,
            'sizes_dict': sizes_dict
        })

    # Create list of dummy tensors for assessing contraction path
    dummy_tensor_list = create_dummy_tensors([data['lhs_str'] for data in term_data], [data['sizes_dict'] for data in term_data])

    # Generate intermediates
    for term_idx, (term_info, dummy_tens) in enumerate(zip(term_data, dummy_tensor_list)):
        _term = term_info['term']
        tensorList = term_info['tensorList']
        prefactor = term_info['prefactor']
        lhs_str = term_info['lhs_str']
        sizes_dict = term_info['sizes_dict']

        # Make einsum string
        ind_str = ind_str or ""
        einsum_string = str(','.join(lhs_str) + '->' + ind_str)

        # Compute most efficient contraction path

        # TODO: WHICH EINSUM? OPT OR NUMPY?
        # Determine if opt_einsum will improve scaling of contraction
        if opt_einsum:
            # using FLOP count to compare
            path_info = oe.contract_path(einsum_string, *dummy_tens, optimize="greedy")
            naive     = path_info[1].naive_cost
            opt       = path_info[1].opt_cost
        else:
            # using scaling to compare
            path_info = np.einsum_path(einsum_string, *dummy_tens, optimize="greedy")
            naive     = int(path_info[1].split('\n')[1].split()[-1])
            opt       = int(path_info[1].split('\n')[2].split()[-1])

        # Print contraction path information
        if options.verbose:
            print('Einsum Contraction Path:')
            print(path_info[1])
            print('Naive Cost: {:}, Optimized Cost: {:}'.format(naive, opt))
            print('')

        # Append terms to modified term list if scaling of contraction cannot be optimized
        if naive <= opt:
            modified_term_list.append(_term)
            continue

        ### DEBUG: Show warning if factor_depth is larger than the number of tensors
        contract_limit = len(path_info[0][1:])
        if contract_limit <= factor_depth:
            print('WARN: The factor_depth requested may produce erroneous intermediates. Proceed with caution...')
        ### DEBUG END

        ### GENERATE CONTRACTION PATH
        # Check if contraction order has been specified
        if custom_path:

            # Initialize intermediates indices list
            int_indices = []

            # Set user-defined path to be the contraction order
            contract_order = custom_path[:]

            # Check that requested contractions are in-range
            for contract in contract_order:

                # Check that contraction path can be performed, skip if it can't
                i_ind, j_ind = contract
                if j_ind >= len(lhs_str):
                    options.print_header("WARNING")
                    print("Not enough tensors for contraction: %s. It will be ignored..." % str(contract))
                    options.print_divider()
                    continue

                # Make list of indices for intermediates
                tens_inds = [lhs_str[i] for i in contract]
                contracted_inds = ''.join(tens_inds)

                # Create string of unique indices 
                int_idx = ''.join(i for i in contracted_inds if contracted_inds.count(i) == 1)

                # Append to list of intermediate indices
                int_indices.append(int_idx)

                # Update lhs_string to include intermediate indices and remove contracted ones
                lhs_str = [s for s in lhs_str if s not in tens_inds] + [int_idx]

        # Use standard contraction path if none has been specified
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

        # Define pre-factor outside of loop
        scale_factor_total = 1.0

        ### FORM INTERMEDIATES
        for num, contract in enumerate(contract_order):

            # Make intermediate name
            tensor_name = 'INT{:02d}'.format(int_name_list.pop(0))

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

            # Append intermediate term
            if not intermediates:
                if options.verbose:
                    print(int_tensor.name)
                    print(int_term)
                    print('')
                
                intermediates.append([int_term, int_tensor])

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

            # Modify 'tensorList' for einsum's contract_path function
            tensorList = [tens for tens in tensorList if tens not in tens_contract]

            # Store INT tensor to 'tensorList'
            tensorList.append(tensor(int_tensor.name, loop_indices, []))

        # Scale terms and store modified einsum expression
        prefactor *= scale_factor_total
        modified_term_list.append(term(prefactor, [], tensorList))

        # Append intermediate indices to 'all_int_indices'
        all_int_indices.append(int_indices[-1])

    if len(intermediates) == 0:
        options.print_header("NO INTERMEDIATES WERE FOUND!")
        return input_terms, None
 
    ## MAKE INDICES EXTERNAL IN MODIFIED TERM LIST
    options.genEinsum.keep_user_defined_dummy_names = True
    for tensor_term, int_indices in zip(modified_term_list, all_int_indices):
        for _tensor in tensor_term.tensors:
            for _index in _tensor.indices:

                index_name = _index.name

                if any(index_name == idx for idx in int_indices):
                    _index.isSummed = False 
                    _index.userDefined = True

    ## RENUMBERING FOR INTERMEDIATES
    renumber_intermediates(modified_term_list, intermediates)

    print("Total intermediates generated: {:}".format(len(intermediates)))
    print("Total modified terms: {:}".format(len(modified_term_list)))
    print("Intermediate generation time :  {:.3f} seconds".format(time.time() - startTime))
    options.print_divider()
    sys.stdout.flush()
    return modified_term_list, intermediates

##TODO: combine rank_sort_term & name_sort_term into one function
def make_canonical(int_term, trans_rdm):

    # Canonicalize any RDM tensors present
    if options.spin_orbital:
        if any(isinstance(t, creDesTensor) for t in int_term.tensors):
            int_term = canonicalize_rdm(int_term, trans_rdm)

    # Create ranking of indices based on the contraction path
    path_rank  = assign_path_rank(int_term)
    space_rank = assign_space_rank(int_term)

    # Store list of tensors that are in canonical order to create a canonicalized term
    canon_tensor_list = []
    scale_factor      = int_term.numConstant

    ### RE-ARRANGE THE INDICES WITHIN EACH TENSOR IN THE TERM
    # Continue canonicalizing after modifying RDM tensors
    for t in int_term.tensors:
        n_inds = len(t.indices)

        # Slice rank info for this tensor
        loop_path_rank  = path_rank[:n_inds]
        loop_space_rank = space_rank[:n_inds]

        # Ensure normal-ordering for RDMs by prioritizing desOp in path rank
        if isinstance(t, creDesTensor):
           for i, op in enumerate(t.ops):
               if isinstance(op, desOp):
                   loop_path_rank[i] += 100

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
    
    # Sort terms by rank and then name

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
        n_inds = len(t.indices)

        # Skip processing non-RDM tensors
        if not isinstance(t, creDesTensor):
            path_rank = path_rank[n_inds:]
            continue

        # Keep track of the path rank for each tensor
        loop_path_rank = path_rank[:n_inds]

        # Extract cre/des operators from creDesTensor
        rdm_ops = t.ops

        # Count how many cre/des operators have external indices
        cre_ext = sum(isinstance(op, creOp) and rank >= 50 for op, rank in zip(rdm_ops, loop_path_rank))
        des_ext = sum(isinstance(op, desOp) and rank >= 50 for op, rank in zip(rdm_ops, loop_path_rank))

        # Reverse order of indices if there are more dummy des operators
        if (cre_ext < des_ext) and (not trans_rdm):
            reversed_rdm_ops = []

            for op in reversed(rdm_ops):
                new_op = creOp(op.indices) if isinstance(op, desOp) else desOp(op.indices)
                reversed_rdm_ops.append(new_op)

            sqa_term.tensors[t_ind] = creDesTensor(reversed_rdm_ops, trans_rdm)

        # Remove used path ranks elements
        path_rank = path_rank[n_inds:]

    return sqa_term


def assign_path_rank(sqa_term):

    # Store whether indices of tensors in term are being contracted over or not
    tensor_indices = []
    all_ind_list   = []

    for t in sqa_term.tensors:
        names = [i.name for i in t.indices]
        tensor_indices.append(''.join(names))
        all_ind_list.extend(''.join(names))

    index_dict = {}
    # Path values for repeating indices
    repeat_ind = list(range(50)) 
    # Path values for unique indices 
    unique_ind = list(range(50,100))

    for inds in tensor_indices:
        for char in inds:

            # Define condition for unique indices
            is_unique = ''.join(tensor_indices).count(char) == 1 and len(inds) > 1

            # Assign path rank for unique indices
            if is_unique:
                if char not in index_dict: 
                    index_dict[char] = unique_ind.pop(0)

            # Assign path rank for repeating indices
            else:
                if char not in index_dict: 
                    index_dict[char] = repeat_ind.pop(0)

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

##TODO: implement contraction signature to avoid false redundancies produced by INT terms when factor_depth > 1
def check_intermediates(interm_list, int_term, int_tensor):

    # Set flag for redundancy check
    isRedundant = False

    if options.verbose:
        print('CHECKING ' + str(int_tensor.name) + '...')

    # Check every intermediate in the existing list
    for list_term, list_tensor in interm_list:

        # Compare the amount of tensors in each term
        if len(list_term.tensors) != len(int_term.tensors):
            continue

        # Compare the tensors via canonical einsum strings
        lhs1 = [''.join(i.name for i in t.indices) for t in list_term.tensors]
        lhs2 = [''.join(i.name for i in t.indices) for t in int_term.tensors]

        rhs1 = ''.join(i.name for i in list_tensor.indices)
        rhs2 = ''.join(i.name for i in int_tensor.indices)
        
        list_canon = canonicalize_einsum_indices(lhs1, rhs1)
        int_canon = canonicalize_einsum_indices(lhs2, rhs2)
        
        if list_canon != int_canon:
            continue

        # Compare the ranks of the tensors that make up each term
        same_path_rank = assign_path_rank(list_term) == assign_path_rank(int_term)
        same_space_rank = assign_space_rank(list_term) == assign_space_rank(int_term)

        # Check orbital subspaces of intermediate tensor
        int_spatial_types = [get_spatial_index_type(ind.indType) for ind in list_tensor.indices]
        tensor_spatial_types = [get_spatial_index_type(ind.indType) for ind in int_tensor.indices]

        same_indices = int_spatial_types == tensor_spatial_types
            
        if same_path_rank and same_space_rank and same_indices:
            isRedundant = True

            # Modify input term and return existing stored intermediate
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

    return new_tensor_list, ext_ind_list, loop_ind_list


def rank_sort_term(sqa_term):

    # Determine rank of tensors in term
    rank_list    = [len(t.indices) for t in sqa_term.tensors]

    # Return copy if all tensors have same rank
    if rank_list.count(rank_list[0]) == len(rank_list):
        rank_sorted_term = sqa_term.copy()

    # Else, sort tensors by decreasing rank
    else:
        rank_order = np.argsort(rank_list)[::-1]
        sorted_tensors = [sqa_term.tensors[i] for i in rank_order]
        rank_sorted_term = term(1.0, [], sorted_tensors)

    return rank_sorted_term, rank_list


def name_sort_term(sqa_term, rank_list):

    # Make list of sqa tensors to create new term sorted by name within each rank
    sorted_tensors = []
    unique_ranks   = sorted(set(rank_list), reverse=True)

    # Iterate through unique ranks
    for rank in unique_ranks:

        # Collect tensors with current rank
        rank_tensors = [tens for tens in sqa_term.tensors if len(tens.indices) == rank]

        # Sort the tensors of current rank by name
        names_to_sort = make_names([tens.name for tens in rank_tensors])
        name_sort     = np.argsort(names_to_sort)[::-1]

        sorted_tensors.extend([rank_tensors[i] for i in name_sort])

    return term(1.0, [], sorted_tensors)


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

def renumber_intermediates(mod_term_list, int_term_list):

    # Initialize intermediate name map
    name_map = {}

    # Reassign INT names to avoid gaps in numbering
    for i, (_term, _tensor) in enumerate(int_term_list):
        new_name = 'INT{:02d}'.format(i+1)
        old_name = _tensor.name

        name_map[old_name] = new_name
        _tensor.name = new_name

    # Update tensor names in modified_term_list to use renumbered intermediates
    for _term in mod_term_list:
        for _tensor in _term.tensors:
            if _tensor.name in name_map:
                _tensor.name = name_map[_tensor.name]


def canonicalize_einsum_indices(lhs_list, rhs_str):

    canon_chars = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
    mapping = {}
    next_char_idx = 0

    # Process indices in lhs and rhs
    for idx_str in lhs_list + [rhs_str]:
        for char in idx_str:
            if char not in mapping:
                mapping[char] = canon_chars[next_char_idx]
                next_char_idx += 1

    # Apply mapping to create canonical forms
    lhs_canon = [''.join(mapping[char] for char in part) for part in lhs_list]
    rhs_canon = ''.join(mapping[char] for char in rhs_str)
    return (lhs_canon, rhs_canon)
   
def build_einsum_string(tensorIndicesList):

    # Initialize einsum string and dictionary of index sizes
    lhs_str = []
    sizes_dict = {}

    # Iterate through lists of indices lists
    for ind_list in tensorIndicesList:

        # Append string of indices for left-hand side expression
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
                    raise ValueError('Index does not belong to a valid orbital subspace')

    return lhs_str, sizes_dict

def create_dummy_tensors(lhs_str_list, sizes_dict_list):

    # Initialize list to hold all dummy tensors
    dummy_batch = []

    # Loop through lhs_str and sizes_dict pairs
    for lhs_str, sizes_dict in zip(lhs_str_list, sizes_dict_list):
        term_tensors = []

        # Create dummy tensors based on sizes_dict
        for tensor_indices in lhs_str:
            shape = tuple(sizes_dict[idx] for idx in tensor_indices)
            term_tensors.append(np.empty(shape, dtype='f4'))

        dummy_batch.append(term_tensors)

    return dummy_batch