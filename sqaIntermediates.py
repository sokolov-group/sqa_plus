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

import sys
import numpy as np
import opt_einsum as oe

from .sqaTensor import tensor, creOp, desOp, kroneckerDelta, creDesTensor
from .sqaTerm import term
from .sqaIndex import is_core_index_type, is_active_index_type, is_virtual_index_type, get_spatial_index_type
from .sqaOptions import options

def genIntermediates(input_terms, ind_str = None, custom_path = None):

    # Import options from sqaOptions class
    trans_rdm = options.genIntermediates.trans_rdm
    factor_depth = options.genIntermediates.factor_depth
    greedy_opt = options.genIntermediates.greedy

    if not input_terms:
        raise ValueError('No input terms provided for intermediate generation.') 

    if factor_depth < 0 or not isinstance(factor_depth, int):
        raise ValueError('Invalid factor depth, must use non-negative integer value.')

    options.print_header('Generating Intermediate Tensors | Factor Depth = {:}'.format(factor_depth))
    sys.stdout.flush()

    # Select optimizing path approach
    if greedy_opt:
        optimizer = 'random-greedy'
    else:
        optimizer = oe.DynamicProgramming(
            minimize='size',    # minimize largest intermediate tensor size
            search_outer=True,  # search through outer products as well
            cost_cap=False,     # don't use cost-capping strategy
        )

    # Create list of integers to form unique names of 'INT'    
    int_name_list = list(np.arange(1, 10000))

    # Initialize container lists
    intermediates = []      # intermediate tensors
    mod_term_list = []      # modified term list that will use intermediate tensors
    all_int_indices = []    # intermediate index string

    # Convert creOp/desOp objects to RDMs
    convert_credes_to_rdm(input_terms)

    # Iterate through every term in list of terms
    for _term_ind, _term in enumerate(input_terms, start=1):

        # Create lists for all tensors
        prefactor = _term.numConstant
        tensor_list = _term.tensors
        tensor_indices_list = [list(t.indices) for t in _term.tensors]

        # Build einsum strings
        lhs_str, einsum_string = build_einsum_string(tensor_indices_list, ind_str)

        # Compute intermediates if order of contracting tensors is specified
        if custom_path:

            # Check that requested contractions are in-range
            contract_order = [(i,j) for i,j in custom_path if max(i,j) < len(lhs_str)]
            if len(contract_order) < len(custom_path):
                print(f"WARNING: Not enough tensors for requested contraction, skipping term {_term_ind}.")

            # Make list of indices for all intermediates
            int_indices = []
            for i_ind, j_ind in contract_order:

                # Get indices from tensors being contracted
                contracted_inds = lhs_str[i_ind] + lhs_str[j_ind]
                
                # Keep only indices appearing once (not contracted over)
                int_ind = ''.join(i for i in contracted_inds if contracted_inds.count(i) == 1)
                int_indices.append(int_ind)
                
                # Remove contracted tensors from lhs_str and add intermediate
                lhs_str.pop(max(i_ind, j_ind))
                lhs_str.pop(min(i_ind, j_ind))

                lhs_str.append(int_ind)

        # If unspecified, use opt_einsum to compute contraction path
        else:

            # Build sizes dictionary
            sizes_dict = build_sizes_dict(tensor_indices_list)

            # Construct dummy tensors for term in order to assess contraction path
            dummy_tens = build_dummy_tensors(lhs_str, sizes_dict)

            # Compute most efficient contraction path
            path, path_info = oe.contract_path(einsum_string, *dummy_tens, optimize=optimizer)
            opt   = path_info.opt_cost
            naive = path_info.naive_cost

            # Append terms to modified term list if scaling of contraction cannot be optimized
            if opt >= naive:
                if options.verbose:
                    print(f"Contraction cannot be optimized, skipping term {_term_ind}.")
                mod_term_list.append(_term)
                continue

            if options.verbose:
                print(f'\n> Intermediate will be created for term {_term_ind}, reducing operations by {path_info.speedup:.2f}x.')

            # Save tuples that indicate optimized order of contracting tensors
            contract_order = [contraction for contraction in path[:factor_depth]]

            # Determine contraction path and indices of intermediates
            int_indices = [contraction[2].split('->')[1] for contraction in path_info.contraction_list][:factor_depth]

        # Define scale outside of loop
        scale_factor_total = 1.0

        # Create intermediate
        for num, contract in enumerate(contract_order):

            # Make intermediate name
            tensor_name = 'INT{:04d}'.format(int_name_list.pop(0))

            # Determine which tensors from tensor_list are being contracted
            tens_contract = [tensor_list[i] for i in contract]

            # Use the external string to modify indexType of the indices in tensors wrt the intermediate term
            new_tensors, def_indices, loop_indices = get_int_indices(tens_contract, int_indices[num])

            # Construct intermediate term w/ updated tensor
            int_term = term(1.0, [], new_tensors)

            # Canonicalize term and tensor representation of intermediate and update scale factor
            if options.verbose:
                print('------------------------------')
                print(f'CANONICALIZING {tensor_name}:')
                print(f'{int_term}')

            int_term, scale_factor = make_canonical(int_term, trans_rdm)
            scale_factor_total *= scale_factor

            # Update indices after canonicalizing term
            new_tensors, def_indices, loop_indices = get_int_indices(int_term.tensors, int_indices[num])

            # Define intermediate tensor wrt to new indices
            int_tensor = tensor(tensor_name, def_indices, [])

            # Check intermediates for redundancy
            if not intermediates:
                intermediates.append([int_term, int_tensor])
            else:
                int_term, int_tensor, isRedundant = check_intermediates(intermediates, int_term, int_tensor)
                if not isRedundant:
                    intermediates.append([int_term, int_tensor])

            # Modify 'tensor_list' for einsum's contract_path function
            tensor_list = [tens for tens in tensor_list if tens not in tens_contract]

            # Append representation of INT tensor with internal/external indices defined wrt to full contraction
            loop_tensor = tensor(int_tensor.name, loop_indices, [])
            tensor_list.append(loop_tensor)

        prefactor *= scale_factor_total
        mod_term_list.append(term(prefactor, [], tensor_list))

        # Append intermediate indices to 'all_int_indices'
        all_int_indices.append(int_indices[-1] if int_indices else [])

    if not intermediates:
        options.print_header("NO INTERMEDIATES WERE FOUND!")
        return input_terms, None

    options.genEinsum.keep_user_defined_dummy_names = True
    finalize_dummy_indices(mod_term_list, all_int_indices)
 
    renumber_intermediates(mod_term_list, intermediates)
    return mod_term_list, intermediates

def convert_credes_to_rdm(_terms_credes, trans_rdm = False):
    """Convert creation/destruction operators to RDMs."""
    for term_credes in _terms_credes:

        ## Append all cre/des operators to list
        credes_ops = [t for t in term_credes.tensors if isinstance(t, (creOp, desOp))]
        if not credes_ops:
            continue
    
        other_tensors = [t for t in term_credes.tensors if not isinstance(t, (creOp, desOp))]

        ## Modify term in list to use creDesTensor object instead of cre/des objects
        term_credes.tensors = other_tensors + [creDesTensor(credes_ops, trans_rdm)]

def build_einsum_string(tensor_indices, ind_str):
    """Build einsum strings."""
    # Append string of indices for left-hand side expression
    inputs = [''.join(i.name for i in ind_list) for ind_list in tensor_indices]

    # Set index string for output indices
    output = "" if not ind_str else ind_str

    # Make einsum string
    einsum_str = "{}->{}".format(",".join(inputs), output)

    return inputs, einsum_str
 
def build_sizes_dict(tensor_indices):
    """Build dictionary of index sizes."""
    # Iterate through lists of tensor indices
    sizes_dict = {}
    for ind_list in tensor_indices:
        for ind in ind_list:
            # Build dictionary of index sizes
            if ind.name in sizes_dict:
                continue

            if is_active_index_type(ind):
                sizes_dict[ind.name] = 2
            elif is_core_index_type(ind):
                sizes_dict[ind.name] = 4
            elif is_virtual_index_type(ind):
                sizes_dict[ind.name] = 6
            else:
                raise ValueError(f"Index {ind} does not belong to a valid orbital subspace")

    return sizes_dict

def build_dummy_tensors(lhs_str, sizes_dict):
    """Build dummy tensors for assessing contraction path."""
    dummy_tensors = []
    for tensor_indices in lhs_str:
        shape = tuple(sizes_dict[idx] for idx in tensor_indices)
        #dummy_tensors.append(np.empty(shape, dtype='f4'))
        dummy_tensors.append(np.random.rand(*shape))

    return dummy_tensors

def make_canonical(int_term, trans_rdm):
    """Canonicalize tensor indices in a term."""
    # Additional canonicalization for RDM tensors
    if not options.spin_adapted:
        if any(isinstance(t, creDesTensor) for t in int_term.tensors):
            int_term = canonicalize_rdm(int_term, trans_rdm)

    # Create ranking of indices based on the contraction path
    path_rank  = assign_path_rank(int_term)
    space_rank = assign_space_rank(int_term)

    # Store list of tensors that are in canonical order to create a canonicalized term
    canon_tensor_list = []
    prefactor = int_term.numConstant

    ### RE-ARRANGE THE INDICES WITHIN EACH TENSOR IN THE TERM
    # Continue canonicalizing after modifying RDM tensors
    tensor_start = 0
    for t in int_term.tensors:
        n_inds = len(t.indices)

        # Slice path and space ranks for the current tensor
        tensor_path = path_rank[tensor_start:tensor_start + n_inds]
        tensor_space = space_rank[tensor_start:tensor_start + n_inds]

        # Ensure normal-ordering for RDMs by prioritizing desOp in path rank
        if isinstance(t, creDesTensor):
           for i, op in enumerate(t.ops):
               if isinstance(op, desOp):
                   tensor_path[i] += 100
        
        # Generate permutation of indices to canonical order
        combined_rank = [path + space for path, space in zip(tensor_path, tensor_space)]
        canon_order = np.argsort(combined_rank).tolist()

        # Get symmetry information from sqa_tensor
        symPerms, symFactors  = t.symPermutes()

        # Only modify tensors with appropriate symmetry
        if canon_order in symPerms:

            # Get index of symmetry for convenience
            ind = symPerms.index(canon_order)

            # Keep track of prefactor
            prefactor *= float(symFactors[ind])

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

        # Update tensor start index for next tensor
        tensor_start += n_inds

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

    return canon_term, prefactor

def canonicalize_rdm(sqa_term, trans_rdm):
    """Canonicalize rdm tensor indices in a term."""
    # Find path rank of indices in term
    path_rank  = assign_path_rank(sqa_term)

    # Check all tensors in term to find RDM
    tensor_start = 0
    for t_ind, t in enumerate(sqa_term.tensors):
        n_inds = len(t.indices)

        # Skip processing non-RDM tensors
        if not isinstance(t, creDesTensor):
            tensor_start += n_inds
            continue

        # Keep track of the path rank for each tensor
        tensor_path_rank = path_rank[tensor_start:tensor_start + n_inds]

        # Extract cre/des operators from creDesTensor
        rdm_ops = t.ops

        # Count how many cre/des operators have external indices
        cre_ext = sum(isinstance(op, creOp) and rank >= 50 for op, rank in zip(rdm_ops, tensor_path_rank))
        des_ext = sum(isinstance(op, desOp) and rank >= 50 for op, rank in zip(rdm_ops, tensor_path_rank))

        # Reverse order of indices if there are more external des operators
        if (cre_ext < des_ext) and (not trans_rdm):
            reversed_rdm_ops = []

            for op in reversed(rdm_ops):
                new_op = creOp(op.indices) if isinstance(op, desOp) else desOp(op.indices)
                reversed_rdm_ops.append(new_op)

            sqa_term.tensors[t_ind] = creDesTensor(reversed_rdm_ops, trans_rdm)

        # Update tensor start index for next tensor
        tensor_start += n_inds

    return sqa_term


def assign_path_rank(sqa_term):
    """Assign rank of indices based on if external or internal(dummy)."""
    # Tensor indices
    tensor_indices = [''.join(i.name for i in t.indices) for t in sqa_term.tensors]
    all_indices = ''.join(tensor_indices)

    # Count occurrences of each index
    index_counts = {char: all_indices.count(char) for char in set(all_indices)}
    
    # Assign ranks: repeated indices get 0-49, unique indices get 50-99
    index_dict = {}
    repeat_rank = 0
    unique_rank = 50
    
    for inds in tensor_indices:
        for char in inds:
            if char not in index_dict:
                # Unique index
                if index_counts[char] == 1 and len(inds) > 1:
                    index_dict[char] = unique_rank
                    unique_rank += 1

                # Repeated index
                else:
                    index_dict[char] = repeat_rank
                    repeat_rank += 1
    
    return [index_dict[char] for char in all_indices]


def assign_space_rank(sqa_term):
    """Assign rank of indices based on subspace."""
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
    """Check for redundacies in intermediate tensors."""
    if options.verbose:
        print(f'\nCHECKING REDUNDANCY OF {int_tensor.name}...')

    int_path_rank = assign_path_rank(int_term)
    int_space_rank  = assign_space_rank(int_term)
    int_tensor_types = [get_spatial_index_type(ind.indType) for ind in int_tensor.indices]

    # Check every intermediate in the existing list
    for list_term, list_tensor in interm_list:

        # Compare the amount of tensors in each term
        if len(list_term.tensors) != len(int_term.tensors):
            continue

        # Compare the names of the tensors that make up each term
        if [t.name for t in list_term.tensors] != [t.name for t in int_term.tensors]:
            continue
    
        # Last check is to compare the indices of the tensors
        list_path_rank = assign_path_rank(list_term)
        list_space_rank = assign_space_rank(list_term)
        list_tensor_types = [get_spatial_index_type(ind.indType) for ind in list_tensor.indices]

        if (list_path_rank == int_path_rank and list_space_rank == int_space_rank and list_tensor_types == int_tensor_types):
            if options.verbose:
                print(f'REDUNDANCY FOUND. {int_tensor.name} IS EQUAL TO {list_tensor.name}.')

            # Modify input term and return existing stored intermediate
            int_term = list_term.copy()
            int_tensor.name = list_tensor.name
            return int_term, int_tensor, True

    if options.verbose:
        print('INTERMEDIATE IS UNIQUE.')

    return int_term, int_tensor, False

def get_int_indices(sqa_tensor_list, ext_string):
    """Set internal(dummy) and external indices for intermediate tensors."""
    # Make copy of tensor list to modify
    new_tensor_list = sqa_tensor_list[:]

    # Create list for index objects that are used in INT tensor definition
    ext_ind_list  = []
    loop_ind_list = []

    # Create sets to avoid duplicates
    ext_names = set()
    loop_names = set()

    # Iterate through tensor indices
    for t in new_tensor_list:
        for ind in t.indices:

            # External indices
            if ind.name in ext_string:
                ind.isSummed = False

                if ind.name not in ext_names:
                    ext_ind_list.append(ind)
                    ext_names.add(ind.name)

                if ind.name not in loop_names:
                    loop_ind_list.append(ind)
                    loop_names.add(ind.name)

            # Dummy indices
            else:
                ind.isSummed = True

    return new_tensor_list, ext_ind_list, loop_ind_list
 

def rank_sort_term(sqa_term):
    """Sort term tensors by rank."""
    # Determine rank of tensors in term
    rank_list = [len(t.indices) for t in sqa_term.tensors]

    if rank_list.count(rank_list[0]) == len(rank_list):
        rank_sorted_term = sqa_term.copy()

    else:
        # Make unique value for each name
        rank_order = np.argsort(rank_list)[::-1]
        rank_sorted_term = term(1.0, [], [sqa_term.tensors[i] for i in rank_order])

    return rank_sorted_term, rank_list


def name_sort_term(sqa_term, rank_list):
    """Sort term tensors by name."""
    # Make list of sqa tensors to create new term sorted by name within each rank
    sorted_tensors = []
    unique_ranks   = sorted(set(rank_list), reverse=True)

    for rank in unique_ranks:

        # Collect tensors with current rank
        rank_tensors = [tens for tens in sqa_term.tensors if len(tens.indices) == rank]

        # Sort by name
        names_to_sort = make_names([tens.name for tens in rank_tensors])
        name_sort     = np.argsort(names_to_sort)[::-1]

        sorted_tensors.extend([rank_tensors[i] for i in name_sort])

    return term(1.0, [], sorted_tensors)

def make_names(name_list):
    """Make unique value for each name."""
    from functools import reduce
    return [reduce(lambda p, c: p * 255 + ord(c), name, 0) for name in name_list]

def finalize_dummy_indices(term_list, indices_list):
    """Finalize intermediate indices in term list as user-defined and not summed."""
    for _term, int_indices_list in zip(term_list, indices_list):
        for _tensor in _term.tensors:
            for _index in _tensor.indices:
                index_name = _index.name
                if index_name in int_indices_list:
                    _index.isSummed = False
                    _index.userDefined = True

def renumber_intermediates(modified_terms, intermediate_terms):
    """ Reassign names of intermediate tensors to avoid numerical gaps."""
    # Intermediate name map
    name_map = {}

    # Assign new sequential names to intermediate tensors
    for index, (_term, _tensor) in enumerate(intermediate_terms, start=1):
        old_name = _tensor.name
        new_name = f"INT{index:02d}"
        name_map[old_name] = new_name
        _tensor.name = new_name

        # If intermediate term references other intermediates (factor_depth > 1), update those names too
        for inner_tensor in _term.tensors:
            if inner_tensor.name in name_map:
                #print(f"<<< Updating Term-Tensor Pair: {_term}, {_tensor}")
                inner_tensor.name = name_map[inner_tensor.name]
                #print(f">>> Updated Term-Tensor Pair: {_term}, {_tensor}")

    # Apply renaming to intermediates in modified_terms
    for _term in modified_terms:
        for _tensor in _term.tensors:
            if _tensor.name in name_map:
                _tensor.name = name_map[_tensor.name]

