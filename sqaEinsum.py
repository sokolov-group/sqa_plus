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
# Author: Koushik Chatterjee <koushikchatterjee7@gmail.com>
#         Ilia Mazin <ilia.mazin@gmail.com>
#         Carlos E. V. de Moura <carlosevmoura@gmail.com>
#

from .sqaIndex import get_spin_index_type, \
                     is_core_index_type, is_active_index_type, is_virtual_index_type, \
                     is_cvs_core_index_type, is_cvs_valence_index_type, \
                     is_alpha_index_type, is_beta_index_type

from .sqaTensor import creOp, desOp, kroneckerDelta, creDesTensor
from .sqaMatrixBlock import dummyLabel
from .sqaOptions import options

from fractions import Fraction

def genEinsum(
    terms,
    lhs_string=None,
    indices_string=None,
    suffix=None,
    trans_indices_string=None,
    intermediate_list=None,
    help=False,
    **tensor_rename,
):
 
    if not terms:
        options.print_header("genEinsum equations")
        print('No terms provided for einsum equations.')
        options.print_divider()
        return

    # Use defaults if options not provided in arguments or set by sqaOptions class
    lhs_string = lhs_string or options.genEinsum.lhs_string or "temp"
    indices_string = indices_string or options.genEinsum.indices_string
    trans_indices_string = trans_indices_string or options.genEinsum.trans_indices_string
    intermediate_list = intermediate_list or options.genEinsum.intermediate_list
    suffix = suffix or options.genEinsum.suffix

    # spin-orbital suffix
    if not suffix and options.spin_orbital:
        suffix = "so"

    # Configuration options
    trans_rdm = options.genEinsum.trans_rdm
    remove_trans_rdm_constant = options.genEinsum.remove_trans_rdm_constant
    remove_core_integrals = options.genEinsum.remove_core_integrals

    opt_einsum_terms = options.genEinsum.opt_einsum_terms
    optimize = options.genEinsum.optimize

    keep_user_defined_dummy_names = options.genEinsum.keep_user_defined_dummy_names
    if not keep_user_defined_dummy_names:
        dummyLabel(terms, keep_user_defined_dummy_names)

    # Store custom names if provided by user
    custom_names = list(tensor_rename.items()) if tensor_rename else []

    ################################################
    # IF PROVIDED, PRINT EINSUMS FOR INT TERMS
    ################################################
    trans_int = []
    removed_int = []
    int_einsum_list = []

    if intermediate_list:

        # Process intermediates with transition RDM contractions
        if trans_rdm:
            trans_int = get_trans_intermediates(intermediate_list)

        # Remove double-counted contributions to core terms
        if remove_core_integrals:
            intermediate_list, removed_int = remove_core_int(intermediate_list, int_terms=True)

        # Generate einsum for each intermediate
        for int_term, int_tensor in intermediate_list:
            int_einsum = _build_intermediate_einsum(
                int_term,
                int_tensor,
                trans_indices_string,
                suffix,
                trans_int,
                custom_names,
                opt_einsum_terms,
                optimize,
            )
            int_einsum_list.append(int_einsum)

    ################################################
    # GENERATE EINSUM EXPRESSIONS FOR PROVIDED TERMS
    ################################################
    # Convert Cre/Des Objects to RDM Objects
    _convert_credes_to_rdm(terms, trans_rdm)

    # Constants terms in CAS blocks are removed by default, print warning
    if trans_rdm and remove_trans_rdm_constant:
        args = (terms, trans_int) if (intermediate_list and trans_int) else (terms,)
        terms, const_terms = remove_trans_rdm_const(*args)

    # If using effective Hamiltonian, remove double-counted contributions to core terms
    if remove_core_integrals:
        args = (terms, removed_int) if (intermediate_list and removed_int) else (terms,)
        terms, core_terms = remove_core_int(*args)
 
    # Generate einsum expressions for each term
    einsum_list = []

    for term_ind, term in enumerate(terms):
        einsum = _build_term_einsum(
            term,
            term_ind,
            lhs_string,
            trans_indices_string,
            indices_string,
            suffix,
            trans_int if (trans_rdm and intermediate_list) else None,
            custom_names,
            opt_einsum_terms,
            optimize,
        )
        einsum_list.append(einsum)

    if intermediate_list:
        options.print_header("genEinsum intermediates")
        for _einsum in int_einsum_list:
            print(_einsum)
        options.print_divider()

    options.print_header("genEinsum equations")
    for _einsum in einsum_list:
        print(_einsum)
    options.print_divider()

    # Modify return for intermediate term definition
    if intermediate_list:
        return int_einsum_list, einsum_list
    else:
        return einsum_list

def _build_intermediate_einsum(
    int_term,
    int_tensor,
    trans_indices_string,
    suffix,
    trans_int,
    custom_names,
    opt_einsum_terms,
    optimize,
):
    """Build einsum expression for an intermediate tensor."""

    # Get tensor information
    int_tensor_inds, int_tensor_names = get_tensor_info(
        int_term.tensors,
        trans_indices_string,
        "".join([i.name for i in int_tensor.indices]),
        suffix,
        trans_int,
        custom_names,
    )

    # Define tensor name for intermediate
    if custom_names and "INT" in [old for old, _ in custom_names]:
        tensor_name = make_custom_name(int_tensor, custom_names)
    else:
        tensor_name = int_tensor.name

    # Build einsum expression
    einsum_func = "einsum(" if opt_einsum_terms else "np.einsum("
    tensor_info = ", ".join([f"'{int_tensor_inds}'"] + int_tensor_names)

    # Add optimize flag
    if optimize:
        opt_flag = ", optimize = einsum_type)" if opt_einsum_terms else ", optimize = True)"
    else:
        opt_flag = ")"

    return f"{tensor_name} = {einsum_func}{tensor_info}{opt_flag}"

def _convert_credes_to_rdm(terms, trans_rdm):
    """Convert creOp/desOp objects to creDesTensor objects."""

    for term_credes in terms:

        # Append all cre/des operators to list
        credes_ops = [tens for tens in term_credes.tensors if isinstance(tens, (creOp, desOp))]

        # Modify term in list to use creDesTensor object instead of cre/des objects
        if credes_ops:
            other_tensors = [tens for tens in term_credes.tensors if tens not in credes_ops]
            term_credes.tensors = other_tensors + [creDesTensor(credes_ops, trans_rdm)]

def _build_term_einsum(
    term,
    term_ind,
    lhs_string,
    trans_indices_string,
    indices_string,
    suffix,
    trans_int,
    custom_names,
    opt_einsum_terms,
    optimize,
):
    """Build einsum expression for a single term."""

    # Set up equals sign for first term and rest of terms
    is_negative = term.numConstant < 0
    if term_ind == 0:
        assign_op = "=- " if is_negative else " = "
    else:
        assign_op = "-= " if is_negative else "+= "

    einsum = f"{lhs_string} {assign_op}"

    # Add scaling factor
    abs_constant = abs(term.numConstant)
    if round(abs_constant, 15) != 1.0:
        frac_constant = Fraction(abs_constant).limit_denominator()
        if round(float(frac_constant), 12) == round(abs_constant, 12):
            einsum += f"{frac_constant} * "
        else:
            einsum += f"{abs_constant} * "

    # Build einsum function call
    einsum_func = "einsum(" if opt_einsum_terms else "np.einsum("
    einsum += einsum_func

    # Get tensor information
    tensor_inds, tensor_names = get_tensor_info(
        term.tensors,
        trans_indices_string,
        indices_string,
        suffix,
        trans_int,
        custom_names,
    )

    # Add tensor information
    tensor_info = ", ".join([f"'{tensor_inds}'"] + tensor_names)
    einsum += tensor_info

    # Add optimize flag
    if optimize:
        opt_flag = ", optimize = einsum_type)" if opt_einsum_terms else ", optimize = True)"
    else:
        opt_flag = ")"

    # Append copy function call for single-tensor terms
    if len(term.tensors) == 1:
        opt_flag += '.copy()'

    return einsum + opt_flag

def get_tensor_info(
    sqa_tensors,
    trans_indices_string,
    indices_string,
    suffix,
    trans_int=None,
    custom_names=None
):

    # Import settings from options class
    spin_integrated_tensors = options.genEinsum.spin_integrated_tensors
    cvs_tensors = options.genEinsum.cvs_tensors

    # Make a list out of the user-provided external indices
    cvs_indices_list = options.genEinsum.cvs_indices_list
    cvs_indices_list = list(cvs_indices_list) if isinstance(cvs_indices_list, str) else cvs_indices_list

    val_indices_list = options.genEinsum.valence_indices_list
    val_indices_list = list(val_indices_list) if isinstance(val_indices_list, str) else val_indices_list

    # Process tensors
    tensor_names = []
    tensor_inds  = []

    for tens in sqa_tensors:

        tensor_name = None

        # Handle special case of kroneckerDelta (kdelta) object
        if isinstance(tens, kroneckerDelta):
            idx0, idx1 = tens.indices[0], tens.indices[1]
            idx0_name, idx1_name = idx0.name, idx1.name

            # Determine orbital space of kdelta
            if is_core_index_type(idx0) and is_core_index_type(idx1):

                if cvs_tensors:
                    if cvs_indices_list and val_indices_list:
                        orb_space = ('ncvs' if idx0_name in cvs_indices_list and idx1_name in cvs_indices_list else
                                     'nval' if idx0_name in val_indices_list and idx1_name in val_indices_list else
                                     'ncore' if (idx0_name not in cvs_indices_list and idx0_name not in val_indices_list and
                                                 idx1_name not in cvs_indices_list and idx1_name not in val_indices_list) else 'none')
                    elif cvs_indices_list:
                        orb_space = ('ncvs' if idx0_name in cvs_indices_list and idx1_name in cvs_indices_list else
                                     'ncore' if idx0_name not in cvs_indices_list and idx1_name not in cvs_indices_list else 'none')
                    elif val_indices_list:
                        orb_space = ('nval' if idx0_name in val_indices_list and idx1_name in val_indices_list else
                                     'ncore' if idx0_name not in val_indices_list and idx1_name not in val_indices_list else 'none')
                    else:
                        orb_space = 'ncvs' if is_cvs_core_index_type(idx0) else 'nval' if is_cvs_valence_index_type(idx0) else 'ncore'
                else:
                    orb_space = 'ncore'

            elif is_active_index_type(idx0) and is_active_index_type(idx1):
                orb_space = 'ncas'

            elif is_virtual_index_type(idx0) and is_virtual_index_type(idx1):
                orb_space = 'nextern'

            else:
                raise TypeError('WARNING: The indices of the kronecker delta term do not belong to the same orbital sub-space')

            orb_space = f"{orb_space}_{suffix}" if suffix else orb_space
            tensor_name = f"np.identity({orb_space})"

            if custom_names and "kdelta" in [old for old, _ in custom_names]:
                tensor_name = make_custom_name(tens, custom_names) + '_'

##                if cvs_tensors:
##                    if cvs_indices_list and val_indices_list:
##                        if (tens.indices[0].name in cvs_indices_list) and (tens.indices[1].name in cvs_indices_list):
##                            orb_space = 'ncvs'
##                        elif (tens.indices[0].name in val_indices_list) and (tens.indices[1].name in val_indices_list):
##                            orb_space = 'nval'
##                        elif (((tens.indices[0].name not in cvs_indices_list) and (tens.indices[0].name not in val_indices_list)) and
##                            ((tens.indices[1].name not in cvs_indices_list) and (tens.indices[1].name not in val_indices_list))):
##                            orb_space = 'ncore'
##                        else:
##                            orb_space = 'none'
##
##                    elif cvs_indices_list:
##                        if (tens.indices[0].name in cvs_indices_list) and (tens.indices[1].name in cvs_indices_list):
##                            orb_space = 'ncvs'
##                        elif (((tens.indices[0].name not in cvs_indices_list)) and ((tens.indices[1].name not in cvs_indices_list))):
##                            orb_space = 'ncore'
##                        else:
##                            orb_space = 'none'
##
##                    elif val_indices_list:
##                        if (tens.indices[0].name in val_indices_list) and (tens.indices[1].name in val_indices_list):
##                            orb_space = 'nval'
##                        elif (((tens.indices[0].name not in val_indices_list)) and ((tens.indices[1].name not in val_indices_list))):
##                            orb_space = 'ncore'
##                        else:
##                            orb_space = 'none'
##
##                    else:
##                        if is_cvs_core_index_type(tens.indices[0]):
##                            orb_space = 'ncvs'
##                        elif is_cvs_valence_index_type(tens.indices[0]):
##                            orb_space = 'nval'
##                        else:
##                            orb_space = 'ncore'
##                else:
##                    orb_space = 'ncore'
##
##            elif (is_active_index_type(tens.indices[0]) and is_active_index_type(tens.indices[1])):
##                orb_space = 'ncas'
##
##            elif (is_virtual_index_type(tens.indices[0]) and is_virtual_index_type(tens.indices[1])):
##                orb_space = 'nextern'
##
##            else:
##                raise TypeError('WARNING: The indices of the kronecker delta term do not belong to the same orbital sub-space')
##
##            if suffix:
##                orb_space += '_' + suffix
##            tensor_name += orb_space + ')'
##
##            # Rename if custom name is provided
##            if custom_names:
##                if ('kdelta') in [x for x,y in custom_names]:
##                    new_name = make_custom_name(tens, custom_names)
##                    tensor_name = new_name + '_'

        # Handle special case of orbital energy vector
        elif len(tens.indices) == 1 and tens.name.lower() == 'e':

            idx, idx_name = tens.indices[0], tens.indices[0].name

            # Determine orbital space
            if is_core_index_type(idx):

                if cvs_tensors:
                    if cvs_indices_list and val_indices_list:
                        orb_space = 'cvs' if idx_name in cvs_indices_list else 'val' if idx_name in val_indices_list else 'core'
                    elif cvs_indices_list:
                        orb_space = 'cvs' if idx_name in cvs_indices_list else 'core'
                    elif val_indices_list:
                        orb_space = 'val' if idx_name in val_indices_list else 'core'
                    else:
                        orb_space = 'cvs' if is_cvs_core_index_type(idx) else 'val' if is_cvs_valence_index_type(idx) else 'core'
                else:
                    orb_space = 'core'

            elif is_virtual_index_type(idx):
                orb_space = 'extern'

            else:
                orb_space = 'active'

            orb_space = f"{orb_space}_{suffix}" if suffix else orb_space
            tensor_name = f"{tens.name.lower()}_{orb_space}"
 
            if custom_names and 'e' in [old.lower() for old, _ in custom_names]:
                tensor_name = make_custom_name(tens, custom_names)

            ### Determine orbital space of energies
            ##if is_core_index_type(tens.indices[0]):
            ##    if cvs_tensors:
            ##        if cvs_indices_list and val_indices_list:
            ##            if tens.indices[0].name in cvs_indices_list:
            ##                orb_space = 'cvs'
            ##            elif tens.indices[0].name in val_indices_list:
            ##                orb_space = 'val'
            ##            else:
            ##                orb_space = 'core'
            ##        elif cvs_indices_list:
            ##            if tens.indices[0].name in cvs_indices_list:
            ##                orb_space = 'cvs'
            ##            else:
            ##                orb_space = 'core'
            ##        elif val_indices_list:
            ##            if tens.indices[0].name in val_indices_list:
            ##                orb_space = 'val'
            ##            else:
            ##                orb_space = 'core'
            ##        else:
            ##            if is_cvs_core_index_type(tens.indices[0]):
            ##                orb_space = 'cvs'
            ##            elif is_cvs_valence_index_type(tens.indices[0]):
            ##                orb_space = 'val'
            ##            else:
            ##                orb_space = 'core'
            ##    else:
            ##        orb_space = 'core'

            ##elif is_virtual_index_type(tens.indices[0]):
            ##    orb_space = 'extern'

            ##if suffix:
            ##    orb_space += '_' + suffix
            ##tensor_name += orb_space

            ### Rename if custom name is provided
            ##if custom_names:
            ##    if ('e' or 'E') in [x for x,y in custom_names]:
            ##        tensor_name = make_custom_name(tens, custom_names)

        # Handle special case of RDM tensor
        elif isinstance(tens, creDesTensor):

            # Modify name of RDM to reflect particle number
            number_suffix = ['c' if isinstance(op, creOp) else 'a' for op in tens.ops]
            tensor_name = tens.name + '_' + ''.join(number_suffix)

            # Append spin-integrated suffix if required
            if spin_integrated_tensors:
                spin_suffix = [('a' if is_alpha_index_type(idx) else 'b')
                             for idx in tens.indices if (is_alpha_index_type(idx) or is_beta_index_type(idx))]
                if spin_suffix:
                    tensor_name += '_' + ''.join(spin_suffix)

            tensor_name += f"_{suffix}" if suffix else ""

            if custom_names and 'rdm' in [old for old, _ in custom_names]:
                tensor_name = make_custom_name(tens, custom_names)

            ### Modify name of RDM to reflect particle number
            ##for op in tens.ops:
            ##    if isinstance(op, creOp):
            ##        tensor_name += 'c'
            ##    elif isinstance(op, desOp):
            ##        tensor_name += 'a'

            ### Append spin-integrated suffix if required
            ##if spin_integrated_tensors:
            ##    spin_suffix = '_'
            ##    for i in range(len(tens.indices)):
            ##        if is_alpha_index_type(tens.indices[i]):
            ##            spin_suffix += 'a'
            ##        elif is_beta_index_type(tens.indices[i]):
            ##            spin_suffix += 'b'
            ##    tensor_name += spin_suffix

            ### Append suffix
            ##if suffix:
            ##    tensor_name += '_' + suffix

            ### Rename if custom name is provided
            ##if custom_names:
            ##    if ('rdm') in [x for x,y in custom_names]:
            ##        tensor_name = make_custom_name(tens, custom_names)

        # Name integrals and amplitudes
        elif tens.name in ('h', 'v', 't1', 't2'):
            tensor_name = tens.name + '_'

            for idx in tens.indices:
                if is_active_index_type(idx):
                    tensor_name += 'a'
                elif is_core_index_type(idx):
                    idx_name = idx.name
                    if cvs_tensors:
                        if cvs_indices_list and val_indices_list:
                            tensor_name += ('x' if idx_name in cvs_indices_list else 
                                            'v' if idx_name in valence_indices_list else 'c')
                        elif cvs_indices_list:
                            tensor_name += 'x' if idx_name in cvs_indices_list else 'c'
                        elif val_indices_list:
                            tensor_name += 'v' if idx_name in valence_indices_list else 'c'
                        else:
                            tensor_name += ('x' if is_cvs_core_index_type(idx) else 
                                            'v' if is_cvs_valence_index_type(idx) else 'c')
                    else:
                        tensor_name += 'c'
                else:
                    tensor_name += 'e'

            # Append spin-integrated suffix if required
            if spin_integrated_tensors:
                spin_suffix = [('a' if is_alpha_index_type(idx) else 'b')
                             for idx in tens.indices if (is_alpha_index_type(idx) or is_beta_index_type(idx))]
                if spin_suffix:
                    tensor_name += '_' + ''.join(spin_suffix)

            # Add suffix for non-amplitude tensors
            if tens.name not in ('t1', 't2') and suffix:
                tensor_name += f"_{suffix}"

            # Rename if custom name is provided
            if custom_names:
                if tens.name in [old for old,_ in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

####
##            # Append letter representing orbital subspace of indices
##            for i in range(len(tens.indices)):
##                if is_active_index_type(tens.indices[i]):
##                    tensor_name += 'a'
##                elif is_core_index_type(tens.indices[i]):
##                    if cvs_tensors:
##                        if cvs_indices_list and val_indices_list:
##                            if (tens.indices[i].name in cvs_indices_list):
##                                tensor_name += 'x'
##                            elif (tens.indices[i].name in val_indices_list):
##                                tensor_name += 'v'
##                            else:
##                                tensor_name += 'c'
##                        elif cvs_indices_list:
##                            if (tens.indices[i].name in cvs_indices_list):
##                                tensor_name += 'x'
##                            else:
##                                tensor_name += 'c'
##                        elif val_indices_list:
##                            if (tens.indices[i].name in val_indices_list):
##                                tensor_name += 'v'
##                            else:
##                                tensor_name += 'c'
##                        else:
##                            if is_cvs_core_index_type(tens.indices[i]):
##                                tensor_name += 'x'
##                            elif is_cvs_valence_index_type(tens.indices[i]):
##                                tensor_name += 'v'
##                            else:
##                                tensor_name += 'c'
##                    else:
##                        tensor_name += 'c'
##                else:
##                    tensor_name += 'e'
##
##            # Append spin-integrated suffix if required
##            if spin_integrated_tensors:
##                spin_suffix = '_'
##                for i in range(len(tens.indices)):
##                    if is_alpha_index_type(tens.indices[i]):
##                        spin_suffix += 'a'
##                    elif is_beta_index_type(tens.indices[i]):
##                        spin_suffix += 'b'
##                tensor_name += spin_suffix
##
##            # Append suffix
##            if not (tens.name == 't1' or tens.name == 't2') and suffix:
##                tensor_name += '_' + suffix
##
##            # Rename if custom name is provided
##            if custom_names:
##                if tens.name in [x for x,y in custom_names]:
##                    tensor_name = make_custom_name(tens, custom_names)

        # Intermediate/custom tensors
        else:

            tensor_name = tens.name

            # Append spin-integrated suffix if required
            if spin_integrated_tensors:
                spin_suffix = [('a' if is_alpha_index_type(idx) else 'b')
                             for idx in tens.indices if (is_alpha_index_type(idx) or is_beta_index_type(idx))]
                if spin_suffix:
                    tensor_name += '_' + ''.join(spin_suffix)

            # Allow to rename intermediate tensors in term definitions
            if custom_names:
                if tens.name[:3] == 'INT' and 'INT' in [old for old,_ in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Create indices of tensor as string
        indices = ''.join(i.name for i in tens.indices)

        # Append 'slices' to appropiate dimensions of spin-integrated tensors
        if options.spin_integrated and not options.genEinsum.spin_integrated_tensors:
            tensor_name = append_spin_integrated_slice(tens, tensor_name, indices)

        # Append 'slices' to appropiate dimensions of tensors w/ CVS core indices
        if cvs_indices_list and not cvs_tensors:
            tensor_name = append_CVS_slice(tens, tensor_name, indices, suffix)

        # Append name of tensor (after and modifications due to special cases)
        tensor_names.append(tensor_name)

        # Append transition state index to appropriate set of indices
        if isinstance(tens, creDesTensor) and tens.trans_rdm:
            indices = trans_indices_string + indices
        elif trans_int and tens.name in trans_int:
            indices = trans_indices_string + indices

        # Append completed index string to list
        tensor_inds.append(indices)

    # Build complete index string with arrow notation
    tensor_inds = ','.join(tensor_inds)
    if trans_indices_string or indices_string:
        tensor_inds += '->' + (trans_indices_string or '') + (indices_string or '')
 
    ### Convert list of indices into one comma-separated string and prepare to append external index string
    ##tensor_inds = ','.join(tensor_inds)

    ### Check if rhs string or transition index is provided before adding arrow
    ##if trans_indices_string or indices_string:
    ##    tensor_inds += '->'

    ### Append transition index first, if present
    ##if trans_indices_string:
    ##    tensor_inds += trans_indices_string

    ### Append rhs string, if provided
    ##if indices_string:
    ##    tensor_inds += indices_string

    return tensor_inds, tensor_names

def remove_core_int(terms, removed_int = None, int_terms = False):

    # Remove terms from standard term list
    if not int_terms:
        options.print_header("WARNING")
        print('Terms with a contraction over repeating dummy core indices of 2e- integrals')
        print('will be removed. Set "remove_core_integrals" flag to FALSE to preserve terms')

        # Create lists to split up SQA terms
        kept_terms = []
        core_terms = []

        # Separate out the terms that have redundant 2e- integral contractions over core space
        for term_ind, term in enumerate(terms):
            coreTerm = False
            for tens_ind, tens in enumerate(term.tensors):
                if tens.name == 'v':
                    if (((options.physicists_notation) and
                         (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[2].name) or
                         (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[3].name) or
                         (terms[term_ind].tensors[tens_ind].indices[1].name) == (terms[term_ind].tensors[tens_ind].indices[2].name) or
                         (terms[term_ind].tensors[tens_ind].indices[1].name) == (terms[term_ind].tensors[tens_ind].indices[3].name))
                        or 
                        ((options.chemists_notation) and
                         (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[1].name) or
                         (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[3].name) or
                         (terms[term_ind].tensors[tens_ind].indices[2].name) == (terms[term_ind].tensors[tens_ind].indices[1].name) or
                         (terms[term_ind].tensors[tens_ind].indices[2].name) == (terms[term_ind].tensors[tens_ind].indices[3].name))
                       ):
                        coreTerm = True
                        break

                elif removed_int and (tens.name in removed_int):
                    coreTerm = True
                    break

            # Append to either list based on coreTerm flag
            if not coreTerm:
                kept_terms.append(terms[term_ind])

            else:
                core_terms.append(terms[term_ind])

        print('')
        print(str(len(core_terms)) + ' terms removed:')
        for term in core_terms:
            print(term)

        options.print_divider()
        print('Remaining terms: ' + str(len(kept_terms)))
        print('')

        return kept_terms, core_terms

    # Filter through intermediate definitions
    else:
        options.print_header("WARNING")
        print('Intermediate tensors defined w/ contractions over repeating dummy core indices of')
        print('2e- integrals will be removed. Set "remove_core_integrals" flag to FALSE to preserve definitions')

        # Track which tensor definitions are removed and kept
        removed_int   = []
        removed_terms = []

        # Determine which intermediate definitions to remove
        for int_ind, (int_term, int_tensor) in enumerate(terms):
            for tens in int_term.tensors:

                # If intermediate is defined w/ 2e- integral
                if tens.name == 'v':
                    if (((options.physicists_notation) and
                         (tens.indices[0].name) == (tens.indices[2].name) or
                         (tens.indices[0].name) == (tens.indices[3].name) or
                         (tens.indices[1].name) == (tens.indices[2].name) or
                         (tens.indices[1].name) == (tens.indices[3].name))
                        or 
                        ((options.chemists_notation) and
                         (tens.indices[0].name) == (tens.indices[1].name) or
                         (tens.indices[0].name) == (tens.indices[3].name) or
                         (tens.indices[2].name) == (tens.indices[1].name) or
                         (tens.indices[2].name) == (tens.indices[3].name))
                        ):
                        removed_int.append(int_tensor.name)
                        removed_terms.append(terms[int_ind])
                        break

                # If intermediate is defined in terms of one of the intermediates to be removed
                elif tens.name in removed_int:
                    removed_int.append(int_tensor.name)
                    removed_terms.append(terms[int_ind])
                    break

        # If some intermediate definitions were removed
        if removed_int:
            print('')
            print(str(len(removed_int)) + ' definitions removed:')
            for tens, term in zip(removed_int, removed_terms):
                print(tens + ": " + str(term[0]))

        options.print_divider()
        print('')

        # Returned shortened intermediate list
        terms = [t for t in terms if t not in removed_terms]

        return terms, removed_int

def remove_trans_rdm_const(terms, trans_int_list = None):

    options.print_header("WARNING")
    print('Terms w/o transRDM tensor in the expression will be removed. Set "remove_trans_rdm_constant"')
    print('flag to FALSE to preserve terms')

    # Create lists to split up SQA terms
    const_terms     = []
    trans_rdm_terms = []

    # Remove terms without tRDM tensors in r.h.s.
    for term_ind, term in enumerate(terms):
        creDes = False

        for tensor in term.tensors:
#            if isinstance(tensor, creOp) or isinstance(tensor, desOp) or isinstance(tensor, creDesTensor):
            if isinstance(tensor, creDesTensor) and tensor.trans_rdm:
                creDes = True
                break

            elif trans_int_list:
                tens_list = [tns.name for tns in term.tensors]
                for trans_int in trans_int_list:
                    if trans_int in tens_list:
                        creDes = True
                        break

        # Append to either list based on creDes flag
        if not creDes:
            const_terms.append(terms[term_ind])

        else:
            trans_rdm_terms.append(terms[term_ind])

    print('')
    print(str(len(const_terms)) + ' terms removed:')
    for term in const_terms:
        print(term)

    options.print_divider()
    print('Remaining terms: ' + str(len(trans_rdm_terms)))
    print('')

    return trans_rdm_terms, const_terms

def get_trans_intermediates(intermediate_list):

    # Store which intermediates are contracted over transition index
    trans_int_list = []

    # Iterate through list of intermediates
    for int_term, int_tensor in intermediate_list:

        # Make list of tensors that define intermediates
        ten_list = [t.name for t in int_term.tensors]

        # Check if one of the tensors is a tRDM
        if 'trdm' in ten_list:
             trans_int_list.append(int_tensor.name)

        # If an intermediate is defined in terms of another intermediate, make sure that intermediate
        # isn't defined w/ a tRDM
        elif trans_int_list:
             for trans_int in trans_int_list:
                 if trans_int in ten_list:
                     trans_int_list.append(int_tensor.name)

    return trans_int_list

def make_custom_name(sqa_tensor, rename_tuple):

    old_name = [old for old, new in rename_tuple]

    if sqa_tensor.name[:3] == 'INT':
        rename_index = old_name.index('INT')
        new_name = rename_tuple[rename_index][1] + sqa_tensor.name[3:]

    else:
        rename_index = old_name.index(sqa_tensor.name)
        new_name = rename_tuple[rename_index][1]

    return new_name

def append_CVS_slice(tens, tens_name, tens_indices, suffix):

    # Make a list out of the user-provided external indices
    cvs_indices_list = options.genEinsum.cvs_indices_list
    if isinstance(cvs_indices_list, str):
        cvs_indices_list = list(cvs_indices_list)

    val_indices_list = options.genEinsum.valence_indices_list
    if isinstance(val_indices_list, str):
        val_indices_list = list(val_indices_list)

    tens_indices = list(tens_indices)

    # Check whether tensor name needs an additional slice
    num_cvs = len([ind for ind in tens_indices if ind in cvs_indices_list])

    # Define ncvs string
    ncvs_string = 'ncvs'
    if suffix is not None:
        ncvs_string += '_' + suffix

    if num_cvs > 0:

        # Special condition for Kronecker delta
        if isinstance(tens, kroneckerDelta):

            tens_name = 'np.identity(' + ncvs_string + ')'

        # Add slices
        else:

            # Make 'starting' string to append to appropriate tensors
            to_append = '['

            # Iterate through all indices of tensor
            for ind in tens_indices:

                # Append slice through CVS indices
                if ind in cvs_indices_list:
                    to_append += ':' + ncvs_string + ','

                elif ind in val_indices_list:
                    to_append +=  ncvs_string + ':,'

                # Ignore non-CVS indices
                else:
                    to_append += ':,'

            # Remove extra comma and append end bracket
            to_append = to_append[:-1] + ']'

            # Append slices to tensor name
            tens_name += to_append

    return tens_name

def append_spin_integrated_slice(tens, tens_name, tens_indices):

    # List of spin index types
    spin_ind_types = [get_spin_index_type(ind) for ind in tens.indices]

    to_append = '['

    # Iterate through all indices of tensor
    for spin_ind_type in spin_ind_types:
        if spin_ind_type == options.alpha_type:
            to_append += '::2,'
        elif spin_ind_type == options.beta_type:
            to_append += '1::2,'

    # Remove extra comma and append end bracket
    to_append = to_append[:-1] + ']'

    # Append slices to tensor name
    tens_name += to_append

    #to_append = '_'
    #for spin_ind_type in spin_ind_types:
    #    if spin_ind_type == options.alpha_type:
    #        to_append += 'a'
    #    elif spin_ind_type == options.beta_type:
    #        to_append += 'b'
    #tens_name += to_append

    return tens_name

def sqalatex(terms, lhs = None, output = None, indbra = False, indket = None, print_default = True):

    texfile = output if output else 'output_default'

    header = f"""  
----------------------- SQA LATEX ----------------------------
    _____ ____    ___   __
   / ___// __ \  /   | / /____  _  __
   \__ \/ / / / / /| |/ __/ _ \| |/_/  Translate to Latex format and generate pdf
  ___/ / /_/ / / ___ / /_/  __/>  <    author:  Koushik Chatterjee
 /____/\___\_\/_/  |_\__/\___/_/|_|    date:  April 28, 2019
                                       VERSION : 1
 Copyright (C) 2018-2020  Koushik Chatterjee (koushikchatterjee7@gmail.com)

 Tex file : {texfile}.tex
 PDF file : {texfile}.pdf
--------------------------------------------------------------
    """
    print(header)

    modifier_tensor = {
        'bold': lambda s: rf'\boldsymbol{{{s}}}',
        'hat': lambda s: rf'\hat{{{s}}}',
        'bra': lambda s: rf'\langle\Psi_{{{s}}}\lvert',
        'ket': lambda s: rf'\rvert\Psi_{{{s}}}\rangle',
        #'gamma': lambda s: r'\Gamma',
        'kdelta': lambda s: r'\delta',
        'cre': lambda s: rf'\hat{{{s}}}^{{\dagger}}',
        'des': lambda s: rf'\hat{{{s}}}',
        'rdm': lambda s: r'\gamma',
    }

    #t_modifier = lambda s: r'\boldsymbol{'+s+r'}'
    t_modifier = lambda s: s

    lhs = 'M={}' if not lhs else f'{t_modifier(lhs)}={{}}'

    tex = []

    for term in terms:

        # Sign of coefficient
        if term.numConstant == 1.0:
            constant = " + "
        elif term.numConstant == -1.0:
            constant = " - "
        elif term.numConstant > 0:
            constant = f" +{term.numConstant} "
        else:
            constant = f" {term.numConstant} "

        credes_ops = ''
        tensor_names = ''
        gamma = ''

        for tens in term.tensors:
            tensor_name = tens.name
            names = [idx.name for idx in tens.indices]

            # Extract superscript and subscript indices based on number of indices
            n_indices = len(tens.indices)
            if n_indices == 1:
                superscripts = ''
                subscripts = names[0]
            elif n_indices == 2:
                superscripts = names[0]
                subscripts = names[1]
            elif n_indices == 4:
                superscripts = ''.join(names[:2])
                subscripts = ''.join(names[2:4])
            elif n_indices == 6:
                superscripts = ''.join(names[:3])
                subscripts = ''.join(names[3:6])
            elif n_indices == 8:
                superscripts = ''.join(names[:4])
                subscripts = ''.join(names[4:8])
            else:
                raise Exception(f"Not implemented: {n_indices}-index {tensor_name}...")
 
            if isinstance(tens, (creOp, desOp)):
                credes_ops += modifier_tensor[tensor_name](subscripts)

            elif tensor_name == 'gamma':
                bra = modifier_tensor['bra']('0')
                ket = modifier_tensor['ket']('0')
                creation_op = modifier_tensor['cre'](superscripts)
                destruction_op = modifier_tensor['des'](subscripts)
                gamma += bra + creation_op + destruction_op + ket + r"\:"

            else:
                if tensor_name in modifier_tensor:
                    tensor_names += modifier_tensor[tensor_name](tensor_name)
                else:
                    tensor_names += t_modifier(tensor_name)

                tensor_names += f"^{{{' '.join(superscripts)}}}"
                tensor_names += f"_{{{' '.join(subscripts)}}}"
                tensor_names += r"\:"

        # Append gamma expression if present
        if gamma:
            tensor_names += gamma

        # Wrap creation/destruction operators in bra-ket notation if present
        if credes_ops:
            ind = '0'
            bra_index = indbra if indbra else ind
            ket_index = indket if indket else ind

            bra = modifier_tensor['bra'](bra_index)
            ket = modifier_tensor['ket'](ket_index)
            tensor_names += bra + credes_ops + ket

        # Combine constant and tensor expression
        tex.append(constant + r'\:' + tensor_names)

    # Print to console if requested
    if print_default:
        print(r'\documentclass{article}')
        print(r'\usepackage{amsmath}')
        print(r'\begin{document}')
        print('')
        print('')
        print(r"\begin{align*}")
        print(lhs)
        for term_latex in tex:
            print(f" & {term_latex}" + r'\\')
        print(r"\end{align*}")
        print('')
        print('')
        print(r'\end{document}')

    # Write to file
    with open(f'{texfile}.tex', "w") as output_file:
        output_file.write(r'\documentclass{article}')
        output_file.write("\n")
        output_file.write(r'\usepackage{amsmath}')
        output_file.write("\n")
        output_file.write(r'\begin{document}')
        output_file.write("\n")
        output_file.write("\n")
        output_file.write(r"\begin{align*}")
        output_file.write("\n")
        output_file.write(lhs)
        output_file.write("\n")

        for term_latex in tex:
            output_file.write(f" & {term_latex}" + r'\\')
            output_file.write("\n")

        output_file.write(r"\end{align*}")
        output_file.write("\n")
        output_file.write("\n")
        output_file.write(r'\end{document}')

    # Compile PDF
    import subprocess
    try:
        result = subprocess.run(
            ['pdflatex', '-interaction=nonstopmode', f'{texfile}.tex'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True, check=False)
        if result.returncode != 0:
            print(f'LaTeX compilation failed with return code {result.returncode}')
    except Exception as e:
        print(f'LaTeX compilation error: {e}')

    return

def einsum_help():
    print("""\n        HELP :: 
        -----------
        terms           : A list of terms
        lhs_string      : Left hand side string (e.g. string 'M' = einsum ..)
        indices_string  : Einsum right side indix string (e.g. -> string 'p')
        transRDM        : Transition RDM True of False
        trans_indices_string   : Transition RDM string if True
        rhs_str         : Extra string for other kind of operation 
                          (e.g. transpose, copy, reshape .. etc)
        optimize        : By default optimization is true
        suffix          : Additional string atachement to the tensor name 
                          ( By default suffix = 'so' for spin orbitlas)
        rdm_str         : RDM string name
        tensor_rename   : Rename tensor if required 
                        (e.g. rename tensor 'X' to 'TEMP': X = 'TEMP'. 
                        For multiple tensors: X = 'TEMP', h = 'Hamiltonian', ..)
--------------------------------------------------------------""")
    yes = {'Yes','yes','y', 'Y', ''}
    no = {'NO','no','n','N'}
    sys.stdout.write("Do you want to continue [y/n] : ")
    choice = raw_input().lower()
    if choice in no:
       exit()
    print("-------------------------------------------------------------- ")
