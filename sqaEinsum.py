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

from sqaIndex import get_spatial_index_type, get_spin_index_type, \
                     is_core_index_type, is_active_index_type, is_virtual_index_type, \
                     is_alpha_index_type, is_beta_index_type

from sqaTensor import tensor, creOp, desOp, kroneckerDelta, creDesTensor

def genEinsum(terms, lhs_str = None, ind_str = None, trans_rdm = False, trans_ind_str = None, suffix = None,
              cvs_ind = None, val_ind = None, use_cvs_tensors = False, rm_trans_rdm_const = False, rm_core_int = False,
              intermediate_list = None, opt_einsum_terms = True, optimize = True, spin_integrated = False, use_spin_integrated_tensors = False,
              help = False, **tensor_rename):

    # Store custom names if provided by user
    custom_names = []
    if tensor_rename:
        for old_name, new_name in tensor_rename.items():
            custom_names.append((old_name, new_name))

    # Default to 'temp' as name of matrix being created w/ einsum function
    if not lhs_str:
        lhs_str = 'temp'

    # Default to spin-orbital suffix if not defined by user
    # if not suffix:
        # suffix = 'so'

    ################################################
    # IF PROVIDED, PRINT EINSUMS FOR INT TERMS
    ################################################
    if intermediate_list:

        # Create empty to list to store einsum expressions for provided intermediate terms
        int_einsum_list = []

        # Determined which intermediates are defined w/ trans_rdm contraction
        trans_int = []
        if trans_rdm:
            trans_int = get_trans_intermediates(intermediate_list)

        # If using effective Hamiltonian, remove double-counted contributions to core terms
        if rm_core_int:
            intermediate_list, removed_int = remove_core_int(intermediate_list, int_terms = True)

        # Iterate through the list of provided intermediates
        for int_ind, (int_term, int_tensor) in enumerate(intermediate_list):

            # Pass tensors of term to function to create string representation of contraction indices and tensor names
            int_tensor_inds, int_tensor_names = get_tensor_info(intermediate_list[int_ind][0].tensors, trans_rdm, trans_ind_str,
                                                                ''.join([i.name for i in intermediate_list[int_ind][1].indices]),
                                                                suffix, cvs_ind, val_ind, use_cvs_tensors,
                                                                spin_integrated, use_spin_integrated_tensors,
                                                                trans_int, custom_names)

            # Rename intermediates in tensor definitions
            if custom_names:
                if 'INT' in [x for x,y in custom_names]:
                    new_name = make_custom_name(int_tensor, custom_names)
                    int_einsum = new_name + ' = '

            else:
                int_einsum = int_tensor.name + ' = '

            # Define term for either optEinsum or built-in Numpy 'einsum' function
            if opt_einsum_terms:
                int_einsum += 'einsum('
            else:
                int_einsum += 'np.einsum('

            # Add contraction and tensor info
            int_tensor_info = (', '.join([str("'") + int_tensor_inds + str("'")] + int_tensor_names))
            int_einsum += int_tensor_info

            # Add optimize flag to einsum if enabled
            if optimize and opt_einsum_terms:
                int_einsum += ', optimize = einsum_type)'
            elif optimize and not opt_einsum_terms:
                int_einsum += ', optimize = True)'
            else:
                int_einsum += ')'

            # Append einsum definition to list, returned at the end of function
            int_einsum_list.append(int_einsum)

    ################################################
    # GENERATE EINSUM EXPRESSIONS FOR PROVIDED TERMS
    ################################################
    ## CONVERT ALL CRE/DES OBJECTS to CREDESTENSOR OBJECT
    for term_ind, term in enumerate(terms):

        # List for storing cre/des operators
        credes_ops = []

        # Append all cre/des operators to list
        for tens in term.tensors:
            if isinstance(tens, creOp) or isinstance(tens, desOp):
                credes_ops.append(tens)

        # Modify term in list to use creDesTensor object instead of cre/des objects
        if credes_ops:
            terms[term_ind].tensors = [ten for ten in terms[term_ind].tensors if ten not in credes_ops]
            terms[term_ind].tensors.append(creDesTensor(credes_ops, trans_rdm))

    # Constants terms in CAS blocks are removed by default, print warning
    if trans_rdm and rm_trans_rdm_const and intermediate_list and trans_int:
        terms, const_terms = remove_trans_rdm_const(terms, trans_int)

    elif trans_rdm and rm_trans_rdm_const:
        terms, const_terms = remove_trans_rdm_const(terms)

    # If using effective Hamiltonian, remove double-counted contributions to core terms
    if intermediate_list and rm_core_int and removed_int:
        terms, core_terms = remove_core_int(terms, removed_int)

    elif rm_core_int:
        terms, core_terms = remove_core_int(terms)

    # Create empty list for storing einsums
    einsum_list = []

    # Iterate through terms and create einsum expressions
    for term_ind, term in enumerate(terms):

        ## PROCEED WITH GENERATING EINSUMS
        # Start to define einsum string
        einsum = lhs_str + ' '

        # Determine sign of term
        pos = True
        if term.numConstant < 0:
            pos = False

        # Set up equals sign for first term and rest of terms
        if term_ind == 0 and pos:
            einsum += ' = '
        elif term_ind == 0 and not pos:
            einsum += '=- '
        elif term_ind != 0 and pos:
            einsum += '+= '
        elif term_ind != 0 and not pos:
            einsum += '-= '

        # Add appropriate scaling factor
        if round(abs(term.numConstant), 15) != 1.0:
            from fractions import Fraction

            frac_constant = Fraction(abs(term.numConstant)).limit_denominator()
            if round(float(frac_constant), 12) == round(abs(term.numConstant), 12):
                einsum = einsum + str(abs(frac_constant)) + ' * '
            else:
                einsum = einsum + str(abs(term.numConstant)) + ' * '

        # Define term for either optEinsum or built-in Numpy 'einsum' function
        if opt_einsum_terms:
            einsum += 'einsum('
        else:
            einsum += 'np.einsum('

        # Pass tensors of term to function to create string representation of contraction indices and tensor names
        if trans_rdm and intermediate_list:
            tensor_inds, tensor_names = get_tensor_info(term.tensors, trans_rdm, trans_ind_str, ind_str, suffix,
                                                        cvs_ind, val_ind, use_cvs_tensors, spin_integrated, use_spin_integrated_tensors,
                                                        trans_int, custom_names)
        else:
            tensor_inds, tensor_names = get_tensor_info(term.tensors, trans_rdm, trans_ind_str, ind_str, suffix,
                                                        cvs_ind, val_ind, use_cvs_tensors, spin_integrated, use_spin_integrated_tensors,
                                                        custom_names)

        # Add contraction and tensor info
        tensor_info = (', '.join([str("'") + tensor_inds + str("'")] + tensor_names))
        einsum += tensor_info

        # Add optimize flag to einsum if enabled
        if optimize and opt_einsum_terms:
            einsum += ', optimize = einsum_type)'
        elif optimize and not opt_einsum_terms:
            einsum += ', optimize = True)'
        else:
            einsum += ')'

        # Append a '.copy()' function call if the term is made up of only one tensor
        if len(term.tensors) == 1:
            einsum += '.copy()'

        # Append completed einsum to list
        if not 'none' in einsum:
            einsum_list.append(einsum)

    # Modify return for intermediate term definition
    if intermediate_list:
        return int_einsum_list, einsum_list
    else:
        return einsum_list

def get_tensor_info(sqa_tensors, trans_rdm, trans_ind_str, ind_str, suffix, cvs_ind = False, val_ind = False, use_cvs_tensors = False,
                    spin_integrated = False, use_spin_integrated_tensors = False, trans_int = None, custom_names = None):

    # Pre-define list of names of tensors used in SQA and make list to store any new tensor 'types'
    tensor_names = []
    tensor_inds  = []

    # Iterate through all the provided tensors
    for tens in sqa_tensors:
        # Handle special case of kroneckerDelta (kdelta) object
        if isinstance(tens, kroneckerDelta):
            tensor_name = 'np.identity('

            # Determine orbital space of kdelta
            if (is_core_index_type(tens.indices[0]) and is_core_index_type(tens.indices[1])):
                if use_cvs_tensors:
                    if cvs_ind and val_ind:
                        if (tens.indices[0].name in cvs_ind) and (tens.indices[1].name in cvs_ind):
                            orb_space = 'ncvs'
                        elif (tens.indices[0].name in val_ind) and (tens.indices[1].name in val_ind):
                            orb_space = 'nval'
                        elif (((tens.indices[0].name not in cvs_ind) and (tens.indices[0].name not in val_ind)) and
                            ((tens.indices[1].name not in cvs_ind) and (tens.indices[1].name not in val_ind))):
                            orb_space = 'ncore'
                        else:
                            orb_space = 'none'
                    elif cvs_ind:
                        if (tens.indices[0].name in cvs_ind) and (tens.indices[1].name in cvs_ind):
                            orb_space = 'ncvs'
                        elif (((tens.indices[0].name not in cvs_ind)) and ((tens.indices[1].name not in cvs_ind))):
                            orb_space = 'ncore'
                        else:
                            orb_space = 'none'
                    elif val_ind:
                        if (tens.indices[0].name in val_ind) and (tens.indices[1].name in val_ind):
                            orb_space = 'nval'
                        elif (((tens.indices[0].name not in val_ind)) and ((tens.indices[1].name not in val_ind))):
                            orb_space = 'ncore'
                        else:
                            orb_space = 'none'
                    else:
                        orb_space = 'ncore'
                else:
                    orb_space = 'ncore'

            elif (is_active_index_type(tens.indices[0]) and is_active_index_type(tens.indices[1])):
                orb_space = 'ncas'

            elif (is_virtual_index_type(tens.indices[0]) and is_virtual_index_type(tens.indices[1])):
                orb_space = 'nextern'

            else:
                raise TypeError('WARNING: The indices of the kronecker delta term do not belong to the same orbital sub-space')

            if suffix:
                orb_space += '_' + suffix
            tensor_name += orb_space + ')'

            # Rename if custom name is provided
            if custom_names:
                if ('kdelta') in [x for x,y in custom_names]:
                    new_name = make_custom_name(tens, custom_names)
                    tensor_name = new_name + '_'

        # Handle special case of orbital energy vector
        elif len(tens.indices) == 1 and (tens.name == 'e' or tens.name == 'E'):

            tensor_name = str(tens.name).lower() + '_'

            # Determine orbital space of energies
            if is_core_index_type(tens.indices[0]):
                if use_cvs_tensors:
                    if cvs_ind and val_ind:
                        if tens.indices[0].name in cvs_ind:
                            orb_space = 'cvs'
                        elif tens.indices[0].name in val_ind:
                            orb_space = 'val'
                        else:
                            orb_space = 'core'
                    elif cvs_ind:
                        if tens.indices[0].name in cvs_ind:
                            orb_space = 'cvs'
                        else:
                            orb_space = 'core'
                    elif val_ind:
                        if tens.indices[0].name in val_ind:
                            orb_space = 'val'
                        else:
                            orb_space = 'core'
                    else:
                        orb_space = 'ncore'
                else:
                    orb_space = 'core'

            elif is_virtual_index_type(tens.indices[0]):
                orb_space = 'extern'

            if suffix:
                orb_space += '_' + suffix
            tensor_name += orb_space

            # Rename if custom name is provided
            if custom_names:
                if ('e' or 'E') in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Handle special case of RDM tensor
        elif isinstance(tens, creDesTensor):
            tensor_name = tens.name + '_'

            # Modify name of RDM to reflect particle number
            for op in tens.ops:
                if isinstance(op, creOp):
                    tensor_name += 'c'
                elif isinstance(op, desOp):
                    tensor_name += 'a'

            # Append suffix
            if suffix:
                tensor_name += '_' + suffix

            # Rename if custom name is provided
            if custom_names:
                if ('rdm') in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Name remaining tensors w/ same convention of orbital space and suffix
        elif tens.name == 'h' or tens.name == 'v' or tens.name == 't1' or tens.name == 't2':
            tensor_name = tens.name + '_'

            # Append letter representing orbital subspace of indices
            for i in range(len(tens.indices)):
                if is_active_index_type(tens.indices[i]):
                    tensor_name += 'a'
                elif is_core_index_type(tens.indices[i]):
                    if use_cvs_tensors:
                        if cvs_ind and val_ind:
                            if (tens.indices[i].name in cvs_ind):
                                tensor_name += 'x'
                            elif (tens.indices[i].name in val_ind):
                                tensor_name += 'v'
                            else:
                                tensor_name += 'c'
                        elif cvs_ind:
                            if (tens.indices[i].name in cvs_ind):
                                tensor_name += 'x'
                            else:
                                tensor_name += 'c'
                        elif val_ind:
                            if (tens.indices[i].name in val_ind):
                                tensor_name += 'v'
                            else:
                                tensor_name += 'c'
                    else:
                        tensor_name += 'c'
                else:
                    tensor_name += 'e'

            # Append spin-integrated suffix if required
            if use_spin_integrated_tensors:
                spin_suffix = '_'
                for i in range(len(tens.indices)):
                    if is_alpha_index_type(tens.indices[i]):
                        spin_suffix += 'a'
                    elif is_beta_index_type(tens.indices[i]):
                        spin_suffix += 'b'
                tensor_name += spin_suffix

            # Append suffix
            if not (tens.name == 't1' or tens.name == 't2') and suffix:
                tensor_name += '_' + suffix

            # Rename if custom name is provided
            if custom_names:
                if tens.name in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Account for intermediate tensors and any custom tensors
        else:
            # Make copy of tensor name
            tensor_name = '%s' % tens.name

            # Append spin-integrated suffix if required
            if use_spin_integrated_tensors:
                spin_suffix = '_'
                for i in range(len(tens.indices)):
                    if is_alpha_index_type(tens.indices[i]):
                        spin_suffix += 'a'
                    elif is_beta_index_type(tens.indices[i]):
                        spin_suffix += 'b'
                tensor_name += spin_suffix

            # Allow to rename intermediate tensors in term definitions
            if custom_names:
                if tens.name[:3] == 'INT' and 'INT' in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Create indices of tensor as string
        indices = ''.join([i.name for i in tens.indices])

        # Append 'slices' to appropiate dimensions of spin-integrated tensors
        if spin_integrated and (not use_spin_integrated_tensors):
            tensor_name = append_spin_integrated_slice(tens, tensor_name, indices)

        # Append 'slices' to appropiate dimensions of tensors w/ CVS core indices
        if (cvs_ind is not None) and (not use_cvs_tensors):
            tensor_name = append_CVS_slice(tens, tensor_name, indices, cvs_ind, val_ind, suffix)

        # Append name of tensor (after and modifications due to special cases)
        tensor_names.append(tensor_name)

        # Append transition state index to appropriate set of indices
        if isinstance(tens, creDesTensor) and tens.trans_rdm:
            indices = trans_ind_str + indices

        elif trans_int and (tens.name in trans_int):
            indices = trans_ind_str + indices

        # Append completed index string to list
        tensor_inds.append(indices)

    # Convert list of indices into one comma-separated string and prepare to append external index string
    tensor_inds = ','.join(tensor_inds)

    # Check if rhs string or transition index is provided before adding arrow
    if trans_ind_str or ind_str:
        tensor_inds += '->'

    # Append transition index first, if present
    if trans_ind_str:
        tensor_inds += trans_ind_str

    # Append rhs string, if provided
    if ind_str:
        tensor_inds += ind_str

    return tensor_inds, tensor_names

def remove_core_int(terms, removed_int = None, int_terms = False):

    # Remove terms from standard term list
    if not int_terms:
        print ('\n--------------------------------- WARNING ---------------------------------')
        print ('Terms with a contraction over repeating dummy core indices of 2e- integrals')
        print ('will be removed. Set "rm_core_int" flag to FALSE to preserve terms')

        # Create lists to split up SQA terms
        kept_terms = []
        core_terms = []

        # Separate out the terms that have redundant 2e- integral contractions over core space
        for term_ind, term in enumerate(terms):
            coreTerm = False
            for tens_ind, tens in enumerate(term.tensors):
                if tens.name == 'v':
                    if (
                        (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[2].name)
                        or (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[3].name)
                        or (terms[term_ind].tensors[tens_ind].indices[1].name) == (terms[term_ind].tensors[tens_ind].indices[2].name)
                        or (terms[term_ind].tensors[tens_ind].indices[1].name) == (terms[term_ind].tensors[tens_ind].indices[3].name)
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

        print ('')
        print (str(len(core_terms)) + ' terms removed:')
        for term in core_terms:
            print (term)

        print ('---------------------------------------------------------------------------')
        print ('Remaining terms: ' + str(len(kept_terms)))
        print ('')

        return kept_terms, core_terms

    # Filter through intermediate definitions
    else:
        print ('--------------------------------- WARNING ---------------------------------')
        print ('Intermediate tensors defined w/ contractions over repeating dummy core indices of')
        print ('2e- integrals will be removed. Set "rm_core_int" flag to FALSE to preserve definitions')

        # Track which tensor definitions are removed and kept
        removed_int   = []
        removed_terms = []

        # Determine which intermediate definitions to remove
        for int_ind, (int_term, int_tensor) in enumerate(terms):
            for tens in int_term.tensors:

                # If intermediate is defined w/ 2e- integral
                if tens.name == 'v':
                    if (
                        (tens.indices[0].name) == (tens.indices[2].name)
                        or (tens.indices[0].name) == (tens.indices[3].name)
                        or (tens.indices[1].name) == (tens.indices[2].name)
                        or (tens.indices[1].name) == (tens.indices[3].name)
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
            print ('')
            print (str(len(removed_int)) + ' definitions removed:')
            for tens, term in zip(removed_int, removed_terms):
                print (tens + ": " + str(term[0]))

        print ('---------------------------------------------------------------------------')
        print ('')

        # Returned shortened intermediate list
        terms = [t for t in terms if t not in removed_terms]

        return terms, removed_int

def remove_trans_rdm_const(terms, trans_int_list = None):

    print ('--------------------------------- WARNING ---------------------------------')
    print ('Terms w/o transRDM tensor in the expression will be removed. Set "rm_trans_rdm_const"')
    print ('flag to FALSE to preserve terms')

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

    print ('')
    print (str(len(const_terms)) + ' terms removed:')
    for term in const_terms:
        print (term)

    print ('---------------------------------------------------------------------------')
    print ('Remaining terms: ' + str(len(trans_rdm_terms)))
    print ('')

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

def append_CVS_slice(tens, tens_name, tens_indices, cvs_ind, val_ind, suffix):

    # Make a list out of the user-provided external indices
    if cvs_ind:
        cvs_ind = list(cvs_ind)
    if val_ind:
        val_ind = list(cvs_ind)

    tens_indices = list(tens_indices)

    # Check whether tensor name needs an additional slice
    num_cvs = len([ind for ind in tens_indices if ind in cvs_ind])

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
                if ind in cvs_ind:
                    to_append += ':' + ncvs_string + ','

                elif ind in val_ind:
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
    spin_ind_types = [get_spin_indType(ind.indType) for ind in tens.indices]

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

    return tens_name

def sqalatex(terms, lhs = None, output = None, indbra = False, indket = None, print_default = True):

 if not output:
  # texfile = r'latex_output.tex'
   texfile = r'output_default'
 else:
   texfile = output

 print("""\n----------------------- SQA LATEX ----------------------------
    _____ ____    ___   __
   / ___// __ \  /   | / /____  _  __
   \__ \/ / / / / /| |/ __/ _ \| |/_/  Translate to Latex format and generate pdf
  ___/ / /_/ / / ___ / /_/  __/>  <    author:  Koushik Chatterjee
 /____/\___\_\/_/  |_\__/\___/_/|_|    date:  April 28, 2019
                                       VERSION : 1
 Copyright (C) 2018-2020  Koushik Chatterjee (koushikchatterjee7@gmail.com)
 
 Tex file : %s
 PDF file : %s
--------------------------------------------------------------""" % (texfile+r'.tex', texfile+r'.pdf'))

 modifier_tensor = {
     'bold': lambda s: r'\boldsymbol{'+s+r'}',
     'hat': lambda s: r'\hat{'+s+r'}',
     'bra': lambda s: r'\langle\Psi_{'+s+r'}\lvert',
     'ket': lambda s: r'\rvert\Psi_{'+s+r'}\rangle',
#     'gamma': lambda s: r'\Gamma',
     'kdelta': lambda s: r'\delta',
     'cre': lambda s: r'\hat{'+s+r'}^{\dagger}',
     'des': lambda s: r'\hat{'+s+r'}',
 }
# t_modifier = lambda s: r'\boldsymbol{'+s+r'}'
 t_modifier = lambda s: s

 if not lhs:
  lhs = 'M={}'
 else:
  lhs = t_modifier(lhs)+r'={}'

 tex = []

 for term in terms:

     constant = ''
     if (term.numConstant == 1.0):
        constant = " + "
     elif (term.numConstant == -1.0):
        constant = " - "
     else:
        constant = " %s " % str(term.numConstant)
        if (term.numConstant > 0):
            
#          constant += " %s " % str(Fraction(Decimal('term.numConstant')))
          constant = " +%s " % str(term.numConstant)
#
     cre_count = 0
     des_count = 0
     credes = ''
     name = ''
     gamma = ''
     for i in range(len(term.tensors)):

         tens = term.tensors[i]
         s = tens.name
     #    credes = None
     #    name = ''

         supers = ''
         subs   = ''
         index = len(tens.indices)
         if (index == 1):
            subs   = tens.indices[0].name
         elif(index == 2):
            supers = tens.indices[0].name
            subs   = tens.indices[1].name
         elif (index == 4):
            supers = tens.indices[0].name+tens.indices[1].name
            subs   = tens.indices[2].name+tens.indices[3].name

         else:
             raise Exception("Not implemented ...")

         if not (isinstance(tens, creOp) or isinstance(tens, desOp)):
            if (s == 'gamma'):
               bra = modifier_tensor['bra']('0')
               ket = modifier_tensor['ket']('0')
               ind1 = modifier_tensor['cre'](supers)
               ind2 = modifier_tensor['des'](subs)
               gamma += bra+ind1+ind2+ket+"\:"
            else:
               if s in modifier_tensor:
                  name += modifier_tensor[s](s)
               else:
                  name += t_modifier(s)

               name += "^{%s}" % " ".join(supers)
               name += "_{%s}" % " ".join(subs)
               name +="\:"
         else:
            if (isinstance(tens, creOp)):
               cre_count += 1
            if (isinstance(tens, desOp)):
               des_count += 1
            credes += modifier_tensor[s](subs)

     if(len(gamma) > 0):
        name += gamma
     if (len(credes) > 0):
         ind = r'0'
         if not indbra:
            indbra = ind
         if not indket:
            indket = ind
         bra = modifier_tensor['bra'](indbra)
         ket = modifier_tensor['ket'](indket)
         name += bra+credes+ket

     tex.append(constant+r'\:'+name)


 if print_default:
    print r'\documentclass{article}'
    print r'\usepackage{amsmath}'
    print r'\begin{document}'
    print ''
    print ''
#    print r"\begin{equation}"
    print r"\begin{align*}"
    print lhs
    for i in tex:
#      print " & "+i+' \\\\'
      print " & "+i+r'\\'
    print r"\end{align*}"
#    print r"\end{equation}"
    print ''
    print ''
    print r'\end{document}'


 ### write to a file ###
# if not output:
#  # texfile = r'latex_output.tex'
#   texfile = r'latex_output'
# else:
#   texfile = output
 output = open(texfile+r'.tex', "w")
 output.write(r'\documentclass{article}')
 output.write("\n")
 output.write(r'\usepackage{amsmath}')
 output.write("\n")
 output.write(r'\begin{document}')
 output.write("\n")
 output.write('')
 output.write("\n")
 output.write('')
 output.write("\n")
# output.write(r"\begin{equation*}")
 output.write(r"\begin{align*}")
 output.write("\n")
 output.write(lhs)
 for i in tex:
#   output.write(" & "+i+' \\\\')
   output.write(" & "+i+r'\\')
   output.write("\n")
 output.write(r"\end{align*}")
 output.write("\n")
# print r"\end{equation*}"
 output.write('')
 output.write("\n")
 output.write('')
 output.write("\n")
 output.write(r'\end{document}')

# os.system("pdflatex latex_output.tex")
 procs = []
 try:
     pread, pwrite = os.pipe()
     cmd = ['pdflatex', texfile+r'.tex']
#     proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
     proc = subprocess.Popen(cmd, stdout=pwrite, stderr=subprocess.STDOUT)
     procs.append(proc)
     os.close(pwrite)
     os.close(pread)

 except OSError as e:
   #  sys.exit()
     print 'Latex compilation error ...'

# pdf()
# proc_cleanup(procs)
 return

def einsum_help():
    print("""\n        HELP :: 
        -----------
        terms         : A list of terms
        lhs_str       : Left hand side string (e.g. string 'M' = einsum ..)
        ind_str       : Einsum right side indix string (e.g. -> string 'p')
        transRDM      : Transition RDM True of False
        trans_ind_str : Transition RDM string if True
        rhs_str       : Extra string for other kind of operation 
                        (e.g. transpose, copy, reshape .. etc)
        optimize      : By default optimization is true
        suffix        : Additional string atachement to the tensor name 
                        ( By default suffix = 'so' for spin orbitlas)
        rdm_str       : RDM string name
        tensor_rename : Rename tensor if required 
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
