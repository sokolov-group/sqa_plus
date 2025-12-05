# Copyright 2009-2022 SecondQuantizationAlgebra Developers. All Rights Reserved.
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
# In addition, any modification or use of this software should
# cite the following paper:
#
#   E. Neuscamman, T. Yanai, and G. K.-L. Chan.
#   J. Chem. Phys. 130, 124102 (2009)
#
# Author: Alexander Yu. Sokolov <alexander.y.sokolov@gmail.com>
# .       Koushik Chatterjee <koushikchatterjee7@gmail.com>
#

import sys, time
from .sqaTensor import kroneckerDelta, sfExOp, creOp, desOp
from .sqaTerm import term, termChop, combineTerms
from .sqaMisc import makeTuples, allDifferent
from .sqaSymmetry import symmetry
from .sqaOptions import options

from .sqaNormalOrder import normalOrder

from .sqaIndex import get_spatial_index_type, get_spin_index_type, \
                     is_core_index_type, is_active_index_type, is_virtual_index_type, \
                     is_cvs_core_index_type, is_cvs_valence_index_type

def matrixBlock(terms):
    "Construct matrix block."

    startTime = time.time()
    options.print_header("SQA Automation")
    sys.stdout.flush()
    fTerms = []

    # Import options from sqaOptions class
    remove_trans_rdm_constant = options.matrixBlock.remove_trans_rdm_constant

    # Make sure that input expression is normal-ordered with respect to physical vacuum
    nterms = []
    for t in terms:
        t_no = normalOrder(t)
        nterms.extend(t_no)
    del(terms)

    # Filter zero terms wrt virtual (note: Filter first for virtual orbitals)
    filterVirtual(nterms)
    termChop(nterms)

    # Normal ordering with respect to core orbitals
    fTerms = normalOrderCore(nterms)
    del(nterms)

    # Evaluate Kronecker delta
    for t in fTerms:
        t.contractDeltaFuncs()

    # Filter zero terms wrt core (note: after filtering virtual, then do for core)
    filterCore(fTerms)
    termChop(fTerms)
    combineTerms(fTerms)

    # Contract delta function for both non-dummy indices
    contractDeltaFuncs_nondummy(fTerms)

    # If (remove_trans_rdm_constant = True) => Remove those constant terms
    if remove_trans_rdm_constant:
        for trm in fTerms:
            # If no creation/destruction operators, set the constant to zero
            if not any(isinstance(t, (creOp, desOp)) for t in trm.tensors):
                trm.numConstant = 0.0
        termChop(fTerms)

    # Reorder tensor indices: (core < active < virtual) order
    reorder_tensor_indices(fTerms)

    # Dummy indices label upate
    dummyLabel(fTerms)

    # Print the final results
    options.print_header("Final results")
    sys.stdout.flush()
    for t in fTerms:
        print(t)

    print("")
    print("Total terms : %s" % len(fTerms))
    print("SQA automation time :  %.3f seconds" % (time.time() - startTime))
    sys.stdout.flush()

    return fTerms

def dummyLabel(_terms, keep_user_defined_dummy_names = True):
    "A function to relabel dummy indices."

    # Import options from sqaOptions class
    user_defined_indices = options.user_defined_indices

    print("Dummy indices relabeling...")
    sys.stdout.flush()

    for _term_ind, _term in enumerate(_terms):
        mymap = {}

        coreInd = list('ijklmnopq')
        actvInd = list('xyzwuvstr')
        virtInd = list('abcdefgh')

        if keep_user_defined_dummy_names:
            coreInd = [c for c in coreInd if c not in user_defined_indices]
            actvInd = [a for a in actvInd if a not in user_defined_indices]
            virtInd = [v for v in virtInd if v not in user_defined_indices]

        if options.verbose:
            _term_unlabeled = _term.copy()

        for _tensor_ind, _tensor in enumerate(_term.tensors):
            for _index_ind, _index in enumerate(_tensor.indices):

                # Decide which new label to assign
                index_type = _index.indType
                index_name = _index.name
                index_summed = _index.isSummed
                index_user_defined = _index.userDefined

                if (index_summed and 
                    ((keep_user_defined_dummy_names and not index_user_defined)
                      or not keep_user_defined_dummy_names)):
                    if index_name not in mymap:
                        if is_core_index_type(index_type):
                            mymap[index_name] = coreInd.pop(0)
                        elif is_active_index_type(index_type):
                            mymap[index_name] = actvInd.pop(0)
                        elif is_virtual_index_type(index_type):
                            mymap[index_name] = virtInd.pop(0)

                    # Update the label
                    _terms[_term_ind].tensors[_tensor_ind].indices[_index_ind].name = mymap[index_name]

        if options.verbose:
            print("{:} ---> {:}".format(_term_unlabeled, _terms[_term_ind]))

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return

def filterVirtual(_terms):
    "A function to calculate expectation value wrt virtual: filter zero terms wrt virtual."

    print("Computing expectation value with respect to virtual ...")
    sys.stdout.flush()

    for t_term in _terms:
        # Collect all indices for cre/desOp in t_term
        cre_des_indices = (
            ind
            for t_tensor in t_term.tensors
            if isinstance(t_tensor, (creOp, desOp))
            for ind in t_tensor.indices
        )
        # Set t_term to 0 if any indices are virtual type
        if any(is_virtual_index_type(ind.indType) for ind in cre_des_indices):
            t_term.numConstant = 0.0

    if options.verbose:
        print("")
        for t_term in _terms:
            print(t_term)
        print("")

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return

def filterCore(_terms):
    "A function to calculate expectation value wrt core: filter zero terms wrt core."

    print("Computing expectation value with respect to core ...")
    sys.stdout.flush()

    for t_term in _terms:
        # Collect all indices for cre/desOp in t_term
        cre_des_indices = (
            ind
            for t_tensor in t_term.tensors
            if isinstance(t_tensor, (creOp, desOp))
            for ind in t_tensor.indices
        )
        # Set t_term to 0 if any indices are core type
        if any(is_core_index_type(ind.indType) for ind in cre_des_indices):
            t_term.numConstant = 0.0

    if options.verbose:
        print("")
        for _term in _terms:
            print(_term)
        print("")

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return

def normalOrderCore(_terms):
    print("Normal ordering with respect to core ...")
    sys.stdout.flush()
    ordered_terms = []

    for _term in _terms:
        ordered_term = normOrderCor(_term)
        ordered_terms.extend(ordered_term)

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return ordered_terms

def normOrderCor(_term):
    "Returns a list of terms resulting from normal ordering the operators in inTerm."
    if options.verbose:
        print('Term = %s' % _term)

    # check if is a term
    if not isinstance(_term, term):
        raise TypeError("Input term must be of class term")

    # determine what types of operators the term contains
    has_creDesOps = False
    has_sfExOps = False
    for _tensor in _term.tensors:
        if isinstance(_tensor, (creOp, desOp)):
            has_creDesOps = True
        elif isinstance(_tensor, sfExOp):
            has_sfExOps = True

    # If term has both creation/destruction operators and spin free excitation operators raise an error
    if has_creDesOps and has_sfExOps:
        raise RuntimeError("Normal ordering not implemented when both creOp/desOp and sfExOp tensors are present")

    # Normal ordering for creOp/desOp
    elif has_creDesOps:

        # Separate the cre/des operators from other tensors
        ops = []
        nonOps = []
        for _tensor in _term.tensors:
            if isinstance(_tensor, (creOp, desOp)):
                ops.append(_tensor.copy())
            else:
                nonOps.append(_tensor.copy())

        # Generate all contraction pairs
        contractionPairs = []
        for i in range(len(ops)):
            iTerm = ops[i]
            iType = iTerm.indices[0].indType

            for j in range(i+1,len(ops)):
                jTerm = ops[j]
                jType = jTerm.indices[0].indType

                if isinstance(iTerm, creOp) and isinstance(jTerm, desOp):
                    if is_core_index_type(iType) or is_core_index_type(jType):
                        contractionPairs.append((i,j));

        # Determine maximum contraction order
        creCount = 0
        maxConOrder = 0
        for i in range(len(ops)-1,-1,-1):
            iTerm = ops[i]
            if isinstance(iTerm, desOp):
                creCount +=1
            elif isinstance(iTerm, creOp) and creCount > 0:
                maxConOrder += 1
                creCount -= 1
        del(creCount,iTerm)

        # Generate all contractions
        contractions = []
        for i in range(maxConOrder+1):
            subCons = makeTuples(i,contractionPairs)
            j = 0
            while j < len(subCons):
                creOpTags = []
                desOpTags = []
                for k in range(i):
                    creOpTags.append(subCons[j][k][1])
                    desOpTags.append(subCons[j][k][0])
                if allDifferent(creOpTags) and allDifferent(desOpTags):
                    j += 1
                else:
                    del(subCons[j])
            for j in range(len(subCons)):
                contractions.append(subCons[j])
        del(subCons,creOpTags,desOpTags,contractionPairs)

        # For each contraction, generate the resulting term
        ordered_terms = []
        for contraction in contractions:
            conSign = 1
            deltaFuncs = []
            subOpString = []
            subOpString.extend(ops)
            for conPair in contraction:
                index1 = ops[conPair[0]].indices[0]
                index2 = ops[conPair[1]].indices[0]
                deltaFuncs.append(kroneckerDelta([index1,index2]))
                subOpString[conPair[0]] = 'contracted'
                subOpString[conPair[1]] = 'contracted'
                for q in subOpString[conPair[0]+1:conPair[1]]:
                    if not (q is 'contracted'):
                        conSign *= -1
            i = 0
            while i < len(subOpString):
                if subOpString[i] is 'contracted':
                    del(subOpString[i])
                else:
                    i += 1
            (sortSign, sortedOps) = sortOpsCore(subOpString)
            totalSign = conSign * sortSign

            ordered_tensors = []
            ordered_tensors.extend(nonOps)
            ordered_tensors.extend(deltaFuncs)
            ordered_tensors.extend(sortedOps)
            ordered_terms.append(term(totalSign * _term.numConstant, _term.constants, ordered_tensors))

    # Normal ordering for sfExOps
    elif has_sfExOps:
        # Make separate lists of the spin free excitation operators and other tensors
        raise Exception('This code does not support for now')

    else:
        ordered_terms = []
        ordered_terms = [_term]

    if options.verbose:
        print("Terms after normal ordering:")
        for ordered_term in ordered_terms:
            print(ordered_term)

    return ordered_terms

def sortOpsCore(_unsorted_ops, returnPermutation = False):
    """
    Sorts a list of creation/destruction operators into normal order and alphabetically, without performing contractions.
    Returns the overall sign resulting from the sort and the sorted operator list.
    """
    sorted_ops = list(_unsorted_ops)
    n_ops = len(sorted_ops)
    sign = 1

    perm = None
    if returnPermutation:
        perm = list(range(n_ops))

    i = 0
    while i < n_ops-1:
        # Bubble core creation operators to the right
        current_op = sorted_ops[i]
        if isinstance(current_op, creOp) and is_core_index_type(current_op.indices[0]):
            j = i
            while j+1 < n_ops and isinstance(sorted_ops[j+1], creOp) and is_core_index_type(sorted_ops[j+1].indices[0]):
                j += 1
            if j+1 >= n_ops:
                break

            sorted_ops[j], sorted_ops[j+1] = sorted_ops[j+1], sorted_ops[j]
            sign *= -1
            if perm is not None:
                perm[j], perm[j+1] = perm[j+1], perm[j]

            i = 0
            continue

        # Bubble core destruction operators to the left
        next_op = sorted_ops[i+1]
        if isinstance(next_op, desOp) and is_core_index_type(next_op.indices[0]):
            sorted_ops[i], sorted_ops[i+1] = next_op, current_op
            sign *= -1
            if perm is not None:
                perm[i], perm[i+1] = perm[i+1], perm[i]

            if sorted_ops[i+1].name == sorted_ops[i].name:
                break

            i = 0
            continue

        i += 1

    return (sign, sorted_ops) if not returnPermutation else (sign, sorted_ops, perm)

def contractDeltaFuncs_nondummy(_terms):
    "Contracts delta function for both non-dummy indices only wrt to orbitals subspaces, otherwise use 'contractDeltaFuncs' function."

    print("Contract delta function for non-dummy indices ...")
    sys.stdout.flush()

    for term in _terms:
        for t in term.tensors:
            if not isinstance(t, kroneckerDelta):
                continue

            i0, i1 = t.indices[0], t.indices[1]
            both_summed = i0.isSummed and i1.isSummed
            spatial_match = get_spatial_index_type(i0) == get_spatial_index_type(i1)
            spin_match = get_spin_index_type(i0) == get_spin_index_type(i1)
            if not both_summed and not (spatial_match and spin_match):
                term.numConstant = 0.0
                break

    termChop(_terms)

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return _terms

def reorder_tensor_indices(_terms):
    print("Reordering indices according to core < active < virtual...")
    sys.stdout.flush()

    # Import options from sqaOptions class
    legacy_ordering = options.legacy_ordering
    chemists_notation = options.chemists_notation

    if chemists_notation:
        v2e_sym_braket = symmetry((1,0,3,2), 1)
    else:
        v2e_sym_braket = symmetry((2,3,0,1), 1)

    for ind_unordered_term, unordered_term in enumerate(_terms):
        for ind_unordered_tensor, unordered_tensor in enumerate(unordered_term.tensors):
            reorder_tensor = ((unordered_tensor.name == 'h' and len(unordered_tensor.indices) == 2) or
                              (unordered_tensor.name == 'v' and len(unordered_tensor.indices) == 4) or
                              ((unordered_tensor.name[0] == 't') and
                               ((len(unordered_tensor.indices) == 2) or len(unordered_tensor.indices) == 4)))

            if (unordered_tensor.name == 'v') and (v2e_sym_braket not in unordered_tensor.symmetries):
                unordered_tensor_symmetries = unordered_tensor.symmetries
                unordered_tensor_symmetries.append(v2e_sym_braket)
                unordered_tensor.symmetries = unordered_tensor_symmetries

            if reorder_tensor:
                original_rank = []
                for ind in unordered_tensor.indices:
                    if is_core_index_type(ind) or is_cvs_core_index_type(ind):
                        original_rank.append(3)
                    elif is_cvs_valence_index_type(ind):
                        original_rank.append(2)
                    elif is_active_index_type(ind):
                        original_rank.append(1)
                    elif is_virtual_index_type(ind):
                        original_rank.append(0)

                permutes_indices, permutes_factors = unordered_tensor.symPermutes()
                permutes_rank = []

                for permute_indices in permutes_indices:
                    permute_rank = []

                    permuted_indices = [unordered_tensor.indices[ind] for ind in permute_indices]
                    for ind in permuted_indices:
                        if is_core_index_type(ind) or is_cvs_core_index_type(ind):
                            permute_rank.append(3)
                        elif is_cvs_valence_index_type(ind):
                            permute_rank.append(2)
                        elif is_active_index_type(ind):
                            permute_rank.append(1)
                        elif is_virtual_index_type(ind):
                            permute_rank.append(0)
                    permutes_rank.append(permute_rank)

                permutes_rank, permutes_indices, permutes_factors = zip(*sorted(zip(permutes_rank, permutes_indices, permutes_factors), reverse=True))
                permutes_rank_ind = [ind for ind, rank in enumerate(permutes_rank) if rank == permutes_rank[0]][-1]

                permuted_indices = [unordered_tensor.indices[ind] for ind in permutes_indices[permutes_rank_ind]]
                ordered_tensor = unordered_tensor.copy()
                ordered_tensor.indices = permuted_indices
                order_factor = permutes_factors[permutes_rank_ind]

                # Legacy ordering to obtain spin-orbital Prism integrals exceptions
                if (ordered_tensor.name in ['t1', 't2']) and legacy_ordering and not chemists_notation:

                    inds_ccae = [options.core_type,   options.core_type,   options.active_type, options.virtual_type]
                    inds_caae = [options.core_type,   options.active_type, options.active_type, options.virtual_type]
                    inds_aaae = [options.active_type, options.active_type, options.active_type, options.virtual_type]
                    exceptions_list = [inds_ccae, inds_caae, inds_aaae]

                    exception_symm = symmetry((0,1,3,2), -1)

                    unordered_tensor_inds = [get_spatial_index_type(ind) for ind in ordered_tensor.indices]

                    if (unordered_tensor_inds in exceptions_list) and (exception_symm in ordered_tensor.symmetries):
                        ordered_tensor.indices = [ordered_tensor.indices[i] for i in [0, 1, 3, 2]]
                        order_factor *= - 1.0

                if options.verbose:
                    print("")
                    print(unordered_term)
                    print(' %s    --->    %s (factor = %s)' % (unordered_tensor, ordered_tensor, order_factor))

                _terms[ind_unordered_term].tensors[ind_unordered_tensor] = ordered_tensor.copy()
                _terms[ind_unordered_term].numConstant *= order_factor

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return _terms

