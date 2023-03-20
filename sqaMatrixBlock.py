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
# Author: Koushik Chatterjee <koushikchatterjee7@gmail.com>
#

import sys, time
from sqaTensor import kroneckerDelta, sfExOp, creOp, desOp
from sqaTerm import term, termChop, combineTerms
from sqaOptions import options
from sqaMisc import makeTuples, allDifferent

from sqaNormalOrder import normalOrder

from sqaIndex import get_spatial_index_type, get_spin_index_type, \
                     is_core_index_type, is_active_index_type, is_virtual_index_type

#####################################
def matrixBlock(terms, transRDM = False):
    "Construct matrix block."

    startTime = time.time()
    print("")
    print("---------------------------------- SQA Automation --------------------------------")
    sys.stdout.flush()
    print("")
    fTerms = []

    # Make sure that input expression is normal-ordered with respect to physical vacuum
    nterms = []
    for t in terms:
        t_no = normalOrder(t)
        nterms.extend(t_no)

    # Filter zero terms wrt virtual (note: Filter first for virtual orbitals)
    filterVirtual(nterms)
    termChop(nterms)

    # Normal ordering with respect to core orbitals
    fTerms = normalOrderCore(nterms)

    # Evaluate Kroneker delta
    for t in fTerms:
        t.contractDeltaFuncs()

    # Filter zero terms wrt core (note: after filtering virtual, then do for core)
    filterCore(fTerms)
    termChop(fTerms)
    combineTerms(fTerms)

    # Contract delta function for both non-dummy indices
    contractDeltaFuncs_nondummy(fTerms)

    # If (transRDM = True) =>  Remove those constant terms
    if (transRDM):
        for trm in fTerms:
            iremove = True
            for i in range(len(trm.tensors)):
                t = trm.tensors[i]
                if (isinstance(t, creOp) or isinstance(t, desOp)):
                   iremove = False
            if (iremove):
                trm.numConstant = 0.0
        termChop(fTerms)

    # Dummy indices label upate
    fTerms = dummyLabel(fTerms)

    # Reorder tensor indices: (core < active < virtual) order
    reorder_tensor_indices(fTerms, reorder_t = True)

    # Print the final results
    print("")
    print("--------------------------------- Final results ----------------------------------\n")
    sys.stdout.flush()
    for t in fTerms:
        index_types = ()
        for t_tensor in t.tensors:
            for t_tensor_index in t_tensor.indices:
                index_types += t_tensor_index.indType[0]
        print(t)

    print("")
    print("Total terms : %s" % len(fTerms))
    print("SQA automation time :  %.3f seconds" % (time.time() - startTime))
    sys.stdout.flush()

    return fTerms

def dummyLabel(_terms):
    "A function to relabel dummy indices."

    print("Dummy indices relabelling...")
    sys.stdout.flush()

    nterms = list(_terms)
    for t in nterms:
        mymap = {}

        coreInd = list('ijklmnopq')
        actvInd = list('xyzwuvstr')
        virtInd = list('abcdefgh')

        reservedInd = []
        for t_tensor in t.tensors:
            for t_tensor_index in range(len(t_tensor.indices)):
                index_name = t_tensor.indices[t_tensor_index].name
                index_user_defined = t_tensor.indices[t_tensor_index].userDefined

                if index_name not in reservedInd and index_user_defined:
                    reservedInd.append(index_name)

        for reserved_index_name in reservedInd:
            if reserved_index_name in coreInd:
                coreInd.remove(reserved_index_name)
            elif reserved_index_name in actvInd:
                actvInd.remove(reserved_index_name)
            elif reserved_index_name in virtInd:
                virtInd.remove(reserved_index_name)

        for t_tensor in t.tensors:
            for t_tensor_index in range(len(t_tensor.indices)):

                # Decide which new label to assign
                index_type = t_tensor.indices[t_tensor_index].indType
                index_name = t_tensor.indices[t_tensor_index].name
                index_user_defined = t_tensor.indices[t_tensor_index].userDefined

                if not index_user_defined:
                        if index_name not in mymap.keys():
                            if is_core_index_type(index_type):
                                mymap[index_name] = coreInd[0]
                                coreInd.pop(0)
                            elif is_active_index_type(index_type):
                                mymap[index_name] = actvInd[0]
                                actvInd.pop(0)
                            elif is_virtual_index_type(index_type):
                                mymap[index_name] = virtInd[0]
                                virtInd.pop(0)

                    # Update the label
                        t_tensor.indices[t_tensor_index].name = mymap[index_name]

    if options.verbose:
        print("")
        for _term in nterms:
            print(_term)
        print("")

    print("Done!")
    print("----------------------------------------------------------------------------------")
    sys.stdout.flush()
    return nterms

def filterVirtual(_terms):
    "A function to calculate expectation value wrt virtual: filter zero terms wrt virtual."

    print("Computing expectation value with respect to virtual:=>")
    sys.stdout.flush()

    for t_term in _terms:
        for t_tensor in t_term.tensors:
            for t_tensor_index in range(len(t_tensor.indices)):

                index_type = t_tensor.indices[t_tensor_index].indType

                if isinstance(t_tensor, desOp):
                    if is_virtual_index_type(index_type):
                        t_term.numConstant = 0.0
                elif isinstance(t_tensor, creOp):
                    if is_virtual_index_type(index_type):
                        t_term.numConstant = 0.0

    if options.verbose:
        print("")
        for t_term in _terms:
            print(t_term)
        print("")

    print("Done!")
    print("----------------------------------------------------------------------------------")
    sys.stdout.flush()
    return

def filterCore(_terms):
    "A function to calculate expectation value wrt core: filter zero terms wrt core."

    print("Computing expectation value with respect to core:=>")
    sys.stdout.flush()

    for t_term in _terms:
        for t_tensor in t_term.tensors:
            for t_tensor_index in range(len(t_tensor.indices)):

                index_type = t_tensor.indices[t_tensor_index].indType

                if isinstance(t_tensor, desOp):
                    if is_core_index_type(index_type):
                        t_term.numConstant = 0.0
                elif isinstance(t_tensor, creOp):
                    if is_core_index_type(index_type):
                        t_term.numConstant = 0.0

    if options.verbose:
        print("")
        for _term in _terms:
            print(_term)
        print("")

    print("Done!")
    print("----------------------------------------------------------------------------------")
    sys.stdout.flush()
    return

def normalOrderCore(_terms):
    print("Normal ordering with respect to core:=>")
    sys.stdout.flush()
    ordered_terms = []

    for _term in _terms:
        ordered_term = normOrderCor(_term)
        ordered_terms.extend(ordered_term)

    print("Done!")
    print("----------------------------------------------------------------------------------")
    sys.stdout.flush()
    return ordered_terms

def normOrderCor(_term):
    "Returns a list of terms resulting from normal ordering the operators in inTerm."
    if options.verbose:
        print('Term = %s' % _term)

    # check if is a term
    if not isinstance(_term, term):
        raise TypeError ("Input term must be of class term")

    # determine what types of operators the term contains
    has_creDesOps = False
    has_sfExOps = False
    for _tensor in _term.tensors:
        if isinstance(_tensor, creOp) or isinstance(_tensor, desOp):
            has_creDesOps = True
        elif isinstance(_tensor, sfExOp):
            has_sfExOps = True

    # If term has both creation/destruction operators and spin free excitation operators raise an error
    if has_creDesOps and has_sfExOps:
        raise RuntimeError ("Normal ordering not implemented when both creOp/desOp and sfExOp tensors are present")

    # Normal ordering for creOp/desOp
    elif has_creDesOps:

        # Separate the cre/des operators from other tensors
        ops = []
        nonOps = []
        for _tensor in _term.tensors:
            if isinstance(_tensor, creOp) or isinstance(_tensor, desOp):
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
        raise Exception ('This code does not support for now')

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
    Sorts a list of creation/destruction operators into normal order and alphabetically.
    Performs no contractions.  Returns the overall sign resulting from the sort and the sorted operator list.
    """
    sorted_ops = _unsorted_ops + []
    i = 0
    sign = 1

    if returnPermutation:
        perm = range(len(_unsorted_ops))

    while i < len(sorted_ops)-1:
        if isinstance(sorted_ops[i], creOp) and is_core_index_type(sorted_ops[i].indices[0]):
            j = i
            for k in range(i,(len(sorted_ops)-1)):
                if isinstance(sorted_ops[k+1], creOp) and is_core_index_type(sorted_ops[k+1].indices[0]):
                    i += 1
                    j = i
                else:
                    i = j
                    break
            if ((i+1) > (len(sorted_ops)-1)):
                break

            temp = sorted_ops[i]
            sorted_ops[i] = sorted_ops[i+1]
            sorted_ops[i+1] = temp

            if returnPermutation:
                temp = perm[i]
                perm[i] = perm[i+1]
                perm[i+1] = temp
            i = 0
            sign *= -1

        elif isinstance(sorted_ops[i+1], desOp) and is_core_index_type(sorted_ops[i+1].indices[0]):
            temp = sorted_ops[i]
            sorted_ops[i] = sorted_ops[i+1]
            sorted_ops[i+1] = temp

            if returnPermutation:
                temp = perm[i+1]
                perm[i+1] = perm[i]
                perm[i] = temp
            i = 0
            sign *= -1
            if (sorted_ops[i+1].name == sorted_ops[i].name):
                break
        else:
            i += 1

    if returnPermutation:
        return (sign, sorted_ops, perm)
    return (sign, sorted_ops)

def print_header():

    print("""
    \n----------------------------------------------------------------------------------
    sqa_plus: Code generator for quasi-particle systems.
    Copyright 2009-2022 SecondQuantizationAlgebra Developers. All Rights Reserved.

    Licensed under the GNU General Public License v3.0;

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

    ----------------------------------------------------------------------------------""")

def contractDeltaFuncs_nondummy(_terms):
    "Contracts delta function for both non-dummy indices only wrt to orbitals subspaces, otherwise use 'contractDeltaFuncs' function."

    print("Contract delta function for non-dummy indices: =>")
    sys.stdout.flush()

    for term in _terms:
        for i in range(len(term.tensors)):
            t = term.tensors[i]

            if isinstance(t, kroneckerDelta):
                i0 = t.indices[0]
                i1 = t.indices[1]
                if not (i0.isSummed and i1.isSummed):
                    if not ((get_spatial_index_type(i0) == get_spatial_index_type(i1)) and
                            (get_spin_index_type(i0) == get_spin_index_type(i1))):
                        term.numConstant = 0.0

    termChop(_terms)

    print("Done!")
    print("----------------------------------------------------------------------------------")
    sys.stdout.flush()
    return _terms

def reorder_tensor_indices(_terms, reorder_t = True):
    print("Reordering indices according to core < active < virtual...")
    sys.stdout.flush()

    for ind_unordered_term, unordered_term in enumerate(_terms):
        for ind_unordered_tensor, unordered_tensor in enumerate(unordered_term.tensors):
            if reorder_t:
                reorder_tensor = ((unordered_tensor.name == 'h' and len(unordered_tensor.indices) == 2) or
                                  (unordered_tensor.name == 'v' and len(unordered_tensor.indices) == 4) or
                                  ((unordered_tensor.name[0] == 't') and
                                   ((len(unordered_tensor.indices) == 2) or len(unordered_tensor.indices) == 4)))
            else:
                reorder_tensor = ((unordered_tensor.name == 'h' and len(unordered_tensor.indices) == 2) or
                                  (unordered_tensor.name == 'v' and len(unordered_tensor.indices) == 4))

            if reorder_tensor:
                original_rank = []
                for ind in unordered_tensor.indices:
                    if is_core_index_type(ind):
                        original_rank.append(3)
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
                        if is_core_index_type(ind):
                            permute_rank.append(3)
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

                _terms[ind_unordered_term].tensors[ind_unordered_tensor] = ordered_tensor.copy()

                order_factor = permutes_factors[permutes_rank_ind]
                _terms[ind_unordered_term].numConstant *= order_factor

                if options.verbose:
                    print("")
                    print(unordered_term)
                    print(' %s    --->    %s (factor = %s)' % (unordered_tensor, ordered_tensor, order_factor))

    print("Done!")
    print("----------------------------------------------------------------------------------")
    sys.stdout.flush()
    return _terms

