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
# Author: Eric Neuscamman <eric.neuscamman@gmail.com>
#

from .sqaTensor import kroneckerDelta, creOp, desOp, sfExOp
from .sqaTerm import term, sortOps
from .sqaMisc import makeTuples, allDifferent, makePermutations


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def normalOrder(inTerm):
    "Returns a list of terms resulting from normal ordering the operators in inTerm."

    # check that inTerm is a term
    if not isinstance(inTerm, term):
        raise TypeError("inTerm must be of class term, not %s" % type(inTerm))

    # determine what types of operators the term contains
    has_creDesOps = any(isinstance(t, (creOp, desOp)) for t in inTerm.tensors)
    has_sfExOps = any(isinstance(t, sfExOp) for t in inTerm.tensors)

    # If term has both creation/destruction operators and spin free excitation operators, raise an error
    if has_creDesOps and has_sfExOps:
        raise RuntimeError("Normal ordering not implemented when both creOp/desOp and sfExOp tensors are present")

    # if the term is already normal ordered, return it unchanged
    if inTerm.isNormalOrdered():
        return [inTerm.copy()]

    # Normal ordering for creOp/desOp
    if has_creDesOps:

        # Separate the cre/des operators from other tensors
        ops = [t.copy() for t in inTerm.tensors if isinstance(t, (creOp, desOp))]
        nonOps = [t.copy() for t in inTerm.tensors if not isinstance(t, (creOp, desOp))]

        # Generate all contraction pairs
        contractionPairs = [
            (i, j)
            for i in range(len(ops))
            for j in range(i+1, len(ops))
            if isinstance(ops[i], desOp) and isinstance(ops[j], creOp)
        ]

        # Determine maximum contraction order
        creCount = 0
        maxConOrder = 0
        for iTerm in reversed(ops):
            if isinstance(iTerm, creOp):
                creCount += 1
            elif isinstance(iTerm, desOp) and creCount > 0:
                maxConOrder += 1
                creCount -= 1

        # Generate all contractions
        contractions = []
        for i in range(maxConOrder + 1):
            subCons = makeTuples(i, contractionPairs)
            # Store only valid contractions (no duplicate tags)
            contractions.extend([
                con for con in subCons
                if allDifferent([con[k][1] for k in range(i)]) and
                   allDifferent([con[k][0] for k in range(i)])
            ])

        # For each contraction, generate the resulting term
        outTerms = []
        for contraction in contractions:
            conSign = 1
            deltaFuncs = []
            subOpString = list(ops)
            for conPair in contraction:
                index1 = ops[conPair[0]].indices[0]
                index2 = ops[conPair[1]].indices[0]
                deltaFuncs.append(kroneckerDelta([index1, index2]))
                subOpString[conPair[0]] = 'contracted'
                subOpString[conPair[1]] = 'contracted'
                # Count sign flips
                conSign *= (-1) ** sum(1 for q in subOpString[conPair[0]+1:conPair[1]] if q != 'contracted')

            # Remove contracted operators
            subOpString = [op for op in subOpString if op != 'contracted']

            # Form outTerms
            sortSign, sortedOps = sortOps(subOpString)
            totalSign = conSign * sortSign
            outTensors = nonOps + deltaFuncs + sortedOps
            outTerms.append(term(totalSign * inTerm.numConstant, inTerm.constants, outTensors))

    # Normal ordering for sfExOps
    elif has_sfExOps:

        # Make separate lists of the spin free excitation operators and other tensors
        sfExOp_list = [t.copy() for t in inTerm.tensors if isinstance(t, sfExOp)]
        other_list = [t.copy() for t in inTerm.tensors if not isinstance(t, sfExOp)]

        # Initialize n, the number of remaining spin free excitation operators
        n = len(sfExOp_list)

        # Set the original term, with all excitation operators moved to the end, as the
        # first iteration's input term
        iter_input_terms = [term(inTerm.numConstant, inTerm.constants, other_list + sfExOp_list)]

        # Successively normal order the last two excitation operators until each term
        # has only one exitation operator left (at which point the term is normal ordered)
        while n > 1:

            # Initialize the list to hold this iteration's output terms
            iter_output_terms = []

            # For each of this iteration's input terms, produce all terms resulting from normal ordering
            # the last two excitation operators
            for t in iter_input_terms:

                # Make a list of the term's tensors that excludes the last two excitation operators
                tensors_except_last_two = t.tensors[:-2]

                # Give short names for the last two excitation operators and their orders
                e1, e2 = t.tensors[-2], t.tensors[-1]
                o1, o2 = e1.order, e2.order

                # Loop over the number of contractions
                for nc in range(min(o1,o2)+1):

                    # Compute the order of excitation operator for the current number of contractions
                    newOrder = o1 + o2 - nc

                    # Compute all nc-tuples of index numbers from e1 and e2, as well as all permutations of 
                    # the order in which the tuples may be combined to form a contraction
                    perms = [0] if nc == 0 else makePermutations(nc)
                    tups1 = makeTuples(nc, range(o1))
                    tups2 = makeTuples(nc, range(o2))

                    # For each contraction, compute the resulting term
                    for perm in perms:
                        for tup1 in tups1:
                            for tup2 in tups2:

                                # Initialize the term's tensor list
                                tensorList = list(tensors_except_last_two)

                                # Compute the pairs of indices to be contracted
                                # Example: (conPairs[0][p], conPairs[1][p]) is the pth contraction pair.
                                conPairs = [
                                    [tup1[perm[i]] for i in range(nc)],
                                    [tup2[i] for i in range(nc)]
                                ]

                                # Initialize the index list for the new excitation operator
                                indexList = [None] * (2*newOrder)

                                # Populate the index list for the new excitation operator
                                # Also, create a kronecker delta function for each contraction pair
                                # Indices from e1
                                for i in range(o1):
                                    indexList[i] = e1.indices[i]
                                    if i in conPairs[0]:
                                        i1 = i+o1
                                        i2 = conPairs[1][conPairs[0].index(i)]
                                        ind1 = e1.indices[i1]
                                        ind2 = e2.indices[i2]
                                        tensorList.insert(0, kroneckerDelta([ind1,ind2]))
                                        indexList[i+newOrder] = e2.indices[o2+i2]
                                    else:
                                        indexList[i+newOrder] = e1.indices[i+o1]
                                # Indices from e2
                                count = 0
                                for i in range(o2):
                                    if not (i in conPairs[1]):
                                        indexList[o1+count] = e2.indices[i]
                                        indexList[o1+count+newOrder] = e2.indices[i+o2]
                                        count += 1

                                # Ensure that all slots in the index list have been filled
                                if any(ind is None for ind in indexList):
                                    raise RuntimeError("There is at least one unassigned index in the new spin free operator.")

                                # Add the new excitation operator to the tensor list
                                tensorList.append(sfExOp(indexList))

                                # Add the resulting term to this iteration's list of output terms
                                iter_output_terms.append(term(inTerm.numConstant, inTerm.constants, tensorList))

            # Set this iteration's list of output terms as the next iteration's input terms
            iter_input_terms = iter_output_terms

            # Decrement the counter for the number of excitation operators
            n -= 1

        # Set the return value as the final iteration's output terms
        outTerms = iter_output_terms

    else:
        raise RuntimeError("Normal ordering function failed to choose what to do.")

#    print "Terms after normal ordering:"
#    for t in outTerms:
#        print t

    return outTerms


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


