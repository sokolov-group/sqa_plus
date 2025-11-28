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
# The term class represents a set of constants and tensors that have been multiplied together.
# The class consists of a numerical constant factor, a list of named constants, and a list
# of tensors.
#
# It is common to encounter a list of term objects, which is often used to represent an expression.
# For example, the Hamiltonian in second quantization can be written as a list of terms.
#

from functools import total_ordering
from multiprocessing import Pool, cpu_count
from collections import deque
from itertools import islice, count

from .sqaIndex import index
from .sqaTensor import tensor, kroneckerDelta, sfExOp, creOp, desOp
from .sqaMisc import makePermutations
from .sqaOptions import options
import time
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


@total_ordering
class term:
    "A class for terms used in operator algebra. Each term is a multiplicative string of constants, tensors, and operators."

    #------------------------------------------------------------------------------------------------

    def __init__(self, numConstant, constList, tensorList, isInCanonicalForm = False):
        self.constants = []
        self.tensors = []
        if isinstance(numConstant, (int, float)):
            self.numConstant = float(numConstant)
        else:
            raise TypeError("numConstant must be given as a float or an int.")

        if not all(isinstance(c, str) for c in constList):
            raise TypeError("constList must be a list of strings.")
        self.constants.extend(constList)

        if not all(isinstance(t, tensor) for t in tensorList):
            raise TypeError("tensorList must be a list of tensor objects.")
        self.tensors.extend(t.copy() for t in tensorList)

        if isinstance(isInCanonicalForm, bool):
            self.isInCanonicalForm = isInCanonicalForm
        else:
            raise TypeError("if specified, isInCanonicalForm must be True or False")

    #------------------------------------------------------------------------------------------------

    def __eq__(self, other):
        if not isinstance(other, term):
            return False

        # sort by number of loose creation operators first
        if self.nCreOps() != other.nCreOps():
            return False
        # next sort by number of loose destruction operators
        if self.nDesOps() != other.nDesOps():
            return False
        # next sort by the orders of the spin free excitation operators
        if self.sfExOp_ranks() != other.sfExOp_ranks():
            return False
        # next sort by the number of constants
        if len(self.constants) != len(other.constants):
            return False
        # next sort by the number of tensors
        if len(self.tensors) != len(other.tensors):
            return False
        # next sort by the tensors' names
        if [t.name for t in self.tensors] != [t.name for t in other.tensors]:
            return False
        # next sort by the tensors
        if self.tensors != other.tensors:
            return False
        # next sort by the constants
        if self.constants != other.constants:
            return False
        # finally compare the numerical constants
        numDiff = self.numConstant - other.numConstant
        return abs(numDiff) < 1e-6

    def __lt__(self, other):
        if not isinstance(other, term):
            raise TypeError("term object can only be compared to other term objects.")

        # sort by number of loose creation operators first
        if self.nCreOps() != other.nCreOps():
            return self.nCreOps() < other.nCreOps()
        # next sort by number of loose destruction operators
        if self.nDesOps() != other.nDesOps():
            return self.nDesOps() < other.nDesOps()
        # next sort by the orders of the spin free excitation operators
        if self.sfExOp_ranks() != other.sfExOp_ranks():
            return self.sfExOp_ranks() < other.sfExOp_ranks()
        # next sort by the number of constants
        if len(self.constants) != len(other.constants):
            return len(self.constants) < len(other.constants)
        # next sort by the number of tensors
        if len(self.tensors) != len(other.tensors):
            return len(self.tensors) < len(other.tensors)
        # next sort by the tensors' names
        self_names = [t.name for t in self.tensors]
        other_names = [t.name for t in other.tensors]
        if self_names != other_names:
            return self_names < other_names
        # next sort by the tensors
        if self.tensors != other.tensors:
            return self.tensors < other.tensors
        # next sort by the constants
        if self.constants != other.constants:
            return self.constants < other.constants
        # finally compare the numerical constants
        numDiff = self.numConstant - other.numConstant
        if abs(numDiff) >= 1e-6:
            return numDiff < 0

        return False

    #------------------------------------------------------------------------------------------------

    def __str__(self):

        retval = [f" ({self.numConstant:10.5f})", " ".join(self.constants), " ".join(str(t) for t in self.tensors)]
        return " ".join(filter(None, retval))

    #------------------------------------------------------------------------------------------------

    def __add__(self,other):
        "Adds the terms self and other. Only works if the constants, tensors, and operators of two terms all match."
        if not self.sameForm(other):
            raise RuntimeError("self and other must have the same constants, tensors, and operators")
        retval = self.copy()
        retval.numConstant += other.numConstant
        return retval

    #------------------------------------------------------------------------------------------------

    def copy(self):
        "Returns a deep copy of the term."
        return term(self.numConstant, self.constants, self.tensors, self.isInCanonicalForm)

    #------------------------------------------------------------------------------------------------

    def nCreOps(self):
        "Returns the number of loose creation operators in the term"
        retval = sum(isinstance(t, creOp) for t in self.tensors)
        return retval

    #------------------------------------------------------------------------------------------------

    def nDesOps(self):
        "Returns the number of loose destruction operators in the term"
        retval = sum(isinstance(t, desOp) for t in self.tensors)
        return retval

    #------------------------------------------------------------------------------------------------

    def sfExOp_ranks(self):
        "Returns a list of the ranks of the spin free excitation operators in the term"
        retval = [len(t.indices)/2 for t in self.tensors if isinstance(t, sfExOp)]
        return retval

    #------------------------------------------------------------------------------------------------

    def scale(self, factor):
        "Multiplies the term by factor."
        if not isinstance(factor, (int, float)):
            raise ValueError("factor must an integer or float")
        self.numConstant *= factor

    #------------------------------------------------------------------------------------------------

    def sameForm(self, other):
        "Determines whether the terms are of the same form."
        self.makeCanonical()
        other.makeCanonical()
        if self.constants != other.constants or self.tensors != other.tensors:
            return False
        return True

    #------------------------------------------------------------------------------------------------

    def contractDeltaFuncs(self):
        "Contracts the indices within any kronecker delta functions in the term."

        i = 0
        while i < len(self.tensors):
            t = self.tensors[i]
            
            # If the term is a delta funciton with a repeated index, remove it
            if isinstance(t, kroneckerDelta) and (t.indices[0] == t.indices[1]) and (t.indices[0].userDefined == t.indices[1].userDefined):
                del(self.tensors[i])

            # If the term is a delta function with contractable indices, contract them and remove the delta func
            elif isinstance(t, kroneckerDelta) and (t.indices[0].isSummed or t.indices[1].isSummed):
                i0 = t.indices[0]
                i1 = t.indices[1]

                # Check that the indices have the same number of type groups
                if len(i0.indType) != len(i1.indType):
                    raise RuntimeError("Cannot contract indices %s, %s.    They have different numbers of type groups." %(str(i0), str(i1)))

                # Determine the type of the new index based on the overlap between
                # the types of the two indices
                typeOverlap = []
                for j in range(len(i0.indType)):
                    typeOverlap.append([])
                    for typeString in i0.indType[j]:
                        if typeString in i1.indType[j]:
                            typeOverlap[-1].append(typeString)

                # If there is no overlap between any of the type groups, the delta function is zero
                if ( len(i0.indType) > 0 or len(i1.indType) > 0 ) and [] in typeOverlap:
                    self.numConstant = 0.0
                    return

                # Create the new index
                if not i0.isSummed:
                    newIndex = index(i0.name, typeOverlap, i0.isSummed, i0.userDefined)
                else:
                    newIndex = index(i1.name, typeOverlap, i1.isSummed, i1.userDefined)

                # Remove the delta function
                del(self.tensors[i])

                # Excecute the index replacement
                for ten in self.tensors:
                    for j in range(len(ten.indices)):
                        if ten.indices[j] in [i0,i1]:
                            ten.indices[j] = newIndex.copy()

            # Otherwise move on to the next tensor
            else:
                i += 1


    #------------------------------------------------------------------------------------------------

    def generateAlphabet(self):
        """ Returns a list of strings that share no elements with the term's index names. """

        # Set of index names used in self
        used_names = {
            (f"user_{idx.userDefined}" if idx.userDefined else idx.name)
            for t in self.tensors
            for idx in t.indices
        }

        # Generate an alphabet with no elements overlapping with used_names
        alphabet = list(
            islice((idx for idx in map(str, count()) if idx not in used_names), len(used_names))
        )

        return alphabet

    #------------------------------------------------------------------------------------------------

    def isNormalOrdered(self):
        """ Returns True if term is normal-ordered, False otherwise. """
        seen_tensor_types = set()
        for t in self.tensors:
            if isinstance(t, (creOp, sfExOp)) and seen_tensor_types.intersection((desOp, sfExOp)):
                return False
            seen_tensor_types.add(type(t))
        return True

    #------------------------------------------------------------------------------------------------

    def makeCanonical(self, rename_user_defined = True):
        """ Converts the term to a unique canonical form. """
        # Use the non recursive function
        self.makeCanonical_non_recursive(rename_user_defined)
        return

    ## DEAD CODE
    def makeCanonical_recursive(self, rename_user_defined = True):
        "Converts the term to a unique canonical form using a recursive algorithm."

        # If the tensor is already in canonical form, do nothing
        if self.isInCanonicalForm:
            return

        # If the term is not normal ordered, raise an error
        if not self.isNormalOrdered():
            raise RuntimeError("A term must have normal ordered operators to be converted to canonical form.")

        # Sort the constants
        self.constants.sort()

        # If there are no tensors in the term then skip the tensor sorting.
        if len(self.tensors) == 0:
            self.isInCanonicalForm = True
            return

        # Create all tensor lists in which the tensors are sorted by type and name.
        # There can be more than one term here if there are multiple tensors with the same name.
        candidateTensorLists = self.getCandidateTensorLists()

        # Generate an alphabet for use in renaming indices
        # This is done to avoid renaming with an index name already in use.
        # Should try using numbers, i.e. '1', '2', '3', etc.
        alphabet = self.generateAlphabet()

        # For each candidate list, rename the indices canonically and record the score for that list.
        bestScore = [-1]
        nTopScore = 0
        for i in range(len(candidateTensorLists)):
            (score,map,sign,newTensorList) = getcim(candidateTensorLists[i],alphabet)
            if score > bestScore:
                (bestScore,bestMap,bestSign,bestTensorList) = (score,map,sign,newTensorList)
                nTopScore = 1
            elif score == bestScore:
                nTopScore += 1

        # Check to see that only one candidate achieved the top score
        if nTopScore > 1:
            raise RuntimeError("%i candidates tied for the top score." %nTopScore)

        # Set the tensors and their indices in the canonical order (the order with the highest score)
        self.tensors = bestTensorList

        # Apply the sign produced from the index sorting
        self.scale(bestSign)

        # Create an index mapping that converts to a canonical alphabet, i.e. a-z
        map = bestMap
        alphabet = list('abcdefghijklmnopqrstuvwxyz')
        if len(alphabet) < len(map):
            raise RuntimeError("Alphabet smaller than number of indices, no more names left!")
        canonMap = {}
        while map.keys():
            minVal = min(map.values())
            for key in map.keys():
                if map[key] == minVal:
                    canonMap[map[key].tup()] = map[key].copy()
                    canonMap[map[key].tup()].name = alphabet.pop(0)
                    del map[key]
                    break
        map = canonMap

        # Rename the indices using the canonical mapping
        for t in self.tensors:
            for i in range(len(t.indices)):
                if t.indices[i].tup() in map.keys():
                    t.indices[i] = map[t.indices[i].tup()].copy()

        # Turn on the canonical form flag to avoid calling this function again unnecessarily
        self.isInCanonicalForm = True

    #------------------------------------------------------------------------------------------------

    def makeCanonical_non_recursive(self, rename_user_defined = True):
        "Converts the term to a unique canonical form using a non-recursive algorithm."

        # If the tensor is already in canonical form, do nothing
        if self.isInCanonicalForm:
            return

        # If the term is not normal ordered, raise an error
        if not self.isNormalOrdered():
            raise RuntimeError("A term must have normal ordered operators to be converted to canonical form.")

        # Sort the constants
        self.constants.sort()

        # If there are no tensors in the term then skip the tensor sorting.
        if not self.tensors:
            self.isInCanonicalForm = True
            return

        # Sort the freely commuting tensors by name
        fcList = [t for t in self.tensors if t.freelyCommutes]
        ncList = [t for t in self.tensors if not t.freelyCommutes]

        nameGroups = {}
        for t in fcList:
            if t.name not in nameGroups:
                nameGroups[t.name] = []
            nameGroups[t.name].append(t)

        # Sort further by tensor subclass and length (use name as tie-breaker)
        ##sort_key = lambda item: (item[1][0].__class__.__name__, len(item[1]), item[0], tuple(ind.indType for ind in item[1][0].indices))
        sort_key = lambda item: (item[1][0].__class__.__name__, len(item[1]), item[0])
        nameGroups = [group for _, group in sorted(nameGroups.items(), key=sort_key)]

        # Add non-commuting tensors as individual groups
        nameGroups.extend([[t] for t in ncList])
 
        # Sort within all groups by index type and name
        sort_key = lambda t: (tuple(ind.indType for ind in t.indices), tuple(str(ind.name) for ind in t.indices))
        ##sort_key = lambda t: (t.__class__.__name__, len(t.indices), t.name, tuple(ind.indType for ind in t.indices), tuple(str(ind.name) for ind in t.indices))
        for group in nameGroups:
            group.sort(key=sort_key)

        # Generate an alphabet for use in renaming indices
        # This is done to avoid renaming with an index name already in use.
        # Should try using numbers, i.e. '1', '2', '3', etc.
        alphabet = self.generateAlphabet()

        # Determine the best ordering and index mapping
        best_tensor_list = None
        best_factor = None
        bestMap = None
        bestScore = None
        nTopScore = 0
        # job format:    (map, gCount, tCount, aCount, gPerms)
        #jobStack = [({},0,0,0,[])]
        jobStack = deque([({}, 0, 0, 0, [])])
        while jobStack:

            # get the next job
            current_map, gCount, tCount, aCount, gPerms = jobStack.pop()

            # If there are no name groups remaining, compute the score
            if gCount == len(nameGroups):

                # Compute a new tensor list in which any dummy indices are given their
                # new names and all indices are sorted.
                # Also compute a list of these ordered indices.
                # Also compute the multiplicitive factor generated by sorting the indices.
                tenList, indexList, factor = self.get_tensor_list(nameGroups, gPerms, current_map)

                # Compute a score based on how alphabetical the indices are
                score = []
                n_ind = len(indexList)
                for i in range(n_ind-1):
                    count = 0
                    for j in range(i+1, n_ind):
                        if indexList[i] < indexList[j]:
                            count += 1
                    score.append(count)

                # If the current score is the best score, save the result
                if (bestScore is None) or (score > bestScore):
                    nTopScore = 1
                    bestScore = score
                    bestMap = current_map
                    best_factor = factor
                    best_tensor_list = tenList

                # If the current score ties for the best, count the number of best scores
                elif score == bestScore:
                    nTopScore += 1

                continue

            # If only cre/des operators remain, sort them and compute the score
            if all([ (len(group) == 1 and isinstance(group[0], (creOp, desOp))) for group in nameGroups[gCount:] ]):
                # Compute a new tensor list in which any dummy indices are given their
                # new names and all indices are sorted.
                # Also compute a list of these ordered indices.
                # Also compute the multiplicitive factor generated by sorting the indices.
                tenList, indexList, factor = self.get_tensor_list(nameGroups[:gCount], gPerms, current_map)

                # Apply the input mapping to the creation/destruction operators
                # Create and apply a new mapping for any new dummy indices
                opList = [nameGroups[i][0].copy() for i in range(gCount,len(nameGroups))]
                nNewMaps = 0
                for op in opList:
                    if op.indices[0].tup() in current_map:
                        op.indices[0] = current_map[op.indices[0].tup()]
                    elif op.indices[0].isSummed:
                        temp_map = current_map.copy()
                        temp_map[op.indices[0].tup()] = index(alphabet[aCount+nNewMaps], op.indices[0].indType, op.indices[0].isSummed, op.indices[0].userDefined)
                        current_map = temp_map
                        nNewMaps += 1
                        op.indices[0] = current_map[op.indices[0].tup()]

                # Sort the operators and apply the resulting sign
                (sign, opList) = sortOps(opList)
                factor *= sign

                # Add the operators' indices to the ordered list of indices.
                # Also add the sorted operators to the new tensor list.
                for op in opList:
                    indexList.append(op.indices[0])
                    tenList.append(op)

                # Compute a score based on how alphabetical the indices are
                score = []
                n_ind = len(indexList)
                for i in range(n_ind-1):
                    count = 0
                    for j in range(i+1, n_ind):
                        if indexList[i] < indexList[j]:
                            count += 1
                    score.append(count)

                # If the current score is the best score, save the result
                if (bestScore is None) or (score > bestScore):
                    nTopScore = 1
                    bestScore = score
                    bestMap = current_map
                    best_factor = factor
                    best_tensor_list = tenList

                # If the current score ties for the best, count the number of best scores
                elif score == bestScore:
                    nTopScore += 1

                continue

            # If only a sfExOp remains, sort its indices and compute the score
            #group_is_last = (gCount == len(nameGroups) - 1)
            if (gCount == len(nameGroups)-1) and (len(nameGroups[gCount]) == 1) and isinstance(nameGroups[gCount][0], sfExOp):

                # Compute a new tensor list in which any dummy indices are given their
                # new names and all indices are sorted.
                # Also compute a list of these ordered indices.
                # Also compute the multiplicitive factor generated by sorting the indices.
                tenList, indexList, factor = self.get_tensor_list(nameGroups[:gCount], gPerms, current_map)

                # Apply the input mapping to the sfExOp
                # Create and apply a new mapping for any new dummy indices
                t = nameGroups[gCount][0].copy()
                nNewMaps = 0
                for i in range(len(t.indices)):
                    if t.indices[i].tup() in current_map:
                        t.indices[i] = current_map[t.indices[i].tup()]
                    elif t.indices[i].isSummed:
                        temp_map = current_map.copy()
                        temp_map[t.indices[i].tup()] = index(alphabet[aCount+nNewMaps], t.indices[i].indType, t.indices[i].isSummed, t.indices[i].userDefined)
                        current_map = temp_map
                        nNewMaps += 1
                        t.indices[i] = current_map[t.indices[i].tup()]

                # Sort the indices of the sfExOp (go go gadget bubble sort!)
                i = 0
                while i < t.order-1:
                    if t.indices[i] > t.indices[i+1]:
                        temp = t.indices[i+1]
                        t.indices[i+1] = t.indices[i]
                        t.indices[i] = temp
                        temp = t.indices[i+t.order+1]
                        t.indices[i+t.order+1] = t.indices[i+t.order]
                        t.indices[i+t.order] = temp
                        i = 0
                    else:
                        i += 1

                # Add the sfExOp to the new tensor list
                tenList.append(t)

                # Add the sfExOp's indices to the ordered index list
                for i in t.indices:
                    indexList.append(i)

                # Compute a score based on how alphabetical the indices are
                score = []
                n_ind = len(indexList)
                for i in range(n_ind-1):
                    count = 0
                    for j in range(i+1, n_ind):
                        if str(indexList[i].name) < str(indexList[j].name):
                            count += 1
                    score.append(count)

                # If the current score is the best score, save the result
                if (bestScore is None) or (score > bestScore):
                    nTopScore = 1
                    bestScore = score
                    bestMap = current_map
                    best_factor = factor
                    best_tensor_list = tenList

                # If the current score ties for the best, count the number of best scores
                elif score == bestScore:
                    nTopScore += 1

                continue

            # If starting a new name group, sort the group's names based on the input mapping
            # and schedule a new job for each permutation of the tensors with no index assignments
            if len(gPerms) <= gCount:
                withMapped = []
                withoutMapped = []
                for i in range(len(nameGroups[gCount])):
                    leastMapped = False
                    for ind in nameGroups[gCount][i].indices:
                        if (ind.tup() in current_map) and ((leastMapped is False) or (current_map[ind.tup()] < leastMapped)):
                            leastMapped = current_map[ind.tup()]
                    if leastMapped is False:
                        withoutMapped.append(i)
                    else:
                        withMapped.append((leastMapped,i))
                withMapped.sort(key=lambda x: x[0])
                withMapped = [i[1] for i in withMapped]
                if len(withoutMapped) <= 1:
                    jobStack.append((current_map,gCount,tCount,aCount,gPerms + [withMapped + withoutMapped]))
                else:
                    for perm in makePermutations(len(withoutMapped)):
                        new_gPerm = withMapped + [withoutMapped[i] for i in perm]
                        jobStack.append((current_map,gCount,tCount,aCount,gPerms + [new_gPerm]))

                continue

            # If continuing an existing group, schedule a new job for each equivelent ordering
            # of the current tensor's indices
            # Give the current tensor a convenient name
            t = nameGroups[gCount][gPerms[gCount][tCount]]

            # Get the tensor's symmetry permutations
            (symPerms,factors) = t.symPermutes()

            # Compute gCount and tCount for the next job
            next_gCount = gCount
            next_tCount = tCount + 1
            if next_tCount == len(nameGroups[gCount]):
                next_gCount += 1
                next_tCount = 0

            # For each of the tensor's symmetry-equivelent index orderings, create mappings for any
            # un-mapped dummy indices and schedule a new job
            for perm in symPerms:
                nNewMaps = 0
                newMap = {}
                newMap.update(current_map)
                for ind in [t.indices[perm[i]] for i in range(len(t.indices))]:
                    if ind.isSummed and ind.tup() not in newMap:
                        newMap[ind.tup()] = index(alphabet[aCount+nNewMaps], ind.indType, ind.isSummed, ind.userDefined)
                        nNewMaps += 1
                jobStack.append((newMap,next_gCount,next_tCount,aCount+nNewMaps,gPerms))

        # Check to see that only one candidate achieved the top score
        #if nTopScore > 1:
        #    print(f'WARNING: {nTopScore} candidates have tied for the top score.')

        if bestMap is None:
            bestMap = current_map

        # Create an index mapping that converts to a canonical alphabet, i.e. a-z
        alphabet = list('abcdefghijklmnopqrstuvwxyz')
        filtered_alphabet = [c for c in alphabet if c not in options.user_defined_indices]

        if len(filtered_alphabet) < len(bestMap):
            raise RuntimeError("Alphabet smaller than number of indices, no more names left!")

        canonMap = {}
        for _, val in sorted(bestMap.items(), key=lambda kv: kv[1].name):
            canonMap[val.tup()] = val.copy()
            canonMap[val.tup()].name = filtered_alphabet.pop(0)

        # Rename the indices using the canonical mapping
        for t in best_tensor_list:
            for i in range(len(t.indices)):
                if t.indices[i].tup() in canonMap:
                    t.indices[i] = canonMap[t.indices[i].tup()].copy()
                if rename_user_defined and t.indices[i].userDefined:
                    t.indices[i].rename()

        # Finalize results
        self.tensors = best_tensor_list 
        self.scale(best_factor)
        self.isInCanonicalForm = True


    def get_tensor_list(self, nameGroups, gPerms, current_map):
        """Compute a new tensor list with renamed dummy indices, a list of the ordered dummy indices,
        and the multiplicative factor generated by sorting the indices."""
        tenList = []
        indexList = []
        factor = 1
        # Build tensor and index lists using group permutations (gPerms)
        for g_ind in range(len(nameGroups)):
            group = nameGroups[g_ind]
            perm = gPerms[g_ind]
            for t_ind in range(len(group)):
                t_src = group[perm[t_ind]].copy()
                # Rename dummy indices using current_map
                for i in range(len(t_src.indices)):
                    if t_src.indices[i].isSummed:
                        t_src.indices[i] = current_map[t_src.indices[i].tup()]
                # Sort tensor indices and accumulate factor
                factor *= t_src.sortIndices()
                # Store sorted indices and tensor
                for ind in t_src.indices:
                    indexList.append(ind)
                tenList.append(t_src)
        return tenList, indexList, factor
    #------------------------------------------------------------------------------------------------

    def getCandidateTensorLists(self):
        """
        Create all tensor lists in which the term's tensors are sorted by name.
        Multiple lists occur when there are multiple freely commuting tensors with the same name.
        """

        # Separate the tensors into two lists:    freely commuting and non-commuting
        fcList = []
        ncList = []
        for t in self.tensors:
            if t.freelyCommutes:
                fcList.append(t.copy())
            else:
                ncList.append(t.copy())

        if len(fcList) > 0:

            # Sort the freely commuting tensors by name
            #fcList.sort(lambda x,y: cmp(x.name,y.name))
            fcList.sort(key=lambda x: x.name)

            # For the freely commuting tensors, compile a list of unique tensor names and the number of times they occur
            uniqueNames = []
            nameCounts = []
            for ten in fcList:
                if ten.name in uniqueNames:
                    nameCounts[uniqueNames.index(ten.name)] += 1
                else:
                    uniqueNames.append(ten.name)
                    nameCounts.append(1)

            # For each unique name, generate ordering permutations
            permutes = []
            for i in nameCounts:
                permutes.append(makePermutations(i))

            # Combine the name perturbations above into all possible lists in which the freely commuting tensors
            # are ordered by name
            tensorsByName = []
            total = 0
            for i in nameCounts:
                tensorsByName.append(fcList[total:total+i])
                total += i
            candidateLists = []
            for perm in permutes[0]:
                candidateLists.append([])
                for p in perm:
                    candidateLists[-1].append(tensorsByName[0][p].copy())
            for i in range(1,len(uniqueNames)):
                nOldVariants = len(candidateLists)
                for perm in permutes[i]:
                    for j in range(nOldVariants):
                        candidateLists.append([])
                        candidateLists[-1].extend(candidateLists[j])
                        for p in perm:
                            candidateLists[-1].append(tensorsByName[i][p].copy())
                del(candidateLists[0:nOldVariants])

        else:
            candidateLists = [[]]

        # Append the non-commuting tensors to each candidate list
        for l in candidateLists:
            for t in ncList:
                l.append(t.copy())

        # Return the list of candidate tensor orders
        return candidateLists

    #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

def process_chunk(terms_chunk):
    for term in terms_chunk:
        term.makeCanonical(rename_user_defined=False)
    return terms_chunk
 
def combineTerms(termList, maxProcesses = None):
    "Combines any like terms in termList"

    if not termList:
        return

    if maxProcesses is None:
        maxProcesses = cpu_count()
    else:
        maxProcesses = max(1, maxProcesses)

    if options.verbose:
        print('')
        print('Combining like terms:')
        print('Converting %i terms to canonical form...' %(len(termList)))
        print('Using max threads %i' %(maxProcesses))

    startTime = time.time()

    # Put the terms in termList into their canonical (unique) forms
    if maxProcesses > 1 and len(termList) > 100:
        # Process in chunks to reduce serialization overhead
        chunk_size = max(1, len(termList) // (maxProcesses * 4))
        chunk_size = min(chunk_size, len(termList))
        chunks = [termList[i:i+chunk_size] for i in range(0, len(termList), chunk_size)]

        # Process in parallel
        with Pool(processes=maxProcesses, maxtasksperchild=1) as pool:
        #with Pool(processes=maxProcesses) as pool:
            processed_chunks = pool.map(process_chunk, chunks)

        # Flatten results
        termList[:] = [term for chunk in processed_chunks for term in chunk]

    else:
        # Convert the terms to canonical form in the main thread
        for i in range(len(termList)):
            if options.verbose:
                print('%6i    %s' %(i,str(termList[i])))
            termList[i].makeCanonical(rename_user_defined = False)

    # Sort the terms
    termList.sort()

    # Combine any terms with the same canonical form
    newTermList = []
    current = termList[0]
    for i in range(1, len(termList)):
        if (current.constants == termList[i].constants) and (current.tensors == termList[i].tensors):
            current.numConstant += termList[i].numConstant
        else:
            newTermList.append(current)
            current = termList[i]

    newTermList.append(current)
    termList[:] = newTermList

    # Rename user defined dummy indices
    for _term in termList:
        for _tensor in _term.tensors:
            for _ind in range(len(_tensor.indices)):
                _tensor.indices[_ind].rename()

    # Remove terms with coefficients of zero
    termChop(termList)

    if options.verbose:
        print("Finished combining terms in %.3f seconds" %(time.time() - startTime))
        print("")

    # Sort the terms
    termList.sort()


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def multiplyTerms(t1,t2):
    if (not isinstance(t1,term)) or (not isinstance(t2,term)):
        raise TypeError("t1 and t2 must be of type term")
    return term(t1.numConstant*t2.numConstant, t1.constants+t2.constants, t1.tensors+t2.tensors)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def termChop(termList, tolerance = 1e-6):
    "Removes any terms with zero constant factors from termList."
    TypeErrorMessage = "termList must be a list of terms"
    if not isinstance(termList, list):
        raise TypeError(TypeErrorMessage)

    if not all(isinstance(t, term) for t in termList):
        raise TypeError(TypeErrorMessage)

    termList[:] = [t for t in termList if abs(t.numConstant) >= tolerance]

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def getcim(tenList, alphabet, tenCount = 0, alphaCount = 0, inputMaps = {}):
    """
    Determines the index mapping necessary to canonicalize the indices in tenList.
    Assumes any creation/destruction operators are in normal order.
    """

    # Determine what types of tensors are left to process
    no_tensors_left    = ( tenCount == len(tenList) )
    only_sfExOp_left = ( (tenCount == len(tenList)-1) and isinstance(tenList[-1], sfExOp) )
    only_creDes_left = ( tenCount < len(tenList) )
    for t in tenList[tenCount:]:
        if ( not isinstance(t, creOp) ) and ( not isinstance(t, desOp) ):
            only_creDes_left = False
            break

    if no_tensors_left or only_creDes_left or only_sfExOp_left:

        # Copy the input mapping into a new map that can be added to
        map = {}
        map.update(inputMaps)

        # For the preceeding tensors, create an ordered list of indeces
        # and a list of tensors with the canonical index sorting
        indexList = []
        newTensorList = []
        sign = 1
        for t in tenList[:tenCount]:
            tcopy = t.copy()
            for j in range(len(tcopy.indices)):
                if tcopy.indices[j].tup() in map.keys():
                    tcopy.indices[j] = map[tcopy.indices[j].tup()]
            # Keep track of the sign produced by sorting the tensor's indices
            sign *= tcopy.sortIndices()
            newTensorList.append(tcopy)
            for ind in tcopy.indices:
                indexList.append(ind)

    if only_creDes_left:

        # Apply the input mapping to the creation/destruction operators
        # Create and apply a new mapping for any new dummy indices
        opList = []
        for op in tenList[tenCount:]:
            opList.append(op.copy())
        nNewMaps = 0
        for op in opList:
            if op.indices[0].tup() in map.keys():
                # Apply input mapping
                op.indices[0] = map[op.indices[0].tup()].copy()
            elif op.indices[0].isSummed:
                # Create new mapping
                map[op.indices[0].tup()] = index(alphabet[alphaCount+nNewMaps], op.indices[0].indType, op.indices[0].isSummed, op.indices[0].userDefined)
                nNewMaps += 1
                op.indices[0] = map[op.indices[0].tup()].copy()

        # Sort the operators and apply the resulting sign
        (s,opList) = sortOps(opList)
        sign *= s

        # Add the operators' indices to the ordered list of indices.
        # Also add the sorted operators to the new tensor list.
        for op in opList:
            indexList.append(op.indices[0])
            newTensorList.append(op)

    if only_sfExOp_left:

        # Apply the input mapping to the sfExOp
        # Create and apply a new mapping for any new dummy indices
        t = tenList[-1].copy()
        nNewMaps = 0
        for i in range(len(t.indices)):
            if t.indices[i].tup() in map.keys():
                # Apply input mapping
                t.indices[i] = map[t.indices[i].tup()].copy()
            elif t.indices[i].isSummed:
                # Create new mapping
                map[t.indices[i].tup()] = index(alphabet[alphaCount+nNewMaps], t.indices[i].indType, t.indices[i].isSummed, t.indices[i].userDefined)
                nNewMaps += 1
                t.indices[i] = map[t.indices[i].tup()].copy()

        # Sort the indices of the sfExOp (go go gadget bubble sort!)
        i = 0
        while i < t.order-1:
            if t.indices[i] > t.indices[i+1]:
                temp = t.indices[i+1]
                t.indices[i+1] = t.indices[i]
                t.indices[i] = temp
                temp = t.indices[i+t.order+1]
                t.indices[i+t.order+1] = t.indices[i+t.order]
                t.indices[i+t.order] = temp
                i = 0
            else:
                i += 1

        # Add the sfExOp to the new tensor list
        newTensorList.append(t)

        # Add the sfExOp's indices to the ordered index list
        for i in t.indices:
            indexList.append(i)
                

    if no_tensors_left or only_creDes_left or only_sfExOp_left:

        # Compute a score based on how alphabetical the ordered list of indices is
        score = []
        for i in range(len(indexList)-1):
            score.append(0)
            for j in range(i+1,len(indexList)):
                if indexList[i] < indexList[j]:
                    score[-1] += 1

        # Return the score, the mapping, the resulting sign, and the canonical tensor list produced by the mapping
        return (score,map,sign,newTensorList)

    # Otherwise, process the next tensor
    else:

        # Get the tensor's symmetry permutations
        (symPerms,factors) = tenList[tenCount].symPermutes()

        # Make a copy of the tensor to be processed
        t = tenList[tenCount].copy()

        # Apply all index maps from previous tensors to the current tensor
        for i in range(len(t.indices)):
            if t.indices[i].tup() in inputMaps:
                t.indices[i] = inputMaps[t.indices[i].tup()]

        # Determine the permutation that maximizes alphabetical order of the mapped indices
        bestPermScore = -1
        for perm in symPerms:
            permScore = 0
            indList = []
            for j in range(len(t.indices)):
                indList.append(t.indices[perm[j]])
            for i in range(len(indList)-1):
                for j in range(i+1,len(indList)):
                    if indList[j] in inputMaps.values() and indList[i] in inputMaps.values() and indList[i] < indList[j]:
                        permScore += 1
            if permScore > bestPermScore:
                bestPermScore = permScore
                bestInitialPerm = perm

        # Does it matter if there are more than one perm with max score?
        # I think it doesn't, one can select any of them

        # Reset the indices and sort them according to bestInitialPerm
        t = tenList[tenCount].copy()
        tcopy = t.copy()
        for i in range(len(t.indices)):
            tcopy.indices[bestInitialPerm[i]] = t.indices[i]
        t = tcopy

        # make a mapping of the next alphabet elements, in order, to the unassigned indices
        # also make any symmetry equivalent mappings
        bestScore = [-1]
        for k in range(len(symPerms)):
            nNewMaps = 0
            map = {}
            map.update(inputMaps)
            for i in range(len(t.indices)):
                p = symPerms[k][i]
                if t.indices[p].isSummed and ( not (t.indices[p].tup() in map) ):
                    map[t.indices[p].tup()] = index(alphabet[alphaCount+nNewMaps], t.indices[p].indType, t.indices[p].isSummed, t.indices[p].userDefined)
                    nNewMaps += 1
            (score,map,sign,newTensorList) = getcim(tenList, alphabet, tenCount+1, alphaCount+nNewMaps, map)
            if score > bestScore: #is it possible to have multiple top scores here? I think so.    Does it matter?
                bestScore = score
                bestMaps = map
                bestSign = sign
                bestNewTensorList = newTensorList

        # return the mapping that gives the best overall score
        return (bestScore,bestMaps,bestSign,bestNewTensorList)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def sortOps(unsortedOps, returnPermutation = False):
    """
    Sorts a list of creation/destruction operators into normal order and alphabetically, without performing contractions.
    Returns the overall sign resulting from the sort and the sorted operator list. Optionally also returns the permutation.
    """
    sortedOps = list(unsortedOps)
    n_ops = len(unsortedOps)

    i = 0
    sign = 1

    if returnPermutation:
        perm = list(range(n_ops))

    while i < n_ops-1:
        if sortedOps[i] <= sortedOps[i+1]:
            i += 1
        else:
            sortedOps[i], sortedOps[i+1] = sortedOps[i+1], sortedOps[i]
            if returnPermutation:
                perm[i], perm[i+1] = perm[i+1], perm[i]
            i = 0
            sign *= -1

    if returnPermutation:
        return sign, sortedOps, perm
    return sign, sortedOps


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def removeCoreOpPairs(inList):
    """
    Removes pairs of core creation and core destruction operators corresponding to the same core index.
    Does not remove a pair if it's creation or destruction operator is repeated.
    Input is a list of terms.
    The terms must be in normal order.
    """

    # prepare input argument
    if not isinstance(inList, list):
        raise TypeError("input must be a list of terms")

    # loop over input terms
    for t in inList:

        # Check that the term is indeed a term
        if not isinstance(t, term):
            raise TypeError("input must be a term or list of terms")

        # Check that the term is normal ordered
        if not t.isNormalOrdered():
            raise ValueError("core index removal function only works for normal ordered terms")

        # Initialize a counter for unremoved creation operators
        creCount = 0

        # Loop over the term's tensors
        i = 0
        while i < len(t.tensors):

            # Initialize flags
            operatorsRemoved = False
            repeatedCreOp = False
            repeatedDesOp = False

            # if the tensor is a core creation operator
            if isinstance(t.tensors[i], creOp) and options.core_type in t.tensors[i].indices[0].indType:

                # Check whether the creation operator is a repeat of an earlier creation operator
                for k in range(i):
                    if t.tensors[k] == t.tensors[i]:
                        repeatedCreOp = True

                # If a repeat, move to the next tensor
                if not repeatedCreOp:

                    # otherwise...
                    
                    # create the matching destruction operator
                    matchingDesOp = desOp(t.tensors[i].indices[0])

                    # search for the matching destruction operator
                    for j in range(i+1,len(t.tensors)):

                        # if a repeat of the creation operator is found, move to the next tensor
                        if t.tensors[j] == t.tensors[i]:
                            break

                        # if the matching destruction operator is found
                        if t.tensors[j] == matchingDesOp:

                            # initialize a counter for the number of operator commutations necessary to move the
                            # matching creation and destruction operators to the begining and end of the term,
                            # respectively
                            commCount = creCount

                            # count the number of destruction operators after the matching destruction operator
                            for k in range(j+1, len(t.tensors)):
                                if isinstance(t.tensors[k], desOp):
                                    commCount += 1

                                # while counting, check whether a repeat of the matching destruction operator is present
                                if t.tensors[k] == t.tensors[j]:
                                    repeatedDesOp = True

                            # if there is a repeat of the matching destruction operator, move to the next tensor
                            if repeatedDesOp:
                                break

                            # scale the term by the factor resulting from commuting the creation operator and
                            # matching destruction operator to the beginning and end of the term, respectively
                            t.scale((-1)**commCount)

                            # delete the creation and matching destruction operators
                            del t.tensors[j]
                            del t.tensors[i]

                            # set the operator removal flag to true
                            operatorsRemoved = True

                            # move to the next tensor
                            break

            # If no operators were removed...
            if not operatorsRemoved:

                # If the current tensor is a creation operator, increase the creation operator count
                if isinstance(t.tensors[i], creOp):
                    creCount += 1

                # increment the index
                i += 1


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def removeCoreOps_sf(inList):
    """
    Removes core indices from inList's terms' spin-free excitation operators.
    This function assumes that the spin-free operators will be converted to density matrices
    by taking their expectation value immediately after this function is finished.
    Terms that have a zero expectation value due to the nature of their spin-free operator's core
    indices are deleted from inList.
    """

    if options.verbose:
        print("removing core creation and destruction operators in preperation for conversion to RDMs by an expectation value...")
        print("")

    # loop repeatedly through the terms until no core indices are left
    hasCore = True
    while hasCore:

        hasCore = False

        # process each term
        t_num = 0
        while t_num < len(inList):

            # use a short name for the current term
            t = inList[t_num]

            # check that t is a term
            if not isinstance(t, term):
                raise TypeError("inList must be a list of term objects")

            # check for normal ordering
            if not t.isNormalOrdered():
                raise ValueError("input terms must be normal ordered")

            # check for spin-orbital creation and destruction operators
            for ten in t.tensors:
                if isinstance(ten, creOp) or isinstance(ten, desOp):
                    raise TypeError("input terms may not contain creOp or desOp objects")

            # find the spin-free excitation operator
            op = None
            for i in range(len(t.tensors)):
                if isinstance(t.tensors[i], sfExOp):
                    op = t.tensors[i]
                    opPos = i

            # if there is no spin-free excitation operator, skip the term
            if op is None:
                t_num += 1
                continue

            # ensure the sfExOp has no indices with multiple type groups or a type group with core and non-core types
            for ind in op.indices:
                if len(ind.indType) > 1:
                    raise ValueError("index %s in term (%s) has more than one type group:    %s" %(ind.name, str(t), str(typeGroup)))
                for typeGroup in ind.indType:
                    if options.core_type[0] in typeGroup and len(typeGroup) > 1:
                        raise ValueError("index %s in term (%s) has a type group including core and non-core types:    %s" %(ind.name, str(t), str(typeGroup)))

            # compute the operator's order
            order = len(op.indices)/2

            # find a core index
            cInd = None
            for i in range(2*order):
                if op.indices[i].indType == (options.core_type,):
                    cInd = op.indices[i]
                    break

            # if there are no core indices, move to the next term
            if cInd is None:
                t_num += 1
                continue

            # if there is a core index, request another loop through the terms because all core indices have not been found
            hasCore = True

            # count the number of times the targeted core index appears among creation and destruction operators
            nCre = 0
            nDes = 0
            for i in range(order):
                if op.indices[i] == cInd:
                    nCre += 1
                if op.indices[order+i] == cInd:
                    nDes += 1

            # if the term is equal to zero, remove it and move to the next term
            if nCre != nDes or nCre > 2 or nDes > 2:
                del inList[t_num]
                continue

            # organize the operator's indices into vertical pairs of cre/des operator indices
            pairs = [ [op.indices[i],op.indices[order+i]] for i in range(order)]

            # print out the initial term
            if options.verbose:
                print("    initial term: ", t)
#                print "verticle pairs: ",
#                for p in pairs:
#                    print " [%s,%s]" %(p[0].name, p[1].name),
#                print ""

            # make sure the number of pairs is equal to the operator's order
            if len(pairs) != order:
                raise ValueError("number of pairs not equal to operator's order")

            # determine the new operator's indices.
            # record how many pairs there were with both elements equal to the targeted core operator
            nMatch = 0
            topUnmatched = []
            botUnmatched = []
            i = order-1
            while i >= 0:
                if pairs[i][0] == cInd and pairs[i][1] == cInd:
                    del pairs[i]
                    nMatch += 1
                elif pairs[i][0] == cInd:
                    botUnmatched.append(pairs.pop(i)[1])
                elif pairs[i][1] == cInd:
                    topUnmatched.append(pairs.pop(i)[0])
                i -= 1
            newIndices = []
            newIndices.extend(topUnmatched)
            for p in pairs:
                newIndices.append(p[0])
            newIndices.extend(botUnmatched)
            for p in pairs:
                newIndices.append(p[1])

            # replace the old operator with the new operator in which the targeted core index has been removed
            if len(newIndices) > 0:
                t.tensors[opPos] = sfExOp(newIndices)
            # if there are no indices left after the core index's removal, remove the old operator
            else:
                del t.tensors[opPos]

            # apply the appropriate constant factor
            if     nCre == 1 and nMatch == 1:
                t.scale(2.0)
            elif nCre == 1 and nMatch == 0:
                t.scale(-1.0)
            elif nCre == 2 and nMatch == 2:
                t.scale(2.0)
            elif nCre == 2 and nMatch == 1:
                t.scale(-1.0)
            elif nCre == 2 and nMatch == 0:
                t.scale(1.0)
            else:
                raise ValueError("unexpected values:    nCre = %i, nMatch = %i" %(nCre, nMatch))

            # print out the final term
            if options.verbose:
#                print "    topUnmatched: ",
#                for p in topUnmatched:
#                    print " %s" %(p.name),
#                print ""
#                print "    botUnmatched: ",
#                for p in botUnmatched:
#                    print " %s" %(p.name),
#                print ""
                print("        final term: ", t)

            # for the special case of two unmatched pairs, the result is a sum of two different operators.
            # the first replaced the original operator, and the second is added here.
            if nCre == 2 and nMatch == 0:
                inList.append(t.copy())
                if len(newIndices) < 4:
                    raise ValueError("expected at least 4 remaining indices for nCre == 2 and nMatch == 0 case, but only %i are present" %len(newIndices))
                (newIndices[0], newIndices[1]) = (newIndices[1], newIndices[0])
                inList[-1].tensors[opPos] = sfExOp(newIndices)
                # print out the additional final term
                if options.verbose:
                    print("2nd final term: ", inList[-1])

            # print a blank line
            if options.verbose:
                print("")

            # increment the index to the next term
            t_num += 1


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def removeVirtOps_sf(inList):
    """
    Removes from inList any terms containing a spin-free operator with a virtual index.
    """

    if options.verbose:
        print("removing terms containing a spin-free operator with a virtual index...")
        print("")

    # loop over the terms in inList
    i = 0
    while i < len(inList):

        # ensure that each element of inList is a term object
        if not isinstance(inList[i], term):
            raise TypeError("inList must be a list of term objects")

        # determine if the term's spin-free excitation operators have any virtual indices
        hasVirtual = False
        for ten in inList[i].tensors:
            if isinstance(ten, sfExOp):
                for ind in ten.indices:
                    if options.virtual_type in ind.indType:
                        hasVirtual = True

        # remove the term if a spin-free excitation operator had a virtual index
        if hasVirtual:
            if options.verbose:
                print(" removing term: ", inList[i])
            del inList[i]

        # otherwise, move to the next term
        else:
            i += 1

    if options.verbose:
        print("")


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
