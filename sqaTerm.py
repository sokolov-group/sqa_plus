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
from collections import deque, defaultdict
from itertools import islice, count

from .sqaIndex import index, ind_type_order
from .sqaTensor import tensor, kroneckerDelta, sfExOp, creOp, desOp
from .sqaMisc import makePermutations
from .sqaOptions import options
import time
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


@total_ordering
class term:
    "A class for terms used in operator algebra. Each term is a multiplicative string of constants, tensors, and operators."
    TOL = 1e-6
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

    def _comparison_key(self):
        """
        Return tuple for comparison:
        (nr creOps, nr desOps, sfExOp ranks, nr constants, nr tensors, tensor names, tensors, constants)
        """
        return (
            self.nCreOps(),
            self.nDesOps(),
            self.sfExOp_ranks(),
            len(self.constants),
            len(self.tensors),
            tuple(t.name for t in self.tensors),
            tuple(self.tensors),
            tuple(self.constants),
        )

    def __eq__(self, other):
        if not isinstance(other, term):
            raise TypeError("term object can only be compared to other term objects.")

        # compare keys and numerical constants
        return (
            self._comparison_key() == other._comparison_key() and
            abs(self.numConstant - other.numConstant) < term.TOL
        )

    def __lt__(self, other):
        if not isinstance(other, term):
            raise TypeError("term object can only be compared to other term objects.")

        if self._comparison_key() != other._comparison_key():
            return self._comparison_key() < other._comparison_key()

        # if keys equal, compare numerical constants
        if abs(self.numConstant - other.numConstant) < term.TOL:
            return False
        return self.numConstant < other.numConstant

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
        #retval = [len(t.indices)/2 for t in self.tensors if isinstance(t, sfExOp)]
        retval = [t.order for t in self.tensors if isinstance(t, sfExOp)]
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

            # Skip non-kronecker delta functions
            if not isinstance(t, kroneckerDelta):
                i += 1
                continue

            # Kronecker delta indices
            i0, i1 = t.indices[0], t.indices[1]

            # Remove delta functions with a repeated index
            if (i0 == i1) and (i0.userDefined == i1.userDefined):
                del self.tensors[i]
                continue

            # Skip delta functions with no contractible indices
            if not (i0.isSummed or i1.isSummed):
                i += 1
                continue

            # Perform contraction for delta functions with contractible indices
            # Validate that indices are compatible
            if len(i0.indType) != len(i1.indType):
                raise RuntimeError(f"Cannot contract indices {i0}, {i1}. They have different numbers of type groups.")

            # Determine new index type based on type overlap
            typeOverlap = [
                [typeString for typeString in i0.indType[j] if typeString in i1.indType[j]]
                for j in range(len(i0.indType))
            ]

            # If no overlap, then delta function is zero
            if (len(i0.indType) > 0 or len(i1.indType) > 0) and [] in typeOverlap:
                self.numConstant = 0.0
                return

            # Create the new index
            if not i0.isSummed:
                newIndex = index(i0.name, typeOverlap, i0.isSummed, i0.userDefined)
            else:
                newIndex = index(i1.name, typeOverlap, i1.isSummed, i1.userDefined)

            # Remove the delta function
            del self.tensors[i]

            # Replace old indices with new index
            for ten in self.tensors:
                ten.indices = [
                    newIndex.copy() if old_idx in (i0, i1) else old_idx
                    for old_idx in ten.indices
                ]

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

    #------------------------------------------------------------------------------------------------

    def makeCanonical_non_recursive(self, rename_user_defined = True):
        """Converts the term to a unique canonical form using a non-recursive algorithm."""

        # Early exit checks
        if self.isInCanonicalForm:
            return

        if not self.tensors:
            self.isInCanonicalForm = True
            return

        if not self.isNormalOrdered():
            raise RuntimeError("A term must have normal ordered operators to be converted to canonical form.")

        # Sort constants
        self.constants.sort()

        # Partition tensors based on commutation properties
        fc_list = [t for t in self.tensors if t.freelyCommutes]
        nc_list = [t for t in self.tensors if not t.freelyCommutes]

        # Sort freely commuting tensors
        fc_groups = defaultdict(list)
        for t in fc_list:
            fc_groups[t.name].append(t)

        # external group sort
        def group_key(g):
            return g[0].__class__.__name__, len(g), g[0].name

        # internal group sort
        def element_key(t):
            score = [ind_type_order(ind.indType) for ind in t.indices]
            return len(t.indices), tuple(sorted(score)), tuple(score)

        fc_groups = [sorted(g, key=element_key) for g in fc_groups.values()]
        fc_groups = sorted(fc_groups, key=group_key)

        # Add non-commuting tensors as individual groups
        name_groups = fc_groups + [[t] for t in nc_list]

        # Generate alphabet for renaming indices
        alphabet = self.generateAlphabet()

        # Initialize best state tracking
        best_state = {
            'tensor_list': None,
            'factor': None,
            'map': None,
            'score': None,
            'index_list': None,
            'count': 0,
        }

        job_stack = deque([({}, 0, 0, 0, [])])

        # Determine the best ordering and index mapping
        while job_stack:
            current_map, g_count, t_count, a_count, g_perms = job_stack.pop()

            # ===== Terminal Case: All Groups Processed =====
            if g_count == len(name_groups):
                ten_list, index_list, factor = self._get_tensor_list(name_groups, g_perms, current_map)
                self._update_best_state(best_state, ten_list, index_list, factor, current_map)
                continue

            # ===== Case: Only creOp/desOp Tensors Left =====
            if all(len(group) == 1 and isinstance(group[0], (creOp, desOp)) for group in name_groups[g_count:]):
                ten_list, index_list, factor = self._get_tensor_list(name_groups[:g_count], g_perms, current_map)

                # Process ops
                op_list = [name_groups[i][0].copy() for i in range(g_count, len(name_groups))]
                new_maps_count = 0

                for op in op_list:
                    if op.indices[0].tup() in current_map:
                        op.indices[0] = current_map[op.indices[0].tup()]
                    elif op.indices[0].isSummed:
                        current_map[op.indices[0].tup()] = index(
                            alphabet[a_count + new_maps_count],
                            op.indices[0].indType,
                            op.indices[0].isSummed,
                            op.indices[0].userDefined
                        )
                        op.indices[0] = current_map[op.indices[0].tup()]
                        new_maps_count += 1

                sign, op_list = sortOps(op_list)
                factor *= sign

                index_list.extend(op.indices[0] for op in op_list)
                ten_list.extend(op_list)

                self._update_best_state(best_state, ten_list, index_list, factor, current_map)
                continue

            # ===== Case: Single sfExOp Tensor Left =====
            if (g_count == len(name_groups) - 1 and len(name_groups[g_count]) == 1 and
                isinstance(name_groups[g_count][0], sfExOp)):

                ten_list, index_list, factor = self._get_tensor_list(name_groups[:g_count], g_perms, current_map)

                t = name_groups[g_count][0].copy()
                new_maps_count = 0

                for i in range(len(t.indices)):
                    if t.indices[i].tup() in current_map:
                        t.indices[i] = current_map[t.indices[i].tup()]
                    elif t.indices[i].isSummed:
                        current_map[t.indices[i].tup()] = index(
                            alphabet[a_count + new_maps_count],
                            t.indices[i].indType,
                            t.indices[i].isSummed,
                            t.indices[i].userDefined
                        )
                        t.indices[i] = current_map[t.indices[i].tup()]
                        new_maps_count += 1

                # Sort sfExOp indices using bubble sort
                i = 0
                while i < t.order-1:
                    if t.indices[i] > t.indices[i+1]:
                        t.indices[i], t.indices[i+1] = t.indices[i+1], t.indices[i]
                        j = i + t.order
                        t.indices[j], t.indices[j+1] = t.indices[j+1], t.indices[j]
                        i = 0
                    else:
                        i += 1

                ten_list.append(t)
                index_list.extend(t.indices)

                self._update_best_state(best_state, ten_list, index_list, factor, current_map)
                continue

            # ===== Start New Name Group =====
            if len(g_perms) <= g_count:
                with_mapped = []
                without_mapped = []

                for i, tensor in enumerate(name_groups[g_count]):
                    least_mapped = False
                    for ind in tensor.indices:
                        if ind.tup() in current_map:
                            mapped_val = current_map[ind.tup()]
                            if least_mapped is False or mapped_val < least_mapped:
                                least_mapped = mapped_val

                    if least_mapped is False:
                        without_mapped.append(i)
                    else:
                        with_mapped.append((least_mapped, i))

                with_mapped.sort(key=lambda x: x[0])

                if len(without_mapped) <= 1:
                    sorted_indices = [i[1] for i in with_mapped] + without_mapped
                    job_stack.append((current_map, g_count, t_count, a_count, g_perms + [sorted_indices]))
                else:
                    for perm in makePermutations(len(without_mapped)):
                        new_perm = [i[1] for i in with_mapped] + [without_mapped[j] for j in perm]
                        job_stack.append((current_map, g_count, t_count, a_count, g_perms + [new_perm]))

                continue

            # ===== Continue Existing Group =====
            t = name_groups[g_count][g_perms[g_count][t_count]]
            sym_perms, _ = t.symPermutes()

            next_g_count = g_count
            next_t_count = t_count + 1
            if next_t_count == len(name_groups[g_count]):
                next_g_count += 1
                next_t_count = 0

            # Schedule new jobs for each symmetry-equivalent index ordering of current tensor, t
            for perm in sym_perms:
                new_map = dict(current_map)
                new_maps_count = 0

                for ind in [t.indices[perm[i]] for i in range(len(t.indices))]:
                    if ind.isSummed and ind.tup() not in new_map:
                        new_map[ind.tup()] = index(
                            alphabet[a_count + new_maps_count],
                            ind.indType,
                            ind.isSummed,
                            ind.userDefined
                        )
                        new_maps_count += 1

                job_stack.append((new_map, next_g_count, next_t_count, a_count + new_maps_count, g_perms))

        # Finalize canonical alphabet mapping
        best_map = best_state['map'] or {}
        alphabet_list = list('abcdefghijklmnopqrstuvwxyz')
        filtered_alphabet = [c for c in alphabet_list if c not in options.user_defined_indices]

        if len(filtered_alphabet) < len(best_map):
            raise RuntimeError("Alphabet smaller than number of indices, no more names left!")

        canon_map = {}
        for _, val in sorted(best_map.items(), key=lambda kv: kv[1].name):
            canon_map[val.tup()] = val.copy()
            canon_map[val.tup()].name = filtered_alphabet.pop(0)

        # Apply canonical mapping to final tensor list
        for t in (best_state['tensor_list'] or []):
            for i in range(len(t.indices)):
                if t.indices[i].tup() in canon_map:
                    t.indices[i] = canon_map[t.indices[i].tup()].copy()
                if rename_user_defined and t.indices[i].userDefined:
                    t.indices[i].rename()

        #if best_state['count'] > 1:
        #    print(f'Multiple states are best (count = {best_state['count']})...')

        # Finalize results
        self.tensors = best_state['tensor_list'] or []
        if best_state['factor'] is not None:
            self.scale(best_state['factor'])
        self.isInCanonicalForm = True

    def _get_tensor_list(self, name_groups, g_perms, current_map):
        """Compute a new tensor list with renamed dummy indices, a list of the ordered dummy indices,
        and the multiplicative factor generated by sorting the indices."""
        ten_list = []
        index_list = []
        factor = 1
        # Build tensor and index lists using group permutations
        for group, perm in zip(name_groups, g_perms):
            for t_ind in perm:
                t_src = group[t_ind].copy()
                # Rename summed (dummy) indices using current_map
                t_src.indices = [
                    current_map[idx.tup()] if idx.isSummed else idx
                    for idx in t_src.indices
                ]
                # Sort indices in-place and accumulate resulting factor
                factor *= t_src.sortIndices()
                index_list.extend(t_src.indices)
                ten_list.append(t_src)
        return ten_list, index_list, factor

    def _update_best_state(self, best_state, ten_list, index_list, factor, current_map):
        """Update best state if current candidate is better."""
        # Compute a score for index_list based on how alphabetical the indices are
        score = []
        for i in range(len(index_list) - 1):
            count = sum(
                1 for j in range(i + 1, len(index_list))
                if index_list[i] < index_list[j]
            )
            score.append(count)

        # Update state tracker 
        if best_state['score'] is None or score > best_state['score']:
            best_state.update({
                'tensor_list': ten_list,
                'factor': factor,
                'map': current_map,
                'score': score,
                'index_list': index_list,
                'count': 1
            })
        elif score == best_state['score']:
            best_state['count'] += 1

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

def process_chunk(terms_chunk):
    for t in terms_chunk:
        t.makeCanonical(rename_user_defined = False)
    return terms_chunk

def combineTerms(term_list, max_processes = None):
    """Combines like terms in list of terms."""

    if not term_list:
        return

    if max_processes is None:
        max_processes = cpu_count()
    else:
        max_processes = max(1, max_processes)

    if options.verbose:
        print('\nCombining like terms:')
        print('Converting %i terms to canonical form...' %(len(term_list)))
        print('Using max threads %i' %(max_processes))

    start_time = time.time()

    # Canonicalize terms
    n_terms = len(term_list)
    if max_processes > 1 and n_terms > 100:
        # Process chunks in parallel to reduce serialization overhead
        chunk_size = max(1, n_terms // (max_processes * 4))
        chunks = [term_list[i:i+chunk_size] for i in range(0, n_terms, chunk_size)]

        with Pool(processes=max_processes, maxtasksperchild=1) as pool:
            processed_chunks = pool.map(process_chunk, chunks)

        # Flatten results
        term_list[:] = [t for chunk in processed_chunks for t in chunk]

    else:
        # Convert terms in serial
        for i, t in enumerate(term_list, start = 1):
            if options.verbose:
                print('%6i    %s' % (i, t))
            t.makeCanonical(rename_user_defined = False)

    # Sort the terms
    term_list.sort()

    # Combine any terms with the same canonical form
    new_term_list = []
    for t in term_list:
        if (new_term_list and
            new_term_list[-1].constants == t.constants and
            new_term_list[-1].tensors == t.tensors
        ):
            new_term_list[-1].numConstant += t.numConstant
        else:
            new_term_list.append(t)
    term_list[:] = new_term_list

    # Rename user defined dummy indices
    for _term in term_list:
        for _tensor in _term.tensors:
            for idx in _tensor.indices:
                idx.rename()

    # Remove terms with coefficients of zero
    termChop(term_list)

    if options.verbose:
        print("Finished combining terms in %.3f seconds\n" %(time.time() - start_time))

    # Sort the terms
    term_list.sort()


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def multiplyTerms(t1,t2):
    if (not isinstance(t1,term)) or (not isinstance(t2,term)):
        raise TypeError("t1 and t2 must be of type term")
    return term(t1.numConstant*t2.numConstant, t1.constants+t2.constants, t1.tensors+t2.tensors)

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def termChop(termList):
    "Removes any terms with zero constant factors from termList."
    TypeErrorMessage = "termList must be a list of terms"
    if not isinstance(termList, list):
        raise TypeError(TypeErrorMessage)

    if not all(isinstance(t, term) for t in termList):
        raise TypeError(TypeErrorMessage)

    termList[:] = [t for t in termList if abs(t.numConstant) >= term.TOL]

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

    from .sqaIndex import is_virtual_index_type

    if options.verbose:
        print("removing terms containing a spin-free operator with a virtual index...\n")

    if not all(isinstance(t, term) for t in inList):
        raise TypeError("inList must be a list of term objects")

    # Filter out sfExOp terms with virtual indices and log removals
    terms_to_remove = [
        t for t in inList
        if any(
            is_virtual_index_type(ind)
            for ten in t.tensors if isinstance(ten, sfExOp)
            for ind in ten.indices)
    ]

    if options.verbose:
        for t in terms_to_remove:
            print(" removing term: ", t)

    inList[:] = [t for t in inList if t not in terms_to_remove]

    # # loop over the terms in inList
    # i = 0
    # while i < len(inList):

    #     # ensure that each element of inList is a term object
    #     if not isinstance(inList[i], term):
    #         raise TypeError("inList must be a list of term objects")

    #     # determine if the term's spin-free excitation operators have any virtual indices
    #     hasVirtual = False
    #     for ten in inList[i].tensors:
    #         if isinstance(ten, sfExOp):
    #             for ind in ten.indices:
    #                 if options.virtual_type in ind.indType:
    #                     hasVirtual = True

    #     # remove the term if a spin-free excitation operator had a virtual index
    #     if hasVirtual:
    #         if options.verbose:
    #             print(" removing term: ", inList[i])
    #         del inList[i]

    #     # otherwise, move to the next term
    #     else:
    #         i += 1

    # if options.verbose:
    #     print("")


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
