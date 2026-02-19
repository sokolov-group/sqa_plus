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
# The tensor class represents an object consisting of a name, an ordered set of indices,
# and possibly a set of symmetry permutations among the indices.
#
# A tensor's name is a string.
#
# The indices are given by a list of objects of the index class.
#
# The symmetries are given by a list of objects of the symmetry class. Note that this
# list does not have to be exhaustive, as the code will attempt to find all possible
# symmetry permutations that can be created from the supplied symmetries.
#
# For example, the list of symmetries
#     [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2),-1), sqa.symmetry((1,0,3,2),1)]
# is redundant, as the third permutation can and will be generated using the first two.
#
#
# There are 4 important children of the tensor class: creOp, desOp, sfExOp, and kroneckerDelta.
#    - creOp and desOp are one-index tensors representing creation and destruction operators.
#    - sfExOp is a 2n-index tensor representing a spin-free excitation operator.
#    - kroneckerDelta is a two-index tensor representing the Kronecker delta function.
#

from functools import total_ordering
from .sqaIndex import index
from .sqaSymmetry import symmetry

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class tensor:
    """A class to represent tensors in operator algebra. Integrals and density matrices are examples."""

    #------------------------------------------------------------------------------------------------

    freelyCommutes = True

    def __init__(self, name, indices=None, symmetries=None):
        """Initialize a tensor with a name, indices, and symmetries."""

        # Initialize permutation
        self.permutations = None
        self.factors = None

        # Process name
        self.name = str(name)

        # Process indices
        if indices is None:
            indices = []

        if not isinstance(indices, list) or not all(isinstance(i, index) for i in indices):
            raise TypeError("indices must be a list of index objects")
        self.indices = [i.copy() for i in indices]

        # Process symmetries
        if symmetries is None:
            symmetries = []

        if not isinstance(symmetries, list) or not all(isinstance(sym, symmetry) for sym in symmetries):
            raise TypeError("symmetries must be a list of symmetry objects")

        self.symmetries = []
        for sym in symmetries:
            if any(s.pattern == sym.pattern for s in self.symmetries):
                raise ValueError("a tensor cannot have two symmetries with the same pattern")
            self.symmetries.append(sym.copy())

    #------------------------------------------------------------------------------------------------

    def _comparison_key(self):
        """Return tuple for comparison: (name, indices, symmetries)."""
        return (self.name, self.indices, self.symmetries)

    def __eq__(self, other):
        if not isinstance(other, tensor):
            return False

        # If other belongs to a tensor subclass, use the subclass's comparison method
        if type(other) is not type(self):
            return other == self

        return self._comparison_key() == other._comparison_key()

    def __lt__(self, other):
        if not isinstance(other, tensor):
            raise TypeError("A tensor may only be compared to another tensor")

        # If other belongs to a tensor subclass, use the subclass's comparison method
        if type(other) is not type(self):
            return self < other

        return self._comparison_key() < other._comparison_key()

    #------------------------------------------------------------------------------------------------

    def __str__(self):
        indices_str = ",".join(index.name for index in self.indices)
        return f"{self.name}({indices_str})"

    def __repr__(self):
        return str(self)

    #------------------------------------------------------------------------------------------------

    def copy(self):
        """Returns a copy of the tensor"""
        retval = tensor(self.name, self.indices, self.symmetries)
        if self.permutations != None and self.factors != None:
            retval.permutations = [ perm + [] for perm in self.permutations ]
            retval.factors = self.factors + []
        return retval

    #------------------------------------------------------------------------------------------------

    def symPermutes(self, force = False):
        """Returns the index permutations and resulting factors allowed by the tensor's symmetry."""

#        # If the result is already known, return it
#        if not force and self.permutations != None and self.factors != None:
#            return (self.permutations,self.factors)

        # Compute the permutations and corresponding factors
        permutations = [list(range(len(self.indices)))]
        factors = [1]
        known_perm = {tuple(permutations[0])}

        idx = 0
        while idx < len(permutations):
            perm = permutations[idx]
            factor = factors[idx]

            # Generate new permutations
            for sym in self.symmetries:
                new_perm = [perm[i] for i in sym.pattern]
                new_perm_tuple = tuple(new_perm)
    
                if new_perm_tuple not in known_perm:
                    known_perm.add(new_perm_tuple)
                    permutations.append(new_perm)
                    factors.append(sym.factor * factor)

            idx += 1

        # Save the results for later so they don't need to be computed again
        self.permutations, self.factors = permutations, factors
        return permutations, factors

    #------------------------------------------------------------------------------------------------

    def sortIndices(self):
        """Sort indices alphabetically within symmetry constraints. Returns the resulting symmetry factor."""

        # If the tensor has no symmetry, do nothing
        if not self.symmetries:
            return 1

        # Get allowed symmetry permutations and corresponding factors
        tuples, factors = map(list, self.symPermutes())
        scores = [0] * len(tuples)
        n_ind = len(self.indices)

        # Score the different permutations and select the winner
        for i in range(n_ind - 1):
            for j in range(i+1, n_ind):
                for p, tup in enumerate(tuples):
                    if self.indices[tup[i]] < self.indices[tup[j]]:
                        scores[p] += 1

            # Determine max score and keep only tuples with that score
            max_score = max(scores)
            keep = [p for p, s in enumerate(scores) if s == max_score]
            tuples = [tuples[p] for p in keep]
            factors = [factors[p] for p in keep]
            scores = [scores[p] for p in keep]

            if not tuples:
                break

        # Raise an error if no unique winner found
        if len(tuples) > 1:
            ref = [self.indices[p] for p in tuples[0]]

            for tup in tuples[1:]:
                compare = [self.indices[p] for p in tup]
                if any(i != j for i,j in zip(compare, ref)):
                    msg = f"No unique winner produced when sorting indices of tensor {self!s}"
                    raise RuntimeError("Scoring system did not produce unique winner.\n" + msg)

        # Set the sorted indices and return factor
        self.indices = [self.indices[p] for p in tuples[0]]
        return factors[0]

    #------------------------------------------------------------------------------------------------

    ## DEAD CODE ##
    ## def hasIndex(self,i):
    ##     """Returns True if i is one of the tensor's indices and False otherwise."""
    ##     if not isinstance(i,index):
    ##         raise TypeError("i must be of the index class")
    ##     return (i in self.indices)

    #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class kroneckerDelta(tensor):
    """A tensor representation of the kronecker delta function."""
    #------------------------------------------------------------------------------------------------

    freelyCommutes = True

    def __init__(self, indices):
        if len(indices) != 2:
            raise ValueError("The kronecker delta function takes exactly two indices")

        # Initialize attributes
        self.permutations = None
        self.factors = None
        self.name = "kdelta"
        self.indices = [i.copy() for i in indices]
        self.symmetries = [symmetry((1, 0), 1)]

    #------------------------------------------------------------------------------------------------

    # def _comparison_key(self):
    #     """kroneckerDelta comparison key - used for ordering"""
    #     return (self.name, self.indices, self.symmetries)

    def __eq__(self, other):
        if isinstance(other, kroneckerDelta):
            return self._comparison_key() == other._comparison_key()
        return False

    def __lt__(self, other):
        if isinstance(other, kroneckerDelta):
            return self._comparison_key() < other._comparison_key()
        elif isinstance(other, (creOp, desOp, creDesTensor, sfExOp)):
            return True
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("A kroneckerDelta may only be compared to another tensor")

    #------------------------------------------------------------------------------------------------

    def copy(self):
        "Returns a copy of the kroneckerDelta object"
        return kroneckerDelta(self.indices)

    #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class sfExOp(tensor):
    """
    A tensor representation of a spin free excitation operator.
    e.g. E(i,j,k,l) = sum over ( sigma, tau )    of ( a+_(i sigma) a+_(j tau) a_(l tau) a_(k sigma) )
    Note that the indices i,j,k,l refer to spacial orbitals and sigma,tau refer to spins.
    """

    #------------------------------------------------------------------------------------------------

    freelyCommutes = False

    #------------------------------------------------------------------------------------------------

    def __init__(self, indices):
        # Check that there are an even number of indices
        if len(indices) % 2 != 0:
            raise ValueError("A spin free excitation operator (the sfExOp class) must have an even number of indices")

        # Initialize attributes
        self.permutations = None
        self.factors = None
        self.order = len(indices) // 2
        self.name = f"E{self.order}"
        self.indices = [i.copy() for i in indices]
        self._init_symmetries()

    def _init_symmetries(self):
        """Generate symmetries for spin-free excitation operators."""
        self.symmetries = []
        for i in range(self.order - 1):
            pattern = [1] if i == 0 else [0]
            for j in range(1, 2 * self.order):
                if j == i:
                    pattern.append(i + 1)
                elif j == i + 1:
                    pattern.append(i)
                elif j == i + self.order:
                    pattern.append(i + 1 + self.order)
                elif j == i + 1 + self.order:
                    pattern.append(i + self.order)
                else:
                    pattern.append(j)
            self.symmetries.append(symmetry(tuple(pattern), 1))

    #------------------------------------------------------------------------------------------------

    # def _comparison_key(self):
    #     """sfExOp comparison key - used for ordering"""
    #     return (self.name, self.indices, self.symmetries)

    def __eq__(self, other):
        if isinstance(other, sfExOp):
            return self._comparison_key() == other._comparison_key()
        return False

    def __lt__(self, other):
        if isinstance(other, sfExOp):
            return self._comparison_key() < other._comparison_key()
        elif isinstance(other, (creOp, desOp, creDesTensor)):
            return True
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("An sfExOp object may only be compared to another tensor")

    #------------------------------------------------------------------------------------------------

    def copy(self):
        return sfExOp(self.indices)

    #------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class creDesTensor(tensor):

    freelyCommutes = False

    def __init__(self, ops, trans_rdm=False, symmetries=None):

        if not isinstance(ops, list) or any(not isinstance(op, (creOp, desOp)) for op in ops):
            raise TypeError("ops must be a normal ordered list of creOp and desOp objects")

        # Initialize attributes
        self.permutations = None
        self.factors = None

        self.trans_rdm = trans_rdm
        self.name = 'trdm' if trans_rdm else 'rdm'

        # Count creation/destruction operators and build index list
        self.ops = ops
        self.nCre = 0
        self.nDes = 0
        self.indices = []
        des_flag = False

        # Ensure normal-ordering
        for op in ops:
            if isinstance(op, desOp):
                self.nDes += 1
                des_flag = True
            elif des_flag:
                raise TypeError("ops must be a normal ordered list of creOp and desOp objects")
            else:
                self.nCre += 1

            self.indices.append(op.indices[0].copy())

        # Create count of total indices
        self.nInd = self.nCre + self.nDes

        # Initialize symmetries
        self._init_symmetries(symmetries)

    def _init_symmetries(self, symmetries):
        """Initialize symmetries for creation/destruction tensor."""
        if symmetries:
            self.symmetries = symmetries
        else:
            self.symmetries = []
            n = len(self.indices)

            if n > 1:
                swap_values = list(range(n - 1))
                if self.nCre > 0:
                    del swap_values[self.nCre - 1]

                for i in swap_values:
                    pattern = list(range(n))
                    pattern[i], pattern[i + 1] = i + 1, i
                    self.symmetries.append(symmetry(tuple(pattern), -1))

            # Add bra/ket symmetries for ground-state RDMs
            if not self.trans_rdm:
                if self.nInd % 2 == 0:
                    reversed_range = tuple(range(self.nInd - 1, -1, -1))
                    self.symmetries.append(symmetry(reversed_range, 1))
                else:
                    print('WARN: trans_rdm set to False with odd number of cre/des operators, switching trans_rdm to True!')
                    self.trans_rdm = True

    #------------------------------------------------------------------------------------------------

    # def _comparison_key(self):
    #     """creDesTensor comparison key - used for ordering"""
    #     return (self.name, self.indices, self.symmetries)

    def __eq__(self, other):
        if isinstance(other, creDesTensor):
            return self._comparison_key() == other._comparison_key()
        return False

    def __lt__(self, other):
        if isinstance(other, creDesTensor):
            return self._comparison_key() < other._comparison_key()
        elif isinstance(other, (creOp, desOp)):
            return True
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("A creDesTensor object may only be compared to another tensor")

    def copy(self):
        ops = []

        for i in range(self.nCre):
            ops.append(creOp(self.indices[i]))

        for i in range(self.nCre,len(self.indices)):
            ops.append(desOp(self.indices[i]))

        return creDesTensor(list(ops), self.trans_rdm, self.symmetries)


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class creOp(tensor):
    """A tensor representation for a creation operator."""

    freelyCommutes = False

    def __init__(self, indices):
        # Process index
        if isinstance(indices, list):
            if len(indices) != 1 or not isinstance(indices[0], index):
                raise TypeError("indices must be an index or a list of indices with length 1")
            input_index = indices[0].copy()
        elif isinstance(indices, index):
            input_index = indices.copy()
        else:
            raise TypeError("indices must be an index or a list of indices with length 1")

        # Initialize attributes
        self.permutations = None
        self.factors = None
        self.name = "cre"
        self.indices = [input_index]
        self.symmetries = []

    #------------------------------------------------------------------------------------------------

    # def _comparison_key(self):
    #     """creOp comparison key - used for ordering"""
    #     return (self.name, self.indices, self.symmetries)

    def __eq__(self, other):
        if isinstance(other, creOp):
            return self._comparison_key() == other._comparison_key()
        return False

    def __lt__(self, other):
        if isinstance(other, creOp):
            return self._comparison_key() < other._comparison_key()
        elif isinstance(other, desOp):
            return True
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("A creOp object may only be compared to another tensor")

    #------------------------------------------------------------------------------------------------

    def copy(self):
        return creOp(self.indices)

    #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class desOp(tensor):
    """A tensor representation for a destruction operator."""

    freelyCommutes = False

    def __init__(self, indices):
        # Process index
        if isinstance(indices, list):
            if len(indices) != 1 or not isinstance(indices[0], index):
                raise TypeError("indices must be an index or a list of indices with length 1")
            input_index = indices[0].copy()
        elif isinstance(indices, index):
            input_index = indices.copy()
        else:
            raise TypeError("indices must be an index or a list of indices with length 1")

        # Initialize attributes
        self.permutations = None
        self.factors = None
        self.name = "des"
        self.indices = [input_index]
        self.symmetries = []

    #------------------------------------------------------------------------------------------------

    # def _comparison_key(self):
    #     """desOp comparison key - used for ordering"""
    #     return (self.name, self.indices, self.symmetries)

    def __eq__(self, other):
        if isinstance(other, desOp):
            return self._comparison_key() == other._comparison_key()
        return False

    def __lt__(self, other):
        if isinstance(other, desOp):
            return self._comparison_key() < other._comparison_key()
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("A desOp object may only be compared to another tensor")

    #------------------------------------------------------------------------------------------------

    def copy(self):
        return desOp(self.indices)

    #------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

