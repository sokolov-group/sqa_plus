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
    "A class to represent tensors in operator algebra. Integrals and density matrices are examples."

    #------------------------------------------------------------------------------------------------

    freelyCommutes = True

    #------------------------------------------------------------------------------------------------

    def __init__(self, name, indices = [], symmetries = []):

        # Initialize data
        (self.permutations,self.factors) = (None,None)
        self.indices = []
        self.symmetries = []

        # Process name
        self.name = str(name)

        # Process indices
        indicesError = "indices must be a list of index objects"
        if not isinstance(indices, list):
            raise TypeError(indicesError)
        for i in indices:
            if not isinstance(i, index):
                raise TypeError(indicesError)
            self.indices.append( i.copy() )

        # Process symmetries
        symmetryError = "symmetries must be a list of symmetry objects"
        if not isinstance(symmetries, list):
            raise TypeError(symmetryError)
        for sym in symmetries:
            if not isinstance(sym, symmetry):
                raise TypeError(symmetryError)
            for s in self.symmetries:
                if s.pattern == sym.pattern:
                    raise ValueError("a tensor cannot have two symmetries with the same pattern")
            self.symmetries.append( sym.copy() )

    #------------------------------------------------------------------------------------------------

    def __eq__(self, other):
        if not isinstance(other, tensor):
            return False

        # If other belongs to a tensor subclass, use the subclass's comparison method
        if isinstance(other, (kroneckerDelta, creOp, desOp, creDesTensor, sfExOp)):
            return other == self

        return (self.name == other.name and self.indices == other.indices and self.symmetries == other.symmetries)

    def __lt__(self, other):
        if not isinstance(other, tensor):
            raise TypeError("A tensor may only be compared to another tensor")

        # If other belongs to a tensor subclass, use the subclass's comparison method
        if isinstance(other, (kroneckerDelta, creOp, desOp, creDesTensor, sfExOp)):
            return self < other

        # compare names
        if self.name != other.name:
            return self.name < other.name
        # compare indices
        if self.indices != other.indices:
            return self.indices < other.indices
        # compare symmetries
        return self.symmetries < other.symmetries

    #------------------------------------------------------------------------------------------------

    def __str__(self):
        retval = self.name + "("
        for i in range(len(self.indices)):
            retval += self.indices[i].name #+ " " + str(self.indices[i].type)
            if i < len(self.indices)-1:
                retval += ","
        retval += ")"
        return retval

    #------------------------------------------------------------------------------------------------

    def copy(self):
        "Returns a copy of the tensor"
        retval = tensor(self.name, self.indices, self.symmetries)
        if self.permutations != None and self.factors != None:
            retval.permutations = [ perm + [] for perm in self.permutations ]
            retval.factors = self.factors + []
        return retval

    #------------------------------------------------------------------------------------------------

    def symPermutes(self, force = False):
        "Returns the index permutations and resulting factors allowed by the tensor's symmetry"

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
        "Sort indices alphabetically within symmetry constraints. Returns the resulting symmetry factor."

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

    def hasIndex(self,i):
        "Returns True if i is one of the tensor's indices and False otherwise."
        if not isinstance(i,index):
            raise TypeError("i must be of the index class")
        return (i in self.indices)

    #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class kroneckerDelta(tensor):
    "A tensor representation of the kronecker delta function."

    #------------------------------------------------------------------------------------------------

    freelyCommutes = True
    symmetries = [symmetry((1,0),1)]
    name = "kdelta"

    #------------------------------------------------------------------------------------------------

    def __init__(self,indices):
        if len(indices) != 2:
            raise ValueError("The kronecker delta function takes exactly two indices")
        self.indices = []
        for i in indices:
            self.indices.append( i.copy() )
        (self.permutations,self.factors) = (None,None)

    #------------------------------------------------------------------------------------------------

    def __eq__(self, other):
        if isinstance(other, kroneckerDelta):
            return (self.name == other.name and self.indices == other.indices and self.symmetries == other.symmetries)
        return False

    # kroneckerDelta class is less than the creOp, desOp, creDesTensor and sfExOp sub classes
    def __lt__(self, other):
        if isinstance(other, kroneckerDelta):
            if self.name != other.name:
                return self.name < other.name
            if self.indices != other.indices:
                return self.indices < other.indices
            return self.symmetries < other.symmetries
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
        if len(indices)/2 != (len(indices)+1)/2:
            raise ValueError("A spin free excitation operator (the sfExOp class) must have an even number of indices")

        # Initialize order
        self.order = len(indices)/2

        # Initialize name
        self.name = "E%i" %self.order

        # Initialize indices
        self.indices = []
        for i in indices:
            self.indices.append( i.copy() )

        # Initialize permutations and factors
        (self.permutations,self.factors) = (None,None)

        # Initialize symmetries
        self.symmetries = []
        for i in range(self.order-1):
            if i == 0:
                temp_tup = (1,)
            else:
                temp_tup = (0,)
            for j in range(1,2*self.order):
                if j == i:
                    temp_tup = temp_tup + (i+1,)
                elif j == i+1:
                    temp_tup = temp_tup + (i,)
                elif j == i+self.order:
                    temp_tup = temp_tup + (i+1+self.order,)
                elif j == i+1+self.order:
                    temp_tup = temp_tup + (i+self.order,)
                else:
                    temp_tup = temp_tup + (j,)
            self.symmetries.append(symmetry(temp_tup, 1))
                    

    #------------------------------------------------------------------------------------------------

    def __eq__(self, other):
        if isinstance(other, sfExOp):
            return (self.name == other.name and self.indices == other.indices and self.symmetries == other.symmetries)
        return False
    
    # sfExOp class is less than the creOp, desOp, and creDesTensor classes
    def __lt__(self, other):
        if isinstance(other, sfExOp):
            if self.name != other.name:
                return self.name < other.name
            if self.indices != other.indices:
                return self.indices < other.indices
            return self.symmetries < other.symmetries
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

    def __init__(self, ops, trans_rdm = False, symmetries = False):

        TypeErrorMessage = "ops must be a normal ordered list of creOp and desOp objects"
        if not isinstance(ops, list):
            raise TypeError(TypeErrorMessage)

        # Initialize list of cre/des operators
        self.ops = ops

        # Initialize name
        self.trans_rdm = trans_rdm
        if trans_rdm:
            self.name = 'trdm'
        else:
            self.name = 'rdm'

        # Initialize permutations and factors
        (self.permutations, self.factors) = (None, None)

        # Count the number of creation/destruction operators
        self.nCre = 0
        self.nDes = 0

        # Build the index list
        self.indices = []
        desFlag = False

        for op in ops:

            # Ensure normal-ordering
            if not isinstance(op, (creOp, desOp)):
                raise TypeError(TypeErrorMessage)

            if isinstance(op, desOp):
                self.nDes += 1
                desFlag = True

            elif isinstance(op, creOp):
                if desFlag:
                    raise TypeError(TypeErrorMessage)
                self.nCre += 1

            self.indices.append(op.indices[0].copy())

        # Create count of total indices
        self.nInd = self.nCre + self.nDes

        # Initialize symmetries
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
            if (len(self.indices) % 2 == 0) and self.trans_rdm == False:
                reversed_range = tuple(range(len(self.indices))[::-1])
                self.symmetries.append(symmetry(reversed_range, 1))

            # Print warning if number of indices is odd and trans_rdm is False
            if (len(self.indices) % 2 != 0) and self.trans_rdm == False:
                print ('trans_rdm flag is set to True, but an ODD number of cre/des operators are present. Switching trans_rdm flag to TRUE !!')
                self.trans_rdm == True

    def __eq__(self, other):
        if isinstance(other, creDesTensor):
            return (self.name == other.name and self.indices == other.indices and self.symmetries == other.symmetries)
        return False
    
    # creDesTensor class is less than the creOp and desOp classes
    def __lt__(self, other):
        if isinstance(other, creDesTensor):
            if self.name != other.name:
                return self.name < other.name
            if self.indices != other.indices:
                return self.indices < other.indices
            return self.symmetries < other.symmetries
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

class creDesTensor_original(tensor):
    """
    A tensor representation of a string of creation/destruction operators
    """

    #------------------------------------------------------------------------------------------------

    freelyCommutes = False
    name = "creDesTensor"

    #------------------------------------------------------------------------------------------------

    def __init__(self, ops):
        
        TypeErrorMessage = "ops must be a normal ordered list of creOp and desOp objects"
        if not type(ops) == type([]):
            raise TypeError(TypeErrorMessage)

        # Initialize permutations and factors
        (self.permutations,self.factors) = (None,None)

        # Build the index list and count the number of creation operators
        self.nCre = 0
        self.indices = []
        desFlag = False
        for op in ops:
            if (not isinstance(op, creOp)) and (not isinstance(op, desOp)):
                raise TypeError(TypeErrorMessage)
            if isinstance(op, desOp):
                desFlag = True
            if isinstance(op, creOp):
                self.nCre += 1
                if desFlag:
                    raise TypeError(TypeErrorMessage)
            self.indices.append(op.indices[0].copy())

        # Initialize symmetries
        self.symmetries = []
        swapValues = range(len(self.indices)-1)
        if self.nCre > 0:
            del(swapValues[self.nCre-1])
        for i in swapValues:
            if i == 0:
                temp_tup = (1,)
            else:
                temp_tup = (0,)
            for j in range(1,len(self.indices)):
                if j == i:
                    temp_tup = temp_tup + (i+1,)
                elif j == i+1:
                    temp_tup = temp_tup + (i,)
                else:
                    temp_tup = temp_tup + (j,)
            self.symmetries.append(symmetry(temp_tup, -1))
                    
    #------------------------------------------------------------------------------------------------

    def __eq__(self, other):
        if isinstance(other, creDesTensor):
            return (self.name == other.name and self.indices == other.indices and self.symmetries == other.symmetries)
        return False
    
    # creDesTensor class is less than the creOp and desOp classes
    def __lt__(self, other):
        if isinstance(other, creDesTensor):
            if self.name != other.name:
                return self.name < other.name
            if self.indices != other.indices:
                return self.indices < other.indices
            return self.symmetries < other.symmetries
        elif isinstance(other, (creOp, desOp)):
            return True
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("A creDesTensor object may only be compared to another tensor")
 
    #------------------------------------------------------------------------------------------------

    def copy(self):
        ops = []
        for i in range(self.nCre):
            ops.append(creOp(self.indices[i]))
        for i in range(self.nCre,len(self.indices)):
            ops.append(desOp(self.indices[i]))
        return creDesTensor(ops)

    #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class creOp(tensor):
    """
    A tensor representation for a creation operator
    """

    #------------------------------------------------------------------------------------------------

    freelyCommutes = False
    name = "cre"
    symmetries = []

    #------------------------------------------------------------------------------------------------

    def __init__(self, indices):

        # Initialize index
        if isinstance(indices, list) and len(indices) == 1 and isinstance(indices[0], index):
            inputIndex = indices[0].copy()
        elif isinstance(indices, index):
            inputIndex = indices.copy()
        else:
            raise TypeError("indices must be an index or a list of indices with length 1")
        self.indices = [inputIndex]

        # Initialize permutations and factors
        (self.permutations,self.factors) = (None,None)

    #------------------------------------------------------------------------------------------------

    def __eq__(self, other):
        if isinstance(other, creOp):
            return (self.name == other.name and self.indices == other.indices and self.symmetries == other.symmetries)
        return False
    
    # creOp class is less than the desOp class
    def __lt__(self, other):
        if isinstance(other, creOp):
            if self.name != other.name:
                return self.name < other.name
            if self.indices != other.indices:
                return self.indices < other.indices
            return self.symmetries < other.symmetries
        elif isinstance(other, desOp):
            return True
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("An creOp object may only be compared to another tensor")
 
    #------------------------------------------------------------------------------------------------

    def copy(self):
        return creOp(self.indices)

    #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@total_ordering
class desOp(tensor):
    """
    A tensor representation for a destruction operator
    """

    #------------------------------------------------------------------------------------------------

    freelyCommutes = False
    name = "des"
    symmetries = []

    #------------------------------------------------------------------------------------------------

    def __init__(self, indices):

        # Initialize index
        if isinstance(indices, list) and len(indices) == 1 and isinstance(indices[0], index):
            inputIndex = indices[0].copy()
        elif isinstance(indices, index):
            inputIndex = indices.copy()
        else:
            raise TypeError("indices must be an index or a list of indices with length 1")
        self.indices = [inputIndex]

        # Initialize permutations and factors
        (self.permutations,self.factors) = (None,None)

    #------------------------------------------------------------------------------------------------

    def __eq__(self, other):
        if isinstance(other, desOp):
            return (self.name == other.name and self.indices == other.indices and self.symmetries == other.symmetries)
        return False
    
    # desOp class is greater than other tensor subclasses
    def __lt__(self, other):
        if isinstance(other, desOp):
            if self.name != other.name:
                return self.name < other.name
            if self.indices != other.indices:
                return self.indices < other.indices
            return self.symmetries < other.symmetries
        elif isinstance(other, tensor):
            return False
        else:
            raise TypeError("An desOp object may only be compared to another tensor")
 
    #------------------------------------------------------------------------------------------------

    def copy(self):
        return desOp(self.indices)

    #------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

