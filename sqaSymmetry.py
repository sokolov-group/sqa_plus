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
# The symmetry class represents a symmetry permutation of a tensor.
#
# The class consists of a permutation pattern and the constant factor
# that is produced when applying the permutation to the tensor's indices.
#
# The permutation, or pattern, for a tensor with n indices is represented
# by a tuple containing the integers 0 through n-1 in some order.
#
# The constant factor is either an int or a float.
#
# For example, the object created by symmetry((0,1,3,2), -1) implies that
# if you swap the tensor's last two indices, the result is equal to -1
# times the original tensor.

class symmetry:
  "A class representing symmetries in a tensor's indices"

  #------------------------------------------------------------------------------------------------

  def __init__(self,pattern,factor):
    patternError = False
    if type(pattern) != type((1,2)):
      patternError = True
    else:
      tempList = []
      for elem in pattern:
        if type(elem) != type(1):
          patternError = True
          break
        else:
          tempList.append(elem)
      if not patternError:
        tempList.sort()
        patternError = ( tempList != range(len(pattern)) )
    if patternError:
      raise ValueError, "pattern must be a tuple of contiguous, non-negative integers including zero. They need not be in order."
    else:
      self.pattern = pattern
    if (factor == 1 or factor == -1) and type(factor) == type(1):
      self.factor = factor
    elif type(factor) == type(1) or type(factor) == type(1.0):
      self.factor = float(factor)
    else:
      raise TypeError, "factor must be a float or an int."

  #------------------------------------------------------------------------------------------------
  
  def __cmp__(self,other):
    if not isinstance(other,symmetry):
      raise TypeError, "can only compare a symmetry object to other symmetry objects"
    return cmp(self.pattern,other.pattern)

  #------------------------------------------------------------------------------------------------
  
  def __str__(self):
    retval = str(self.pattern)
    retval += "%15.5e" %self.factor
    return retval

  #------------------------------------------------------------------------------------------------

  def copy(self):
    "Returns a copy of the symmetry"
    return symmetry(self.pattern, self.factor)

  #------------------------------------------------------------------------------------------------


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


