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
# The index class is used to represent a tensor index.
# An index consists of three parts: a name, a list of types, and a bool.
#
# The index's name may be any string.
#
# The bool indicates whether the index is summed over or otherwise is a dummy index whose name may be changed.
# A value of True means that the various operator algebra functions are allowed the change the index's name.
# A value of False means that the index's name may not be changed.
#
# The list of types is actually a list of lists of strings. Each list of strings represents a type group.
# Examples of type groups are the index's spin type or whether the index is core, active, virtual, etc.
# The reason for this format is to evaluate Kronecker delta functions,
# which will evaluate to zero if there is no overlap in any of the two indices' type groups.
#
# For example, if index1 has indexType = [['alpha'], ['core', 'active']] and
# index2 has indexType = [['beta'], ['active']], the Kronecker delta between them will be zero
# because the first type group has no matching strings.
#
# If the first type group had been omitted, then the delta function would not be
# zero because both indices have 'active' as one of the types in the second type group.
# Note that while the type groups can be inputted as a list of lists of strings,
# they are actually stored as a tuple of tuples of strings.

from sqaOptions import options

class index:
	"A class for tensor and operator indices."

	def __init__(self, name, indexType = (), isSummed = False, userDefined = True):
		# Initialize index name
		self.name = str(name)

		# Initialize index types
		indType = []
		for l in indexType:
			if not ( type(l) in [type([]), type(())] ):
				raise TypeError, "indexType must be a list or tuple of lists or tuples of strings"
			indType.append([])
			indType[-1].extend(l)
			indType[-1].sort()
            # indType.sort()
		self.indType = []
		for l in indType:
			self.indType.append(())
			for s in l:
				if type(s) != type('a'):
					raise TypeError, "indexType must be a list or tuple of lists or tuples of strings"
				self.indType[-1] = self.indType[-1] + (s,)
		self.indType = tuple(self.indType)

		# Initialize flag for whether the index is summed over
		if type(isSummed) != type(True):
			raise TypeError, "isSummed must be True or False"
		self.isSummed = isSummed

		if type(userDefined) != type(True):
			raise TypeError, "userDefined must be True or False"
		self.userDefined = userDefined

	def __cmp__(self,other):
		if (not isinstance(other,index)):
			raise ValueError, "can only compare index class with other index class objects."
		retval = cmp(self.isSummed, other.isSummed)
		if retval != 0:
			return retval
		retval = cmp(self.name, other.name)
		if retval != 0:
			return retval
		retval = cmp(self.userDefined, other.userDefined)
		if retval != 0:
			return retval
		return cmp(self.indType, other.indType)

	def tup(self):
		"Returns a tuple representation of the index. The return object in unmutable and thus can be used as a dictionary key."
		return (self.name, self.indType, self.isSummed, self.userDefined)

	def copy(self):
		"Returns a deep copy of the index"
		return index(self.name, self.indType, self.isSummed, self.userDefined)

# SecondQuantizationAlgebra Plus
#
# Functions implemented to automate test index types
#
# Author: Carlos E. V. de Moura <carlosevmoura@gmail.com>

def is_spin_integrated_index_type(indice_types):
    spin_integrated = False

    if isinstance(indice_types, index):
        indice_types = indice_types.indType

    for index_type in indice_types:
        if index_type in (options.alpha_type, options.beta_type):
            spin_integrated = True

    return spin_integrated

def get_spin_index_type(indice_types):
    spin_index_types = ''

    if isinstance(indice_types, index):
        indice_types = indice_types.indType

    for index_type in indice_types:
        if index_type in (options.alpha_type, options.beta_type):
            spin_index_types = index_type

    return spin_index_types

def get_spatial_index_type(indice_types):
    spatial_index_types = ''

    if isinstance(indice_types, index):
        indice_types = indice_types.indType

    for index_type in indice_types:
        if not index_type in (options.alpha_type, options.beta_type):
            spatial_index_types = index_type

    return spatial_index_types

def is_index_type(indice_types, sqa_index_type):
    is_type = False

    for index_type in indice_types:
      if index_type in sqa_index_type:
          is_type = True

    return is_type

def is_core_index_type(index_type):
    spatial_index_type = get_spatial_index_type(index_type)
    return is_index_type(spatial_index_type, options.core_type)

def is_active_index_type(index_type):
    spatial_index_type = get_spatial_index_type(index_type)
    return is_index_type(spatial_index_type, options.active_type)

def is_virtual_index_type(index_type):
    spatial_index_type = get_spatial_index_type(index_type)
    return is_index_type(spatial_index_type, options.virtual_type)

def is_alpha_index_type(index_type):
    spin_index_type = get_spin_index_type(index_type)
    return is_index_type(spin_index_type, options.alpha_type)

def is_beta_index_type(index_type):
    spin_index_type = get_spin_index_type(index_type)
    return is_index_type(spin_index_type, options.beta_type)
