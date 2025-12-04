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

from functools import total_ordering
from .sqaOptions import options

@total_ordering
class index:
    "A class for tensor and operator indices."

    def __init__(self, name, indexType = (), isSummed = False, userDefined = True):
        # Initialize name
        self.name = str(name)

        # Initialize flag for whether the index is summed over (dummy index)
        if not isinstance(isSummed, bool):
            raise TypeError("isSummed must be a boolean")
        self.isSummed = isSummed

        # Initialize if index is user-defined
        if not isinstance(userDefined, (bool, str)):
            raise TypeError("userDefined must be a boolean or a string")

        if userDefined is True:
            self.userDefined = self.name
            options.add_user_defined_index(self.name)
        else:
            self.userDefined = userDefined

        # Initialize index types
        indexError = "indexType must be a list/tuple of lists/tuples of strings"
        indType = []
        for group in indexType:
            if not isinstance(group, (list, tuple)) or not all(isinstance(s, str) for s in group):
                raise TypeError(indexError)
            indType.append(tuple(sorted(group)))
        self.indType = tuple(indType)

    def __eq__(self, other):
        if not isinstance(other,index):
            raise False
        return (self.isSummed == other.isSummed and self.name == other.name and self.indType == other.indType)
 
    def __lt__(self, other):
        if not isinstance(other,index):
            raise ValueError("can only compare index class with other index class objects.")
        if self.isSummed != other.isSummed:
            return self.isSummed < other.isSummed
        if self.name != other.name:
            return self.name < other.name
        return self.indType < other.indType

    def tup(self):
        "Returns a tuple representation of the index. The return object in unmutable and thus can be used as a dictionary key."
        return (self.name, self.indType, self.isSummed, self.userDefined)

    def copy(self):
        "Returns a deep copy of the index"
        return index(self.name, self.indType, self.isSummed, self.userDefined)

    def rename(self):
        "Rename index according to user defined name."
        if isinstance(self.userDefined, str):
            self.name = self.userDefined

# SecondQuantizationAlgebra Plus
#
# Functions implemented to automate test index types
#
# Author: Carlos E. V. de Moura <carlosevmoura@gmail.com>

def is_spin_integrated_index_type(indice_types):
    """Returns True if indices contain spin-integrated (alpha or beta) types."""
    if isinstance(indice_types, index):
        indice_types = indice_types.indType
    return any(ind in (options.alpha_type, options.beta_type) for ind in indice_types)

def get_spin_index_type(indice_types):
    """Returns spin index type of indices."""
    if isinstance(indice_types, index):
        indice_types = indice_types.indType
    return next((ind for ind in indice_types if ind in (options.alpha_type, options.beta_type)), '')

def get_spatial_index_type(indice_types):
    """Returns spatial index type of indices."""
    if isinstance(indice_types, index):
        indice_types = indice_types.indType
    return next((ind for ind in indice_types if ind not in (options.alpha_type, options.beta_type)), '')

def is_index_type(indice_types, target_index_types):
    """Returns True of any index matches one in target_index_types."""
    return any(ind in target_index_types for ind in indice_types)

def is_core_index_type(index_type):
    """Returns True of any index is a core type."""
    spatial_index_type = get_spatial_index_type(index_type)
    core_types = (options.core_type, options.cvs_core_type, options.cvs_valence_type)
    return any(is_index_type(spatial_index_type, ct) for ct in core_types)

def is_cvs_index_type(index_type):
    """Returns True of any index is a cvs type."""
    spatial_index_type = get_spatial_index_type(index_type)
    cvs_types = (options.cvs_core_type, options.cvs_valence_type)
    return any(is_index_type(spatial_index_type, ct) for ct in cvs_types)

def is_cvs_core_index_type(index_type):
    """Returns True of any index is a cvs core type."""
    spatial_index_type = get_spatial_index_type(index_type)
    return is_index_type(spatial_index_type, options.cvs_core_type)

def is_cvs_valence_index_type(index_type):
    """Returns True of any index is a cvs valence type."""
    spatial_index_type = get_spatial_index_type(index_type)
    return is_index_type(spatial_index_type, options.cvs_valence_type)

def is_active_index_type(index_type):
    """Returns True of any index is an active type."""
    spatial_index_type = get_spatial_index_type(index_type)
    return is_index_type(spatial_index_type, options.active_type)

def is_virtual_index_type(index_type):
    """Returns True of any index is a virtual type."""
    spatial_index_type = get_spatial_index_type(index_type)
    return is_index_type(spatial_index_type, options.virtual_type)

def is_alpha_index_type(index_type):
    """Returns True of any index is an alpha-spin type."""
    spin_index_type = get_spin_index_type(index_type)
    return is_index_type(spin_index_type, options.alpha_type)

def is_beta_index_type(index_type):
    """Returns True of any index is a beta-spin type."""
    spin_index_type = get_spin_index_type(index_type)
    return is_index_type(spin_index_type, options.beta_type)
