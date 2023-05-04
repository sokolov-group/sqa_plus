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

class sqaSpinBasis(object):
    "A class to define types of spin-basis"

    def __init__(self, value):
        self.spin_basis = value        

    def __get__(self, instance, owner):
        return (instance.spin_basis == self.spin_basis)

    def __set__(self, instance, value):
        if type(value) == type(True):
            if value:
                instance.spin_basis = self.spin_basis
        else:
            raise Exception("Spin-basis should be 'True' or 'False' ...")

class sqaOptions(object):
    "A class to hold options for the SecondQuantizationAlgebra+."

    spin_orbital = sqaSpinBasis('spin_orbital')
    spin_adapted = sqaSpinBasis('spin_adapted')
    spin_integrated = sqaSpinBasis('spin_integrated')

    def __init__(self):
        # Set default options
        self.verbose = False

        # Indices types
        self.alpha_type = ('alpha',)
        self.beta_type = ('beta',)

        self.core_type = ('core',)
        self.active_type = ('active',)
        self.virtual_type = ('virtual',)

        self.cvs_core_type = ('cvs',)
        self.cvs_valence_type = ('valence',)

        # Spin basis
        self.spin_basis = 'spin_orbital'
        self.explicit_spin_cases = False

        # Core-Valence Separation approach
        self.cvs_approach = False

# Create an object of the options class
options = sqaOptions()
