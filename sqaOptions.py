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
            if value == 'spin_adapted':
                instance.reorder_amplitudes = True
                instance.legacy_ordering = False
        else:
            raise Exception("Spin-basis must be True or False ...")

class sqaTensorNotation(object):
    "A class to define types of tensors notations"

    def __init__(self, value):
        self.tensors_notation = value        

    def __get__(self, instance, owner):
        return (instance.tensors_notation == self.tensors_notation)

    def __set__(self, instance, value):
        if type(value) == type(True):
            if value:
                instance.tensors_notation = self.tensors_notation
        else:
            raise Exception("Tensors notation must be True or False ...")

class sqaOptions(object):
    "A class to hold options for the SecondQuantizationAlgebra+."

    spin_orbital = sqaSpinBasis('spin_orbital')
    spin_adapted = sqaSpinBasis('spin_adapted')
    spin_integrated = sqaSpinBasis('spin_integrated')

    chemists_notation = sqaTensorNotation('chemists')
    physicists_notation = sqaTensorNotation('physicists')

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

        # Tensors notations
        self.tensors_notation = 'physicists'

        # Core-Valence Separation approach
        self.cvs_approach = False

        # List of indices names defined by the user
        self.user_defined_indices = []

        # Options about custom reorder of tensor indices
        self.reorder_amplitudes = False
        self.legacy_ordering = True

        # MatrixBlock options
        self.matrixBlock = lambda:None
        self.matrixBlock.remove_trans_rdm_constant = False

        # genEinsum options
        self.genEinsum = lambda:None
        self.genEinsum.lhs_string = None
        self.genEinsum.indices_string = None
        self.genEinsum.spin_integrated_tensors = True

        self.genEinsum.trans_rdm = False
        self.genEinsum.trans_indices_string = None
        self.genEinsum.suffix = None

        self.genEinsum.cvs_indices_list = None
        self.genEinsum.valence_indices_list = None
        self.genEinsum.cvs_tensors = False

        self.genEinsum.remove_trans_rdm_constant = False
        self.genEinsum.remove_core_integrals = False
        self.genEinsum.intermediate_list = None

        self.genEinsum.opt_einsum_terms = True
        self.genEinsum.optimize = True

        # genIntermediates options
        self.genIntermediates = lambda:None
        self.genIntermediates.trans_rdm = False
        self.genIntermediates.factor_depth = 1

    def add_user_defined_index(self, name):
        if name not in self.user_defined_indices:
            self.user_defined_indices.append(name)

    def print_header(self, string):
        print("\n" + " {:} ".format(string).center(120, "-") + "\n")

    def print_divider(self):
        print("\n" + "-" * 120 + "\n")

# Create an object of the options class
options = sqaOptions()
