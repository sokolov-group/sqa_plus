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
# Author: Carlos E. V. de Moura <carlosevmoura@gmail.com>
#

from sqaIndex import index
from sqaOptions import options

class indexList:
    "A list of indices created on demand."

    def __init__(self, index_name, index_type, dummy):
        self.index_name = index_name
        self.index_type = index_type
        self.index_reserved = 0
        self.dummy = dummy
        self.user_defined = False

    def new_indices(self, index_external):
        new_index_list = []
        for index_number in range(self.index_reserved, self.index_reserved + index_external):
            new_index_list.append(index(self.index_name % index_number, self.index_type, self.dummy, self.user_defined))
        self.index_reserved += index_external

        return new_index_list

    def new_index(self):
        index_number = self.index_reserved
        self.index_reserved += 1

        return index(self.index_name % index_number, self.index_type, self.dummy, self.user_defined)

def create_dummy_indices_list(spin_integrated = False):
    "A function which returns lists of dummy indices of all classes."

    def create_dummy_indices_spin_orbital():
        # Define operator types
        tg_cor = options.core_type
        tg_act = options.active_type
        tg_vir = options.virtual_type

        dummy = True

        # Core dummy indices
        cor_inds = indexList('c%ia', [tg_cor], dummy)

        # Active dummy indices
        act_inds = indexList('a%ia', [tg_act], dummy)

        # Virtual dummy indices
        vir_inds = indexList('v%ia', [tg_vir], dummy)

        return (cor_inds, act_inds, vir_inds)

    def create_dummy_indices_spin_integrated():
        # Define operator types
        tg_cor = options.core_type
        tg_act = options.active_type
        tg_vir = options.virtual_type

        tg_alpha = options.alpha_type
        tg_beta = options.beta_type

        dummy = True

        # Core dummy indices
        cor_alpha_inds = indexList('c%ia', [tg_alpha, tg_cor], dummy)
        cor_beta_inds  = indexList('c%ib', [tg_beta,  tg_cor], dummy)

        # Active dummy indices
        act_alpha_inds = indexList('a%ia', [tg_alpha, tg_act], dummy)
        act_beta_inds  = indexList('a%ib', [tg_beta,  tg_act], dummy)

        # Virtual dummy indices
        vir_alpha_inds = indexList('v%ia', [tg_alpha, tg_vir], dummy)
        vir_beta_inds  = indexList('v%ib', [tg_beta,  tg_vir], dummy)

        return (cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds)

    if spin_integrated:
        dummy_indices = create_dummy_indices_spin_integrated()
    else:
        dummy_indices = create_dummy_indices_spin_orbital()

    return dummy_indices
