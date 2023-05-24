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

class dummyIndexList:
    "A list of indices created on demand."

    def __init__(self, index_name, index_type, sqa_default_index):

        if self.check_index_name_conflicts(index_name, sqa_default_index):
            self.index_name = index_name
        else:
            raise Exception("Index name {:} is reserved to SQA+ internal dummy indices.")

        self.index_type = index_type
        self.index_reserved = 0
        self.index_dummy = True
        self.index_user_defined = False

    def check_index_name_conflicts(self, index_name, sqa_default_index):
        sqa_default_index_names = ['cc%i', 'cc%ia', 'cc%ib',
                                   'xx%i', 'xx%ia', 'xx%ib',
                                   'vv%i', 'vv%ia', 'vv%ib',
                                   'aa%i', 'aa%ia', 'aa%ib',
                                   'ee%i', 'ee%ia', 'ee%ib']

        return sqa_default_index or ((not sqa_default_index) and (index_name not in sqa_default_index_names))

    def new_indices(self, index_external):
        new_index_list = []
        for index_number in range(self.index_reserved, self.index_reserved + index_external):
            new_index_list.append(index(self.index_name % index_number, self.index_type, self.index_dummy, self.index_user_defined))
        self.index_reserved += index_external

        return new_index_list

    def new_index(self):
        index_number = self.index_reserved
        self.index_reserved += 1

        return index(self.index_name % index_number, self.index_type, self.index_dummy, self.index_user_defined)

class sqaIndexLists:
    "A class to define all types of required dummy indices in actual SQA+ implementation."

    def __init__(self):

        self.sqa_default_index = True

        # Define spin-orbital indices
        ## Core dummy indices
        self.core = dummyIndexList('cc%i', [options.core_type], self.sqa_default_index)

        ## CVS core dummy indices
        self.cvs_core = dummyIndexList('xx%i', [options.cvs_core_type], self.sqa_default_index)

        self.cvs_valence = dummyIndexList('vv%i', [options.cvs_valence_type], self.sqa_default_index)

        ## Active dummy indices
        self.active = dummyIndexList('aa%i', [options.active_type], self.sqa_default_index)

        ## Virtual dummy indices
        self.virtual = dummyIndexList('ee%i', [options.virtual_type], self.sqa_default_index)

        # Define spin-integrated indices
        ## Core dummy indices
        self.core_alpha = dummyIndexList('cc%ia', [options.alpha_type, options.core_type], self.sqa_default_index)
        self.core_beta  = dummyIndexList('cc%ib', [options.beta_type,  options.core_type], self.sqa_default_index)

        ## CVS core dummy indices
        self.cvs_core_alpha = dummyIndexList('xx%ia', [options.alpha_type, options.cvs_core_type], self.sqa_default_index)
        self.cvs_core_beta  = dummyIndexList('xx%ib', [options.beta_type,  options.cvs_core_type], self.sqa_default_index)

        self.cvs_valence_alpha = dummyIndexList('vv%ia', [options.alpha_type, options.cvs_valence_type], self.sqa_default_index)
        self.cvs_valence_beta  = dummyIndexList('vv%ib', [options.beta_type,  options.cvs_valence_type], self.sqa_default_index)
 
        ## Active dummy indices
        self.active_alpha = dummyIndexList('aa%ia', [options.alpha_type, options.active_type], self.sqa_default_index)
        self.active_beta  = dummyIndexList('aa%ib', [options.beta_type,  options.active_type], self.sqa_default_index)

        ## Virtual dummy indices
        self.virtual_alpha = dummyIndexList('ee%ia', [options.alpha_type, options.virtual_type], self.sqa_default_index)
        self.virtual_beta  = dummyIndexList('ee%ib', [options.beta_type,  options.virtual_type], self.sqa_default_index)

# Create an object of the indexLists class
indexLists = sqaIndexLists()
