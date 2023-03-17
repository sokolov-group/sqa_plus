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

class sqaOptions:
  "A class to hold options for the second quantization algebra program."

  def __init__(self):
    # Set default options
    self.verbose = False
    self.alpha_type = ('alpha',)
    self.beta_type = ('beta',)
    self.core_type = ('core',)
    self.active_type = ('active',)
    self.virtual_type = ('virtual',)

# Create an object of the options class
options = sqaOptions()
