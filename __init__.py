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

import sys
print("""
\n----------------------------------------------------------------------------------
sqa_plus: Code generator for quasi-particle systems.
Copyright 2009-2022 SecondQuantizationAlgebra Developers. All Rights Reserved.
Available at https://github.com/sokolov-group/sqa_plus

Licensed under the GNU General Public License v3.0;

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

----------------------------------------------------------------------------------""")
sys.stdout.flush()

from sqaOptions import \
  options

from sqaIndex import \
  index

from sqaSymmetry import \
  symmetry

from sqaTensor import \
  creOp,              \
  desOp,              \
  kroneckerDelta,     \
  sfExOp,             \
  creDesTensor,       \
  tensor

from sqaTerm import  \
  combineTerms,      \
  multiplyTerms,     \
  removeCoreOpPairs, \
  removeCoreOps_sf,  \
  removeVirtOps_sf,  \
  sortOps,           \
  term,              \
  termChop

from sqaNormalOrder import \
  normalOrder

from sqaIntermediates import \
  genIntermediates

from sqaCommutator import \
  commutator

from sqaMisc import       \
  allDifferent,           \
  get_num_perms,          \
  makePermutations,       \
  makeTuples,             \
  assign_rdm_types,       \
  combine_transpose,      \
  convert_ops_to_rdms_so

from sqaDecomposition_sf import  \
  decomp_3op_to_2op_2rdm_sf,     \
  decomp_3ops_to_2ops_2rdms_sf,  \
  decomp_3rdm_to_2rdm_sf,        \
  decomp_3rdms_to_2rdms_sf,      \
  decomp_4op_to_2op_2rdm_sf,     \
  decomp_4ops_to_2ops_2rdms_sf,  \
  decomp_4rdm_to_2rdm_sf,        \
  decomp_4rdm_to_3rdm_sf,        \
  decomp_4rdms_to_2rdms_sf,      \
  decomp_4rdms_to_3rdms_sf

from sqaDecomposition_so import  \
  decomp_3op_to_2op_2rdm_so,     \
  decomp_3op_to_2op_3rdm_so,     \
  decomp_3ops_to_2ops_2rdms_so,  \
  decomp_3ops_to_2ops_3rdms_so,  \
  decomp_3rdm_to_2rdm_so,        \
  decomp_3rdms_to_2rdms_so,      \
  decomp_4op_to_2op_2rdm_so,     \
  decomp_4op_to_2op_3rdm_so,     \
  decomp_4ops_to_2ops_2rdms_so,  \
  decomp_4ops_to_2ops_3rdms_so,  \
  decomp_4rdm_to_2rdm_so,        \
  decomp_4rdm_to_3rdm_so,        \
  decomp_4rdms_to_2rdms_so,      \
  decomp_4rdms_to_3rdms_so

from sqaMatrixBlock import    \
  matrixBlock,                \
  dummyLabel,                 \
  filterVirtual,              \
  filterCore,                 \
  normalOrderCore,            \
  sortOpsCore,                \
  contractDeltaFuncs_nondummy

from sqaHeff import         \
  Heff,                     \
  dyallH,                   \
  dyallH_act,               \
  Tamplitude,               \
  Tamplitude_excitation,    \
  Tamplitude_deexcitation,  \
  Vperturbation,            \
  getT,                     \
  getV                      \

from sqaEinsum import genEinsum

from sqaIndexList import    \
  indexList,                \
  create_dummy_indices_list

from sqaSpinAdapted import convertSpinIntegratedToAdapted
