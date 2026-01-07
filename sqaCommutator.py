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
# Computes the commutator, including normal ordering,
# between terms or lists of terms.
#

from .sqaTerm import term, combineTerms, multiplyTerms, termChop
from .sqaNormalOrder import normalOrder

def commutator(leftInput, rightInput, contract = True, combine = True):

    # Convert inputs that are terms into lists of terms
    leftTerms  = [leftInput]  if isinstance(leftInput, term)  else leftInput
    rightTerms = [rightInput] if isinstance(rightInput, term) else rightInput

    # Check input integrity
    TypeErrorMessage = "commutator inputs must be terms or lists of terms"
    if not isinstance(leftTerms, list) or not isinstance(rightTerms, list):
            raise TypeError(TypeErrorMessage)
    if not all(isinstance(t, term) for t in leftTerms + rightTerms):
        raise TypeError(TypeErrorMessage)

    # Construct terms resulting from the commutator (before normal ordering)
    preNOTerms = []
    for l in leftTerms:
        for r in rightTerms:
            preNOTerms.append(multiplyTerms(l, r))
            pdt = multiplyTerms(r, l)
            pdt.scale(-1)
            preNOTerms.append(pdt)

    # For each term, apply Wick's theorem to convert it to normal order
    noTerms = []
    for t in preNOTerms:
        noTerms.extend(normalOrder(t))
    del(preNOTerms)

    # Contract any delta functions resulting from the normal ordering
    if contract:
        for t in noTerms:
            t.contractDeltaFuncs()

    # Remove any terms that are zero
    termChop(noTerms)

    # Combine any like terms unless told not to
    if combine:
        combineTerms(noTerms)

    # Return result
    return noTerms
