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

import secondQuantizationAlgebra as sqa

# Define short names for type groups
ta = sqa.options.alpha_type
tb = sqa.options.beta_type
tc = sqa.options.core_type
tt = sqa.options.active_type
tv = sqa.options.virtual_type

# Define core indices
ac = [sqa.index('c%i' %i, [ta,tc], False) for i in range(10)]
bc = [sqa.index('c%i' %i, [tb,tc], False) for i in range(10)]

# Define active indices
at = [sqa.index('a%i' %i, [ta,tt], True) for i in range(10)]
bt = [sqa.index('a%i' %i, [tb,tt], True) for i in range(10)]

terms = []
terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.desOp(ac[0])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.creOp(ac[1]), sqa.desOp(ac[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[0]), sqa.tensor('r0', [at[0], at[1]], []), sqa.creOp(ac[1]), \
                                sqa.desOp(ac[0]), sqa.desOp(at[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[0]), sqa.creOp(ac[1]), \
                                sqa.desOp(ac[0]), sqa.desOp(at[1]), sqa.desOp(ac[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.creOp(ac[1]), sqa.creOp(ac[0]), \
                                sqa.desOp(ac[0]), sqa.desOp(at[1]), sqa.desOp(ac[1])]))
terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[1]), sqa.creOp(ac[0]), \
                                sqa.desOp(ac[1]), sqa.desOp(ac[0]), sqa.desOp(ac[1])]))

print ""
print "Testing removeCoreOpPairs function."
print ""
print "Core indices labeled with c."
print "Active indices labeled with a."

print ""
print "initial terms:"
for t in terms:
  print t

sqa.removeCoreOpPairs(terms)

output = ""
for t in terms:
  output += str(t) + "\n"

correctOutput = \
" (   1.00000) \n" + \
" (  -1.00000) cre(c0) \n" + \
" (   1.00000) cre(a0) r0(a0,a1) cre(c1) des(a1) \n" + \
" (   1.00000) cre(a0) des(a1) \n" + \
" (  -1.00000) cre(c0) cre(c0) des(c0) des(a1) \n" + \
" (  -1.00000) cre(a0) cre(c1) des(c1) des(c1) \n"

print ""
print "terms after core operator pair removal:"
print output

if output == correctOutput:
  print "Test passed!"
else:
  print "Test failed.  Correct output is:"
  print correctOutput
print ""
