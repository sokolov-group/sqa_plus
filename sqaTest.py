# Copyright 2018-2022 SecondQuantizationAlgebra Developers. All Rights Reserved.
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
# Authors: Alexander Yu. Sokolov <sokolov.8@osu.edu>
#          Koushik Chatterjee <koushikchatterjee7@gmail.com>
#          Ilia Mazin <ilia.mazin@gmail.com>
#          Carlos E. V. de Moura <carlosevmoura@gmail.com>
#

import secondQuantizationAlgebra as sqa
import time

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# Test 1: Double commutator evaluation for a first order contribution to the  M_00 block of the effective Hamiltonian matrix
#
print ""
print "Starting test 1"
print ""
startTime = time.time()

# Create left hand side operator list
i = sqa.index('I', [tg_c])
x = sqa.index('X', [tg_a])

l_op  = [sqa.creOp(i), sqa.desOp(x)] # CA

# Create right hand side operator list
j = sqa.index('J', [tg_c])
y = sqa.index('Y', [tg_a])

r_op  = [sqa.creOp(y), sqa.desOp(j)] # CA

# Define order of the effective Hamiltonian
effH = []
effH = sqa.Heff(1)

# Define terms from operator lists
term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

print ("First Commutator")
term3 = sqa.commutator(effH, term1)

print ("Second Commutator")
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

test1_string_output = ''
for t in term5:
    test1_string_output += str(t) + '\n'

test1_correct_answer = " (  -1.00000) v(I,Y,J,X) \n" + \
                       " (   1.00000) kdelta(X,Y) v(I,x,J,y) rdm(x,y) \n" + \
                       " (   1.00000) v(I,Y,J,x) cre(x) des(X) \n" + \
                       " (   1.00000) v(I,x,J,X) cre(Y) des(x) \n" + \
                       " (  -1.00000) kdelta(X,Y) v(I,x,J,y) cre(y) des(x) \n" + \
                       " (  -1.00000) v(I,x,J,y) rdm(x,y) cre(Y) des(X) \n" + \
                       " (  -1.00000) v(I,x,J,y) cre(Y) cre(y) des(X) des(x) \n"

print "Test 1 output:"
print test1_string_output

if test1_string_output == test1_correct_answer:
    print "Test 1 passed",
else:
    print "Test 1 failed",

print "(%.3f seconds) \n" % (time.time() - startTime)


# Test 2: Construction of the overlap matrix for an M_01 sector of the effective Hamiltonian matrix
#
print ""
print "Starting test 2"
print ""
startTime = time.time()

# Create left hand side operator list
i = sqa.index('I', [tg_c])
x = sqa.index('X', [tg_a])

l_op  = [sqa.creOp(i), sqa.desOp(x)]

# Create right hand side operator list
j = sqa.index('J', [tg_c])
y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])
u = sqa.index('U', [tg_a])

r_op  = [sqa.creOp(z), sqa.creOp(u), sqa.desOp(y), sqa.desOp(j)]

# Define terms from operator lists
term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], r_op)

# Perform commutator
term3 = sqa.commutator(term1, term2)
term4 = sqa.matrixBlock(term3)

test2_string_output = ''
for t in term4:
    test2_string_output += str(t) + '\n'

test2_correct_answer = " (  -1.00000) kdelta(I,J) kdelta(U,X) cre(Z) des(Y) \n" + \
                       " (   1.00000) kdelta(I,J) kdelta(X,Z) cre(U) des(Y) \n" + \
                       " (  -1.00000) kdelta(I,J) cre(U) cre(Z) des(X) des(Y) \n"

print "Test 2 output:"
print test2_string_output

if test2_string_output == test2_correct_answer:
    print "Test 2 passed",
else:
    print "Test 2 failed",

print "(%.3f seconds) \n" % (time.time() - startTime)


# Test 3:
#
