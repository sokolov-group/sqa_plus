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
# Authors: Alexander Yu. Sokolov <alexander.y.sokolov@gmail.com>
#          Koushik Chatterjee <koushikchatterjee7@gmail.com>
#          Ilia Mazin <ilia.mazin@gmail.com>
#          Carlos E. V. de Moura <carlosevmoura@gmail.com>
#          Donna Odhiambo <donna.odhiambo@proton.me>
#

import pytest
import sqa_plus as sqa

# Helper Function
def format_term_output(terms):
    """Format list of terms into a string output."""
    return ''.join(str(term) + ' \n' for term in terms)

# ============================================================================
# Test 1: Double Commutator Evaluation (M_00 Block)
# ============================================================================

@pytest.fixture
def test_m00_op():
    """
    LHS: cre(I) des(X)
    RHS: cre(Y) des(J)
    """
    tg_c = sqa.options.core_type
    tg_a = sqa.options.active_type

    i = sqa.index('I', [tg_c])
    x = sqa.index('X', [tg_a])
    j = sqa.index('J', [tg_c])
    y = sqa.index('Y', [tg_a])

    term_left  = sqa.term(1.0, [], [sqa.creOp(i), sqa.desOp(x)])
    term_right = sqa.term(1.0, [], [sqa.creOp(y), sqa.desOp(j)])

    return term_left, term_right

expected_output_m00 = (
    " (  -1.00000) v(I,Y,J,X) \n"
    " (   1.00000) kdelta(X,Y) v(I,x,J,y) rdm(x,y) \n"
    " (   1.00000) v(I,Y,J,x) cre(x) des(X) \n"
    " (   1.00000) v(I,x,J,X) cre(Y) des(x) \n"
    " (  -1.00000) kdelta(X,Y) v(I,x,J,y) cre(y) des(x) \n"
    " (  -1.00000) v(I,x,J,y) rdm(x,y) cre(Y) des(X) \n"
    " (  -1.00000) v(I,x,J,y) cre(Y) cre(y) des(X) des(x) \n"
)

def test_double_commutator_m00_block(test_m00_op):
    """
    Test 1: Double commutator evaluation for a first order contribution
    to the M_00 block of the effective Hamiltonian matrix.
    """
    term_left, term_right = test_m00_op

    # Define order of the effective Hamiltonian
    effH = sqa.Heff(1)

    # Perform first commutator
    inner_commutator = sqa.commutator(effH, term_right)

    # Perform second commutator
    outer_commutator = sqa.commutator(term_left, inner_commutator)

    result = sqa.matrixBlock(outer_commutator)
    output = format_term_output(result)

    assert output == expected_output_m00, (
        f"\nExpected:\n{expected_output_m00}\n\n"
        f"Got:\n{output}"
    )


# ============================================================================
# Test 2: Overlap Matrix (M_01 Block)
# ============================================================================

@pytest.fixture
def test_m01_op():
    """
    LHS: cre(I) des(X)
    RHS: cre(Z) cre(U) des(Y) des(J)
    """
    tg_c = sqa.options.core_type
    tg_a = sqa.options.active_type

    i = sqa.index('I', [tg_c])
    x = sqa.index('X', [tg_a])
    j = sqa.index('J', [tg_c])
    y = sqa.index('Y', [tg_a])
    z = sqa.index('Z', [tg_a])
    u = sqa.index('U', [tg_a])

    term_left  = sqa.term(1.0, [], [sqa.creOp(i), sqa.desOp(x)])
    term_right = sqa.term(1.0, [], [sqa.creOp(z), sqa.creOp(u), sqa.desOp(y), sqa.desOp(j)])

    return term_left, term_right

expected_output_m01 = (
    " (  -1.00000) kdelta(I,J) kdelta(U,X) cre(Z) des(Y) \n"
    " (   1.00000) kdelta(I,J) kdelta(X,Z) cre(U) des(Y) \n"
    " (  -1.00000) kdelta(I,J) cre(U) cre(Z) des(X) des(Y) \n"
)

def test_overlap_matrix_m01_sector(test_m01_op):
    """
    Test construction of the overlap matrix for a
    M_01 sector of the effective Hamiltonian matrix.
    """
    term_left, term_right = test_m01_op

    # Perform commutator
    result_commutator = sqa.commutator(term_left, term_right)

    result = sqa.matrixBlock(result_commutator)
    output = format_term_output(result)

    assert output == expected_output_m01, (
        f"\nExpected:\n{expected_output_m01}\n\n"
        f"Got:\n{output}"
    )

# ============================================================================
# Test 3: Operator Multiplication
# ============================================================================

@pytest.fixture
def test_V_op():
    """
    LHS: cre(X) des(A)
    RHS: cre(B) des(J)
    """
    tg_c = sqa.options.core_type
    tg_a = sqa.options.active_type
    tg_v = sqa.options.virtual_type

    x = sqa.index('X', [tg_a])
    a = sqa.index('A', [tg_v])
    j = sqa.index('J', [tg_c])
    b = sqa.index('B', [tg_v])

    term_left  = sqa.term(1.0, [], [sqa.creOp(x), sqa.desOp(a)])
    term_right = sqa.term(1.0, [], [sqa.creOp(b), sqa.desOp(j)])

    return term_left, term_right

expected_output_V = (
    " (  -1.00000) v(J,A,x,B) cre(X) des(x) \n"
    " (  -1.00000) kdelta(A,B) h(J,x) cre(X) des(x) \n"
    " (   1.00000) kdelta(A,B) v(J,i,i,x) cre(X) des(x) \n"
    " (   0.50000) kdelta(A,B) v(J,x,y,z) cre(X) cre(x) des(y) des(z) \n"
)

def test_perturbation_operator_multiplication(test_V_op):
    """
    Test multiplication of perturbation operator by
    single excitation operators from either side.
    """
    term_left, term_right = test_V_op

    V = sqa.Vperturbation()

    # Perform multiplication: LHS*V*RHS
    product = [
        sqa.multiplyTerms(sqa.multiplyTerms(term_left, term_v), term_right)
        for term_v in V
    ]

    result = sqa.matrixBlock(product)
    output = format_term_output(result)

    assert output == expected_output_V, (
        f"\nExpected:\n{expected_output_V}\n\n"
        f"Got:\n{output}"
    )
