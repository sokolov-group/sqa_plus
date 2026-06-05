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

##TODO: use indexList class for dummy indices

import pytest
import sqa_plus as sqa

# Helper Function
def format_term_output(terms):
    """Format list of terms into a string output."""
    return ''.join(str(term) + ' \n' for term in terms)

# ============================================================================
# Fixtures for Symmetries and Indices
# ============================================================================

@pytest.fixture(scope="session")
def symmetries():
    """Define tensor symmetries."""
    # Density matrix symmetries
    d1sym = [sqa.symmetry((1, 0), 1)]
    d2sym_aaaa = [sqa.symmetry((1, 0, 2, 3), -1), sqa.symmetry((2, 3, 0, 1), 1)]
    d2sym_abab = [sqa.symmetry((2, 3, 0, 1), 1)]
    d3sym_aaaaaa = [sqa.symmetry((1, 0, 2, 3, 4, 5), -1),
                    sqa.symmetry((0, 2, 1, 3, 4, 5), -1),
                    sqa.symmetry((3, 4, 5, 0, 1, 2), 1)]
    d3sym_aabaab = [sqa.symmetry((1, 0, 2, 3, 4, 5), -1),
                    sqa.symmetry((3, 4, 5, 0, 1, 2), 1)]

    # Amplitude symmetries
    a1sym = []
    a2sym_aaaa = [sqa.symmetry((1, 0, 2, 3), -1), sqa.symmetry((0, 1, 3, 2), -1)]
    a2sym_abab = []
    a2sym_sf = [sqa.symmetry((1, 0, 3, 2), 1)]

    # Integral symmetries
    h1sym = d1sym
    h2sym_aaaa = d2sym_aaaa
    h2sym_abab = d2sym_abab
    h2sym_sf = [sqa.symmetry((2, 3, 0, 1), 1), sqa.symmetry((1, 0, 3, 2), 1)]

    return {
        'd1sym': d1sym,
        'd2sym_aaaa': d2sym_aaaa,
        'd2sym_abab': d2sym_abab,
        'd3sym_aaaaaa': d3sym_aaaaaa,
        'd3sym_aabaab': d3sym_aabaab,
        'a1sym': a1sym,
        'a2sym_aaaa': a2sym_aaaa,
        'a2sym_abab': a2sym_abab,
        'a2sym_sf': a2sym_sf,
        'h1sym': h1sym,
        'h2sym_aaaa': h2sym_aaaa,
        'h2sym_abab': h2sym_abab,
        'h2sym_sf': h2sym_sf,
    }

@pytest.fixture(scope="session")
def dummy_indices():
    """Define dummy indices for testing."""
    return [sqa.index('i%i' % i, [], True, False) for i in range(20)]

@pytest.fixture(scope="session")
def fixed_indices():
    """Define non-dummy indices for testing."""
    return [sqa.index('j%i' % i, [], False) for i in range(20)]

# ============================================================================
# Test 1: Normal Ordering
# ============================================================================

expected_output_test1 = (
    " (  -1.00000) kdelta(j1,j3) kdelta(j2,j4) \n"
    " (   1.00000) kdelta(j1,j4) kdelta(j2,j3) \n"
    " (   1.00000) kdelta(j1,j3) cre(j4) des(j2) \n"
    " (  -1.00000) kdelta(j1,j4) cre(j3) des(j2) \n"
    " (  -1.00000) kdelta(j2,j3) cre(j4) des(j1) \n"
    " (   1.00000) kdelta(j2,j4) cre(j3) des(j1) \n"
    " (   1.00000) cre(j3) cre(j4) des(j1) des(j2) \n"
)

def test_normal_ordering(fixed_indices):
    """
    Test 1: Normal ordering of a(j1) a(j2) a+(j3) a+(j4).
    """
    test1_tensors = []
    test1_tensors.append(sqa.desOp(fixed_indices[1]))
    test1_tensors.append(sqa.desOp(fixed_indices[2]))
    test1_tensors.append(sqa.creOp(fixed_indices[3]))
    test1_tensors.append(sqa.creOp(fixed_indices[4]))
    test1_term = sqa.term(1.0, [], test1_tensors)

    test1_output = sqa.normalOrder(test1_term)
    sqa.combineTerms(test1_output)

    output = format_term_output(test1_output)

    assert output == expected_output_test1, (
        f"\nExpected:\n{expected_output_test1}\n\n"
        f"Got:\n{output}"
    )


# ============================================================================
# Test 2: Commutator with Index Symmetries
# ============================================================================

expected_output_test2 = (
    " (   4.00000) a2(a,b,c,d) h1(a,e) cre(c) cre(d) des(b) des(e) \n"
    " (  -4.00000) a2(a,b,c,d) h1(c,e) cre(d) cre(e) des(a) des(b) \n"
)

def test_commutator_1e_hamiltonian_with_2e_amplitude(dummy_indices, symmetries):
    """
    Test 2: Commutator of all alpha 1e hamiltonian with all alpha 2 electron amplitude.
    The result is transpose equated, so that symmetrization is required in the resulting Hamiltonian.
    Without index types.
    """
    h1 = sqa.tensor('h1', dummy_indices[0:2], symmetries['h1sym'])
    a2 = sqa.tensor('a2', dummy_indices[2:6], symmetries['a2sym_aaaa'])
    h1_term = sqa.term(1.0, [], [h1, sqa.creOp(dummy_indices[0]), sqa.desOp(dummy_indices[1])])

    a2_terms = []
    a2_terms.append(sqa.term( 1.0, [], [a2, sqa.creOp(dummy_indices[2]),
                                            sqa.creOp(dummy_indices[3]),
                                            sqa.desOp(dummy_indices[5]),
                                            sqa.desOp(dummy_indices[4])]))

    a2_terms.append(sqa.term(-1.0, [], [a2, sqa.creOp(dummy_indices[4]),
                                             sqa.creOp(dummy_indices[5]),
                                             sqa.desOp(dummy_indices[3]),
                                             sqa.desOp(dummy_indices[2])]))

    test2_output = sqa.commutator(h1_term, a2_terms)
    sqa.combine_transpose(test2_output)

    output = format_term_output(test2_output)

    assert output == expected_output_test2, (
        f"\nExpected:\n{expected_output_test2}\n\n"
        f"Got:\n{output}"
    )


# ============================================================================
# Test 3: Commutator with Index Types
# ============================================================================

expected_output_test3 = (
    " (   2.00000) a2_abab(a,b,c,d) h2_aaaa(a,e,f,g) cre(b) cre(f) cre(g) des(c) des(d) des(e) \n"
    "index a types = (('alpha',), ('active', 'virtual'))\n"
    "index b types = (('beta',), ('active', 'virtual'))\n"
    "index c types = (('alpha',), ('active', 'core'))\n"
    "index d types = (('beta',), ('active', 'core'))\n"
    "index e types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index f types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index g types = (('alpha',), ('active', 'core', 'virtual'))\n"
    " (   2.00000) a2_abab(a,b,c,d) h2_aaaa(a,e,f,g) cre(c) cre(d) cre(e) des(b) des(f) des(g) \n"
    "index a types = (('alpha',), ('active', 'virtual'))\n"
    "index b types = (('beta',), ('active', 'virtual'))\n"
    "index c types = (('alpha',), ('active', 'core'))\n"
    "index d types = (('beta',), ('active', 'core'))\n"
    "index e types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index f types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index g types = (('alpha',), ('active', 'core', 'virtual'))\n"
    " (  -2.00000) a2_abab(a,b,c,d) h2_aaaa(c,e,f,g) cre(a) cre(b) cre(e) des(d) des(f) des(g) \n"
    "index a types = (('alpha',), ('active', 'virtual'))\n"
    "index b types = (('beta',), ('active', 'virtual'))\n"
    "index c types = (('alpha',), ('active', 'core'))\n"
    "index d types = (('beta',), ('active', 'core'))\n"
    "index e types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index f types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index g types = (('alpha',), ('active', 'core', 'virtual'))\n"
    " (  -2.00000) a2_abab(a,b,c,d) h2_aaaa(c,e,f,g) cre(d) cre(f) cre(g) des(a) des(b) des(e) \n"
    "index a types = (('alpha',), ('active', 'virtual'))\n"
    "index b types = (('beta',), ('active', 'virtual'))\n"
    "index c types = (('alpha',), ('active', 'core'))\n"
    "index d types = (('beta',), ('active', 'core'))\n"
    "index e types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index f types = (('alpha',), ('active', 'core', 'virtual'))\n"
    "index g types = (('alpha',), ('active', 'core', 'virtual'))\n"
)

@pytest.fixture
def test3_setup():
    test3_indices = []
    test3_indices.append(sqa.index('i0', [['alpha'], ['core', 'active', 'virtual']], True, False))
    test3_indices.append(sqa.index('i1', [['alpha'], ['core', 'active', 'virtual']], True, False))
    test3_indices.append(sqa.index('i2', [['alpha'], ['core', 'active', 'virtual']], True, False))
    test3_indices.append(sqa.index('i3', [['alpha'], ['core', 'active', 'virtual']], True, False))
    test3_indices.append(sqa.index('i4', [['alpha'], ['active', 'virtual']], True, False))
    test3_indices.append(sqa.index('i5', [['beta'],  ['active', 'virtual']], True, False))
    test3_indices.append(sqa.index('i6', [['alpha'], ['core', 'active']], True, False))
    test3_indices.append(sqa.index('i7', [['beta'],  ['core', 'active']], True, False))
    return test3_indices

def test_commutator_2e_hamiltonian_abab_amplitude(test3_setup, symmetries):
    """
    Test 3: Commutator of all alpha 2e hamiltonian with alpha/beta 2 electron amplitude.
    With index types.
    """
    test3_indices = test3_setup

    h2_aaaa = sqa.tensor('h2_aaaa', test3_indices[0:4], symmetries['h2sym_aaaa'])
    a2_abab = sqa.tensor('a2_abab', test3_indices[4:8], symmetries['a2sym_abab'])

    h2_aaaa_term = sqa.term(1.0, [], [h2_aaaa, sqa.creOp(test3_indices[0]),
                                            sqa.creOp(test3_indices[1]),
                                            sqa.desOp(test3_indices[3]),
                                            sqa.desOp(test3_indices[2])])

    a2_abab_terms = []
    a2_abab_terms.append(sqa.term( 1.0, [], [a2_abab, sqa.creOp(test3_indices[4]),
                                                           sqa.creOp(test3_indices[5]),
                                                           sqa.desOp(test3_indices[7]),
                                                           sqa.desOp(test3_indices[6])]))

    a2_abab_terms.append(sqa.term(-1.0, [], [a2_abab, sqa.creOp(test3_indices[6]),
                                                           sqa.creOp(test3_indices[7]),
                                                           sqa.desOp(test3_indices[5]),
                                                           sqa.desOp(test3_indices[4])]))

    test3_output = sqa.commutator(h2_aaaa_term, a2_abab_terms)

    test3_string_output = ""
    for t in test3_output:
        test3_string_output += str(t) + ' \n'
        test3_ind_list = []
        for ten in t.tensors:
            for ind in ten.indices:
                if not (ind in test3_ind_list):
                    test3_ind_list.append(ind.copy())
                    test3_string_output += 'index %s types = %s\n' % (ind.name, str(ind.indType))

    assert test3_string_output == expected_output_test3, (
        f"\nExpected:\n{expected_output_test3}\n\n"
        f"Got:\n{test3_string_output}"
    )


# ============================================================================
# Test 4: Normal Ordering of Spin-Free Operators
# ============================================================================

expected_output_test4 = (
    " (   1.00000) kdelta(i1,i2) kdelta(i4,i6) kdelta(i5,i7) E2(i0,i3,i8,i9) \n"
    " (   1.00000) kdelta(i1,i2) kdelta(i4,i7) kdelta(i5,i6) E2(i0,i3,i9,i8) \n"
    " (   1.00000) kdelta(i1,i3) kdelta(i4,i6) kdelta(i5,i7) E2(i0,i2,i9,i8) \n"
    " (   1.00000) kdelta(i1,i3) kdelta(i4,i7) kdelta(i5,i6) E2(i0,i2,i8,i9) \n"
    " (   1.00000) kdelta(i1,i2) kdelta(i4,i6) E3(i0,i3,i7,i8,i5,i9) \n"
    " (   1.00000) kdelta(i1,i2) kdelta(i4,i7) E3(i0,i3,i6,i9,i5,i8) \n"
    " (   1.00000) kdelta(i1,i2) kdelta(i5,i6) E3(i0,i3,i7,i4,i8,i9) \n"
    " (   1.00000) kdelta(i1,i2) kdelta(i5,i7) E3(i0,i3,i6,i4,i9,i8) \n"
    " (   1.00000) kdelta(i1,i3) kdelta(i4,i6) E3(i0,i2,i7,i5,i8,i9) \n"
    " (   1.00000) kdelta(i1,i3) kdelta(i4,i7) E3(i0,i2,i6,i5,i9,i8) \n"
    " (   1.00000) kdelta(i1,i3) kdelta(i5,i6) E3(i0,i2,i7,i8,i4,i9) \n"
    " (   1.00000) kdelta(i1,i3) kdelta(i5,i7) E3(i0,i2,i6,i9,i4,i8) \n"
    " (   1.00000) kdelta(i1,i6) kdelta(i4,i7) E3(i0,i2,i3,i8,i9,i5) \n"
    " (   1.00000) kdelta(i1,i6) kdelta(i5,i7) E3(i0,i2,i3,i8,i4,i9) \n"
    " (   1.00000) kdelta(i1,i7) kdelta(i4,i6) E3(i0,i2,i3,i9,i8,i5) \n"
    " (   1.00000) kdelta(i1,i7) kdelta(i5,i6) E3(i0,i2,i3,i9,i4,i8) \n"
    " (   1.00000) kdelta(i4,i6) kdelta(i5,i7) E3(i0,i2,i3,i1,i8,i9) \n"
    " (   1.00000) kdelta(i4,i7) kdelta(i5,i6) E3(i0,i2,i3,i1,i9,i8) \n"
    " (   1.00000) kdelta(i1,i2) E4(i0,i3,i6,i7,i4,i5,i8,i9) \n"
    " (   1.00000) kdelta(i1,i3) E4(i0,i2,i6,i7,i5,i4,i8,i9) \n"
    " (   1.00000) kdelta(i1,i6) E4(i0,i2,i3,i7,i8,i4,i5,i9) \n"
    " (   1.00000) kdelta(i1,i7) E4(i0,i2,i3,i6,i9,i4,i5,i8) \n"
    " (   1.00000) kdelta(i4,i6) E4(i0,i2,i3,i7,i1,i8,i5,i9) \n"
    " (   1.00000) kdelta(i4,i7) E4(i0,i2,i3,i6,i1,i9,i5,i8) \n"
    " (   1.00000) kdelta(i5,i6) E4(i0,i2,i3,i7,i1,i4,i8,i9) \n"
    " (   1.00000) kdelta(i5,i7) E4(i0,i2,i3,i6,i1,i4,i9,i8) \n"
    " (   1.00000) E5(i0,i2,i3,i6,i7,i1,i4,i5,i8,i9) \n"
)

def test_normal_ordering_spin_free_operators():
    """
    Test 4: Normal ordering of the spin free operator string E(i0,i1) E(i2,i3,i4,i5) E(i6,i7,i8,i9).
    """
    test4_indices = [sqa.index('i%i' % i) for i in range(10)]
    test4_term = sqa.term(1.0, [], [sqa.sfExOp(test4_indices[0:2]),
                                         sqa.sfExOp(test4_indices[2:6]),
                                         sqa.sfExOp(test4_indices[6:10])])

    test4_output = sqa.normalOrder(test4_term)
    sqa.combineTerms(test4_output)

    output = format_term_output(test4_output)

    assert output == expected_output_test4, (
        f"\nExpected:\n{expected_output_test4}\n\n"
        f"Got:\n{output}"
    )

# ============================================================================
# Test 5: Commutator of Spin-Free Hamiltonian and Amplitudes
# ============================================================================

expected_output_test5 = (
    " (   2.00000) a2(a,b,c,d) h2(a,b,e,f) E2(c,d,e,f) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    " (   2.00000) a2(a,b,c,d) h2(a,b,e,f) E2(e,f,c,d) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    " (  -2.00000) a2(a,b,c,d) h2(c,d,e,f) E2(a,b,e,f) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    " (  -2.00000) a2(a,b,c,d) h2(c,d,e,f) E2(e,f,a,b) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    " (   4.00000) a2(a,b,c,d) h2(a,e,f,g) E3(b,f,g,d,c,e) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    "index g types = (('active', 'core', 'virtual'),)\n"
    " (   4.00000) a2(a,b,c,d) h2(a,e,f,g) E3(c,d,e,f,b,g) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    "index g types = (('active', 'core', 'virtual'),)\n"
    " (  -4.00000) a2(a,b,c,d) h2(c,e,f,g) E3(a,b,e,f,d,g) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    "index g types = (('active', 'core', 'virtual'),)\n"
    " (  -4.00000) a2(a,b,c,d) h2(c,e,f,g) E3(d,f,g,b,a,e) \n"
    "index a types = (('active', 'virtual'),)\n"
    "index b types = (('active', 'virtual'),)\n"
    "index c types = (('active', 'core'),)\n"
    "index d types = (('active', 'core'),)\n"
    "index e types = (('active', 'core', 'virtual'),)\n"
    "index f types = (('active', 'core', 'virtual'),)\n"
    "index g types = (('active', 'core', 'virtual'),)\n"
)

@pytest.fixture
def test5_setup():
    test5_indices = []
    test5_indices.append(sqa.index('i0', [['core', 'active', 'virtual']], True, False))
    test5_indices.append(sqa.index('i1', [['core', 'active', 'virtual']], True, False))
    test5_indices.append(sqa.index('i2', [['core', 'active', 'virtual']], True, False))
    test5_indices.append(sqa.index('i3', [['core', 'active', 'virtual']], True, False))
    test5_indices.append(sqa.index('i4', [['active', 'virtual']], True, False))
    test5_indices.append(sqa.index('i5', [['active', 'virtual']], True, False))
    test5_indices.append(sqa.index('i6', [['core', 'active']], True, False))
    test5_indices.append(sqa.index('i7', [['core', 'active']], True, False))
    return test5_indices

def test_commutator_spin_free_hamiltonian_amplitudes(test5_setup, symmetries):
    """
    Test 5: Commutator of 2e spin free hamiltonian with 2e spin free amplitudes.
    """
    test5_indices = test5_setup

    test5_h2 = sqa.tensor('h2', test5_indices[0:4], symmetries['h2sym_sf'])
    test5_hterm = sqa.term(1.0, [], [test5_h2, sqa.sfExOp(test5_indices[0:4])])

    test5_a2 = sqa.tensor('a2', test5_indices[4:8], symmetries['a2sym_sf'])
    test5_aterms = []
    test5_aterms.append(sqa.term( 1.0, [], [test5_a2, sqa.sfExOp(test5_indices[4:6] + test5_indices[6:8])]))
    test5_aterms.append(sqa.term(-1.0, [], [test5_a2, sqa.sfExOp(test5_indices[6:8] + test5_indices[4:6])]))

    test5_output = sqa.commutator(test5_hterm, test5_aterms)

    test5_string_output = ""
    for t in test5_output:
        test5_string_output += str(t) + ' \n'
        test5_ind_list = []
        for ten in t.tensors:
            for ind in ten.indices:
                if not (ind in test5_ind_list):
                    test5_ind_list.append(ind.copy())
                    test5_string_output += 'index %s types = %s\n' % (ind.name, str(ind.indType))

    assert test5_string_output == expected_output_test5, (
        f"\nExpected:\n{expected_output_test5}\n\n"
        f"Got:\n{test5_string_output}"
    )

# ============================================================================
# Test 6: Core Operator Pair Removal
# ============================================================================

# old
#expected_output_test6 = (
#    " (   1.00000) \n"
#    " (  -1.00000) cre(c0) \n"
#    " (   1.00000) cre(a0) r0(a0,a1) cre(c1) des(a1) \n"
#    " (   1.00000) cre(a0) des(a1) \n"
#    " (  -1.00000) cre(c0) cre(c0) des(c0) des(a1) \n"
#    " (  -1.00000) cre(a0) cre(c1) des(c1) des(c1) \n"
#)

expected_output_test6 = (
    " (   1.00000) \n"
    " (  -1.00000) cre(c0) \n"
    " (  -1.00000) cre(a0) r0(a0,a1) cre(c1) des(a1) \n"
    " (   1.00000) cre(a0) des(a1) \n"
    " (  -1.00000) cre(c0) cre(c0) des(c0) des(a1) \n"
    " (   1.00000) cre(a0) cre(c1) des(c1) des(c1) \n"
)


@pytest.fixture
def test6_setup():
    ta = sqa.options.alpha_type
    tc = sqa.options.core_type
    tt = sqa.options.active_type

    # Define core indices
    ac = [sqa.index('c%i' % i, [ta, tc], False) for i in range(10)]

    # Define active indices
    at = [sqa.index('a%i' % i, [ta, tt], True) for i in range(10)]

    return ac, at

def test_remove_core_operator_pairs(test6_setup):
    """
    Test 6: Test of removeCoreOpPairs function.
    """
    ac, at = test6_setup

    test6_terms = []
    test6_terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.desOp(ac[0])]))

    test6_terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.creOp(ac[1]),
                                               sqa.desOp(ac[1])]))

    test6_terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[0]),
                                               sqa.tensor('r0', [at[0], at[1]], []),
                                               sqa.creOp(ac[1]),
                                               sqa.desOp(ac[0]), sqa.desOp(at[1])]))

    test6_terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[0]),
                                               sqa.creOp(ac[1]),
                                               sqa.desOp(ac[0]), sqa.desOp(at[1]),
                                               sqa.desOp(ac[1])]))

    test6_terms.append(sqa.term(1.0, [], [sqa.creOp(ac[0]), sqa.creOp(ac[1]),
                                               sqa.creOp(ac[0]),
                                               sqa.desOp(ac[0]), sqa.desOp(at[1]),
                                               sqa.desOp(ac[1])]))

    test6_terms.append(sqa.term(1.0, [], [sqa.creOp(at[0]), sqa.creOp(ac[1]),
                                               sqa.creOp(ac[0]),
                                               sqa.desOp(ac[1]), sqa.desOp(ac[0]),
                                               sqa.desOp(ac[1])]))

    sqa.removeCoreOpPairs(test6_terms)

    output = format_term_output(test6_terms)

    assert output == expected_output_test6, (
        f"\nExpected:\n{expected_output_test6}\n\n"
        f"Got:\n{output}"
    )
