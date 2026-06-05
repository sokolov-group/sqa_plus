import os
import pytest
import sqa_plus as sqa

# Helper Function
def compare_einsum_output(test_einsum, ref_path):
    """Compare test einsum terms with reference terms."""
    with open(ref_path) as f:
        ref_einsum =  [line.rstrip("\n") for line in f]

    test_einsum = [line.rstrip("\n") for line in test_einsum]

    if len(ref_einsum) != len(test_einsum):
        raise AssertionError(
            f"Number of einsum terms do not match: expected {len(ref_einsum)}, "
            f"got {len(test_einsum)}."
        )

    differences = [
        f"Line {i}:\n  Expected:  {ref}\n  Got:       {test}"
        for i, (ref, test) in enumerate(zip(ref_einsum, test_einsum), start=1)
        if ref != test
    ]

    if differences:
        raise AssertionError("Differences found in einsum terms:\n" + "\n".join(differences))

# ============================================================================
# Test: Double Commutator Evaluation (EE M_00 Block)
# ============================================================================

@pytest.fixture
def test_m00_op():
    """
    LHS: cre(I) des(X)
    RHS: cre(B) des(Y)
    """
    tg_c = sqa.options.core_type
    tg_a = sqa.options.active_type
    tg_v = sqa.options.virtual_type

    i = sqa.index('I', [tg_c])
    x = sqa.index('X', [tg_a])
    y = sqa.index('Y', [tg_a])
    b = sqa.index('B', [tg_v])

    term_left = sqa.term(1.0, [], [sqa.creOp(i), sqa.desOp(x)])
    term_right = sqa.term(1.0, [], [sqa.creOp(b), sqa.desOp(y)])

    return term_left, term_right

def test_ee_h0h0_einsum(test_m00_op):
    """
    Test: Double commutator evaluation for a second order contribution
    to the M_00 block of the effective Hamiltonian matrix.
    """
    term_left, term_right = test_m00_op

    # Define order of the effective Hamiltonian
    effH = sqa.Heff(2)

    # Perform first commutator
    inner_commutator = sqa.commutator(effH, term_right)

    # Perform second commutator
    outer_commutator = sqa.commutator(term_left, inner_commutator)

    result = sqa.matrixBlock(outer_commutator)

    # Configure einsum options
    sqa.options.genEinsum.remove_core_integrals = False
    sqa.options.genEinsum.opt_einsum_terms = False
    sqa.options.genEinsum.indices_string = "IXYB"

    # Generate einsum terms
    einsumlist = sqa.genEinsum(result)

    # Load reference file from same directory
    test_dir = os.path.dirname(os.path.abspath(__file__))
    ref_path = os.path.join(test_dir, 'EE_h0h0_test.ref')

    # Compare with reference
    compare_einsum_output(einsumlist, ref_path)

