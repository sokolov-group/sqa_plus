import os
import pytest
import sqa_plus as sqa

# Helper Function
def compare_einsum_output(test_einsum, ref_path):
    """Compare test einsum terms with reference terms."""
    with open(ref_path) as f:
        ref_einsum = [line.rstrip("\n") for line in f]
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
# Test: Double Commutator Evaluation (EE h0h1 Block)
# ============================================================================

@pytest.fixture
def test_h0h1_op():
    """
    LHS: cre(X) des(A)
    RHS: cre(B) cre(U) des(Z) des(Y)
    """
    tg_a = sqa.options.active_type
    tg_v = sqa.options.virtual_type

    x = sqa.index('X', [tg_a])
    y = sqa.index('Y', [tg_a])
    z = sqa.index('Z', [tg_a])
    u = sqa.index('U', [tg_a])
    a = sqa.index('A', [tg_v])
    b = sqa.index('B', [tg_v])

    term_left  = sqa.term(1.0, [], [sqa.creOp(x), sqa.desOp(a)])
    term_right = sqa.term(1.0, [], [sqa.creOp(b), sqa.creOp(u), sqa.desOp(z), sqa.desOp(y)])

    return term_left, term_right


def test_ee_h0h1_effH_einsum(test_h0h1_op):

    term_left, term_right = test_h0h1_op

    # Define order of the effective Hamiltonian
    effH = sqa.Heff(0)

    # Perform first commutator
    inner_commutator = sqa.commutator(effH, term_right)

    # Perform second commutator
    outer_commutator = sqa.commutator(term_left, inner_commutator)

    result = sqa.matrixBlock(outer_commutator)

    # Configure einsum options
    sqa.options.genEinsum.remove_core_integrals = False
    sqa.options.genEinsum.opt_einsum_terms = False
    sqa.options.genEinsum.indices_string = 'XAYZBU'

    # Generate einsum terms
    einsumlist = sqa.genEinsum(result)

    # Load reference file from same directory
    test_dir = os.path.dirname(os.path.abspath(__file__))
    ref_path = os.path.join(test_dir, 'EE_h0h1_effH_test.ref')

    # Compare with reference
    compare_einsum_output(einsumlist, ref_path)

