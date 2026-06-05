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
# Test: Single Commutator Evaluation (EE h1h1 Overlap Block)
# ============================================================================
@pytest.fixture
def test_h1h1_overlap_op():
    """
    LHS: cre(X) cre(Y) des(Z) des(A)
    RHS: cre(D) des(W) des(V) des(U)
    """
    tg_a = sqa.options.active_type
    tg_v = sqa.options.virtual_type

    x = sqa.index('X', [tg_a])
    y = sqa.index('Y', [tg_a])
    z = sqa.index('Z', [tg_a])
    u = sqa.index('U', [tg_a])
    v = sqa.index('V', [tg_a])
    w = sqa.index('W', [tg_a])
    a = sqa.index('A', [tg_v])
    d = sqa.index('D', [tg_v])

    term_left  = sqa.term(1.0, [], [sqa.creOp(x), sqa.creOp(y), sqa.desOp(z), sqa.desOp(a)])
    term_right = sqa.term(1.0, [], [sqa.creOp(d), sqa.creOp(w), sqa.desOp(v), sqa.desOp(u)])

    return term_left, term_right

def test_ee_h1h1_overlap_einsum(test_h1h1_overlap_op):
    term_left, term_right = test_h1h1_overlap_op

    # Perform commutaEE_h1h1_overlap_test.pytor
    commutator = sqa.commutator(term_left, term_right)

    result = sqa.matrixBlock(commutator)

    # Configure einsum options
    sqa.options.genEinsum.remove_core_integrals = False
    sqa.options.genEinsum.opt_einsum_terms = False
    sqa.options.genEinsum.indices_string = 'XYAZUVDW'

    # Generate einsum terms
    einsumlist = sqa.genEinsum(result)

    # Load reference file from same directory
    test_dir = os.path.dirname(os.path.abspath(__file__))
    ref_path = os.path.join(test_dir, 'EE_h1h1_overlap_test.ref')

    # Compare with reference
    compare_einsum_output(einsumlist, ref_path)
