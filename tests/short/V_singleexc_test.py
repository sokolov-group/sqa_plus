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
# Test: First Order Perturbation (V Single Excitation Block)
# ============================================================================
@pytest.fixture
def test_v_singleexc_op():
    """
    LHS: cre(I) des(X)
    RHS: cre(Y) des(J)
    """
    tg_c = sqa.options.core_type
    tg_a = sqa.options.active_type

    i = sqa.index('I', [tg_c])
    j = sqa.index('J', [tg_c])
    x = sqa.index('X', [tg_a])
    y = sqa.index('Y', [tg_a])

    term_left  = sqa.term(1.0, [], [sqa.creOp(i), sqa.desOp(x)])
    term_right = sqa.term(1.0, [], [sqa.creOp(y), sqa.desOp(j)])

    return term_left, term_right

def test_v_singleexc_einsum(test_v_singleexc_op):

    term_left, term_right = test_v_singleexc_op

    # Define perturbation operators
    V = sqa.Vperturbation()

    # Compute LHS * V * RHS
    pdt = [
        sqa.multiplyTerms(sqa.multiplyTerms(term_left, term_v), term_right)
        for term_v in V
    ]

    result = sqa.matrixBlock(pdt)

    # Configure einsum options
    sqa.options.genEinsum.remove_core_integrals = False
    sqa.options.genEinsum.opt_einsum_terms = False
    sqa.options.genEinsum.indices_string = 'IXJY'

    # Generate einsum terms
    einsumlist = sqa.genEinsum(result)

    # Load reference file from same directory
    test_dir = os.path.dirname(os.path.abspath(__file__))
    ref_path = os.path.join(test_dir, 'V_singleexc_test.ref')

    # Compare with reference
    compare_einsum_output(einsumlist, ref_path)

