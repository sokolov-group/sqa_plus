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
# Test: Single Commutator Evaluation (EE M1 CAS Block)
# ============================================================================

@pytest.fixture
def test_m1_op():
    """
    RHS: cre(A) des(X)
    """
    tg_a = sqa.options.active_type
    tg_v = sqa.options.virtual_type

    x = sqa.index('X', [tg_a])
    a = sqa.index('A', [tg_v])

    term = sqa.term(1.0, [], [sqa.creOp(a), sqa.desOp(x)])

    return term

def test_ee_m1_cas_einsum(test_m1_op):

    term = test_m1_op

    # Define order of the effective Hamiltonian
    effH = sqa.Heff(1)

    # Perform commutator
    commutator = sqa.commutator(effH, term)
    result = sqa.matrixBlock(commutator)

    # Configure einsum options
    sqa.options.genEinsum.remove_core_integrals = False
    sqa.options.genEinsum.opt_einsum_terms = False
    sqa.options.genEinsum.trans_rdm = True
    sqa.options.genEinsum.indices_string = 'XA'
    sqa.options.genEinsum.trans_indices_string = 'I'

    # Generate einsum terms
    einsumlist = sqa.genEinsum(result)

    # Load reference file from same directory
    test_dir = os.path.dirname(os.path.abspath(__file__))
    ref_path = os.path.join(test_dir, 'EE_M1_CAS_ae_test.ref')

    # Compare with reference
    compare_einsum_output(einsumlist, ref_path)

