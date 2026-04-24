import os
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
# Test: Matrix Block
# ============================================================================

def test_matrixblock_einsum():

    # Define perturbation operators
    V = sqa.Vperturbation()

    # Define amplitude operators
    T = sqa.Tamplitude(1)

    # Compute T * V
    term = [sqa.multiplyTerms(t_v, t_t) for t_v in V for t_t in T]

    Heff = sqa.matrixBlock(term)

    # Configure einsum options
    sqa.options.genEinsum.remove_core_integrals = False
    sqa.options.genEinsum.opt_einsum_terms = False
    sqa.options.genEinsum.trans_rdm = True
    sqa.options.genEinsum.trans_indices_string = "P"

    # Generate einsum terms
    einsumlist = sqa.genEinsum(Heff)

    # Load reference file from same directory
    test_dir = os.path.dirname(os.path.abspath(__file__))
    ref_path = os.path.join(test_dir, 'matrixblock_test.ref')

    # Compare with reference
    compare_einsum_output(einsumlist, ref_path)

