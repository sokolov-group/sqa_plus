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
# Test: Second Order Terms
# ============================================================================

def test_T_einsum():

    # Define perturbation operators
    V_1 = sqa.Vperturbation()
    V_2 = sqa.Vperturbation()

    # Define amplitude operators
    T_1 = sqa.Tamplitude(1)
    T_2 = sqa.Tamplitude(1)

    # Compute 1/2 * V * T
    term_1 = [sqa.multiplyTerms(t_v, t_t) for t_v in V_1 for t_t in T_1]
    for t in term_1:
        t.scale(0.5)

    term_2 = [sqa.multiplyTerms(t_t, t_v) for t_v in V_2 for t_t in T_2]
    for t in term_2:
        t.scale(-0.5)

    # Add two terms together
    term = term_1 + term_2

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
    ref_path = os.path.join(test_dir, 'T_extend_test.ref')

    # Compare with reference
    compare_einsum_output(einsumlist, ref_path)

