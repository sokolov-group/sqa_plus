import sqa_plus

# Define operator types
tg_c = sqa_plus.options.core_type
tg_a = sqa_plus.options.active_type
tg_v = sqa_plus.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# External indices
x = sqa_plus.index('X', [tg_a])
a = sqa_plus.index('A', [tg_v])

r_op = [sqa_plus.creOp(a), sqa_plus.desOp(x)]

effH = []
effH = sqa_plus.Heff(1)

term1 = sqa_plus.term(1.0, [], r_op)

term2 = sqa_plus.commutator(effH, term1)

term3 = sqa_plus.matrixBlock(term2)

sqa_plus.options.genEinsum.remove_core_integrals = False
sqa_plus.options.genEinsum.opt_einsum_terms = False
sqa_plus.options.genEinsum.trans_rdm = True
einsumlist = sqa_plus.genEinsum(term3, 'temp', 'XA', trans_indices_string="I")

## Read the reference file:
rfile = open("EE_M1_CAS_ae_test.ref","r")
ref_file = rfile.readlines()

def compare_ref(test_einsum, ref_einsum):

    print("\nCompare and Test with respect to the Reference set")
  
    test_einsum = [ i.strip("\n") for i in test_einsum ]
    ref_einsum   = [ i.strip("\n") for i in ref_einsum ]
 
    check = True
    if(len(ref_einsum) != len(test_einsum)):

        print("Error: The number of einsum terms don't match with the reference")
        check = False

    else:

        for l in range(len(ref_einsum)):
            if(ref_einsum[l] != test_einsum[l]):
                check = False
                print("Error: Differences in einsum terms in line: "+str(l+1)) 
                print("Test:      " + test_einsum[l])
                print("Reference: " + ref_einsum[l] + "\n")
                
    if(check == False):
        print("Test Failed")
    else: 
        print("Test Passed")
  
    return check

run_test = compare_ref(einsumlist, ref_file)
            
