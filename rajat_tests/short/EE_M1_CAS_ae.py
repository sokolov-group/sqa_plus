import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = False

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# External indices
x = sqa.index('X', [tg_a])
a = sqa.index('A', [tg_v])

r_op = [sqa.creOp(a), sqa.desOp(x)]

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(1.0, [], r_op)

term2 = sqa.commutator(effH, term1)

term3 = sqa.matrixBlock(term2)

einsumlist = sqa.generateEinsum(term3, 'temp', 'XA', transRDM=True, trans_ind_str="I")


## Read the reference file:

rfile = open("EE_M1_CAS_ae.ref","r")
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
            
