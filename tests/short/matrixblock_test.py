import sqa_plus

# Second order terms
print("Running 2nd order terms")
term3 = []
term4 = []
term  = []   

V = sqa_plus.Vperturbation()
T = sqa_plus.Tamplitude(1)

for t_v in V:
    for t_t in T:
        term3.append(sqa_plus.multiplyTerms(t_v, t_t))

Heff = sqa_plus.matrixBlock(term3)

sqa_plus.options.genEinsum.remove_core_integrals = False
sqa_plus.options.genEinsum.opt_einsum_terms = False
sqa_plus.options.genEinsum.trans_rdm = True
einsumlist = sqa_plus.genEinsum(Heff, 'temp', '', trans_indices_string="P")

## Read the reference file:
rfile = open("matrixblock_test.ref","r")
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
 
