import sqa_plus

# Define operator types
tg_c = sqa_plus.options.core_type
tg_a = sqa_plus.options.active_type
tg_v = sqa_plus.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# External indices
i = sqa_plus.index('I', [tg_c])
j = sqa_plus.index('J', [tg_c])
k = sqa_plus.index('K', [tg_c])
l = sqa_plus.index('L', [tg_c])

x = sqa_plus.index('X', [tg_a])
y = sqa_plus.index('Y', [tg_a])
z = sqa_plus.index('Z', [tg_a])
u = sqa_plus.index('U', [tg_a])

a = sqa_plus.index('A', [tg_v])
b = sqa_plus.index('B', [tg_v])
c = sqa_plus.index('C', [tg_v])
d = sqa_plus.index('D', [tg_v])

# Operations on Left and right reference wavefunctions.
# LHS
l_op  = [sqa_plus.creOp(i), sqa_plus.desOp(x)] # CA
l_ind = 'IX'                                                         

# RHS
r_op  = [sqa_plus.creOp(y), sqa_plus.desOp(j)] # CA
r_ind = 'JY'                                                         

term1 = sqa_plus.term(1.0, [], r_op)
term2 = sqa_plus.term(1.0, [], l_op)

print("Running 1st order terms")
term3 = []
term4 = []
V = sqa_plus.Vperturbation()

for v in V:
    term3.append(sqa_plus.multiplyTerms(term2, v))
for t in term3:
    term4.append(sqa_plus.multiplyTerms(t, term1))

Heff = sqa_plus.matrixBlock(term4)

sqa_plus.options.genEinsum.remove_core_integrals = False
sqa_plus.options.genEinsum.opt_einsum_terms = False
einsumlist = sqa_plus.genEinsum(Heff, 'temp', str(l_ind + r_ind), "")

## Read the reference file:
rfile = open("V_singleexc_test.ref","r")
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
 
