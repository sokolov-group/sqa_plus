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
m = sqa_plus.index('M', [tg_c])
n = sqa_plus.index('N', [tg_c])

x = sqa_plus.index('X', [tg_a])
y = sqa_plus.index('Y', [tg_a])
z = sqa_plus.index('Z', [tg_a])
u = sqa_plus.index('U', [tg_a])
v = sqa_plus.index('V', [tg_a])
w = sqa_plus.index('W', [tg_a])

a = sqa_plus.index('A', [tg_v])
b = sqa_plus.index('B', [tg_v])
c = sqa_plus.index('C', [tg_v])
d = sqa_plus.index('D', [tg_v])
e = sqa_plus.index('E', [tg_v])
f = sqa_plus.index('F', [tg_v])

# LHS
l_op  = [sqa_plus.creOp(x), sqa_plus.creOp(y), sqa_plus.desOp(z), sqa_plus.desOp(a)] # AAEA
l_ind = 'XYAZ'

# RHS
r_op  = [sqa_plus.creOp(d), sqa_plus.creOp(w), sqa_plus.desOp(v), sqa_plus.desOp(u)] # AAEA
r_ind = 'UVDW'

term1 = sqa_plus.term(1.0, [], l_op)
term2 = sqa_plus.term(1.0, [], r_op)

term3 = sqa_plus.commutator(term1, term2)
term4 = sqa_plus.matrixBlock(term3)

sqa_plus.options.genEinsum.remove_core_integrals = False
sqa_plus.options.genEinsum.opt_einsum_terms = False
einsumlist = sqa_plus.genEinsum(term4, 'temp', str(l_ind + r_ind), "")

## Read the reference file:
rfile = open("EE_h1h1_overlap_test.ref","r")
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
 
