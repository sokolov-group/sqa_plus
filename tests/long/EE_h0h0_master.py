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
i = sqa.index('I', [tg_c])
j = sqa.index('J', [tg_c])
k = sqa.index('K', [tg_c])
l = sqa.index('L', [tg_c])

x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])
u = sqa.index('U', [tg_a])

a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
c = sqa.index('C', [tg_v])
d = sqa.index('D', [tg_v])


# LHS
l_op  = [sqa.creOp(i), sqa.desOp(x)] # CA
l_ind = 'IX'                                                         

#l_op  = [sqa.creOp(i), sqa.desOp(a)] # CE
#l_ind = 'IA'                                                         

#l_op  = [sqa.creOp(x), sqa.desOp(a)] # AE
#l_ind = 'XA'                                                         


# RHS
#r_op  = [sqa.creOp(y), sqa.desOp(j)] # CA
#r_ind = 'JY'                                                         

#r_op  = [sqa.creOp(b), sqa.desOp(j)] # CE
#r_ind = 'JB'                                                         

r_op  = [sqa.creOp(b), sqa.desOp(y)] # AE
r_ind = 'YB'                                                         


# Define Hamiltonian
effH = []
effH = sqa.Heff(2)


term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)


term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

einsumlist = sqa.generateEinsum(term5, 'temp', str(l_ind + r_ind), "")

## Read the reference file:

rfile = open("EE_h0h0_master.ref","r")
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
 
