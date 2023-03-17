import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = False

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# Define the LHS and RHS subspaces
#lhs = 'ca'
#lhs = 'ce'
lhs = 'ae'

#rhs = 'ca'
#rhs = 'ce'
rhs = 'ae'

# Create LHS and external indices
if lhs == 'ca':
    i = sqa.index('I', [tg_c])
    x = sqa.index('X', [tg_a])

    l_op  = [sqa.creOp(i), sqa.desOp(x)] # CA
    l_ind = 'IX'                                                         

elif lhs == 'ce':
    i = sqa.index('I', [tg_c])
    a = sqa.index('A', [tg_v])

    l_op  = [sqa.creOp(i), sqa.desOp(a)] # CE
    l_ind = 'IA'                                                         

elif lhs == 'ae':
    x = sqa.index('X', [tg_a])
    a = sqa.index('A', [tg_v])

    l_op  = [sqa.creOp(x), sqa.desOp(a)] # AE
    l_ind = 'XA'                                                         

else:
    print 'Please define valid orbital subspaces for the LHS'

# Create RHS and external indices
if rhs == 'ca':
    j = sqa.index('J', [tg_c])
    y = sqa.index('Y', [tg_a])

    r_op  = [sqa.creOp(y), sqa.desOp(j)] # CA
    r_ind = 'JY'                                                         

elif rhs == 'ce':
    j = sqa.index('J', [tg_c])
    b = sqa.index('B', [tg_v])

    r_op  = [sqa.creOp(b), sqa.desOp(j)] # CE
    r_ind = 'JB'                                                         

elif rhs == 'ae':
    y = sqa.index('Y', [tg_a])
    b = sqa.index('B', [tg_v])

    r_op  = [sqa.creOp(b), sqa.desOp(y)] # AE
    r_ind = 'YB'                                                         

else:
    print 'Please define valid orbital subspaces for the RHS'

# Define order of the effective Hamiltonian
effH = []
#effH = sqa.Heff(1)
effH = sqa.Heff(2)

print ("Effective Hamiltonian")
for t in effH:
  print (t)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

print ("First Commutator")
term3 = sqa.commutator(effH, term1)

print ("Second Commutator")
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

# Convert result to einsum expressions
einsum_list = sqa.genEinsum(term5, lhs_str = 'temp', ind_str = l_ind + r_ind, trans_rdm = False, trans_ind_str = None, rm_core_int = True, intermediate_list = None, opt_einsum_terms = False, optimize = True)

# Print out terms in einsum format
for term in einsum_list:
   print (term)
