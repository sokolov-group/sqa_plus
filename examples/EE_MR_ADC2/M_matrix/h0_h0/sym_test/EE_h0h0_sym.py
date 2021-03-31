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


l_op_term = sqa.term(1.0, [], l_op)
r_op_term = sqa.term(1.0, [], r_op)

# Define Hamiltonian
effH = []
effH = sqa.Heff(2)

term1 = []
for t in effH:
    tt = sqa.multiplyTerms(l_op_term, t)
    term1.append(tt)

term2 = []
for t in term1:
    tt = sqa.multiplyTerms(t, r_op_term)
    term2.append(tt)

sym_term = []
for t in term2:
    tt = sqa.normalOrder(t)
    sym_term.extend(tt)

term4 = sqa.matrixBlock(sym_term)

einsum_list = sqa.genEinsum(term4, 'temp', l_ind + r_ind, rm_core_int = True)
for einsum in einsum_list:
    print (einsum)

