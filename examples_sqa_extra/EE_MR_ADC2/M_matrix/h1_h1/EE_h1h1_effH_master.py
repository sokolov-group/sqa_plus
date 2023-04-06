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
m = sqa.index('M', [tg_c])
n = sqa.index('N', [tg_c])

x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])
u = sqa.index('U', [tg_a])
v = sqa.index('V', [tg_a])
w = sqa.index('W', [tg_a])

a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
c = sqa.index('C', [tg_v])
d = sqa.index('D', [tg_v])
e = sqa.index('E', [tg_v])
f = sqa.index('F', [tg_v])


# LHS
#l_op  = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(y), sqa.desOp(x)] # CCAA
#l_ind = 'IJXY'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(x), sqa.desOp(a)] # CCEA
#l_ind = 'IJAX'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(b), sqa.desOp(a)] # CCEE
#l_ind = 'IJAB'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(b), sqa.desOp(a)] # CAEE
#l_ind = 'IXAB'                                                         
                                                                       
#l_op  = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(b), sqa.desOp(a)] # AAEE
#l_ind = 'XYAB'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(z), sqa.desOp(y)] # CAAA
#l_ind = 'IXYZ'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(y), sqa.desOp(a)] # CAEA
#l_ind = 'IXAY'                                                         
                                                                       
l_op  = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(z), sqa.desOp(a)] # AAEA
l_ind = 'XYAZ'


# RHS
#r_op  = [sqa.creOp(u), sqa.creOp(v), sqa.desOp(m), sqa.desOp(l)] # CCAA  
#r_ind = 'LMUV'

#r_op  = [sqa.creOp(d), sqa.creOp(u), sqa.desOp(m), sqa.desOp(l)] # CCEA
#r_ind = 'LMDU'

#r_op  = [sqa.creOp(d), sqa.creOp(e), sqa.desOp(m), sqa.desOp(l)] # CCEE
#r_ind = 'LMDE'

#r_op  = [sqa.creOp(d), sqa.creOp(e), sqa.desOp(u), sqa.desOp(l)] # CAEE
#r_ind = 'LUDE'

#r_op  = [sqa.creOp(d), sqa.creOp(e), sqa.desOp(v), sqa.desOp(u)] # AAEE
#r_ind = 'UVDE'

#r_op  = [sqa.creOp(v), sqa.creOp(w), sqa.desOp(u), sqa.desOp(l)] # CAAA
#r_ind = 'LUVW'

#r_op  = [sqa.creOp(d), sqa.creOp(v), sqa.desOp(u), sqa.desOp(l)] # CAEA
#r_ind = 'LUDV'

r_op  = [sqa.creOp(d), sqa.creOp(w), sqa.desOp(v), sqa.desOp(u)] # AAEA
r_ind = 'UVDW'

# Define Hamiltonian
effH = []
effH = sqa.Heff(0)

for t in effH:
  print(t)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

print("First Commutator")

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', str(l_ind + r_ind), "")
