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
l = sqa.index('l', [tg_c], dummy)
m = sqa.index('m', [tg_c], dummy)
u = sqa.index('u', [tg_a], dummy)
v = sqa.index('v', [tg_a], dummy)
Xsym_1 = sqa.symmetry((1,0,2,3),-1)
Xsym_2 = sqa.symmetry((0,1,3,2),-1)
X = [sqa.tensor('X', [l,m,u,v], [Xsym_1, Xsym_2])]
l_op  = [sqa.creOp(l), sqa.creOp(m), sqa.desOp(v), sqa.desOp(u)] # CCAA  
term1 = sqa.term(0.25, [], X + l_op)


#l = sqa.index('l', [tg_c], dummy)
#m = sqa.index('m', [tg_c], dummy)
#d = sqa.index('dd', [tg_v], dummy)
#u = sqa.index('u', [tg_a], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#X = [sqa.tensor('X', [l,m,d,u], [Xsym_1])]
#l_op  = [sqa.creOp(l), sqa.creOp(m), sqa.desOp(u), sqa.desOp(d)] # CCEA
#term1 = sqa.term(0.50, [], X + l_op)


#l = sqa.index('l', [tg_c], dummy)
#m = sqa.index('m', [tg_c], dummy)
#d = sqa.index('dd', [tg_v], dummy)
#e = sqa.index('ee', [tg_v], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [l,m,d,e], [Xsym_1, Xsym_2])]
#l_op  = [sqa.creOp(l), sqa.creOp(m), sqa.desOp(e), sqa.desOp(d)] # CCEE
#term1 = sqa.term(0.25, [], X + l_op)


#l = sqa.index('l', [tg_c], dummy)
#u = sqa.index('u', [tg_a], dummy)
#d = sqa.index('dd', [tg_v], dummy)
#e = sqa.index('ee', [tg_v], dummy)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [l,u,d,e], [Xsym_2])]
#l_op  = [sqa.creOp(l), sqa.creOp(u), sqa.desOp(e), sqa.desOp(d)] # CAEE
#term1 = sqa.term(0.50, [], X + l_op)


#u = sqa.index('u', [tg_a], dummy)
#v = sqa.index('v', [tg_a], dummy)
#d = sqa.index('dd', [tg_v], dummy)
#e = sqa.index('ee', [tg_v], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [u,v,d,e], [Xsym_1, Xsym_2])]
#l_op  = [sqa.creOp(u), sqa.creOp(v), sqa.desOp(e), sqa.desOp(d)] # AAEE
#term1 = sqa.term(0.25, [], X + l_op)


#l = sqa.index('l', [tg_c], dummy)
#u = sqa.index('u', [tg_a], dummy)
#v = sqa.index('v', [tg_a], dummy)
#w = sqa.index('w', [tg_a], dummy)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [l,u,v,w], [Xsym_2])]
#l_op  = [sqa.creOp(l), sqa.creOp(u), sqa.desOp(w), sqa.desOp(v)] # CAAA
#term1 = sqa.term(0.50, [], X + l_op)


#l = sqa.index('l', [tg_c], dummy)
#u = sqa.index('u', [tg_a], dummy)
#d = sqa.index('dd', [tg_v], dummy)
#v = sqa.index('v', [tg_a], dummy)
#X = [sqa.tensor('X', [l,u,d,v])]
#l_op  = [sqa.creOp(l), sqa.creOp(u), sqa.desOp(v), sqa.desOp(d)] # CAEA
#term1 = sqa.term(1.0, [], X + l_op)


#u = sqa.index('u', [tg_a], dummy)
#v = sqa.index('v', [tg_a], dummy)
#d = sqa.index('dd', [tg_v], dummy)
#w = sqa.index('w', [tg_a], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#X = [sqa.tensor('X', [u,v,d,w], [Xsym_1])]
#l_op  = [sqa.creOp(u), sqa.creOp(v), sqa.desOp(w), sqa.desOp(d)] # AAEA
#term1 = sqa.term(0.50, [], X + l_op)



# RHS
#r_op  = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(j), sqa.desOp(i)] # CCAA
#r_ind = 'IJXY'                                                         

#r_op  = [sqa.creOp(a), sqa.creOp(x), sqa.desOp(j), sqa.desOp(i)] # CCEA
#r_ind = 'IJAX'                                                         

#r_op  = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(j), sqa.desOp(i)] # CCEE
#r_ind = 'IJAB'                                                         

#r_op  = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(x), sqa.desOp(i)] # CAEE
#r_ind = 'IXAB'                                                         

#r_op  = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(y), sqa.desOp(x)] # AAEE
#r_ind = 'XYAB'                                                         

r_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(x), sqa.desOp(i)] # CAAA
r_ind = 'IXYZ'                                                         

#r_op  = [sqa.creOp(a), sqa.creOp(y), sqa.desOp(x), sqa.desOp(i)] # CAEA
#r_ind = 'IXAY'                                                         

#r_op  = [sqa.creOp(a), sqa.creOp(z), sqa.desOp(y), sqa.desOp(x)] # AAEA
#r_ind = 'XYAZ'


# Define Hamiltonian
effH = []
effH = sqa.Heff(1)

term2 = sqa.term(1.0, [], r_op)

print("First Commutator")

term3 = sqa.commutator(effH, term2)
term4 = sqa.commutator(term1, term3)

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', r_ind, "")
