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
k = sqa.index('K', [tg_c])
l = sqa.index('L', [tg_c])
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])

# Dummy indices
j = sqa.index('j', [tg_c], dummy)
c = sqa.index('c', [tg_v], dummy)
d = sqa.index('d', [tg_v], dummy)
w = sqa.index('w', [tg_a], dummy)
z = sqa.index('z', [tg_a], dummy)

#
# L.h.s.
#
l_op = [sqa.creOp(x), sqa.desOp(y), sqa.desOp(a)] # XYA
ext_indices = "XYA"

#l_op = [sqa.creOp(i), sqa.desOp(x), sqa.desOp(a)] # IXA
#ext_indices = "IXA"

#l_op = [sqa.creOp(i), sqa.desOp(x), sqa.desOp(y)] # IXY
#ext_indices = "IXY"

#l_op = [sqa.creOp(x), sqa.desOp(a), sqa.desOp(b)] # XAB
#ext_indices = "XAB"

#
# R.h.s.
#
r_op = [sqa.creOp(c), sqa.creOp(z), sqa.desOp(w)] # WZC
X = [sqa.tensor('X', [w, z, c], [])]
prefactor = 1.0

#r_op = [sqa.creOp(c), sqa.creOp(w), sqa.desOp(j)] # JWC
#X = [sqa.tensor('X', [j, w, c], [])]
#prefactor = 1.0

#r_op = [sqa.creOp(z), sqa.creOp(w), sqa.desOp(j)] # JWZ
#Xsym = [sqa.symmetry((0,2,1),-1)]
#X = [sqa.tensor('X', [j, w, z], Xsym)]
#prefactor = 0.5

#r_op = [sqa.creOp(d), sqa.creOp(c), sqa.desOp(w)] # WCD
#Xsym = [sqa.symmetry((0,2,1),-1)]
#X = [sqa.tensor('X', [w, c, d], Xsym)]
#prefactor = 0.5


effH = []
effH = sqa.Heff(0)

term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(prefactor, [], X + r_op)

print "First Commutator"
term3 = sqa.commutator(effH,term2)

print "Second Commutator"
term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(term1,t))
    term4.append(sqa.multiplyTerms(t,term1))

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', ext_indices)

