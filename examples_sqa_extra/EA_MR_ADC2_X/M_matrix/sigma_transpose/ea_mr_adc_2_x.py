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
j = sqa.index('jj', [tg_c], dummy)
c = sqa.index('cc', [tg_v], dummy)
d = sqa.index('dd', [tg_v], dummy)
w = sqa.index('ww', [tg_a], dummy)
z = sqa.index('zz', [tg_a], dummy)

#
# L.h.s.
#
#l_op = [sqa.creOp(w), sqa.desOp(z), sqa.desOp(c)] # WZC (AAE)
#X = [sqa.tensor('X', [w, z, c], [])]
#prefactor = 1.0

#l_op = [sqa.creOp(j), sqa.desOp(c), sqa.desOp(d)] # JCD (CEE)
#Xsym = [sqa.symmetry((0,2,1),-1)]
#X = [sqa.tensor('X', [j, c, d], Xsym)]
#prefactor = 0.5

l_op = [sqa.creOp(j), sqa.desOp(w), sqa.desOp(c)] # JWC (CAE)
X = [sqa.tensor('X', [j, w, c], [])]
prefactor = 1.0

#l_op = [sqa.creOp(w), sqa.desOp(c), sqa.desOp(d)] # WCD (AEE)
#Xsym = [sqa.symmetry((0,2,1),-1)]
#X = [sqa.tensor('X', [w, c, d], Xsym)]
#prefactor = 0.5

#l_op = [sqa.creOp(j), sqa.desOp(w), sqa.desOp(z)] # JWZ (CAA)
#Xsym = [sqa.symmetry((0,2,1),-1)]
#X = [sqa.tensor('X', [j, w, z], Xsym)]
#prefactor = 0.5

#
# R.h.s.
#
#r_op = [sqa.creOp(a), sqa.creOp(y), sqa.desOp(x)] # XYA (AAE)
#ext_indices = "XYA"

#r_op = [sqa.creOp(b), sqa.creOp(a), sqa.desOp(i)] # IAB (CEE)
#ext_indices = "IAB"

#r_op = [sqa.creOp(a), sqa.creOp(x), sqa.desOp(i)] # IXA (CAE)
#ext_indices = "IXA"

#r_op = [sqa.creOp(b), sqa.creOp(a), sqa.desOp(x)] # XAB (AEE)
#ext_indices = "XAB"

r_op = [sqa.creOp(y), sqa.creOp(x), sqa.desOp(i)] # IXY (CAA)
ext_indices = "IXY"

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(prefactor, [], X + l_op)
term2 = sqa.term(1.0, [], r_op)

print "First Commutator"
term3 = sqa.commutator(effH,term2)

print "Second Commutator"
term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(term1,t))
    term4.append(sqa.multiplyTerms(t,term1))

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', ext_indices)

