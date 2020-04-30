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
a = sqa.index('A', [tg_v])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])

# Dummy indices
k = sqa.index('kk', [tg_c], dummy)
l = sqa.index('ll', [tg_c], dummy)
b = sqa.index('bb', [tg_v], dummy)
w = sqa.index('ww', [tg_a], dummy)
z = sqa.index('zz', [tg_a], dummy)

Xsym = [sqa.symmetry((1, 0, 2), -1)]

#
# L.h.s.
#
#l_op = [sqa.creOp(y), sqa.desOp(x), sqa.desOp(i)] # CAA (IXY)
#ext_indices = "IXY"

#l_op = [sqa.creOp(x), sqa.desOp(j), sqa.desOp(i)] # CCA (IJX)
#ext_indices = "IJX"

#l_op = [sqa.creOp(a), sqa.desOp(j), sqa.desOp(i)] # CCE (IJA)
#ext_indices = "IJA"

#l_op = [sqa.creOp(a), sqa.desOp(y), sqa.desOp(x)] # AAE (XYA)
#ext_indices = "XYA"

l_op = [sqa.creOp(a), sqa.desOp(x), sqa.desOp(i)] # CAE (IXA)
ext_indices = "IXA"

#
# R.h.s.
#
#r_op = [sqa.creOp(k), sqa.creOp(w), sqa.desOp(z)] # CAA (KWZ)
#X = [sqa.tensor('X', [k, w, z], [])]   # CAA
#prefactor = 1.0

#r_op = [sqa.creOp(k), sqa.creOp(l), sqa.desOp(w)] # CCA (KLW)
#X = [sqa.tensor('X', [k, l, w], Xsym)] # CCA
#prefactor = 0.5

#r_op = [sqa.creOp(k), sqa.creOp(l), sqa.desOp(b)] # CCE (KLB)
#X = [sqa.tensor('X', [k, l, b], Xsym)] # CCE
#prefactor = 0.5

#r_op = [sqa.creOp(w), sqa.creOp(z), sqa.desOp(b)] # AAE (WZB)
#X = [sqa.tensor('X', [w, z, b], Xsym)] # AAE
#prefactor = 0.5

r_op = [sqa.creOp(k), sqa.creOp(w), sqa.desOp(b)] # CAE (KWB)
X = [sqa.tensor('X', [k, w, b], [])]   # CAE
prefactor = 1.0

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(prefactor, [], X + r_op)

print "First Commutator"
term3 = sqa.commutator(effH, term2)
for t in term3:
    print t

term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(t,term1))
    term4.append(sqa.multiplyTerms(term1,t))

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', ext_indices)
#
#
#
##################################################################
