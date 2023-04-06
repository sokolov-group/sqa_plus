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
#l_op = [sqa.creOp(z), sqa.desOp(w), sqa.desOp(k)] # CAA (KWZ)
#X = [sqa.tensor('X', [k, w, z], [])]   # CAA
#prefactor = 1.0

#l_op = [sqa.creOp(w), sqa.desOp(l), sqa.desOp(k)] # CCA (KLW)
#X = [sqa.tensor('X', [k, l, w], Xsym)] # CCA
#prefactor = 0.5

#l_op = [sqa.creOp(b), sqa.desOp(l), sqa.desOp(k)] # CCE (KLB)
#X = [sqa.tensor('X', [k, l, b], Xsym)] # CCE
#prefactor = 0.5

#l_op = [sqa.creOp(b), sqa.desOp(z), sqa.desOp(w)] # AAE (WZB)
#X = [sqa.tensor('X', [w, z, b], Xsym)] # AAE
#prefactor = 0.5

l_op = [sqa.creOp(b), sqa.desOp(w), sqa.desOp(k)] # CAE (KWB)
X = [sqa.tensor('X', [k, w, b], [])]   # CAE
prefactor = 1.0

#
# R.h.s.
#
#r_op =  [sqa.creOp(i), sqa.creOp(x), sqa.desOp(y)] # CAA (IXY)
#ext_indices = "IXY"

#r_op =  [sqa.creOp(i), sqa.creOp(j), sqa.desOp(x)] # CCA (IJX)
#ext_indices = "IJX"

#r_op =  [sqa.creOp(i), sqa.creOp(j), sqa.desOp(a)] # CCE (IJA)
#ext_indices = "IJA"

r_op =  [sqa.creOp(x), sqa.creOp(y), sqa.desOp(a)] # AAE (XYA)
ext_indices = "XYA"

#r_op =  [sqa.creOp(i), sqa.creOp(x), sqa.desOp(a)] # CAE (IXA)
#ext_indices = "IXA"

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(prefactor, [], X + l_op)
term2 = sqa.term(1.0, [], r_op)

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
