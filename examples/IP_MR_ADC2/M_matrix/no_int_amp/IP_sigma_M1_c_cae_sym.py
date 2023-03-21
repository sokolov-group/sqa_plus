import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True
i = sqa.index('i', [tg_c], dummy)
#j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#l = sqa.index('l', [tg_c], dummy)
a = sqa.index('aa', [tg_v], dummy)
#b = sqa.index('bb', [tg_v], dummy)
x = sqa.index('x', [tg_a], dummy)
#y = sqa.index('y', [tg_a], dummy)
#w = sqa.index('w', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
X = [sqa.tensor('X', [i,x,a], [])]

# External indices
#i = sqa.index('I', [tg_c])
#j = sqa.index('J', [tg_c])
k = sqa.index('K', [tg_c])
#l = sqa.index('L', [tg_c])
#a = sqa.index('A', [tg_v])
#b = sqa.index('B', [tg_v])
#x = sqa.index('X', [tg_a])
#y = sqa.index('Y', [tg_a])
#w = sqa.index('W', [tg_a])
#z = sqa.index('Z', [tg_a])

r_op = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(a)]
l_op = [sqa.desOp(k)]

effH = []
effH = sqa.Heff(1)
for t in effH:
  print t

term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], X + r_op)

term3 = []
for t in effH:
    term3.append(sqa.multiplyTerms(term2, t))

term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(t,term1))

term5 = []
for t in term4:
    t.scale(-1.0)
    t_no = sqa.normalOrder(t)
    term5.extend(t_no)
    print t

term6 = sqa.matrixBlock(term5)

sqa.generateEinsum(term6, 'sigma[s_c:f_c]', 'K')
#
#
#
##################################################################
