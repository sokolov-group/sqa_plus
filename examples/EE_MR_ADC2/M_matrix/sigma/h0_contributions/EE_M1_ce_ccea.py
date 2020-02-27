import sqa_extra.secondQuantizationAlgebra as sqa

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True
#i = sqa.index('i', [tg_c], dummy)
j = sqa.index('j', [tg_c], dummy)
k = sqa.index('k', [tg_c], dummy)
#l = sqa.index('l', [tg_c], dummy)
#a = sqa.index('aa', [tg_v], dummy)
b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#d = sqa.index('dd', [tg_v], dummy)
x = sqa.index('x', [tg_a], dummy)
#y = sqa.index('y', [tg_a], dummy)
#w = sqa.index('w', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)

# External indices
i = sqa.index('I', [tg_c])
a = sqa.index('A', [tg_v])

#j = sqa.index('J', [tg_c])
#k = sqa.index('K', [tg_c])
#b = sqa.index('B', [tg_v])
#x = sqa.index('X', [tg_a])

Xsym_1 = sqa.symmetry((1,0,2,3),-1)

X = [sqa.tensor('X', [j,k,b,x], [Xsym_1])]

l_op = [sqa.creOp(i), sqa.desOp(a)]
r_op = [sqa.creOp(b), sqa.creOp(x), sqa.desOp(k), sqa.desOp(j)]

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(0.5, [], X + r_op)
term2 = sqa.term(1.0, [], l_op)

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)
term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, '    temp', 'IA', "")