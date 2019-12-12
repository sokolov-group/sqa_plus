import sqa_extra.secondQuantizationAlgebra as sqa

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True
#i = sqa.index('i', [tg_c], dummy)
#j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#l = sqa.index('l', [tg_c], dummy)
#a = sqa.index('aa', [tg_v], dummy)
#b = sqa.index('bb', [tg_v], dummy)
c = sqa.index('cc', [tg_v], dummy)
#d = sqa.index('dd', [tg_v], dummy)
x = sqa.index('x', [tg_a], dummy)
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#w = sqa.index('w', [tg_a], dummy)

X = [sqa.tensor('X', [x,c])]

# External indices
#x = sqa.index('X', [tg_a])
#c = sqa.index('C', [tg_v])

i = sqa.index('I', [tg_c])
j = sqa.index('J', [tg_c])
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])

l_op = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(b), sqa.desOp(a)]
r_op = [sqa.creOp(c), sqa.desOp(x)]

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(1.0, [], X + r_op)
term2 = sqa.term(1.0, [], l_op)

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)
term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, '    temp', 'IJAB', "")
