import sqa_extra.secondQuantizationAlgebra as sqa

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# External indices
i = sqa.index('I', [tg_c])
x = sqa.index('X', [tg_a])
a = sqa.index('A', [tg_v])
y = sqa.index('Y', [tg_a])

j = sqa.index('J', [tg_c])
z = sqa.index('Z', [tg_a])
b = sqa.index('B', [tg_v])
w = sqa.index('W', [tg_a])

X = [sqa.tensor('X', [j,z,b,w])]

l_op = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(y), sqa.desOp(a)]
r_op = [sqa.creOp(b), sqa.creOp(w), sqa.desOp(z), sqa.desOp(j)]

effH = []
effH = sqa.Heff(0)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, '    temp', 'IXAYJZBW', "")
