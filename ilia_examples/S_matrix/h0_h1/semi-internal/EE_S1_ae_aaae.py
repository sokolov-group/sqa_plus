import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# External indices
x = sqa.index('X', [tg_a])
a = sqa.index('A', [tg_v])

y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])
b = sqa.index('B', [tg_v])
w = sqa.index('W', [tg_a])

l_op = [sqa.creOp(x), sqa.desOp(a)]
r_op = [sqa.creOp(b), sqa.creOp(w), sqa.desOp(z), sqa.desOp(y)]

term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], r_op)

term3 = sqa.commutator(term1, term2)
term4 = sqa.matrixBlock(term3)

sqa.generateEinsum(term4, '        temp', 'XAYZBW', "")
