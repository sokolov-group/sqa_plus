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
i = sqa.index('I', [tg_c])
a = sqa.index('A', [tg_v])

j = sqa.index('J', [tg_c])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])

l_op = [sqa.creOp(i), sqa.desOp(a)]
r_op = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(x), sqa.desOp(j)]

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)
term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', 'IAJXYZ', "")
