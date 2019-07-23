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
y = sqa.index('Y', [tg_a])
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])

z = sqa.index('Z', [tg_a])  
w = sqa.index('W', [tg_a])
c = sqa.index('C', [tg_v])
d = sqa.index('D', [tg_v])

l_op = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(y), sqa.desOp(x)]
r_op = [sqa.creOp(z), sqa.creOp(w), sqa.desOp(d), sqa.desOp(c)]


term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], r_op)

term3 = sqa.commutator(term1, term2)
term4 = sqa.matrixBlock(term3)

sqa.generateEinsum(term4, 'temp', 'XYABZWCD',"")
