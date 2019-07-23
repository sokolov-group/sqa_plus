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
j = sqa.index('J', [tg_c])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])

k = sqa.index('K', [tg_c])  
l = sqa.index('L', [tg_c])
z = sqa.index('Z', [tg_a])
w = sqa.index('W', [tg_a])

l_op = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(j), sqa.desOp(i)]
r_op = [sqa.creOp(k), sqa.creOp(l), sqa.desOp(w), sqa.desOp(z)]


term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], r_op)

term3 = sqa.commutator(term1, term2)
term4 = sqa.matrixBlock(term3)

sqa.generateEinsum(term4, 'temp', 'IJXYKLZW',"")
