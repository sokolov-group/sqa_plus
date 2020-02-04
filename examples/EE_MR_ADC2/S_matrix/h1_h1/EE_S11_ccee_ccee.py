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
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])

k = sqa.index('K', [tg_c])  
l = sqa.index('L', [tg_c])
c = sqa.index('C', [tg_v])
d = sqa.index('D', [tg_v])

l_op = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(b), sqa.desOp(a)]
r_op = [sqa.creOp(c), sqa.creOp(d), sqa.desOp(l), sqa.desOp(k)]


term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], r_op)

term3 = sqa.commutator(term1, term2)
term4 = sqa.matrixBlock(term3)

sqa.generateEinsum(term4, 'temp', 'IJABKLCD',"")
