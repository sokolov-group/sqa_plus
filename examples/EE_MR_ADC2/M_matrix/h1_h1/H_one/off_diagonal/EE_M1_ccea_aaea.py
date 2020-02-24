import sqa_extra.secondQuantizationAlgebra as sqa

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
z = sqa.index('Z', [tg_a])
b = sqa.index('B', [tg_v])
u = sqa.index('U', [tg_a])

l_op = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(x), sqa.desOp(a)]
r_op = [sqa.creOp(b), sqa.creOp(u), sqa.desOp(z), sqa.desOp(y)]

effH = []
effH = sqa.Heff(1)

for t in effH:
  print (t)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

print ("First Commutator")

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', 'IJAXYZBU', "")
