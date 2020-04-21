import sqa_extra.secondQuantizationAlgebra as sqa

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# External indices
j = sqa.index('J', [tg_c])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])

l_op = [sqa.creOp(j), sqa.creOp(x), sqa.desOp(z), sqa.desOp(y)]

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(1.0, [], l_op)

term2 = []

for t in effH:
    term2.append(sqa.multiplyTerms(term1, t))

term3 = []

for t in term2:
    tt = sqa.normalOrder(t)
    term3.extend(tt)

term4 = sqa.matrixBlock(term3)

sqa.generateEinsum(term4, '    temp', 'JXYZ', transRDM=True, trans_ind_str='I')
