import sqa_extra.secondQuantizationAlgebra as sqa

#sqa.options.verbose = True

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
k = sqa.index('K', [tg_c])
l = sqa.index('L', [tg_c])
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
w = sqa.index('W', [tg_a])
z = sqa.index('Z', [tg_a])

l_op = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(b), sqa.desOp(a)]

effH = []
effH = sqa.Heff(2)
for t in effH:
  print t

term1 = sqa.term(1.0, [], l_op)

term2 = []
for t in effH:
    term2.append(sqa.multiplyTerms(term1,t))

term3 = []
for t in term2:
    t_no = sqa.normalOrder(t)
    term3.extend(t_no)
    print t

term4 = sqa.matrixBlock(term3)

sqa.generateEinsum(term4, 'temp', 'XYAB', "")
