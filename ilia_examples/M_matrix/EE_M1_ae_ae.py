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
x = sqa.index('I', [tg_a])
a = sqa.index('X', [tg_v])
y = sqa.index('J', [tg_a])
b = sqa.index('Y', [tg_v])

r_op = [sqa.creOp(b), sqa.desOp(y)]
l_op = [sqa.creOp(x), sqa.desOp(a)]

effH = []
effH = sqa.Heff(1)

for t in effH:
  print (t)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

print ("First Commutator")

term3 = sqa.commutator(effH, term1)
for t in term3:
    print (t)

term4 = sqa.commutator(term2, term3)
for t in term4:
    print (t)

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'M[s_ae:f_ae, s_ae:f_ae]', 'XAYB',"")
