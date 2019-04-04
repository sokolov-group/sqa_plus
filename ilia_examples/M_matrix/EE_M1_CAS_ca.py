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
x = sqa.index('X', [tg_a])

r_op = [sqa.creOp(x), sqa.desOp(i)]

effH = []
effH = sqa.Heff(1)

for t in effH:
  print (t)

term1 = sqa.term(1.0, [], r_op)

print ("First Commutator")

term2 = sqa.commutator(effH, term1)
for t in term2:
    print (t)

term3 = sqa.matrixBlock(term2)

sqa.generateEinsum(term3, 'M[s_casci:f_casci, s_ca:f_ca]', 'JX', transRDM=True, trans_ind_str="I")
