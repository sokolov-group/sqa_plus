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

i_op = [sqa.desOp(i)]
j_op = [sqa.creOp(j)]

effH = []
effH = sqa.Heff(1)
for t in effH:
  print t

term1 = sqa.term(1.0, [], i_op)
term2 = sqa.term(1.0, [], j_op)

print "First Commutator"
term3 = sqa.commutator(effH,term2)
for t in term3:
    print t

term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(t,term1))

print "Second Commutator"
for t in term4:
    print t

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'M[s_c:f_c, s_c:f_c]', 'IJ', "")
#
#
#
##################################################################
