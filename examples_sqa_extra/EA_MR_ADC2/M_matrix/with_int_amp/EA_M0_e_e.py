import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = False

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# External indices
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])

a_op = [sqa.desOp(a)]
b_op = [sqa.creOp(b)]

effH = []
effH = sqa.Heff(0)

term1 = sqa.term(1.0, [], a_op)
term2 = sqa.term(1.0, [], b_op)

print "First Commutator"
term3 = sqa.commutator(effH,term2)

print "Second Commutator"
term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(term1,t))
    term4.append(sqa.multiplyTerms(t,term1))

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'M[s_e:f_e, s_e:f_e]', 'AB')

#
#
#
##################################################################
