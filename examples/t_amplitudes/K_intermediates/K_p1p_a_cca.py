import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# External indices
u = sqa.index('U', [tg_a])
v = sqa.index('V', [tg_a])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
w = sqa.index('W', [tg_a])
z = sqa.index('Z', [tg_a])

l_op = [sqa.desOp(x)]
r_op = [sqa.creOp(y), sqa.creOp(w), sqa.desOp(z)]

effH = []
effH = sqa.Heff(0)
for t in effH:
  print t

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

print "First Commutator"
term3 = sqa.commutator(effH,term1)
for t in term3:
    print t

term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(term2,t))

term5 = []
print "Second Commutator"
for t in term4:
    t_no = sqa.normalOrder(t)
    term5.extend(t_no)
    print t

term6 = sqa.matrixBlock(term5)

sqa.generateEinsum(term6, 'K12', 'XYWZ')
#
#
#
##################################################################
