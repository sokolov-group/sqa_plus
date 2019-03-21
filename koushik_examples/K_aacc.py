import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

z = sqa.index('Z', [tg_a])
x = sqa.index('X', [tg_a])
w = sqa.index('W', [tg_a])
y = sqa.index('Y', [tg_a])

zw_op = [sqa.desOp(z), sqa.desOp(w)]
xy_op = [sqa.creOp(x), sqa.creOp(y)]

effH = []
effH = sqa.Heff(0)

term1 = sqa.commutator(effH, sqa.term(1.0, [], xy_op))
print "First Commutator"
for t in term1:
    print t

term2 = []
print "Multiply by the left operator"
for t in term1:
    op_t = sqa.multiplyTerms(sqa.term(1.0, [], zw_op), t)
    op_t_no = sqa.normalOrder(op_t)
    term2.extend(op_t_no)

for t in term2:
    t.contractDeltaFuncs()
#    print 'term2 = ', t
sqa.combineTerms(term2)
for t in term2:
    print t

sqa.dummyLabel(term2)
sqa.generateEinsum(term2, "K_aacc_so", "ZWXY")
exit()

term3 = sqa.matrixBlock(term2)
for t in term3:
    print 'term3 = ',t
#sqa.generateEinsum(term3, "K_aacc_so", "ZWXY")
