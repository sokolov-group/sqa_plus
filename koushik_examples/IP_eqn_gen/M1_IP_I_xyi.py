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
k = sqa.index('K', [tg_c])
a = sqa.index('A', [tg_v])

x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])

#
effH = []
effH = sqa.Heff(0)
for t in effH:
  print t
term1 = sqa.term(1.0, [], [sqa.desOp(i)])
term2 = sqa.term(1.0, [], [sqa.creOp(j),sqa.creOp(y),sqa.desOp(x)])
print term1, term2
#
print "First Commutator"
term3 = sqa.commutator(effH,term2)
for t in term3:
    print t
#
print ''
#
term4 = []
for t in term3:
  term4.append(sqa.multiplyTerms(term1,t))
  term4.append(sqa.multiplyTerms(t,term1))

term5 = []
for t in term4:
  t_no = sqa.normalOrder(t)
  term5.extend(t_no)
  print t
#
for t in term5:
 print 'normalorder = ',
#
#sqa.combineTerms(term2)
#
term6 = sqa.matrixBlock(term5)

sqa.generateEinsum(term6, 'M[si:fi, sxyi:fxyi]', 'IXYJ')
#
##################################################################
