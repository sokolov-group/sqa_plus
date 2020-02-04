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
effH = sqa.Heff(1)
for t in effH:
  print t
term2 = sqa.term(1.0, [], [sqa.creOp(j),sqa.creOp(k),sqa.desOp(x)])
print term2
#
print "First Commutator"
term3 = sqa.commutator(effH,term2)
for t in term3:
    print t
#
term6 = sqa.matrixBlock(term3)

print 'I, term2 =',term2
sqa.generateEinsum(term6, 'M[sI:fI, sxji:fxji]', 'XKJ', transRDM = True, trans_ind_str = 'I')
#
##################################################################
