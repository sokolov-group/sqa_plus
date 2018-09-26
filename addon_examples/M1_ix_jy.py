import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

## Core dummy indices
#cc = [sqa.index('c%i' %p, [tg_c], dummy) for p in range(10)]
## Active dummy indices
#aa = [sqa.index('a%i' %p, [tg_a], dummy) for p in range(10)]
## Virtual dummy indices
#vv = [sqa.index('v%i' %p, [tg_v], dummy) for p in range(10)]

# External indices
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
i = sqa.index('I', [tg_c])
j = sqa.index('J', [tg_c])

# Operators with external indices
# Note that we define i^\dag as q-annihilation operator with index i, etc
ix_op = [sqa.creOp(i), sqa.desOp(x)]
yj_op = [sqa.creOp(y), sqa.desOp(j)]
#
effH = []
effH = sqa.Heff(1)
#for t in effH:
#  print t
#
term1 = sqa.commutator(effH,sqa.term(1.0, [], yj_op))
print "First Commutator"
for t in term1:
    print t
#
term2 = sqa.commutator(sqa.term(1.0, [], ix_op),term1)
print "Second Commutator"
for t in term2:
    print t
#
#sqa.combineTerms(term2)
#
sqa.addon(term2)
#
exit()
#
##################################################################
