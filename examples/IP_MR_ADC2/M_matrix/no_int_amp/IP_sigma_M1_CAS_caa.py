import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True
i = sqa.index('i', [tg_c], dummy)
x = sqa.index('x', [tg_a], dummy)
y = sqa.index('y', [tg_a], dummy)
X = [sqa.tensor('X', [i,x,y], [])]

# External indices
j = sqa.index('J', [tg_c])
k = sqa.index('K', [tg_c])
l = sqa.index('L', [tg_c])
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
w = sqa.index('W', [tg_a])
z = sqa.index('Z', [tg_a])

r_op = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(y)]

effH = []
effH = sqa.Heff(1)
for t in effH:
  print t

term2 = sqa.term(1.0, [], X + r_op)

print "First Commutator"
term3 = sqa.commutator(effH,term2)
for t in term3:
    print t

term6 = sqa.matrixBlock(term3)

sqa.generateEinsum(term6, 'sigma[s_casci:f_casci]', '', transRDM=True, trans_ind_str="P")
#
#
#
##################################################################
