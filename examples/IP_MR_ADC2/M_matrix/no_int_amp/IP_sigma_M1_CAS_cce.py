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
j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#l = sqa.index('l', [tg_c], dummy)
a = sqa.index('aa', [tg_v], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#x = sqa.index('x', [tg_a], dummy)
#y = sqa.index('y', [tg_a], dummy)
#w = sqa.index('w', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
Xsym = [sqa.symmetry((1,0,2),-1)]
#Xsym = []
X = [sqa.tensor('X', [i,j,a], Xsym)]

# External indices
#i = sqa.index('I', [tg_c])
#j = sqa.index('J', [tg_c])
#k = sqa.index('K', [tg_c])
#l = sqa.index('L', [tg_c])
#a = sqa.index('A', [tg_v])
#b = sqa.index('B', [tg_v])
#x = sqa.index('X', [tg_a])
#y = sqa.index('Y', [tg_a])
#w = sqa.index('W', [tg_a])
#z = sqa.index('Z', [tg_a])

r_op = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(a)]

effH = []
effH = sqa.Heff(1)
for t in effH:
  print t

term2 = sqa.term(0.5, [], X + r_op)

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
