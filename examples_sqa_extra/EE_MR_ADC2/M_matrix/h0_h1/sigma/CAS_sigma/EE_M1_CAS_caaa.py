import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True
#i = sqa.index('i', [tg_c], dummy)
j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#l = sqa.index('l', [tg_c], dummy)
#a = sqa.index('aa', [tg_v], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#d = sqa.index('dd', [tg_v], dummy)
x = sqa.index('x', [tg_a], dummy)
y = sqa.index('y', [tg_a], dummy)
z = sqa.index('z', [tg_a], dummy)
#w = sqa.index('w', [tg_a], dummy)

# External indices
#j = sqa.index('J', [tg_c])
#x = sqa.index('X', [tg_a])
#y = sqa.index('Y', [tg_a])
#z = sqa.index('Z', [tg_a])

Xsym_1 = sqa.symmetry((0,1,3,2),-1)

X = [sqa.tensor('X', [j,x,y,z], [Xsym_1])]

r_op = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(x), sqa.desOp(j)]

effH = []
effH = sqa.Heff(1)

term1 = sqa.term(0.5, [], X + r_op)
term2 = sqa.commutator(effH, term1)
term3 = sqa.matrixBlock(term2)

sqa.generateEinsum(term3, '    temp', '', transRDM=True, trans_ind_str='I')
