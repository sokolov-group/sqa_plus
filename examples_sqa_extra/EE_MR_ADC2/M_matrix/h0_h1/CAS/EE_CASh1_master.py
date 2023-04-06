import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = False

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v


# RHS
#j = sqa.index('J', [tg_c])
#k = sqa.index('K', [tg_c])
#x = sqa.index('X', [tg_a])
#y = sqa.index('Y', [tg_a])
#r_op  = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(k), sqa.desOp(j)]   # CCAA
#r_ind = 'JKXY'


#j = sqa.index('J', [tg_c])
#k = sqa.index('K', [tg_c])
#a = sqa.index('A', [tg_v])
#x = sqa.index('X', [tg_a])
#r_op  = [sqa.creOp(a), sqa.creOp(x), sqa.desOp(k), sqa.desOp(j)]   # CCEA
#r_ind = 'JKAX'


#j = sqa.index('J', [tg_c])
#k = sqa.index('K', [tg_c])
#a = sqa.index('A', [tg_v])
#b = sqa.index('B', [tg_v])
#r_op  = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(k), sqa.desOp(j)]   # CCEE
#r_ind = 'JKAB'


#j = sqa.index('J', [tg_c])
#x = sqa.index('X', [tg_a])
#a = sqa.index('A', [tg_v])
#b = sqa.index('B', [tg_v])
#r_op  = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(x), sqa.desOp(j)]   # CAEE
#r_ind = 'JXAB'


#x = sqa.index('X', [tg_a])
#y = sqa.index('Y', [tg_a])
#a = sqa.index('A', [tg_v])
#b = sqa.index('B', [tg_v])
#r_op  = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(y), sqa.desOp(x)]   # AAEE
#r_ind = 'XYAB'


#j = sqa.index('J', [tg_c])
#x = sqa.index('X', [tg_a])
#y = sqa.index('Y', [tg_a])
#z = sqa.index('Z', [tg_a])
#r_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(x), sqa.desOp(j)]   # CAAA
#r_ind = 'JXYZ'


#j = sqa.index('J', [tg_c])
#x = sqa.index('X', [tg_a])
#a = sqa.index('A', [tg_v])
#y = sqa.index('Y', [tg_a])
#r_op  = [sqa.creOp(a), sqa.creOp(y), sqa.desOp(x), sqa.desOp(j)]   # CAEA
#r_ind = 'JXAY'


x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
a = sqa.index('A', [tg_v])
z = sqa.index('Z', [tg_a])
r_op  = [sqa.creOp(a), sqa.creOp(z), sqa.desOp(y), sqa.desOp(x)]   # AAEA
r_ind = 'XYAZ'


#############################
#############################


effH = []
effH = sqa.Heff(1)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.commutator(effH, term1)
term3 = sqa.matrixBlock(term2)

sqa.generateEinsum(term3, 'temp', r_ind, transRDM=True, trans_ind_str='I')
