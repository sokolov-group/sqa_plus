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
i = sqa.index('I', [tg_c])
j = sqa.index('J', [tg_c])
k = sqa.index('K', [tg_c])
l = sqa.index('L', [tg_c])
a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
c = sqa.index('C', [tg_v])
d = sqa.index('D', [tg_v])
x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
w = sqa.index('W', [tg_a])
z = sqa.index('Z', [tg_a])

#l_op = [sqa.desOp(a)] # A
#l_op = [sqa.creOp(x), sqa.desOp(y), sqa.desOp(a)] # XYA
#l_op = [sqa.creOp(i), sqa.desOp(a), sqa.desOp(b)] # IAB
#l_op = [sqa.creOp(i), sqa.desOp(x), sqa.desOp(a)] # IXA
#l_op = [sqa.creOp(i), sqa.desOp(x), sqa.desOp(y)] # IXY
#l_op = [sqa.creOp(x), sqa.desOp(a), sqa.desOp(b)] # XAB
#r_op = [sqa.creOp(c)] # C
#r_op = [sqa.creOp(c), sqa.creOp(z), sqa.desOp(w)] # WZC
r_op = [sqa.creOp(d), sqa.creOp(c), sqa.desOp(j)] # JCD
#r_op = [sqa.creOp(c), sqa.creOp(w), sqa.desOp(j)] # JWC
#r_op = [sqa.creOp(z), sqa.creOp(w), sqa.desOp(j)] # JWZ
#r_op = [sqa.creOp(d), sqa.creOp(c), sqa.desOp(w)] # WCD

effH = []
effH = sqa.Heff(1)

#term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], r_op)

print "First Commutator"
term3 = sqa.commutator(effH,term2)

#print "Second Commutator"
#term4 = []
#for t in term3:
#    term4.append(sqa.multiplyTerms(term1,t))
#    term4.append(sqa.multiplyTerms(t,term1))

term5 = sqa.matrixBlock(term3)

sqa.generateEinsum(term5, 'temp', 'JCD', transRDM=True, trans_ind_str="I")

