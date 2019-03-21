import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# Core dummy indices
cc = [sqa.index('c%i' %p, [tg_c], dummy) for p in range(50)]
# Active dummy indices
aa = [sqa.index('a%i' %p, [tg_a], dummy) for p in range(50)]
# Virtual dummy indices
vv = [sqa.index('v%i' %p, [tg_v], dummy) for p in range(50)]

# External indices
a = sqa.index('A', [tg_v], False)
b = sqa.index('B', [tg_v], False)
x = sqa.index('X', [tg_a], False)
y = sqa.index('Y', [tg_a], False)

c = sqa.index('C', [tg_v], False)
d = sqa.index('D', [tg_v], False)
z = sqa.index('Z', [tg_a], False)
w = sqa.index('W', [tg_a], False)

i = sqa.index('I', [tg_c], False)
j = sqa.index('J', [tg_c], False)
k = sqa.index('K', [tg_c], False)
l = sqa.index('L', [tg_c], False)
p = sqa.index('P', [tg_v], False)
#p = sqa.index('P', [tg_c + tg_a])


# Operators with external indices
# Note that we define i^\dag as q-annihilation operator with index i, etc
p_op = [sqa.desOp(p)]
i_op = [sqa.creOp(i)]
j_op = [sqa.creOp(j)]

#
term1 = sqa.term(1.0, [], p_op)
term2 = sqa.term(1.0, [], j_op)
#term2 = sqa.term(1.0, [], [sqa.creOp(i), sqa.creOp(j), sqa.desOp(x)])

T1 = []
cc1 = []
aa1 = []
vv1 = []
for i in range(4):
    cc1.append(cc.pop(0))
    aa1.append(aa.pop(0))
    vv1.append(vv.pop(0))
T1.extend(sqa.Tamplitude(T1, 1, cc1, aa1, vv1))

T2 = []
for i in range(4):
    cc1.append(cc.pop(0))
    aa1.append(aa.pop(0))
    vv1.append(vv.pop(0))
T2.extend(sqa.Tamplitude(T2, 2, cc1, aa1, vv1))

q0_T2 = []
q0_T2.extend(sqa.commutator(term1,T2))

q0_T1 = []

q0_T1.extend(sqa.commutator(term1,T1))

q0_T1_T1 = []
q0_T1_T1.extend(sqa.commutator(q0_T1,T1))

for t in q0_T1_T1:
  t.scale(0.5)

q0_T2.extend(q0_T1_T1)


## q^1
#q_terms = q0_T1

## q^2
q_terms = q0_T2

for t in q_terms:
 print t



fterms = []

for t in q_terms:
   fterms.append(sqa.multiplyTerms(t,term2))
   fterms.append(sqa.multiplyTerms(term2,t))

#fterms.append(sqa.multiplyTerms(term1,term2))
#fterms.append(sqa.multiplyTerms(term2,term1))
#
term4 = []
for t in q_terms:
#for t in fterms:
  tt = sqa.normalOrder(t)
  term4.extend(tt)
#
term5 = sqa.matrixBlock(term4)
#term5 = sqa.matrixBlock(fterms)
#term5 = sqa.matrixBlock(q_terms)
#term5 = sqa.matrixBlock(T1)

sqa.generateEinsum(term5, lhs_str = 'T1[nocc_so:, si:fi]', ind_str = 'P', transRDM = True, trans_ind_str = 'I')
#sqa.generateEinsum(term5, lhs_str = 'T1[:ncore_so, si:fi]', ind_str = 'PJ')

#
#
##################################################################
