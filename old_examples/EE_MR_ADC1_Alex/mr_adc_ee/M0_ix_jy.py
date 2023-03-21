import sqa_extra.secondQuantizationAlgebra as sqa
#import sqa_extra.sqaAddon as addonkc
#from sqaAddon import addon,normalOrderCore

sqa.options.verbose = True
#
# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# Core dummy indices
cc = [sqa.index('c%i' %p, [tg_c], dummy) for p in range(10)]
# Active dummy indices
aa = [sqa.index('a%i' %p, [tg_a], dummy) for p in range(10)]
# Virtual dummy indices
vv = [sqa.index('v%i' %p, [tg_v], dummy) for p in range(10)]

# External indices
x = sqa.index('xx', [tg_a])
y = sqa.index('yy', [tg_a])
i = sqa.index('ii', [tg_c])
j = sqa.index('jj', [tg_c])

# Operators with external indices
# Note that we define i^\dag as q-annihilation operator with index i, etc
ix_op = [sqa.creOp(i), sqa.desOp(x)]
yj_op = [sqa.creOp(y), sqa.desOp(j)]
#
# Hamiltonian operator : 
# H^0 = E_fc - SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a} + H_act
Hamil = []
#
Hamil.append(sqa.term(1.0, ['E_fc'], []))
cor = cc.pop(0)
vir = vv.pop(0)
Hamil.append(sqa.term(1.0, ['e_i'],[sqa.creOp(cor),sqa.desOp(cor)]))
#Hamil.append(sqa.term(-1.0, ['e_i'],[sqa.desOp(cor),sqa.creOp(cor)]))
Hamil.append(sqa.term(1.0, ['e_v'],[sqa.creOp(vir),sqa.desOp(vir)]))
#act = aa.pop(0)
#cd = [sqa.creOp(act),sqa.desOp(act)]
#credes = sqa.creDesTensor(cd)
#Hamil.append(sqa.term(1.0, [], [credes]))
Hact = sqa.tensor('Hact', [], [])
Hamil.append(sqa.term(1.0, [], [Hact]))
#
for t in Hamil:
 print t
#
#term1 = sqa.commutator(sqa.term(1.0, [], ix_op),sqa.commutator(Hamil,sqa.term(1.0, [], yj_op)))
term1 = sqa.commutator(Hamil,sqa.term(1.0, [], yj_op))
#
#sqa.combine_transpose(term1)
#
# Print first commutator
print "First Commutator"
for t in term1:
 print t
#
term2 = sqa.commutator(sqa.term(1.0, [], ix_op),term1)
#
#sqa.combine_transpose(term2)
#
# Print second commutator
print "Second Commutator"
for t in term2:
 print t
#
sqa.addon(term2)


effH = []
sqa.Heff(1)
exit()
















# V tensor. Note that we permute first two indices to compensate the minus sign introduced in the V operator (see below)
v_symmetry = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2), -1)]
v_t = [sqa.tensor('v', [aa[0], cc[0], cc[1], aa[1]], v_symmetry)]

# V operator
v_op = [sqa.creOp(cc[1]), sqa.desOp(cc[0]), sqa.creOp(aa[0]), sqa.desOp(aa[1])]

# V <gamma_pq> tensor
v_gamma_symmetry = [sqa.symmetry((1,0),1)]
v_gamma_t = [sqa.tensor('v_gamma', [cc[3], cc[2]], v_gamma_symmetry)]

# V <gamma_pq> operator
v_gamma_op = [sqa.creOp(cc[3]), sqa.desOp(cc[2])]

# Define terms. Note the sign change in the third and fourth terms due to permutation {k^\dag l} = - {l k^\dag}
terms = [sqa.term(1.0, [], v_t + ix_op + v_op + yj_op)]
terms += [sqa.term(-1.0, [], v_t + ix_op + yj_op + v_op)]
terms += [sqa.term(1.0, [], v_gamma_t + ix_op + v_gamma_op + yj_op)]
terms += [sqa.term(-1.0, [], v_gamma_t + ix_op + yj_op + v_gamma_op)]
#for t in terms:
#    print t

# Normal order terms
noTerms = []

for t in terms:
    t_no = sqa.normalOrder(t)
    noTerms.extend(t_no)

for t in noTerms:
    t.contractDeltaFuncs()

# Remove zero terms and combine like terms
sqa.termChop(noTerms)
sqa.combineTerms(noTerms)

# Print the final results
print "Final results:"
for t in noTerms:
    index_types = ()
    for t_tensor in t.tensors:
        for t_tensor_index in t_tensor.indices:
            index_types += t_tensor_index.indType[0]
    print t
    print index_types
#
#
#i############## Koushik Chatterjee ##############
print ""
print ""
print " koushik check"
#for t in noTerms:
print noTerms[-1]

print noTerms[-1].tensors[0]
print noTerms[-1].tensors[0].indices[2].name
print noTerms[-1].tensors[0].indices[2].indType
print vv[0].indType
print vv[0].indType[0]
print noTerms[-1].tensors[0].indices[2].isSummed

print ""
print ""
print " koushik calling"
sqa.addon(noTerms)
exit()
coreterm = []
for t in noTerms:
  print t
  tno = sqa.normalOrderCore(t)

exit()

#  coreterm.extend(tno)
#
for t in coreterm:
    t.contractDeltaFuncs()
#
#
sqa.termChop(coreterm)
sqa.combineTerms(coreterm)
#
for t in coreterm:
    print t
#addon(noTerms)
################################################
#
#













