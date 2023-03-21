import sqa_extra.secondQuantizationAlgebra as sqa
#import sqa_extra.sqaAddon as addonkc
#from sqaAddon import addon,normalOrderCore

sqa.options.verbose = True

# Derivation of the following term:
# <\Psi_0| i^\dag x [\hat{V}, y^\dag j] |\Psi_0>
# = \sum_{klwz} v_{kz}^{lw} [<\Psi_0| i^\dag x {k^\dag l} z^\dag w y^\dag j|\Psi_0> - <\Psi_0| i^\dag x y^\dag j {k^\dag l} z^\dag w] |\Psi_0>]
# - \sum_{klwz} v_{kz}^{lw} <\Psi_0| z^\dag w \Psi_0> [<\Psi_0| i^\dag x {k^\dag l} y^\dag j |\Psi_0> - <\Psi_0| i^\dag x y^\dag j {k^\dag l} |\Psi_0>]

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# Core dummy indices
cc = [sqa.index('c%i' %p, [tg_c], dummy) for p in range(6)]
# Active dummy indices
aa = [sqa.index('a%i' %p, [tg_a], dummy) for p in range(6)]
# Virtual dummy indices
vv = [sqa.index('v%i' %p, [tg_v], dummy) for p in range(6)]

# External indices
x = sqa.index('xx', [tg_a])
y = sqa.index('yy', [tg_a])
i = sqa.index('ii', [tg_c])
j = sqa.index('jj', [tg_c])

# Operators with external indices
# Note that we define i^\dag as q-annihilation operator with index i, etc
ix_op = [sqa.desOp(i), sqa.desOp(x)]
yj_op = [sqa.creOp(y), sqa.creOp(j)]

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













