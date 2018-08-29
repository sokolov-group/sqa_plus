import sqa.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Derivation of the following term:
# <\Psi_0| x^\dag a [\hat{V}, b^\dag y] |\Psi_0>
# = \sum_{cdwz} v_{cz}^{dw}[<\Psi_0| x^\dag a {c^\dag d} z^\dag w b^\dag y |\Psi_0> - <\Psi_0| x^\dag a b^\dag y {c^\dag d} z^\dag w |\Psi_0>]
# - \sum_{cdwz} v_{cz}^{dw} <\Psi_0| z^\dag w \Psi_0> [<\Psi_0| x^\dag a {c^\dag d} b^\dag y |\Psi_0> - <\Psi_0| x^\dag a b^\dag y {c^\dag d} |\Psi_0>]

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
a = sqa.index('aa', [tg_v])
b = sqa.index('bb', [tg_v])

# Operators with external indices
xa_op = [sqa.creOp(x), sqa.desOp(a)]
by_op = [sqa.creOp(b), sqa.desOp(y)]

# V tensor
v_symmetry = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2), -1)]
v_t = [sqa.tensor('v', [vv[0], aa[0], vv[1], aa[1]], v_symmetry)]

# V operator
v_op = [sqa.creOp(vv[0]), sqa.desOp(vv[1]), sqa.creOp(aa[0]), sqa.desOp(aa[1])]

# V <gamma_pq> tensor
v_gamma_symmetry = [sqa.symmetry((1,0),1)]
v_gamma_t = [sqa.tensor('v_gamma', [vv[2], vv[3]], v_gamma_symmetry)]

# V <gamma_pq> operator
v_gamma_op = [sqa.creOp(vv[2]), sqa.desOp(vv[3])]

# Define terms
terms = [sqa.term(1.0, [], v_t + xa_op + v_op + by_op)]
terms += [sqa.term(-1.0, [], v_t + xa_op + by_op + v_op)]
terms += [sqa.term(-1.0, [], v_gamma_t + xa_op + v_gamma_op + by_op)]
terms += [sqa.term(1.0, [], v_gamma_t + xa_op + by_op + v_gamma_op)]

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

