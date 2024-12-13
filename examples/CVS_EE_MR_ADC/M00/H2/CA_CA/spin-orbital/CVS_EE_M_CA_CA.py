import sqa_plus
sqa_plus.options.spin_orbital = True
sqa_plus.options.cvs_approach = True
order_Heff = 2

import time
start = time.time()

indices_string = 'ca_ca'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: M00 H{:} {:}".format(order_Heff, indices_string.upper()))

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

# Generating Term 
## External Indices
i = sqa_plus.index('I', [tg_cvs_cor])
a = sqa_plus.index('A', [tg_vir])
j = sqa_plus.index('J', [tg_cvs_cor])
b = sqa_plus.index('B', [tg_vir])
x = sqa_plus.index('X', [tg_act])
y = sqa_plus.index('Y', [tg_act])

## Define operators
cre_i = sqa_plus.creOp(i)
des_a = sqa_plus.desOp(a)
des_j = sqa_plus.desOp(j)
cre_b = sqa_plus.creOp(b)
des_x = sqa_plus.desOp(x)
cre_y = sqa_plus.creOp(y)

print("\n## Generating Term a_I^\dag a_X ... a_Y\dag a_J ...")
final_indices_string = 'IXJY'

l_op = sqa_plus.term(1.0, [], [cre_i, des_x])
r_op = sqa_plus.term(1.0, [], [cre_y, des_j])

print("\n## Left operator terms:")
print(l_op)

print("\n## Right operator terms:")
print(r_op)

# Spin-Orbital H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculate the commutators
print("\n## Calculate the commutator ... [H(2), r_op] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, r_op)

print("\n## Calculate the commutator [l_op, [H(2), r_op]] ...")
terms_EE_M0 = sqa_plus.commutator(l_op, terms_commutator)

# Expected value of Spin-Orbital EE M00 H2 CA-CA
expected_EE_M0 = sqa_plus.matrixBlock(terms_EE_M0)

# Generate Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M0, 'temp', final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))

