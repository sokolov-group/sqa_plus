import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_integrated = True

import time
start = time.time()

sqa_plus.options.print_header("Spin-Integrated V(I,X): T-Slice")

# Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices and Operators
print("# Defining a^\dag_p a_q...")
i_b = sqa_plus.index('I', [tg_beta, tg_cor])
x_b = sqa_plus.index('X', [tg_beta, tg_act])

term_l_op = sqa_plus.term(-0.5, [], [sqa_plus.creOp(i_b), sqa_plus.desOp(x_b)])

## Spin-Integrated T - T^\dag
print("\n## Creating first set of spin-integrated T - T^\dag...")
terms_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-integrated A terms created.".format(len(terms_T)))

print("\n## Creating second set of spin-integrated T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-integrated A terms created.".format(len(terms2_T)))

# First Term
print("\n# Creating spin-integrated [2V + [H^(0), T - T^\dag], T - T^\dag]:")

# Spin-Integrated V
print("\n## Creating spin-integrated V...")
terms_V = sqa_plus.Vperturbation()
for term_V in terms_V:
    term_V.scale(2.0)
print("{:} spin-integrated V terms created.".format(len(terms_V)))

## Spin-Integrated Dyall H
print("\n## Creating spin-integrated H^(0)...")
terms_Dyall_H = sqa_plus.dyallH()
print("{:} spin-integrated H^(0) terms created.".format(len(terms_Dyall_H)))

## Commutator [H0, T - T^\dag]
print("\n## Calculating [H^(0), T - T^\dag]:")
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, terms_T)
print(">>> {:} spin-integrated terms created.".format(len(terms_commutator_H0_T)))

## Term V + H1
print("\n## Calculating V + H^(1):")
terms_V_H1 = terms_V + terms_commutator_H0_T
print(">>> {:} spin-integrated terms created.".format(len(terms_V_H1)))

## Commutator [V + H1, T - T^\dag]
print("\n## Calculating [V + H^(1), T - T^\dag]:")
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, terms2_T)
print(">>> {:} spin-integrated terms created.".format(len(terms_commutator_H0_T)))

## Term a^\dag_p a_q [V + H1, T - T^\dag]
print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
terms_V1 = []
for term_commutator in terms_commutator_V_H1_T:
    terms_V1.append(sqa_plus.multiplyTerms(term_l_op, term_commutator))
print(">>> {:} spin-integrated a^\dag_p a_q [V + H^(1), T - T^\dag] terms created.\n".format(len(terms_V1)))

sqa_plus.options.print_divider()

## Solving for V
print("\n## Calculating matrix element of V1:")
## Expected value of Spin-Integrated V1
expected_terms_V1 = sqa_plus.matrixBlock(terms_V1)
del terms_V1
print(">>> {:} spin-integrated a^\dag_p a_q [V + H^(1), A^(1)] terms created.".format(len(expected_terms_V1)))

# Spin-Adaptation of V1
expected_V1_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_terms_V1)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_V1_sa, 'V1', 'IX')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
