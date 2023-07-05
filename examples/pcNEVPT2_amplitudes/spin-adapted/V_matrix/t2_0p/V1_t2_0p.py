import sqa_plus
sqa_plus.options.spin_integrated = True

import time
start = time.time()

sqa_plus.options.print_header("Spin-Adapted V(I,A)")

# Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices and Operators
i_a = sqa_plus.index('I', [tg_alpha, tg_cor])
a_a = sqa_plus.index('A', [tg_alpha, tg_vir])

terms_l_op = sqa_plus.term(- 0.5, [], [sqa_plus.creOp(i_a), sqa_plus.desOp(a_a)])

## Spin-Integrated Dyall H
print("## Creating spin-integrated  H^(0)...")
terms_Dyall_H = sqa_plus.dyallH()
print("{:} spin-integrated H^(0) terms created.".format(len(terms_Dyall_H)))

## Spin-Integrated T - T^\dag
print("\n## Creating spin-integrated T - T^\dag...")
terms_T = sqa_plus.Tamplitude(1)
print("{:} spin-integrated A terms created.".format(len(terms_T)))

## Commutator [H^(0), T - T^\dag]
print("\n## Calculating [H^(0), T - T^\dag]:")
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, terms_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_H0_T)))

## Spin-Integrated V
print("\n## Creating spin-integrated V...")
terms_V = sqa_plus.Vperturbation()
print("{:} spin-integrated V terms created.".format(len(terms_V)))

for term_V in terms_V:
    term_V.scale(2.0)

## Spin-Integrated H^(1)
terms_V_H1 = terms_V + terms_commutator_H0_T
print("{:} spin-integrated H^(1) terms created.".format(len(terms_V_H1)))

# Calculating [V + H^(1), T - T^\dag]
print("\n# Creating [V + H^(1), T - T^\dag]")

## Spin-Integrated T - T^\dag
print("\n## Creating spin-integrated T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(1)
print("{:} spin-integrated A terms created.".format(len(terms2_T)))

## Commutator [V + H^(1), T - T^\dag]
print("\n## Calculating [V + H^(1), T - T^\dag]:")
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, terms2_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_V_H1_T)))
 
print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
terms_V1 = []
for _term_Heff in terms_commutator_V_H1_T:
    terms_V1.append(sqa_plus.multiplyTerms(terms_l_op, _term_Heff))

## Expected value of Spin-Integrated V1
expected_terms_V1 = sqa_plus.matrixBlock(terms_V1)
print(">>> {:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(len(expected_terms_V1)))

sqa_plus.options.print_divider()

# Spin-Adaptation of V1
expected_V1_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_terms_V1)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_V1_sa, 'V1', 'IA')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
