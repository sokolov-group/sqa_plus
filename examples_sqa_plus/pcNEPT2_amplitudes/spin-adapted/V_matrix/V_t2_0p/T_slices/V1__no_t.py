import sqa_plus
spin_integrated = True

print ("# Spin-Adapted V_{(I,A)}:\n")

# Define indices
dummy = True

tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

dummy_indices_list = sqa_plus.create_dummy_indices_list(spin_integrated)

# Spin-Integrated H(1)
print("# Creating spin-integrated  H^(1) = V + [H^(0), T - T^\dag]:")

## Spin-Integrated Dyall H
print("## Creating spin-integrated  H^(0)...")
terms_Dyall_H = sqa_plus.dyallH(dummy_indices_list, spin_integrated)
print("{:} spin-integrated H^(0) terms created.".format(len(terms_Dyall_H)))

## Spin-Integrated T - T^\dag
print("\n## Creating spin-integrated T - T^\dag...")
terms_T = sqa_plus.Tamplitude(dummy_indices_list, 1, spin_integrated)
print("{:} spin-integrated A terms created.".format(len(terms_T)))

# Spin-Integrated V
print("\n## Creating spin-integrated V...")
terms_V = sqa_plus.Vperturbation(dummy_indices_list, spin_integrated)
print("{:} spin-integrated V terms created.".format(len(terms_V)))

## Spin-Integrated H^(1)
terms_H1 = terms_V

## Commutator [H^(0), T - T^\dag]
print("\n## Calculating [H^(0), T - T^\dag]:")
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, terms_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_H0_T)))

## Spin-Integrated V + H^(1)
terms_V_H1 = terms_H1
print("{:} spin-integrated V + H^(1) terms created.".format(len(terms_V_H1)))
for term_V in terms_V_H1:
    term_V.scale(2.0)

## Spin-Integrated T - T^\dag
print("\n## Creating spin-integrated T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(dummy_indices_list, 1, spin_integrated)
print("{:} spin-integrated A terms created.".format(len(terms2_T)))

## Commutator [V + H^(1), T - T^\dag]
print("\n## Calculating [V + H^(1), T - T^\dag]:")
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, terms2_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_V_H1_T)))

# Spin-Integrated a^{\dag}_p a_q [V + H^(1), T - T^\dag]
print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
print("# Defining a^\dag_p a_q...")

## External Indices and Operators
i_a = sqa_plus.index('I', [tg_alpha, tg_cor])
a_a = sqa_plus.index('A', [tg_alpha, tg_vir])

terms_l_op = [sqa_plus.term(0.5, [], [sqa_plus.creOp(i_a), sqa_plus.desOp(a_a)])]

print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
terms_V1 = []
for term_commutator in terms_commutator_V_H1_T:
    for term_l in terms_l_op:
        terms_V1.append(sqa_plus.multiplyTerms(term_l, term_commutator))
print("{:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(len(terms_V1)))

## Expected value of Spin-Integrated V1
expected_V1 = sqa_plus.matrixBlock(terms_V1)

# Spin-Adaptation of V1
expected_V1_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_V1)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_V1_sa, 'V1', 'IA', rm_core_int = True, suffix = '')

print ("\n# Equations using genEinsum:")
for item in result:
    print(item)

