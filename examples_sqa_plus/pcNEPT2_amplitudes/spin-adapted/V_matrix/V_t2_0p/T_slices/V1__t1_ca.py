import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
spin_integrated = True

selected_amplitude = 't1_ca'

print ("# Spin-Adapted V_{(I,A)}:\n")

# Define indices
dummy = True

tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

dummy_indices_list = sqa_plus.create_dummy_indices_list(spin_integrated)

if selected_amplitude == 't1_ae':
    ind_list_select = [sqa_plus.options.active_type, sqa_plus.options.virtual_type]

elif selected_amplitude == 't1_ca':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.active_type]

elif selected_amplitude == 't1_ce':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.virtual_type]

elif selected_amplitude == 't2_aaea':
    ind_list_select = [sqa_plus.options.active_type, sqa_plus.options.active_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.active_type]

elif selected_amplitude == 't2_aaee':
    ind_list_select = [sqa_plus.options.active_type, sqa_plus.options.active_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]

elif selected_amplitude == 't2_caaa':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.active_type,
                       sqa_plus.options.active_type, sqa_plus.options.active_type]

elif selected_amplitude == 't2_caea':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.active_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.active_type]

elif selected_amplitude == 't2_caee':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.active_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]

elif selected_amplitude == 't2_ccaa':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.core_type,
                       sqa_plus.options.active_type, sqa_plus.options.active_type]

elif selected_amplitude == 't2_ccea':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.core_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.active_type]

elif selected_amplitude == 't2_ccee':
    ind_list_select = [sqa_plus.options.core_type, sqa_plus.options.core_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]

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

print("\n## Selecting spin-integrated T - T^\dag...")
selected_terms_T = []

for term_T in terms_T:
    select_term = False
    for tensor_T in term_T.tensors:
        ind_types = []
        for ind_T in tensor_T.indices:
            ind_types.append(get_spatial_index_type(ind_T.indType))
        if ind_types == ind_list_select:
            select_term = True
    if select_term:
        selected_terms_T.append(term_T)

terms_T = selected_terms_T
for term_T in selected_terms_T:
    print(term_T)
print("{:} spin-integrated A terms selected.".format(len(terms_T)))

## Commutator [H^(0), T - T^\dag]
print("\n## Calculating [H^(0), T - T^\dag]:")
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, terms_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_H0_T)))

## Spin-Integrated H^(1)
terms_H1 = terms_commutator_H0_T

## Spin-Integrated V + H^(1)
terms_V_H1 = terms_H1
print("{:} spin-integrated V + H^(1) terms created.".format(len(terms_V_H1)))

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

