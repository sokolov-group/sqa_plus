import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_integrated = True

import time
start = time.time()

selected_amplitude = 't2_ccaa'

sqa_plus.options.print_header("Spin-Adapted V(X,Y): T-Slice {:}".format(selected_amplitude))

# Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices and Operators
x_b = sqa_plus.index('X', [tg_beta, tg_act])
y_b = sqa_plus.index('Y', [tg_beta, tg_act])

terms_l_op = sqa_plus.term(-1.0, [], [sqa_plus.creOp(x_b), sqa_plus.desOp(y_b)])

## Spin-Integrated Dyall H
print("## Creating spin-integrated  H^(0)...")
terms_Dyall_H = sqa_plus.dyallH()
print("{:} spin-integrated H^(0) terms created.".format(len(terms_Dyall_H)))

## Spin-Integrated T - T^\dag
print("\n## Creating spin-integrated T - T^\dag...")
terms_T = sqa_plus.Tamplitude(1)
terms2_T = sqa_plus.Tamplitude(1)
terms3_T = sqa_plus.Tamplitude(1)
print("{:} spin-integrated A terms created.".format(len(terms_T)))

## Solving for T in slices
print("\n## Selecting spin-integrated slice of T - T^\dag: {:}...".format(selected_amplitude))
selected_terms_T = []
selected_terms2_T = []
selected_terms3_T = []

if selected_amplitude == 'ce_caea':
    ind_t1 = [sqa_plus.options.core_type, sqa_plus.options.virtual_type]
    ind_t2 = [sqa_plus.options.core_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type]
    ind_list_select = [ind_t1, ind_t2]

elif selected_amplitude == 'ae_aaea':
    ind_t1 = [sqa_plus.options.active_type, sqa_plus.options.virtual_type]
    ind_t2 = [sqa_plus.options.active_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type]
    ind_list_select = [ind_t1, ind_t2]

elif selected_amplitude == 'ca_caaa':
    ind_t1 = [sqa_plus.options.core_type, sqa_plus.options.active_type]
    ind_t2 = [sqa_plus.options.core_type, sqa_plus.options.active_type, sqa_plus.options.active_type, sqa_plus.options.active_type]
    ind_list_select = [ind_t1, ind_t2]

elif selected_amplitude == 't2_aaee':
     ind_list_select = [[sqa_plus.options.active_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

elif selected_amplitude == 't2_caee':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

elif selected_amplitude == 't2_ccaa':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type, sqa_plus.options.active_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't2_ccea':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't2_ccee':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

print("\n# Amplitude Set...")
for term_T in terms_T:
    select_term = False
    for tensor_T in term_T.tensors:
        ind_types = []
        for ind_T in tensor_T.indices:
            ind_types.append(get_spatial_index_type(ind_T.indType))
        if ind_types in ind_list_select:
            select_term = True
    if select_term:
        selected_terms_T.append(term_T)

for term_T in selected_terms_T:
    print(term_T)
print(">>> Slice {:}: {:} spin-integrated A terms selected.".format(selected_amplitude, len(selected_terms_T)))

for term_T in terms2_T:
    select_term = False
    for tensor_T in term_T.tensors:
        ind_types = []
        for ind_T in tensor_T.indices:
            ind_types.append(get_spatial_index_type(ind_T.indType))
        if ind_types in ind_list_select:
            select_term = True
    if select_term:
        selected_terms2_T.append(term_T)

for term_T in selected_terms2_T:
    print(term_T)
print(">>> Slice {:}: {:} spin-integrated A terms selected.".format(selected_amplitude, len(selected_terms2_T)))

for term_T in terms3_T:
    select_term = False
    for tensor_T in term_T.tensors:
        ind_types = []
        for ind_T in tensor_T.indices:
            ind_types.append(get_spatial_index_type(ind_T.indType))
        if ind_types in ind_list_select:
            select_term = True
    if select_term:
        selected_terms3_T.append(term_T)

for term_T in selected_terms3_T:
    print(term_T)
print(">>> Slice {:}: {:} spin-integrated A terms selected.".format(selected_amplitude, len(selected_terms3_T)))

## Spin-Integrated V
print("\n## Creating spin-integrated V...")
terms_V = sqa_plus.Vperturbation()
print("{:} spin-integrated V terms created.".format(len(terms_V)))

## Commutator [V, T - T^\dag]
print("\n## Calculating [V, T - T^\dag]:")
terms_commutator_V_T = sqa_plus.commutator(terms_V, selected_terms_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_V_T)))

print("\n## Calculating a^\dag_p a_q [V, T - T^\dag]:")
terms_V1 = []
for _term_Heff in terms_commutator_V_T:
    terms_V1.append(sqa_plus.multiplyTerms(terms_l_op, _term_Heff))

## Expected value of Spin-Integrated V1
expected_terms_V1 = sqa_plus.matrixBlock(terms_V1)
print(">>> {:} spin-integrated a^\dag_p a_q [V, A1] terms created.\n".format(len(expected_terms_V1)))

###################################################

## Commutator [H^(0), T - T^\dag]
print("\n## Calculating [H^(0), T - T^\dag]:")
terms_commutator_H_T = sqa_plus.commutator(terms_Dyall_H, selected_terms2_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_H_T)))

## Commutator [[H^(0), T - T^\dag], T - T^\dag]
print("\n## Calculating [[H^(0), T - T^\dag], T - T^\dag]:")
terms_commutator_H_T_T = sqa_plus.commutator(terms_commutator_H_T, selected_terms3_T)
print("{:} spin-integrated terms created.".format(len(terms_commutator_H_T_T)))

print("\n## Calculating a^\dag_p a_q [[H^(0), T - T^\dag], T - T^\dag]:")
terms_V1 = []
for _term_Heff in terms_commutator_H_T_T:
    terms_V1.append(sqa_plus.multiplyTerms(terms_l_op, _term_Heff))

## Expected value of Spin-Integrated V1
hold = sqa_plus.matrixBlock(terms_V1)
print(">>> {:} spin-integrated a^\dag_p a_q [[H0, A1], A1] terms created.\n".format(len(hold)))

for t in hold:
    t.scale(0.5)
expected_terms_V1.extend(hold)
sqa_plus.combineTerms(expected_terms_V1)
print(">>> {:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(len(expected_terms_V1)))

sqa_plus.options.print_divider()

# Spin-Adaptation of V1
expected_V1_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_terms_V1, remove_5rdms = True)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_V1_sa, 'V1', 'XY')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
