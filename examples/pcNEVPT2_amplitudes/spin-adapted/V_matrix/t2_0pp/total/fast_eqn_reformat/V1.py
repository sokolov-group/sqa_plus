import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_integrated = True

import time
start = time.time()

selected_amplitude = 'ca_caaa'

sqa_plus.options.print_header("Spin-Adapted V(X,Y): T-Slice {:}".format(selected_amplitude))

# Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices and Operators
print("# Defining a^\dag_p a_q...")
x_a = sqa_plus.index('X', [tg_alpha, tg_act])
y_a = sqa_plus.index('Y', [tg_alpha, tg_act])

term_l_op = sqa_plus.term(-0.5, [], [sqa_plus.creOp(x_a), sqa_plus.desOp(y_a)])

## Spin-Integrated T - T^\dag
print("\n## Creating first set of spin-integrated T - T^\dag...")
terms_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-integrated A terms created.".format(len(terms_T)))

print("\n## Creating second set of spin-integrated T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-integrated A terms created.".format(len(terms2_T)))

## Solving for T in slices
print("\n## Selecting spin-integrated T - T^\dag slice: {:}...".format(selected_amplitude))
selected_terms_T = []
selected_terms2_T = []

if selected_amplitude == 'ce_caea':
    ind_t1 = [sqa_plus.options.core_type, sqa_plus.options.virtual_type]
    ind_t2 = [sqa_plus.options.core_type, sqa_plus.options.active_type,
              sqa_plus.options.virtual_type, sqa_plus.options.active_type]
    ind_list_select = [ind_t1, ind_t2]

elif selected_amplitude == 'ae_aaea':
    ind_t1 = [sqa_plus.options.active_type, sqa_plus.options.virtual_type]
    ind_t2 = [sqa_plus.options.active_type, sqa_plus.options.active_type,
              sqa_plus.options.virtual_type, sqa_plus.options.active_type]
    ind_list_select = [ind_t1, ind_t2]

elif selected_amplitude == 'ca_caaa':
    ind_t1 = [sqa_plus.options.core_type, sqa_plus.options.active_type]
    ind_t2 = [sqa_plus.options.core_type, sqa_plus.options.active_type,
                       sqa_plus.options.active_type, sqa_plus.options.active_type]
    ind_list_select = [ind_t1, ind_t2]

elif selected_amplitude == 't2_aaee':
     ind_list_select = [[sqa_plus.options.active_type, sqa_plus.options.active_type,
                     sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

elif selected_amplitude == 't2_caee':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.active_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

elif selected_amplitude == 't2_ccaa':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type,
                       sqa_plus.options.active_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't2_ccea':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't2_ccee':
    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type,
                       sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

print("\n# First Amplitude Set...")
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

print("\n# Second Amplitude Set...")
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
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, selected_terms_T)
print(">>> Slice {:}: {:} spin-integrated terms created.".format(selected_amplitude, len(terms_commutator_H0_T)))

## Term V + H1
print("\n## Calculating V + H^(1):")
terms_V_H1 = terms_V + terms_commutator_H0_T
print(">>> Slice {:}: {:} spin-integrated terms created.".format(selected_amplitude, len(terms_V_H1)))

## Commutator [V + H1, T - T^\dag]
print("\n## Calculating [V + H^(1), T - T^\dag]:")
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, selected_terms2_T)
print(">>> Slice {:}: {:} spin-integrated terms created.".format(selected_amplitude, len(terms_commutator_H0_T)))

## Term a^\dag_p a_q [V + H1, T - T^\dag]
print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
terms_V1 = []
for term_commutator in terms_commutator_V_H1_T:
    terms_V1.append(sqa_plus.multiplyTerms(term_l_op, term_commutator))
print(">>> Slice {:}: {:} spin-integrated a^\dag_p a_q [V + H^(1), T - T^\dag] terms created.\n".format(selected_amplitude, len(terms_V1)))

sqa_plus.options.print_divider()

## Solving for V
print("\n## Calculating matrix element of {:} slice of V1:".format(selected_amplitude))
## Expected value of Spin-Integrated V1
expected_terms_V1 = sqa_plus.matrixBlock(terms_V1)
del terms_V1
print(">>> Slice V1: {:} spin-integrated a^\dag_p a_q [V + H^(1), A^(1)] terms created.".format(len(expected_terms_V1)))

# Spin-Adaptation of V1
expected_terms_V1_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_terms_V1, remove_5rdms = True)
print(">>> Slice V1: {:} spin-adapted a^\dag_p a_q [V + H^(1), A^(1)] terms created.\n".format(len(expected_terms_V1_sa)))

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_terms_V1_sa, 'V1', 'XY', remove_core_integrals = True)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
