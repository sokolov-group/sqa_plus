import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_integrated = True
#sqa_plus.options.verbose = True

import time
start = time.time()

amplitudes_list = ['t1_ae', 't1_ca', 't1_ce',
                   't2_aaea', 't2_aaee', 't2_caaa', 't2_caea',
                   't2_caee', 't2_ccaa', 't2_ccea', 't2_ccee']

sqa_plus.options.print_header("Spin-Adapted V(X,Y)")

# Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices and Operators
x_a = sqa_plus.index('X', [tg_alpha, tg_act])
y_a = sqa_plus.index('Y', [tg_alpha, tg_act])

x_b = sqa_plus.index('X', [tg_beta, tg_act])
y_b = sqa_plus.index('Y', [tg_beta, tg_act])

terms_l_op = [sqa_plus.term(-0.5, [], [sqa_plus.creOp(x_a), sqa_plus.desOp(y_a)])]
##terms_l_op = [sqa_plus.term(-0.5, [], [sqa_plus.creOp(x_b), sqa_plus.desOp(y_b)])]

# Spin-Integrated H(1)
print("# Creating spin-integrated  H^(1) = V + [H^(0), T - T^\dag]:")

## Spin-Integrated Dyall H
print("## Creating spin-integrated  H^(0)...")
terms_Dyall_H = sqa_plus.dyallH()

print("{:} spin-integrated H^(0) terms created.".format(len(terms_Dyall_H)))

# Spin-Integrated V
print("\n## Creating spin-integrated V...")
terms_V = sqa_plus.Vperturbation()
print("{:} spin-integrated V terms created.".format(len(terms_V)))

## Spin-Integrated T - T^\dag
print("\n## Creating spin-integrated T - T^\dag...")
terms_T = sqa_plus.Tamplitude(order = 1)

print("{:} spin-integrated A terms created.".format(len(terms_T)))

## Spin-Integrated T - T^\dag
print("\n## Creating spin-integrated T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-integrated A terms created.".format(len(terms2_T)))

# Spin-Integrated a^{\dag}_p a_q [V + H^(1), T - T^\dag]
print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
print("# Defining a^\dag_p a_q...")

## Solving for T in slices
expected_V1 = []

print("\n## Selecting spin-integrated T - T^\dag...")
for selected_amplitude in amplitudes_list:
    selected_terms_T = []

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

    print("\n# Selecting {:}...".format(selected_amplitude))
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

    for term_T in selected_terms_T:
        print(term_T)
    print(">>> Slice {:}: {:} spin-integrated A terms selected.".format(selected_amplitude, len(selected_terms_T)))

    ## Commutator [H^(0), T - T^\dag]
    print("\n## Calculating [H^(0), T - T^\dag]:")
    terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, selected_terms_T)
    print(">>> Slice {:}: {:} spin-integrated terms created.".format(selected_amplitude, len(terms_commutator_H0_T)))

    ## Spin-Integrated H^(1)
    terms_H1 = terms_commutator_H0_T

    ## Spin-Integrated V + H^(1)
    terms_V_H1 = terms_H1
    print(">>> Slice {:}: {:} spin-integrated V + H^(1) terms created.".format(selected_amplitude, len(terms_V_H1)))

    ## Commutator [V + H^(1), T - T^\dag]
    print("\n## Calculating [V + H^(1), T - T^\dag]:")
    terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, terms2_T)
    print(">>> Slice {:}: {:} spin-integrated terms created.".format(selected_amplitude, len(terms_commutator_V_H1_T)))

    print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
    terms_V1 = []
    for term_commutator in terms_commutator_V_H1_T:
        for term_l in terms_l_op:
            terms_V1.append(sqa_plus.multiplyTerms(term_l, term_commutator))
    print(">>> Slice {:}: {:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(selected_amplitude, len(terms_V1)))

    ## Expected value of Spin-Integrated V1
    expected_terms_V1 = sqa_plus.matrixBlock(terms_V1)
    print(">>> Slice {:}: {:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(selected_amplitude, len(expected_terms_V1)))
    for term_V1 in expected_terms_V1:
        expected_V1.append(term_V1)

    print("----------------------------------------------------------------------------------\n\n")

## Solving for V
## Spin-Integrated V + H^(1)
print(">>> Slice V: {:} spin-integrated V terms created.".format(len(terms_V)))
for term_V in terms_V:
    term_V.scale(2.0)

## Commutator [V + H^(1), T - T^\dag]
print("\n## Calculating [V + H^(1), T - T^\dag]:")
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V, terms2_T)
print(">>> Slice V: {:} spin-integrated terms created.".format(len(terms_commutator_V_H1_T)))

print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
terms_V1 = []
for term_commutator in terms_commutator_V_H1_T:
    for term_l in terms_l_op:
        terms_V1.append(sqa_plus.multiplyTerms(term_l, term_commutator))
print(">>> Slice V: {:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(len(terms_V1)))

## Expected value of Spin-Integrated V1
expected_terms_V1 = sqa_plus.matrixBlock(terms_V1)
print(">>> Slice V: {:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(len(expected_terms_V1)))
expected_V1.extend(expected_terms_V1)

for term_V1 in expected_terms_V1:
    expected_V1.append(term_V1)

print("----------------------------------------------------------------------------------\n\n")

# Final combineTerms
sqa_plus.combineTerms(expected_V1)
print("### Total: {:} spin-integrated a^\dag_p a_q [V + H1, A1] terms created.\n".format(len(expected_V1)))

sqa_plus.options.print_divider()

# Spin-Adaptation of V1
#from custom_remove_rdm import remove_high_order_rdms
#sqa_plus.options.add_spin_adaptation_custom_function(remove_high_order_rdms)
##made function in sqaSpinAdapted to handle this

#expected_V1_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_V1, remove_5rdms = True)
code_to_profile = '''sqa_plus.convertSpinIntegratedToAdapted(expected_V1, remove_5rdms = True)'''
cProfile.run(code_to_profile)

## Generating Numpy einsum equations
#result = sqa_plus.genEinsum(expected_V1_sa, 'V1', 'XY', remove_core_integrals = True)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
