import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_orbital = True
sqa_plus.options.cvs_approach = True

import time
start = time.time()

selected_amplitude = 't1_m1p'

sqa_plus.options.print_header("Spin-Orbital CVS-EE: M00 H2 CA_CA | T-Slice {:}".format(selected_amplitude))

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

# Generating Term 
## External Indices and Operators
i = sqa_plus.index('I', [tg_cvs_cor])
j = sqa_plus.index('J', [tg_cvs_cor])
x = sqa_plus.index('X', [tg_act])
y = sqa_plus.index('Y', [tg_act])

## Define operators
cre_i = sqa_plus.creOp(i)
des_j = sqa_plus.desOp(j)
des_x = sqa_plus.desOp(x)
cre_y = sqa_plus.creOp(y)

## Define terms
print("\n## Generating Term a_I^\dag a_X ... a_Y^\dag a_J ...")
final_indices_string = 'IXJY'

l_op = sqa_plus.term(1.0, [], [cre_i, des_x])
r_op = sqa_plus.term(1.0, [], [cre_y, des_j])

print("\n## Left operator terms:")
print(l_op)

print("\n## Right operator terms:")
print(r_op)

# Spin-Orbital H_eff
terms_Heff = []
## Spin-Orbital Dyall H
print("\n## Creating spin-orbital H^(0)...")
terms_Dyall_H = sqa_plus.dyallH()
print("{:} spin-orbital H^(0) terms created.".format(len(terms_Dyall_H)))

## Spin-Orbital T - T^\dag
print("\n## Creating first set of spin-orbital T - T^\dag...")
terms_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-orbital A^(1) terms created.".format(len(terms_T)))

print("\n## Creating second set of spin-orbital T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-orbital A^(1) terms created.".format(len(terms2_T)))

## Solving for T in slices
print("\n## Selecting spin-orbital T - T^\dag slice: {:}...".format(selected_amplitude))
selected_terms_T = []
selected_terms2_T = []

if selected_amplitude == 't1_0p':
    ind_list_select = [[sqa_plus.options.cvs_core_type, sqa_plus.options.virtual_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.virtual_type],
                       [sqa_plus.options.cvs_core_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't1_m1p':
    ind_list_select = [[sqa_plus.options.active_type, sqa_plus.options.virtual_type],
                       [sqa_plus.options.active_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't1_p1p':
    ind_list_select = [[sqa_plus.options.cvs_core_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_core_type, sqa_plus.options.active_type, sqa_plus.options.active_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.active_type, sqa_plus.options.active_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't1_m2':
     ind_list_select = [[sqa_plus.options.active_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

elif selected_amplitude == 't1_m1':
    ind_list_select = [[sqa_plus.options.cvs_core_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.active_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

elif selected_amplitude == 't1_p2':
    ind_list_select = [[sqa_plus.options.cvs_core_type, sqa_plus.options.cvs_core_type, sqa_plus.options.active_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_core_type, sqa_plus.options.cvs_valence_type, sqa_plus.options.active_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.cvs_core_type, sqa_plus.options.active_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.cvs_valence_type, sqa_plus.options.active_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't1_p1':
    ind_list_select = [[sqa_plus.options.cvs_core_type, sqa_plus.options.cvs_core_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_core_type, sqa_plus.options.cvs_valence_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.cvs_core_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.cvs_valence_type, sqa_plus.options.virtual_type, sqa_plus.options.active_type]]

elif selected_amplitude == 't1_0':
    ind_list_select = [[sqa_plus.options.cvs_core_type, sqa_plus.options.cvs_core_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type],
                       [sqa_plus.options.cvs_core_type, sqa_plus.options.cvs_valence_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.cvs_core_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type],
                       [sqa_plus.options.cvs_valence_type, sqa_plus.options.cvs_valence_type, sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]

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
print(">>> Slice {:}: {:} spin-orbital A terms selected.".format(selected_amplitude, len(selected_terms_T)))

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
print(">>> Slice {:}: {:} spin-orbital A terms selected.".format(selected_amplitude, len(selected_terms2_T)))

# Second Commutator in H2
print("\n# Creating spin-orbital 1/2 * [V + H^(1), T - T^\dag]:")

# Spin-Orbital V
print("\n## Creating spin-orbital V...")
terms_V = sqa_plus.Vperturbation()
for term_V in terms_V:
    term_V.scale(2.0)
print("{:} spin-orbital V terms created.".format(len(terms_V)))

## Commutator [H0, T - T^\dag]
print("\n## Calculating [H^(0), T - T^\dag]:")
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, selected_terms_T)
print(">>> Slice {:}: {:} spin-orbital terms created.".format(selected_amplitude, len(terms_commutator_H0_T)))

## Term V + H1
print("\n## Calculating V + H^(1):")
terms_V_H1 = terms_V + terms_commutator_H0_T
print(">>> Slice {:}: {:} spin-orbital terms created.".format(selected_amplitude, len(terms_V_H1)))

## Commutator [V + H1, T - T^\dag]
print("\n## Calculating 1/2 * [V + H^(1), T - T^\dag]:")
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, selected_terms2_T)
for _term in terms_commutator_V_H1_T:
    _term.scale(0.5)
print(">>> Slice {:}: {:} spin-orbital terms created.".format(selected_amplitude, len(terms_commutator_H0_T)))

print("\n## Calculating H^(2):")
terms_Heff.extend(terms_commutator_V_H1_T)
del(terms_commutator_V_H1_T, terms_T, terms2_T)
print(">>> Slice {:}: {:} spin-orbital H^(2) terms created.".format(selected_amplitude, len(terms_Heff)))

sqa_plus.options.print_divider()

## Calculate the commutators
print("\n## Calculating [H^(2), a_Y^\dag a_J] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, r_op)

print("\n## Calculating [a_I^\dag a_X, [H^(2), a_Y^\dag a_J]] ...")
terms_EE_M0 = sqa_plus.commutator(l_op, terms_commutator)

# Expected value of Spin-Orbital EE M2 CA-CA
expected_EE_M0 = sqa_plus.matrixBlock(terms_EE_M0)

# Generate Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M0, 'temp', final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
