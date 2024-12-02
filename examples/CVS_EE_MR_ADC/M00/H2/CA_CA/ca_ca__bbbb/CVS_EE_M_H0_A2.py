import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True

import time
start = time.time()

spin_indices_string = 'bb_bb'
selected_amplitude = 't2_0pp'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: M00 H2 CA_CA ({:}) | T-Slice {:}".format(spin_indices_string, selected_amplitude))

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

# Generating Term 
## External Indices and Operators
#Alpha
i_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_cor])
j_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_cor])
x_alpha = sqa_plus.index('X', [tg_alpha, tg_act])
y_alpha = sqa_plus.index('Y', [tg_alpha, tg_act])

#Beta
i_beta = sqa_plus.index('I', [tg_beta, tg_cvs_cor])
j_beta = sqa_plus.index('J', [tg_beta, tg_cvs_cor])
x_beta = sqa_plus.index('X', [tg_beta, tg_act])
y_beta = sqa_plus.index('Y', [tg_beta, tg_act])


## Define operators
#Alpha
cre_i_alpha = sqa_plus.creOp(i_alpha)
des_j_alpha = sqa_plus.desOp(j_alpha)
des_x_alpha = sqa_plus.desOp(x_alpha)
cre_y_alpha = sqa_plus.creOp(y_alpha)

#Beta
cre_i_beta = sqa_plus.creOp(i_beta)
des_j_beta = sqa_plus.desOp(j_beta)
des_x_beta = sqa_plus.desOp(x_beta)
cre_y_beta = sqa_plus.creOp(y_beta)

## Define terms
print("\n## Generating Term a_I^\dag a_X ... a_Y^\dag a_J ...")
final_indices_string = 'IXJY'

## aa - aa
#l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_x_alpha])
#r_op = sqa_plus.term(1.0, [], [cre_y_alpha, des_j_alpha])

## bb - bb
l_op = sqa_plus.term(1.0, [], [cre_i_beta, des_x_beta])
r_op = sqa_plus.term(1.0, [], [cre_y_beta, des_j_beta])

print("\n## Left operator terms:")
print(l_op)

print("\n## Right operator terms:")
print(r_op)

# Spin-Adapted H_eff
terms_Heff = []
## Spin-Integrated Dyall H
print("\n## Creating spin-integrated H^(0)...")
terms_Dyall_H = sqa_plus.dyallH()
print("{:} spin-integrated H^(0) terms created.".format(len(terms_Dyall_H)))

# First Commutator in H2
## Spin-Integrated T - T^\dag
print("\n## Creating set of spin-integrated T^{(2)} - T^{(2)\dag}...")
terms_T = sqa_plus.Tamplitude(order = 2)
print("{:} spin-integrated A^(2) terms created.".format(len(terms_T)))

## Select t2_0pp from A(2) set:
print("\n## Selecting spin-integrated T^{(2)} - T{(2)^\dag} slice: t2_0pp...")
ind_list_select = [sqa_plus.options.active_type, sqa_plus.options.active_type]
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

for term_T in selected_terms_T:
    print(term_T)
print(">>> Slice {:}: {:} spin-integrated A terms selected.".format(selected_amplitude, len(selected_terms_T)))

## Commutator [H0, T2 - T^2\dag]
print("\n## Calculating [H^(0), T^{(2)} - T^{(2)\dag}]:")
terms_commutator_H0_T2 = sqa_plus.commutator(terms_Dyall_H, selected_terms_T)

terms_Heff.extend(terms_commutator_H0_T2)
del(terms_commutator_H0_T2, terms_T)

print(">>> Slice {:}: {:} spin-integrated H^(2) terms created.".format(selected_amplitude, len(terms_Heff)))

sqa_plus.options.print_divider()

## Calculate the commutators
print("\n## Calculating [H^(2), a_Y^\dag a_J] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, r_op)

print("\n## Calculating [a_I^\dag a_X, [H^(2), a_Y^\dag a_J]] ...")
terms_EE_M0 = sqa_plus.commutator(l_op, terms_commutator)

# Expected value of Spin-Integrated EE M2 CA-CA
expected_EE_M0 = sqa_plus.matrixBlock(terms_EE_M0)

# Spin-Adaptation of EE M2 CA-CA
expected_EE_M0_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M0)

# Generate Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M0_sa, 'temp', final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
