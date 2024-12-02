import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_orbital = True

import time
start = time.time()

selected_amplitude = ''

sqa_plus.options.print_header("Spin-Orbital V(I,X): T-Slice {:}".format(selected_amplitude))

# Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

## External Indices and Operators
print("# Defining a^\dag_p a_q...")
i = sqa_plus.index('I', [tg_cor])
x = sqa_plus.index('X', [tg_act])

term_l_op = sqa_plus.term(-0.5, [], [sqa_plus.creOp(i), sqa_plus.desOp(x)])

## Spin-Integrated T - T^\dag
print("\n## Creating first set of spin-orbital T - T^\dag...")
terms_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-orbital A terms created.".format(len(terms_T)))

print("\n## Creating second set of spin-orbital T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(order = 1)
print("{:} spin-orbital A terms created.".format(len(terms2_T)))

#### Solving for T in slices
##print("\n## Selecting spin-orbital T - T^\dag slice: {:}...".format(selected_amplitude))
##selected_terms_T = []
##
##if selected_amplitude == 't1_0p':
##    ind_t1 = [sqa_plus.options.core_type, sqa_plus.options.virtual_type]
##    ind_t2 = [sqa_plus.options.core_type, sqa_plus.options.active_type,
##              sqa_plus.options.virtual_type, sqa_plus.options.active_type]
##    ind_list_select = [ind_t1, ind_t2]
##
##elif selected_amplitude == 't1_m1p':
##    ind_t1 = [sqa_plus.options.active_type, sqa_plus.options.virtual_type]
##    ind_t2 = [sqa_plus.options.active_type, sqa_plus.options.active_type,
##              sqa_plus.options.virtual_type, sqa_plus.options.active_type]
##    ind_list_select = [ind_t1, ind_t2]
##
##elif selected_amplitude == 't1_p1p':
##    ind_t1 = [sqa_plus.options.core_type, sqa_plus.options.active_type]
##    ind_t2 = [sqa_plus.options.core_type, sqa_plus.options.active_type,
##              sqa_plus.options.active_type, sqa_plus.options.active_type]
##    ind_list_select = [ind_t1, ind_t2]
##
##elif selected_amplitude == 't1_m2':
##     ind_list_select = [[sqa_plus.options.active_type, sqa_plus.options.active_type,
##                         sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]
##
##elif selected_amplitude == 't1_m1':
##    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.active_type,
##                       sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]
##
##elif selected_amplitude == 't1_p2':
##    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type,
##                       sqa_plus.options.active_type, sqa_plus.options.active_type]]
##
##elif selected_amplitude == 't1_p1':
##    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type,
##                       sqa_plus.options.virtual_type, sqa_plus.options.active_type]]
##
##elif selected_amplitude == 't1_0':
##    ind_list_select = [[sqa_plus.options.core_type, sqa_plus.options.core_type,
##                       sqa_plus.options.virtual_type, sqa_plus.options.virtual_type]]
##
##print("\n# First Amplitude Set...")
##for term_T in terms_T:
##    select_term = False
##    for tensor_T in term_T.tensors:
##        ind_types = []
##        for ind_T in tensor_T.indices:
##            ind_types.append(get_spatial_index_type(ind_T.indType))
##        if ind_types in ind_list_select:
##            select_term = True
##    if select_term:
##        selected_terms_T.append(term_T)
##
##for term_T in selected_terms_T:
##    print(term_T)
##print(">>> Slice {:}: {:} spin-orbital A terms selected.".format(selected_amplitude, len(selected_terms_T)))

# First Term
print("\n# Creating spin-orbital [2V + [H^(0), T - T^\dag], T - T^\dag]:")

# Spin-Integrated V
print("\n## Creating spin-orbital V...")
terms_V = sqa_plus.Vperturbation()
for term_V in terms_V:
    term_V.scale(2.0)
print("{:} spin-orbital V terms created.".format(len(terms_V)))

## Spin-Integrated Dyall H
print("\n## Creating spin-orbital H^(0)...")
terms_Dyall_H = sqa_plus.dyallH()
print("{:} spin-orbital H^(0) terms created.".format(len(terms_Dyall_H)))

## Commutator [H0, T - T^\dag]
print("\n## Calculating [H^(0), T - T^\dag]:")
#terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, selected_terms_T) ## outer A
#terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, terms2_T) ## inner A
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, terms_T) ## full
print(">>> Slice {:}: {:} spin-orbital terms created.".format(selected_amplitude, len(terms_commutator_H0_T)))

## Term V + H1
print("\n## Calculating V + H^(1):")
terms_V_H1 = terms_V + terms_commutator_H0_T
print(">>> Slice {:}: {:} spin-orbital terms created.".format(selected_amplitude, len(terms_V_H1)))

## Commutator [V + H1, T - T^\dag]
print("\n## Calculating [V + H^(1), T - T^\dag]:")
#terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, terms2_T) ## outer A
#terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, selected_terms_T) ## inner A
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, terms2_T) ## full
print(">>> Slice {:}: {:} spin-orbital terms created.".format(selected_amplitude, len(terms_commutator_H0_T)))

## Term a^\dag_p a_q [V + H1, T - T^\dag]
print("\n## Calculating a^\dag_p a_q [V + H^(1), T - T^\dag]:")
terms_V1 = []
for term_commutator in terms_commutator_V_H1_T:
    terms_V1.append(sqa_plus.multiplyTerms(term_l_op, term_commutator))
print(">>> Slice {:}: {:} spin-orbital a^\dag_p a_q [V + H^(1), T - T^\dag] terms created.\n".format(selected_amplitude, len(terms_V1)))

sqa_plus.options.print_divider()

## Solving for V
print("\n## Calculating matrix element of {:} slice of V1:".format(selected_amplitude))
## Expected value of Spin-Integrated V1
expected_terms_V1 = sqa_plus.matrixBlock(terms_V1)
del terms_V1
print(">>> Slice V1: {:} spin-orbital a^\dag_p a_q [V + H^(1), A^(1)] terms created.".format(len(expected_terms_V1)))

# Generating Numpy einsum equations
#result = sqa_plus.genEinsum(expected_terms_V1, 'V1__'+selected_amplitude, 'XY', remove_core_integrals = True)
result = sqa_plus.genEinsum(expected_terms_V1, 'V1', 'XY', remove_core_integrals = True)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
