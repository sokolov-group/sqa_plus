import sqa_plus
from sqa_plus.sqaIndex import get_spatial_index_type
sqa_plus.options.spin_orbital = True
sqa_plus.options.cvs_approach = True

import time
start = time.time()

selected_amplitude = 't2_0pp'

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

# First Commutator in H2
## Spin-Orbital T - T^\dag
print("\n## Creating set of spin-orbital T^{(2)} - T^{(2)\dag}...")
terms_T = sqa_plus.Tamplitude(order = 2)
print("{:} spin-orbital A^(2) terms created.".format(len(terms_T)))

## Select t2_0pp from A(2) set:
print("\n## Selecting spin-orbital T^{(2)} - T{(2)^\dag} slice: t2_0pp...")
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
print(">>> Slice {:}: {:} spin-orbital A terms selected.".format(selected_amplitude, len(selected_terms_T)))

## Commutator [H0, T2 - T^2\dag]
print("\n## Calculating [H^(0), T^{(2)} - T^{(2)\dag}]:")
terms_commutator_H0_T2 = sqa_plus.commutator(terms_Dyall_H, selected_terms_T)

terms_Heff.extend(terms_commutator_H0_T2)
del(terms_commutator_H0_T2, terms_T)

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
