import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 0

import time
start = time.time()

indices_string = 'ca_caaa'
spin_indices_string = 'abba'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: Sigma H{:} {:} ({:})".format(order_Heff, indices_string.upper(), spin_indices_string))

# Generating operators
print("\n## Generating operators ...\n")

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices
b_alpha = sqa_plus.index('B', [tg_alpha, tg_vir])
b_beta  = sqa_plus.index('B', [tg_beta,  tg_vir])

c_alpha = sqa_plus.index('C', [tg_alpha, tg_vir])
c_beta  = sqa_plus.index('C', [tg_beta,  tg_vir])

d_alpha = sqa_plus.index('D', [tg_alpha, tg_vir])
d_beta  = sqa_plus.index('D', [tg_beta,  tg_vir])

k_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_cor])

k_val_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_val])
k_val_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_val])

l_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_cor])
l_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_cor])

l_val_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_val])
l_val_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_val])

w_alpha = sqa_plus.index('W', [tg_alpha, tg_act])
w_beta  = sqa_plus.index('W', [tg_beta,  tg_act])

y_alpha = sqa_plus.index('Y', [tg_alpha, tg_act])
y_beta  = sqa_plus.index('Y', [tg_beta,  tg_act])

z_alpha = sqa_plus.index('Z', [tg_alpha, tg_act])
z_beta  = sqa_plus.index('Z', [tg_beta,  tg_act])

## Dummy Indices
dummy = True

a_alpha = sqa_plus.index('a', [tg_alpha, tg_vir], dummy)
a_beta  = sqa_plus.index('a', [tg_beta,  tg_vir], dummy)

i_alpha = sqa_plus.index('i', [tg_alpha, tg_cvs_cor], dummy)
i_beta  = sqa_plus.index('i', [tg_beta,  tg_cvs_cor], dummy)

x_alpha = sqa_plus.index('x', [tg_alpha, tg_act], dummy)
x_beta  = sqa_plus.index('x', [tg_beta,  tg_act], dummy)

## Define terms
Xsym_1 = sqa_plus.symmetry((1,0,2,3),-1) # Xpqrs = Xqprs
Xsym_2 = sqa_plus.symmetry((0,1,3,2),-1) # Xpqrs = Xqpsr

indices_string_left  = indices_string.split('_')[0]
indices_string_right = indices_string.split('_')[1]

#CAAA
if indices_string_right in ['caaa']:
    final_indices_string = 'KWYZ'

    if spin_indices_string == 'aaaa':
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'abba':
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'baba':
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_beta)])

#CE
if indices_string_left in ['ce']:
    X_tensor_aa = [sqa_plus.tensor('X_aa', [i_alpha, a_alpha], [])]
    X_tensor_bb = [sqa_plus.tensor('X_aa', [i_beta,  a_beta], [])]

    terms_left = [sqa_plus.term(1.0, [], X_tensor_aa + [sqa_plus.creOp(i_alpha), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_bb + [sqa_plus.creOp(i_beta), sqa_plus.desOp(a_beta)])] 

#CA
elif indices_string_left in ['ca']:
    X_tensor_aa = [sqa_plus.tensor('X_aa', [i_alpha, x_alpha], [])]
    X_tensor_bb = [sqa_plus.tensor('X_aa', [i_beta,  x_beta], [])]

    terms_left = [sqa_plus.term(1.0, [], X_tensor_aa + [sqa_plus.creOp(i_alpha), sqa_plus.desOp(x_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_bb + [sqa_plus.creOp(i_beta), sqa_plus.desOp(x_beta)])] 


indices_string = final_indices_string + '_' + spin_indices_string

print("\n## Right operator terms:")
print(term_right)

print("\n## Left operator terms:")
for _term in terms_left:
    print(_term)

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculating the commutator
print("## Calculating the commutator [H(1), h(1)^\dag] ...")
##combine=True?
terms_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating [h^(0), [H(1), h(1)^\dag]] ...")
terms_EE_M01 = sqa_plus.commutator(terms_left, terms_commutator)

# Expected value of Spin-Adapted EE M11
expected_EE_M01 = sqa_plus.matrixBlock(terms_EE_M01)

################
### the following skips spin adaptation, and only uses spin integration!
##expected_EE_M11.sort()
##
### Generating Numpy einsum equations
##result = sqa_plus.genEinsum(expected_EE_M11, 'sigma_' + indices_string, final_indices_string, use_cvs_tensors = True, rm_core_int = True)
##
##end = time.time()
##print("> Total elapsed time: {:.2f} seconds.".format(end - start))
##exit()
#################

# Spin-Adaptation of EE M01
from custom_si_to_sa_functions import convert_X_si_to_sa
expected_EE_M01 = convert_X_si_to_sa(expected_EE_M01, sqa_plus.options)

expected_EE_M01_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M01)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M01_sa, 'sigma_' + indices_string, final_indices_string, use_cvs_tensors = True, rm_core_int = True)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
