import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 0

import time
start = time.time()

indices_string = 'cae'
spin_indices_string = 'bab'

sqa_plus.options.print_header("Spin-Adapted CVS-IP: Sigma H0 {:} ({:})".format(indices_string.upper(), spin_indices_string))

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

z_alpha = sqa_plus.index('Z', [tg_alpha, tg_act])
z_beta  = sqa_plus.index('Z', [tg_beta,  tg_act])

## Dummy Indices
dummy = True

a_alpha = sqa_plus.index('a', [tg_alpha, tg_vir], dummy)
a_beta  = sqa_plus.index('a', [tg_beta,  tg_vir], dummy)

i_alpha = sqa_plus.index('i', [tg_alpha, tg_cvs_cor], dummy)
i_beta  = sqa_plus.index('i', [tg_beta,  tg_cvs_cor], dummy)

i_val_alpha = sqa_plus.index('i', [tg_alpha, tg_cvs_val], dummy)
i_val_beta  = sqa_plus.index('i', [tg_beta,  tg_cvs_val], dummy)

j_alpha = sqa_plus.index('j', [tg_alpha, tg_cvs_cor], dummy)
j_beta  = sqa_plus.index('j', [tg_beta,  tg_cvs_cor], dummy)

j_val_alpha = sqa_plus.index('j', [tg_alpha, tg_cvs_val], dummy)
j_val_beta  = sqa_plus.index('j', [tg_beta,  tg_cvs_val], dummy)

x_alpha = sqa_plus.index('x', [tg_alpha, tg_act], dummy)
x_beta  = sqa_plus.index('x', [tg_beta,  tg_act], dummy)

y_alpha = sqa_plus.index('y', [tg_alpha, tg_act], dummy)
y_beta  = sqa_plus.index('y', [tg_beta,  tg_act], dummy)

## Define terms
### * <- CAA
symmetry_X_ppq = [sqa_plus.symmetry((1,0,2), -1)]

if indices_string in ['caa']:
    X_tensor_aaa = [sqa_plus.tensor('X_aaa', [i_alpha, x_alpha, y_alpha], [])]
    X_tensor_abb = [sqa_plus.tensor('X_abb', [i_alpha, x_beta,  y_beta],  [])]
    X_tensor_bab = [sqa_plus.tensor('X_bab', [i_beta,  x_alpha, y_beta],  [])]

    terms_right = [sqa_plus.term(1.0, [], X_tensor_aaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abb + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta),  sqa_plus.desOp(y_beta)]),
                   sqa_plus.term(1.0, [], X_tensor_bab + [sqa_plus.creOp(i_beta),  sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_beta)])]

    if spin_indices_string == 'aaa':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(z_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'abb':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(z_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'bab':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(z_beta), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_beta)])

    final_indices_string = 'KWZ'

elif indices_string in ['cca']:
    X_tensor_aaa = [sqa_plus.tensor('X_aaa', [i_alpha, j_alpha, x_alpha], symmetry_X_ppq)]
    X_tensor_abb = [sqa_plus.tensor('X_abb', [i_alpha, j_beta,  x_beta],  [])]
    X_tensor_bab = [sqa_plus.tensor('X_bab', [i_beta,  j_alpha, x_beta],  [])]

    terms_right = [sqa_plus.term(0.5, [], X_tensor_aaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(x_alpha)]),
                   sqa_plus.term(0.5, [], X_tensor_abb + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta),  sqa_plus.desOp(x_beta)]),
                   sqa_plus.term(0.5, [], X_tensor_bab + [sqa_plus.creOp(i_beta),  sqa_plus.creOp(j_alpha), sqa_plus.desOp(x_beta)])]

    if spin_indices_string == 'aaa':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'abb':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_beta), sqa_plus.desOp(l_beta), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'bab':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_beta), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_beta)])

    final_indices_string = 'KLW'

elif indices_string in ['cva']:
    X_tensor_aaa = [sqa_plus.tensor('X_aaa', [i_alpha, j_val_alpha, x_alpha], [])]
    X_tensor_abb = [sqa_plus.tensor('X_abb', [i_alpha, j_val_beta,  x_beta],  [])]
    X_tensor_bab = [sqa_plus.tensor('X_bab', [i_beta,  j_val_alpha, x_beta],  [])]

    terms_right = [sqa_plus.term(1.0, [], X_tensor_aaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(x_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abb + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta),  sqa_plus.desOp(x_beta)]),
                   sqa_plus.term(1.0, [], X_tensor_bab + [sqa_plus.creOp(i_beta),  sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(x_beta)])]

    if spin_indices_string == 'aaa':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'abb':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_beta), sqa_plus.desOp(l_val_beta), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'bab':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_beta), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_beta)])

    final_indices_string = 'KLW'

elif indices_string in ['cce']:
    X_tensor_aaa = [sqa_plus.tensor('X_aaa', [i_alpha, j_alpha, a_alpha], symmetry_X_ppq)]
    X_tensor_abb = [sqa_plus.tensor('X_abb', [i_alpha, j_beta,  a_beta],  [])]
    X_tensor_bab = [sqa_plus.tensor('X_bab', [i_beta,  j_alpha, a_beta],  [])]

    terms_right = [sqa_plus.term(0.5, [], X_tensor_aaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(a_alpha)]),
                   sqa_plus.term(0.5, [], X_tensor_abb + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta),  sqa_plus.desOp(a_beta)]),
                   sqa_plus.term(0.5, [], X_tensor_bab + [sqa_plus.creOp(i_beta),  sqa_plus.creOp(j_alpha), sqa_plus.desOp(a_beta)])]

    if spin_indices_string == 'aaa':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'abb':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(l_beta), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'bab':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_beta)])

    final_indices_string = 'KLB'

elif indices_string in ['cve']:
    X_tensor_aaa = [sqa_plus.tensor('X_aaa', [i_alpha, j_val_alpha, a_alpha], [])]
    X_tensor_abb = [sqa_plus.tensor('X_abb', [i_alpha, j_val_beta,  a_beta],  [])]
    X_tensor_bab = [sqa_plus.tensor('X_bab', [i_beta,  j_val_alpha, a_beta],  [])]

    terms_right = [sqa_plus.term(1.0, [], X_tensor_aaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(a_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abb + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta),  sqa_plus.desOp(a_beta)]),
                   sqa_plus.term(1.0, [], X_tensor_bab + [sqa_plus.creOp(i_beta),  sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(a_beta)])]

    if spin_indices_string == 'aaa':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_alpha), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'abb':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(l_val_beta), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'bab':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_beta)])

    final_indices_string = 'KLB'

elif indices_string in ['cae']:
    X_tensor_aaa = [sqa_plus.tensor('X_aaa', [i_alpha, x_alpha, a_alpha], [])]
    X_tensor_abb = [sqa_plus.tensor('X_abb', [i_alpha, x_beta,  a_beta],  [])]
    X_tensor_bab = [sqa_plus.tensor('X_bab', [i_beta,  x_alpha, a_beta],  [])]

    terms_right = [sqa_plus.term(1.0, [], X_tensor_aaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(a_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abb + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta),  sqa_plus.desOp(a_beta)]),
                   sqa_plus.term(1.0, [], X_tensor_bab + [sqa_plus.creOp(i_beta),  sqa_plus.creOp(x_alpha), sqa_plus.desOp(a_beta)])]

    if spin_indices_string == 'aaa':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'abb':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(k_alpha)])

    elif spin_indices_string == 'bab':
        term_left   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(z_alpha), sqa_plus.desOp(k_beta)])

    final_indices_string = 'KZB'

indices_string = indices_string[:3]

print("\n## Right operator terms:")
for _term in terms_right:
    print(_term)

print("\n## Left operator terms:")
print(term_left)

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculating the commutator
print("## Calculating the commutator [H(0), a_S^\dag a_T^\dag a_U] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, terms_right)

print("\n## Calculating a_P^\dag a_Q a_R [H(0), a_S^\dag a_T^\dag a_U] ...")
terms_IP_M11 = []
for term_commutator in terms_commutator:
    terms_IP_M11.append(sqa_plus.multiplyTerms(term_commutator, term_left))
    terms_IP_M11.append(sqa_plus.multiplyTerms(term_left, term_commutator))

# Expected value of Spin-Adapted IP M11
expected_IP_M11 = sqa_plus.matrixBlock(terms_IP_M11)

# Spin-Adaptation of IP M11
from custom_si_to_sa_functions import convert_X_si_to_sa
sqa_plus.options.add_spin_adaptation_custom_function(convert_X_si_to_sa)

expected_IP_M11_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_M11)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_IP_M11_sa, 'sigma_' + indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))

