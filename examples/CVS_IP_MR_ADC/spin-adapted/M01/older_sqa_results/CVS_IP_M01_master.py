import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 1

import time
start = time.time()

indices_string = 'c_caa'
spin_indices_string = 'a_bab'

sqa_plus.options.print_header("Spin-Adapted CVS-IP: M01")

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
a_alpha = sqa_plus.index('A', [tg_alpha, tg_vir])
a_beta  = sqa_plus.index('A', [tg_beta,  tg_vir])

b_alpha = sqa_plus.index('B', [tg_alpha, tg_vir])
b_beta  = sqa_plus.index('B', [tg_beta,  tg_vir])

i_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_cor])
i_beta  = sqa_plus.index('I', [tg_beta,  tg_cvs_cor])

j_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_cor])
j_beta  = sqa_plus.index('J', [tg_beta,  tg_cvs_cor])

k_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_cor])

l_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_cor])
l_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_cor])

l_val_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_val])
l_val_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_val])

x_alpha = sqa_plus.index('X', [tg_alpha, tg_act])
x_beta  = sqa_plus.index('X', [tg_beta,  tg_act])

y_alpha = sqa_plus.index('Y', [tg_alpha, tg_act])
y_beta  = sqa_plus.index('Y', [tg_beta,  tg_act])

w_alpha = sqa_plus.index('W', [tg_alpha, tg_act])
w_beta  = sqa_plus.index('W', [tg_beta,  tg_act])

z_alpha = sqa_plus.index('Z', [tg_alpha, tg_act])
z_beta  = sqa_plus.index('Z', [tg_beta,  tg_act])

## Define terms
if indices_string in ['c_caa']:
    if spin_indices_string == 'a_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_alpha)])

    elif spin_indices_string == 'a_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(z_beta)])

    elif spin_indices_string == 'a_bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_beta)])

    final_indices_string = 'IKWZ'
    indices_string = indices_string + '_' + spin_indices_string

elif indices_string in ['c_cce']:
    if spin_indices_string == 'a_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IKLB'
    indices_string = indices_string + '_' + spin_indices_string

elif indices_string in ['c_cve']:
    if spin_indices_string == 'a_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IKLB'
    indices_string = indices_string + '_' + spin_indices_string

elif indices_string in ['c_cca']:
    if spin_indices_string == 'a_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'IKLY'
    indices_string = indices_string + '_' + spin_indices_string

elif indices_string in ['c_cva']:
    if spin_indices_string == 'a_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'IKLY'
    indices_string = indices_string + '_' + spin_indices_string

elif indices_string in ['c_cae']:
    if spin_indices_string == 'a_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IKYB'
    indices_string = indices_string + '_' + spin_indices_string

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculating the commutator
print("## Calculating the commutator [H(0), a_S^\dag a_T^\dag a_U] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating a_Q [H(0), a_S^\dag a_T^\dag a_U] ...")
terms_IP_M01 = []
for term_commutator in terms_commutator:
    terms_IP_M01.append(sqa_plus.multiplyTerms(term_commutator, term_left))
    terms_IP_M01.append(sqa_plus.multiplyTerms(term_left, term_commutator))

# Expected value of Spin-Adapted IP M01
expected_IP_M01 = sqa_plus.matrixBlock(terms_IP_M01)

# Spin-Adaptation of IP M01
expected_IP_M01_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_M01)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_IP_M01_sa, 'M_' + indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
