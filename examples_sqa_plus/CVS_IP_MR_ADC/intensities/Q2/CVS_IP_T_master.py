import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
q_order = 2

import time
start = time.time()

indices_string = 'a_cae'
spin_indices_string = 'a_bab'

sqa_plus.options.print_header("Spin-Adapted CVS-IP: T Q{:} {:} ({:})".format(q_order, indices_string.upper(), spin_indices_string))

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

i_val_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_val])
i_val_beta  = sqa_plus.index('I', [tg_beta,  tg_cvs_val])

j_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_cor])
j_beta  = sqa_plus.index('J', [tg_beta,  tg_cvs_cor])

j_val_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_val])
j_val_beta  = sqa_plus.index('J', [tg_beta,  tg_cvs_val])

k_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_cor])

k_val_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_val])
k_val_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_val])

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
if indices_string in ['c_c']:
    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])

    final_indices_string = 'IJ'

elif indices_string in ['v_c']:
    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])

    final_indices_string = 'IJ'

elif indices_string in ['a_c']:
    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])

    final_indices_string = 'XJ'

elif indices_string in ['e_c']:
    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])

    final_indices_string = 'AJ'

elif indices_string in ['c_caa']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])

    final_indices_string = 'IJYZ'

elif indices_string in ['v_caa']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])

    final_indices_string = 'IJYZ'

elif indices_string in ['a_caa']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])

    final_indices_string = 'XJYZ'

elif indices_string in ['e_caa']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])

    final_indices_string = 'AJYZ'

elif indices_string in ['c_cce']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IJKB'

elif indices_string in ['v_cce']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IJKB'

elif indices_string in ['c_cve']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IJKB'

elif indices_string in ['v_cve']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IJKB'

elif indices_string in ['a_cce']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'XJKB'

elif indices_string in ['a_cve']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'XJKB'

elif indices_string in ['e_cce']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'AJKB'

elif indices_string in ['e_cve']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'AJKB'

elif indices_string in ['c_cae']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IJYB'

elif indices_string in ['v_cae']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'IJYB'

elif indices_string in ['a_cae']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'XJYB'

elif indices_string in ['e_cae']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])

    final_indices_string = 'AJYB'

elif indices_string in ['c_cca']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'IJKY'

elif indices_string in ['v_cca']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'IJKY'

elif indices_string in ['c_cva']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'IJKY'

elif indices_string in ['v_cva']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'IJKY'

elif indices_string in ['a_cca']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'XJKY'

elif indices_string in ['a_cva']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'XJKY'

elif indices_string in ['e_cca']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'AJKY'

elif indices_string in ['e_cva']:
    if spin_indices_string == 'a_aaa':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'a_bab':
        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'AJKY'

## Create Dummy Indices List
if q_order == 0:
    print("## Calculating the excitation operator q^(0) ...")
    # Zeroth-order excitation operator
    terms_q = [term_q0]

elif q_order == 1:
    # First-order excitation operator
    print("## Calculating the excitation operator [q^(0), T - T^\dag] ...")
    terms_T = sqa_plus.Tamplitude(1)

    ## Calculating the commutator
    print("## Calculating the commutator...")
    terms_q = sqa_plus.commutator(term_q0, terms_T)

elif q_order == 2:
    # First-order excitation operator
    print("## Calculating the excitation operator [q^(0), T - T^\dag] ...")
    terms_T1 = sqa_plus.Tamplitude(1)

    print("## Calculating the commutator...")
    terms_q = sqa_plus.commutator(term_q0, terms_T1)

    print("## Calculating the excitation operator 1/2 * [[q^(0), T - T^\dag], T - T^\dag] ...")
    terms_T1_2 = sqa_plus.Tamplitude(1)

    print("## Calculating the commutator...")
    terms_q = sqa_plus.commutator(terms_q, terms_T1_2)

    for term_q in terms_q:
      term_q.scale(0.5)

    # First-order excitation operator
    print("## Calculating the excitation operator [q^(0), T^(2) - T^(2)^\dag] ...")
    terms_T2 = sqa_plus.Tamplitude(2)

    print("## Calculating the commutator...")
    terms_q_2 = sqa_plus.commutator(term_q0, terms_T2)

    terms_q.extend(terms_q_2)

print("\n## Calculating h [q, T - T^\dag] ...")
terms_IP_T = []
for term_q in terms_q:
    terms_IP_T.append(sqa_plus.multiplyTerms(term_q, term_h))
    terms_IP_T.append(sqa_plus.multiplyTerms(term_h, term_q))

# Expected value of Spin-Adapted IP T
expected_IP_T = sqa_plus.matrixBlock(terms_IP_T)

# Spin-Adaptation of IP T
expected_IP_T_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_T)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_IP_T_sa, 'T_' + indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
