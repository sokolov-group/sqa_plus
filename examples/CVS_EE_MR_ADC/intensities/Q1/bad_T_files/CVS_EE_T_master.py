import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
q_order = 1

import time
start = time.time()

indices_string = 'q_ce'
spin_indices_string = 'aa_aa'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: T Q{:} {:} ({:})".format(q_order, indices_string.upper(), spin_indices_string))

# Generating operators
print("\n## Generating operators ...\n")

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

dummy = True

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
term_q0 = []
dSym = [sqa_plus.symmetry((1,0),1)]

# CC
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_cor], dummy)

dTen = sqa_plus.tensor('d_cc', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
#dTen = sqa_plus.tensor('d_cc', [p_alpha, q_alpha], dSym)
#term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# VV
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_val], dummy)

dTen = sqa_plus.tensor('d_vv', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))

# CV
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_val], dummy)

dTen = sqa_plus.tensor('d_cv', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_cv', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# AA
p_alpha = sqa_plus.index('P', [tg_alpha, tg_act], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act], dummy)

dTen = sqa_plus.tensor('d_aa', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
#dTen = sqa_plus.tensor('d_aa', [p_alpha, q_alpha], dSym)
#term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# EE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_vir], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir], dummy)

dTen = sqa_plus.tensor('d_ee', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
#dTen = sqa_plus.tensor('d_ee', [p_alpha, q_alpha], dSym)
#term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# CA
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act], dummy)

dTen = sqa_plus.tensor('d_ca', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ca', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# VA
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act], dummy)

dTen = sqa_plus.tensor('d_va', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_va', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# CE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir], dummy)

dTen = sqa_plus.tensor('d_ce', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ce', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# VE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir], dummy)

dTen = sqa_plus.tensor('d_ve', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ve', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# AE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_act], dummy)
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir], dummy)

dTen = sqa_plus.tensor('d_ae', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ae', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

if indices_string in ['q_ca']:

    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(x_alpha), sqa_plus.desOp(i_alpha)])
    
    final_indices_string = 'IX'

elif indices_string in ['q_ce']:

    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.desOp(i_alpha)])
    
    final_indices_string = 'IA'

#elif indices_string in ['c_c']:
#    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])
#
#    final_indices_string = 'IJ'
#
#
#elif indices_string in ['v_c']:
#    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])
#
#    final_indices_string = 'IJ'
#
#elif indices_string in ['a_c']:
#    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])
#
#    final_indices_string = 'XJ'
#
#elif indices_string in ['e_c']:
#    term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#    term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha)])
#
#    final_indices_string = 'AJ'
#
#elif indices_string in ['c_caa']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])
#
#    final_indices_string = 'IJYZ'
#
#elif indices_string in ['v_caa']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])
#
#    final_indices_string = 'IJYZ'
#
#elif indices_string in ['a_caa']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])
#
#    final_indices_string = 'XJYZ'
#
#elif indices_string in ['e_caa']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(z_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(z_beta)])
#
#    final_indices_string = 'AJYZ'
#
#elif indices_string in ['c_cce']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'IJKB'
#
#elif indices_string in ['v_cce']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'IJKB'
#
#elif indices_string in ['c_cve']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'IJKB'
#
#elif indices_string in ['v_cve']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'IJKB'
#
#elif indices_string in ['a_cce']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'XJKB'
#
#elif indices_string in ['a_cve']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'XJKB'
#
#elif indices_string in ['e_cce']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'AJKB'
#
#elif indices_string in ['e_cve']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'AJKB'
#
#elif indices_string in ['c_cae']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'IJYB'
#
#elif indices_string in ['v_cae']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'IJYB'
#
#elif indices_string in ['a_cae']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'XJYB'
#
#elif indices_string in ['e_cae']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(b_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(b_beta)])
#
#    final_indices_string = 'AJYB'
#
#elif indices_string in ['c_cca']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'IJKY'
#
#elif indices_string in ['v_cca']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'IJKY'
#
#elif indices_string in ['c_cva']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'IJKY'
#
#elif indices_string in ['v_cva']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(i_val_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'IJKY'
#
#elif indices_string in ['a_cca']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'XJKY'
#
#elif indices_string in ['a_cva']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(x_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'XJKY'
#
#elif indices_string in ['e_cca']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'AJKY'
#
#elif indices_string in ['e_cva']:
#    if spin_indices_string == 'a_aaa':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_alpha)])
#
#    elif spin_indices_string == 'a_abb':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_alpha), sqa_plus.creOp(k_val_beta), sqa_plus.desOp(y_beta)])
#
#    elif spin_indices_string == 'a_bab':
#        term_q0 = sqa_plus.term(1.0, [], [sqa_plus.desOp(a_alpha)])
#        term_h  = sqa_plus.term(1.0, [], [sqa_plus.creOp(j_beta), sqa_plus.creOp(k_val_alpha), sqa_plus.desOp(y_beta)])
#
#    final_indices_string = 'AJKY'

## Create Dummy Indices List
if q_order == 0:
    print("## Calculating the excitation operator q^(0) ...")
    # Zeroth-order excitation operator
    terms_q = term_q0

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

terms_EE_T = sqa_plus.commutator(terms_q, term_h)

# Expected value of Spin-Adapted EE T
expected_EE_T = sqa_plus.matrixBlock(terms_EE_T)

# Spin-Adaptation of EE T
expected_EE_T_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_T)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_T_sa, 'T_' + indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
