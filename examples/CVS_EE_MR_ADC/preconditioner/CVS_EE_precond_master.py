import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True

import time
start = time.time()

diagonal_indices_string = 'caea'
spin_indices_string = 'aaaa_abab'

sqa_plus.options.print_header("Spin-Adapted Preconditioner {:} {:}".format(diagonal_indices_string, spin_indices_string))

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

c_alpha = sqa_plus.index('C', [tg_alpha, tg_vir])
c_beta  = sqa_plus.index('C', [tg_beta,  tg_vir])

d_alpha = sqa_plus.index('D', [tg_alpha, tg_vir])
d_beta  = sqa_plus.index('D', [tg_beta,  tg_vir])

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

u_alpha = sqa_plus.index('U', [tg_alpha, tg_act])
u_beta  = sqa_plus.index('U', [tg_beta,  tg_act])

v_alpha = sqa_plus.index('V', [tg_alpha, tg_act])
v_beta  = sqa_plus.index('V', [tg_beta,  tg_act])

## Define terms
if diagonal_indices_string in ['ccee']:
    if spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)]) #LHS = KLCD
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_beta), sqa_plus.desOp(l_beta), sqa_plus.desOp(k_alpha)]) #RHS = IJAB

    final_indices_string = 'IJAB'
    diagonal_pairs_dict = {'L': 'J', 'K': 'I', 'C': 'A', 'D':'B'}

elif diagonal_indices_string in ['cvee']:

    if spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)]) #LHS = KLCD
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_beta), sqa_plus.desOp(l_val_beta), sqa_plus.desOp(k_alpha)]) #RHS = IJAB

    final_indices_string = 'IJAB'
    diagonal_pairs_dict = {'L': 'J', 'K': 'I', 'C': 'A', 'D':'B'}

elif diagonal_indices_string in ['ccaa']:
    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_alpha)]) #LHS = KLWZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)]) #RHS = IJXY

    elif spin_indices_string == 'abba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]) #LHS = KLWZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(l_beta), sqa_plus.desOp(k_alpha)]) #RHS = IJXY

    elif spin_indices_string == 'baba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]) #LHS = KLWZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_beta)]) #RHS = IJXY

    final_indices_string = 'IJXYWZ'
    diagonal_pairs_dict = {'K': 'I', 'L': 'J'}

elif diagonal_indices_string in ['cvaa']:
    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_alpha)]) #LHS = KLWZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_alpha)]) #RHS = IJXY

    elif spin_indices_string == 'abba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]) #LHS = KLWZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(l_val_beta), sqa_plus.desOp(k_alpha)]) #RHS = IJXY

    elif spin_indices_string == 'baba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]) #LHS = KLWZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_beta)]) #RHS = IJXY

    final_indices_string = 'IJXYWZ'
    diagonal_pairs_dict = {'K': 'I', 'L': 'J'}

elif diagonal_indices_string in ['ccea']:
    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(b_alpha)]) #LHS = KLBW
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IJAX

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(b_alpha)]) #LHS = KLBW
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(j_beta), sqa_plus.desOp(i_alpha)]) #RHS = IJAX
   
    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(w_beta), sqa_plus.desOp(b_alpha)]) #LHS = KLBW
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_beta)]) #RHS = IJAX

    final_indices_string = 'IJAXW'
    diagonal_pairs_dict = {'B': 'A', 'L': 'J', 'K': 'I'}

elif diagonal_indices_string in ['cvea']:
    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(b_alpha)]) #LHS = KLBW
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(j_val_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IJAX

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(a_alpha)]) #LHS = IJAX
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(l_val_beta), sqa_plus.desOp(k_alpha)]) #RHS = KLBW
   
    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(w_beta), sqa_plus.desOp(b_alpha)]) #LHS = KLBW
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(j_val_alpha), sqa_plus.desOp(i_beta)]) #RHS = IJAX

    final_indices_string = 'IJAXW'
    diagonal_pairs_dict = {'B': 'A', 'L': 'J', 'K': 'I'}

elif diagonal_indices_string in ['caee']:
    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(d_alpha), sqa_plus.desOp(c_alpha)]) #LHS = KWCD
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IXAB

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)]) 
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(k_alpha)]) 

    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)]) 
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)]) 

    final_indices_string = 'IABXW'
    diagonal_pairs_dict = {'K': 'I', 'C': 'A', 'D': 'B'}

elif diagonal_indices_string in ['caea']:
    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(b_alpha)]) #LHS = KWBZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IXAY

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(b_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])

    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_beta), sqa_plus.desOp(b_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)])

    elif spin_indices_string == 'aaaa_abab': ##equal to aaaa_abab
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(b_alpha)]) #LHS = KWBZ
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])

    final_indices_string = 'IAXYWZ'
    diagonal_pairs_dict = {'K': 'I', 'B': 'A'}

elif diagonal_indices_string in ['ce_caea']:
    if spin_indices_string == 'aa_aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)]) #LHS = KB
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IXAY

    elif spin_indices_string == 'aa_abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])

    elif spin_indices_string == 'aa_baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)])

    elif spin_indices_string == 'bb_aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)]) #LHS = KB
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IXAY

    elif spin_indices_string == 'bb_abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])

    elif spin_indices_string == 'bb_baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.desOp(b_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)])

    elif spin_indices_string == 'aa_abba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(b_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(a_beta), sqa_plus.creOp(y_alpha), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])

    final_indices_string = 'IXAY'
    diagonal_indices_string = diagonal_indices_string + '__' + spin_indices_string
    diagonal_pairs_dict = {'K': 'I', 'B': 'A'}

elif diagonal_indices_string in ['ca_caaa']:
    if spin_indices_string == 'aa_aaaa': ##aa - aaaa
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(w_alpha)]) #LHS = KW
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IXYZ

    elif spin_indices_string == 'aa_abba': ##aa - abab
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(w_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])

    elif spin_indices_string == 'aa_baba': ##aa - baab
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(w_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)])

#    elif spin_indices_string == 'aa_abab':
#        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(w_alpha)])
#        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])
#
#    elif spin_indices_string == 'aa_baab':
#        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.desOp(w_alpha)])
#        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)])

    final_indices_string = 'IXYZW'
    diagonal_indices_string = diagonal_indices_string + '__' + spin_indices_string
    diagonal_pairs_dict = {'I': 'K'}

elif diagonal_indices_string in ['caaa']:
    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(v_alpha), sqa_plus.desOp(u_alpha)]) #LHS = KWUV
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IXYZ

    elif spin_indices_string == 'abba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(v_beta), sqa_plus.desOp(u_alpha)]) #LHS = KWUV
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)]) #RHS = IXYZ

    elif spin_indices_string == 'baba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(v_beta), sqa_plus.desOp(u_alpha)]) #LHS = KWUV
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)]) #RHS = IXYZ

    elif spin_indices_string == 'aaaa_abba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(v_alpha), sqa_plus.desOp(u_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)]) #RHS = IXYZ

    elif spin_indices_string == 'abba_aaaa':
         term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(v_beta), sqa_plus.desOp(u_alpha)]) #LHS = KWUV
         term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]) #RHS = IXYZ

    elif spin_indices_string == 'baba_abba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(v_beta), sqa_plus.desOp(u_alpha)]) #LHS = KWUV
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)]) #RHS = IXYZ

    final_indices_string = 'IWVUXZY'
    diagonal_pairs_dict = {'K': 'I'}
    diagonal_indices_string = diagonal_indices_string + '__' + spin_indices_string

# Spin-Integrated H_eff
if 'ce_caea' in diagonal_indices_string or 'ca_caaa' in diagonal_indices_string:
    terms_Heff = sqa_plus.Heff(0)
    terms_Heff.extend(sqa_plus.Heff(1))
else:
    terms_Heff = sqa_plus.Heff(0)

## Calculating the commutator
print("## Calculating the commutator [H(0), h_NU^\dag] ...")
first_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating [h_MU, [H(0), h_NU^\dag]] ...")
terms_EE_M11 = sqa_plus.commutator(term_left, first_commutator)

# Expected value of Spin-Adapted EE M11
expected_EE_M11 = sqa_plus.matrixBlock(terms_EE_M11)

# Spin-Adaptation of EE M11
expected_EE_M11_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M11)

from sqa_plus.sqaEinsum import remove_core_int
expected_EE_M11_sa, removed_core = remove_core_int(expected_EE_M11_sa)
sqa_plus.options.genEinsum.remove_core_integrals = False
sqa_plus.options.genEinsum.keep_user_defined_dummy_names = True

# Replacing indices to obtain the diagonal
for term_ind, term_sa in enumerate(expected_EE_M11_sa):
    for tensor_ind, tensor_sa in enumerate(term_sa.tensors):
        for index_ind, index_sa in enumerate(tensor_sa.indices):
            if index_sa.name in diagonal_pairs_dict.keys():
                expected_EE_M11_sa[term_ind].tensors[tensor_ind].indices[index_ind].userDefined = True
                expected_EE_M11_sa[term_ind].tensors[tensor_ind].indices[index_ind].name = diagonal_pairs_dict[index_sa.name]

from sqa_plus.sqaMatrixBlock import contractDeltaFuncs_nondummy
expected_EE_M11_sa = contractDeltaFuncs_nondummy(expected_EE_M11_sa)

sqa_plus.combineTerms(expected_EE_M11_sa)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M11_sa, 'precond_' + diagonal_indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
