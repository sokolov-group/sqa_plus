import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 0

import time
start = time.time()

indices_string = 'ccea_ccea'
spin_indices_string = 'abab'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: Sigma H0 {:} ({:})".format(indices_string.upper(), spin_indices_string))

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

u_alpha = sqa_plus.index('U', [tg_alpha, tg_act])
u_beta  = sqa_plus.index('U', [tg_beta,  tg_act])

v_alpha = sqa_plus.index('V', [tg_alpha, tg_act])
v_beta  = sqa_plus.index('V', [tg_beta,  tg_act])

## Dummy Indices
dummy = True

a_alpha = sqa_plus.index('a', [tg_alpha, tg_vir], dummy)
a_beta  = sqa_plus.index('a', [tg_beta,  tg_vir], dummy)

b_alpha = sqa_plus.index('b', [tg_alpha, tg_vir], dummy)
b_beta  = sqa_plus.index('b', [tg_beta,  tg_vir], dummy)

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

z_alpha = sqa_plus.index('z', [tg_alpha, tg_act], dummy)
z_beta  = sqa_plus.index('z', [tg_beta,  tg_act], dummy)

## Define terms
# a^{PQ}_{RS} = a_P^\dag a_Q^\dag a_S a_R

X_sym_1 = [sqa_plus.symmetry((1,0,2,3), -1)] #symmetry: pprs
X_sym_2 = [sqa_plus.symmetry((0,1,3,2), -1)] #symmetry: pqrr
X_sym_3 = X_sym_1 + X_sym_2 #symmetry: pprr

indices_string_left  = indices_string.split('_')[0]
indices_string_right = indices_string.split('_')[1]

## Left Op: h^(1)

if indices_string_left in ['caaa']:
    final_indices_string = 'KWUV'

#### amplitudes, overlap, precond: kuxv, iwyz --> form used here: kwuv, ixyz
# precond  
#    if spin_indices_string == 'aaaa':
#        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(u_alpha), sqa_plus.desOp(v_alpha), sqa_plus.desOp(x_alpha)]) #lhs = iuxv
#        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)]) #rhs = kwyz

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(v_alpha), sqa_plus.desOp(u_alpha)])

    elif spin_indices_string == 'abab':    
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(v_beta), sqa_plus.desOp(u_alpha)]) 

#    elif spin_indices_string == 'baab': 
#        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(v_beta), sqa_plus.desOp(u_alpha)])

    elif spin_indices_string == 'baba': ##baba
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(v_alpha), sqa_plus.desOp(u_beta)])

    elif spin_indices_string == 'bbbb': ##bbbb
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_beta), sqa_plus.desOp(v_beta), sqa_plus.desOp(u_beta)])

elif indices_string_left in ['caea']:
    final_indices_string = 'KWCU'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(u_alpha), sqa_plus.desOp(c_alpha)]) #KWCU

    elif spin_indices_string == 'abab':    
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(u_beta), sqa_plus.desOp(c_alpha)]) 

    elif spin_indices_string == 'baab': 
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(u_beta), sqa_plus.desOp(c_alpha)])
 
#CCAA
elif indices_string_left in ['ccaa']:
    final_indices_string = 'KLWU'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(u_alpha), sqa_plus.desOp(w_alpha)])

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(u_beta), sqa_plus.desOp(w_alpha)]) 

    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(u_beta), sqa_plus.desOp(w_alpha)]) 

#CVAA
elif indices_string_left in ['cvaa']:
    final_indices_string = 'KLWU'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(u_alpha), sqa_plus.desOp(w_alpha)])

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_beta), sqa_plus.desOp(u_beta), sqa_plus.desOp(w_alpha)]) 

    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(u_beta), sqa_plus.desOp(w_alpha)]) 

#CCEA
elif indices_string_left in ['ccea']:
    final_indices_string = 'KLCW'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(c_alpha)])

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(c_alpha)]) 

    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(w_beta), sqa_plus.desOp(c_alpha)]) 

#CVEA
elif indices_string_left in ['cvea']:
    final_indices_string = 'KLCW'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(c_alpha)])

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(c_alpha)]) 

    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(w_beta), sqa_plus.desOp(c_alpha)]) 

#CCEE
elif indices_string_left in ['ccee']:
    final_indices_string = 'KLCD'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(d_alpha), sqa_plus.desOp(c_alpha)]) 

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)]) 

    elif spin_indices_string == 'bbbb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_beta), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_beta)]) 

    elif spin_indices_string == 'baba':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(d_alpha), sqa_plus.desOp(c_beta)]) 

#CVEE
elif indices_string_left in ['cvee']:
    final_indices_string = 'KLCD'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_alpha), sqa_plus.desOp(d_alpha), sqa_plus.desOp(c_alpha)]) 

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_val_beta), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)]) 

#CAEE
elif indices_string_left in ['caee']:
    final_indices_string = 'KWCD'

    if spin_indices_string == 'aaaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(d_alpha), sqa_plus.desOp(c_alpha)])

    elif spin_indices_string == 'abab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)])

    elif spin_indices_string == 'baab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)])


## Right Op: h^{(1)\dag}
if indices_string_right in ['caaa']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, x_alpha, y_alpha, z_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, x_beta,  y_alpha, z_beta],  [])]
#    X_tensor_abba = [sqa_plus.tensor('X_abba', [i_alpha, x_beta,  y_beta, z_alpha], [])]
    X_tensor_abba = [sqa_plus.tensor('X_abab', [i_alpha, x_beta, z_alpha, y_beta], [])]

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  x_beta,  y_beta, z_beta],  X_sym_2)]
    X_tensor_baba = [sqa_plus.tensor('X_baba', [i_beta,  x_alpha, y_beta, z_alpha], [])]
#    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  x_alpha, y_alpha, z_beta],  [])]
    X_tensor_baab = [sqa_plus.tensor('X_baba', [i_beta,  x_alpha, z_beta, y_alpha],  [])]

    terms_right = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta),  sqa_plus.desOp(x_beta),  sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(-1.0, [], X_tensor_baab + [sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_beta),  sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)]),

                   sqa_plus.term(0.5, [], X_tensor_bbbb + [sqa_plus.creOp(y_beta), sqa_plus.creOp(z_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_beta)]),
                   sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(y_beta), sqa_plus.creOp(z_alpha),  sqa_plus.desOp(x_alpha),  sqa_plus.desOp(i_beta)]),
                   sqa_plus.term(-1.0, [], X_tensor_abba + [sqa_plus.creOp(y_beta), sqa_plus.creOp(z_alpha),  sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)])]

elif indices_string_right in ['caea']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, x_alpha, a_alpha, y_alpha], [])]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, x_beta,  a_alpha, y_beta],  [])]
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  x_alpha, a_alpha, y_beta],  [])]

    X_tensor_bbbb = [sqa_plus.tensor('X_aaaa', [i_beta,  x_beta,  a_beta, y_beta],  [])]
    X_tensor_baba = [sqa_plus.tensor('X_abab', [i_beta,  x_alpha, a_beta, y_alpha], [])]
    X_tensor_abba = [sqa_plus.tensor('X_baab', [i_alpha, x_beta,  a_beta, y_alpha], [])]

    terms_right = [sqa_plus.term(1.0, [], X_tensor_aaaa + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta),  sqa_plus.desOp(x_beta),  sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta),  sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)]),

                   sqa_plus.term(1.0, [], X_tensor_aaaa + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta),  sqa_plus.desOp(x_beta),  sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_beta),  sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)])]


elif indices_string_right in ['ccaa']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_alpha, x_alpha, y_alpha], X_sym_3)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_beta,  x_alpha, y_beta],  [])]
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  j_alpha, x_alpha, y_beta],  [])]

    terms_right = [sqa_plus.term(0.25, [], X_tensor_aaaa + [sqa_plus.creOp(x_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(x_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(j_beta), sqa_plus.desOp(i_alpha)])]
                   #sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(x_alpha), sqa_plus.creOp(y_beta), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_beta)])]

elif indices_string_right in ['cvaa']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_val_alpha, x_alpha, y_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_val_beta,  x_alpha, y_beta],  [])]
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  j_val_alpha, x_alpha, y_beta],  [])]

    terms_right = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(x_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(j_val_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(x_alpha), sqa_plus.creOp(y_beta),  sqa_plus.desOp(j_val_beta),  sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(x_alpha), sqa_plus.creOp(y_beta),  sqa_plus.desOp(j_val_alpha), sqa_plus.desOp(i_beta)])]

elif indices_string_right in ['ccea']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_alpha, a_alpha, x_alpha], X_sym_1)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_beta,  a_alpha, x_beta],  [])]
    X_tensor_baba = [sqa_plus.tensor('X_baba', [i_beta,  j_alpha, a_beta, x_alpha],  [])]
    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  j_beta,  a_beta, x_beta],   X_sym_1)]

    terms_right = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_beta),  sqa_plus.desOp(j_beta),  sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(a_beta), sqa_plus.creOp(x_alpha),  sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_beta)]),
                   sqa_plus.term(1.0, [], X_tensor_bbbb + [sqa_plus.creOp(a_beta), sqa_plus.creOp(x_beta),  sqa_plus.desOp(j_beta), sqa_plus.desOp(i_beta)])]


elif indices_string_right in ['cvea']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_val_alpha, a_alpha, x_alpha], [])]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_val_beta,  a_alpha, x_beta],  [])]
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  j_val_alpha, a_alpha, x_beta],  [])]

    terms_right = [sqa_plus.term(1.0, [], X_tensor_aaaa + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(j_val_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_beta),  sqa_plus.desOp(j_val_beta),  sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(x_beta),  sqa_plus.desOp(j_val_alpha), sqa_plus.desOp(i_beta)])]

elif indices_string_right in ['ccee']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_alpha, a_alpha, b_alpha], X_sym_3)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_beta,  a_alpha, b_beta], [])] 
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  j_alpha, a_alpha, b_beta], [])] 

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  j_beta,  a_beta, b_beta],  X_sym_3)]
#    X_tensor_baba = [sqa_plus.tensor('X_baba', [i_beta,  j_alpha, a_beta, b_alpha], [])]
    X_tensor_abba = [sqa_plus.tensor('X_abba', [i_alpha, j_beta,  a_beta, b_alpha], [])]


    ##X_baba as transpose
    X_tensor_baba = [sqa_plus.tensor('X_abab', [j_alpha, i_beta, b_alpha, a_beta], [])]
 
    terms_right = [sqa_plus.term(0.25, [], X_tensor_aaaa + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_alpha), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0,  [], X_tensor_abab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_beta), sqa_plus.desOp(j_beta), sqa_plus.desOp(i_alpha)]),
                   #sqa_plus.term(1.0,  [], X_tensor_baab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_beta), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_beta)]),

                   sqa_plus.term(0.25, [], X_tensor_bbbb + [sqa_plus.creOp(a_beta), sqa_plus.creOp(b_beta), sqa_plus.desOp(j_beta), sqa_plus.desOp(i_beta)]),
                   sqa_plus.term(1.0,  [], X_tensor_baba + [sqa_plus.creOp(a_beta), sqa_plus.creOp(b_alpha), sqa_plus.desOp(j_alpha), sqa_plus.desOp(i_beta)])]
                   #sqa_plus.term(1.0,  [], X_tensor_abba + [sqa_plus.creOp(a_beta), sqa_plus.creOp(b_alpha), sqa_plus.desOp(j_beta), sqa_plus.desOp(i_alpha)])]

elif indices_string_right in ['cvee']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_val_alpha, a_alpha, b_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_val_beta,  a_alpha, b_beta],  [])]

    terms_right = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_alpha), sqa_plus.desOp(j_val_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0,  [], X_tensor_abab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_beta),  sqa_plus.desOp(j_val_beta),  sqa_plus.desOp(i_alpha)])]

elif indices_string_right in ['caee']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, x_alpha, a_alpha, b_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, x_beta,  a_alpha, b_beta],  [])]
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  x_alpha, a_alpha, b_beta],  [])]

    terms_right = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(i_alpha)]),
                   sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(a_alpha), sqa_plus.creOp(b_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_beta)])]

indices_string = final_indices_string + '_' + spin_indices_string

print("\n## Right operator terms:")
for _term in terms_right:
    print(_term)

print("\n## Left operator terms:")
print(term_left)

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculating the commutator
print("## Calculating the commutator [H(0), h(1)^\dag] ...")
#can set combine = False
terms_commutator = sqa_plus.commutator(terms_Heff, terms_right)

print("\n## Calculating [h^(1), [H(0), h(1)^\dag]] ...")
terms_EE_M11 = sqa_plus.commutator(term_left, terms_commutator)

# Expected value of Spin-Adapted EE M11
expected_EE_M11 = sqa_plus.matrixBlock(terms_EE_M11)

# Spin-Adaptation of EE M11
from custom_si_to_sa_functions import convert_X_si_to_sa
#sqa_plus.options.add_spin_adaptation_custom_function(convert_X_si_to_sa, sqa_plus.options)
expected_EE_M11 = convert_X_si_to_sa(expected_EE_M11, sqa_plus.options)
expected_EE_M11_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M11)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M11_sa, 'sigma_' + indices_string, final_indices_string, use_cvs_tensors = True, rm_core_int = True)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
