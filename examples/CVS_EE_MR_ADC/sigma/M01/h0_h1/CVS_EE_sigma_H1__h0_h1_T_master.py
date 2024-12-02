import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 1

import time
start = time.time()

indices_string = 'ccea_ca'
spin_indices_string = 'bb'

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
k_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_cor])

w_alpha = sqa_plus.index('W', [tg_alpha, tg_act])
w_beta  = sqa_plus.index('W', [tg_beta,  tg_act])

c_alpha = sqa_plus.index('C', [tg_alpha, tg_vir])
c_beta  = sqa_plus.index('C', [tg_beta,  tg_vir])

## Dummy Indices
dummy = True

i_alpha = sqa_plus.index('i', [tg_alpha, tg_cvs_cor], dummy)
i_beta  = sqa_plus.index('i', [tg_beta,  tg_cvs_cor], dummy)

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

a_alpha = sqa_plus.index('a', [tg_alpha, tg_vir], dummy)
a_beta  = sqa_plus.index('a', [tg_beta,  tg_vir], dummy)

b_alpha = sqa_plus.index('b', [tg_alpha, tg_vir], dummy)
b_beta  = sqa_plus.index('b', [tg_beta,  tg_vir], dummy)

## Define terms
# a^{PQ}_{RS} = a_P^\dag a_Q^\dag a_S a_R

X_sym_1 = [sqa_plus.symmetry((1,0,2,3), -1)] #symmetry: pprs
X_sym_2 = [sqa_plus.symmetry((0,1,3,2), -1)] #symmetry: pqrr
X_sym_3 = X_sym_1 + X_sym_2 #symmetry: pprr

indices_string_left = indices_string.split('_')[0] 
indices_string_right  = indices_string.split('_')[1]

##Dummy Indices
from sqa_plus.sqaIndexList import indexLists
cvs_alpha_inds = indexLists.cvs_core_alpha
val_alpha_inds = indexLists.cvs_valence_alpha
act_alpha_inds = indexLists.active_alpha
vir_alpha_inds = indexLists.virtual_alpha

cvs_beta_inds = indexLists.cvs_core_beta
val_beta_inds = indexLists.cvs_valence_beta
act_beta_inds = indexLists.active_beta
vir_beta_inds = indexLists.virtual_beta


if indices_string_left in ['caaa']:

    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, x_alpha, y_alpha, z_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, x_beta,  y_alpha, z_beta], [])]
    X_tensor_abba = [sqa_plus.tensor('X_abab', [i_alpha, x_beta,  y_beta, z_alpha], [])]

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  x_beta,  y_beta, z_beta],  X_sym_2)]
    X_tensor_baba = [sqa_plus.tensor('X_baba', [i_beta,  x_alpha, y_beta, z_alpha], [])]
    X_tensor_baab = [sqa_plus.tensor('X_baba', [i_beta,  x_alpha, y_alpha, z_beta], [])]

    terms_left = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(y_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(y_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(z_alpha), sqa_plus.desOp(y_beta)]),
                  sqa_plus.term(0.5, [], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(y_beta)]),
                  sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(y_beta)]),
                  sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(z_beta), sqa_plus.desOp(y_alpha)])]

elif indices_string_left in ['caea']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, x_alpha, a_alpha, y_alpha])]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, x_beta,  a_alpha, y_beta])]
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  x_alpha, a_alpha, y_beta])]

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  x_beta,  a_beta, y_beta])]
    X_tensor_baba = [sqa_plus.tensor('X_baba', [i_beta,  x_alpha, a_beta, y_alpha])]
    X_tensor_abba = [sqa_plus.tensor('X_abba', [i_alpha, x_beta,  a_beta, y_alpha])]
    
    terms_left = [sqa_plus.term(1.0, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_beta), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(a_beta)]),
                  sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(a_beta)]),
                  sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(y_alpha), sqa_plus.desOp(a_beta)])]

elif indices_string_left in ['ccaa']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_alpha, x_alpha, y_alpha], X_sym_3)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_beta,  x_alpha, y_beta], [])]
    X_tensor_baab = [sqa_plus.tensor('X_abab', [i_beta,  j_alpha, x_alpha, y_beta], [])]

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  j_beta,  x_beta, y_beta],  X_sym_3)]
    X_tensor_baba = [sqa_plus.tensor('X_abab', [i_beta,  j_alpha, x_beta, y_alpha], [])]
    X_tensor_abba = [sqa_plus.tensor('X_abab', [i_alpha, j_beta,  x_beta, y_alpha], [])]

    terms_left = [sqa_plus.term(0.25,[], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]),
                  #sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]),

                  sqa_plus.term(0.25,[], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_beta)])]
                  #sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_beta)]),
                  #sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_beta)])]
 
elif indices_string_left in ['cvaa']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_val_alpha, x_alpha, y_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_val_beta,  x_alpha, y_beta],  [])]
    #X_tensor_baab = [sqa_plus.tensor('X_abab', [i_beta,  j_val_alpha, y_beta,  x_alpha], X_sym_2)]

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  j_val_beta,  x_beta, y_beta],  X_sym_2)]
    #X_tensor_baba = [sqa_plus.tensor('X_abab', [i_beta,  j_val_alpha, x_beta, y_alpha], [])]
    #X_tensor_abba = [sqa_plus.tensor('X_baab', [i_alpha, j_val_beta,  x_beta, y_alpha], [])]

    terms_left = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]),
                  #sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_alpha)]),
     
                  sqa_plus.term(0.5, [], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(y_beta), sqa_plus.desOp(x_beta)])]
                  #sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_beta)])]
                  #sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_beta)])]

elif indices_string_left in ['ccea']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_alpha, a_alpha, x_alpha], X_sym_1)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_beta,  a_alpha, x_beta],  [])]
    #X_tensor_baab = [sqa_plus.tensor('X_abab', [i_beta,  j_alpha, a_alpha, x_beta],  [])]

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  j_beta,  a_beta, x_beta],  X_sym_1)]
    #X_tensor_baba = [sqa_plus.tensor('X_abab', [i_beta,  j_alpha, a_beta, x_alpha], [])]
    X_tensor_abba = [sqa_plus.tensor('X_abba', [i_alpha, j_beta,  a_beta, x_alpha], [])]

    terms_left = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(a_alpha)]),
                  #sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(x_beta), sqa_plus.desOp(a_alpha)]),
     
                  sqa_plus.term(0.5, [], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(a_beta)]),
                  #sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(a_beta)])]
                  sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(a_beta)])]

elif indices_string_left in ['cvea']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_val_alpha, a_alpha, x_alpha], X_sym_3)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_val_beta,  a_alpha, x_beta],  X_sym_3)]
    X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  j_val_alpha, a_alpha, x_beta],  X_sym_3)]

    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  j_val_beta,  a_beta, x_beta],  X_sym_3)]
    X_tensor_baba = [sqa_plus.tensor('X_abab', [i_beta,  j_val_alpha, a_beta, x_alpha], X_sym_3)]
    X_tensor_abba = [sqa_plus.tensor('X_baab', [i_alpha, j_val_beta,  a_beta, x_alpha], X_sym_3)]

    terms_left = [sqa_plus.term(1.0, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(x_beta), sqa_plus.desOp(a_alpha)]),
     
                  sqa_plus.term(1.0, [], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(x_beta), sqa_plus.desOp(a_beta)]),
                  sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(a_beta)]),
                  sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(x_alpha), sqa_plus.desOp(a_beta)])]

elif indices_string_left in ['ccee']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_alpha, a_alpha, b_alpha], X_sym_3)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_beta,  a_alpha, b_beta], [])] 
    #X_tensor_baab = [sqa_plus.tensor('X_baab', [i_beta,  j_alpha,  a_alpha, b_beta], [])] 
    X_tensor_bbbb = [sqa_plus.tensor('X_bbbb', [i_beta,  j_beta,  a_beta, b_beta],  X_sym_3)]
    #X_tensor_baba = [sqa_plus.tensor('X_baba', [i_beta,  j_alpha, a_beta, b_alpha], [])]
    #X_tensor_abba = [sqa_plus.tensor('X_abba', [i_alpha, j_beta,  a_beta, b_alpha], [])]

    terms_left = [sqa_plus.term(0.25, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(0.25,[], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_beta)])]
                  
                  #sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)]),
                  #sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_beta)]),
                  #sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_beta)])]

elif indices_string_left in ['cvee']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, j_val_alpha, a_alpha, b_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, j_val_beta,  a_alpha, b_beta], [])]
    X_tensor_bbbb = [sqa_plus.tensor('X_aaaa', [i_beta,  j_val_beta,  a_beta, b_beta],  X_sym_2)]
    #X_tensor_baba = [sqa_plus.tensor('X_abab', [i_beta,  j_val_alpha, a_beta, b_alpha], X_sym_3)]
 
    terms_left = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(0.5, [], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_beta)])]
                  #sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_beta)])]
  
elif indices_string_left in ['caee']:
    X_tensor_aaaa = [sqa_plus.tensor('X_aaaa', [i_alpha, x_alpha, a_alpha, b_alpha], X_sym_2)]
    X_tensor_abab = [sqa_plus.tensor('X_abab', [i_alpha, x_beta,  a_alpha, b_beta],  [])]
    #X_tensor_baab = [sqa_plus.tensor('X_abab', [i_beta,  x_alpha, a_alpha, b_beta],  [])]

    X_tensor_bbbb = [sqa_plus.tensor('X_aaaa', [i_beta,  x_beta,  a_beta, b_beta],  X_sym_2)]
    X_tensor_baba = [sqa_plus.tensor('X_abab', [i_beta,  x_alpha, a_beta, b_alpha], [])]
    #X_tensor_abba = [sqa_plus.tensor('X_baab', [i_alpha, x_beta,  a_beta, b_alpha], [])]

    terms_left = [sqa_plus.term(0.5, [], X_tensor_aaaa + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_alpha)]),
                  sqa_plus.term(1.0, [], X_tensor_abab + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)]),
                  #sqa_plus.term(1.0, [], X_tensor_baab + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_alpha)])]
                  sqa_plus.term(0.5, [], X_tensor_bbbb + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_beta), sqa_plus.desOp(b_beta), sqa_plus.desOp(a_beta)]),
                  sqa_plus.term(1.0, [], X_tensor_baba + [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_beta)])]
                  #sqa_plus.term(1.0, [], X_tensor_abba + [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(b_alpha), sqa_plus.desOp(a_beta)])]
                
if indices_string_right in ['ce']:
    final_indices_string = 'KC'
    if spin_indices_string == 'aa':
        term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(c_alpha), sqa_plus.desOp(k_alpha)])
    if spin_indices_string == 'bb':
        term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(c_beta), sqa_plus.desOp(k_beta)])

elif indices_string_right in ['ca']:
    final_indices_string = 'KW'
    if spin_indices_string == 'aa':
        term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_alpha), sqa_plus.desOp(k_alpha)])
    if spin_indices_string == 'bb':
        term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(w_beta), sqa_plus.desOp(k_beta)])

indices_string = final_indices_string + '_' + spin_indices_string

print("\n## Left operator terms:")
for _term in terms_left:
    print(_term)

print("\n## Right operator terms:")
print(term_right)

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculating the commutator
print("## Calculating the commutator [H({:}), h(0)^\dag] ...".format(order_Heff))
terms_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating [h(1), [H({:}), h(0)^\dag]] ...".format(order_Heff))
terms_EE_M01 = sqa_plus.commutator(terms_left, terms_commutator)

# Expected value of Spin-Adapted EE M01
expected_EE_M01 = sqa_plus.matrixBlock(terms_EE_M01)

# Spin-Adaptation of EE M01
from custom_si_to_sa_functions import convert_X_si_to_sa
sqa_plus.options.add_spin_adaptation_custom_function(convert_X_si_to_sa)
expected_EE_M01_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M01)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M01_sa, 'sigma_' + indices_string, final_indices_string, use_cvs_tensors = True, rm_core_int = True)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
