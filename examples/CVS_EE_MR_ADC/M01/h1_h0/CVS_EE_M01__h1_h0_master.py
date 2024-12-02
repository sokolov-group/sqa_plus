import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 1

import time
start = time.time()

indices_string = 'ccaa_ca'
spin_indices_string = 'bbbb_aa'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: M01 H{:} {:} ({:})".format(order_Heff, indices_string.upper(), spin_indices_string))

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
indices_string_left  = indices_string.split('_')[0]
indices_string_right = indices_string.split('_')[1]

#LHS
if indices_string_left in ['ccaa']:
    l_ind = 'KLWZ'
    ##sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(y_alpha), sqa_plus.desOp(x_alpha)]) #LHS = KLWZ
 
    if spin_indices_string == 'aaaa_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(w_alpha)])

    if spin_indices_string == 'abba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(w_alpha)]) 

    if spin_indices_string == 'baba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(z_beta), sqa_plus.desOp(w_alpha)]) 

    ### blocks to check bbbb, baab, abab 
    if spin_indices_string == 'bbbb_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(w_beta)])

    if spin_indices_string == 'baab_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(w_beta)])

    if spin_indices_string == 'abab_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(z_alpha), sqa_plus.desOp(w_beta)])

if indices_string_left in ['ccee']:
    l_ind = 'KLCD'
     
    if spin_indices_string == 'aaaa_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(d_alpha), sqa_plus.desOp(c_alpha)])

    if spin_indices_string == 'abba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)]) 

    if spin_indices_string == 'baba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)]) 

if indices_string_left in ['ccea']:
    l_ind = 'KLBW'
     
    if spin_indices_string == 'aaaa_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(b_alpha)])

    if spin_indices_string == 'abba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(l_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(b_alpha)])

    if spin_indices_string == 'baba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(l_alpha), sqa_plus.desOp(w_beta), sqa_plus.desOp(b_alpha)])

if indices_string_left in ['caee']:
    l_ind = 'KWCD'
 
    if spin_indices_string == 'aaaa_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(d_alpha), sqa_plus.desOp(c_alpha)])

    if spin_indices_string == 'abba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)])

    if spin_indices_string == 'baba_aa': 
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(d_beta), sqa_plus.desOp(c_alpha)])
 
if indices_string_left in ['caea']:
    l_ind = 'KWBZ'

    if spin_indices_string == 'aaaa_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(b_alpha)])

    elif spin_indices_string == 'abba_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(b_alpha)]) 

    elif spin_indices_string == 'baba_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_beta), sqa_plus.desOp(b_alpha)]) 

if indices_string_left in ['caaa']:
    l_ind = 'KWYZ'

    if spin_indices_string == 'aaaa_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'abba_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_alpha), sqa_plus.creOp(w_beta), sqa_plus.desOp(z_beta), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'baba_aa':
        term_left = sqa_plus.term(1.0, [], [sqa_plus.creOp(k_beta), sqa_plus.creOp(w_alpha), sqa_plus.desOp(z_beta), sqa_plus.desOp(y_alpha)])

#RHS
if indices_string_right in ['ce']:
    r_ind = 'IA'
    term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(a_alpha), sqa_plus.desOp(i_alpha)])
#    term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(a_beta), sqa_plus.desOp(i_beta)])


elif indices_string_right in ['ca']:
    r_ind = 'IX'
    term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(x_alpha), sqa_plus.desOp(i_alpha)])
#    term_right   =  sqa_plus.term(1.0, [], [sqa_plus.creOp(x_beta), sqa_plus.desOp(i_beta)])

#spin_indices_string = spin_indices_string.split('_')[0] + '_' + 'bb'

indices_string = indices_string + '_' + spin_indices_string

print("\n## Left operator terms:")
print(term_left)

print("\n## Right operator terms:")
print(term_right)

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculating the commutator
print("## Calculating the commutator [H({:}), h^(0)\dag] ...".format(order_Heff))
first_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating [h^(1), [H({:}), h^(0)\dag]] ...".format(order_Heff))
terms_EE_M01 = sqa_plus.commutator(term_left, first_commutator)

# Expected value of Spin-Adapted EE M11
expected_EE_M01 = sqa_plus.matrixBlock(terms_EE_M01)

## Spin-Adaptation of EE M11
expected_EE_M01_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M01)

# Generating Numpy einsum equations
final_indices_string = r_ind + l_ind
result = sqa_plus.genEinsum(expected_EE_M01_sa, 'M_' + indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
