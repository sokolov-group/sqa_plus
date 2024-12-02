import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
#order_Heff = 0
order_Heff = 1
#order_Heff = 2

import time
start = time.time()

indices_string = 'ce_ca'
spin_indices_string = 'aa_bb'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: M00 H{:} {:} ({:})".format(order_Heff, indices_string.upper(), spin_indices_string))

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

# Generating Term 

## External Indices
#Alpha
i_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_cor])
a_alpha = sqa_plus.index('A', [tg_alpha, tg_vir])
j_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_cor])
b_alpha = sqa_plus.index('B', [tg_alpha, tg_vir])
x_alpha = sqa_plus.index('X', [tg_alpha, tg_act])
y_alpha = sqa_plus.index('Y', [tg_alpha, tg_act])

#Beta
i_beta = sqa_plus.index('I', [tg_beta, tg_cvs_cor])
a_beta = sqa_plus.index('A', [tg_beta, tg_vir])
j_beta = sqa_plus.index('J', [tg_beta, tg_cvs_cor])
b_beta = sqa_plus.index('B', [tg_beta, tg_vir])
x_beta = sqa_plus.index('X', [tg_beta, tg_act])
y_beta = sqa_plus.index('Y', [tg_beta, tg_act])


## Define operators
#Alpha
cre_i_alpha = sqa_plus.creOp(i_alpha)
des_j_alpha = sqa_plus.desOp(j_alpha)

cre_x_alpha = sqa_plus.creOp(x_alpha)
des_x_alpha = sqa_plus.desOp(x_alpha)
cre_y_alpha = sqa_plus.creOp(y_alpha)

cre_a_alpha = sqa_plus.creOp(a_alpha)
des_a_alpha = sqa_plus.desOp(a_alpha)
cre_b_alpha = sqa_plus.creOp(b_alpha)

#Beta
cre_i_beta = sqa_plus.creOp(i_beta)
des_j_beta = sqa_plus.desOp(j_beta)

cre_x_beta = sqa_plus.creOp(x_beta)
des_x_beta = sqa_plus.desOp(x_beta)
cre_y_beta = sqa_plus.creOp(y_beta)

cre_a_beta = sqa_plus.creOp(a_beta)
des_a_beta = sqa_plus.desOp(a_beta)
cre_b_beta = sqa_plus.creOp(b_beta)

## Define terms
if indices_string in ['ce_ce']:
    print("\n## Generating Term a_I^\dag a_A ... a_B\dag a_J ...")   
    final_indices_string = 'IAJB'

    if spin_indices_string in ['aa_aa']:
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_a_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_b_alpha, des_j_alpha])
    
    elif spin_indices_string in ['bb_bb']:
        l_op = sqa_plus.term(1.0, [], [cre_i_beta, des_a_beta])
        r_op = sqa_plus.term(1.0, [], [cre_b_beta, des_j_beta])
    
    elif spin_indices_string in ['aa_bb']:
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_a_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_b_beta, des_j_beta])
    
    elif spin_indices_string in ['bb_aa']:
        l_op = sqa_plus.term(1.0, [], [cre_i_beta, des_a_beta])
        r_op = sqa_plus.term(1.0, [], [cre_b_alpha, des_j_alpha])

elif indices_string in ['ca_ca']:
    print("\n## Generating Term a_I^\dag a_X ... a_Y\dag a_J ...")
    final_indices_string = 'IXJY'

    if spin_indices_string in ['aa_aa']:
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_x_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_y_alpha, des_j_alpha])

    elif spin_indices_string in ['bb_bb']: 
        l_op = sqa_plus.term(1.0, [], [cre_i_beta, des_x_beta])
        r_op = sqa_plus.term(1.0, [], [cre_y_beta, des_j_beta])

    elif spin_indices_string in ['aa_bb']: 
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_x_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_y_beta, des_j_beta])

    elif spin_indices_string in ['bb_aa']: 
        l_op = sqa_plus.term(1.0, [], [cre_i_beta, des_x_beta])
        r_op = sqa_plus.term(1.0, [], [cre_y_alpha, des_j_alpha])

elif indices_string in ['ca_ce']:
    print("\n## Generating Term a_I^\dag a_X ... a_A\dag a_J ...")
    final_indices_string = 'IXJA'

    if spin_indices_string in ['aa_aa']:
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_x_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_a_alpha, des_j_alpha])
 
    elif spin_indices_string in ['aa_bb']: 
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_x_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_a_beta, des_j_beta])

elif indices_string in ['ce_ca']:
    print("\n## Generating Term a_I^\dag a_A ... a_X\dag a_J ...\n")
    final_indices_string = 'IAJX'

    if spin_indices_string in ['aa_aa']:
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_a_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_x_alpha, des_j_alpha])
 
    elif spin_indices_string in ['aa_bb']: 
        l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_a_alpha])
        r_op = sqa_plus.term(1.0, [], [cre_x_beta, des_j_beta])
 
print("\n## Left operator terms:")
print(l_op)

print("\n## Right operator terms:")
print(r_op)

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculate the commutators
print("\n## Calculate the commutator ... [H({:}), r_op] ...".format(order_Heff))
terms_commutator = sqa_plus.commutator(terms_Heff, r_op)

print("\n## Calculate the commutator [l_op, [H({:}), r_op]] ...".format(order_Heff))
terms_EE_M0 = sqa_plus.commutator(l_op, terms_commutator)

# Expected value of Spin-Integrated EE M0 h0-h0
expected_EE_M0 = sqa_plus.matrixBlock(terms_EE_M0)

# Spin-Adaptation of EE M0 h0-h0
expected_EE_M0_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M0)

# Generate Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M0_sa, 'temp', final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
