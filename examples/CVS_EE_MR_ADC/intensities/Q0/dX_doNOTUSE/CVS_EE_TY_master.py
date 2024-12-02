import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
from sqa_plus.sqaIndexList import indexLists
q_order = 0

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
tg_g = tg_cvs_cor + tg_cvs_val + tg_act + tg_vir
r  = sqa_plus.index('R', [tg_g])

##Dummy Indices
#cor_alpha_inds = indexLists.core_alpha
cvs_alpha_inds = indexLists.cvs_core_alpha
val_alpha_inds = indexLists.cvs_valence_alpha
act_alpha_inds = indexLists.active_alpha
vir_alpha_inds = indexLists.virtual_alpha

#cor_beta_inds = indexLists.core_beta
cvs_beta_inds = indexLists.cvs_core_beta
val_beta_inds = indexLists.cvs_valence_beta
act_beta_inds = indexLists.active_beta
vir_beta_inds = indexLists.virtual_beta

## Define terms
term_q0 = []
dSym = [sqa_plus.symmetry((1,0),1)]

# CC
p_alpha = cvs_alpha_inds.new_index()
q_alpha = cvs_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_cc', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))

# VV
p_alpha = val_alpha_inds.new_index()
q_alpha = val_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_vv', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))

# AA
p_alpha = act_alpha_inds.new_index()
q_alpha = act_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_aa', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))

# EE
p_alpha = vir_alpha_inds.new_index()
q_alpha = vir_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_ee', [p_alpha, q_alpha], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))

# CA
p_alpha = cvs_alpha_inds.new_index()
q_alpha = act_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_ca', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ca', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# VA
p_alpha = val_alpha_inds.new_index()
q_alpha = act_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_va', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_va', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# CE
p_alpha = cvs_alpha_inds.new_index()
q_alpha = vir_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_ce', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ce', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# VE
p_alpha = val_alpha_inds.new_index()
q_alpha = vir_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_ve', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ve', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

# AE
p_alpha = act_alpha_inds.new_index()
q_alpha = vir_alpha_inds.new_index()

dTen = sqa_plus.tensor('d_ae', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
dTen = sqa_plus.tensor('d_ae', [p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))

if indices_string in ['ca']:

    i_alpha = cvs_alpha_inds.new_index()
    x_alpha = act_alpha_inds.new_index()
 
    Y = sqa_plus.tensor('Y_KW__aa', [i_alpha, x_alpha])
    term_h  = sqa_plus.term(1.0, [], [Y, sqa_plus.creOp(x_alpha), sqa_plus.desOp(i_alpha)])
    
    #final_indices_string = 'IX'

elif indices_string in ['caaa']:

    i_alpha = cvs_alpha_inds.new_index()
    x_alpha = act_alpha_inds.new_index()
    y_alpha = act_alpha_inds.new_index()
    z_alpha = act_alpha_inds.new_index()
 
    Y = sqa_plus.tensor('Y_KWYZ__aa', [i_alpha, x_alpha, y_alpha, z_alpha])
    term_h = sqa_plus.term(1.0, [], [Y, sqa_plus.creOp(y_alpha), sqa_plus.creOp(z_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)])

elif indices_string in ['ce']:
    i_alpha = cvs_alpha_inds.new_index()
    a_alpha = vir_alpha_inds.new_index()

    Y = sqa_plus.tensor('Y_KC__aa', [i_alpha, a_alpha]) 
    term_h  = sqa_plus.term(1.0, [], [Y, sqa_plus.creOp(a_alpha), sqa_plus.desOp(i_alpha)])
 
elif indices_string in ['ce_caea']:

    i_alpha = cvs_alpha_inds.new_index()
    x_alpha = act_alpha_inds.new_index()
    a_alpha = vir_alpha_inds.new_index()
    y_alpha = act_alpha_inds.new_index()
 
    Y_1 = sqa_plus.tensor('Y_KC', [i_alpha, a_alpha])
    Y_2 = sqa_plus.tensor('Y_KWCU', [i_alpha, x_alpha, a_alpha, y_alpha])

    op_Y1 = sqa_plus.term(1.0, [], [Y_1, sqa_plus.creOp(a_alpha), sqa_plus.desOp(i_alpha)]) 
    op_Y2 = sqa_plus.term(1.0, [], [Y_2, sqa_plus.creOp(a_alpha), sqa_plus.creOp(y_alpha), sqa_plus.desOp(x_alpha), sqa_plus.desOp(i_alpha)])

    term_h = [op_Y1, op_Y2]
 
    #final_indices_string = 'IA'

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

terms_EE_TY = sqa_plus.commutator(terms_q, term_h, combine = True)

# Expected value of Spin-Adapted EE T
expected_EE_TY = sqa_plus.matrixBlock(terms_EE_TY)

# Spin-Adaptation of EE T
expected_EE_TY_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_TY)

# Generating Numpy einsum equations
#result = sqa_plus.genEinsum(expected_EE_TY_sa, 'TY_' + indices_string, final_indices_string)
result = sqa_plus.genEinsum(expected_EE_TY_sa, 'TY_' + indices_string, 'R')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
