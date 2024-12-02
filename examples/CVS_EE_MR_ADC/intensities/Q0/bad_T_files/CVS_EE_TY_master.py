import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
from sqa_plus.sqaIndexList import indexLists, dummyIndexList
q_order = 0

import time
start = time.time()

indices_string = 'q_ce'
spin_indices_string = 'aa'

sqa_plus.options.print_header("Spin-Adapted CVS-EE: TY Q{:} {:} ({:})".format(q_order, indices_string.upper(), spin_indices_string))

# Generating operators
print("\n## Generating operators ...\n")

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type
tg_g = tg_cvs_cor + tg_cvs_val + tg_act + tg_vir

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices
#k_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_cor])
#k_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_cor])

#w_alpha = sqa_plus.index('W', [tg_alpha, tg_act])
#w_beta  = sqa_plus.index('W', [tg_beta,  tg_act])

#c_alpha = sqa_plus.index('C', [tg_alpha, tg_vir])
#c_beta  = sqa_plus.index('C', [tg_beta,  tg_vir])

##Dummy Indices
cvs_alpha_inds = indexLists.cvs_core_alpha
val_alpha_inds = indexLists.cvs_valence_alpha
act_alpha_inds = indexLists.active_alpha
vir_alpha_inds = indexLists.virtual_alpha

cvs_beta_inds = indexLists.cvs_core_beta
val_beta_inds = indexLists.cvs_valence_beta
act_beta_inds = indexLists.active_beta
vir_beta_inds = indexLists.virtual_beta

## Define terms
term_q0 = []
#dSym = [sqa_plus.symmetry((1,0),1)]
#dSym = [sqa_plus.symmetry((0,2,1),1)]

final_indices_list = []
TY_index = []

# CC
#p_alpha = cvs_alpha_inds.new_index()
#q_alpha = cvs_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_cc', [r, p_alpha, q_alpha], dSym)

p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_cor])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('IJ')
TY_index.append('[s_c:f_c, s_c:f_c]')

# VV
#p_alpha = val_alpha_inds.new_index()
#q_alpha = val_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_vv', [r, p_alpha, q_alpha], dSym)

p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_val])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('IJ')
TY_index.append('[s_v:f_v, s_v:f_v]')

# AA
#p_alpha = act_alpha_inds.new_index()
#q_alpha = act_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_aa', [r, p_alpha, q_alpha], dSym)

p_alpha = sqa_plus.index('P', [tg_alpha, tg_act])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('XY')
TY_index.append('[s_a:f_a, s_a:f_a]')

# EE
#p_alpha = vir_alpha_inds.new_index()
#q_alpha = vir_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_ee', [r, p_alpha, q_alpha], dSym)

p_alpha = sqa_plus.index('P', [tg_alpha, tg_vir])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('AB')
TY_index.append('[s_e:f_e, s_e:f_e]')

# CV
#p_alpha = cvs_alpha_inds.new_index()
#q_alpha = val_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_cv', [r, p_alpha, q_alpha])

p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_val])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('IJ')
TY_index.append('[s_c:f_c, s_v:f_v]')

#dTen = sqa_plus.tensor('d_cv', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))
final_indices_list.append('JI')
TY_index.append('[s_v:f_v, s_c:f_c]')

# CA
#p_alpha = cvs_alpha_inds.new_index()
#q_alpha = act_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_ca', [r, p_alpha, q_alpha])

p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('IX')
TY_index.append('[s_c:f_c, s_a:f_a]')

#dTen = sqa_plus.tensor('d_ca', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))
final_indices_list.append('XI')
TY_index.append('[s_a:f_a, s_c:f_c]')

# VA
#p_alpha = val_alpha_inds.new_index()
#q_alpha = act_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_va', [r, p_alpha, q_alpha])

p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('IX')
TY_index.append('[s_v:f_v, s_a:f_a]')

#dTen = sqa_plus.tensor('d_va', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))
final_indices_list.append('XI')
TY_index.append('[s_a:f_a, s_v:f_v]')

# CE
#p_alpha = cvs_alpha_inds.new_index()
#q_alpha = vir_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_ce', [r, p_alpha, q_alpha])

p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('IA')
TY_index.append('[s_c:f_c, s_e:f_e]')

#dTen = sqa_plus.tensor('d_ce', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))
final_indices_list.append('AI')
TY_index.append('[s_e:f_e, s_c:f_c]')

# VE
#p_alpha = val_alpha_inds.new_index()
#q_alpha = vir_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_ve', [r, p_alpha, q_alpha])

p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('IA')
TY_index.append('[s_v:f_v, s_e:f_e]')

#dTen = sqa_plus.tensor('d_ve', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))
final_indices_list.append('AI')
TY_index.append('[s_e:f_e, s_v:f_v]')

# AE
#p_alpha = act_alpha_inds.new_index()
#q_alpha = vir_alpha_inds.new_index()
#dTen = sqa_plus.tensor('d_ae', [r, p_alpha, q_alpha])

p_alpha = sqa_plus.index('P', [tg_alpha, tg_act])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('XA')
TY_index.append('[s_a:f_a, s_e:f_e]')

#dTen = sqa_plus.tensor('d_ae', [r, p_alpha, q_alpha])
term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)]))
final_indices_list.append('AX')
TY_index.append('[s_e:f_e, s_a:f_a]')

if indices_string in ['q_ca']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    w_alpha = act_alpha_inds.new_index()
    w_beta  = act_beta_inds.new_index()

    Y_aa = sqa_plus.tensor('Y_KW__aa', [k_alpha, w_alpha])
    Y_bb = sqa_plus.tensor('Y_KW__bb', [k_beta, w_beta])

    term_h  = [sqa_plus.term(1.0, [], [Y_aa, sqa_plus.creOp(w_alpha), sqa_plus.desOp(k_alpha)]),
               sqa_plus.term(1.0, [], [Y_bb, sqa_plus.creOp(w_beta), sqa_plus.desOp(k_beta)])]
    
    #final_indices_string = 'KW'

elif indices_string in ['q_ce']:

    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    c_alpha = vir_alpha_inds.new_index()
    c_beta  = vir_beta_inds.new_index()

    Y_aa = sqa_plus.tensor('Y_KC__aa', [k_alpha, c_alpha]) 
    Y_bb = sqa_plus.tensor('Y_KC__bb', [k_beta, c_beta]) 

    term_h  = [sqa_plus.term(1.0, [], [Y_aa, sqa_plus.creOp(c_alpha), sqa_plus.desOp(k_alpha)]),
               sqa_plus.term(1.0, [], [Y_bb, sqa_plus.creOp(c_beta), sqa_plus.desOp(k_beta)])]
    
    #final_indices_string = 'KC'

## Evaluate transition moments contributions

TY_set = zip(term_q0, TY_index)
for term_q0, index in TY_set:

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
    
        terms_q = sqa_plus.commutator(term_q0, terms_T, combine = True)
    
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
    result = sqa_plus.genEinsum(expected_EE_TY_sa, 'TY' + index, 'PQ')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
