import sqa_plus
sqa_plus.options.cvs_approach = True
sqa_plus.options.spin_integrated = True

from sqa_plus.sqaIndexList import indexLists
q_order = 1

import time
start = time.time()

indices_string = 'q_caea'
sqa_plus.options.print_header("Spin-Adapted CVS-EE: TY Q{:} {:}".format(q_order, indices_string.upper()))

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

r = sqa_plus.index('R', [tg_g])

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
final_indices_list = []
TY_index = []

# CC
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_cor])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('CC (aa)')
TY_index.append('[:, s_c:f_c, s_c:f_c]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_cvs_cor])
q_beta = sqa_plus.index('Q', [tg_beta, tg_cvs_cor])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]))
final_indices_list.append('CC (bb)')
TY_index.append('[:, s_c:f_c, s_c:f_c]')

# VV
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_val])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('VV (aa)')
TY_index.append('[:, s_v:f_v, s_v:f_v]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_cvs_val])
q_beta = sqa_plus.index('Q', [tg_beta, tg_cvs_val])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]))
final_indices_list.append('VV (bb)')
TY_index.append('[:, s_v:f_v, s_v:f_v]')

# AA
p_alpha = sqa_plus.index('P', [tg_alpha, tg_act])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('AA (aa)')
TY_index.append('[:, s_a:f_a, s_a:f_a]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_act])
q_beta = sqa_plus.index('Q', [tg_beta, tg_act])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]))
final_indices_list.append('AA (bb)')
TY_index.append('[:, s_a:f_a, s_a:f_a]')

# EE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_vir])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]))
final_indices_list.append('EE (aa)')
TY_index.append('[:, s_e:f_e, s_e:f_e]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_vir])
q_beta = sqa_plus.index('Q', [tg_beta, tg_vir])

term_q0.append(sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]))
final_indices_list.append('EE (bb)')
TY_index.append('[:, s_e:f_e, s_e:f_e]')

# CV
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_cvs_val])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)])]

term_q0.append(temp)
final_indices_list.append('CV (aa)')
TY_index.append('[:, s_c:f_c, s_v:f_v]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_cvs_cor])
q_beta = sqa_plus.index('Q', [tg_beta, tg_cvs_val])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_beta), sqa_plus.desOp(p_beta)])]

term_q0.append(temp)
final_indices_list.append('CV (bb)')
TY_index.append('[:, s_c:f_c, s_v:f_v]')

# CA
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)])]

term_q0.append(temp)
final_indices_list.append('CA (aa)')
TY_index.append('[:, s_c:f_c, s_a:f_a]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_cvs_cor])
q_beta = sqa_plus.index('Q', [tg_beta, tg_act])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_beta), sqa_plus.desOp(p_beta)])]

term_q0.append(temp)
final_indices_list.append('CA (bb)')
TY_index.append('[:, s_c:f_c, s_a:f_a]')

# VA
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_act])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)])]

term_q0.append(temp)
final_indices_list.append('VA (aa)')
TY_index.append('[:, s_v:f_v, s_a:f_a]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_cvs_val])
q_beta = sqa_plus.index('Q', [tg_beta, tg_act])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_beta), sqa_plus.desOp(p_beta)])]

term_q0.append(temp)
final_indices_list.append('VA (bb)')
TY_index.append('[:, s_v:f_v, s_a:f_a]')

# CE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_cor])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)])]

term_q0.append(temp)
final_indices_list.append('CE (aa)')
TY_index.append('[:, s_c:f_c, s_e:f_e]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_cvs_cor])
q_beta = sqa_plus.index('Q', [tg_beta, tg_vir])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_beta), sqa_plus.desOp(p_beta)])]

term_q0.append(temp)
final_indices_list.append('CE (bb)')
TY_index.append('[:, s_c:f_c, s_e:f_e]')

# VE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_cvs_val])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)])]

term_q0.append(temp)
final_indices_list.append('VE (aa)')
TY_index.append('[:, s_v:f_v, s_e:f_e]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_cvs_val])
q_beta = sqa_plus.index('Q', [tg_beta, tg_vir])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_beta), sqa_plus.desOp(p_beta)])]

term_q0.append(temp)
final_indices_list.append('VE (bb)')
TY_index.append('[:, s_v:f_v, s_e:f_e]')

# AE
p_alpha = sqa_plus.index('P', [tg_alpha, tg_act])
q_alpha = sqa_plus.index('Q', [tg_alpha, tg_vir])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_alpha), sqa_plus.desOp(q_alpha)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_alpha), sqa_plus.desOp(p_alpha)])]

term_q0.append(temp)
final_indices_list.append('AE (aa)')
TY_index.append('[:, s_a:f_a, s_e:f_e]')

## beta
p_beta = sqa_plus.index('P', [tg_beta, tg_act])
q_beta = sqa_plus.index('Q', [tg_beta, tg_vir])

temp = [sqa_plus.term(1.0, [], [sqa_plus.creOp(p_beta), sqa_plus.desOp(q_beta)]), 
        sqa_plus.term(1.0, [], [sqa_plus.creOp(q_beta), sqa_plus.desOp(p_beta)])]

term_q0.append(temp)
final_indices_list.append('AE (bb)')
TY_index.append('[:, s_a:f_a, s_e:f_e]')

if indices_string in ['q_ca']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    w_alpha = act_alpha_inds.new_index()
    w_beta  = act_beta_inds.new_index()

    Y_aa = sqa_plus.tensor('Y_KW__aa', [r, k_alpha, w_alpha])
    Y_bb = sqa_plus.tensor('Y_KW__bb', [r, k_beta, w_beta])

    term_h  = [sqa_plus.term(1.0, [], [Y_aa, sqa_plus.creOp(w_alpha), sqa_plus.desOp(k_alpha)]),
               sqa_plus.term(1.0, [], [Y_bb, sqa_plus.creOp(w_beta), sqa_plus.desOp(k_beta)])]

elif indices_string in ['q_ce']:

    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    c_alpha = vir_alpha_inds.new_index()
    c_beta  = vir_beta_inds.new_index()

    Y_aa = sqa_plus.tensor('Y_KC__aa', [r, k_alpha, c_alpha]) 
    Y_bb = sqa_plus.tensor('Y_KC__bb', [r, k_beta, c_beta]) 

    term_h  = [sqa_plus.term(1.0, [], [Y_aa, sqa_plus.creOp(c_alpha), sqa_plus.desOp(k_alpha)]),
               sqa_plus.term(1.0, [], [Y_bb, sqa_plus.creOp(c_beta), sqa_plus.desOp(k_beta)])]

elif indices_string in ['q_caea']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    w_alpha = act_alpha_inds.new_index()
    w_beta  = act_beta_inds.new_index()

    c_alpha = vir_alpha_inds.new_index()
    c_beta  = vir_beta_inds.new_index()

    u_alpha = act_alpha_inds.new_index()
    u_beta  = act_beta_inds.new_index()

    Y_aaaa = sqa_plus.tensor('Y_KWCU__aaaa', [r, k_alpha, w_alpha, c_alpha, u_alpha])
    Y_abab = sqa_plus.tensor('Y_KWCU__abab', [r, k_alpha, w_beta,  c_alpha, u_beta])
    Y_baab = sqa_plus.tensor('Y_KWCU__baab', [r, k_beta,  w_beta,  c_alpha, u_beta])

    #Y_bbbb = sqa_plus.tensor('Y_KWCU__bbbb', [r, k_beta,  w_beta,  c_beta, u_beta])
    #Y_baba = sqa_plus.tensor('Y_KWCU__baba', [r, k_beta,  w_alpha, c_beta, u_alpha])
    #Y_abba = sqa_plus.tensor('Y_KWCU__abba', [r, k_alpha, w_beta,  c_beta, u_alpha])

    term_h = [sqa_plus.term(1.0, [], [Y_aaaa, sqa_plus.creOp(c_alpha), sqa_plus.creOp(u_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_abab, sqa_plus.creOp(c_alpha), sqa_plus.creOp(u_beta),  sqa_plus.desOp(w_beta),  sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_baab, sqa_plus.creOp(c_alpha), sqa_plus.creOp(u_beta),  sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_beta)])]
              #sqa_plus.term(1.0, [], [Y_bbbb, sqa_plus.creOp(c_beta), sqa_plus.creOp(u_beta), sqa_plus.desOp(w_beta), sqa_plus.desOp(k_beta)]),
              #sqa_plus.term(1.0, [], [Y_baba, sqa_plus.creOp(c_beta), sqa_plus.creOp(u_alpha),  sqa_plus.desOp(w_alpha),  sqa_plus.desOp(k_beta)]),
              #sqa_plus.term(1.0, [], [Y_abba, sqa_plus.creOp(c_beta), sqa_plus.creOp(u_alpha),  sqa_plus.desOp(w_beta), sqa_plus.desOp(k_alpha)])]

elif indices_string in ['q_caaa']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    w_alpha = act_alpha_inds.new_index()
    w_beta  = act_beta_inds.new_index()

    u_alpha = act_alpha_inds.new_index()
    u_beta  = act_beta_inds.new_index()

    v_alpha = act_alpha_inds.new_index()
    v_beta  = act_beta_inds.new_index()

    Y_aaaa = sqa_plus.tensor('Y_KWUV__aaaa', [r, k_alpha, w_alpha, u_alpha, v_alpha])
    Y_abab = sqa_plus.tensor('Y_KWUV__abab', [r, k_alpha, w_beta, u_alpha, v_beta])

    term_h = [sqa_plus.term(1.0, [], [Y_aaaa, sqa_plus.creOp(u_alpha), sqa_plus.creOp(v_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_abab, sqa_plus.creOp(u_alpha), sqa_plus.creOp(v_beta),  sqa_plus.desOp(w_beta),  sqa_plus.desOp(k_alpha)])]

elif indices_string in ['q_ccaa']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    l_alpha = cvs_alpha_inds.new_index()
    l_beta  = cvs_beta_inds.new_index()

    w_alpha = act_alpha_inds.new_index()
    w_beta  = act_beta_inds.new_index()

    u_alpha = act_alpha_inds.new_index()
    u_beta  = act_beta_inds.new_index()


    Y_aaaa = sqa_plus.tensor('Y_KLWU__aaaa', [r, k_alpha, l_alpha, w_alpha, u_alpha])
    Y_abab = sqa_plus.tensor('Y_KLWU__abab', [r, k_alpha, l_beta,  w_alpha, u_beta])

    term_h = [sqa_plus.term(1.0, [], [Y_aaaa, sqa_plus.creOp(w_alpha), sqa_plus.creOp(u_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_abab, sqa_plus.creOp(w_alpha), sqa_plus.creOp(u_beta),  sqa_plus.desOp(l_beta),  sqa_plus.desOp(k_alpha)])]

elif indices_string in ['q_ccee']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    l_alpha = cvs_alpha_inds.new_index()
    l_beta  = cvs_beta_inds.new_index()

    c_alpha = vir_alpha_inds.new_index()
    c_beta  = vir_beta_inds.new_index()

    d_alpha = vir_alpha_inds.new_index()
    d_beta  = vir_beta_inds.new_index()


    Y_aaaa = sqa_plus.tensor('Y_KLCD__aaaa', [r, k_alpha, l_alpha, c_alpha, d_alpha])
    Y_abab = sqa_plus.tensor('Y_KLCD__abab', [r, k_alpha, l_beta,  c_alpha, d_beta])

    term_h = [sqa_plus.term(1.0, [], [Y_aaaa, sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_abab, sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_beta),  sqa_plus.desOp(l_beta),  sqa_plus.desOp(k_alpha)])]

elif indices_string in ['q_ccea']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    l_alpha = cvs_alpha_inds.new_index()
    l_beta  = cvs_beta_inds.new_index()

    c_alpha = vir_alpha_inds.new_index()
    c_beta  = vir_beta_inds.new_index()

    w_alpha = act_alpha_inds.new_index()
    w_beta  = act_beta_inds.new_index()

    Y_aaaa = sqa_plus.tensor('Y_KLCW__aaaa', [r, k_alpha, l_alpha, c_alpha, w_alpha])
    Y_abab = sqa_plus.tensor('Y_KLCW__abab', [r, k_alpha, l_beta,  c_alpha, w_beta])
    Y_abba = sqa_plus.tensor('Y_KLCW__abba', [r, k_alpha, l_beta,  c_beta,  w_alpha])


    term_h = [sqa_plus.term(1.0, [], [Y_aaaa, sqa_plus.creOp(c_alpha), sqa_plus.creOp(w_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_abab, sqa_plus.creOp(c_alpha), sqa_plus.creOp(w_beta),  sqa_plus.desOp(l_beta),  sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_abba, sqa_plus.creOp(c_beta), sqa_plus.creOp(w_alpha),  sqa_plus.desOp(l_beta),  sqa_plus.desOp(k_alpha)])]

elif indices_string in ['q_caee']:
   
    k_alpha = cvs_alpha_inds.new_index()
    k_beta  = cvs_beta_inds.new_index()

    w_alpha = act_alpha_inds.new_index()
    w_beta  = act_beta_inds.new_index()

    c_alpha = vir_alpha_inds.new_index()
    c_beta  = vir_beta_inds.new_index()

    d_alpha = vir_alpha_inds.new_index()
    d_beta  = vir_beta_inds.new_index()


    Y_aaaa = sqa_plus.tensor('Y_KWCD__aaaa', [r, k_alpha, w_alpha, c_alpha, d_alpha])
    Y_abab = sqa_plus.tensor('Y_KWCD__abab', [r, k_alpha, w_beta,  c_alpha, d_beta])
    Y_baab = sqa_plus.tensor('Y_KWCD__baab', [r, k_beta,  w_alpha,  c_alpha, d_beta])

    term_h = [sqa_plus.term(1.0, [], [Y_aaaa, sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_abab, sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_beta),  sqa_plus.desOp(w_beta),  sqa_plus.desOp(k_alpha)]),
              sqa_plus.term(1.0, [], [Y_baab, sqa_plus.creOp(c_alpha), sqa_plus.creOp(d_beta),  sqa_plus.desOp(w_alpha),  sqa_plus.desOp(k_beta)])]


## Evaluate transition moments contributions

TY_set = zip(term_q0, TY_index, final_indices_list)
for term_q0, index, final in TY_set:

    sqa_plus.options.print_header(final)
    sqa_plus.options.spin_integrated = True

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
    result = sqa_plus.genEinsum(expected_EE_TY_sa, 'TY' + index, 'RPQ')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
