import sqa_plus as sqa 
import time
from sqa_plus import sqaEinsum
start = time.time()

# Generating operators
print("\n## Generating operators ...\n")

## Define indices
dummy = True

tg_cor = sqa.options.core_type
tg_vir = sqa.options.virtual_type

tg_alpha = sqa.options.alpha_type
tg_beta = sqa.options.beta_type

# Core dummy indices
cor_alpha_inds = sqa.dummyIndexList('cc%ia', [tg_alpha, tg_cor], dummy)
cor_beta_inds  = sqa.dummyIndexList('cc%ib', [tg_beta,  tg_cor], dummy)

# Virtual dummy indices
vir_alpha_inds = sqa.dummyIndexList('vv%ia', [tg_alpha, tg_vir], dummy)
vir_beta_inds  = sqa.dummyIndexList('vv%ib', [tg_beta,  tg_vir], dummy)

### External Indices
I_alpha =  sqa.index('I', [tg_alpha, tg_cor])
i_beta  =  sqa.index('i', [tg_beta,  tg_cor])

J_alpha =  sqa.index('J', [tg_alpha, tg_cor])
j_beta  =  sqa.index('j', [tg_beta,  tg_cor])

L_alpha =  sqa.index('L', [tg_alpha, tg_cor])
l_beta  =  sqa.index('l', [tg_beta,  tg_cor])

K_alpha =  sqa.index('K', [tg_alpha, tg_cor])
k_beta  =  sqa.index('k', [tg_beta,  tg_cor])

A_alpha =  sqa.index('A', [tg_alpha, tg_vir])
a_beta  =  sqa.index('a', [tg_beta,  tg_vir])

B_alpha =  sqa.index('B', [tg_alpha, tg_vir])
b_beta  =  sqa.index('b', [tg_beta,  tg_vir])

D_alpha =  sqa.index('D', [tg_alpha, tg_vir])
d_beta  =  sqa.index('d', [tg_beta,  tg_vir])

C_alpha =  sqa.index('C', [tg_alpha, tg_vir])
c_beta  =  sqa.index('c', [tg_beta,  tg_vir])

## Define O operators
cre_I_alpha = sqa.creOp(I_alpha)
cre_J_alpha = sqa.creOp(J_alpha)
des_L_alpha = sqa.desOp(L_alpha)
des_K_alpha = sqa.desOp(K_alpha)

cre_i_beta = sqa.creOp(i_beta)
cre_j_beta = sqa.creOp(j_beta)
des_l_beta = sqa.desOp(l_beta)
des_k_beta = sqa.desOp(k_beta)


cre_A_alpha = sqa.creOp(A_alpha)
cre_B_alpha = sqa.creOp(B_alpha)
des_D_alpha = sqa.desOp(D_alpha)
des_C_alpha = sqa.desOp(C_alpha)

cre_a_beta = sqa.creOp(a_beta)
cre_b_beta = sqa.creOp(b_beta)
des_d_beta = sqa.desOp(d_beta)
des_c_beta = sqa.desOp(c_beta)

# Eigenvectors definition
## Define Y0 matrices

### Internal Indices
e_alpha = vir_alpha_inds.new_index()
r_alpha = vir_alpha_inds.new_index()
f_alpha = vir_alpha_inds.new_index()
s_alpha = vir_alpha_inds.new_index()
m_alpha = cor_alpha_inds.new_index()
p_alpha = cor_alpha_inds.new_index()
q_alpha = cor_alpha_inds.new_index()
n_alpha = cor_alpha_inds.new_index()

e_beta = vir_beta_inds.new_index()
r_beta = vir_beta_inds.new_index()
f_beta = vir_beta_inds.new_index()
s_beta = vir_beta_inds.new_index()
m_beta = cor_beta_inds.new_index()
p_beta = cor_beta_inds.new_index()
q_beta = cor_beta_inds.new_index()
n_beta = cor_beta_inds.new_index()

### Internal Indices operators
cre_e_alpha = sqa.creOp(e_alpha)
cre_f_alpha = sqa.creOp(f_alpha)
cre_q_alpha = sqa.creOp(q_alpha)
cre_p_alpha = sqa.creOp(p_alpha)
des_m_alpha = sqa.desOp(m_alpha)
des_n_alpha = sqa.desOp(n_alpha)
des_s_alpha = sqa.desOp(s_alpha)
des_r_alpha = sqa.desOp(r_alpha)

cre_e_beta = sqa.creOp(e_beta)
cre_f_beta = sqa.creOp(f_beta)
cre_q_beta = sqa.creOp(q_beta)
cre_p_beta = sqa.creOp(p_beta)
des_m_beta = sqa.desOp(m_beta)
des_n_beta = sqa.desOp(n_beta)
des_s_beta = sqa.desOp(s_beta)
des_r_beta = sqa.desOp(r_beta)


## Defining M operators

a_id_dag_aa = [cre_I_alpha, des_D_alpha]
a_id_dag_aa_terms = sqa.term(1.0, [], a_id_dag_aa)

a_id_dag_bb = [cre_i_beta, des_d_beta]
a_id_dag_bb_terms = sqa.term(1.0, [], a_id_dag_bb)

a_la_aa = [cre_A_alpha, des_L_alpha]
a_la_aa_terms = sqa.term(1.0, [], a_la_aa)

a_la_bb = [cre_a_beta, des_l_beta]
a_la_bb_terms = sqa.term(1.0, [], a_la_bb)

a_ijcd_dag_aa = [cre_I_alpha, cre_J_alpha, des_D_alpha, des_C_alpha]
a_ijcd_dag_aa_terms = [sqa.term(1.0, [], a_ijcd_dag_aa)]

a_ijcd_dag_ab = [cre_I_alpha, cre_j_beta, des_d_beta, des_C_alpha]
a_ijcd_dag_ab_terms = sqa.term(1.0, [], a_ijcd_dag_ab)

a_ijcd_dag_bb = [cre_i_beta, cre_j_beta, des_d_beta, des_c_beta]
a_ijcd_dag_bb_terms = sqa.term(1.0, [], a_ijcd_dag_bb)


Y0_x_aa = [sqa.tensor('Y', [m_alpha, e_alpha], []), cre_e_alpha, des_m_alpha]
Y0_x_bb = [sqa.tensor('Y', [m_beta,  e_beta], []), cre_e_beta, des_m_beta]

Y0_x_aa_terms = [sqa.term(1.0, [], Y0_x_aa)]
Y0_x_bb_terms = sqa.term(1.0, [], Y0_x_bb)

Y_symm = [sqa.symmetry((1,0,2,3), -1), sqa.symmetry((0,1,3,2), -1)]

Y1_x_aa = [sqa.tensor('Y', [ m_alpha, n_alpha, e_alpha, f_alpha], Y_symm), cre_e_alpha, cre_f_alpha ,des_n_alpha, des_m_alpha]
Y1_x_aa_terms = sqa.term(0.25, [], Y1_x_aa)

Y1_x_ab = [sqa.tensor('Y', [ m_alpha, n_beta, e_alpha, f_beta],[]), cre_e_alpha, cre_f_beta ,des_n_beta, des_m_alpha]
Y1_x_ab_terms = sqa.term(1.0, [], Y1_x_ab)

Y1_x_bb = [sqa.tensor('Y', [m_beta, n_beta, e_beta, f_beta], Y_symm), cre_e_beta, cre_f_beta ,des_n_beta, des_m_beta]
Y1_x_bb_terms = sqa.term(0.25, [], Y1_x_bb)

O_op = [[cre_I_alpha, des_L_alpha],
        [cre_A_alpha, des_C_alpha],
        [cre_I_alpha, des_C_alpha],
        [cre_A_alpha, des_L_alpha]]

#O_op =  [[cre_i_beta, des_l_beta],
#        [cre_a_beta, des_c_beta],
#        [cre_i_beta, des_c_beta],
#        [cre_a_beta, des_l_beta]]
num = 0
order = 1
for op in O_op:
    
    termm = sqa.term(1.0, [], op)
    from sqaHeff import Heff
    TY_3 =  Heff(order, termm)
    #TY_3 =  sqa.Heff(order, termm)
    
    term_ = sqa.commutator(TY_3, Y0_x_aa_terms,  combine = True)

        # Expected value
    expected_G = sqa.matrixBlock(term_)
    print("\n------------------------------- matrixBlock equations -------------------------------\n")
    for item in expected_G:
        print(item)
    print("\n-------------------------------------------------------------------------------------\n")
    
    # Generating Numpy einsum equations
    external_TY = ["IL","AC","IC","AL"]
#    external_TY = ["il","ac","ic","al"]
    result = sqa.genEinsum(expected_G, 'TY_a', external_TY[num], spin_integrated = True , use_spin_integrated_tensors = True, rm_core_int = True)
    num += 1    
    print("\n-------------------------------- genEinsum equations --------------------------------\n")
    for item in result:
        print(item)
    print("\n-------------------------------------------------------------------------------------\n")
    
    end = time.time()
    print("> Total elapsed time: {:.2f} seconds.".format(end - start))
    
#    term__ = sqa.commutator(TY_3, Y0_x_bb_terms,  combine = True)
#
#        # Expected value
#    expected_G = sqa.matrixBlock(term__)
#    print("\n------------------------------- matrixBlock equations -------------------------------\n")
#    for item in expected_G:
#        print(item)
#    print("\n-------------------------------------------------------------------------------------\n")
#    
#    # Generating Numpy einsum equations
##    external_TY = ["IL","AC","IC","AL"]
#    external_TY = ["il","ac","ic","al"]
#    result = sqa.genEinsum(expected_G, 'TY_b', external_TY[num], spin_integrated = True, use_spin_integrated_tensors = True, rm_core_int = True)
#    num += 1
#    
#    print("\n-------------------------------- genEinsum equations --------------------------------\n")
#    for item in result:
#        print(item)
#    print("\n-------------------------------------------------------------------------------------\n")
#    
#    end = time.time()
#    print("> Total elapsed time: {:.2f} seconds.".format(end - start))
