import sqa_plus
sqa_plus.options.spin_integrated = True

amplitude_string = 't1_p1p'
k_block_string = 'K22'
spin_cases_string = 'bba_bba'

import time
start = time.time()

sqa_plus.options.print_header("Spin-Adapted {:} {:} {:}".format(amplitude_string, k_block_string, spin_cases_string))

# Create spin-integrated amplitude operator
## Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

u_alpha = sqa_plus.index('U', [tg_alpha, tg_act])
u_beta  = sqa_plus.index('U', [tg_beta,  tg_act])

v_alpha = sqa_plus.index('V', [tg_alpha, tg_act])
v_beta  = sqa_plus.index('V', [tg_beta,  tg_act])

x_alpha = sqa_plus.index('X', [tg_alpha, tg_act])
x_beta  = sqa_plus.index('X', [tg_beta,  tg_act])

y_alpha = sqa_plus.index('Y', [tg_alpha, tg_act])
y_beta  = sqa_plus.index('Y', [tg_beta,  tg_act])

w_alpha = sqa_plus.index('W', [tg_alpha, tg_act])
w_beta  = sqa_plus.index('W', [tg_beta,  tg_act])

z_alpha = sqa_plus.index('Z', [tg_alpha, tg_act])
z_beta  = sqa_plus.index('Z', [tg_beta,  tg_act])

## Define operator types
cre_u_alpha = sqa_plus.creOp(u_alpha)
cre_u_beta  = sqa_plus.creOp(u_beta)
des_u_alpha = sqa_plus.desOp(u_alpha)
des_u_beta  = sqa_plus.desOp(u_beta)

cre_v_alpha = sqa_plus.creOp(v_alpha)
cre_v_beta  = sqa_plus.creOp(v_beta)
des_v_alpha = sqa_plus.desOp(v_alpha)
des_v_beta  = sqa_plus.desOp(v_beta)

cre_x_alpha = sqa_plus.creOp(x_alpha)
cre_x_beta  = sqa_plus.creOp(x_beta)
des_x_alpha = sqa_plus.desOp(x_alpha)
des_x_beta  = sqa_plus.desOp(x_beta)

cre_y_alpha = sqa_plus.creOp(y_alpha)
cre_y_beta  = sqa_plus.creOp(y_beta)
des_y_alpha = sqa_plus.desOp(y_alpha)
des_y_beta  = sqa_plus.desOp(y_beta)

cre_w_alpha = sqa_plus.creOp(w_alpha)
cre_w_beta  = sqa_plus.creOp(w_beta)
des_w_alpha = sqa_plus.desOp(w_alpha)
des_w_beta  = sqa_plus.desOp(w_beta)

cre_z_alpha = sqa_plus.creOp(z_alpha)
cre_z_beta  = sqa_plus.creOp(z_beta)
des_z_alpha = sqa_plus.desOp(z_alpha)
des_z_beta  = sqa_plus.desOp(z_beta)

# Create spin-integrated Dyall Hamiltonian
print("# Create spin-integrated Dyall Hamiltonian ...")
terms_Heff = sqa_plus.dyallH_act()

## Define terms
if amplitude_string == 't1_0p':
    if spin_cases_string == 'aa_aa':
        print("## Create K_caca: a_X^\dag a_Y [H_{act}, a_Z^\dag a_W] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_x_alpha, des_y_alpha])
        term_r_op = sqa_plus.term(1.0, [], [cre_z_alpha, des_w_alpha])

        k_string = 'K_caca_' + spin_cases_string
        final_indices_string = 'XYWZ'

    elif spin_cases_string == 'aa_bb':
        print("## Create K_caca: a_X^\dag a_Y [H_{act}, a_Z^\dag a_W] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_x_alpha, des_y_alpha])
        term_r_op = sqa_plus.term(1.0, [], [cre_z_beta, des_w_beta])

        k_string = 'K_caca_' + spin_cases_string
        final_indices_string = 'XYWZ'

elif amplitude_string == 't1_m1':
    print("## Create K_ca: a_X^\dag [H_{act}, a_Y] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [cre_x_alpha])
    term_r_op = sqa_plus.term(1.0, [], [des_y_alpha])

    k_string = 'K_ca'
    final_indices_string = 'XY'

elif amplitude_string == 't1_p1':
    print("## Create K_ac: a_X [H_{act}, a_Y^\dag] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [des_x_alpha])
    term_r_op = sqa_plus.term(1.0, [], [cre_y_alpha])

    k_string = 'K_ac'
    final_indices_string = 'XY'

elif amplitude_string == 't1_m2':
    print("## Create K_ccaa: a_X^\dag a_Y^\dag [H_{act}, a_Z a_W] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [cre_x_alpha, cre_y_beta])
    term_r_op = sqa_plus.term(1.0, [], [des_z_beta,  des_w_alpha])

    k_string = 'K_ccaa'
    final_indices_string = 'XYWZ'

elif amplitude_string == 't1_p2':
    print("## Create K_aacc: a_X a_Y [H_{act}, a_Z^\dag a_W^\dag] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [des_x_alpha, des_y_beta])
    term_r_op = sqa_plus.term(1.0, [], [cre_z_beta,  cre_w_alpha])

    k_string = 'K_aacc'
    final_indices_string = 'XYWZ'

elif amplitude_string == 't1_m1p':
    if k_block_string == 'K11':
        print("## Create K11: a_X^\dag [H_{act}, a_Y] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_x_alpha])
        term_r_op = sqa_plus.term(1.0, [], [des_y_alpha])

        k_string = 'K11_a_a'
        final_indices_string = 'XY'

    elif k_block_string == 'K12':
        print("## Create K12: a_X^\dag [H_{act}, a_Y^\dag a_Z a_W] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_x_alpha])
        term_r_op = sqa_plus.term(1.0, [], [cre_y_beta, des_z_beta, des_w_alpha])

        k_string = 'K12_a_abb'
        final_indices_string = 'XWZY'

    elif k_block_string == 'K22':
        if spin_cases_string == 'aaa_aaa':
            print("## Create K22: a_U^\dag a_V^\dag a_X [H_{act}, a_Y^\dag a_Z a_W] ...\n")

            term_l_op = sqa_plus.term(1.0, [], [cre_u_alpha, cre_v_alpha, des_x_alpha])
            term_r_op = sqa_plus.term(1.0, [], [cre_y_alpha, des_z_alpha, des_w_alpha])

            k_string = k_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

        elif spin_cases_string == 'abb_abb':
            print("## Create K22: a_U^\dag a_V^\dag a_X [H_{act}, a_Y^\dag a_Z a_W] ...\n")

            term_l_op = sqa_plus.term(1.0, [], [cre_u_alpha, cre_v_beta,  des_x_beta])
            term_r_op = sqa_plus.term(1.0, [], [cre_y_beta,  des_z_beta,  des_w_alpha])

            k_string = k_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

elif amplitude_string == 't1_p1p':
    if k_block_string == 'K11':
        print("## Create K11: a_X [H_{act}, a_Y^\dag] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [des_x_alpha])
        term_r_op = sqa_plus.term(1.0, [], [cre_y_alpha])

        k_string = 'K11_a_a'
        final_indices_string = 'XY'

    elif k_block_string == 'K12':
        print("## Create K12: a_X [H_{act}, a_Y^\dag a_Z^\dag a_W] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [des_x_alpha])
        term_r_op = sqa_plus.term(1.0, [], [cre_y_alpha, cre_z_beta, des_w_beta])

        k_string = 'K12_a_bba'
        final_indices_string = 'XWZY'

    elif k_block_string == 'K22':
        if spin_cases_string == 'aaa_aaa':
            print("## Create K22: a_U^\dag a_V a_X [H_{act}, a_Y^\dag a_Z^\dag a_W] ...\n")

            term_l_op = sqa_plus.term(1.0, [], [cre_u_alpha, des_v_alpha, des_x_alpha])
            term_r_op = sqa_plus.term(1.0, [], [cre_y_alpha, cre_z_alpha, des_w_alpha])

            k_string = k_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

        elif spin_cases_string == 'bba_bba':
            print("## Create K22: a_U^\dag a_V a_X [H_{act}, a_Y^\dag a_Z^\dag a_W] ...\n")

            term_l_op = sqa_plus.term(1.0, [], [cre_u_beta,  des_v_beta, des_x_alpha])
            term_r_op = sqa_plus.term(1.0, [], [cre_y_alpha, cre_z_beta, des_w_beta])

            k_string = k_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

print(term_l_op)
print(term_r_op)

## Compute the commutator
print("\n## Calculating the commutator ...")
terms_r_com = sqa_plus.commutator(terms_Heff, term_r_op)

print("\n## Multiply ...")
terms_K = []
for term_r_com in terms_r_com:
    terms_K.append(sqa_plus.multiplyTerms(term_l_op, term_r_com))

# Compute expected value of spin-integrated K matrix
print("## Compute expected value of spin-integrated K matrix ...")
terms_expected_K = sqa_plus.matrixBlock(terms_K)

# Spin-Adaptating K matrix expression
terms_expected_K_sa = sqa_plus.convertSpinIntegratedToAdapted(terms_expected_K)

# Create Numpy einsum equations
terms_expected_K_sa.sort()
result = sqa_plus.genEinsum(terms_expected_K_sa, k_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))