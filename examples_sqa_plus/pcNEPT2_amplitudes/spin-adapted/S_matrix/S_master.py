import sqa_plus
spin_integrated = True
explicit_spin_cases = True

amplitude_string = 't1_p1p'
s_block_string = 'S22'
spin_cases_string = 'bba_bba'

import time
start = time.time()

print("\n----------------------------------------------------------------------------------")
print("Spin-Adapted {:} S".format(amplitude_string).center(82))
print("----------------------------------------------------------------------------------\n")

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

## Define terms
if amplitude_string == 't1_0p':
    if s_block_string == 'S12':
        print("## Create S12: a_X^\dag a_Y ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_x_alpha, des_y_alpha])]

        s_string = 'S12_a_a'
        final_indices_string = 'XY'

    elif s_block_string == 'S22':
        if spin_cases_string == 'aa_aa':
            print("## Create S22: a_X^\dag a_Y^\dag a_Z a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [cre_x_alpha, des_y_alpha, cre_z_alpha, des_w_alpha])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'XYWZ'

        elif spin_cases_string == 'aa_bb':
            print("## Create S22: a_X^\dag a_Y^\dag a_Z a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [cre_x_alpha, des_y_alpha, cre_z_beta, des_w_beta])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'XYWZ'

elif amplitude_string == 't1_m1':
    print("## Create S11: a_X^\dag a_Y ...\n")

    terms_S = [sqa_plus.term(1.0, [], [cre_x_alpha, des_y_alpha])]

    s_string = 'S_m1'
    final_indices_string = 'XY'

elif amplitude_string == 't1_p1':
    print("## Create S11: a_X a_Y^\dag ...\n")

    terms_S = [sqa_plus.term(1.0, [], [des_x_alpha, cre_y_alpha])]

    s_string = 'S_p1'
    final_indices_string = 'XY'

elif amplitude_string == 't1_m2':
    print("## Create S: a_X^\dag a_Y^\dag a_Z a_W ...\n")

    terms_S = [sqa_plus.term(1.0, [], [cre_x_alpha, cre_y_beta, des_z_beta, des_w_alpha])]

    s_string = 'S_m2'
    final_indices_string = 'XYWZ'

elif amplitude_string == 't1_p2':
    print("## Create S: a_X a_Y a_Z^\dag a_W^\dag ...\n")

    terms_S = [sqa_plus.term(1.0, [], [des_x_alpha, des_y_beta, cre_z_beta, cre_w_alpha])]

    s_string = 'S_p2'
    final_indices_string = 'XYWZ'

elif amplitude_string == 't1_m1p':
    if s_block_string == 'S11':
        print("## Create S11: a_X^\dag a_Y ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_x_alpha, des_y_alpha])]

        s_string = 'S11_a_a'
        final_indices_string = 'XY'

    elif s_block_string == 'S12':
        if spin_cases_string == 'a_abb':
            print("## Create S12: a_X^\dag a_Y^\dag a_Z a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [cre_x_alpha, cre_y_beta, des_z_beta, des_w_alpha])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'XWZY'

    elif s_block_string == 'S22':
        if spin_cases_string == 'aaa_aaa':
            print("## Create S22: a_U^\dag a_V^\dag a_X a_Y^\dag a_Z a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [cre_u_alpha, cre_v_alpha, des_x_alpha, cre_y_alpha, des_z_alpha, des_w_alpha])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

        elif spin_cases_string == 'abb_abb':
            print("## Create S22: a_U^\dag a_V^\dag a_X a_Y^\dag a_Z a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [cre_u_alpha, cre_v_beta, des_x_beta, cre_y_beta, des_z_beta, des_w_alpha])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

elif amplitude_string == 't1_p1p':
    if s_block_string == 'S11':
        print("## Create S11: a_X a_Y^\dag ...\n")

        terms_S = [sqa_plus.term(1.0, [], [des_x_alpha, cre_y_alpha])]

        s_string = 'S11_a_a'
        final_indices_string = 'XY'

    elif s_block_string == 'S12':
        if spin_cases_string == 'a_bba':
            print("## Create S12: a_X a_Y^\dag a_Z a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [des_x_alpha, cre_y_alpha, cre_z_beta, des_w_beta])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'XWZY'

    elif s_block_string == 'S22':
        if spin_cases_string == 'aaa_aaa':
            print("## Create S22: a_U^\dag a_V a_X a_Y^\dag a_Z^\dag a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [cre_u_alpha, des_v_alpha, des_x_alpha, cre_y_alpha, cre_z_alpha, des_w_alpha])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

        elif spin_cases_string == 'bba_bba':
            print("## Create S22: a_U^\dag a_V^\dag a_X a_Y^\dag a_Z a_W ...\n")

            terms_S = [sqa_plus.term(1.0, [], [cre_u_beta, des_v_beta, des_x_alpha, cre_y_alpha, cre_z_beta, des_w_beta])]

            s_string = s_block_string + "_" + spin_cases_string
            final_indices_string = 'UVXWZY'

for term_S in terms_S:
    print(term_S)

# Compute expected value of spin-integrated S matrix
print("## Compute expected value of spin-integrated S matrix ...")
terms_expected_S = sqa_plus.matrixBlock(terms_S)

# Spin-Adaptating S matrix expression
terms_expected_S_sa = sqa_plus.convertSpinIntegratedToAdapted(terms_expected_S)

# Create Numpy einsum equations
terms_expected_S_sa.sort()
result = sqa_plus.genEinsum(terms_expected_S_sa, s_string, final_indices_string, rm_core_int = True, suffix = '')

print("\n-------------------------------- genEinsum equations --------------------------------\n")
for item in result:
    print(item)
print("\n-------------------------------------------------------------------------------------\n")

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))