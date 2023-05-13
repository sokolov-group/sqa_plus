import sqa_plus
sqa_plus.options.spin_integrated = True

amplitude_string = 't1_p1p'
v_block_string = 'V2'
spin_cases_string = 'ab_ba'

import time
start = time.time()

sqa_plus.options.print_header("Spin-Adapted {:} {:} {:}".format(amplitude_string, v_block_string, spin_cases_string))

# Create spin-integrated amplitude operator
## Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

i_alpha = sqa_plus.index('I', [tg_alpha, tg_cor])
i_beta  = sqa_plus.index('I', [tg_beta,  tg_cor])

j_alpha = sqa_plus.index('J', [tg_alpha, tg_cor])
j_beta  = sqa_plus.index('J', [tg_beta,  tg_cor])

x_alpha = sqa_plus.index('X', [tg_alpha, tg_act])
x_beta  = sqa_plus.index('X', [tg_beta,  tg_act])

y_alpha = sqa_plus.index('Y', [tg_alpha, tg_act])
y_beta  = sqa_plus.index('Y', [tg_beta,  tg_act])

z_alpha = sqa_plus.index('Z', [tg_alpha, tg_act])
z_beta  = sqa_plus.index('Z', [tg_beta,  tg_act])

a_alpha = sqa_plus.index('A', [tg_alpha, tg_vir])
a_beta  = sqa_plus.index('A', [tg_beta,  tg_vir])

b_alpha = sqa_plus.index('B', [tg_alpha, tg_vir])
b_beta  = sqa_plus.index('B', [tg_beta,  tg_vir])

## Define operator types
cre_i_alpha = sqa_plus.creOp(i_alpha)
cre_i_beta  = sqa_plus.creOp(i_beta)
des_i_alpha = sqa_plus.desOp(i_alpha)
des_i_beta  = sqa_plus.desOp(i_beta)

cre_j_alpha = sqa_plus.creOp(j_alpha)
cre_j_beta  = sqa_plus.creOp(j_beta)
des_j_alpha = sqa_plus.desOp(j_alpha)
des_j_beta  = sqa_plus.desOp(j_beta)

cre_x_alpha = sqa_plus.creOp(x_alpha)
cre_x_beta  = sqa_plus.creOp(x_beta)
des_x_alpha = sqa_plus.desOp(x_alpha)
des_x_beta  = sqa_plus.desOp(x_beta)

cre_y_alpha = sqa_plus.creOp(y_alpha)
cre_y_beta  = sqa_plus.creOp(y_beta)
des_y_alpha = sqa_plus.desOp(y_alpha)
des_y_beta  = sqa_plus.desOp(y_beta)

cre_z_alpha = sqa_plus.creOp(z_alpha)
cre_z_beta  = sqa_plus.creOp(z_beta)
des_z_alpha = sqa_plus.desOp(z_alpha)
des_z_beta  = sqa_plus.desOp(z_beta)

cre_a_alpha = sqa_plus.creOp(a_alpha)
cre_a_beta  = sqa_plus.creOp(a_beta)
des_a_alpha = sqa_plus.desOp(a_alpha)
des_a_beta  = sqa_plus.desOp(a_beta)

cre_b_alpha = sqa_plus.creOp(b_alpha)
cre_b_beta  = sqa_plus.creOp(b_beta)
des_b_alpha = sqa_plus.desOp(b_alpha)
des_b_beta  = sqa_plus.desOp(b_beta)

# Create spin-integrated V operator
print("# Create spin-integrated V operator ...")
terms_V = sqa_plus.Vperturbation()

## Define terms
if amplitude_string == 't1_0':
    print("## Create V: -1.0 * a_I^\dag a_J^\dag a_B a_A * V...\n")

    terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_j_beta, des_b_beta, des_a_alpha])]

    v_string = 'V1_0'
    final_indices_string = 'IJAB'

elif amplitude_string == 't1_0p':
    if v_block_string == 'V1':
        print("## Create V: -1.0 * a_I^\dag a_A * V...\n")

        terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, des_a_alpha])]

        v_string = 'V1_a_a'
        final_indices_string = 'IA'

    elif v_block_string == 'V2':
        if spin_cases_string == 'aa_aa':
            print("## Create V: -1.0 * a_I^\dag a_X^\dag a_Y a_A * V...\n")

            terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_x_alpha, des_y_alpha, des_a_alpha])]

            v_string = v_block_string + "_" + spin_cases_string
            final_indices_string = 'IAXY'

        elif spin_cases_string == 'aa_bb':
            print("## Create V: -1.0 * a_I^\dag a_X^\dag a_Y a_A * V...\n")

            terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_x_beta, des_y_beta, des_a_alpha])]

            v_string = v_block_string + "_" + spin_cases_string
            final_indices_string = 'IAXY'

elif amplitude_string == 't1_m1':
    print("## Create V: -1.0 * a_I^\dag a_X^\dag a_B a_A * V...\n")

    terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_x_beta, des_b_beta, des_a_alpha])]

    v_string = 'V1_m1'
    final_indices_string = 'IXAB'

elif amplitude_string == 't1_p1':
    print("## Create V: -1.0 * a_I^\dag a_J^\dag a_X a_A * V...\n")

    terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_j_beta, des_x_beta, des_a_alpha])]

    v_string = 'V1_p1'
    final_indices_string = 'IJAX'

elif amplitude_string == 't1_m2':
    print("## Create V: -1.0 * a_X^\dag a_Y^\dag a_B a_A * V...\n")

    terms_op = [sqa_plus.term(-1.0, [], [cre_x_alpha, cre_y_beta, des_b_beta, des_a_alpha])]

    v_string = 'V1_m2'
    final_indices_string = 'XYAB'

elif amplitude_string == 't1_p2':
    print("## Create V: -1.0 * a_I^\dag a_J^\dag a_Y a_X * V...\n")

    terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_j_beta, des_y_beta, des_x_alpha])]

    v_string = 'V1_p2'
    final_indices_string = 'IJXY'

elif amplitude_string == 't1_m1p':
    if v_block_string == 'V1':
        print("## Create V: -1.0 * a_X^\dag a_A * V...\n")

        terms_op = [sqa_plus.term(-1.0, [], [cre_x_alpha, des_a_alpha])]

        v_string = 'V1_a_a'
        final_indices_string = 'XA'

    elif v_block_string == 'V2':
        if spin_cases_string == 'aa_aa':
            print("## Create V: -1.0 * a_X^\dag a_Y^\dag a_Z a_A * V...\n")

            terms_op = [sqa_plus.term(-1.0, [], [cre_x_alpha, cre_y_alpha, des_z_alpha, des_a_alpha])]

            v_string = v_block_string + "_" + spin_cases_string
            final_indices_string = 'XYZA'

        elif spin_cases_string == 'ab_ba':
            print("## Create V: -1.0 * a_X^\dag a_Y^\dag a_Z a_A * V...\n")

            terms_op = [sqa_plus.term(-1.0, [], [cre_x_alpha, cre_y_beta, des_z_beta, des_a_alpha])]

            v_string = v_block_string + "_" + spin_cases_string
            final_indices_string = 'XYZA'

elif amplitude_string == 't1_p1p':
    if v_block_string == 'V1':
        print("## Create V: -1.0 * a_I^\dag a_X * V...\n")

        terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, des_x_alpha])]

        v_string = 'V1_a_a'
        final_indices_string = 'IX'

    elif v_block_string == 'V2':
        if spin_cases_string == 'aa_aa':
            print("## Create V: -1.0 * a_I^\dag a_X^\dag a_Y a_Z * V...\n")

            terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_x_alpha, des_y_alpha, des_z_alpha])]

            v_string = v_block_string + "_" + spin_cases_string
            final_indices_string = 'IXYZ'

        elif spin_cases_string == 'ab_ba':
            print("## Create V: -1.0 * a_I^\dag a_X^\dag a_Y a_Z * V...\n")

            terms_op = [sqa_plus.term(-1.0, [], [cre_i_alpha, cre_x_beta, des_y_beta, des_z_alpha])]

            v_string = v_block_string + "_" + spin_cases_string
            final_indices_string = 'IXYZ'

for term_op in terms_op:
    print(term_op)

# Multiply
print("\n## Multiply...")
terms_op_V = []
for term_op in terms_op:
    for term_V in terms_V:
        terms_op_V.append(sqa_plus.multiplyTerms(term_op, term_V))

# Compute expected value of spin-integrated V matrix
print("## Compute expected value of spin-integrated V matrix ...")
terms_expected_V = sqa_plus.matrixBlock(terms_op_V)

# Spin-Adaptating V matrix expression
terms_expected_V_sa = sqa_plus.convertSpinIntegratedToAdapted(terms_expected_V)

# Create Numpy einsum equations
terms_expected_V_sa.sort()
result = sqa_plus.genEinsum(terms_expected_V_sa, v_string, final_indices_string, rm_core_int = True, suffix = '')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
