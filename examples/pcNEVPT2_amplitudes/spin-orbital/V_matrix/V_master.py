import sqa_plus

amplitude_string = 't1_m1p'
v_block_string = 'V2'

import time
start = time.time()

sqa_plus.options.print_header("Spin-Orbital {:} V".format(amplitude_string))

# Create spin-orbital amplitude operator
## Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

i = sqa_plus.index('I', [tg_cor])
j = sqa_plus.index('J', [tg_cor])

x = sqa_plus.index('X', [tg_act])
y = sqa_plus.index('Y', [tg_act])
z = sqa_plus.index('Z', [tg_act])

a = sqa_plus.index('A', [tg_vir])
b = sqa_plus.index('B', [tg_vir])

## Define operator types
cre_i = sqa_plus.creOp(i)
des_i = sqa_plus.desOp(i)

cre_j = sqa_plus.creOp(j)
des_j = sqa_plus.desOp(j)

cre_x = sqa_plus.creOp(x)
des_x = sqa_plus.desOp(x)

cre_y = sqa_plus.creOp(y)
des_y = sqa_plus.desOp(y)

cre_z = sqa_plus.creOp(z)
des_z = sqa_plus.desOp(z)

cre_a = sqa_plus.creOp(a)
des_a = sqa_plus.desOp(a)

cre_b = sqa_plus.creOp(b)
des_b = sqa_plus.desOp(b)

# Create spin-orbital V operator
print("# Create spin-orbital V operator ...")
terms_V = sqa_plus.Vperturbation()

## Define terms
if amplitude_string == 't1_0':
    print("## Create V: a_I^\dag a_J^\dag a_B a_A * V...\n")

    terms_op = [sqa_plus.term(1.0, [], [cre_i, cre_j, des_b, des_a])]

    v_string = 'V1_0'
    final_indices_string = 'IJAB'

elif amplitude_string == 't1_0p':
    if v_block_string == 'V1':
        print("## Create V: - 1.0 * a_I^\dag a_A * V...\n")

        terms_op = [sqa_plus.term(- 1.0, [], [cre_i, des_a])]

        v_string = 'V0p'
        final_indices_string = 'IA'

    elif v_block_string == 'V2':
        print("## Create V: - 1.0 * a_I^\dag a_X^\dag a_Y a_A * V...\n")

        terms_op = [sqa_plus.term(- 1.0, [], [cre_i, cre_x, des_y, des_a])]

        v_string = 'V0p'
        final_indices_string = 'IAXY'

elif amplitude_string == 't1_m1':
    print("## Create V: a_I^\dag a_X^\dag a_B a_A * V...\n")

    terms_op = [sqa_plus.term(1.0, [], [cre_i, cre_x, des_b, des_a])]

    v_string = 'V1_m1'
    final_indices_string = 'IXAB'

elif amplitude_string == 't1_p1':
    print("## Create V: a_I^\dag a_J^\dag a_X a_A * V...\n")

    terms_op = [sqa_plus.term(1.0, [], [cre_i, cre_j, des_x, des_a])]

    v_string = 'V1_p1'
    final_indices_string = 'IJAX'

elif amplitude_string == 't1_m2':
    print("## Create V: a_X^\dag a_Y^\dag a_B a_A * V...\n")

    terms_op = [sqa_plus.term(1.0, [], [cre_x, cre_y, des_b, des_a])]

    v_string = 'V1_m2'
    final_indices_string = 'XYAB'

elif amplitude_string == 't1_p2':
    print("## Create V: a_I^\dag a_J^\dag a_Y a_X * V...\n")

    terms_op = [sqa_plus.term(1.0, [], [cre_i, cre_j, des_y, des_x])]

    v_string = 'V1_p2'
    final_indices_string = 'IJXY'

elif amplitude_string == 't1_m1p':
    if v_block_string == 'V1':
        print("## Create V: - 1.0 * a_X^\dag a_A * V...\n")

        terms_op = [sqa_plus.term(- 1.0, [], [cre_x, des_a])]

        v_string = v_block_string
        final_indices_string = 'XA'

    elif v_block_string == 'V2':
        print("## Create V: - 1.0 * a_X^\dag a_Y^\dag a_Z a_A * V...\n")

        terms_op = [sqa_plus.term(- 1.0, [], [cre_x, cre_y, des_z, des_a])]

        v_string = v_block_string
        final_indices_string = 'XYZA'

elif amplitude_string == 't1_p1p':
    if v_block_string == 'V1':
        print("## Create V: - 1.0 * a_I^\dag a_X * V...\n")

        terms_op = [sqa_plus.term(- 1.0, [], [cre_i, des_x])]

        v_string = v_block_string
        final_indices_string = 'IX'

    elif v_block_string == 'V2':
        print("## Create V: - 1.0 * a_I^\dag a_X^\dag a_Y a_Z * V...\n")

        terms_op = [sqa_plus.term(- 1.0, [], [cre_i, cre_x, des_y, des_z])]

        v_string = v_block_string
        final_indices_string = 'IXYZ'

for term_op in terms_op:
    print(term_op)

# Multiply
print("\n## Multiply...")
terms_op_V = []
for term_op in terms_op:
    for term_V in terms_V:
        terms_op_V.append(sqa_plus.multiplyTerms(term_op, term_V))

# Compute expected value of spin-orbital V matrix
print("## Compute expected value of spin-orbital V matrix ...")
terms_expected_V = sqa_plus.matrixBlock(terms_op_V)

# Create Numpy einsum equations
terms_expected_V.sort()
result = sqa_plus.genEinsum(terms_expected_V, v_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
