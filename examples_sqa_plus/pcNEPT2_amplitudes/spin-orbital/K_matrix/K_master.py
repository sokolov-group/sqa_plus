import sqa_plus
spin_integrated = False
explicit_spin_cases = False
use_legacy_order = True

amplitude_string = 't1_p1p'
k_block_string = 'K22'

import time
start = time.time()

print("\n----------------------------------------------------------------------------------")
print("Spin-Orbital {:} K".format(amplitude_string).center(82))
print("----------------------------------------------------------------------------------\n")

# Create spin-orbital Dyall Hamiltonian
print("# Create spin-orbital Dyall Hamiltonian ...")
indices_lists = sqa_plus.create_dummy_indices_list(spin_integrated)
terms_Heff = sqa_plus.dyallH_act(indices_lists, spin_integrated, explicit_spin_cases)

# Create spin-orbital amplitude operator
## Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

u = sqa_plus.index('U', [tg_act])
v = sqa_plus.index('V', [tg_act])
x = sqa_plus.index('X', [tg_act])
y = sqa_plus.index('Y', [tg_act])
w = sqa_plus.index('W', [tg_act])
z = sqa_plus.index('Z', [tg_act])

## Define operator types
cre_u = sqa_plus.creOp(u)
des_u = sqa_plus.desOp(u)

cre_v = sqa_plus.creOp(v)
des_v = sqa_plus.desOp(v)

cre_x = sqa_plus.creOp(x)
des_x = sqa_plus.desOp(x)

cre_y = sqa_plus.creOp(y)
des_y = sqa_plus.desOp(y)

cre_w = sqa_plus.creOp(w)
des_w = sqa_plus.desOp(w)

cre_z = sqa_plus.creOp(z)
des_z = sqa_plus.desOp(z)

## Define terms
if amplitude_string == 't1_0p':
    print("## Create K_caca: a_X^\dag a_Y [H_{act}, a_Z^\dag a_W] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [cre_x, des_y])
    term_r_op = sqa_plus.term(1.0, [], [cre_z, des_w])

    k_string = 'K_caca'
    final_indices_string = 'XYWZ'

elif amplitude_string == 't1_m1':
    print("## Create K_ca: a_X^\dag [H_{act}, a_Y] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [cre_x])
    term_r_op = sqa_plus.term(1.0, [], [des_y])

    k_string = 'K_ca'
    final_indices_string = 'XY'

elif amplitude_string == 't1_p1':
    print("## Create K_ac: a_X [H_{act}, a_Y^\dag] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [des_x])
    term_r_op = sqa_plus.term(1.0, [], [cre_y])

    k_string = 'K_ac'
    final_indices_string = 'XY'

elif amplitude_string == 't1_m2':
    print("## Create K_ccaa: a_X^\dag a_Y^\dag [H_{act}, a_W a_Z] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [cre_x, cre_y])
    term_r_op = sqa_plus.term(1.0, [], [des_w, des_z])

    k_string = 'K_ccaa'
    final_indices_string = 'XYZW'

elif amplitude_string == 't1_p2':
    print("## Create K_aacc: a_Z a_W [H_{act}, a_X^\dag a_Y^\dag] ...\n")

    term_l_op = sqa_plus.term(1.0, [], [des_z, des_w])
    term_r_op = sqa_plus.term(1.0, [], [cre_x, cre_y])

    k_string = 'K_aacc'
    final_indices_string = 'ZWXY'

elif amplitude_string == 't1_m1p':
    if k_block_string == 'K11':
        print("## Create K11: a_X^\dag [H_{act}, a_Y] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_x])
        term_r_op = sqa_plus.term(1.0, [], [des_y])

        k_string = k_block_string
        final_indices_string = 'XY'

    elif k_block_string == 'K12':
        print("## Create K12: a_X^\dag [H_{act}, a_Z^\dag a_W a_Y] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_x])
        term_r_op = sqa_plus.term(1.0, [], [cre_z, des_w, des_y])

        k_string = k_block_string
        final_indices_string = 'XZWY'

    elif k_block_string == 'K22':
        print("## Create K22: a_X^\dag a_Z^\dag a_W [H_{act}, a_U^\dag a_V a_Y] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_x, cre_z, des_w])
        term_r_op = sqa_plus.term(1.0, [], [cre_u, des_v, des_y])

        k_string = k_block_string
        final_indices_string = 'XZWUVY'

elif amplitude_string == 't1_p1p':
    if k_block_string == 'K11':
        print("## Create K11: a_X [H_{act}, a_Y^\dag] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [des_x])
        term_r_op = sqa_plus.term(1.0, [], [cre_y])

        k_string = k_block_string
        final_indices_string = 'XY'

    elif k_block_string == 'K12':
        print("## Create K12: a_X [H_{act}, a_Y^\dag a_W^\dag a_Z] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [des_x])
        term_r_op = sqa_plus.term(1.0, [], [cre_y, cre_w, des_z])

        k_string = k_block_string
        final_indices_string = 'XYWZ'

    elif k_block_string == 'K22':
        print("## Create K22: a_U^\dag a_V a_X [H_{act}, a_Y^\dag a_W^\dag a_Z] ...\n")

        term_l_op = sqa_plus.term(1.0, [], [cre_u, des_v, des_x])
        term_r_op = sqa_plus.term(1.0, [], [cre_y, cre_w, des_z])

        k_string = k_block_string
        final_indices_string = 'UVXYWZ'

print(term_l_op)
print(term_r_op)

## Compute the commutator
print("\n## Calculating the commutator ...")
terms_r_com = sqa_plus.commutator(terms_Heff, term_r_op)

print("\n## Multiply ...")
terms_K = []
for term_r_com in terms_r_com:
    terms_K.append(sqa_plus.multiplyTerms(term_l_op, term_r_com))

# Compute expected value of spin-orbital K matrix
print("## Compute expected value of spin-orbital K matrix ...")
terms_expected_K = sqa_plus.matrixBlock(terms_K, transRDM = False, legacy_order = use_legacy_order)

# Create Numpy einsum equations
terms_expected_K.sort()
result = sqa_plus.genEinsum(terms_expected_K, k_string, final_indices_string, rm_core_int = True, suffix = 'so')

print("\n-------------------------------- genEinsum equations --------------------------------\n")
for item in result:
    print(item)
print("\n-------------------------------------------------------------------------------------\n")

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))