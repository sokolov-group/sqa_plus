import sqa_plus
spin_integrated = False
explicit_spin_cases = False
use_legacy_order = True

amplitude_string = 't1_p1p'
s_block_string = 'S11'

import time
start = time.time()

print("\n----------------------------------------------------------------------------------")
print("Spin-Orbital {:} S".format(amplitude_string).center(82))
print("----------------------------------------------------------------------------------\n")

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
    if s_block_string == 'S12':
        print("## Create S12: a_X^\dag a_Y ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_x, des_y])]

        s_string = s_block_string
        final_indices_string = 'XY'

    elif s_block_string == 'S22':
        print("## Create S22: a_X^\dag a_Y^\dag a_Z a_W ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_x, des_y, cre_z, des_w])]

        s_string = s_block_string
        final_indices_string = 'XYWZ'

elif amplitude_string == 't1_m1':
    print("## Create S11: a_X^\dag a_Y ...\n")

    terms_S = [sqa_plus.term(1.0, [], [cre_x, des_y])]

    s_string = 'S_m1'
    final_indices_string = 'XY'

elif amplitude_string == 't1_p1':
    print("## Create S11: a_X a_Y^\dag ...\n")

    terms_S = [sqa_plus.term(1.0, [], [des_x, cre_y])]

    s_string = 'S_p1'
    final_indices_string = 'XY'

elif amplitude_string == 't1_m2':
    print("## Create S: a_X^\dag a_Y^\dag a_Z a_W ...\n")

    terms_S = [sqa_plus.term(1.0, [], [cre_x, cre_y, des_z, des_w])]

    s_string = 'S_m2'
    final_indices_string = 'XYWZ'

elif amplitude_string == 't1_p2':
    print("## Create S: a_X a_Y a_Z^\dag a_W^\dag ...\n")

    terms_S = [sqa_plus.term(1.0, [], [des_x, des_y, cre_z, cre_w])]

    s_string = 'S_p2'
    final_indices_string = 'XYWZ'

elif amplitude_string == 't1_m1p':
    if s_block_string == 'S11':
        print("## Create S11: a_X^\dag a_Y ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_x, des_y])]

        s_string = s_block_string
        final_indices_string = 'XY'

    elif s_block_string == 'S12':
        print("## Create S12: a_X^\dag a_Y^\dag a_Z a_W ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_x, cre_y, des_z, des_w])]

        s_string = s_block_string
        final_indices_string = 'XWZY'

    elif s_block_string == 'S22':
        print("## Create S22: a_U^\dag a_V^\dag a_X a_Y^\dag a_Z a_W ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_u, cre_v, des_x, cre_y, des_z, des_w])]

        s_string = s_block_string
        final_indices_string = 'UVXWZY'

elif amplitude_string == 't1_p1p':
    if s_block_string == 'S11':
        print("## Create S11: a_X a_Y^\dag ...\n")

        terms_S = [sqa_plus.term(1.0, [], [des_x, cre_y])]

        s_string = s_block_string
        final_indices_string = 'XY'

    elif s_block_string == 'S12':
        print("## Create S12: a_X a_Y^\dag a_Z a_W ...\n")

        terms_S = [sqa_plus.term(1.0, [], [des_x, cre_y, cre_z, des_w])]

        s_string = s_block_string
        final_indices_string = 'XWZY'

    elif s_block_string == 'S22':
        print("## Create S22: a_U^\dag a_V a_X a_Y^\dag a_Z^\dag a_W ...\n")

        terms_S = [sqa_plus.term(1.0, [], [cre_u, des_v, des_x, cre_y, cre_z, des_w])]

        s_string = s_block_string
        final_indices_string = 'UVXWZY'

for term_S in terms_S:
    print(term_S)

# Compute expected value of spin-orbital S matrix
print("## Compute expected value of spin-orbital S matrix ...")
terms_expected_S = sqa_plus.matrixBlock(terms_S, transRDM = False, legacy_order = use_legacy_order)

# Create Numpy einsum equations
terms_expected_S.sort()
result = sqa_plus.genEinsum(terms_expected_S, s_string, final_indices_string, rm_core_int = True, suffix = 'so')

print("\n-------------------------------- genEinsum equations --------------------------------\n")
for item in result:
    print(item)
print("\n-------------------------------------------------------------------------------------\n")

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))