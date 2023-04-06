import sqa_plus
spin_integrated = False
explicit_spin_cases = False
use_legacy_order = True

amplitude_string = 't1_p1p'

import time
start = time.time()

print("\n----------------------------------------------------------------------------------")
print("Spin-Adapted {:} correlation energy".format(amplitude_string).center(82))
print("----------------------------------------------------------------------------------\n")

# Create spin-orbital V operator
print("# Create spin-orbital V operator ...")
terms_V = sqa_plus.getV(spin_integrated, explicit_spin_cases)

# Create spin-orbital amplitude operator
## Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

dummy = True

i = sqa_plus.index('i', [tg_cor], dummy)
j = sqa_plus.index('j', [tg_cor], dummy)

x = sqa_plus.index('x', [tg_act], dummy)
y = sqa_plus.index('y', [tg_act], dummy)
z = sqa_plus.index('z', [tg_act], dummy)

a = sqa_plus.index('a', [tg_vir], dummy)
b = sqa_plus.index('b', [tg_vir], dummy)

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

## Define terms
if amplitude_string == 't1_0':
    print("## Create T: t1_{ij}^{ab} a_a^\dag a_b^\dag a_j a_i ...\n")
    t1_symm_ppqq = [sqa_plus.symmetry((0,1,3,2), -1), sqa_plus.symmetry((1,0,3,2), -1)]

    t1 = sqa_plus.tensor("t1", [i, j, a, b], t1_symm_ppqq)

    terms_t = [sqa_plus.term(0.25, [], [t1, cre_a, cre_b, des_j, des_i])]

    e_string = 'e_0'

elif amplitude_string == 't1_0p':
    print("## Create T: t1_{ix}^{ay} a_a^\dag a_y^\dag a_x a_i ...\n")
    t1_symm = [sqa_plus.symmetry((1,0), 1)]
    t1_symm_pqrs = []

    t1 = sqa_plus.tensor("t1", [i, a], t1_symm)
    t2 = sqa_plus.tensor("t1", [i, x, a, y], t1_symm_pqrs)

    op_t1 = [t1, cre_a, des_i]
    op_t2 = [t2, cre_a, cre_y, des_x, des_i]

    terms_t = [sqa_plus.term(1.00, [], op_t1),
               sqa_plus.term(1.00, [], op_t2)]

    e_string = 'e_0p'

elif amplitude_string == 't1_m1':
    print("## Create T: t1_{ix}^{ab} a_a^\dag a_b^\dag a_x a_i ...\n")
    t1_symm_pqrr = [sqa_plus.symmetry((0,1,3,2), -1)]

    t1 = sqa_plus.tensor("t1", [i, x, a, b], t1_symm_pqrr)

    terms_t = [sqa_plus.term(0.5, [], [t1, cre_a, cre_b, des_x, des_i])]

    e_string = 'e_m1'

elif amplitude_string == 't1_p1':
    print("## Create T: t1_{ij}^{ax} a_a^\dag a_x^\dag a_j a_i ...\n")
    t1_symm_ppqr = [sqa_plus.symmetry((1,0,2,3), -1)]

    t1 = sqa_plus.tensor("t1", [i, j, a, x], t1_symm_ppqr)

    terms_t = [sqa_plus.term(0.5, [], [t1, cre_a, cre_x, des_j, des_i])]

    e_string = 'e_p1'

elif amplitude_string == 't1_m2':
    print("## Create T: t1_{xy}^{ab} a_a^\dag a_b^\dag a_y a_x ...\n")
    t1_symm_ppqq = [sqa_plus.symmetry((0,1,3,2), -1), sqa_plus.symmetry((1,0,2,3), -1)]

    t1 = sqa_plus.tensor("t1", [x, y, a, b], t1_symm_ppqq)

    terms_t = [sqa_plus.term(0.25, [], [t1, cre_a, cre_b, des_y, des_x])]

    e_string = 'e_m2'

elif amplitude_string == 't1_p2':
    print("## Create T: t1_{ij}^{xy} a_x^\dag a_y^\dag a_j a_i ...\n")
    t1_symm_ppqq = [sqa_plus.symmetry((0,1,3,2), -1), sqa_plus.symmetry((1,0,2,3), -1)]

    t1 = sqa_plus.tensor("t1", [i, j, x, y], t1_symm_ppqq)

    terms_t = [sqa_plus.term(0.25, [], [t1, cre_x, cre_y, des_j, des_i])]

    e_string = 'e_p2'

elif amplitude_string == 't1_m1p':
    print("## Create T: t1_{x}^{a} a_a^\dag a_x ...\n")
    print("## Create T: t1_{xy}^{az} a_a^\dag a_z^\dag a_y a_x ...\n")
    t1_symm = [sqa_plus.symmetry((1,0), 1)]
    t1_symm_ppqr = [sqa_plus.symmetry((1,0,2,3), -1)]

    t1 = sqa_plus.tensor("t1", [x, a])
    t2 = sqa_plus.tensor("t1", [x, y, a, z], t1_symm_ppqr)

    op_t1 = [t1, cre_a, des_x]
    op_t2 = [t2, cre_a, cre_z, des_y, des_x]

    terms_t = [sqa_plus.term( 1.00, [], op_t1),
               sqa_plus.term( 0.50, [], op_t2)]

    e_string = 'e_m1p'

elif amplitude_string == 't1_p1p':
    print("## Create T: t1_{i}^{x} a_x^\dag a_i ...\n")
    print("## Create T: t1_{iz}^{xy} a_x^\dag a_y^\dag a_z a_i ...\n")
    t1_symm = [sqa_plus.symmetry((1,0), 1)]
    t1_symm_pqrr = [sqa_plus.symmetry((0,1,3,2), -1)]

    t1 = sqa_plus.tensor("t1", [i, x], t1_symm)
    t2 = sqa_plus.tensor("t1", [i, z, x, y], t1_symm_pqrr)

    op_t1 = [t1, cre_x, des_i]
    op_t2 = [t2, cre_x, cre_y, des_z, des_i]

    terms_t = [sqa_plus.term( 1.00, [], op_t1),
               sqa_plus.term( 0.50, [], op_t2)]

    e_string = 'e_p1p'

for term_t in terms_t:
    print(term_t)

# Multiply V * T
print("\n## Multiply V * T ...")
terms_e = []
for term_t in terms_t:
    for term_V in terms_V:
        terms_e.append(sqa_plus.multiplyTerms(term_V, term_t))

# Compute expected value of spin-orbital pc-NEVPT2 energy
print("## Compute expected value of spin-orbital pc-NEVPT2 energy ...")
terms_expected_e = sqa_plus.matrixBlock(terms_e, transRDM = False, legacy_order = use_legacy_order)

# Create Numpy einsum equations
terms_expected_e.sort()
result = sqa_plus.genEinsum(terms_expected_e, e_string, '', rm_core_int = True, suffix = '')

print("\n-------------------------------- genEinsum equations --------------------------------\n")
for item in result:
    print(item)
print("\n-------------------------------------------------------------------------------------\n")

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))