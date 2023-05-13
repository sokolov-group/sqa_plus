import sqa_plus
sqa_plus.options.spin_integrated = True

amplitude_string = 't1_m1p'

import time
start = time.time()

sqa_plus.options.print_header("Spin-Adapted {:} correlation energy".format(amplitude_string))

# Create spin-integrated amplitude operator
## Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

dummy = True

i_alpha = sqa_plus.index('i', [tg_alpha, tg_cor], dummy)
i_beta  = sqa_plus.index('i', [tg_beta,  tg_cor], dummy)

j_alpha = sqa_plus.index('j', [tg_alpha, tg_cor], dummy)
j_beta  = sqa_plus.index('j', [tg_beta,  tg_cor], dummy)

x_alpha = sqa_plus.index('x', [tg_alpha, tg_act], dummy)
x_beta  = sqa_plus.index('x', [tg_beta,  tg_act], dummy)

y_alpha = sqa_plus.index('y', [tg_alpha, tg_act], dummy)
y_beta  = sqa_plus.index('y', [tg_beta,  tg_act], dummy)

z_alpha = sqa_plus.index('z', [tg_alpha, tg_act], dummy)
z_beta  = sqa_plus.index('z', [tg_beta,  tg_act], dummy)

a_alpha = sqa_plus.index('a', [tg_alpha, tg_vir], dummy)
a_beta  = sqa_plus.index('a', [tg_beta,  tg_vir], dummy)

b_alpha = sqa_plus.index('b', [tg_alpha, tg_vir], dummy)
b_beta  = sqa_plus.index('b', [tg_beta,  tg_vir], dummy)

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

t1_symm = [sqa_plus.symmetry((1,0), 1)]
t2_symm = [sqa_plus.symmetry((1,0,2,3), -1), sqa_plus.symmetry((0,1,3,2), -1)]

## Define terms
if amplitude_string == 't1_0':
    print("## Create T: t1_{ij}^{ab} a_a^\dag a_b^\dag a_j a_i ...\n")

    t1_aaaa = sqa_plus.tensor("t1", [i_alpha, j_alpha, a_alpha, b_alpha], t2_symm)
    t1_abab = sqa_plus.tensor("t1", [i_alpha, j_beta,  a_alpha, b_beta],  t2_symm)
    t1_bbbb = sqa_plus.tensor("t1", [i_beta,  j_beta,  a_beta,  b_beta],  t2_symm)

    terms_t = [sqa_plus.term(0.25, [], [t1_aaaa, cre_a_alpha, cre_b_alpha, des_j_alpha, des_i_alpha]),
               sqa_plus.term(1.00, [], [t1_abab, cre_a_alpha, cre_b_beta,  des_j_beta,  des_i_alpha]),
               sqa_plus.term(0.25, [], [t1_bbbb, cre_a_beta,  cre_b_beta,  des_j_beta,  des_i_beta])]

    e_string = 'e_0'

elif amplitude_string == 't1_0p':
    print("## Create T: t1_{ix}^{ay} a_a^\dag a_y^\dag a_x a_i ...\n")

    t1_aa = sqa_plus.tensor("t1", [i_alpha, a_alpha], t1_symm)
    t1_bb = sqa_plus.tensor("t1", [i_beta,  a_beta],  t1_symm)

    t1_aaaa = sqa_plus.tensor("t1", [i_alpha, x_alpha, a_alpha, y_alpha], t2_symm)
    t1_abab = sqa_plus.tensor("t1", [i_alpha, x_beta,  a_alpha, y_beta],  t2_symm)
    t1_abba = sqa_plus.tensor("t1", [i_alpha, x_beta,  a_beta,  y_alpha], t2_symm)
    t1_baba = sqa_plus.tensor("t1", [i_beta,  x_alpha, a_beta,  y_alpha], t2_symm)
    t1_baab = sqa_plus.tensor("t1", [i_beta,  x_alpha, a_alpha, y_beta],  t2_symm)
    t1_bbbb = sqa_plus.tensor("t1", [i_beta,  x_beta,  a_beta,  y_beta],  t2_symm)

    op_t1_aa = [t1_aa, cre_a_alpha, des_i_alpha]
    op_t1_bb = [t1_bb, cre_a_beta,  des_i_beta]

    op_t1_aaaa = [t1_aaaa, cre_a_alpha, cre_y_alpha, des_x_alpha, des_i_alpha]
    op_t1_abab = [t1_abab, cre_a_alpha, cre_y_beta,  des_x_beta,  des_i_alpha]
    op_t1_abba = [t1_abba, cre_a_beta,  cre_y_alpha, des_x_beta,  des_i_alpha]
    op_t1_baba = [t1_baba, cre_a_beta,  cre_y_alpha, des_x_alpha, des_i_beta]
    op_t1_baab = [t1_baab, cre_a_alpha, cre_y_beta,  des_x_alpha, des_i_beta]
    op_t1_bbbb = [t1_bbbb, cre_a_beta,  cre_y_beta,  des_x_beta,  des_i_beta]

    terms_t = [sqa_plus.term(1.00, [], op_t1_aa),
               sqa_plus.term(1.00, [], op_t1_bb),
               sqa_plus.term(1.00, [], op_t1_aaaa),
               sqa_plus.term(1.00, [], op_t1_abab),
               sqa_plus.term(1.00, [], op_t1_abba),
               sqa_plus.term(1.00, [], op_t1_baba),
               sqa_plus.term(1.00, [], op_t1_baab),
               sqa_plus.term(1.00, [], op_t1_bbbb)]

    e_string = 'e_0p'

elif amplitude_string == 't1_m1':
    print("## Create T: t1_{ix}^{ab} a_a^\dag a_b^\dag a_x a_i ...\n")

    t1_aaaa = sqa_plus.tensor("t1", [i_alpha, x_alpha, a_alpha, b_alpha], t2_symm)
    t1_abab = sqa_plus.tensor("t1", [i_alpha, x_beta,  a_alpha, b_beta],  t2_symm)
    t1_baba = sqa_plus.tensor("t1", [i_beta,  x_alpha, a_beta,  b_alpha], t2_symm)
    t1_bbbb = sqa_plus.tensor("t1", [i_beta,  x_beta,  a_beta,  b_beta],  t2_symm)

    terms_t = [sqa_plus.term(0.5, [], [t1_aaaa, cre_a_alpha, cre_b_alpha, des_x_alpha, des_i_alpha]),
               sqa_plus.term(1.0, [], [t1_abab, cre_a_alpha, cre_b_beta,  des_x_beta,  des_i_alpha]),
               sqa_plus.term(1.0, [], [t1_baba, cre_a_beta,  cre_b_alpha, des_x_alpha, des_i_beta]),
               sqa_plus.term(0.5, [], [t1_bbbb, cre_a_beta,  cre_b_beta,  des_x_beta,  des_i_beta])]

    e_string = 'e_m1'

elif amplitude_string == 't1_p1':
    print("## Create T: t1_{ij}^{ax} a_a^\dag a_x^\dag a_j a_i ...\n")

    t1_aaaa = sqa_plus.tensor("t1", [i_alpha, j_alpha, a_alpha, x_alpha], t2_symm)
    t1_abab = sqa_plus.tensor("t1", [i_alpha, j_beta,  a_alpha, x_beta],  t2_symm)
    t1_abba = sqa_plus.tensor("t1", [i_alpha, j_beta,  a_beta,  x_alpha], t2_symm)
    t1_bbbb = sqa_plus.tensor("t1", [i_beta,  j_beta,  a_beta,  x_beta],  t2_symm)

    terms_t = [sqa_plus.term(0.5, [], [t1_aaaa, cre_a_alpha, cre_x_alpha, des_j_alpha, des_i_alpha]),
               sqa_plus.term(1.0, [], [t1_abab, cre_a_alpha, cre_x_beta,  des_j_beta,  des_i_alpha]),
               sqa_plus.term(1.0, [], [t1_abba, cre_a_beta,  cre_x_alpha, des_j_beta,  des_i_alpha]),
               sqa_plus.term(0.5, [], [t1_bbbb, cre_a_beta,  cre_x_beta,  des_j_beta,  des_i_beta])]

    e_string = 'e_p1'

elif amplitude_string == 't1_m2':
    print("## Create T: t1_{xy}^{ab} a_a^\dag a_b^\dag a_y a_x ...\n")

    t1_aaaa = sqa_plus.tensor("t1", [x_alpha, y_alpha, a_alpha, b_alpha], t2_symm)
    t1_abab = sqa_plus.tensor("t1", [x_alpha, y_beta,  a_alpha, b_beta],  t2_symm)
    t1_baba = sqa_plus.tensor("t1", [x_beta,  y_alpha, a_beta,  b_alpha], t2_symm)
    t1_bbbb = sqa_plus.tensor("t1", [x_beta,  y_beta,  a_beta,  b_beta],  t2_symm)

    terms_t = [sqa_plus.term(0.25, [], [t1_aaaa, cre_a_alpha, cre_b_alpha, des_y_alpha, des_x_alpha]),
               sqa_plus.term(0.50, [], [t1_abab, cre_a_alpha, cre_b_beta,  des_y_beta,  des_x_alpha]),
               sqa_plus.term(0.50, [], [t1_baba, cre_a_beta,  cre_b_alpha, des_y_alpha, des_x_beta]),
               sqa_plus.term(0.25, [], [t1_bbbb, cre_a_beta,  cre_b_beta,  des_y_beta,  des_x_beta])]

    e_string = 'e_m2'

elif amplitude_string == 't1_p2':
    print("## Create T: t1_{ij}^{xy} a_x^\dag a_y^\dag a_j a_i ...\n")

    t1_aaaa = sqa_plus.tensor("t1", [i_alpha, j_alpha, x_alpha, y_alpha], t2_symm)
    t1_abab = sqa_plus.tensor("t1", [i_alpha, j_beta,  x_alpha, y_beta],  t2_symm)
    t1_bbbb = sqa_plus.tensor("t1", [i_beta,  j_beta,  x_beta,  y_beta],  t2_symm)

    terms_t = [sqa_plus.term(0.25, [], [t1_aaaa, cre_x_alpha, cre_y_alpha, des_j_alpha, des_i_alpha]),
               sqa_plus.term(1.00, [], [t1_abab, cre_x_alpha, cre_y_beta,  des_j_beta,  des_i_alpha]),
               sqa_plus.term(0.25, [], [t1_bbbb, cre_x_beta,  cre_y_beta,  des_j_beta,  des_i_beta])]

    e_string = 'e_p2'

elif amplitude_string == 't1_m1p':
    print("## Create T: t1_{x}^{a} a_a^\dag a_x ...\n")
    print("## Create T: t1_{xy}^{az} a_a^\dag a_z^\dag a_y a_x ...\n")

    t1_aa = sqa_plus.tensor("t1", [x_alpha, a_alpha])
    t1_bb = sqa_plus.tensor("t1", [x_beta,  a_beta])

    t2_aaaa = sqa_plus.tensor("t1", [x_alpha, y_alpha, a_alpha, z_alpha], t2_symm)
    t2_abab = sqa_plus.tensor("t1", [x_alpha, y_beta,  a_alpha, z_beta],  t2_symm)
    t2_abba = sqa_plus.tensor("t1", [x_alpha, y_beta,  a_beta,  z_alpha], t2_symm)
    t2_bbbb = sqa_plus.tensor("t1", [x_beta,  y_beta,  a_beta,  z_beta],  t2_symm)

    op_t1_aa = [t1_aa, cre_a_alpha, des_x_alpha]
    op_t1_bb = [t1_bb, cre_a_beta,  des_x_beta]

    op_t2_aaaa = [t2_aaaa, cre_a_alpha, cre_z_alpha, des_y_alpha, des_x_alpha]
    op_t2_abab = [t2_abab, cre_a_alpha, cre_z_beta,  des_y_beta,  des_x_alpha]
    op_t2_abba = [t2_abba, cre_a_beta,  cre_z_alpha, des_y_beta,  des_x_alpha]
    op_t2_bbbb = [t2_bbbb, cre_a_beta,  cre_z_beta,  des_y_beta,  des_x_beta]

    terms_t = [sqa_plus.term( 1.00, [], op_t1_aa),
               sqa_plus.term( 1.00, [], op_t1_bb),
               sqa_plus.term( 0.50, [], op_t2_aaaa),
               sqa_plus.term( 1.00, [], op_t2_abab),
               sqa_plus.term( 1.00, [], op_t2_abba),
               sqa_plus.term( 0.50, [], op_t2_bbbb)]

    e_string = 'e_m1p'

elif amplitude_string == 't1_p1p':
    print("## Create T: t1_{i}^{x} a_x^\dag a_i ...\n")
    print("## Create T: t1_{iz}^{xy} a_x^\dag a_y^\dag a_z a_i ...\n")

    t1_aa = sqa_plus.tensor("t1", [i_alpha, x_alpha], t1_symm)
    t1_bb = sqa_plus.tensor("t1", [i_beta,  x_beta],  t1_symm)

    t1_aaaa = sqa_plus.tensor("t1", [i_alpha, z_alpha, x_alpha, y_alpha], t2_symm)
    t1_abab = sqa_plus.tensor("t1", [i_alpha, z_beta,  x_alpha, y_beta],  t2_symm)
    t1_baab = sqa_plus.tensor("t1", [i_beta,  z_alpha, x_alpha, y_beta],  t2_symm)
    t1_bbbb = sqa_plus.tensor("t1", [i_beta,  z_beta,  x_beta,  y_beta],  t2_symm)

    op_t1_aa = [t1_aa, cre_x_alpha, des_i_alpha]
    op_t1_bb = [t1_bb, cre_x_beta,  des_i_beta]

    op_t1_aaaa = [t1_aaaa, cre_x_alpha, cre_y_alpha, des_z_alpha, des_i_alpha]
    op_t1_abab = [t1_abab, cre_x_alpha, cre_y_beta,  des_z_beta,  des_i_alpha]
    op_t1_baab = [t1_baab, cre_x_alpha, cre_y_beta,  des_z_alpha, des_i_beta]
    op_t1_bbbb = [t1_bbbb, cre_x_beta,  cre_y_beta,  des_z_beta,  des_i_beta]

    terms_t = [sqa_plus.term( 1.00, [], op_t1_aa),
               sqa_plus.term( 1.00, [], op_t1_bb),
               sqa_plus.term( 0.50, [], op_t1_aaaa),
               sqa_plus.term( 1.00, [], op_t1_abab),
               sqa_plus.term( 1.00, [], op_t1_baab),
               sqa_plus.term( 0.50, [], op_t1_bbbb)]

    e_string = 'e_p1p'

for term_t in terms_t:
    print(term_t)

# Multiply V * T
print("\n## Multiply V * T ...")
terms_e = []
for term_t in terms_t:
    for term_V in terms_V:
        terms_e.append(sqa_plus.multiplyTerms(term_V, term_t))

# Compute expected value of spin-integrated pc-NEVPT2 energy
print("## Compute expected value of spin-integrated pc-NEVPT2 energy ...")
terms_expected_e = sqa_plus.matrixBlock(terms_e)

# Spin-Adaptating pc-NEVPT2 energy expression
terms_expected_e_sa = sqa_plus.convertSpinIntegratedToAdapted(terms_expected_e)

# Create Numpy einsum equations
result = sqa_plus.genEinsum(terms_expected_e_sa, e_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))