import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True

import time
start = time.time()

diagonal_indices_string = 'c_caa'
spin_indices_string = 'b_bab'

sqa_plus.options.print_header("Spin-Adapted Preconditioner {:} {:}".format(diagonal_indices_string, spin_indices_string))

# Generating operators
print("\n## Generating operators ...\n")

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices
a_alpha = sqa_plus.index('A', [tg_alpha, tg_vir])
a_beta  = sqa_plus.index('A', [tg_beta,  tg_vir])

b_alpha = sqa_plus.index('B', [tg_alpha, tg_vir])
b_beta  = sqa_plus.index('B', [tg_beta,  tg_vir])

i_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_cor])
i_beta  = sqa_plus.index('I', [tg_beta,  tg_cvs_cor])

i_val_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_val])
i_val_beta  = sqa_plus.index('I', [tg_beta,  tg_cvs_val])

j_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_cor])
j_beta  = sqa_plus.index('J', [tg_beta,  tg_cvs_cor])

j_val_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_val])
j_val_beta  = sqa_plus.index('J', [tg_beta,  tg_cvs_val])

k_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_cor])

k_val_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_val])
k_val_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_val])

l_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_cor])
l_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_cor])

l_val_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_val])
l_val_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_val])

x_alpha = sqa_plus.index('X', [tg_alpha, tg_act])
x_beta  = sqa_plus.index('X', [tg_beta,  tg_act])

y_alpha = sqa_plus.index('Y', [tg_alpha, tg_act])
y_beta  = sqa_plus.index('Y', [tg_beta,  tg_act])

w_alpha = sqa_plus.index('W', [tg_alpha, tg_act])
w_beta  = sqa_plus.index('W', [tg_beta,  tg_act])

z_alpha = sqa_plus.index('Z', [tg_alpha, tg_act])
z_beta  = sqa_plus.index('Z', [tg_beta,  tg_act])

## Define terms
if diagonal_indices_string in ['cce']:
    if spin_indices_string == 'aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(a_alpha)])

    elif spin_indices_string == 'bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(a_beta)])

    elif spin_indices_string == 'abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta),  sqa_plus.desOp(l_beta), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(a_beta)])

    final_indices_string = 'IJA'
    diagonal_pairs_dict = {'L': 'J', 'K': 'I', 'B': 'A'}

elif diagonal_indices_string in ['cve']:
    if spin_indices_string == 'aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_alpha), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(a_alpha)])

    elif spin_indices_string == 'bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(a_beta)])

    elif spin_indices_string == 'abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta),  sqa_plus.desOp(l_val_beta), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(a_beta)])

    final_indices_string = 'IJA'
    diagonal_pairs_dict = {'L': 'J', 'K': 'I', 'B': 'A'}

elif diagonal_indices_string in ['cca']:
    if spin_indices_string == 'aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_alpha), sqa_plus.desOp(x_alpha)])

    elif spin_indices_string == 'bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_beta), sqa_plus.desOp(l_alpha), sqa_plus.desOp(k_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_alpha), sqa_plus.desOp(x_beta)])

    elif spin_indices_string == 'abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_beta),  sqa_plus.desOp(l_beta), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_beta), sqa_plus.desOp(x_beta)])

    final_indices_string = 'IJXY'
    diagonal_pairs_dict = {'L': 'J', 'K': 'I'}

elif diagonal_indices_string in ['cva']:
    if spin_indices_string == 'aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_alpha), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(x_alpha)])

    elif spin_indices_string == 'bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_beta), sqa_plus.desOp(l_val_alpha), sqa_plus.desOp(k_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(j_val_alpha), sqa_plus.desOp(x_beta)])

    elif spin_indices_string == 'abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(y_beta),  sqa_plus.desOp(l_val_beta), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(j_val_beta), sqa_plus.desOp(x_beta)])

    final_indices_string = 'IJXY'
    diagonal_pairs_dict = {'L': 'J', 'K': 'I'}

elif diagonal_indices_string in ['cae']:
    if spin_indices_string == 'aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(a_alpha)])

    elif spin_indices_string == 'abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta),  sqa_plus.desOp(w_beta), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(a_beta)])

    elif spin_indices_string == 'bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(b_beta), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(a_beta)])

    final_indices_string = 'IAXY'
    diagonal_pairs_dict = {'K': 'I', 'B': 'A'}

elif diagonal_indices_string in ['c_caa']:
    if spin_indices_string == 'a_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'a_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'b_bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(k_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_beta)])

    final_indices_string = 'IXY'
    diagonal_pairs_dict = {'K': 'I'}

elif diagonal_indices_string in ['caa']:
    if spin_indices_string == 'aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(z_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_alpha)])

    elif spin_indices_string == 'abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(z_beta),  sqa_plus.desOp(w_beta), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'bab':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(z_beta), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_beta)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_beta), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'aaa_abb':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(z_alpha), sqa_plus.desOp(w_alpha), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_beta), sqa_plus.desOp(y_beta)])

    elif spin_indices_string == 'abb_aaa':
        term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(z_beta),  sqa_plus.desOp(w_beta), sqa_plus.desOp(k_alpha)])
        term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i_alpha), sqa_plus.creOp(x_alpha), sqa_plus.desOp(y_alpha)])

    final_indices_string = 'IWZXY'
    diagonal_pairs_dict = {'K': 'I'}

# Spin-Integrated H_eff
if diagonal_indices_string in ['c_caa']:
    terms_Heff = sqa_plus.Heff(0)
    terms_Heff.extend(sqa_plus.Heff(1))
else:
    terms_Heff = sqa_plus.Heff(0)

## Calculating the commutator
print("## Calculating the commutator [H(0), a_S^\dag a_T^\dag a_U] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating a_P^\dag a_Q a_R [H(0), a_S^\dag a_T^\dag a_U] ...")
terms_IP_M11 = []
for term_commutator in terms_commutator:
    terms_IP_M11.append(sqa_plus.multiplyTerms(term_commutator, term_left))
    terms_IP_M11.append(sqa_plus.multiplyTerms(term_left, term_commutator))

# Expected value of Spin-Adapted IP M11
expected_IP_M11 = sqa_plus.matrixBlock(terms_IP_M11)

# Spin-Adaptation of IP M11
expected_IP_M11_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_M11)

from sqa_plus.sqaEinsum import remove_core_int
expected_IP_M11_sa, removed_core = remove_core_int(expected_IP_M11_sa)
sqa_plus.options.genEinsum.remove_core_integrals = False
sqa_plus.options.genEinsum.keep_user_defined_dummy_names = True

# Replacing indices to obtain the diagonal
for term_ind, term_sa in enumerate(expected_IP_M11_sa):
    for tensor_ind, tensor_sa in enumerate(term_sa.tensors):
        for index_ind, index_sa in enumerate(tensor_sa.indices):
            if index_sa.name in diagonal_pairs_dict.keys():
                expected_IP_M11_sa[term_ind].tensors[tensor_ind].indices[index_ind].userDefined = True
                expected_IP_M11_sa[term_ind].tensors[tensor_ind].indices[index_ind].name = diagonal_pairs_dict[index_sa.name]

from sqa_plus.sqaMatrixBlock import contractDeltaFuncs_nondummy
expected_IP_M11_sa = contractDeltaFuncs_nondummy(expected_IP_M11_sa)

sqa_plus.combineTerms(expected_IP_M11_sa)

# Generating Numpy einsum equations
result = sqa_plus.genEinsum(expected_IP_M11_sa, 'precond_' + diagonal_indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
