import time
import sqa_plus
from sqa_plus import options

options.cvs_approach = True
options.spin_integrated = True

diagonal_indices_string = 'cva'
#spin_indices_string = 'aaa'
#spin_indices_string = 'abb'
spin_indices_string = 'bab'

start = time.time()
options.print_header("Spin-Adapted Preconditioner {:} {:}".format(diagonal_indices_string, spin_indices_string))

# Generating operators
print("\n## Generating operators ...\n")

## Define indices
tg_cvs_cor = options.cvs_core_type
tg_cvs_val = options.cvs_valence_type
tg_act = options.active_type
tg_vir = options.virtual_type

tg_a = options.alpha_type
tg_b = options.beta_type

## External Indices
a_alpha = sqa_plus.index('A', [tg_a, tg_vir])
a_beta  = sqa_plus.index('A', [tg_b, tg_vir])

b_alpha = sqa_plus.index('B', [tg_a, tg_vir])
b_beta  = sqa_plus.index('B', [tg_b, tg_vir])

i_alpha = sqa_plus.index('I', [tg_a, tg_cvs_cor])
i_beta  = sqa_plus.index('I', [tg_b, tg_cvs_cor])

j_alpha = sqa_plus.index('J', [tg_a, tg_cvs_cor])
j_beta  = sqa_plus.index('J', [tg_b, tg_cvs_cor])

j_val_alpha = sqa_plus.index('J', [tg_a, tg_cvs_val])
j_val_beta  = sqa_plus.index('J', [tg_b, tg_cvs_val])

k_alpha = sqa_plus.index('K', [tg_a, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_b, tg_cvs_cor])

l_alpha = sqa_plus.index('L', [tg_a, tg_cvs_cor])
l_beta  = sqa_plus.index('L', [tg_b, tg_cvs_cor])

l_val_alpha = sqa_plus.index('L', [tg_a, tg_cvs_val])
l_val_beta  = sqa_plus.index('L', [tg_b, tg_cvs_val])

x_alpha = sqa_plus.index('X', [tg_a, tg_act])
x_beta  = sqa_plus.index('X', [tg_b, tg_act])

y_alpha = sqa_plus.index('Y', [tg_a, tg_act])
y_beta  = sqa_plus.index('Y', [tg_b, tg_act])

w_alpha = sqa_plus.index('W', [tg_a, tg_act])
w_beta  = sqa_plus.index('W', [tg_b, tg_act])

z_alpha = sqa_plus.index('Z', [tg_a, tg_act])
z_beta  = sqa_plus.index('Z', [tg_b, tg_act])

a     = {'a': a_alpha,     'b': a_beta}
b     = {'a': b_alpha,     'b': b_beta}
i     = {'a': i_alpha,     'b': i_beta}
j     = {'a': j_alpha,     'b': j_beta}
j_val = {'a': j_val_alpha, 'b': j_val_beta}
k     = {'a': k_alpha,     'b': k_beta}
l     = {'a': l_alpha,     'b': l_beta}
l_val = {'a': l_val_alpha, 'b': l_val_beta}
w     = {'a': w_alpha,     'b': w_beta}
x     = {'a': x_alpha,     'b': x_beta}
y     = {'a': y_alpha,     'b': y_beta}
z     = {'a': z_alpha,     'b': z_beta}

## Define terms
spins = {
    'aaa':     ('a', 'a', 'a'),
    'bab':     ('b', 'a', 'b'),
    'abb':     ('a', 'b', 'b'),
}

if diagonal_indices_string == 'c_caa':
    ls1 = spin_indices_string[0]
    rs1, rs2, rs3 = spins[spin_indices_string[-3:]]
    term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(k[ls1])])
    term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i[rs1]), sqa_plus.creOp(x[rs2]), sqa_plus.desOp(y[rs3])])
    final_indices_string = 'IXY'
    diagonal_pairs_dict  = {'K': 'I'}

elif diagonal_indices_string == 'caa':
    ls1, ls2, ls3 = spins[spin_indices_string[:3]]
    rs1, rs2, rs3 = spins[spin_indices_string[-3:]]
    term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(z[ls3]), sqa_plus.desOp(w[ls2]), sqa_plus.desOp(k[ls1])])
    term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(i[rs1]), sqa_plus.creOp(x[rs2]), sqa_plus.desOp(y[rs3])])
    final_indices_string = 'IWZXY'
    diagonal_pairs_dict  = {'K': 'I'}

else:
    ops = {
        'cce': (k, l,     b, i, j,     a, 'IJA',  {'L': 'J', 'K': 'I', 'B': 'A'}),
        'cve': (k, l_val, b, i, j_val, a, 'IJA',  {'L': 'J', 'K': 'I', 'B': 'A'}),
        'cca': (k, l,     y, i, j,     x, 'IJXY', {'L': 'J', 'K': 'I'}),
        'cva': (k, l_val, y, i, j_val, x, 'IJXY', {'L': 'J', 'K': 'I'}),
        'cae': (k, w,     b, i, x,     a, 'IAXY', {'K': 'I', 'B': 'A'}),
 
    }
    s1, s2, s3 = spins[spin_indices_string]

    l1, l2, l3, r1, r2, r3, final_indices_string, diagonal_pairs_dict = ops[diagonal_indices_string]

    term_left  = sqa_plus.term(1.0, [], [sqa_plus.creOp(l3[s3]), sqa_plus.desOp(l2[s2]), sqa_plus.desOp(l1[s1])])
    term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(r1[s1]), sqa_plus.creOp(r2[s2]), sqa_plus.desOp(r3[s3])])

diagonal_indices_string = diagonal_indices_string + '__' + spin_indices_string

print("##  Left Op: %s" % term_left)
print("## Right Op: %s" % term_right)

# Spin-Integrated H_eff
if 'c_caa' in diagonal_indices_string:
    terms_Heff = sqa_plus.Heff(0)
    terms_Heff.extend(sqa_plus.Heff(1))
else:
    terms_Heff = sqa_plus.Heff(0)

## Calculating the commutator
print("## Calculating the commutator [H(0), a_S^\\dag a_T^\\dag a_U] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating a_P^\\dag a_Q a_R [H(0), a_S^\\dag a_T^\\dag a_U] ...")
terms_IP_M11 = []
for term_comm in terms_commutator:
    terms_IP_M11.append(sqa_plus.multiplyTerms(term_comm, term_left))
    terms_IP_M11.append(sqa_plus.multiplyTerms(term_left, term_comm))

# Expected value of Spin-Adapted IP M11
expected_IP_M11 = sqa_plus.matrixBlock(terms_IP_M11)

# Spin-Adaptation of IP M11
expected_IP_M11_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_M11)

from sqa_plus.sqaEinsum import remove_core_int
expected_IP_M11_sa, removed_core = remove_core_int(expected_IP_M11_sa)

# Replacing indices to obtain the diagonal
for term_sa in expected_IP_M11_sa:
    for tensor_sa in term_sa.tensors:
        for index_sa in tensor_sa.indices:
            if index_sa.name in diagonal_pairs_dict:
                index_sa.userDefined = True
                index_sa.name = diagonal_pairs_dict[index_sa.name]

from sqa_plus.sqaMatrixBlock import contractDeltaFuncs_nondummy
expected_IP_M11_sa = contractDeltaFuncs_nondummy(expected_IP_M11_sa)

sqa_plus.combineTerms(expected_IP_M11_sa)

# Generating Numpy einsum equations
options.genEinsum.remove_core_integrals = False
options.genEinsum.keep_user_defined_dummy_names = True
options.genEinsum.lhs_string = 'precond_' + diagonal_indices_string
options.genEinsum.indices_string = final_indices_string
result = sqa_plus.genEinsum(expected_IP_M11_sa)
print("> Total elapsed time: {:.2f} seconds.".format(time.time() - start))
