import time
import sqa_plus
from sqa_plus import options

options.cvs_approach = True
options.spin_integrated = True

#order_Heff = 0
order_Heff = 1

indices_string = 'c_cva'
spin_indices_string = 'a_aaa'

start = time.time()
options.print_header("Spin-Adapted CVS-IP: M01 H{:}".format(order_Heff))

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
b_alpha = sqa_plus.index('B', [tg_a, tg_vir])
b_beta  = sqa_plus.index('B', [tg_b, tg_vir])

i_alpha = sqa_plus.index('I', [tg_a, tg_cvs_cor])
i_beta  = sqa_plus.index('I', [tg_b, tg_cvs_cor])

k_alpha = sqa_plus.index('K', [tg_a, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_b, tg_cvs_cor])

l_alpha = sqa_plus.index('L', [tg_a, tg_cvs_cor])
l_beta  = sqa_plus.index('L', [tg_b, tg_cvs_cor])

l_val_alpha = sqa_plus.index('L', [tg_a, tg_cvs_val])
l_val_beta  = sqa_plus.index('L', [tg_b, tg_cvs_val])

w_alpha = sqa_plus.index('W', [tg_a, tg_act])
w_beta  = sqa_plus.index('W', [tg_b, tg_act])

y_alpha = sqa_plus.index('Y', [tg_a, tg_act])
y_beta  = sqa_plus.index('Y', [tg_b, tg_act])

z_alpha = sqa_plus.index('Z', [tg_a, tg_act])
z_beta  = sqa_plus.index('Z', [tg_b, tg_act])

i     = {'a': i_alpha,     'b': i_beta}
k     = {'a': k_alpha,     'b': k_beta}
l     = {'a': l_alpha,     'b': l_beta}
l_val = {'a': l_val_alpha, 'b': l_val_beta}
w     = {'a': w_alpha,     'b': w_beta}
y     = {'a': y_alpha,     'b': y_beta}
z     = {'a': z_alpha,     'b': z_beta}
b     = {'a': b_alpha,     'b': b_beta}

## Define terms
ops = {
    'c_caa': (i, k, w,     z, 'IKWZ'),
    'c_cce': (i, k, l,     b, 'IKLB'),
    'c_cve': (i, k, l_val, b, 'IKLB'),
    'c_cca': (i, k, l,     y, 'IKLY'),
    'c_cva': (i, k, l_val, y, 'IKLY'),
    'c_cae': (i, k, y,     b, 'IKYB'),
}

spins = {
    'a_aaa': ('a', 'a', 'a', 'a'),
    'a_abb': ('a', 'a', 'b', 'b'),
    'a_bab': ('a', 'b', 'a', 'b'),
}

l_ind, *r_ind = ops[indices_string]
s1, s2, s3, s4 = spins[spin_indices_string]

term_left  = sqa_plus.term(1.0, [], [sqa_plus.desOp(l_ind[s1])])
term_right = sqa_plus.term(1.0, [], [sqa_plus.creOp(r_ind[0][s2]), sqa_plus.creOp(r_ind[1][s3]), sqa_plus.desOp(r_ind[2][s4])])
final_indices_string = r_ind[3]
indices_string = indices_string + '_' + spin_indices_string

print("##  Left Op: %s" % term_left)
print("## Right Op: %s" % term_right)

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculating the commutator
print("## Calculating the commutator [H(0), a_S^\\dag a_T^\\dag a_U] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, term_right)

print("\n## Calculating a_Q [H(0), a_S^\\dag a_T^\\dag a_U] ...")
terms_IP_M01 = []
for term_comm in terms_commutator:
    terms_IP_M01.append(sqa_plus.multiplyTerms(term_comm, term_left))
    terms_IP_M01.append(sqa_plus.multiplyTerms(term_left, term_comm))

# Expected value of Spin-Adapted IP M01
expected_IP_M01 = sqa_plus.matrixBlock(terms_IP_M01)

# Spin-Adaptation of IP M01
expected_IP_M01_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_M01)

# Generating Numpy einsum equations
options.genEinsum.lhs_string = 'M_' + indices_string
options.genEinsum.indices_string = final_indices_string
result = sqa_plus.genEinsum(expected_IP_M01_sa)

print("> Total elapsed time: {:.2f} seconds.".format(time.time() - start))

