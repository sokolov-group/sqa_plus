import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 1

import time
start = time.time()

sqa_plus.options.print_header("Spin-Adapted CVS-IP: M00 H{:}".format(order_Heff))

# Generating Term a_I^\dag a_J
print("\n## Generating Term a_I^\dag a_J ...\n")

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices
i_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_cor])
j_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_cor])

## Define operators
des_i_alpha = sqa_plus.desOp(i_alpha)
cre_j_alpha = sqa_plus.creOp(j_alpha)

term_des_i = sqa_plus.term(1.0, [], [des_i_alpha])
term_cre_j = sqa_plus.term(1.0, [], [cre_j_alpha])

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculate the commutator
print("## Calculate the commutator [H(2), a_J^\dag] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, term_cre_j)

print("\n## Calculate a_I [H(2), a_J^\dag] ...")
terms_IP_M0_c_c = []
for term_commutator in terms_commutator:
    terms_IP_M0_c_c.append(sqa_plus.multiplyTerms(term_commutator, term_des_i))

# Expected value of Spin-Integrated IP M2 C-C
expected_IP_M0_c_c = sqa_plus.matrixBlock(terms_IP_M0_c_c)

# Spin-Adaptation of IP M2 C-C
expected_IP_M0_c_c_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_M0_c_c)

# Generate Numpy einsum equations
result = sqa_plus.genEinsum(expected_IP_M0_c_c_sa, 'M00', 'IJ')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))

