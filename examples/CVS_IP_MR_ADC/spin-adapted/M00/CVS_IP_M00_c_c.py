import time
import sqa_plus

from sqa_plus import options
options.cvs_approach = True
options.spin_integrated = True

#order_Heff = 0
#order_Heff = 1
order_Heff = 2

start = time.time()

options.print_header("Spin-Adapted CVS-IP: M00 H{:}".format(order_Heff))

# Generating Term a_I^\dag a_J
print("\n## Generating Term a_I^\\dag a_J ...\n")

## Define indices
tg_cvs_cor = options.cvs_core_type
tg_cvs_val = options.cvs_valence_type
tg_act = options.active_type
tg_vir = options.virtual_type

tg_a = options.alpha_type
tg_b = options.beta_type

## External Indices
i_alpha = sqa_plus.index('I', [tg_a, tg_cvs_cor])
j_alpha = sqa_plus.index('J', [tg_a, tg_cvs_cor])

## Define operators
des_i_alpha = sqa_plus.desOp(i_alpha)
cre_j_alpha = sqa_plus.creOp(j_alpha)

term_des_i = sqa_plus.term(1.0, [], [des_i_alpha])
term_cre_j = sqa_plus.term(1.0, [], [cre_j_alpha])

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculate the commutator
print(f"## Calculate the commutator [H({order_Heff}), a_J^\\dag] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, term_cre_j)

print(f"\n## Calculate a_I [H({order_Heff}), a_J^\\dag] ...")
terms_IP_M00 = [sqa_plus.multiplyTerms(term_comm, term_des_i) for term_comm in terms_commutator]

# Expected value of Spin-Integrated IP M00 C-C
expected_IP_M00 = sqa_plus.matrixBlock(terms_IP_M00)

# Spin-Adaptation of IP M00 C-C
expected_IP_M00_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_IP_M00)

# Generate Numpy einsum equations
options.genEinsum.lhs_string = 'M00'
options.genEinsum.indices_string = 'IJ'
result = sqa_plus.genEinsum(expected_IP_M00_sa)

print("> Total elapsed time: {:.2f} seconds.".format(time.time() - start))

