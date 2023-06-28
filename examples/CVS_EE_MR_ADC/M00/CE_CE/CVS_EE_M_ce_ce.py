import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
#order_Heff = 0
#order_Heff = 1
order_Heff = 2

import time
start = time.time()

sqa_plus.options.print_header("Spin-Adapted CVS-EE: M00 H{:}".format(order_Heff))

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

# Generating Term 
print("\n## Generating Term a_I^\dag a_A ... a_B\dag a_J ...\n")

### ALPHA
#### External Indices
i_alpha = sqa_plus.index('I', [tg_alpha, tg_cvs_cor])
a_alpha = sqa_plus.index('A', [tg_alpha, tg_vir])
j_alpha = sqa_plus.index('J', [tg_alpha, tg_cvs_cor])
b_alpha = sqa_plus.index('B', [tg_alpha, tg_vir])

#### Define operators
cre_i_alpha = sqa_plus.creOp(i_alpha)
des_a_alpha = sqa_plus.desOp(a_alpha)
des_j_alpha = sqa_plus.desOp(j_alpha)
cre_b_alpha = sqa_plus.creOp(b_alpha)

#### Define left and right hand side
l_op = sqa_plus.term(1.0, [], [cre_i_alpha, des_a_alpha])
r_op = sqa_plus.term(1.0, [], [cre_b_alpha, des_j_alpha])

## BETA
### External Indices
#i_beta = sqa_plus.index('I', [tg_beta, tg_cvs_cor])
#a_beta = sqa_plus.index('A', [tg_beta, tg_vir])
#j_beta = sqa_plus.index('J', [tg_beta, tg_cvs_cor])
#b_beta = sqa_plus.index('B', [tg_beta, tg_vir])

### Define operators
#cre_i_beta = sqa_plus.creOp(i_beta)
#des_a_beta = sqa_plus.desOp(a_beta)
#des_j_beta = sqa_plus.desOp(j_beta)
#cre_b_beta = sqa_plus.creOp(b_beta)

### Define left and right hand side
#l_op = sqa_plus.term(1.0, [], [cre_i_beta, des_a_beta])
#r_op = sqa_plus.term(1.0, [], [cre_b_beta, des_j_beta])

# Spin-Adapted H_eff
terms_Heff = sqa_plus.Heff(order_Heff)

## Calculate the commutators
print("## Calculate the commutator ... [H(2), r_op] ...")
terms_commutator = sqa_plus.commutator(terms_Heff, r_op)

print("\n## Calculate the commutator [l_op, [H(2), r_op]] ...")
terms_EE_M0_ce_ce = sqa_plus.commutator(l_op, terms_commutator)

# Expected value of Spin-Integrated EE M2 CE-CE
expected_EE_M0_ce_ce = sqa_plus.matrixBlock(terms_EE_M0_ce_ce)

# Spin-Adaptation of EE M2 CE-CE
expected_EE_M0_ce_ce_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_M0_ce_ce)

# Generate Numpy einsum equations
result = sqa_plus.genEinsum(expected_EE_M0_ce_ce_sa, 'temp', 'IAJB')

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
