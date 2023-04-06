import sqa_plus
use_legacy_order = True

import time
start = time.time()

print("\n----------------------------------------------------------------------------------")
print("Spin-Orbital {:} t2_0p V".center(82))
print("----------------------------------------------------------------------------------\n")

# Create spin-orbital amplitude operator
## Define indices
tg_cor = sqa_plus.options.core_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

dummy_indices_list = sqa_plus.create_dummy_indices_list()

## Define terms
print("## Create a_I^\dag a_A ...\n")

i = sqa_plus.index('I', [tg_cor])
a = sqa_plus.index('A', [tg_vir])

## Define operator types
cre_i = sqa_plus.creOp(i)
des_a = sqa_plus.desOp(a)

terms_l_op = [sqa_plus.term(0.5, [], [sqa_plus.creOp(i), sqa_plus.desOp(a)])]

# Spin-Orbital H(1)
print("# Create spin-orbital  H^(1) = V + [H^(0), T - T^\dag]:")

# Spin-Orbital Dyall H
print("## Create spin-orbital  H^(0)...")
terms_Dyall_H = sqa_plus.dyallH(dummy_indices_list)
print("{:} spin-orbital H terms created.".format(len(terms_Dyall_H)))

## Spin-Orbital T - T^\dag
print("\n## Create spin-orbital T - T^\dag...")
terms_T = sqa_plus.Tamplitude(1, dummy_indices_list)
print("{:} spin-orbital A terms created.".format(len(terms_T)))

## Commutator [H^(0), T - T^\dag]
print("\n## Calculate [H^(0), T - T^\dag]:")
terms_commutator_H0_T = sqa_plus.commutator(terms_Dyall_H, terms_T)
print("{:} spin-orbital terms created.".format(len(terms_commutator_H0_T)))

# Spin-Orbital V
print("\n## Create spin-orbital V...")
terms_V = sqa_plus.Vperturbation(dummy_indices_list)
print("{:} spin-orbital V terms created.\n".format(len(terms_V)))

## Spin-Orbital H^(1)
terms_H1 = terms_V + terms_commutator_H0_T
print("{:} spin-orbital H^(1) terms created.".format(len(terms_H1)))

# Calculate [V + H^(1), T - T^\dag]
print("\n# Create [V + H^(1), T - T^\dag]")

# Spin-Orbital V
print("\n## Create spin-orbital V...")
terms2_V = sqa_plus.Vperturbation(dummy_indices_list)
print("{:} spin-orbital V terms created.\n".format(len(terms2_V)))

## Spin-Orbital V + H^(1)
terms_V_H1 = terms2_V + terms_H1
print("{:} spin-orbital V + H^(1) terms created.".format(len(terms_V_H1)))

## Spin-Orbital T - T^\dag
print("\n## Create spin-orbital T - T^\dag...")
terms2_T = sqa_plus.Tamplitude(1, dummy_indices_list)
print("{:} spin-orbital A terms created.".format(len(terms2_T)))

## Commutator [V + H^(1), T - T^\dag]
print("\n## Calculate [V + H^(1), T - T^\dag]:")
terms_commutator_V_H1_T = sqa_plus.commutator(terms_V_H1, terms2_T)
print("{:} spin-orbital terms created.".format(len(terms_commutator_V_H1_T)))

print("\n## Calculating a^\dag_i a_a [V + H^(1), T - T^\dag]:")
terms_V1 = []
for term_commutator in terms_commutator_V_H1_T:
    for term_l in terms_l_op:
        terms_V1.append(sqa_plus.multiplyTerms(term_l, term_commutator))
print("{:} spin-orbital a^i a_a [V + H1, A1] terms created.\n".format(len(terms_V1)))

# Compute expected value of spin-orbital V matrix
print("## Compute expected value of spin-orbital V matrix ...")
terms_expected_V = sqa_plus.matrixBlock(terms_V1, transRDM = False, legacy_order = use_legacy_order)

# Create Numpy einsum equations
terms_expected_V.sort()
result = sqa_plus.genEinsum(terms_expected_V, 'V1', 'IA', rm_core_int = True, suffix = 'so')

print("\n-------------------------------- genEinsum equations --------------------------------\n")
for item in result:
    print(item)
print("\n-------------------------------------------------------------------------------------\n")

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))