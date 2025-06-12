import sqa_plus
#sqa_plus.options.verbose = True
sqa_plus.options.spin_integrated = True
sqa_plus.options.genEinsum.remove_trans_rdm_constant = True
sqa_plus.options.genEinsum.remove_core_integrals = True
sqa_plus.options.genEinsum.trans_rdm = True
sqa_plus.options.genEinsum.trans_indices_string = ""

import time
start = time.time()

# Create spin-integrated V and T operators
print("# Create spin-integrated V operator ...")
terms_V = sqa_plus.Vperturbation()
terms_T = sqa_plus.Tamplitude(1, only_excitations = True)

# Multiply V * T
print("\n## Calculate 0.5 * <|V * T|> ...")
terms_VT = []
for term_T in terms_T:
    for term_V in terms_V:
        terms_VT.append(sqa_plus.multiplyTerms(term_V, term_T))

for t in terms_VT:
    t.scale(0.5)

# Evaluate term in spin-integrated form
terms_VT_si = sqa_plus.matrixBlock(terms_VT)

# Spin-adapt the resulting expression
terms_VT_sa = sqa_plus.convertSpinIntegratedToAdapted(terms_VT_si)

# Generate numpy code 
result = sqa_plus.genEinsum(terms_VT_sa, "")

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))