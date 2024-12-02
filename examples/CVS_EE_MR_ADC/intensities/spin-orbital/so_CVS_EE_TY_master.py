import sqa_plus
sqa_plus.options.cvs_approach = True
from sqa_plus.sqaIndexList import indexLists
q_order = 1

import time
start = time.time()

indices_string = 'q_ca'

sqa_plus.options.print_header("Spin-Orbital CVS-EE: TY Q{:} {:}".format(q_order, indices_string.upper()))

# Generating operators
print("\n## Generating operators ...\n")

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type
tg_g = tg_cvs_cor + tg_cvs_val + tg_act + tg_vir

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

dummy = True

## External Indices
r  = sqa_plus.index('R', [tg_g])

#i = sqa_plus.index('I', [tg_cvs_cor])
#x = sqa_plus.index('X', [tg_act])
#a = sqa_plus.index('A', [tg_vir])

##Dummy Indices
cvs_inds = indexLists.cvs_core
val_inds = indexLists.cvs_valence
act_inds = indexLists.active
vir_inds = indexLists.virtual

## Define terms
term_q0 = []
#dSym = [sqa_plus.symmetry((1,0),1)]
dSym = [sqa_plus.symmetry((0,2,1),1)]

# CC
p = cvs_inds.new_index()
q = cvs_inds.new_index()

dTen = sqa_plus.tensor('d_cc_so', [r, p, q], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))

# VV
p = val_inds.new_index()
q = val_inds.new_index()

dTen = sqa_plus.tensor('d_vv_so', [r, p, q], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))

# AA
p = act_inds.new_index()
q = act_inds.new_index()

dTen = sqa_plus.tensor('d_aa_so', [r, p, q], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))

# EE
p = vir_inds.new_index()
q = vir_inds.new_index()

dTen = sqa_plus.tensor('d_ee_so', [r, p, q], dSym)
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))

# CV
p = cvs_inds.new_index()
q = val_inds.new_index()

dTen = sqa_plus.tensor('d_cv_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))
dTen = sqa_plus.tensor('d_cv_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q), sqa_plus.desOp(p)]))

# CA
p = cvs_inds.new_index()
q = act_inds.new_index()

dTen = sqa_plus.tensor('d_ca_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))
dTen = sqa_plus.tensor('d_ca_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q), sqa_plus.desOp(p)]))

# VA
p = val_inds.new_index()
q = act_inds.new_index()

dTen = sqa_plus.tensor('d_va_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))
dTen = sqa_plus.tensor('d_va_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q), sqa_plus.desOp(p)]))

# CE
p = cvs_inds.new_index()
q = vir_inds.new_index()

dTen = sqa_plus.tensor('d_ce_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))
dTen = sqa_plus.tensor('d_ce_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q), sqa_plus.desOp(p)]))

# VE
p = val_inds.new_index()
q = vir_inds.new_index()

dTen = sqa_plus.tensor('d_ve_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))
dTen = sqa_plus.tensor('d_ve_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q), sqa_plus.desOp(p)]))

# AE
p = act_inds.new_index()
q = vir_inds.new_index()

dTen = sqa_plus.tensor('d_ae_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(p), sqa_plus.desOp(q)]))
dTen = sqa_plus.tensor('d_ae_so', [r, p, q])
term_q0.append(sqa_plus.term(1.0, [], [dTen, sqa_plus.creOp(q), sqa_plus.desOp(p)]))

if indices_string in ['q_ca']:

    i = cvs_inds.new_index()
    x = act_inds.new_index()

    Y = sqa_plus.tensor('Y_KW', [i, x])
    term_h  = sqa_plus.term(1.0, [], [Y, sqa_plus.creOp(x), sqa_plus.desOp(i)])
    
    final_indices_string = 'R'

elif indices_string in ['q_ce']:

    i = cvs_inds.new_index()
    a = vir_inds.new_index()
 
    Y = sqa_plus.tensor('Y_KC', [i, a]) 
    term_h  = sqa_plus.term(1.0, [], [Y, sqa_plus.creOp(a), sqa_plus.desOp(i)])
    
    final_indices_string = 'R'

## Create Dummy Indices List
if q_order == 0:
    print("## Calculating the excitation operator q^(0) ...")
    # Zeroth-order excitation operator
    terms_q = term_q0

elif q_order == 1:
    # First-order excitation operator
    print("## Calculating the excitation operator [q^(0), T - T^\dag] ...")
    terms_T = sqa_plus.Tamplitude(1)

    ## Calculating the commutator
    print("## Calculating the commutator...")

    terms_q = sqa_plus.commutator(term_q0, terms_T, combine = True)

elif q_order == 2:
    # First-order excitation operator
    print("## Calculating the excitation operator [q^(0), T - T^\dag] ...")
    terms_T1 = sqa_plus.Tamplitude(1)

    print("## Calculating the commutator...")
    terms_q = sqa_plus.commutator(term_q0, terms_T1)

    print("## Calculating the excitation operator 1/2 * [[q^(0), T - T^\dag], T - T^\dag] ...")
    terms_T1_2 = sqa_plus.Tamplitude(1)

    print("## Calculating the commutator...")
    terms_q = sqa_plus.commutator(terms_q, terms_T1_2)

    for term_q in terms_q:
      term_q.scale(0.5)

    # First-order excitation operator
    print("## Calculating the excitation operator [q^(0), T^(2) - T^(2)^\dag] ...")
    terms_T2 = sqa_plus.Tamplitude(2)

    print("## Calculating the commutator...")
    terms_q_2 = sqa_plus.commutator(term_q0, terms_T2)

    terms_q.extend(terms_q_2)

print("\n## Calculating h [q, T - T^\dag] ...")

terms_EE_TY = sqa_plus.commutator(terms_q, term_h, combine = True)

# Expected value of Spin-Adapted EE T
expected_EE_TY = sqa_plus.matrixBlock(terms_EE_TY)

#### Spin-Adaptation of EE T
###expected_EE_TY_sa = sqa_plus.convertSpinIntegratedToAdapted(expected_EE_TY)

# Generating Numpy einsum equations
#result = sqa_plus.genEinsum(expected_EE_TY_sa, 'TY_' + indices_string, final_indices_string)
result = sqa_plus.genEinsum(expected_EE_TY, 'TY_' + indices_string, final_indices_string)

end = time.time()
print("> Total elapsed time: {:.2f} seconds.".format(end - start))
