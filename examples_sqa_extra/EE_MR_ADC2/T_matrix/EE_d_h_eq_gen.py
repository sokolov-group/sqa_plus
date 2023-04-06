import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = False

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# Core dummy indices
cc = [sqa.index('c%i' %p, [tg_c], dummy) for p in range(500)]

# Active dummy indices
aa = [sqa.index('a%i' %p, [tg_a], dummy) for p in range(500)]

# Virtual dummy indices
vv = [sqa.index('v%i' %p, [tg_v], dummy) for p in range(500)]

# External indices
i = sqa.index('I', [tg_c], False)
j = sqa.index('J', [tg_c], False)
k = sqa.index('K', [tg_c], False)
l = sqa.index('L', [tg_c], False)

x = sqa.index('X', [tg_a], False)
y = sqa.index('Y', [tg_a], False)
z = sqa.index('Z', [tg_a], False)
u = sqa.index('U', [tg_a], False)

a = sqa.index('A', [tg_v], False)
b = sqa.index('B', [tg_v], False)
c = sqa.index('C', [tg_v], False)
d = sqa.index('D', [tg_v], False)

##########################
# Input parameters
##########################

# Define d order
#######################
d_order = 0
#d_order = 1
#d_order = 2


# h^(0) operators:
#######################
# CA
#h_op = [sqa.creOp(x), sqa.desOp(i)]
#target_ind = 'IX'

# CE
#h_op = [sqa.creOp(a), sqa.desOp(i)]
#target_ind = 'IA'

# AE
#h_op = [sqa.creOp(a), sqa.desOp(x)]
#target_ind = 'XA'


# h^(1) operators:
#######################
# CCAA
#h_op = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(j), sqa.desOp(i)]
#target_ind = 'IJXY'

# CCEE
#h_op = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(j), sqa.desOp(i)]
#target_ind = 'IJAB'

# CCEA
#h_op = [sqa.creOp(a), sqa.creOp(x), sqa.desOp(j), sqa.desOp(i)]
#target_ind = 'IJAX'

# CAEE
#h_op = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(x), sqa.desOp(i)]
#target_ind = 'IXAB'

# AAEE
#h_op = [sqa.creOp(a), sqa.creOp(b), sqa.desOp(y), sqa.desOp(x)]
#target_ind = 'XYAB'

# CAAA
#h_op = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(x), sqa.desOp(i)]
#target_ind = 'IXYZ'

# CAEA
#h_op = [sqa.creOp(a), sqa.creOp(y), sqa.desOp(x), sqa.desOp(i)]
#target_ind = 'IXAY'

# AAEA
h_op = [sqa.creOp(a), sqa.creOp(z), sqa.desOp(y), sqa.desOp(x)]
target_ind = 'XYAZ'


##########################

# Define bra-ket symmetry
dSym = sqa.symmetry((1,0),1)

# Construct d operator
d_op = []

# CC
p = cc.pop(0)
q = cc.pop(0)
dTen = sqa.tensor('d_cc_so', [p,q], [dSym])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(p), sqa.desOp(q)]))

# AA
p = aa.pop(0)
q = aa.pop(0)
dTen = sqa.tensor('d_aa_so', [p,q], [dSym])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(p), sqa.desOp(q)]))

# EE
p = vv.pop(0)
q = vv.pop(0)
dTen = sqa.tensor('d_ee_so', [p,q], [dSym])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(p), sqa.desOp(q)]))

# CA
p = cc.pop(0)
q = aa.pop(0)
dTen = sqa.tensor('d_ca_so', [p,q])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(p), sqa.desOp(q)]))
p = cc.pop(0)
q = aa.pop(0)
dTen = sqa.tensor('d_ca_so', [p,q])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(q), sqa.desOp(p)]))

# CE
p = cc.pop(0)
q = vv.pop(0)
dTen = sqa.tensor('d_ce_so', [p,q])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(p), sqa.desOp(q)]))
p = cc.pop(0)
q = vv.pop(0)
dTen = sqa.tensor('d_ce_so', [p,q])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(q), sqa.desOp(p)]))

# AE
p = aa.pop(0)
q = vv.pop(0)
dTen = sqa.tensor('d_ae_so', [p,q])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(p), sqa.desOp(q)]))
p = aa.pop(0)
q = vv.pop(0)
dTen = sqa.tensor('d_ae_so', [p,q])
d_op.append(sqa.term(1.0, [], [dTen, sqa.creOp(q), sqa.desOp(p)]))

##########################

# Evaluating expression
h_term = sqa.term(1.0, [], h_op)
d0 = d_op

# Compute d tilde operators
d_terms = None

# Zeroth-order
if d_order == 0:
    d_terms = d0

# First-order
elif d_order == 1:

    # First-order excitation operator
    T1 = []
    cc1 = []
    aa1 = []
    vv1 = []

    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))

    T1.extend(sqa.Tamplitude(1, cc1, aa1, vv1))

    d1 = []
    d1.extend(sqa.commutator(d0,T1))
    d_terms = d1

# Second-order
else: 

    # First-order excitation operator
    T1 = []
    cc1 = []
    aa1 = []
    vv1 = []

    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))

    T1.extend(sqa.Tamplitude(1, cc1, aa1, vv1))

    d1 = []
    d1.extend(sqa.commutator(d0,T1))

    # First-order excitation operator (2)
    T1_ = []
    cc1 = []
    aa1 = []
    vv1 = []

    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))

    T1_.extend(sqa.Tamplitude(1, cc1, aa1, vv1))

    d2 = []
    d2.extend(sqa.commutator(d1,T1_))

    for t in d2:
      t.scale(0.5)

    # Second-order excitation operator
    T2 = []
    cc1 = []
    aa1 = []
    vv1 = []

    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T2.extend(sqa.Tamplitude(2, cc1, aa1, vv1))

    d0_T2 = []
    d0_T2.extend(sqa.commutator(d0,T2))
    d2.extend(d0_T2)

    d_terms = d2

#######################

# Commute 'd' and 'h' terms
f_terms = sqa.commutator(d_terms, h_term)

# Prepare final result
final_result = sqa.matrixBlock(f_terms)

tensor_name = 'temp'
einsum_list = sqa.genEinsum(final_result, tensor_name, target_ind, rm_core_int = True)

for einsum in einsum_list:
    print(einsum)
