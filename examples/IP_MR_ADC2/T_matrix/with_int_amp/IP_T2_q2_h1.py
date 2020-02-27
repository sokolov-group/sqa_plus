#import sqa_extra_old.secondQuantizationAlgebra as sqa
import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

# Core dummy indices
cc = [sqa.index('c%i' %p, [tg_c], dummy) for p in range(50)]
# Active dummy indices
aa = [sqa.index('a%i' %p, [tg_a], dummy) for p in range(50)]
# Virtual dummy indices
vv = [sqa.index('v%i' %p, [tg_v], dummy) for p in range(50)]

# External indices
b = sqa.index('B', [tg_v], False)
c = sqa.index('C', [tg_v], False)
d = sqa.index('D', [tg_v], False)
w = sqa.index('W', [tg_a], False)
y = sqa.index('Y', [tg_a], False)
z = sqa.index('Z', [tg_a], False)
j = sqa.index('J', [tg_c], False)
k = sqa.index('K', [tg_c], False)
l = sqa.index('L', [tg_c], False)

# Reserved for the q operator
i = sqa.index('I', [tg_c], False)
x = sqa.index('X', [tg_a], False)
a = sqa.index('A', [tg_v], False)

##########################
# Input parameters
##########################
tensor_name = "temp"

# h^(0) operators:
#h_op = [sqa.creOp(j)]
# h^(1) operators:
h_op = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(z)]
#h_op = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(b)]
#h_op = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(b)]
#h_op = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(b)]
#h_op = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(y)]

# q operators
q_op = [sqa.desOp(i)]
#q_op = [sqa.desOp(x)]
#q_op = [sqa.desOp(a)]
q_order = 2

#target_ind = "IJ"
#target_ind = "XJ"
#target_ind = "AJ"

#target_ind = "IJYZ"
#target_ind = "XJYZ"
#target_ind = "AJYZ"

#target_ind = "IJKB"
#target_ind = "XJKB"
#target_ind = "AJKB"

#target_ind = "IJYB"
#target_ind = "XJYB"
#target_ind = "AJYB"

#target_ind = "IYZB"
#target_ind = "XYZB"
#target_ind = "AYZB"

#target_ind = "IJKY"
#target_ind = "XJKY"
target_ind = "AJKY"
##########################

# Evaluating expression
h_term = sqa.term(1.0, [], h_op)
q0 = sqa.term(1.0, [], q_op)

# Compute q tilde operators
q_terms = None
# Zeroth-order
if q_order == 0:
    q_terms = [q0]
# First-order
elif q_order == 1:
    # First-order excitation operator
    T1 = []
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(4):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1.extend(sqa.Tamplitude(T1, 1, cc1, aa1, vv1))

    q1 = []
    q1.extend(sqa.commutator(q0,T1))

    q_terms = q1
# Second-order
else: 
    # First-order excitation operator
    T1 = []
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(4):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1.extend(sqa.Tamplitude(T1, 1, cc1, aa1, vv1))

    q1 = []
    q1.extend(sqa.commutator(q0,T1))

    # First-order excitation operator (2)
    T1_ = []
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(4):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1_.extend(sqa.Tamplitude(T1_, 1, cc1, aa1, vv1))

    q2 = []
    q2.extend(sqa.commutator(q1,T1_))
    for t in q2:
      t.scale(0.5)

    # Second-order excitation operator
    T2 = []
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(4):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T2.extend(sqa.Tamplitude(T2, 2, cc1, aa1, vv1))

    q0_T2 = []
    q0_T2.extend(sqa.commutator(q0,T2))
    q2.extend(q0_T2)

    q_terms = q2

for t in q_terms:
 print t

f_terms = []

for t in q_terms:
    f_terms.append(sqa.multiplyTerms(t,h_term))
    f_terms.append(sqa.multiplyTerms(h_term,t))

t_terms = []
for t in f_terms:
    tt = sqa.normalOrder(t)
    t_terms.extend(tt)

final_result = sqa.matrixBlock(t_terms)

sqa.generateEinsum(final_result, tensor_name, target_ind)
#sqa.generateEinsum(final_result, 'T[:ncore_so, si:fi]', 'P', transRDM=True, trans_ind_str='P')

