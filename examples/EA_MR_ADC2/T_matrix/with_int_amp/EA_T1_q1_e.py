import sqa_extra.secondQuantizationAlgebra as sqa

#sqa.options.verbose = True

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
# h^(0) operators
h_op = [sqa.creOp(b)]

# q operators
q_op = [sqa.desOp(i)]
#q_op = [sqa.desOp(x)]
#q_op = [sqa.desOp(a)]

q_order = 1

tensor_name = "temp"

target_ind = "IB"
#target_ind = "XB"
#target_ind = "AB"
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
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1.extend(sqa.Tamplitude(1, cc1, aa1, vv1))

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
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1.extend(sqa.Tamplitude(1, cc1, aa1, vv1))

    q1 = []
    q1.extend(sqa.commutator(q0,T1))

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

    q2 = []
    q2.extend(sqa.commutator(q1,T1_))
    for t in q2:
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

