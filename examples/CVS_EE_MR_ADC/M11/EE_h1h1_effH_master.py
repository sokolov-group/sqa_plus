#import sqa_extra.secondQuantizationAlgebra as sqa

#sqa.options.verbose = False

import sqa_plus
sqa_plus.options.spin_integrated = True
sqa_plus.options.cvs_approach = True
order_Heff = 0

import time
start = time.time()

sqa_plus.options.print_header("Spin-Adapted CVS-EE: M11 H{:}".format(order_Heff))

## Define indices
tg_cvs_cor = sqa_plus.options.cvs_core_type
tg_cvs_val = sqa_plus.options.cvs_valence_type
tg_act = sqa_plus.options.active_type
tg_vir = sqa_plus.options.virtual_type

tg_alpha = sqa_plus.options.alpha_type
tg_beta = sqa_plus.options.beta_type

## External Indices
c_alpha = sqa_plus.index('C', [tg_alpha, tg_vir])
c_beta  = sqa_plus.index('C', [tg_beta,  tg_vir])

d_alpha = sqa_plus.index('D', [tg_alpha, tg_vir])
d_beta  = sqa_plus.index('D', [tg_beta,  tg_vir])

k_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_cor])
k_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_cor])

k_val_alpha = sqa_plus.index('K', [tg_alpha, tg_cvs_val])
k_val_beta  = sqa_plus.index('K', [tg_beta,  tg_cvs_val])

l_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_cor])
l_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_cor])

l_val_alpha = sqa_plus.index('L', [tg_alpha, tg_cvs_val])
l_val_beta  = sqa_plus.index('L', [tg_beta,  tg_cvs_val])

w_alpha = sqa_plus.index('W', [tg_alpha, tg_act])
w_beta  = sqa_plus.index('W', [tg_beta,  tg_act])

z_alpha = sqa_plus.index('Z', [tg_alpha, tg_act])
z_beta  = sqa_plus.index('Z', [tg_beta,  tg_act])

a_alpha = sqa_plus.index('a', [tg_alpha, tg_vir])
a_beta  = sqa_plus.index('a', [tg_beta,  tg_vir])

b_alpha = sqa_plus.index('b', [tg_alpha, tg_vir])
b_beta  = sqa_plus.index('b', [tg_beta,  tg_vir])

i_alpha = sqa_plus.index('i', [tg_alpha, tg_cvs_cor])
i_beta  = sqa_plus.index('i', [tg_beta,  tg_cvs_cor])

i_val_alpha = sqa_plus.index('i', [tg_alpha, tg_cvs_val])
i_val_beta  = sqa_plus.index('i', [tg_beta,  tg_cvs_val])

j_alpha = sqa_plus.index('j', [tg_alpha, tg_cvs_cor])
j_beta  = sqa_plus.index('j', [tg_beta,  tg_cvs_cor])

j_val_alpha = sqa_plus.index('j', [tg_alpha, tg_cvs_val])
j_val_beta  = sqa_plus.index('j', [tg_beta,  tg_cvs_val])

x_alpha = sqa_plus.index('x', [tg_alpha, tg_act])
x_beta  = sqa_plus.index('x', [tg_beta,  tg_act])

y_alpha = sqa_plus.index('y', [tg_alpha, tg_act])
y_beta  = sqa_plus.index('y', [tg_beta,  tg_act])

# External indices
i = sqa.index('I', [tg_c])
j = sqa.index('J', [tg_c])
k = sqa.index('K', [tg_c])
l = sqa.index('L', [tg_c])
m = sqa.index('M', [tg_c])
n = sqa.index('N', [tg_c])

x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])
u = sqa.index('U', [tg_a])
v = sqa.index('V', [tg_a])
w = sqa.index('W', [tg_a])

a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
c = sqa.index('C', [tg_v])
d = sqa.index('D', [tg_v])
e = sqa.index('E', [tg_v])
f = sqa.index('F', [tg_v])


# LHS
#l_op  = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(y), sqa.desOp(x)] # CCAA
#l_ind = 'IJXY'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(x), sqa.desOp(a)] # CCEA
#l_ind = 'IJAX'                                                         
                                                                       
l_op  = [sqa.creOp(i), sqa.creOp(j), sqa.desOp(b), sqa.desOp(a)] # CCEE
l_ind = 'IJAB'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(b), sqa.desOp(a)] # CAEE
#l_ind = 'IXAB'                                                         
                                                                       
#l_op  = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(b), sqa.desOp(a)] # AAEE
#l_ind = 'XYAB'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(z), sqa.desOp(y)] # CAAA
#l_ind = 'IXYZ'                                                         
                                                                       
#l_op  = [sqa.creOp(i), sqa.creOp(x), sqa.desOp(y), sqa.desOp(a)] # CAEA
#l_ind = 'IXAY'                                                         
                                                                       
#l_op  = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(z), sqa.desOp(a)] # AAEA
#l_ind = 'XYAZ'


# RHS
#r_op  = [sqa.creOp(u), sqa.creOp(v), sqa.desOp(m), sqa.desOp(l)] # CCAA  
#r_ind = 'LMUV'

#r_op  = [sqa.creOp(d), sqa.creOp(u), sqa.desOp(m), sqa.desOp(l)] # CCEA
#r_ind = 'LMDU'

r_op  = [sqa.creOp(d), sqa.creOp(e), sqa.desOp(m), sqa.desOp(l)] # CCEE
r_ind = 'LMDE'

#r_op  = [sqa.creOp(d), sqa.creOp(e), sqa.desOp(u), sqa.desOp(l)] # CAEE
#r_ind = 'LUDE'

#r_op  = [sqa.creOp(d), sqa.creOp(e), sqa.desOp(v), sqa.desOp(u)] # AAEE
#r_ind = 'UVDE'

#r_op  = [sqa.creOp(v), sqa.creOp(w), sqa.desOp(u), sqa.desOp(l)] # CAAA
#r_ind = 'LUVW'

#r_op  = [sqa.creOp(d), sqa.creOp(v), sqa.desOp(u), sqa.desOp(l)] # CAEA
#r_ind = 'LUDV'

r_op  = [sqa.creOp(d), sqa.creOp(w), sqa.desOp(v), sqa.desOp(u)] # AAEA
r_ind = 'UVDW'

# Define Hamiltonian
effH = []
effH = sqa.Heff(0)

for t in effH:
  print (t)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(1.0, [], l_op)

print ("First Commutator")

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

sqa.generateEinsum(term5, 'temp', str(l_ind + r_ind), "")
