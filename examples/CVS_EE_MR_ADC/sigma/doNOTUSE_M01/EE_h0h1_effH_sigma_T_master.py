#import sqa_extra.secondQuantizationAlgebra as sqa
import sqa_plus as sqa
sqa.options.cvs_approach = True

sqa.options.verbose = False

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type

tg_cvs_cor = sqa.options.cvs_core_type
tg_cvs_val = sqa.options.cvs_valence_type

# Define indices
dummy = True

# External indices
i = sqa.index('I', [tg_cvs_cor])
j = sqa.index('J', [tg_cvs_cor])
k = sqa.index('K', [tg_cvs_cor])
l = sqa.index('L', [tg_cvs_cor])
m = sqa.index('M', [tg_cvs_cor])
n = sqa.index('N', [tg_cvs_cor])

i_val = sqa.index('I', [tg_cvs_val])
j_val = sqa.index('J', [tg_cvs_val])
k_val = sqa.index('K', [tg_cvs_val])
l_val = sqa.index('L', [tg_cvs_val])
m_val = sqa.index('M', [tg_cvs_val])
n_val = sqa.index('N', [tg_cvs_val])

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
#i = sqa.index('i', [tg_c], dummy)
#x = sqa.index('x', [tg_a], dummy)
#X = [sqa.tensor('X', [i,x])]
#l_op  = [sqa.creOp(i), sqa.desOp(x)] # CA
#l_ind = 'IX'                                                         
#scaling_factor = 1.0

#i = sqa.index('i', [tg_c], dummy)
#a = sqa.index('aa', [tg_v], dummy)
#X = [sqa.tensor('X', [i,a])]
#l_op  = [sqa.creOp(i), sqa.desOp(a)] # CE
#l_ind = 'IA'                                                         
#scaling_factor = 1.0

#x = sqa.index('x', [tg_a], dummy)
#a = sqa.index('aa', [tg_v], dummy)
#X = [sqa.tensor('X', [x,a])]
#l_op  = [sqa.creOp(x), sqa.desOp(a)] # AE
#l_ind = 'XA'                                                         
#scaling_factor = 1.0

#j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [j,k,y,z], [Xsym_1, Xsym_2])]
#l_op  = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(z), sqa.desOp(y)] # CCAA
#l_ind = 'JKYZ'                                                         
#scaling_factor = 0.25
                                                                     
#j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#y = sqa.index('y', [tg_a], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#X = [sqa.tensor('X', [j,k,b,y], [Xsym_1])]
#l_op  = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(y), sqa.desOp(b)] # CCEA
#l_ind = 'JKBY'                                                         
#scaling_factor = 0.5
                                                                     
#j = sqa.index('j', [tg_cvs_cor], dummy)
#k = sqa.index('k', [tg_cvs_cor], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#Xsym_3 = sqa.symmetry((1,0,3,2),1)
#X = [sqa.tensor('X', [j,k,b,c], [Xsym_1, Xsym_2, Xsym_3])]
#l_op  = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(c), sqa.desOp(b)] # CCEE
#l_ind = 'JKBC'                                                         
#scaling_factor = 0.25
                                                              
#j = sqa.index('j', [tg_c], dummy)
#y = sqa.index('y', [tg_a], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [j,y,b,c], [Xsym_2])]
#l_op  = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(c), sqa.desOp(b)] # CAEE
#l_ind = 'JYBC'                                                         
#scaling_factor = 0.5
                                                                     
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [y,z,b,c], [Xsym_1, Xsym_2])]
#l_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(c), sqa.desOp(b)] # AAEE
#l_ind = 'YZBC'                                                         
#scaling_factor = 0.25
                                                                     
#j = sqa.index('j', [tg_c], dummy)
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#u = sqa.index('u', [tg_a], dummy)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [j,y,z,u], [Xsym_2])]
#l_op  = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(u), sqa.desOp(z)] # CAAA
#l_ind = 'JYZU'                                                         
#scaling_factor = 0.5
                                                                    
j = sqa.index('j', [tg_c], dummy)
y = sqa.index('y', [tg_a], dummy)
b = sqa.index('bb', [tg_v], dummy)
z = sqa.index('z', [tg_a], dummy)
X = [sqa.tensor('X', [j,y,b,z])]
l_op  = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(z), sqa.desOp(b)] # CAEA
l_ind = 'JYBZ'                                                         
scaling_factor = 1.0
                                                                    
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#u = sqa.index('u', [tg_a], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#X = [sqa.tensor('X', [y,z,b,u], [Xsym_1])]
#l_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(u), sqa.desOp(b)] # AAEA
#l_ind = 'YZBU'
#scaling_factor = 0.5

# RHS
#r_op  = [sqa.creOp(x), sqa.desOp(i)] # CA
#r_ind = 'IX'                                                         

r_op  = [sqa.creOp(a), sqa.desOp(i)] # CE
r_ind = 'IA'                                                         

#r_op  = [sqa.creOp(a), sqa.desOp(x)] # AE
#r_ind = 'XA'                                                         

#r_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(k), sqa.desOp(j)] # CCAA
#r_ind = 'JKYZ'                                                         
                                                                      
#r_op  = [sqa.creOp(b), sqa.creOp(y), sqa.desOp(k), sqa.desOp(j)] # CCEA
#r_ind = 'JKBY'                                                         

#r_op  = [sqa.creOp(b), sqa.creOp(c), sqa.desOp(k), sqa.desOp(j)] # CCEE
#r_ind = 'JKBC'                                                         
                                                               
#r_op  = [sqa.creOp(b), sqa.creOp(c), sqa.desOp(y), sqa.desOp(j)] # CAEE
#r_ind = 'JYBC'                                                         
                                                                      
#r_op  = [sqa.creOp(b), sqa.creOp(c), sqa.desOp(z), sqa.desOp(y)] # AAEE
#r_ind = 'YZBC'                                                         
                                                                      
#r_op  = [sqa.creOp(z), sqa.creOp(u), sqa.desOp(y), sqa.desOp(j)] # CAAA
#r_ind = 'JYZU'                                                         
                                                                      
#r_op  = [sqa.creOp(b), sqa.creOp(z), sqa.desOp(y), sqa.desOp(j)] # CAEA
#r_ind = 'JYBZ'                                                         
                                                                      
#r_op  = [sqa.creOp(b), sqa.creOp(u), sqa.desOp(z), sqa.desOp(y)] # AAEA
#r_ind = 'YZBU'


# Define Hamiltonian
effH = []
effH = sqa.Heff(1)

for t in effH:
  print (t)

term1 = sqa.term(1.0, [], r_op)
term2 = sqa.term(scaling_factor * 1.0, [], l_op + X)

print ("First Commutator")

term3 = sqa.commutator(effH, term1)
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)

sqa.genEinsum(term5, 'temp', r_ind, "")
