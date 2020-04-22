import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = False

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

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
#l_op  = [sqa.creOp(i), sqa.desOp(x)] # CA
#l_ind = 'IX'                                                         

#l_op  = [sqa.creOp(i), sqa.desOp(a)] # CE
#l_ind = 'IA'                                                         

l_op  = [sqa.creOp(x), sqa.desOp(a)] # AE
l_ind = 'XA'                                                         

#l_op  = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(z), sqa.desOp(y)] # CCAA
#l_ind = 'JKYZ'                                                         
                                                                      
#l_op  = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(y), sqa.desOp(b)] # CCEA
#l_ind = 'JKBY'                                                         
                                                                      
#l_op  = [sqa.creOp(j), sqa.creOp(k), sqa.desOp(c), sqa.desOp(b)] # CCEE
#l_ind = 'JKBC'                                                         
                                                               
#l_op  = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(c), sqa.desOp(b)] # CAEE
#l_ind = 'JYBC'                                                         
                                                                      
#l_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(c), sqa.desOp(b)] # AAEE
#l_ind = 'YZBC'                                                         
                                                                      
#l_op  = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(u), sqa.desOp(z)] # CAAA
#l_ind = 'JYZU'                                                         
                                                                      
#l_op  = [sqa.creOp(j), sqa.creOp(y), sqa.desOp(z), sqa.desOp(b)] # CAEA
#l_ind = 'JYBZ'                                                         
                                                                      
#l_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(u), sqa.desOp(b)] # AAEA
#l_ind = 'YZBU'


# RHS
#i = sqa.index('i', [tg_c], dummy)
#x = sqa.index('x', [tg_a], dummy)
#X = [sqa.tensor('X', [i,x])]
#r_op  = [sqa.creOp(x), sqa.desOp(i)] # CA
#r_ind = 'IX'                                                         

#i = sqa.index('i', [tg_c], dummy)
#a = sqa.index('aa', [tg_v], dummy)
#X = [sqa.tensor('X', [i,a])]
#r_op  = [sqa.creOp(a), sqa.desOp(i)] # CE
#r_ind = 'IA'                                                         

#x = sqa.index('x', [tg_a], dummy)
#a = sqa.index('aa', [tg_v], dummy)
#X = [sqa.tensor('X', [x,a])]
#r_op  = [sqa.creOp(a), sqa.desOp(x)] # AE
#r_ind = 'XA'                                                         

#j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [j,k,y,z], [Xsym_1, Xsym_2])]
#r_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(k), sqa.desOp(j)] # CCAA
#r_ind = 'JKYZ'                                                         
#scaling_factor = 0.25                                  
                                                                     
#j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#y = sqa.index('y', [tg_a], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#X = [sqa.tensor('X', [j,k,b,y], [Xsym_1])]
#r_op  = [sqa.creOp(b), sqa.creOp(y), sqa.desOp(k), sqa.desOp(j)] # CCEA
#r_ind = 'JKBY'                                                         
#scaling_factor = 0.5                                  
                                                                     
#j = sqa.index('j', [tg_c], dummy)
#k = sqa.index('k', [tg_c], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [j,k,b,c], [Xsym_1, Xsym_2])]
#r_op  = [sqa.creOp(b), sqa.creOp(c), sqa.desOp(k), sqa.desOp(j)] # CCEE
#r_ind = 'JKBC'                                                         
#scaling_factor = 0.25                                  
                                                              
#j = sqa.index('j', [tg_c], dummy)
#y = sqa.index('y', [tg_a], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [j,y,b,c], [Xsym_2])]
#r_op  = [sqa.creOp(b), sqa.creOp(c), sqa.desOp(y), sqa.desOp(j)] # CAEE
#r_ind = 'JYBC'                                                         
#scaling_factor = 0.5                                  
                                                                     
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#c = sqa.index('cc', [tg_v], dummy)
#Xsym_1 = sqa.symmetry((1,0,2,3),-1)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [y,z,b,c], [Xsym_1, Xsym_2])]
#r_op  = [sqa.creOp(b), sqa.creOp(c), sqa.desOp(z), sqa.desOp(y)] # AAEE
#r_ind = 'YZBC'                                                         
#scaling_factor = 0.25                                  
                                   
#j = sqa.index('j', [tg_c], dummy)
#y = sqa.index('y', [tg_a], dummy)
#z = sqa.index('z', [tg_a], dummy)
#u = sqa.index('u', [tg_a], dummy)
#Xsym_2 = sqa.symmetry((0,1,3,2),-1)
#X = [sqa.tensor('X', [j,y,z,u], [Xsym_2])]
#r_op  = [sqa.creOp(z), sqa.creOp(u), sqa.desOp(y), sqa.desOp(j)] # CAAA
#r_ind = 'JYZU'                                                         
#scaling_factor = 0.5                                  
                                                                    
#j = sqa.index('j', [tg_c], dummy)
#y = sqa.index('y', [tg_a], dummy)
#b = sqa.index('bb', [tg_v], dummy)
#z = sqa.index('z', [tg_a], dummy)
#X = [sqa.tensor('X', [j,y,b,z])]
#r_op  = [sqa.creOp(b), sqa.creOp(z), sqa.desOp(y), sqa.desOp(j)] # CAEA
#r_ind = 'JYBZ'                                                         
#scaling_factor = 1.00                                  
                                                                    
y = sqa.index('y', [tg_a], dummy)
z = sqa.index('z', [tg_a], dummy)
b = sqa.index('bb', [tg_v], dummy)
u = sqa.index('u', [tg_a], dummy)
Xsym_1 = sqa.symmetry((1,0,2,3),-1)
X = [sqa.tensor('X', [y,z,b,u], [Xsym_1])]
r_op  = [sqa.creOp(b), sqa.creOp(u), sqa.desOp(z), sqa.desOp(y)] # AAEA
r_ind = 'YZBU'
scaling_factor = 0.5                                


# Define Hamiltonian
effH = []
effH = sqa.Heff(1)

print ("Terms of Hamiltonian")
for t in effH:
  print (t)

term1 = sqa.term(scaling_factor * 1.0, [], X + r_op)
term2 = sqa.term(1.0, [], l_op)

print ("First Commutator")
term3 = sqa.commutator(effH, term1)

print ("Second Commutator")
term4 = sqa.commutator(term2, term3)

term5 = sqa.matrixBlock(term4)
sqa.generateEinsum(term5, 'temp', l_ind, "")
