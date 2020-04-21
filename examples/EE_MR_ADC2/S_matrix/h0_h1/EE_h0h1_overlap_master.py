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

#l_op  = [sqa.creOp(x), sqa.desOp(a)] # AE
#l_ind = 'XA'                                                         

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
                                                                      
l_op  = [sqa.creOp(y), sqa.creOp(z), sqa.desOp(u), sqa.desOp(b)] # AAEA
l_ind = 'YZBU'


# RHS
#r_op  = [sqa.creOp(x), sqa.desOp(i)] # CA
#r_ind = 'IX'                                                         

#r_op  = [sqa.creOp(a), sqa.desOp(i)] # CE
#r_ind = 'IA'                                                         

r_op  = [sqa.creOp(a), sqa.desOp(x)] # AE
r_ind = 'XA'                                                         

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


term1 = sqa.term(1.0, [], l_op)
term2 = sqa.term(1.0, [], r_op)


term3 = sqa.commutator(term1, term2)
term4 = sqa.matrixBlock(term3)

sqa.generateEinsum(term4, 'temp', str(l_ind + r_ind), "")
