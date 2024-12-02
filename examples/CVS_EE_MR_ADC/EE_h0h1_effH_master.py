#import sqa_extra.secondQuantizationAlgebra as sqa
import sqa_plus as sqa

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

l_op  = [sqa.creOp(k), sqa.creOp(l), sqa.desOp(w), sqa.desOp(u)] # CCAA
l_ind = 'KLUW'                                                         

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
#r_op  = [sqa.creOp(x), sqa.desOp(i)] # CA
#r_ind = 'IX'                                                         

#r_op  = [sqa.creOp(a), sqa.desOp(i)] # CE
#r_ind = 'IA'                                                         

#r_op  = [sqa.creOp(a), sqa.desOp(x)] # AE
#r_ind = 'XA'                                                         

r_op  = [sqa.creOp(x), sqa.creOp(y), sqa.desOp(j), sqa.desOp(i)] # CCAA
r_ind = 'IJXY'                                                         

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
                                                                     
#r_op  = [sqa.creOp(c), sqa.creOp(u), sqa.desOp(w), sqa.desOp(k)] # CAEA
#r_ind = 'KWCU'                                                     
                                                                     
#r_op  = [sqa.creOp(b), sqa.creOp(u), sqa.desOp(z), sqa.desOp(y)] # AAEA
#r_ind = 'YZBU'

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

from sqa_plus.sqaEinsum import remove_core_int
term6, removed_core = remove_core_int(term5)
sqa.options.genEinsum.remove_core_integrals = False
sqa.options.genEinsum.keep_user_defined_dummy_names = True

diagonal_pairs_dict = {'I': 'K', 'J': 'L'}

# Replacing indices to obtain the diagonal
for term_ind, term_sa in enumerate(term6):
    for tensor_ind, tensor_sa in enumerate(term_sa.tensors):
        for index_ind, index_sa in enumerate(tensor_sa.indices):
            if index_sa.name in diagonal_pairs_dict.keys():
                term6[term_ind].tensors[tensor_ind].indices[index_ind].userDefined = True
                term6[term_ind].tensors[tensor_ind].indices[index_ind].name = diagonal_pairs_dict[index_sa.name]

from sqa_plus.sqaMatrixBlock import contractDeltaFuncs_nondummy
term6 = contractDeltaFuncs_nondummy(term6)

sqa.combineTerms(term6)

sqa.genEinsum(term6, 'temp', 'KLUW')
