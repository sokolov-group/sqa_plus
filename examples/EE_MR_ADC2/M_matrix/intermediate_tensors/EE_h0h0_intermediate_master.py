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

x = sqa.index('X', [tg_a])
y = sqa.index('Y', [tg_a])
z = sqa.index('Z', [tg_a])
u = sqa.index('U', [tg_a])

a = sqa.index('A', [tg_v])
b = sqa.index('B', [tg_v])
c = sqa.index('C', [tg_v])
d = sqa.index('D', [tg_v])

# LHS
#l_op  = [sqa.creOp(i), sqa.desOp(x)] # CA
#l_ind = 'IX'                                                         

#l_op  = [sqa.creOp(i), sqa.desOp(a)] # CE
#l_ind = 'IA'                                                         

l_op  = [sqa.creOp(x), sqa.desOp(a)] # AE
l_ind = 'XA'                                                         


# RHS
#r_op  = [sqa.creOp(y), sqa.desOp(j)] # CA
#r_ind = 'JY'                                                         

#r_op  = [sqa.creOp(b), sqa.desOp(j)] # CE
#r_ind = 'JB'                                                         

r_op  = [sqa.creOp(b), sqa.desOp(y)] # AE
r_ind = 'YB'                                                         

# Define terms w/ desired operators
l_op_term = sqa.term(1.0, [], l_op)
r_op_term = sqa.term(1.0, [], r_op)

# Define Hamiltonian
effH = []
effH = sqa.Heff(1)

# Perform multiplication of terms and normal-ordering
term1 = []
for t in effH:
    tt = sqa.multiplyTerms(l_op_term, t)
    term1.append(tt)

term2 = []
for t in term1:
    tt = sqa.multiplyTerms(t, r_op_term)
    term2.append(tt)

term3 = []
for t in term2:
    tt = sqa.normalOrder(t)
    term3.extend(tt)

# Filter through all terms keep only terms that involve 4-rdm
term4 = []
for term in term3:
    cre_des = 0
    for tens in term.tensors:
        if isinstance(tens, sqa.creOp) or isinstance(tens, sqa.desOp):
            cre_des += 1 
    if cre_des == 2:
        term4.append(term)

term5 = sqa.matrixBlock(term4)

# Original terms (no intermediates)
einsum_list_1 = sqa.genEinsum(term5, 'temp', l_ind + r_ind, rm_core_int = True)

for einsum in einsum_list_1:
    print (einsum)
print ('')

##############################################
# Scan for intermediates w/ factor_depth = 1
# INT terms will only be defined in terms of tensors in original equations
term6, intermediates = sqa.genIntermediates(term5, l_ind + r_ind, factor_depth = 1)

# Generate einsum for intermediates and new terms
int_einsum_list, einsum_list_2 = sqa.genEinsum(term6, 'temp', l_ind + r_ind, rm_core_int = True, intermediate_list = intermediates)

# Print intermediate term definitions
for int_def in int_einsum_list:
    print (int_def)
print ('')

# Print terms w/ INT terms substitution
for einsum in einsum_list_2:
    print (einsum)

##############################################
# Scan for intermediates w/ factor_depth = 2
# INT terms can be defined in terms of one other INT terms
term6, intermediates = sqa.genIntermediates(term5, l_ind + r_ind, factor_depth = 2)

# Generate einsum for intermediates and new terms
int_einsum_list, einsum_list_2 = sqa.genEinsum(term6, 'temp', l_ind + r_ind, rm_core_int = True, intermediate_list = intermediates)

# Print intermediate term definitions
for int_def in int_einsum_list:
    print (int_def)
print ('')

# Print terms w/ INT terms substitution
for einsum in einsum_list_2:
    print (einsum)
