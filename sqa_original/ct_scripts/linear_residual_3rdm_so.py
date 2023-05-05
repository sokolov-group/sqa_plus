#    file:  linear_residual_3rdm_so.py
#  author:  Eric Neuscamman
#    date:  March 31, 2009
# summary:  Computes tensor contractions for <[H,r]>
#
# (c) 2008-2009 Eric Neuscamman (eric.neuscamman@gmail.com)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# In addition, any modification or use of this software should
# cite the following paper:
#
#   E. Neuscamman, T. Yanai, and G. K.-L. Chan.
#   J. Chem. Phys. 130, 124102 (2009)

import secondQuantizationAlgebra as sqa
import time

startTime = time.time()

# Define indices
alphaIndices = [sqa.index('p%i' %i, [sqa.options.alpha_type], True) for i in range(20)]
betaIndices  = [sqa.index('q%i' %i, [sqa.options.beta_type], True) for i in range(20)]

# Define density matrix symmetries
d1sym = [sqa.symmetry((1,0),1)]
d2sym_aaaa = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((1,0,3,2),1), sqa.symmetry((2,3,0,1),1)]
d2sym_baba = [sqa.symmetry((2,3,0,1),1)]
d3sym_aaaaaa = [sqa.symmetry((1,0,2,3,4,5),-1), sqa.symmetry((0,2,1,3,4,5),-1), sqa.symmetry((3,4,5,0,1,2),1)]
d3sym_baabaa = [sqa.symmetry((0,2,1,3,4,5),-1), sqa.symmetry((3,4,5,0,1,2),1)]

# Define amplitude symmetries
a1sym = []
a2sym_aaaa = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((1,0,3,2),1)]
a2sym_baba = []

# Define density matrices
d1aa = sqa.tensor('d1aa', alphaIndices[0:2], d1sym)
d1bb = sqa.tensor('d1bb', betaIndices[0:2], d1sym)
d2aaaa = sqa.tensor('d2aaaa', alphaIndices[0:4], d2sym_aaaa)
d2bbbb = sqa.tensor('d2bbbb', betaIndices[0:4], d2sym_aaaa)
d2baba = sqa.tensor('d2baba', betaIndices[0:1] + alphaIndices[1:2] + betaIndices[2:3] + alphaIndices[3:4], d2sym_baba)
d3aaaaaa = sqa.tensor('d3aaaaaa', alphaIndices[0:6], d3sym_aaaaaa)
d3bbbbbb = sqa.tensor('d3bbbbbb',  betaIndices[0:6], d3sym_aaaaaa)
d3baabaa = sqa.tensor('d3baabaa', betaIndices[0:1] + alphaIndices[1:3] + betaIndices[3:4] + alphaIndices[4:6], d3sym_baabaa)
d3abbabb = sqa.tensor('d3abbabb', alphaIndices[0:1] + betaIndices[1:3] + alphaIndices[3:4] + betaIndices[4:6], d3sym_baabaa)

# Define Hamiltonian
h1aa = sqa.term(1.0, [], [sqa.tensor('h1aa', alphaIndices[0:2], d1sym), sqa.creOp(alphaIndices[0]), sqa.desOp(alphaIndices[1])])
h1bb = sqa.term(1.0, [], [sqa.tensor('h1bb', betaIndices[0:2], d1sym), sqa.creOp(betaIndices[0]), sqa.desOp(betaIndices[1])])
h2aaaa = sqa.term(0.25, [], [sqa.tensor('h2aaaa', alphaIndices[2:6], d2sym_aaaa), \
                             sqa.creOp(alphaIndices[2]), sqa.creOp(alphaIndices[3]), \
                             sqa.desOp(alphaIndices[5]), sqa.desOp(alphaIndices[4])])
h2bbbb = sqa.term(0.25, [], [sqa.tensor('h2bbbb', betaIndices[2:6], d2sym_aaaa), \
                             sqa.creOp(betaIndices[2]), sqa.creOp(betaIndices[3]), \
                             sqa.desOp(betaIndices[5]), sqa.desOp(betaIndices[4])])
h2baba = sqa.term(1.00, [], [sqa.tensor('h2baba', [betaIndices[6],alphaIndices[6],betaIndices[7],alphaIndices[7]], d2sym_baba), \
                             sqa.creOp(betaIndices[6]), sqa.creOp(alphaIndices[6]), \
                             sqa.desOp(alphaIndices[7]),  sqa.desOp(betaIndices[7])])
h = [h1aa,h1bb,h2aaaa,h2bbbb,h2baba]

# Define Operators
r1aa = [sqa.term( 1.0, [], [sqa.tensor('r1aa', alphaIndices[8:10], a1sym), sqa.creOp(alphaIndices[8]), sqa.desOp(alphaIndices[9])]), \
        sqa.term(-1.0, [], [sqa.tensor('r1aa', alphaIndices[8:10], a1sym), sqa.creOp(alphaIndices[9]), sqa.desOp(alphaIndices[8])])]
r1bb = [sqa.term( 1.0, [], [sqa.tensor('r1bb', betaIndices[8:10], a1sym), sqa.creOp(betaIndices[8]), sqa.desOp(betaIndices[9])]), \
        sqa.term(-1.0, [], [sqa.tensor('r1bb', betaIndices[8:10], a1sym), sqa.creOp(betaIndices[9]), sqa.desOp(betaIndices[8])])]
r2aaaa = [sqa.term( 1.00, [], [sqa.tensor('r2aaaa', alphaIndices[10:14], a2sym_aaaa), \
                               sqa.creOp(alphaIndices[10]), sqa.creOp(alphaIndices[11]), \
                               sqa.desOp(alphaIndices[13]), sqa.desOp(alphaIndices[12])]), \
          sqa.term(-1.00, [], [sqa.tensor('r2aaaa', alphaIndices[10:14], a2sym_aaaa), \
                               sqa.creOp(alphaIndices[12]), sqa.creOp(alphaIndices[13]), \
                               sqa.desOp(alphaIndices[11]), sqa.desOp(alphaIndices[10])])]
r2bbbb = [sqa.term( 1.00, [], [sqa.tensor('r2bbbb', betaIndices[10:14], a2sym_aaaa), \
                               sqa.creOp(betaIndices[10]), sqa.creOp(betaIndices[11]), \
                               sqa.desOp(betaIndices[13]), sqa.desOp(betaIndices[12])]), \
          sqa.term(-1.00, [], [sqa.tensor('r2bbbb', betaIndices[10:14], a2sym_aaaa), \
                               sqa.creOp(betaIndices[12]), sqa.creOp(betaIndices[13]), \
                               sqa.desOp(betaIndices[11]), sqa.desOp(betaIndices[10])])]
r2baba = [sqa.term( 1.00, [], [sqa.tensor('r2baba', [betaIndices[14],alphaIndices[14],betaIndices[15],alphaIndices[15]], a2sym_baba), \
                               sqa.creOp(betaIndices[14]), sqa.creOp(alphaIndices[14]), \
                               sqa.desOp(alphaIndices[15]),  sqa.desOp(betaIndices[15])]), \
          sqa.term(-1.00, [], [sqa.tensor('r2baba', [betaIndices[14],alphaIndices[14],betaIndices[15],alphaIndices[15]], a2sym_baba), \
                               sqa.creOp(betaIndices[15]),  sqa.creOp(alphaIndices[15]), \
                               sqa.desOp(alphaIndices[14]), sqa.desOp(betaIndices[14])])]
r = r1aa + r1bb + r2aaaa + r2bbbb + r2baba

# Compute commutator
result = sqa.commutator(h,r)
print 'commutator complete at %.3f' %(time.time() - startTime)
print ''

# Combine any terms differing only by transpose of cre/des operators
sqa.combine_transpose(result)
print ''
print 'transpose combination complete at %.3f' %(time.time() - startTime)
print ''

# Convert cre/des ops to density matrices
rdm_dummy_name = 'rdm_dummy'
sqa.convert_ops_to_rdms_so(result, rdm_dummy_name)
sqa.assign_rdm_types(result, rdm_dummy_name, [d1aa, d1bb, d2aaaa, d2bbbb, d2baba, d3aaaaaa, d3bbbbbb, d3baabaa, d3abbabb])
print 'density matrix conversion complete at %.3f' %(time.time() - startTime)
print ''

# Combine any like terms
sqa.combineTerms(result)
print 'term combination complete at %.3f' %(time.time() - startTime)
print ''

# Print result
print 'result (%i terms):' %(len(result))
for t in result:
  print t
print ''

## Compile a list of index names
#index_name_list = []
#for t in result:
#  for ten in t.tensors:
#    for ind in ten.indices:
#      if not (ind.name in index_name_list):
#        index_name_list.append(ind.name)
#
## Write contractor input code
#con_code = ''
#con_code += 'header\n'
#con_code += '#include "fct.fh"\n'
#con_code += '\n'
#con_code += '!Author:  Eric Neuscamman\n'
#con_code += '\n'
#con_code += 'endheader\n'
#con_code += '\n'
#con_code += 'subroutine f_ct_res_3rdm_opalg_so\n'
#con_code += '\n'
#for ind in index_name_list:
#  con_code += 'index %s\n' %ind
#con_code += '\n'
#con_code += 'input nn       integer\n'
#con_code += 'input nc       integer\n'
#con_code += 'input na       integer\n'
#con_code += 'input nv       integer\n'
#con_code += 'input h1aa     real(kind=8)  1:nn 1:nn\n'
#con_code += 'input h1bb     real(kind=8)  1:nn 1:nn\n'
#con_code += 'input h2aaaa   real(kind=8)  1:nn 1:nn 1:nn 1:nn\n'
#con_code += 'input h2bbbb   real(kind=8)  1:nn 1:nn 1:nn 1:nn\n'
#con_code += 'input h2baba   real(kind=8)  1:nn 1:nn 1:nn 1:nn\n'
#con_code += 'input d1aa     real(kind=8)  IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d1bb     real(kind=8)  IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d2aaaa   real(kind=8)  IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d2bbbb   real(kind=8)  IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d2baba   real(kind=8)  IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d3aaaaaa real(kind=8)  IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d3bbbbbb real(kind=8)  IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d3baabaa real(kind=8)  IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'input d3abbabb real(kind=8)  IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += '\n'
#con_code += 'inout r1aa     real(kind=8)  ABEGIN:AEND IBEGIN:IEND\n'
#con_code += 'inout r1bb     real(kind=8)  ABEGIN:AEND IBEGIN:IEND\n'
#con_code += 'inout r2aaaa   real(kind=8)  ABEGIN:AEND ABEGIN:AEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'inout r2bbbb   real(kind=8)  ABEGIN:AEND ABEGIN:AEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += 'inout r2baba   real(kind=8)  ABEGIN:AEND ABEGIN:AEND IBEGIN:IEND IBEGIN:IEND\n'
#con_code += '\n'
#for t in result:
#  con_code += 'contraction\n'
#  target_string = None
#  for ten in t.tensors:
#    if ten.name[0] != 'r':
#      con_code += 'tensor %6s  ' %(ten.name)
#      for ind in ten.indices:
#       con_code += ' %2s' %(ind.name)
#      con_code += '\n'
#    else:
#      target_string = 'target %6s  ' %(ten.name)
#      for ind in ten.indices:
#        target_string += ' %2s' %(ind.name)
#      target_string += '\n'
#  if target_string is None:
#    raise RuntimeError("term does not have a target tensor")
#  con_code += target_string
#  con_code += 'alpha 1.00d+00\n'
#  c = '%13.6e' %(t.numConstant)
#  con_code += 'beta %sd%s\n' %(c[0:9],c[10:])
#  con_code += 'endcontraction\n'
#  con_code += '\n'
#
## Output the contractor input code
#print con_code
