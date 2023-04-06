import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

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

i_op = [sqa.desOp(i)]
j_op = [sqa.creOp(j)]

effH = []
effH = sqa.Heff(2)
for t in effH:
  print t

term5 = sqa.matrixBlock(effH)

sqa.generateEinsum(term5, 'temp', '')
#
#
#
##################################################################
