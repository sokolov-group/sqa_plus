import sqa_extra.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

# Define operator types
tg_c = sqa.options.core_type
tg_a = sqa.options.active_type
tg_v = sqa.options.virtual_type
tg_g = tg_c + tg_a + tg_v

# Define indices
dummy = True

effH = []
effH = sqa.Heff(2)

term1 = sqa.matrixBlock(effH)

sqa.generateEinsum(term1, '        temp', '', transRDM = True, trans_ind_str = 'P')
