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
i = sqa.index('I', [tg_c], False)
j = sqa.index('J', [tg_c], False)
x = sqa.index('X', [tg_a], False)
y = sqa.index('Y', [tg_a], False)

k = sqa.index('K', [tg_c], False)
l = sqa.index('L', [tg_c], False)
z = sqa.index('Z', [tg_a], False)
w = sqa.index('W', [tg_a], False)

#
effH = []
effH = sqa.Heff(0)
for t in effH:
  print t
#
term1 = sqa.term(1.0, [], [sqa.creOp(k), sqa.creOp(l), sqa.desOp(z), sqa.desOp(w)])
term2 = sqa.term(1.0, [], [sqa.creOp(x), sqa.creOp(y), sqa.desOp(j), sqa.desOp(i)])

print term1, term2

#
print "First Commutator"
term3 = sqa.commutator(effH,term2)
#
print "Multiply ..."
term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(term1,t))
#
term5 = []
for t in term4:
   t_no = sqa.normalOrder(t)
   term5.extend(t_no)
#
for t in term5:
 print t
#
term6 = sqa.matrixBlock(term5)

#sqa.generateEinsum(term6, 'K[:,:,:,:,:,:,:,:]', 'IJAXKLBY', "")
sqa.generateEinsum(term6, 'K', 'YXWZ', '')
#
#exit()
####################################################################
# for overlap S[+1]
termS = sqa.multiplyTerms(term1,term2)
term7 = []
tt = sqa.normalOrder(termS)
term7.extend(tt)

print ''
for t in term7:
  print t
#
term8 = sqa.matrixBlock(term7)
#
print "overlap S[+1]"
for t in term8:
   print t
#
sqa.generateEinsum(term8, 'S_p2', 'YXWZ', '')
#
exit()
####################################################################
print "V"
#
dummy = True
cc = [sqa.index('c%i' %p, [tg_c], dummy) for p in range(50)]
aa = [sqa.index('a%i' %p, [tg_a], dummy) for p in range(50)]
vv = [sqa.index('v%i' %p, [tg_v], dummy) for p in range(50)]
#
V = []
#
V = sqa.Vperturbation_type(cc, aa, vv)
#
termV = []
for t in V:
  tV = sqa.multiplyTerms(term1,t)
  tVnorm = sqa.normalOrder(tV)
  termV.extend(tVnorm)
#
term9 = sqa.matrixBlock(termV)
#for t in term9:
#  print t
sqa.generateEinsum(term9, 'Vp2', 'KLWZ', "")

####################################################################
