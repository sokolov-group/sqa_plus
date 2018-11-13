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
a = sqa.index('A', [tg_v], False)
b = sqa.index('B', [tg_v], False)

# Dummy indices
k = sqa.index('K', [tg_c], False)
l = sqa.index('L', [tg_c], False)
c = sqa.index('C', [tg_v], False)
d = sqa.index('D', [tg_v], False)

# symmetry
h1sym = [sqa.symmetry((1,0),1)]
v2sym = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2), -1)]
#
#
effH = []
effH = sqa.Heff(0)
for t in effH:
  print t
#
term1 = sqa.term(1.0, [], [sqa.creOp(i), sqa.creOp(j), sqa.desOp(b), sqa.desOp(a)])
term2 = sqa.term(1.0, [], [sqa.creOp(c), sqa.creOp(d), sqa.desOp(l), sqa.desOp(k)])
print term1, term2
#
print "First Commutator"
term3 = sqa.commutator(effH,term2)
for t in term3:
    print t
#
print "Multiply ..."
term4 = []
for t in term3:
    term4.append(sqa.multiplyTerms(term1,t))
#
term5 = []
for t in term4:
#
   t_no = sqa.normalOrder(t)
   term5.extend(t_no)
#
#for t in term4:
# print t

term6 = sqa.matrixBlock(term5)

sqa.generateEinsum(term6, 'M[sia:fia, sia:fia]', 'IAJB', ".reshape(nia_so, nia_so).copy()")
#
exit()
##################################################################
# for overlap S[+1]
termS = sqa.multiplyTerms(term1,term2)
term7 = []
tt = sqa.normalOrder(termS)
term7.extend(tt)
#
term8 = sqa.matrixBlock(term7)
#
print "overlap S[0']"
for t in term8:
   print t
#
exit()
##################################################################
print "V"
#
dummy = True
cc = [sqa.index('c%i' %p, [tg_c], False) for p in range(50)]
aa = [sqa.index('a%i' %p, [tg_a], False) for p in range(50)]
vv = [sqa.index('v%i' %p, [tg_v], False) for p in range(50)]
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
##################################################################
