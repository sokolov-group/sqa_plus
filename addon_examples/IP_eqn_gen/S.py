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
a = sqa.index('A', [tg_v], False)
b = sqa.index('B', [tg_v], False)
x = sqa.index('X', [tg_a], False)
y = sqa.index('Y', [tg_a], False)

c = sqa.index('C', [tg_v], False)
d = sqa.index('D', [tg_v], False)
z = sqa.index('Z', [tg_a], False)
w = sqa.index('W', [tg_a], False)

i = sqa.index('I', [tg_c], False)
j = sqa.index('J', [tg_c], False)
k = sqa.index('K', [tg_c], False)
l = sqa.index('L', [tg_c], False)

#
effH = []
effH = sqa.Heff(0)
for t in effH:
  print t
#
term1 = sqa.term(1.0, [], [sqa.creOp(a), sqa.desOp(j), sqa.desOp(i)])
term2 = sqa.term(1.0, [], [sqa.creOp(k), sqa.creOp(l), sqa.desOp(b)])

print term1, term2

# for overlap S[]
termS = []
termcom0 = sqa.multiplyTerms(term1,term2)
termS.append(termcom0)

termcom1 = sqa.multiplyTerms(term2,term1)
termS.append(termcom1)
print termS
term7 = []
for t in termS:
  tt = sqa.normalOrder(t)
  term7.extend(tt)

print ''
for t in term7:
  print t
#exit()
#
term8 = sqa.matrixBlock(term7)
#exit()
#
print 'terms= ',term1, term2
print "overlap S[+1]"
for t in term8:
   print t
#
#sqa.generateEinsum(term8)
#sqa.generateEinsum(term8, M_str = 'S_p2')
sqa.generateEinsum(term8, 'S_p2', 'YXWZ', transRDM = True, trans_ind_str = 'll')
####################################################################
