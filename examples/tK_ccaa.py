import sqa.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

dummy = True

p = [sqa.index('p%i' %j, [], dummy) for j in range(6)]

x = sqa.index('x')
y = sqa.index('y')
w = sqa.index('w')
z = sqa.index('z')
ywz_op = [sqa.creOp(y), sqa.creOp(w), sqa.desOp(z)]
x_op = [sqa.desOp(x)]

h_symmetry = sqa.symmetry((1,0),1)
h_t = [sqa.tensor('h', [p[0], p[1]], [h_symmetry])]
h_op = [sqa.creOp(p[0]), sqa.desOp(p[1])]

v_symmetry = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2), -1)]
v_t = [sqa.tensor('v', [p[2], p[3], p[4], p[5]], v_symmetry)]
v_op = [sqa.creOp(p[2]), sqa.creOp(p[3]), sqa.desOp(p[5]), sqa.desOp(p[4])]

terms = []
terms += [sqa.term(1.0, [], h_t + ywz_op + h_op + x_op)]
terms += [sqa.term(-1.0, [], h_t + h_op + ywz_op + x_op)]
terms += [sqa.term(0.25, [], v_t + ywz_op + v_op + x_op)]
terms += [sqa.term(-0.25, [], v_t + v_op + ywz_op + x_op)]

noTerms = []

for t in terms:
    t_no = sqa.normalOrder(t)
    noTerms.extend(t_no)

for t in noTerms:
    t.contractDeltaFuncs()

sqa.termChop(noTerms)
sqa.combineTerms(noTerms)

print "Final results:"
for t in noTerms:
    print t
