import sqa.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

dummy = True

p = [sqa.index('p%i' %j, [], dummy) for j in range(6)]

x = sqa.index('x')
y = sqa.index('y')
z = sqa.index('z')
w = sqa.index('w')
x_op = [sqa.creOp(x)]
zwy_op = [sqa.creOp(z), sqa.desOp(w), sqa.desOp(y)]

h_symmetry = sqa.symmetry((1,0),1)
h_t = [sqa.tensor('h', [p[0], p[1]], [h_symmetry])]
h_op = [sqa.creOp(p[0]), sqa.desOp(p[1])]

v_symmetry = [sqa.symmetry((1,0,2,3),-1), sqa.symmetry((0,1,3,2), -1)]
v_t = [sqa.tensor('v', [p[2], p[3], p[4], p[5]], v_symmetry)]
v_op = [sqa.creOp(p[2]), sqa.creOp(p[3]), sqa.desOp(p[5]), sqa.desOp(p[4])]

terms = []
terms += [sqa.term(1.0, [], h_t + h_op + zwy_op + x_op)]
terms += [sqa.term(-1.0, [], h_t + zwy_op + h_op + x_op)]
terms += [sqa.term(0.25, [], v_t + v_op + zwy_op + x_op)]
terms += [sqa.term(-0.25, [], v_t + zwy_op + v_op + x_op)]
terms += [sqa.term(1.0, [], h_t + x_op + h_op + zwy_op)]
terms += [sqa.term(-1.0, [], h_t + x_op + zwy_op + h_op)]
terms += [sqa.term(0.25, [], v_t + x_op + v_op + zwy_op)]
terms += [sqa.term(-0.25, [], v_t + x_op + zwy_op + v_op)]

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

