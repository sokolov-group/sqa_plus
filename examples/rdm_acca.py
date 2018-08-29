import sqa.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

dummy = True

p = [sqa.index('p%i' %j, [], dummy) for j in range(6)]

x = sqa.index('x')
y = sqa.index('y')
w = sqa.index('w')
z = sqa.index('z')
x_op = sqa.desOp(x)
z_op = sqa.creOp(z)
w_op = sqa.creOp(w)
y_op = sqa.desOp(y)

terms = [sqa.term(1.0, [], [x_op, z_op, w_op, y_op])]

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

