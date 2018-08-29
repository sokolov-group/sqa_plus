import sqa.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

dummy = True

p = [sqa.index('p%i' %j, [], dummy) for j in range(6)]

x = sqa.index('x')
y = sqa.index('y')
w = sqa.index('w')
z = sqa.index('z')
x_op = sqa.creOp(x)
y_op = sqa.desOp(y)
w_op = sqa.desOp(w)
z_op = sqa.creOp(z)

terms = [sqa.term(1.0, [], [z_op, w_op, y_op, x_op])]

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

