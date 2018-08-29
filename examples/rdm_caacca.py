import sqa.secondQuantizationAlgebra as sqa

sqa.options.verbose = True

dummy = True

p = [sqa.index('p%i' %j, [], dummy) for j in range(6)]

z = sqa.index('z')
w = sqa.index('w')
x = sqa.index('x')
y = sqa.index('y')
wp = sqa.index('wp')
zp = sqa.index('zp')
z_op = sqa.creOp(z)
w_op = sqa.desOp(w)
x_op = sqa.desOp(x)
y_op = sqa.creOp(y)
wp_op = sqa.creOp(wp)
zp_op = sqa.desOp(zp)

terms = [sqa.term(1.0, [], [z_op, w_op, x_op, y_op, wp_op, zp_op])]

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

