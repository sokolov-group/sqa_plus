#    file:  sqaHeff.py
#  author:  Koushik Chatterjee
#    date:  September 21, 2018
# summary:  Heff : Construct effective Hamiltonian(L) ..
#
# (c) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#
#import sqa_extra.secondQuantizationAlgebra as sqa

from sqaIndex import index
from sqaCommutator import commutator
from sqaTerm import term
from sqaTensor import tensor, creOp, desOp
from sqaOptions import options
from sqaSymmetry import symmetry

#####################################
#
def Heff(order):
 "Construct effective Hamiltonian(L)."
#   order = 0 : L(0) = H(0)
#   order = 1 : L(1) = V + [H(0),T(1) - T'(1)]
#   order = 2 : L(2) = [H(0),(T(2) - T'(2))] + 1/2[(V + L(1)), (T(1) - T'(1))] 
#
 print "Construct effective Hamiltonian:=>"
#
 print ''
#
# Define operator types
 tg_c = options.core_type
 tg_a = options.active_type
 tg_v = options.virtual_type
 tg_g = tg_c + tg_a + tg_v
 dummy = True
# Core dummy indices
 cc = [index('c%i' %p, [tg_c], dummy) for p in range(10)]
# Active dummy indices
 aa = [index('a%i' %p, [tg_a], dummy) for p in range(10)]
# Virtual dummy indices
 vv = [index('v%i' %p, [tg_v], dummy) for p in range(10)]
#
 effH = []
 Hamil = []
# E_fc :
 Hamil.append( term(1.0, ['E_fc'], []))
#
 cor = cc.pop(0)
 vir = vv.pop(0) 
# core and vitual part : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
 Hamil.append( term(-1.0, ['e_i'],[ desOp(cor), creOp(cor)]))
 Hamil.append( term(1.0, ['e_v'],[ creOp(vir), desOp(vir)])) 
# active part : H_act
 Hact = []
 act1 = aa.pop(0)
 act2 = aa.pop(0)
 cor = cc.pop(0)
 act3 = aa.pop(0)
 act4 = aa.pop(0)
# symmetry
 h1sym = [ symmetry((1,0),1)]
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
#
 h1 =  tensor('h',[act2, act1], h1sym)
 v1 =  tensor('v', [act2, cor, act1, cor], v2sym)
 v2 =  tensor('v', [act3, act4, act1, act2], v2sym)
#
 Hact.append( term(1.0, [], [h1,  creOp(act1), desOp(act2)]))
 Hact.append( term(1.0, [], [v1,  creOp(act1), desOp(act2)]))
 Hact.append( term(1.0, ['1/4'], [v2,  creOp(act1), creOp(act2), desOp(act4), desOp(act3)]))
#
# for t in Hact:
#    print 'Hact=', t
#
 Hamil.extend(Hact)
# for t in Hamil:
#    print 'Hamiltonian(0)=', t
#
 if (order == 0):
    effH.extend(Hamil)
#
 elif (order == 1):
# L(1) = V + [H(0),T(1) - T'(1)]
    V = []
    T = []
    effH.extend(pertbV(V, 'V[n=0]', cc, aa, vv))
#
    T.extend(ampT(T,'C-A', cc, aa, vv))
#
    com1 =  commutator(Hamil, T)
#
    effH.extend(com1)
#
# for t in effH:
#    print 'Effec. Hamil=', t
#
 return effH
#
#####################################
#
def pertbV(V, vtype, cc, aa, vv):
#
 "Construct perturbation operator V according to str(vtype)"
#
 V = []
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
 d1sym = [ symmetry((1,0),1)]
#
 if (vtype == 'V[n=0]'):
#
    cor1 = cc.pop(0)
    cor2 = cc.pop(0)
    cor3 = cc.pop(0)
    cor4 = cc.pop(0)
#
    act1 = aa.pop(0)
    act2 = aa.pop(0)
#
    vir1 = vv.pop(0)
    vir2 = vv.pop(0)
    vir3 = vv.pop(0)
    vir4 = vv.pop(0)
#
    ten1 =  tensor('V', [cor3, cor4, cor1, cor2], v2sym)
#    V.append( term(1.0, ['1/4'], [ten1, creOp(cor1), creOp(cor2), desOp(cor4), desOp(cor3)]))
    V.append( term(1.0, ['1/4'], [ten1, desOp(cor4), desOp(cor3),creOp(cor1), creOp(cor2)]))
#
    ten2 =  tensor('V', [vir3, vir4, vir1, vir2], v2sym)
    V.append( term(1.0, ['1/4'], [ten2, creOp(vir1), creOp(vir2), desOp(vir4), desOp(vir3)]))
#
    ten3 =  tensor('V', [vir1, cor2, cor1, vir2], v2sym)
#    V.append( term(1.0, [], [ten3, creOp(cor1), creOp(vir2), desOp(cor2), desOp(vir1)]))
    V.append( term(1.0, [], [ten3, desOp(cor2), creOp(cor1), creOp(vir2), desOp(vir1)]))
#
    ten4 =  tensor('V', [cor2, act2, cor1, act1], v2sym)
#    V.append( term(1.0, [], [ten4, creOp(cor1), desOp(cor2), creOp(act1), desOp(act2)]))
    V.append( term(-1.0, [], [ten4, desOp(cor2), creOp(cor1), creOp(act1), desOp(act2)]))
#
    ten5 = ten4
    ten6 =  tensor('gma', [act2, act1], d1sym)
#    V.append( term(-1.0, [], [ten5, ten6,  creOp(cor1), desOp(cor2)]))
    V.append( term(1.0, [], [ten5, ten6,  desOp(cor2), creOp(cor1)]))
#
    ten7 =  tensor('V', [vir2, act2, vir1, act1], v2sym)
    V.append( term(1.0, [], [ten7, creOp(vir1), desOp(vir2), creOp(act1), desOp(act2)]))
    ten8 = ten7
    ten9 = ten6
    V.append( term(-1.0, [], [ten8, ten9,  creOp(vir1), desOp(vir2)]))
#
# for t in V:
#    print 'perterbative =', t
 return V
#####################################
#
def ampT(T,ttype,cc, aa, vv):
# Cluster operator : T - T^dag, Where T = T1 + T2
# Single excitatio : T1
 t1_sym = [ symmetry((1,0),1)]
 t2_sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
#
 T = []
 Tex = []
 Tdex = []
 if (ttype == 'C-E'):
# Core-External
     ind1 = cc.pop(0)
     ind2 = vv.pop(0)
     ind3 = aa.pop(0)
     ind4 = aa.pop(0)
 elif (ttype == 'C-A'):
# Core-Active
     ind1 = cc.pop(0)
     ind2 = aa.pop(0)
     ind3 = aa.pop(0)
     ind4 = aa.pop(0)
 elif (ttype == 'A-E'):
# Active-External
     ind1 = aa.pop(0)
     ind2 = vv.pop(0)
     ind3 = aa.pop(0)
     ind4 = aa.pop(0)
#
 t1_tens =  tensor('t', [ind2,ind1],t1_sym)
 t2_tens =  tensor('t', [ind2,ind3,ind1,ind4],t2_sym)
 T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
 T2_ex =  term(1.0, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])
#
 T1_dex =  term(-1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])
 T2_dex =  term(-1.0, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])
#
 Tex.append(T1_ex)
 Tex.append(T2_ex)
#
 Tdex.append(T1_dex)
 Tdex.append(T2_dex)
# T = T-T^dag
 T.append(T1_ex)
 T.append(T2_ex)
 T.append(T1_dex)
 T.append(T2_dex)
# 
# for t in T:
#     print 't ampli=', t
 return T
#############

