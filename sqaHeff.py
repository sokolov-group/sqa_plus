# file:  sqaHeff.py
# author:  Koushik Chatterjee
# date:  September 28, 2018
# summary:  
#           Heff : Construct effective Hamiltonian(L^N) of order 'N'.
#
# Copyright (C) 2018-2020 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
#
#import sqa_extra.secondQuantizationAlgebra as sqa
import sys, os, time
import subprocess
import collections

from sqaIndex import index
from sqaCommutator import commutator
from sqaTerm import term, combineTerms
from sqaTensor import tensor, creOp, desOp, kroneckerDelta
from sqaOptions import options
from sqaSymmetry import symmetry

#####################################
#
def Heff(order, internal_excit = True):
# print_header()
 "Construct effective Hamiltonian(L)."
#
 print ("------------------------ Hamiltonian(%s) ----------------------" % order)
 sys.stdout.flush()
#   order = 0 : L(0) = H(0)
#   order = 1 : L(1) = V + [H(0),T(1) - T'(1)]
#   order = 2 : L(2) = [H(0),(T(2) - T'(2))] + 1/2[(V + L(1)), (T(1) - T'(1))] 
#
# Define operator types
 tg_c = options.core_type
 tg_a = options.active_type
 tg_v = options.virtual_type
 tg_g = tg_c + tg_a + tg_v
 dummy = True
# Core dummy indices
 cc = [index('c%i' %p, [tg_c], dummy) for p in range(800)]
# Active dummy indices
 aa = [index('a%i' %p, [tg_a], dummy) for p in range(800)]
# Virtual dummy indices
 vv = [index('v%i' %p, [tg_v], dummy) for p in range(800)]
#
 if (order == 0):
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(10):
       cc1.append(cc.pop(0))
       aa1.append(aa.pop(0))
       vv1.append(vv.pop(0))
#
    L = dyallH(cc1, aa1, vv1)
#
 elif (order == 1):   # L(1) = V + [H(0),T(1) - T'(1)]
    L = []

    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(10):
       cc1.append(cc.pop(0))
       aa1.append(aa.pop(0))
       vv1.append(vv.pop(0))
#
    effH = dyallH(cc1, aa1, vv1)
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(200):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
#
    V = Vperturbation_type(cc1, aa1, vv1)
#
    L.extend(V)
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1 = Tamplitude(1, cc1, aa1, vv1)
#
    com1 = commutator(effH, T1)
    print "Commutation: Done ..."
    sys.stdout.flush()
#
    L.extend(com1)
#
 elif (order == 2):    #   L(2) = [H(0),T(2) - T'(2)]+ 1/2 [V + L(1),T(1) - T'(1)]
    L = []
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(10):
       cc1.append(cc.pop(0))
       aa1.append(aa.pop(0))
       vv1.append(vv.pop(0))
#
    effH = dyallH(cc1, aa1, vv1)
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T2 = Tamplitude(2, cc1, aa1, vv1, internal_excit)
    com1 = commutator(effH, T2)
    print "First Commutation: Done ..."
    sys.stdout.flush()
#
    L.extend(com1)              # 1st Commutator
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(200):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
#
    V = Vperturbation_type(cc1, aa1, vv1)
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1 = Tamplitude(1, cc1, aa1, vv1)
    com2 = commutator(V, T1)
    print "Second Commutation: Done ..."
    sys.stdout.flush()
#
    L.extend(com2)              # 2nd Commutator
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(10):
       cc1.append(cc.pop(0))
       aa1.append(aa.pop(0))
       vv1.append(vv.pop(0))
#
    effH = dyallH(cc1, aa1, vv1)
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1_new1 = Tamplitude(1, cc1, aa1, vv1)
#
    com3 = commutator(effH, T1_new1)
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
    T1_new2 = Tamplitude(1, cc1, aa1, vv1)
#
    com4 = commutator(com3, T1_new2)
    print "Third Commutation: Done ..."
    sys.stdout.flush()
#
    for t in com4:
       t.scale(0.5)
    L.extend(com4)              # 3rd Commutator
#
 else:
    raise Exception('Unknown type of effective Hamiltonian of order = %s' % (order))
#
 print "Done ..."
 print("""--------------------------------------------------------------""")
 sys.stdout.flush()
 return L
#
def dyallH(cc, aa, vv):
 Hamil = []
# E_fc :
 #c = index('Const.', [], dummy = True)
# c = index('Const.', [], True)
# Efc = tensor('E_fc',[c], [])
 Efc = tensor('E_fc', [], [])
# Hamil.append( term(1.0, ['E_fc'], []))
 Hamil.append( term(1.0, [], [Efc]))
#
# core and vitual part : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
 e_core = tensor('e', [cc[1]], [])
 e_virt = tensor('e', [vv[1]], [])
# Hamil.append( term(-1.0, ['e_i'],[ desOp(cor), creOp(cor)]))
# Hamil.append( term(1.0, ['e_v'],[ creOp(vir), desOp(vir)])) 
 Hamil.append( term(-1.0, [],[e_core, desOp(cc[1]), creOp(cc[1])]))
 Hamil.append( term(1.0, [],[e_virt, creOp(vv[1]), desOp(vv[1])]))
#
# active part : H_act
 Hact = []
# symmetry
 h1sym = [ symmetry((1,0),1)]
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
#
 h1 =  tensor('h',[aa[1], aa[0]], h1sym)
 v1 =  tensor('v', [aa[3], cc[2], aa[2], cc[2]], v2sym)
 v2 =  tensor('v', [aa[6], aa[7], aa[4], aa[5]], v2sym)
#
 Hact.append( term(1.0, [], [h1,  creOp(aa[0]), desOp(aa[1])]))
 Hact.append( term(1.0, [], [v1,  creOp(aa[3]), desOp(aa[2])]))
 Hact.append( term(0.25, [], [v2,  creOp(aa[4]), creOp(aa[5]), desOp(aa[7]), desOp(aa[6])]))
#
 Hamil.extend(Hact)
#
 return Hamil

def Vperturbation_type(cc, aa, vv, vtype = None):
#
 "Construct perturbation operator V according to excitation rank."
#
 V = []
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
 d1sym = [ symmetry((1,0),1)]
#
# if (vtype == 'full') or (vtype == 'Full') or (vtype == ''):
 if not (vtype):
#   Default V includes all type of perturbation rank.
    print "Perturbation(V) type = All "
#    print ""
    V.extend(Vperturbation(cc, aa, vv))
    print "Done ..."
    sys.stdout.flush()
#
 else:
    if (vtype == 'V[n=0]'):
#   
       print "Perturbation(V) type = ",vtype
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
       ten1 =  tensor('v', [cor3, cor4, cor1, cor2], v2sym)
#       V.append( term(0.25, [], [ten1, creOp(cor1), creOp(cor2), desOp(cor4), desOp(cor3)]))
       V.append( term(0.25, [], [ten1, desOp(cor4), desOp(cor3),creOp(cor1), creOp(cor2)]))
#   
       ten2 =  tensor('v', [vir3, vir4, vir1, vir2], v2sym)
       V.append( term(0.25, [], [ten2, creOp(vir1), creOp(vir2), desOp(vir4), desOp(vir3)]))
#   
       ten3 =  tensor('v', [vir1, cor2, cor1, vir2], v2sym)
#       V.append( term(1.0, [], [ten3, creOp(cor1), creOp(vir2), desOp(cor2), desOp(vir1)]))
       V.append( term(1.0, [], [ten3, desOp(cor2), creOp(cor1), creOp(vir2), desOp(vir1)]))
#   
       ten4 =  tensor('v', [cor2, act2, cor1, act1], v2sym)
#       V.append( term(1.0, [], [ten4, creOp(cor1), desOp(cor2), creOp(act1), desOp(act2)]))
       V.append( term(-1.0, [], [ten4, desOp(cor2), creOp(cor1), creOp(act1), desOp(act2)]))
#   
       ten5 = ten4
       ten6 =  tensor('gamma', [act2, act1], d1sym)
#       V.append( term(-1.0, [], [ten5, ten6,  creOp(cor1), desOp(cor2)]))
       V.append( term(1.0, [], [ten5, ten6,  desOp(cor2), creOp(cor1)]))
#   
       ten7 =  tensor('v', [vir2, act2, vir1, act1], v2sym)
       V.append( term(1.0, [], [ten7, creOp(vir1), desOp(vir2), creOp(act1), desOp(act2)]))
       ten8 = ten7
       ten9 = ten6
       V.append( term(-1.0, [], [ten8, ten9,  creOp(vir1), desOp(vir2)]))
#   
    else:
       raise Exception('Unknown type of V operator ...')
#
# for t in V:
#    print 'perterbative =', t
 return V
#####################################
#
def Tamplitude(order, cc1, aa1, vv1, internal_excit = True):
# Cluster operator : T - T^dag, Where T = T1 + T2
# Single excitatio : T1
#
# Define t amplitude according to their order
 if (order == 1):
     tname = 't1'
 elif (order == 2):
     tname = 't2'

 t1_sym = [ symmetry((1,0),1)]

 T = []
 cc = list(cc1)
 aa = list(aa1)
 vv = list(vv1)
# Core-External

 ind1 = cc.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens =  tensor(tname, [ind1,ind2],t1_sym)
 T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
 T1_dex =  term(-1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])

 ind1 = cc.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_ex =  term(1.0, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])
 T2_dex =  term(-1.0, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])

 T.append(T1_ex)
 T.append(T2_ex)
 T.append(T1_dex)
 T.append(T2_dex)

# Core-Active
 ind1 = cc.pop(0)
 ind2 = aa.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens = tensor(tname, [ind1,ind2],t1_sym)
 T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
 T1_dex =  term(-1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])

 ind1 = cc.pop(0)
 ind2 = aa.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_ex =  term(0.5, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])
 T2_dex =  term(-0.5, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])

 T.append(T1_ex)
 T.append(T2_ex)
 T.append(T1_dex)
 T.append(T2_dex)

# Active-External
 ind1 = aa.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens =  tensor(tname, [ind1,ind2],t1_sym)
 T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
 T1_dex =  term(-1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])

 ind1 = aa.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_ex =  term(0.5, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])
 T2_dex =  term(-0.5, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])

 T.append(T1_ex)
 T.append(T2_ex)
 T.append(T1_dex)
 T.append(T2_dex)

 # Doubles
 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens1 = custom_tensor(tname, ind1,ind2,ind5,ind6)
 T2_ex = term(0.25, [], [t2_tens1,  creOp(ind5), creOp(ind6), desOp(ind2), desOp(ind1)])
 T2_dex = term(-0.25, [], [t2_tens1,  creOp(ind1), creOp(ind2), desOp(ind6), desOp(ind5)])
 T.append(T2_ex)
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens2 = custom_tensor(tname, ind1,ind2,ind5,ind3)
 T2_ex = term(0.5, [], [t2_tens2,  creOp(ind5), creOp(ind3), desOp(ind2), desOp(ind1)])
 T2_dex = term(-0.5, [], [t2_tens2,  creOp(ind1), creOp(ind2), desOp(ind3), desOp(ind5)])
 T.append(T2_ex)
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens3 = custom_tensor(tname, ind1,ind3,ind5,ind6)
 T2_ex = term(0.5, [], [t2_tens3,  creOp(ind5), creOp(ind6), desOp(ind3), desOp(ind1)])
 T2_dex = term(-0.5, [], [t2_tens3,  creOp(ind1), creOp(ind3), desOp(ind6), desOp(ind5)])
 T.append(T2_ex)
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens4 = custom_tensor(tname, ind1,ind2,ind3,ind4)
 T2_ex = term(0.25, [], [t2_tens4,  creOp(ind3), creOp(ind4), desOp(ind2), desOp(ind1)])
 T2_dex = term(-0.25, [], [t2_tens4,  creOp(ind1), creOp(ind2), desOp(ind4), desOp(ind3)])
 T.append(T2_ex)
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens5 = custom_tensor(tname, ind4,ind3,ind5,ind6)
 T2_ex = term(0.25, [], [t2_tens5,  creOp(ind5), creOp(ind6), desOp(ind3), desOp(ind4)])
 T2_dex = term(-0.25, [], [t2_tens5,  creOp(ind4), creOp(ind3), desOp(ind6), desOp(ind5)])
 T.append(T2_ex)
 T.append(T2_dex)

 if (order == 2 and internal_excit):
     t1_asym = [ symmetry((1,0), -1)]
     ind1 = aa.pop(0)
     ind2 = aa.pop(0)
     t1_tens =  tensor(tname, [ind1,ind2],t1_asym)
     T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
     T.append(T1_ex)

 return T


def Tamplitude_excitation(order, cc1, aa1, vv1):
# Cluster operator : T - T^dag, Where T = T1 + T2
# Single excitatio : T1
#
# Define t amplitude according to their order
 if (order == 1):
     tname = 't1'
 elif (order == 2):
     tname = 't2'

 t1_sym = [ symmetry((1,0),1)]

 T = []
 cc = list(cc1)
 aa = list(aa1)
 vv = list(vv1)
# Core-External

 ind1 = cc.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens =  tensor(tname, [ind1,ind2],t1_sym)
 T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])

 ind1 = cc.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_ex =  term(1.0, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])

 T.append(T1_ex)
 T.append(T2_ex)

# Core-Active
 ind1 = cc.pop(0)
 ind2 = aa.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens = tensor(tname, [ind1,ind2],t1_sym)
 T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])

 ind1 = cc.pop(0)
 ind2 = aa.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_ex =  term(0.5, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])

 T.append(T1_ex)
 T.append(T2_ex)

# Active-External
 ind1 = aa.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens =  tensor(tname, [ind1,ind2],t1_sym)
 T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])

 ind1 = aa.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_ex =  term(0.5, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])

 T.append(T1_ex)
 T.append(T2_ex)

 # Doubles
 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens1 = custom_tensor(tname, ind1,ind2,ind5,ind6)
 T2_ex = term(0.25, [], [t2_tens1,  creOp(ind5), creOp(ind6), desOp(ind2), desOp(ind1)])
 T.append(T2_ex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens2 = custom_tensor(tname, ind1,ind2,ind5,ind3)
 T2_ex = term(0.5, [], [t2_tens2,  creOp(ind5), creOp(ind3), desOp(ind2), desOp(ind1)])
 T.append(T2_ex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens3 = custom_tensor(tname, ind1,ind3,ind5,ind6)
 T2_ex = term(0.5, [], [t2_tens3,  creOp(ind5), creOp(ind6), desOp(ind3), desOp(ind1)])
 T.append(T2_ex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens4 = custom_tensor(tname, ind1,ind2,ind3,ind4)
 T2_ex = term(0.25, [], [t2_tens4,  creOp(ind3), creOp(ind4), desOp(ind2), desOp(ind1)])
 T.append(T2_ex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens5 = custom_tensor(tname, ind4,ind3,ind5,ind6)
 T2_ex = term(0.25, [], [t2_tens5,  creOp(ind5), creOp(ind6), desOp(ind3), desOp(ind4)])
 T.append(T2_ex)

 return T


def Tamplitude_deexcitation(order, cc1, aa1, vv1, internal_excit = True):
# Cluster operator : T - T^dag, Where T = T1 + T2
# Single excitatio : T1
#
# Define t amplitude according to their order
 if (order == 1):
     tname = 't1'
 elif (order == 2):
     tname = 't2'

 t1_sym = [ symmetry((1,0),1)]

 T = []
 cc = list(cc1)
 aa = list(aa1)
 vv = list(vv1)
# Core-External

 ind1 = cc.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens =  tensor(tname, [ind1,ind2],t1_sym)
 T1_dex =  term(1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])

 ind1 = cc.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_dex =  term(1.0, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])

 T.append(T1_dex)
 T.append(T2_dex)

# Core-Active
 ind1 = cc.pop(0)
 ind2 = aa.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens = tensor(tname, [ind1,ind2],t1_sym)
 T1_dex =  term(1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])

 ind1 = cc.pop(0)
 ind2 = aa.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_dex =  term(0.5, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])

 T.append(T1_dex)
 T.append(T2_dex)

# Active-External
 ind1 = aa.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t1_tens =  tensor(tname, [ind1,ind2],t1_sym)
 T1_dex =  term(1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])

 ind1 = aa.pop(0)
 ind2 = vv.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 t2_tens = custom_tensor(tname, ind1,ind4,ind2,ind3)
 T2_dex =  term(0.5, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])

 T.append(T1_dex)
 T.append(T2_dex)

 # Doubles
 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens1 = custom_tensor(tname, ind1,ind2,ind5,ind6)
 T2_dex = term(0.25, [], [t2_tens1,  creOp(ind1), creOp(ind2), desOp(ind6), desOp(ind5)])
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens2 = custom_tensor(tname, ind1,ind2,ind5,ind3)
 T2_dex = term(0.5, [], [t2_tens2,  creOp(ind1), creOp(ind2), desOp(ind3), desOp(ind5)])
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens3 = custom_tensor(tname, ind1,ind3,ind5,ind6)
 T2_dex = term(0.5, [], [t2_tens3,  creOp(ind1), creOp(ind3), desOp(ind6), desOp(ind5)])
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens4 = custom_tensor(tname, ind1,ind2,ind3,ind4)
 T2_dex = term(0.25, [], [t2_tens4,  creOp(ind1), creOp(ind2), desOp(ind4), desOp(ind3)])
 T.append(T2_dex)

 ind1 = cc.pop(0)
 ind2 = cc.pop(0)
 ind3 = aa.pop(0)
 ind4 = aa.pop(0)
 ind5 = vv.pop(0)
 ind6 = vv.pop(0)
 t2_tens5 = custom_tensor(tname, ind4,ind3,ind5,ind6)
 T2_dex = term(0.25, [], [t2_tens5,  creOp(ind4), creOp(ind3), desOp(ind6), desOp(ind5)])
 T.append(T2_dex)

 return T


def Vperturbation(cc, aa, vv):
 from sqaAddon import matrixBlock, dummyLabel
#
 "Construct general perturbation operator V full (default) include all types of V."
#
 V = []
 V81 = []
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
 h1sym = [ symmetry((1,0),1)]
 d1sym = [ symmetry((1,0),1)]

 cc1 = list(cc)
 aa1 = list(aa)
 vv1 = list(vv)

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten1 =  tensor('h', [cor1, act1], h1sym)
 V.append( term(1.0, [], [ten1, creOp(act1), desOp(cor1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten2 =  tensor('h', [act1, cor1], h1sym)
 V.append( term(1.0, [], [ten2, creOp(cor1), desOp(act1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten3 =  tensor('h', [act1, vir1], h1sym)
 V.append( term(1.0, [], [ten3, creOp(vir1), desOp(act1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten4 =  tensor('h', [vir1, act1], h1sym)
 V.append( term(1.0, [], [ten4, creOp(act1), desOp(vir1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten5 =  tensor('h', [cor1, vir1], h1sym)
 V.append( term(1.0, [], [ten5, creOp(vir1), desOp(cor1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten6 =  tensor('h', [vir1, cor1], h1sym)
 V.append( term(1.0, [], [ten6, creOp(cor1), desOp(vir1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten7 =  tensor('v', [cor2, act2, cor1, act1], v2sym)
 ten8 =  tensor('gamma', [act2, act1], d1sym)
 V.append( term(1.0, [], [ten7, ten8,  desOp(cor2), creOp(cor1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 ten8 =  tensor('gamma', [act2, act1], d1sym)
 ten9 =  tensor('v', [vir2, act2, vir1, act1], v2sym)
# V.append( term(-1.0, [], [ten9, creOp(vir2), desOp(vir1), creOp(act2), desOp(act1)]))
 V.append( term(-1.0, [], [ten9, ten8,  creOp(vir1), desOp(vir2)]))
#
 for ityp1 in range(3):
        if (ityp1 == 0):                              # ityp1 = 0 => Core
               p = cc.pop(0)
        elif (ityp1 == 1):                            #         1 => Active
               p = aa.pop(0)
        else:                                         #         2 => Virtual
               p = vv.pop(0)
        for ityp2 in range(3):
               if (ityp2 == 0):
                      q = cc.pop(0)
               elif (ityp2 == 1):
                      q = aa.pop(0)
               else:
                      q = vv.pop(0)
#
               for ityp3 in range(3):
                      if (ityp3 == 0):
                             s = cc.pop(0)
                      elif (ityp3 == 1):
                             s = aa.pop(0)
                      else:
                             s = vv.pop(0)
#
                      for ityp4 in range(3):
                             if (ityp4 == 0):
                                   r = cc.pop(0)
                             elif (ityp4 == 1):
                                   r = aa.pop(0)
                             else:
                                   r = vv.pop(0)
#
                             if not (p.indType[0][0]=='active' and q.indType[0][0]=='active' and r.indType[0][0]=='active' and s.indType[0][0]=='active'):
#
                                   if (p.indType[0][0]=='core' and q.indType[0][0]=='core' and r.indType[0][0]=='core' and s.indType[0][0]=='core'):
#                                       vTen = tensor('v', [r,s,p,q], v2sym)
                                       vTen = tensor('v', [r,s,p,q], v2sym)
                                       V81.append(term(-0.25, [], [vTen, desOp(r), desOp(s), creOp(p), creOp(q)]))
#
                                   else:
                                       vTen = tensor('v', [r,s,p,q], v2sym)
#                                       vTen = tensor('v', [p,q,r,s], v2sym)
                                       V81.append(term(0.25, [], [vTen, creOp(p), creOp(q),desOp(s), desOp(r)]))
#
 
 p = cc.pop(0)
 q = vv.pop(0)
 s = cc.pop(0)
 r = vv.pop(0)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(-1.0, [], [vTen, creOp(p), creOp(q), desOp(s), desOp(r)]))

 p = cc.pop(0)
 q = aa.pop(0)
 s = cc.pop(0)
 r = aa.pop(0)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(-1.0, [], [vTen, creOp(p), creOp(q), desOp(s), desOp(r)]))

 p = cc.pop(0)
 q = vv.pop(0)
 s = cc.pop(0)
 r = vv.pop(0)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(1.0, [], [vTen, creOp(q), desOp(r), desOp(s), creOp(p)]))

 p = cc.pop(0)
 q = aa.pop(0)
 s = cc.pop(0)
 r = aa.pop(0)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(1.0, [], [vTen, desOp(s), creOp(p), creOp(q), desOp(r)]))

 V.extend(V81)

 return V


def print_header():

    print("""\n--------------------------------------------------------------
    SQA_extra: Construct effective Hamiltonian(L)
    author:  Koushik Chatterjee
    date:  August 31, 2018

    Copyright (C) 2018-2020  Koushik Chatterjee (koushikchatterjee7@gmail.com)

    This program is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU General Public License
    for more details.
--------------------------------------------------------------""")
#
#
def generateEinsum_old(terms, lhs_str = None, ind_str = None, tens_ext = None, transRDM = False, trans_ind_str = None, rhs_str = None, optimize = True, h_str = None, v_str = None, e_str = None, t_str = None, rdm_str = None, delta_str = None, suffix = None):
#
# summary: Generate Einsum structures for each term. 
#          terms   : A list of all terms.
#          ind_str : Indices of the matrix (string).
#
# Copyright (C) 2018-2020 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more detai
#
# print "################ Construct Einsum ################"

 print("""\n--------------------------------------------------------------
 Einsum generator: Transform into einsum ...
 author:  Koushik Chatterjee
 date:  August 31, 2018

 Copyright (C) 2018-2020  Koushik Chatterjee (koushikchatterjee7@gmail.com)
--------------------------------------------------------------""")

# print ""
#
 if not (lhs_str):
   lhs_str = 'M'
#
 term1st = 0
 for term in terms:
     outputS = []
     outputF = []
     OpsList = []
     tens_name = []
#
     icopy = '.copy()'
     if not (len(term.tensors) == 1):
        icopy = ''
#
     if not (rdm_str):
        OpsindStr = 'rdm_'
        if (transRDM):
           OpsindStr = 'trdm_'
     else:
        OpsindStr = rdm_str
#
     for i in range(len(term.tensors)):
         TensStr = ''
         tens = term.tensors[i]
         tensor_extern_name = None
         tens_name.append(tens.name)
         if not (isinstance(tens, creOp) or isinstance(tens,desOp) or isinstance(tens, kroneckerDelta)):
#
#            tens_name.append(tens.name)
#
            if ((tens.name == 'v') or (tens.name == 'V')):
               if not (v_str):
#                  indStr = 'v_'
                  indStr = str(tens.name)+'_'
               else:
                  indStr = v_str
            elif (tens.name == 'h'):
               if not (h_str):
#                  indStr = 'h_'
                  indStr = str(tens.name)+'_'
               else:
                  indStr = h_str
            elif (tens.name == 't1') or (tens.name == 't2'):
               if not (t_str):
#                  indStr = 't_'
                  indStr = str(tens.name)+'_'
               else:
                  indStr = t_str
#
            elif ((tens.name == 'E') or (tens.name == 'e')):
                  indStr = tens.name
#
            else:
               tensor_extern_name = tens.name
               if (tens_ext == None):
                  indStr = tens.name
               else:
                  indStr = tens_ext
#                  ii = 0
#                  for i in range(len(tens_name)):
#                     if (tens_name[i] == tens_ext):
#                        ii += 1
#                  if (ii > 0):
#                     indStr = tens_ext+str(ii)
#                  else: 
#                     indStr = tens_ext
#
            for tens_ind in range(len(tens.indices)): 
               TensStr += str(tens.indices[tens_ind].name)
#
               if ((tens.name == 'E') or (tens.name == 'e')):
                  if (tens.indices[tens_ind].indType[0][0] == 'core'):
                     if not (e_str):
#                        indStr = 'e_core_so'
                        indStr = str(tens.name)+'_core'
                     else:
                        indStr = e_str
#
                  elif (tens.indices[tens_ind].indType[0][0] == 'virtual'):
                     if not (e_str):
#                        indStr = 'e_extern_so'
                        indStr = str(tens.name)+'_extern'
                     else:
                        indStr = e_str
#
                  else:
                     raise Exception('Unknown active orbitals energy.')
#
#
               elif (tens.name == 'gamma'):
#                  if not (g_str):
                     indStr = 'rdm_ca'
#                     indStr = str(tens.name)+'_'
#                  else:
#                     indStr = g_str
#
               else:
                  if not (tens.name == tensor_extern_name):
                     if not (tens.indices[tens_ind].indType[0][0] == 'virtual'):
                        indStr += str(tens.indices[tens_ind].indType[0][0][0])
                     else:
                        indStr += 'e'
#
            if not ((tens.name == 't1') or (tens.name == 't2') or (tens.name == tensor_extern_name)) :
               if not (suffix):
                  indStr += '_so'
               else:
                  indStr += suffix
#
            outputF.append(indStr)   
#
            outputS.append(TensStr)
#
         elif (isinstance(tens, kroneckerDelta)):
#               tens_name.append(tens.name)
#
               for tens_ind in range(len(tens.indices)):
                   TensStr += str(tens.indices[tens_ind].name)
#
               indStr = 'np.identity'
#
               if not (delta_str):
                  if (tens.indices[0].indType[0][0] == 'core' and (tens.indices[1].indType[0][0] == 'core')):
                      delstr = 'ncore'
                  elif (tens.indices[0].indType[0][0] == 'active' and (tens.indices[1].indType[0][0] == 'active')):
                      delstr = 'ncas'
                  else:
                      delstr = 'nextern'
               else:
                  delstr = str(delta_str)
#
               if not (suffix):
                  delstr += "_so"
               else:
                  delstr += suffix
#
               indStr += str('(')+delstr+str(')')
#
               outputS.append(TensStr)
               outputF.append(indStr)
#
         else:
               OpsList.append(tens.indices[0].name)
               if (tens.name == 'cre'): 
                   OpsindStr += 'c'
               elif (tens.name == 'des'):
                   OpsindStr += 'a'
# 
     if (len(OpsList)>0):
            if (transRDM):
               if (trans_ind_str == None):
                  raise Exception("Defined 'trans_ind_str' and run again...")
               else:
                  OpsStr = trans_ind_str
            else:
               OpsStr = ''
            for i in OpsList:
                OpsStr += str(i)
            outputS.append(OpsStr)
            if not (suffix):
               OpsindStr += "_so"
            else:
               OpsindStr += suffix
            outputF.append(OpsindStr)
#
########################
#     if not (ind_str):
#        rhs_ind_str = ''
#     else:
#        if (transRDM):
#          rhs_ind_str = '->'+trans_ind_str+ind_str
#        else:
#          rhs_ind_str = "->"+ind_str
     if not (ind_str):
        if (transRDM):
          rhs_ind_str = '->'+trans_ind_str
        else:
          rhs_ind_str = ''
     else:
        if (transRDM):
          rhs_ind_str = '->'+trans_ind_str+ind_str
        else:
          rhs_ind_str = "->"+ind_str
########################
#
     sign = ''
     if not (term1st == 0):
        sign = '+'
        if (term.numConstant < 0.0):
           sign = '-'
        if not (abs(term.numConstant) == 1.0):
           cons = ' '+str(abs(term.numConstant))+' *'
        else:
           cons = ''
     else:
        cons = ' '
        if not (term.numConstant == 1.0):
           if (term.numConstant == -1.0):
               cons = '-'
           else:
               cons = ' '+str(term.numConstant)+' *'
#
     IOpt = ' optimize = True'
     if not (optimize):
        IOpt = ' optimize = False'
#
     Icomnd = ''
     if (rhs_str):
        Icomnd = rhs_str
#
##     print "M[s"+ind1+":f"+ind1+", s"+ind2+":f"+ind2+"] "+sign+"="+cons+" np.einsum('"+str(outputS).translate(None, "'")[1:-1]+"->"+ind_str+"', "+str(outputF).translate(None, "'")[1:-1]+","+IOpt+")"
#     print lhs_str+" "+sign+"="+cons+" np.einsum('"+str(outputS).translate(None, "'")[1:-1]+rhs_ind_str+"', "+str(outputF).translate(None, "'")[1:-1]+","+IOpt+")"+Icomnd+icopy
#
####
     lhs_ind_str = str(outputS).translate(None, "'")[1:-1]
     tensList = str(outputF).translate(None, "'")[1:-1]
     if (len(term.tensors)==0):
        if (transRDM):
          print lhs_str+" "+sign+"="+cons+term.constants[0]+" * "+'np.identity('+trans_ind_str+')'
        else:
          lhs_ind_str = term.constants[0][0]
          print lhs_str+" "+sign+"="+cons+term.constants[0]
     else:
        print lhs_str+" "+sign+"="+cons+" np.einsum('"+lhs_ind_str+rhs_ind_str+"', "+tensList+","+IOpt+")"+Icomnd+icopy
####
#
     term1st += 1
#
#####################################

def generateEinsum(terms, lhs_str = None, ind_str = None, transRDM = False, trans_ind_str = None, rhs_str = None, optimize = True, suffix = None, rdm_str = None, help = False, **tensor_rename):
#
# summary: Generate Einsum structures for each term. 
#          terms   : A list of all terms.
#          ind_str : Indices of the matrix (string).
#
# Copyright (C) 2018-2020 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more detai
#
# print "################ Construct Einsum ################"

 print("""\n----------------------- SQA EINSUM ---------------------------
  ___  _  _  _  __  _ _  _   _ 
 | __|| || \| |/ _|| | || \_/ | Einsum generator: Transform into einsum.
 | _| | || \\  |\_ \| U || \_/ | author:  Koushik Chatterjee
 |___||_||_|\_||__/|___||_| |_| date:  August 31, 2018
                                VERSION : 1
 Copyright (C) 2018-2020  Koushik Chatterjee (koushikchatterjee7@gmail.com)

 For Help :: help = True
--------------------------------------------------------------""")

# print ""
 if help:
   einsum_help()
#
 if not (lhs_str):
   lhs_str = 'Matrix'
#
 internal_tensor = ['e', 'E', 't1', 't2', 'h', 'v', 'cre', 'des', 'rdm', 'trdm', 'gamma']
 external_tensor = []
#
 einsum_list = []
#
 if not (suffix):
    suffix = 'so'       # suffix by default
#
 term1st = 0
 for term in terms:
     tensorlist = []
     tensor_indices_list = []
     credes_list = []
     credes_indices_list = []
#
#     name_default = {}                        # Regular Dictonary
     name_default = collections.OrderedDict()  # Ordered dictionary
#
     icopy = '.copy()'
     if not (len(term.tensors) == 1):
        icopy = ''
#
     for i in range(len(term.tensors)):
         tens = term.tensors[i]
#
         tensor_name, tensr_ind_str, tensr_indtype_str = tensor_name_indices(tens, suffix)
#
         if (isinstance(tens, creOp) or isinstance(tens,desOp)):
            credes_list.append(tensor_name)            
            credes_indices_list.append(tensr_ind_str)
         else:
            tensorlist.append((tensor_name, tensr_ind_str, tensr_indtype_str))
            tensor_indices_list.append(tensr_ind_str)
#
     if (len(credes_list) > 0) :
        cre_count = credes_list.count('cre')
        des_count = credes_list.count('des')
        if not (rdm_str):
           set_rdm = 'rdm'
           if (transRDM):
              set_rdm = 'trdm'
 #          for i in range(cre_count):
 #              set_rdm += 'c'
 #          for i in range(des_count):
 #              set_rdm += 'a'
        else:
           set_rdm = rdm_str
           internal_tensor.append(set_rdm)
#           if (transRDM):
#              set_rdm = 'trans'+rdm_str
#
        rdm_indtype_str = ''
        for i in range(cre_count):
            rdm_indtype_str += 'c'
        for i in range(des_count):
            rdm_indtype_str += 'a'
#
        rdm_ind_str = ''
        for i in credes_indices_list:
            rdm_ind_str += str(i)
#
        if (transRDM):
           if (trans_ind_str == None):
              raise Exception("Defined 'trans_ind_str' and run again...")
           else:
              rdm_ind_str = trans_ind_str+rdm_ind_str
#
        tensorlist.append((set_rdm, rdm_ind_str, rdm_indtype_str))
#
     for i in range(len(tensorlist)):
        if (tensorlist[i][0] == 'gamma'):
           key_tensr = tensorlist[i][0]
           value_tensr = tensorlist[i]
           gamma_tuple = list(tensorlist[i])
           gamma_tuple[0] = 'rdm'
           if (rdm_str):
              gamma_tuple[0] = rdm_str
           tensorlist[i] = tuple(gamma_tuple)
           value_tensr = tensorlist[i]
        else:
           key_tensr = tensorlist[i][0]
           value_tensr = tensorlist[i]
        name_default.setdefault(key_tensr, []).append(value_tensr)
#        name_default.update({key_tensr : value_tensr})
##        name_default.update({tensorlist[i][0] : tensorlist[i]})
#
        if ((tensorlist[i][0] not in internal_tensor) and (tensorlist[i][0] not in external_tensor)):
           external_tensor.append(tensorlist[i][0])

     checkey = 0
#     if (tensor_name):
     name_default_bak = dict(name_default)
     for key0, value0 in tensor_rename.items():
         if (key0 not in internal_tensor) and (key0 not in external_tensor):
            checkey += 1
            if (checkey > 0):
               raise Exception("Unknown tensor key: '%s' to rename ..." % key0)
#
         for key1, value1 in name_default.items():
            if (key0 == key1):
               for i in range(len(value1)):
                   list_tuple = list(value1[i])
                   list_tuple[0] = value0
                   external_tensor.append(value0)
                   value1[i] = tuple(list_tuple)
                  # name_default[key1] = value1
                  ## del name_default[key1]
                  ## name_default.update({value0 : value1})
#
     LHS_ind = ''
     RHS_tensr = ''
     lhs_ein_ind = []
     rhs_ein_ten = []

     for key, value in name_default.items():
         for i in range(len(value)):
             lhs_ein_ind.append(value[i][1])
#
             if (value[i][0] in external_tensor):
                ein_tens_name = value[i][0]
             else:
                if ((value[i][0] == 't1') or (value[i][0] == 't2')):
                   ein_tens_name = value[i][0]+'_'+value[i][2]
                else:
                   ein_tens_name = value[i][0]+'_'+value[i][2]+'_'+suffix
#
#             rhs_ein_ten.append(value[0]+'_'+value[2])
             rhs_ein_ten.append(ein_tens_name)


#     for i in range(len(tensorlist)):
#     #   LHS_ind += tensorlist[i][1]
#        lhs_ein_ind.append(tensorlist[i][1])
#     #   RHS_tensr += tensorlist[i][0]+tensorlist[i][2]
#        rhs_ein_ten.append(tensorlist[i][0]+'_'+tensorlist[i][2])
#
     if not (ind_str):
        if (transRDM):
          rhs_ind_str = '->'+trans_ind_str
        else:
          rhs_ind_str = ''
     else:
        if (transRDM):
          rhs_ind_str = '->'+trans_ind_str+ind_str
        else:
          rhs_ind_str = '->'+ind_str
#
     sign = ''
     if not (term1st == 0):
        sign = '+'
        if (term.numConstant < 0.0):
           sign = '-'
        if not (abs(term.numConstant) == 1.0):
           cons = ' '+str(abs(term.numConstant))+' *'
        else:
           cons = ''
     else:
        cons = ' '
        if not (term.numConstant == 1.0):
           if (term.numConstant == -1.0):
               cons = '-'
           else:
               cons = ' '+str(term.numConstant)+' *'
#
     IOpt = ' optimize = True'
     if not (optimize):
        IOpt = ' optimize = False'
#
     Icomnd = ''
     if (rhs_str):
        Icomnd = rhs_str


     lhs_ind_str = str(lhs_ein_ind).translate(None, "'")[1:-1]
     rhs_ten_str = str(rhs_ein_ten).translate(None, "'")[1:-1]

## Changes made by Rajat to change print statement to list:
     if (len(term.tensors)==0):
        if (transRDM):
          print lhs_str+" "+sign+"="+cons+term.constants[0]+" * "+'np.identity('+trans_ind_str+')'
          einsum_list.append(lhs_str+" "+sign+"="+cons+term.constants[0]+" * "+'np.identity('+trans_ind_str+')\n')
        else:
          lhs_ind_str = term.constants[0][0]
          print lhs_str+" "+sign+"="+cons+term.constants[0]
          einsum_list.append(lhs_str+" "+sign+"="+cons+term.constants[0]+"\n") 
     else:
        print lhs_str+" "+sign+"="+cons+" np.einsum('"+lhs_ind_str+rhs_ind_str+"', "+rhs_ten_str+","+IOpt+")"+Icomnd+icopy
        einsum_list.append(lhs_str+" "+sign+"="+cons+" np.einsum('"+lhs_ind_str+rhs_ind_str+"', "+rhs_ten_str+","+IOpt+")"+Icomnd+icopy+"\n")
          
     term1st += 1
# 
 return einsum_list             ## Changes made by Rajat

def einsum_help():
    print("""\n        HELP :: 
        -----------
        terms         : A list of terms
        lhs_str       : Left hand side string (e.g. string 'M' = einsum ..)
        ind_str       : Einsum right side indix string (e.g. -> string 'p')
        transRDM      : Transition RDM True of False
        trans_ind_str : Transition RDM string if True
        rhs_str       : Extra string for other kind of operation 
                        (e.g. transpose, copy, reshape .. etc)
        optimize      : By default optimization is true
        suffix        : Additional string atachement to the tensor name 
                        ( By default suffix = 'so' for spin orbitlas)
        rdm_str       : RDM string name
        tensor_rename : Rename tensor if required 
                        (e.g. rename tensor 'X' to 'TEMP': X = 'TEMP'. 
                        For multiple tensors: X = 'TEMP', h = 'Hamiltonian', ..)
--------------------------------------------------------------""")
    yes = {'Yes','yes','y', 'Y', ''}
    no = {'NO','no','n','N'}
    sys.stdout.write("Do you want to continue [y/n] : ")
    choice = raw_input().lower()
    if choice in no:
       exit()
    print("-------------------------------------------------------------- ")

def tensor_name_indices(tensr, suffix = ''):
# Return Tensor name and indices as string
 tensr_ind_str = ''
 tensr_indtype_str = ''
 tensor_name = str(tensr.name)
#
 for ind in range(len(tensr.indices)):
     tensr_ind_str += str(tensr.indices[ind].name)
#
     if not (tensr.indices[ind].indType[0][0] == 'virtual'):
        tensr_indtype_str += str(tensr.indices[ind].indType[0][0][0])
     else:
       tensr_indtype_str += 'e'
#
#     if (isinstance(tens, creOp) or isinstance(tens,desOp)):
#
     if (isinstance(tensr, kroneckerDelta)):
       tensor_name = 'np.identity'
       if (tensr.indices[0].indType[0][0] == 'core' and (tensr.indices[1].indType[0][0] == 'core')):
          dstr = 'ncore'
       elif (tensr.indices[0].indType[0][0] == 'active' and (tensr.indices[1].indType[0][0] == 'active')):
          dstr = 'ncas'
       else:
          dstr = 'nextern'
       tensor_name += str('(')+dstr+'_'+suffix+str(')')
#
     elif (((tensr.name == 'E') or (tensr.name == 'e')) and (len(tensr.indices) == 1)):
       tensr_indtype_str = tensr.indices[0].indType[0][0]
       if (tensr.indices[ind].indType[0][0] == 'virtual'):
          tensr_indtype_str = 'extern'     # Change name 'virtual' to 'extern'
     else:
       if (tensr.name == 'gamma'):
          tensr_indtype_str = 'ca'
#      
 return tensor_name, tensr_ind_str, tensr_indtype_str

def custom_tensor(tname, *tup):
# Switch symmetry either one of (bra / ket) or both
# Implemented only for T2 tensors
 ind = list(tup)
 if (len(ind) != 4):
    raise Exception("Implemented only for 4 indices, not for %s ." % (len(ind)))
 else:
    if (ind[0].indType == ind[1].indType) and (ind[2].indType == ind[3].indType):
       symm  =  [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
    elif (ind[0].indType == ind[1].indType) and (ind[2].indType != ind[3].indType):
       symm  =  [ symmetry((1,0,2,3),-1)]
    elif (ind[0].indType != ind[1].indType) and (ind[2].indType == ind[3].indType):
       symm  =  [ symmetry((0,1,3,2), -1)]
    else:
       symm  = []

 tname_tensor= tensor(tname, ind, symm)

 return tname_tensor


def sqalatex(terms, lhs = None, output = None, indbra = False, indket = None, print_default = True):

 if not output:
  # texfile = r'latex_output.tex'
   texfile = r'output_default'
 else:
   texfile = output

 print("""\n----------------------- SQA LATEX ----------------------------
    _____ ____    ___   __
   / ___// __ \  /   | / /____  _  __
   \__ \/ / / / / /| |/ __/ _ \| |/_/  Translate to Latex format and generate pdf
  ___/ / /_/ / / ___ / /_/  __/>  <    author:  Koushik Chatterjee
 /____/\___\_\/_/  |_\__/\___/_/|_|    date:  April 28, 2019
                                       VERSION : 1
 Copyright (C) 2018-2020  Koushik Chatterjee (koushikchatterjee7@gmail.com)
 
 Tex file : %s
 PDF file : %s
--------------------------------------------------------------""" % (texfile+r'.tex', texfile+r'.pdf'))

 modifier_tensor = {
     'bold': lambda s: r'\boldsymbol{'+s+r'}',
     'hat': lambda s: r'\hat{'+s+r'}',
     'bra': lambda s: r'\langle\Psi_{'+s+r'}\lvert',
     'ket': lambda s: r'\rvert\Psi_{'+s+r'}\rangle',
#     'gamma': lambda s: r'\Gamma',
     'kdelta': lambda s: r'\delta',
     'cre': lambda s: r'\hat{'+s+r'}^{\dagger}',
     'des': lambda s: r'\hat{'+s+r'}',
 }
# t_modifier = lambda s: r'\boldsymbol{'+s+r'}'
 t_modifier = lambda s: s

 if not lhs:
  lhs = 'M={}'
 else:
  lhs = t_modifier(lhs)+r'={}'

 tex = []

 for term in terms:

     constant = ''
     if (term.numConstant == 1.0):
        constant = " + "
     elif (term.numConstant == -1.0):
        constant = " - "
     else:
        constant = " %s " % str(term.numConstant)
        if (term.numConstant > 0):
            
#          constant += " %s " % str(Fraction(Decimal('term.numConstant')))
          constant = " +%s " % str(term.numConstant)
#
     cre_count = 0
     des_count = 0
     credes = ''
     name = ''
     gamma = ''
     for i in range(len(term.tensors)):

         tens = term.tensors[i]
         s = tens.name
     #    credes = None
     #    name = ''

         supers = ''
         subs   = ''
         index = len(tens.indices)
         if (index == 1):
            subs   = tens.indices[0].name
         elif(index == 2):
            supers = tens.indices[0].name
            subs   = tens.indices[1].name
         elif (index == 4):
            supers = tens.indices[0].name+tens.indices[1].name
            subs   = tens.indices[2].name+tens.indices[3].name

         else:
             raise Exception("Not implemented ...")

         if not (isinstance(tens, creOp) or isinstance(tens, desOp)):
            if (s == 'gamma'):
               bra = modifier_tensor['bra']('0')
               ket = modifier_tensor['ket']('0')
               ind1 = modifier_tensor['cre'](supers)
               ind2 = modifier_tensor['des'](subs)
               gamma += bra+ind1+ind2+ket+"\:"
            else:
               if s in modifier_tensor:
                  name += modifier_tensor[s](s)
               else:
                  name += t_modifier(s)

               name += "^{%s}" % " ".join(supers)
               name += "_{%s}" % " ".join(subs)
               name +="\:"
         else:
            if (isinstance(tens, creOp)):
               cre_count += 1
            if (isinstance(tens, desOp)):
               des_count += 1
            credes += modifier_tensor[s](subs)

     if(len(gamma) > 0):
        name += gamma
     if (len(credes) > 0):
         ind = r'0'
         if not indbra:
            indbra = ind
         if not indket:
            indket = ind
         bra = modifier_tensor['bra'](indbra)
         ket = modifier_tensor['ket'](indket)
         name += bra+credes+ket

     tex.append(constant+r'\:'+name)


 if print_default:
    print r'\documentclass{article}'
    print r'\usepackage{amsmath}'
    print r'\begin{document}'
    print ''
    print ''
#    print r"\begin{equation}"
    print r"\begin{align*}"
    print lhs
    for i in tex:
#      print " & "+i+' \\\\'
      print " & "+i+r'\\'
    print r"\end{align*}"
#    print r"\end{equation}"
    print ''
    print ''
    print r'\end{document}'


 ### write to a file ###
# if not output:
#  # texfile = r'latex_output.tex'
#   texfile = r'latex_output'
# else:
#   texfile = output
 output = open(texfile+r'.tex', "w")
 output.write(r'\documentclass{article}')
 output.write("\n")
 output.write(r'\usepackage{amsmath}')
 output.write("\n")
 output.write(r'\begin{document}')
 output.write("\n")
 output.write('')
 output.write("\n")
 output.write('')
 output.write("\n")
# output.write(r"\begin{equation*}")
 output.write(r"\begin{align*}")
 output.write("\n")
 output.write(lhs)
 for i in tex:
#   output.write(" & "+i+' \\\\')
   output.write(" & "+i+r'\\')
   output.write("\n")
 output.write(r"\end{align*}")
 output.write("\n")
# print r"\end{equation*}"
 output.write('')
 output.write("\n")
 output.write('')
 output.write("\n")
 output.write(r'\end{document}')

# os.system("pdflatex latex_output.tex")
 procs = []
 try:
     pread, pwrite = os.pipe()
     cmd = ['pdflatex', texfile+r'.tex']
#     proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
     proc = subprocess.Popen(cmd, stdout=pwrite, stderr=subprocess.STDOUT)
     procs.append(proc)
     os.close(pwrite)
     os.close(pread)

 except OSError as e:
   #  sys.exit()
     print 'Latex compilation error ...'

# pdf()
# proc_cleanup(procs)
 return


#
def getT(order = 1, cc = None, aa = None, vv = None):
 "Get T amplitudes (order)."
#
 tg_c = options.core_type
 tg_a = options.active_type
 tg_v = options.virtual_type
 tg_g = tg_c + tg_a + tg_v
 dummy = True
#
 if not cc:
   cc = [index('c%i' %p, [tg_c], dummy) for p in range(200,230)]
 if not aa:
   aa = [index('a%i' %p, [tg_a], dummy) for p in range(200,230)]
 if not vv:
   vv = [index('v%i' %p, [tg_v], dummy) for p in range(200,230)]
#
 T = Tamplitude(order, cc, aa, vv)
#
 return T


def getT_excitation(order = 1, cc = None, aa = None, vv = None):
 "Get T amplitudes (order)."
#
 tg_c = options.core_type
 tg_a = options.active_type
 tg_v = options.virtual_type
 tg_g = tg_c + tg_a + tg_v
 dummy = True
#
 if not cc:
   cc = [index('c%i' %p, [tg_c], dummy) for p in range(200,230)]
 if not aa:
   aa = [index('a%i' %p, [tg_a], dummy) for p in range(200,230)]
 if not vv:
   vv = [index('v%i' %p, [tg_v], dummy) for p in range(200,230)]
#
 T = Tamplitude_excitation(order, cc, aa, vv)
#
 return T


def getT_deexcitation(order = 1, cc = None, aa = None, vv = None):
 "Get T amplitudes (order)."
#
 tg_c = options.core_type
 tg_a = options.active_type
 tg_v = options.virtual_type
 tg_g = tg_c + tg_a + tg_v
 dummy = True
#
 if not cc:
   cc = [index('c%i' %p, [tg_c], dummy) for p in range(200,230)]
 if not aa:
   aa = [index('a%i' %p, [tg_a], dummy) for p in range(200,230)]
 if not vv:
   vv = [index('v%i' %p, [tg_v], dummy) for p in range(200,230)]
#
 T = Tamplitude_deexcitation(order, cc, aa, vv)
#
 return T


def getV(cc = None, aa = None, vv = None):
 "Get V pertubation terms."
#
 tg_c = options.core_type
 tg_a = options.active_type
 tg_v = options.virtual_type
 tg_g = tg_c + tg_a + tg_v
 dummy = True
#
 if not cc:
   cc = [index('c%i' %p, [tg_c], dummy) for p in range(200)]
 if not aa:
   aa = [index('a%i' %p, [tg_a], dummy) for p in range(200)]
 if not vv:
   vv = [index('v%i' %p, [tg_v], dummy) for p in range(200)]
#
 V = Vperturbation(cc, aa, vv)
#
 return V
#
