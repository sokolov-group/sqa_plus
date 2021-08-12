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
from sqaTensor2 import creDesTensor
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
 c = index('Const.', [], True)
 Efc = tensor('E_fc',[c], [])
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
# v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1),  symmetry((2,3,0,1), 1)]
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
       ten6 = creDesTensor([creOp(act2), desOp(act1)])
#       ten6 =  tensor('gamma', [act2, act1], d1sym)
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


def Vperturbation(cc, aa, vv):
 from sqaAddon import matrixBlock, dummyLabel
#
 "Construct general perturbation operator V full (default) include all types of V."
#
 V = []
 V81 = []
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
 v2sym_braket = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1),  symmetry((2,3,0,1), 1)]
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
 ten8 =  creDesTensor([creOp(act2), desOp(act1)])
 #ten8 =  tensor('gamma', [act2, act1], d1sym)
 V.append( term(1.0, [], [ten7, ten8,  desOp(cor2), creOp(cor1)]))

 cor1 = cc.pop(0)
 vir1 = vv.pop(0)
 act1 = aa.pop(0)
 cor2 = cc.pop(0)
 vir2 = vv.pop(0)
 act2 = aa.pop(0)

 #ten8 =  tensor('gamma', [act2, act1], d1sym)
 ten8 =  creDesTensor([creOp(act2), desOp(act1)])
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
                                       vTen = tensor('v', [r,s,p,q], v2sym_braket)
                                       V81.append(term(-0.25, [], [vTen, desOp(r), desOp(s), creOp(p), creOp(q)]))
#
                                   elif (p.indType[0][0]=='virtual' and q.indType[0][0]=='virtual' and r.indType[0][0]=='virtual' and s.indType[0][0]=='virtual'):
                                       vTen = tensor('v', [r,s,p,q], v2sym_braket)
                                       V81.append(term(0.25, [], [vTen, creOp(p), creOp(q),desOp(s), desOp(r)]))
#
                                   else:
                                       vTen = tensor('v', [r,s,p,q], v2sym)
                                       V81.append(term(0.25, [], [vTen, creOp(p), creOp(q),desOp(s), desOp(r)]))
 
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


#####################################

def genEinsum(terms, lhs_str = None, ind_str = None, trans_rdm = False, trans_ind_str = None, suffix = None, rm_trans_rdm_const = False, rm_core_int = False, intermediate_list = None, opt_einsum_terms = True, optimize = True, help = False, **tensor_rename):

    # Store custom names if provided by user
    custom_names = []
    if tensor_rename:
        for old_name, new_name in tensor_rename.items():
            custom_names.append((old_name, new_name))

    # Default to 'temp' as name of matrix being created w/ einsum function
    if not lhs_str:
        lhs_str = 'temp'

    # Default to spin-orbital suffix if not defined by user
    if not suffix:
        suffix = 'so'

    ################################################    
    # IF PROVIDED, PRINT EINSUMS FOR INT TERMS
    ################################################    
    if intermediate_list:

        # Create empty to list to store einsum expressions for provided intermediate terms
        int_einsum_list = []

        # Determined which intermediates are defined w/ trans_rdm contraction
        trans_int = []
        if trans_rdm:
            trans_int = get_trans_intermediates(intermediate_list)

        # If using effective Hamiltonian, remove double-counted contributions to core terms
        if rm_core_int:
            intermediate_list, removed_int = remove_core_int(intermediate_list, int_terms = True)

        # Iterate through the list of provided intermediates
        for int_ind, (int_term, int_tensor) in enumerate(intermediate_list):
 
            # Pass tensors of term to function to create string representation of contraction indices and tensor names
            int_tensor_inds, int_tensor_names = get_tensor_info(intermediate_list[int_ind][0].tensors, trans_rdm, trans_ind_str, ''.join([i.name for i in intermediate_list[int_ind][1].indices]), suffix, trans_int, custom_names)

            # Rename intermediates in tensor definitions
            if custom_names:
                if 'INT' in [x for x,y in custom_names]:
                    new_name = make_custom_name(int_tensor, custom_names)
                    int_einsum = new_name + ' = '

            else:
                int_einsum = int_tensor.name + ' = '
 
            # Define term for either optEinsum or built-in Numpy 'einsum' function
            if opt_einsum_terms:
                int_einsum += 'einsum('
            else:
                int_einsum += 'np.einsum('
 
            # Add contraction and tensor info
            int_tensor_info = (', '.join([str("'") + int_tensor_inds + str("'")] + int_tensor_names))
            int_einsum += int_tensor_info
 
            # Add optimize flag to einsum if enabled
            if optimize and opt_einsum_terms:
                int_einsum += ', optimize = einsum_type)'
            elif optimize and not opt_einsum_terms:
                int_einsum += ', optimize = True)'
            else:
                int_einsum += ')'
 
            # Append einsum definition to list, returned at the end of function
            int_einsum_list.append(int_einsum)

    ################################################    
    # GENERATE EINSUM EXPRESSIONS FOR PROVIDED TERMS
    ################################################    

    # Constants terms in CAS blocks are removed by default, print warning
    if trans_rdm and rm_trans_rdm_const and intermediate_list and trans_int:
        terms, const_terms = remove_trans_rdm_const(terms, trans_int)

    elif trans_rdm and rm_trans_rdm_const:
        terms, const_terms = remove_trans_rdm_const(terms)

    # If using effective Hamiltonian, remove double-counted contributions to core terms
    if intermediate_list and rm_core_int and removed_int:
        terms, core_terms = remove_core_int(terms, removed_int)

    elif rm_core_int:
        terms, core_terms = remove_core_int(terms)
    
    # Create empty list for storing einsums
    einsum_list = []

    # Iterate through terms and create einsum expressions
    for term_ind, term in enumerate(terms):

        ## CONVERT ALL CRE/DES OBJECTS to CREDESTENSOR OBJECT
        # List for storing cre/des operators
        credes_ops = []

        # Append all cre/des operators to list
        for tens in term.tensors:
            if isinstance(tens, creOp) or isinstance(tens, desOp):
                credes_ops.append(tens)

        # Modify term in list to use creDesTensor object instead of cre/des objects
        if credes_ops:
            terms[term_ind].tensors = [ten for ten in terms[term_ind].tensors if ten not in credes_ops]
            terms[term_ind].tensors.append(creDesTensor(credes_ops, trans_rdm))

        ## PROCEED WITH GENERATING EINSUMS
        # Start to define einsum string
        einsum = lhs_str + ' '

        # Determine sign of term
        pos = True
        if term.numConstant < 0:
            pos = False

        # Set up equals sign for first term and rest of terms
        if term_ind == 0 and pos:
            einsum += ' = '
        elif term_ind == 0 and not pos:
            einsum += '=- '
        elif term_ind != 0 and pos:
            einsum += '+= '
        elif term_ind != 0 and not pos:
            einsum += '-= '

        # Add appropriate scaling factor
        if abs(term.numConstant) != 1:
            einsum = einsum + str(abs(term.numConstant)) + ' * '

        # Define term for either optEinsum or built-in Numpy 'einsum' function
        if opt_einsum_terms:
            einsum += 'einsum('
        else:
            einsum += 'np.einsum('

        # Pass tensors of term to function to create string representation of contraction indices and tensor names
        if trans_rdm and intermediate_list:
            tensor_inds, tensor_names = get_tensor_info(term.tensors, trans_rdm, trans_ind_str, ind_str, suffix, trans_int, custom_names)
        else:
            tensor_inds, tensor_names = get_tensor_info(term.tensors, trans_rdm, trans_ind_str, ind_str, suffix, custom_names)

        # Add contraction and tensor info
        tensor_info = (', '.join([str("'") + tensor_inds + str("'")] + tensor_names))
        einsum += tensor_info

        # Add optimize flag to einsum if enabled
        if optimize and opt_einsum_terms:
            einsum += ', optimize = einsum_type)'
        elif optimize and not opt_einsum_terms:
            einsum += ', optimize = True)'
        else:
            einsum += ')'

        # Append a '.copy()' function call if the term is made up of only one tensor
        if len(term.tensors) == 1:
            einsum += '.copy()'

        # Append completed einsum to list
        einsum_list.append(einsum)

    # Modify return for intermediate term definition
    if intermediate_list:
        return int_einsum_list, einsum_list
    else:
        return einsum_list


def get_tensor_info(sqa_tensors, trans_rdm, trans_ind_str, ind_str, suffix, trans_int = None, custom_names = None):

    # Pre-define list of names of tensors used in SQA and make list to store any new tensor 'types'
    tensor_names = []
    tensor_inds  = []

    # Iterate through all the provided tensors
    for tens in sqa_tensors:

        # Handle special case of kroneckerDelta (kdelta) object
        if isinstance(tens, kroneckerDelta):

            tensor_name = 'np.identity('

            # Determine orbital space of kdelta
            if (tens.indices[0].indType[0][0][0] and tens.indices[1].indType[0][0][0]) == 'c':
                orb_space = 'ncore'

            elif (tens.indices[0].indType[0][0][0] and tens.indices[1].indType[0][0][0]) == 'a':
                orb_space = 'ncas'

            elif (tens.indices[0].indType[0][0][0] and tens.indices[1].indType[0][0][0]) == 'v':
                orb_space = 'nextern'

            else:
                raise TypeError('WARNING: The indices of the kronecker delta term do not belong to the same orbital sub-space')

            if suffix:
                orb_space += '_' + suffix
            tensor_name += orb_space + ')'

            # Rename if custom name is provided
            if custom_names:
                if ('kdelta') in [x for x,y in custom_names]:
                    new_name = make_custom_name(tens, custom_names)
                    tensor_name = new_name + '_'

        # Handle special case of orbital energy vector
        elif len(tens.indices) == 1 and (tens.name == 'e' or tens.name == 'E'):

            tensor_name = str(tens.name).lower() + '_'

            # Determine orbital space of energies
            if (tens.indices[0].indType[0][0][0]) == 'c':
                orb_space = 'core'

            elif (tens.indices[0].indType[0][0][0]) == 'v':
                orb_space = 'extern'

            if suffix:
                orb_space += '_' + suffix
            tensor_name += orb_space

            # Rename if custom name is provided
            if custom_names:
                if ('e' or 'E') in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Handle special case of RDM tensor
        elif isinstance(tens, creDesTensor):

            tensor_name = tens.name + '_'
            
            # Modify name of RDM to reflect particle number
            for op in tens.ops:
                if isinstance(op, creOp):
                    tensor_name += 'c'
                elif isinstance(op, desOp):
                    tensor_name += 'a'

            # Append suffix
            if suffix:
                tensor_name += '_' + suffix

            # Rename if custom name is provided
            if custom_names:
                if ('rdm') in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Name remaining tensors w/ same convention of orbital space and suffix
        elif tens.name == 'h' or tens.name == 'v' or tens.name == 't1' or tens.name == 't2':
            
            tensor_name = tens.name + '_'

            # Append letter representing orbital subspace of indices
            for i in range(len(tens.indices)):
                if tens.indices[i].indType[0][0][0] != 'v':
                    tensor_name += tens.indices[i].indType[0][0][0]
                else:
                    tensor_name += 'e'

            # Append suffix
            if not (tens.name == 't1' or tens.name == 't2') and suffix:
                tensor_name += '_' + suffix

            # Rename if custom name is provided
            if custom_names:
                if tens.name in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Account for intermediate tensors and any custom tensors
        else:
            # Make copy of tensor name
            tensor_name = '%s' % tens.name

            # Allow to rename intermediate tensors in term definitions             
            if custom_names:
                if tens.name[:3] == 'INT' and 'INT' in [x for x,y in custom_names]:
                    tensor_name = make_custom_name(tens, custom_names)

        # Append name of tensor (after and modifications due to special cases)
        tensor_names.append(tensor_name)

        # Create indices of tensor as string
        indices = ''.join([i.name for i in tens.indices])

        # Append transition state index to appropriate set of indices
        if isinstance(tens, creDesTensor) and trans_rdm:
            indices = trans_ind_str + indices
            ind_str = trans_ind_str + ind_str            

        elif trans_int and (tens.name in trans_int):
            indices = trans_ind_str + indices
            ind_str = trans_ind_str + ind_str            

        # Append completed index string to list
        tensor_inds.append(indices)

    # Convert list of indices into one comma-separated string and prepare to append external index string
    tensor_inds = ','.join(tensor_inds)
    tensor_inds += '->'
    tensor_inds += ind_str

    return tensor_inds, tensor_names


def remove_core_int(terms, removed_int = None, int_terms = False):

    # Remove terms from standard term list
    if not int_terms:
        print ('--------------------------------- WARNING ---------------------------------')
        print ('Terms with a contraction over repeating dummy core indices of 2e- integrals')
        print ('will be removed. Set "rm_core_int" flag to FALSE to preserve terms')
    
        # Create lists to split up SQA terms
        kept_terms = []
        core_terms = []
    
        # Separate out the terms that have redundant 2e- integral contractions over core space
        for term_ind, term in enumerate(terms):
            coreTerm = False
            for tens_ind, tens in enumerate(term.tensors):
                if tens.name == 'v':
                    if (
                        (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[2].name)
                        or (terms[term_ind].tensors[tens_ind].indices[0].name) == (terms[term_ind].tensors[tens_ind].indices[3].name)
                        or (terms[term_ind].tensors[tens_ind].indices[1].name) == (terms[term_ind].tensors[tens_ind].indices[2].name)
                        or (terms[term_ind].tensors[tens_ind].indices[1].name) == (terms[term_ind].tensors[tens_ind].indices[3].name)
                    ):
                        coreTerm = True
                        break

                elif removed_int and (tens.name in removed_int):
                    coreTerm = True
                    break
           
            # Append to either list based on coreTerm flag
            if not coreTerm:
                kept_terms.append(terms[term_ind])
    
            else:
                core_terms.append(terms[term_ind])
    
        print ('')
        print (str(len(core_terms)) + ' terms removed:')
        for term in core_terms:
            print (term)
        
        print ('---------------------------------------------------------------------------')
        print ('Remaining terms: ' + str(len(kept_terms)))
        print ('')
    
        return kept_terms, core_terms

    # Filter through intermediate definitions
    else:
        print ('--------------------------------- WARNING ---------------------------------')
        print ('Intermediate tensors defined w/ contractions over repeating dummy core indices of')
        print ('2e- integrals will be removed. Set "rm_core_int" flag to FALSE to preserve definitions')

        # Track which tensor definitions are removed and kept
        removed_int   = []
        removed_terms = []

        # Determine which intermediate definitions to remove
        for int_ind, (int_term, int_tensor) in enumerate(terms):
            for tens in int_term.tensors:

                # If intermediate is defined w/ 2e- integral
                if tens.name == 'v':
                    if (
                        (tens.indices[0].name) == (tens.indices[2].name)
                        or (tens.indices[0].name) == (tens.indices[3].name)
                        or (tens.indices[1].name) == (tens.indices[2].name)
                        or (tens.indices[1].name) == (tens.indices[3].name)
                    ):
                        removed_int.append(int_tensor.name)
                        removed_terms.append(terms[int_ind])
                        break

                # If intermediate is defined in terms of one of the intermediates to be removed
                elif tens.name in removed_int:
                    removed_int.append(int_tensor.name)
                    removed_terms.append(terms[int_ind])
                    break

        # If some intermediate definitions were removed
        if removed_int:
            print ('')
            print (str(len(removed_int)) + ' definitions removed:')
            for tens, term in zip(removed_int, removed_terms):
                print (tens + ": " + str(term[0])) 
        
        print ('---------------------------------------------------------------------------')
        print ('')

        # Returned shortened intermediate list
        terms = [t for t in terms if t not in removed_terms]

        return terms, removed_int


def remove_trans_rdm_const(terms, trans_int_list = None):

    print ('--------------------------------- WARNING ---------------------------------')
    print ('Terms w/o transRDM tensor in the expression will be removed. Set "rm_trans_rdm_const"')
    print ('flag to FALSE to preserve terms')

    # Create lists to split up SQA terms
    const_terms     = []
    trans_rdm_terms = []

    # Remove terms without tRDM tensors in r.h.s.
    for term_ind, term in enumerate(terms):
        creDes = False

        for tensor in term.tensors:
            if isinstance(tensor, creOp) or isinstance(tensor, desOp) or isinstance(tensor, creDesTensor):
                creDes = True
                break       

            elif trans_int_list:
                tens_list = [tns.name for tns in term.tensors]
                for trans_int in trans_int_list:
                    if trans_int in tens_list:
                        creDes = True
                        break
 
        # Append to either list based on creDes flag
        if not creDes:
            const_terms.append(terms[term_ind])

        else:
            trans_rdm_terms.append(terms[term_ind])

    print ('')
    print (str(len(const_terms)) + ' terms removed:')
    for term in const_terms:
        print (term)
    
    print ('---------------------------------------------------------------------------')
    print ('Remaining terms: ' + str(len(trans_rdm_terms)))
    print ('')

    return trans_rdm_terms, const_terms


def get_trans_intermediates(intermediate_list):

    # Store which intermediates are contracted over transition index
    trans_int_list = []

    # Iterate through list of intermediates
    for int_term, int_tensor in intermediate_list:

        # Make list of tensors that define intermediates
        ten_list = [t.name for t in int_term.tensors]

        # Check if one of the tensors is a tRDM
        if 'trdm' in ten_list:
             trans_int_list.append(int_tensor.name)

        # If an intermediate is defined in terms of another intermediate, make sure that intermediate
        # isn't defined w/ a tRDM
        elif trans_int_list:
             for trans_int in trans_int_list:
                 if trans_int in ten_list:
                     trans_int_list.append(int_tensor.name)

    return trans_int_list


def make_custom_name(sqa_tensor, rename_tuple):

    old_name = [old for old, new in rename_tuple]

    if sqa_tensor.name[:3] == 'INT':
        rename_index = old_name.index('INT')
        new_name = rename_tuple[rename_index][1] + sqa_tensor.name[3:]

    else:
        rename_index = old_name.index(sqa_tensor.name)
        new_name = rename_tuple[rename_index][1]

    return new_name


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
   cc = [index('c%i' %p, [tg_c], dummy) for p in range(30)]
 if not aa:
   aa = [index('a%i' %p, [tg_a], dummy) for p in range(30)]
 if not vv:
   vv = [index('v%i' %p, [tg_v], dummy) for p in range(30)]
#
 T = Tamplitude(order, cc, aa, vv)
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


#############################
def analyzeTerm(input_terms, max_act_ind):

    # Initialize counters
    max_dim = 0
    max_dim_term = 0
    removed_terms = []

    for t_num, term in enumerate(input_terms, 1):
        print ('##### Term #{} #####'.format(t_num))
        print (term)
        print ('num of tensors: ', len(term.tensors))

        core_cre = 0
        core_des = 0
        act_cre = 0
        act_des = 0
        ext_cre = 0
        ext_des = 0


        # Iterate through every tensor in term's contraction
        for tensor in term.tensors:

            # Track largest contraction
            if max_dim < len(term.tensors):
                max_dim = len(term.tensors)
                max_dim_term = int(t_num)

            # Counting particle and hole operators in subspaces
            if isinstance(tensor, creOp):

                index = tensor.indices[0]  

                if index.indType[0][0][0] == 'c':
                    core_cre += 1

                if index.indType[0][0][0] == 'a':
                    act_cre += 1

                if index.indType[0][0][0] == 'v':
                    ext_cre += 1

            elif isinstance(tensor, desOp):
                
                index = tensor.indices[0]  
                
                if index.indType[0][0][0] == 'c':
                    core_des += 1

                if index.indType[0][0][0] == 'a':
                    act_des += 1

                if index.indType[0][0][0] == 'v':
                    ext_des += 1

            else:

                 print (tensor)

        print ('core_cre: ', core_cre)
        print ('core_des: ', core_des)
        print ('act_cre: ' , act_cre)
        print ('act_des: ' , act_des)
        print ('ext_cre: ' , ext_cre)
        print ('ext_des: ' , ext_des)

        # Filter out terms with more than request active indices
        if (act_cre + act_des) > max_act_ind:

            print ('####################################') 
	    print ('Term exceeds the maximum amount of active indices requested. Removing...') 
            print ('####################################') 

            removed_terms.append(input_terms[t_num])
            input_terms.pop(t_num)

        print ('\n')

    print ('Term #{} has largest tensor contraction'.format(max_dim_term))
    print (input_terms[max_dim_term - 1])
    print (str(max_dim) + ' tensors in contraction')

    print (str(len(removed_terms)) + ' terms removed')


    return input_terms, removed_terms



