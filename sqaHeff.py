# file:  sqaHeff.py
# author:  Koushik Chatterjee
# date:  September 28, 2018
# summary:  
#           Heff : Construct effective Hamiltonian(L^N) of order 'N'.
#
# Copyright (C) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
#
#import sqa_extra.secondQuantizationAlgebra as sqa
import sys
import collections

from sqaIndex import index
from sqaCommutator import commutator
from sqaTerm import term, combineTerms
from sqaTensor import tensor, creOp, desOp, kroneckerDelta
from sqaOptions import options
from sqaSymmetry import symmetry

#####################################
#
def Heff(order):
# print_header()
# "Construct effective Hamiltonian(L)."
 print("""\n--------------------------------------------------------------
 Effective Hamiltonian: Construct effective Hamiltonian ...
 author:  Koushik Chatterjee
 date:  August 31, 2018

 Copyright (C) 2018  Koushik Chatterjee (koushikchatterjee7@gmail.com)
 --------------------------------------------------------------""")
#   order = 0 : L(0) = H(0)
#   order = 1 : L(1) = V + [H(0),T(1) - T'(1)]
#   order = 2 : L(2) = [H(0),(T(2) - T'(2))] + 1/2[(V + L(1)), (T(1) - T'(1))] 
#
# print "Construct effective Hamiltonian:=>"
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
 cc = [index('c%i' %p, [tg_c], dummy) for p in range(50)]
# Active dummy indices
 aa = [index('a%i' %p, [tg_a], dummy) for p in range(50)]
# Virtual dummy indices
 vv = [index('v%i' %p, [tg_v], dummy) for p in range(50)]
#
 effH = []
 Hamil = []
# E_fc :
 c = index('Const.', [], dummy)
 Efc = tensor('E_fc',[c], [])
 Hamil.append( term(1.0, ['E_fc'], []))
# Hamil.append( term(1.0, [], [Efc]))
#
 cor = cc.pop(0)
 vir = vv.pop(0) 
# core and vitual part : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
 e_core = tensor('e', [cor], [])
 e_virt = tensor('e', [vir], [])
# Hamil.append( term(-1.0, ['e_i'],[ desOp(cor), creOp(cor)]))
# Hamil.append( term(1.0, ['e_v'],[ creOp(vir), desOp(vir)])) 
 Hamil.append( term(-1.0, [],[e_core, desOp(cor), creOp(cor)]))
 Hamil.append( term(1.0, [],[e_virt, creOp(vir), desOp(vir)]))
#
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
 Hact.append( term(0.25, [], [v2,  creOp(act1), creOp(act2), desOp(act4), desOp(act3)]))
#
# for t in Hact:
#    print 'Hact=', t
#
 Hamil.extend(Hact)
# for t in Hamil:
#    print 'Hamiltonian(0)=', t
#
 if (order == 0):
    L0 = []
    L0.extend(Hamil)
    effH = L0
    return effH
#
############################
#
 elif (order >= 1):
# L(1) = V + [H(0),T(1) - T'(1)]
    L1 = []
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(30):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
#
    V = []
#
#    vtype = 'V[n=0]'
#    V = Vperturbation_type(V, cc1, aa1, vv1, vtype = 'V[n=0]')  # Example
    V = Vperturbation_type(cc1, aa1, vv1)
    L1.extend(V)
#
    cc1 = []
    aa1 = []
    vv1 = []
    for i in range(4):
        cc1.append(cc.pop(0))
        aa1.append(aa.pop(0))
        vv1.append(vv.pop(0))
#    ttype = 'full'
#    ttype = 'C-A'
    T1 = []
#    T1.extend(Tamplitude(T1, 1, cc1, aa1, vv1, ttype = 'full'))
    T1.extend(Tamplitude(T1, 1, cc1, aa1, vv1))
#
    com1 = commutator(Hamil, T1)
#
    L1.extend(com1)
#
    if (order == 1):
#   L(1) = V + [H(0),T(1) - T'(1)]
       effH = L1
       return effH
#
############################
#
    elif (order == 2):
#   L(2) = [H(0),T(2) - T'(2)]+ 1/2 [V + L(1),T(1) - T'(1)]
#
       L2 = []
#
       cc1 = []
       aa1 = []
       vv1 = []
       for i in range(4):
           cc1.append(cc.pop(0))
           aa1.append(aa.pop(0))
           vv1.append(vv.pop(0))
#
       T2 = []
#       T2.extend(Tamplitude(T2, 2, cc1, aa1, vv1, ttype = 'full'))
       T2.extend(Tamplitude(T2, 2, cc1, aa1, vv1))
#
       com2 = commutator(Hamil, T2)
       L2.extend(com2)
######################
# for checking purpose for T2
#       effH = L2
#       return effH
######################
#
# (V + L1)
       VL1 = []
       VL1.extend(V)
       VL1.extend(L1)    
#
       for t in VL1:
          t.scale(0.5)
#
       com3 = commutator(VL1, T1)
#       for term in com3:
#          term.numConstant = 0.5 * term.numConstant
       L2.extend(com3)
#
       effH = L2
#
       return effH
#
############################
    else:
       raise Exception('Unknown type of effective Hamiltonian of order = %s' % (order))
#
#####################################
#
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
    print "Perturbation(V) type = All types of V"
    print ""
    V.extend(Vperturbation(cc, aa, vv))
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
def Tamplitude(T, order, cc1, aa1, vv1, ttype = None):
# Cluster operator : T - T^dag, Where T = T1 + T2
# Single excitatio : T1
#
# Define t amplitude according to their order
 if (order == 1):
    tname = 't1'
 elif (order == 2):
    tname = 't2'
#
 t1_sym = [ symmetry((1,0),1)]
 t2_sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
#
 T = []
 #Tex = []
 #Tdex = []
 T_CE = []
 T_CA = []
 T_AE = []
 T_othr = []
 for typ in range(4):
     cc = list(cc1)
     aa = list(aa1)
     vv = list(vv1)
     if (typ == 0):                # Core-External
         ind1 = cc.pop(0)
         ind2 = vv.pop(0)
         ind3 = aa.pop(0)
         ind4 = aa.pop(0)
         t1_tens =  tensor(tname, [ind2,ind1],t1_sym)
         t2_tens =  tensor(tname, [ind2,ind3,ind1,ind4],t2_sym)
         T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
         T2_ex =  term(1.0, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])
         T1_dex =  term(-1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])
         T2_dex =  term(-1.0, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])
#
         T_CE.append(T1_ex)
         T_CE.append(T2_ex)
         T_CE.append(T1_dex)
         T_CE.append(T2_dex)
     elif (typ == 1):              # Core-Active
         ind1 = cc.pop(0)
         ind2 = aa.pop(0)
         ind3 = aa.pop(0)
         ind4 = aa.pop(0)
         t1_tens =  tensor(tname, [ind2,ind1],t1_sym)
         t2_tens =  tensor(tname, [ind2,ind3,ind1,ind4],t2_sym)
         T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
         T2_ex =  term(0.5, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])
         T1_dex =  term(-1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])
         T2_dex =  term(-0.5, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])
#
         T_CA.append(T1_ex)
         T_CA.append(T2_ex)
         T_CA.append(T1_dex)
         T_CA.append(T2_dex)

     elif (typ == 2):              # Active-External
         ind1 = aa.pop(0)
         ind2 = vv.pop(0)
         ind3 = aa.pop(0)
         ind4 = aa.pop(0)
         t1_tens =  tensor(tname, [ind2,ind1],t1_sym)
         t2_tens =  tensor(tname, [ind2,ind3,ind1,ind4],t2_sym)
         T1_ex =  term(1.0, [], [t1_tens,  creOp(ind2), desOp(ind1)])
         T2_ex =  term(0.5, [], [t2_tens,  creOp(ind2), creOp(ind3), desOp(ind4), desOp(ind1)])
         T1_dex =  term(-1.0, [], [t1_tens,  creOp(ind1), desOp(ind2)])
         T2_dex =  term(-0.5, [], [t2_tens,  creOp(ind1), creOp(ind4), desOp(ind3), desOp(ind2)])
#
         T_AE.append(T1_ex)
         T_AE.append(T2_ex)
         T_AE.append(T1_dex)
         T_AE.append(T2_dex)
     else:                           # All types
         ind1 = cc.pop(0)
         ind2 = cc.pop(0)
         ind3 = aa.pop(0)
         ind4 = aa.pop(0)
         ind5 = vv.pop(0)
         ind6 = vv.pop(0)
# For other type of T2 excitations and de-excitations
         t2_tens1 =  tensor(tname, [ind5,ind6,ind1,ind2],t2_sym)
         T2_ex = term(0.25, [], [t2_tens1,  creOp(ind5), creOp(ind6), desOp(ind2), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens2 = tensor(tname, [ind5,ind3,ind1,ind2],t2_sym)
         T2_ex = term(0.5, [], [t2_tens2,  creOp(ind5), creOp(ind3), desOp(ind2), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens3 = tensor(tname, [ind5,ind6,ind1,ind3],t2_sym)
         T2_ex = term(0.5, [], [t2_tens3,  creOp(ind5), creOp(ind6), desOp(ind3), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens4 = tensor(tname, [ind3,ind4,ind1,ind2],t2_sym)
         T2_ex = term(0.25, [], [t2_tens4,  creOp(ind3), creOp(ind4), desOp(ind2), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens5 = tensor(tname, [ind5,ind6,ind4,ind3],t2_sym)
         T2_ex = term(0.25, [], [t2_tens5,  creOp(ind5), creOp(ind6), desOp(ind3), desOp(ind4)])
         T_othr.append(T2_ex)
#
         T2_dex = term(-0.25, [], [t2_tens1,  creOp(ind1), creOp(ind2), desOp(ind6), desOp(ind5)])
         T_othr.append(T2_dex)
         T2_dex = term(-0.5, [], [t2_tens2,  creOp(ind1), creOp(ind2), desOp(ind3), desOp(ind5)])
         T_othr.append(T2_dex)
         T2_dex = term(-0.5, [], [t2_tens3,  creOp(ind1), creOp(ind3), desOp(ind6), desOp(ind5)])
         T_othr.append(T2_dex)
         T2_dex = term(-0.25, [], [t2_tens4,  creOp(ind1), creOp(ind2), desOp(ind4), desOp(ind3)])
         T_othr.append(T2_dex)
         T2_dex = term(-0.25, [], [t2_tens5,  creOp(ind4), creOp(ind3), desOp(ind6), desOp(ind5)])
         T_othr.append(T2_dex)
#
# T = T-T^dag
 if not (ttype):
# Include all type of excitations and de-excitations
     T.extend(T_CE)
     T.extend(T_CA)
     T.extend(T_AE)
     T.extend(T_othr)
 else:
     if (ttype == 'C-E'):
# Core-External
         T.extend(T_CE)
     elif (ttype == 'C-A'):
# Core-Active
         T.extend(T_CA)
     elif (ttype == 'A-E'):
# Active-External
         T.extend(T_AE)
#
     else:
         raise Exception('Unknown type of T operator...')
#
# for t in T:
#     print 't ampli=', t
#
 return T
#
#####################################
#
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
#
# Define operator types
 tg_c = options.core_type
 tg_a = options.active_type
 tg_v = options.virtual_type
 tg_g = tg_c + tg_a + tg_v
#
# Define indices
 dummy = True
#
 p = index('p', [], dummy)
 q = index('q', [], dummy)
 r = index('r', [], dummy)
 s = index('s', [], dummy)
#
############################
# list1 = ['c%i' %p for p in range(30)]
# list2 = ['a%i' %p for p in range(30)]
# list3 = ['v%i' %p for p in range(30)]
 list1 = []
 list2 = []
 list3 = []
 for i in range(len(cc)):
    list1.append(cc[i].name)
 for i in range(len(cc)):
    list2.append(aa[i].name)
 for i in range(len(cc)):
    list3.append(vv[i].name)
############################
#
 indc1 = list1.pop(0)
 indc2 = list1.pop(0)
 inda1 = list2.pop(0)
 inda2 = list2.pop(0)
 indv1 = list3.pop(0)
 indv2 = list3.pop(0)
#
 cor1 = index(indc1, [tg_c], dummy)
 cor2 = index(indc2, [tg_c], dummy)
 act1 = index(inda1, [tg_a], dummy)
 act2 = index(inda2, [tg_a], dummy)
 vir1 = index(indv1, [tg_v], dummy)
 vir2 = index(indv2, [tg_v], dummy)
#
 ten1 =  tensor('h', [cor1, act1], h1sym)
 V.append( term(1.0, [], [ten1, creOp(act1), desOp(cor1)]))
#
 ten2 =  tensor('h', [act1, cor1], h1sym)
 V.append( term(1.0, [], [ten2, creOp(cor1), desOp(act1)]))
#
#
 ten3 =  tensor('h', [act1, vir1], h1sym)
 V.append( term(1.0, [], [ten3, creOp(vir1), desOp(act1)]))
#
#
 ten4 =  tensor('h', [vir1, act1], h1sym)
 V.append( term(1.0, [], [ten4, creOp(act1), desOp(vir1)]))
#
#
 ten5 =  tensor('h', [cor1, vir1], h1sym)
 V.append( term(1.0, [], [ten5, creOp(vir1), desOp(cor1)]))
#
#
 ten6 =  tensor('h', [vir1, cor1], h1sym)
 V.append( term(1.0, [], [ten6, creOp(cor1), desOp(vir1)]))
#
#
#
 ten7 =  tensor('v', [cor2, act2, cor1, act1], v2sym)
# V.append( term(1.0, [], [ten7, desOp(cor1), creOp(cor2), creOp(act2), desOp(act1)]))
#
 ten8 =  tensor('gamma', [act2, act1], d1sym)
 V.append( term(1.0, [], [ten7, ten8,  desOp(cor2), creOp(cor1)]))
#
 ten9 =  tensor('v', [vir2, act2, vir1, act1], v2sym)
# V.append( term(-1.0, [], [ten9, creOp(vir2), desOp(vir1), creOp(act2), desOp(act1)]))
 V.append( term(-1.0, [], [ten9, ten8,  creOp(vir1), desOp(vir2)]))
#
 for ityp1 in range(3):
#        list1 = list(coreInd)
#        list2 = list(actvInd)
#        list3 = list(virtInd)
        if (ityp1 == 0):                              # ityp1 = 0 => Core
               ind = list1.pop(0)
               p = index(ind, [tg_c], dummy)
               list1.append(ind)
        elif (ityp1 == 1):                            #         1 => Active
               ind = list2.pop(0)
               p = index(ind, [tg_a], dummy)
               list2.append(ind)
        else:                                         #         2 => Virtual
               ind = list3.pop(0)
               p = index(ind, [tg_v], dummy)
               list3.append(ind)
#
        for ityp2 in range(3):
               if (ityp2 == 0):
                      ind = list1.pop(0)
                      q = index(ind, [tg_c], dummy)
                      list1.append(ind)
               elif (ityp2 == 1):
                      ind = list2.pop(0)
                      q = index(ind, [tg_a], dummy)
                      list2.append(ind)
               else:
                      ind = list3.pop(0)
                      q = index(ind, [tg_v], dummy)
                      list3.append(ind)
#
               for ityp3 in range(3):
                      if (ityp3 == 0):
                             ind = list1.pop(0)
                             s = index(ind, [tg_c], dummy)
                             list1.append(ind)
                      elif (ityp3 == 1):
                             ind = list2.pop(0)
                             s = index(ind, [tg_a], dummy)
                             list2.append(ind)
                      else:
                             ind = list3.pop(0)
                             s = index(ind, [tg_v], dummy)
                             list3.append(ind)
#
                      for ityp4 in range(3):
                             if (ityp4 == 0):
                                  ind = list1.pop(0)
                                  r = index(ind, [tg_c], dummy)
                                  list1.append(ind)
                             elif (ityp4 == 1):
                                  ind = list2.pop(0)
                                  r = index(ind, [tg_a], dummy)
                                  list2.append(ind)
                             else:
                                  ind = list3.pop(0)
                                  r = index(ind, [tg_v], dummy)
                                  list3.append(ind)
#
                             if not (p.indType[0][0]=='active' and q.indType[0][0]=='active' and r.indType[0][0]=='active' and s.indType[0][0]=='active'):
#
                                   if (p.indType[0][0]=='core' and q.indType[0][0]=='core' and r.indType[0][0]=='core' and s.indType[0][0]=='core'):
#                                       vTen = tensor('v', [r,s,p,q], v2sym)
                                       vTen = tensor('v', [r,s,p,q], v2sym)
                                       V81.append(term(-0.25, [], [vTen,desOp(r), desOp(s), creOp(p), creOp(q)]))
#
                                   else:
                                       vTen = tensor('v', [r,s,p,q], v2sym)
#                                       vTen = tensor('v', [p,q,r,s], v2sym)
                                       V81.append(term(0.25, [], [vTen, creOp(p), creOp(q),desOp(s), desOp(r)]))
#
 
 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list3.pop(0), [tg_v], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list3.pop(0), [tg_v], dummy)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(-1.0, [], [vTen, creOp(p), creOp(q), desOp(s), desOp(r)]))

 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list2.pop(0), [tg_a], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list2.pop(0), [tg_a], dummy)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(-1.0, [], [vTen, creOp(p), creOp(q), desOp(s), desOp(r)]))

 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list3.pop(0), [tg_v], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list3.pop(0), [tg_v], dummy)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(1.0, [], [vTen, creOp(q), desOp(r), desOp(s), creOp(p)]))

 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list2.pop(0), [tg_a], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list2.pop(0), [tg_a], dummy)
 vTen = tensor('v', [r,s,p,q], v2sym)
 V81.append(term(1.0, [], [vTen, desOp(s), creOp(p), creOp(q), desOp(r)]))

 V.extend(V81)
# Dummy indices label upate
# dummyLabel(V)
# print len(V81)
 return V
#
#####################################
def print_header():

    print("""\n--------------------------------------------------------------
    SQA_extra: Construct effective Hamiltonian(L)
    author:  Koushik Chatterjee
    date:  August 31, 2018

    Copyright (C) 2018  Koushik Chatterjee (koushikchatterjee7@gmail.com)

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
# Copyright (C) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
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

 Copyright (C) 2018  Koushik Chatterjee (koushikchatterjee7@gmail.com)
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
# Copyright (C) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
# This program is distributed in the hope that it will
# be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License
# for more detai
#
# print "################ Construct Einsum ################"

 print("""\n--------------------------------------------------------------
 VERSION : 1
 Einsum generator: Transform into einsum ...
 author:  Koushik Chatterjee
 date:  August 31, 2018

 Copyright (C) 2018  Koushik Chatterjee (koushikchatterjee7@gmail.com)


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
         tensor_name, tensr_ind_str, tensr_indtype_str = tensor_name_indices(tens)
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
        name_default.update({key_tensr : value_tensr})
#        name_default.update({tensorlist[i][0] : tensorlist[i]})
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
               list_tuple = list(value1)
               list_tuple[0] = value0
               external_tensor.append(value0)
               value1 = tuple(list_tuple)
               name_default[key1] = value1
              # del name_default[key1]
              # name_default.update({value0 : value1})
#
     LHS_ind = ''
     RHS_tensr = ''
     lhs_ein_ind = []
     rhs_ein_ten = []

     for key, value in name_default.items():
         lhs_ein_ind.append(value[1])
#
         if (value[0] in external_tensor):
            ein_tens_name = value[0]
         else:
            if ((value[0] == 't1') or (value[0] == 't2')):
               ein_tens_name = value[0]+'_'+value[2]
            else:
               ein_tens_name = value[0]+'_'+value[2]+'_'+suffix
#
#         rhs_ein_ten.append(value[0]+'_'+value[2])
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

     if (len(term.tensors)==0):
        if (transRDM):
          print lhs_str+" "+sign+"="+cons+term.constants[0]+" * "+'np.identity('+trans_ind_str+')'
        else:
          lhs_ind_str = term.constants[0][0]
          print lhs_str+" "+sign+"="+cons+term.constants[0]
     else:
        print lhs_str+" "+sign+"="+cons+" np.einsum('"+lhs_ind_str+rhs_ind_str+"', "+rhs_ten_str+","+IOpt+")"+Icomnd+icopy
#
     term1st += 1
# 
 return

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

def tensor_name_indices(tensr):
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
       tensor_name += str('(')+dstr+str(')')
#
     elif (((tensr.name == 'E') or (tensr.name == 'e')) and (len(tensr.indices) == 1)):
       tensr_indtype_str = tensr.indices[0].indType[0][0]
     else:
       if (tensr.name == 'gamma'):
          tensr_indtype_str = 'ca'
#      
 return tensor_name, tensr_ind_str, tensr_indtype_str
