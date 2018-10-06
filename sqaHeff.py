#    file:  sqaHeff.py
#  author:  Koushik Chatterjee
#    date:  September 28, 2018
# summary:  
#           Heff : Construct effective Hamiltonian(L^N) of order 'N'.
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
from sqaTerm import term, combineTerms
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
 e_core = tensor('E', [cor], [])
 e_virt = tensor('E', [vir], [])
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
    effH.extend(Hamil)
#
 elif (order == 1):
# L(1) = V + [H(0),T(1) - T'(1)]
    V = []
    T = []
#
#    vtype = 'V[n=0]'
    effH.extend(Vperturbation_type(V, cc, aa, vv, vtype = 'full'))
#
#    ttype = 'full'
#    ttype = 'C-A'
    T.extend(Tamplitude(T, cc, aa, vv, ttype = 'full'))
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
def Vperturbation_type(V, cc, aa, vv, vtype):
#
 "Construct perturbation operator V according to excitation rank."
#
 V = []
 v2sym = [ symmetry((1,0,2,3),-1),  symmetry((0,1,3,2), -1)]
 d1sym = [ symmetry((1,0),1)]
#
 if (vtype == 'full') or (vtype == 'Full') or (vtype == ''):
#   Default V includes all type of perturbation rank.
    print "Perturbation(V) type = All types of V"
    print ""
    V.extend(Vperturbation(V))
#
 elif (vtype == 'V[n=0]'):
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
    ten1 =  tensor('V', [cor3, cor4, cor1, cor2], v2sym)
#    V.append( term(0.25, [], [ten1, creOp(cor1), creOp(cor2), desOp(cor4), desOp(cor3)]))
    V.append( term(0.25, [], [ten1, desOp(cor4), desOp(cor3),creOp(cor1), creOp(cor2)]))
#
    ten2 =  tensor('V', [vir3, vir4, vir1, vir2], v2sym)
    V.append( term(0.25, [], [ten2, creOp(vir1), creOp(vir2), desOp(vir4), desOp(vir3)]))
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
    ten6 =  tensor('gamma', [act2, act1], d1sym)
#    V.append( term(-1.0, [], [ten5, ten6,  creOp(cor1), desOp(cor2)]))
    V.append( term(1.0, [], [ten5, ten6,  desOp(cor2), creOp(cor1)]))
#
    ten7 =  tensor('V', [vir2, act2, vir1, act1], v2sym)
    V.append( term(1.0, [], [ten7, creOp(vir1), desOp(vir2), creOp(act1), desOp(act2)]))
    ten8 = ten7
    ten9 = ten6
    V.append( term(-1.0, [], [ten8, ten9,  creOp(vir1), desOp(vir2)]))
 else:
    raise Exception('Unknown type of V..')
#
# for t in V:
#    print 'perterbative =', t
 return V
#####################################
#
def Tamplitude(T, cc1, aa1, vv1, ttype):
# Cluster operator : T - T^dag, Where T = T1 + T2
# Single excitatio : T1
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
         t1_tens =  tensor('t', [ind2,ind1],t1_sym)
         t2_tens =  tensor('t', [ind2,ind3,ind1,ind4],t2_sym)
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
         t1_tens =  tensor('t', [ind2,ind1],t1_sym)
         t2_tens =  tensor('t', [ind2,ind3,ind1,ind4],t2_sym)
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
         t1_tens =  tensor('t', [ind2,ind1],t1_sym)
         t2_tens =  tensor('t', [ind2,ind3,ind1,ind4],t2_sym)
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
         t2_tens1 =  tensor('t', [ind5,ind6,ind1,ind2],t2_sym)
         T2_ex = term(0.25, [], [t2_tens1,  creOp(ind5), creOp(ind6), desOp(ind2), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens2 = tensor('t', [ind5,ind3,ind1,ind2],t2_sym)
         T2_ex = term(0.5, [], [t2_tens2,  creOp(ind5), creOp(ind3), desOp(ind2), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens3 = tensor('t', [ind5,ind6,ind1,ind3],t2_sym)
         T2_ex = term(0.5, [], [t2_tens3,  creOp(ind5), creOp(ind6), desOp(ind3), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens4 = tensor('t', [ind3,ind4,ind1,ind2],t2_sym)
         T2_ex = term(0.25, [], [t2_tens4,  creOp(ind3), creOp(ind4), desOp(ind2), desOp(ind1)])
         T_othr.append(T2_ex)
         t2_tens5 = tensor('t', [ind5,ind6,ind4,ind3],t2_sym)
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
# Include all type of excitations and de-excitations
     T.extend(T_CE)
     T.extend(T_CA)
     T.extend(T_AE)
     T.extend(T_othr)
#
# for t in T:
#     print 't ampli=', t
#
 return T
#
#####################################
#
def Vperturbation(V):
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
# coreInd = list('ijklmn')
# actvInd = list('xyzwuv')
# virtInd = list('abcdef')
#
 list1 = ['c%i' %p for p in range(30)]
 list2 = ['a%i' %p for p in range(30)]
 list3 = ['v%i' %p for p in range(30)]
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
 ten7 =  tensor('V', [cor2, act2, cor1, act1], v2sym)
# V.append( term(1.0, [], [ten7, desOp(cor1), creOp(cor2), creOp(act2), desOp(act1)]))
#
 ten8 =  tensor('gamma', [act2, act1], d1sym)
 V.append( term(1.0, [], [ten7, ten8,  desOp(cor2), creOp(cor1)]))
#
 ten9 =  tensor('V', [vir2, act2, vir1, act1], v2sym)
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
#                                       vTen = tensor('V', [r,s,p,q], v2sym)
                                       vTen = tensor('V', [r,s,p,q], v2sym)
                                       V81.append(term(-0.25, [], [vTen,desOp(r), desOp(s), creOp(p), creOp(q)]))
#
                                   else:
                                       vTen = tensor('V', [r,s,p,q], v2sym)
#                                       vTen = tensor('V', [p,q,r,s], v2sym)
                                       V81.append(term(0.25, [], [vTen, creOp(p), creOp(q),desOp(s), desOp(r)]))
#
 
 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list3.pop(0), [tg_v], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list3.pop(0), [tg_v], dummy)
 vTen = tensor('V', [r,s,p,q], v2sym)
 V81.append(term(-1.0, [], [vTen, creOp(p), creOp(q), desOp(s), desOp(r)]))

 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list2.pop(0), [tg_a], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list2.pop(0), [tg_a], dummy)
 vTen = tensor('V', [r,s,p,q], v2sym)
 V81.append(term(-1.0, [], [vTen, creOp(p), creOp(q), desOp(s), desOp(r)]))

 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list3.pop(0), [tg_v], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list3.pop(0), [tg_v], dummy)
 vTen = tensor('V', [r,s,p,q], v2sym)
 V81.append(term(1.0, [], [vTen, creOp(q), desOp(r), desOp(s), creOp(p)]))

 p = index(list1.pop(0), [tg_c], dummy)
 q = index(list2.pop(0), [tg_a], dummy)
 s = index(list1.pop(0), [tg_c], dummy)
 r = index(list2.pop(0), [tg_a], dummy)
 vTen = tensor('V', [r,s,p,q], v2sym)
 V81.append(term(1.0, [], [vTen, desOp(s), creOp(p), creOp(q), desOp(r)]))

 V.extend(V81)
# Dummy indices label upate
# dummyLabel(V)
# print len(V81)
 return V
#
#####################################
#
def generateEinsum(terms, ind_str, idn_trans = False, optimize = True):
#
# summary: Generate Einsum structures for each term. 
#          terms   : A list of all terms.
#          ind_str : Indices of the matrix (string).
#
# (c) 2018-2019 Koushik Chatterjee (koushikchatterjee7@gmail.com)
#
 ind1 = ''
 ind2 = ''
 indList = list(ind_str)
 for i in range(len(ind_str)):
     if (i <= float(len(indList)%2)):
        ind0 = indList.pop(0)
        ind1 += ind0
     else:
        ind0 = indList.pop(0)
        ind2 += ind0
#
 if not (len(ind_str)%2 == 0):
    idn_trans = True
#
 term1st = 0
 for term in terms:
     outputS = []
     outputF = []
     OpsList = []
     tens_name = []
#
     OpsindStr = 'rdm_'
     if (idn_trans):
        OpsindStr = 'trdm_'
#
     for i in range(len(term.tensors)):
         TensStr = ''
         tens = term.tensors[i]
         if not (isinstance(tens, creOp) or isinstance(tens,desOp)):
#
            tens_name.append(tens.name)
#
            if ((tens.name == 'v') or (tens.name == 'V')):
               indStr = 'v_'
            elif (tens.name == 'h'):
               indStr = 'h_'
            elif (tens.name == 't'):
               indStr = 't_'
#
            for tens_ind in range(len(tens.indices)): 
               TensStr += str(tens.indices[tens_ind].name)
#
               if not (tens.indices[tens_ind].indType[0][0] == 'virtual'):
                  indStr += str(tens.indices[tens_ind].indType[0][0][0])
               else:
                  indStr += 'e'
#
               if ((tens.name == 'E') or (tens.name == 'e')):
                  if (tens.indices[tens_ind].indType[0][0] == 'core'):
                     indStr = 'e_core_so'
                  elif (tens.indices[tens_ind].indType[0][0] == 'virtual'):
                     indStr = 'e_extern_so'
                  else:
                     raise Exception('Unknown active orbitals energy.')
            if not (tens.name == 't'):
               indStr += '_so'
#
            outputF.append(indStr)   
#
            outputS.append(TensStr)
         else:
               OpsList.append(tens.indices[0].name)
               if (tens.name == 'cre'): 
                   OpsindStr += 'c'
               elif (tens.name == 'des'):
                   OpsindStr += 'a'
# 
     if (len(OpsList)>0):
            OpsStr = ''
            for i in OpsList:
                OpsStr += str(i)
            outputS.append(OpsStr)
            OpsindStr += "_so[1,:]"
            outputF.append(OpsindStr)
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
        cons = ''
        if not (term.numConstant == 1.0):
           if (term.numConstant == -1.0):
               cons = '-'
           else:
               cons = ' '+str(term.numConstant)+' *'
#
     IOpt = 'optimize = True'
     if not (optimize):
        IOpt = 'optimize = False'
#
     print "M[s"+ind1+":f"+ind1+", s"+ind2+":f"+ind2+"] "+sign+"="+cons+" np.einsum('"+str(outputS).translate(None, "'")[1:-1]+"->"+ind_str+"', "+str(outputF).translate(None, "'")[1:-1]+","+IOpt+")"
#
     term1st += 1
#
#####################################
