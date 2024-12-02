#import sys
#from sqaIndex import index, get_spatial_index_type
#from sqaTensor import tensor, creOp, desOp, creDesTensor
#from sqaTerm import term
#from sqaOptions import options
#from sqaSymmetry import symmetry
#from sqaCommutator import commutator
#from sqaIndexList import create_dummy_indices_list

import sys
import sqa_plus
from sqa_plus.sqaIndex import index, get_spatial_index_type, is_core_index_type, is_active_index_type, is_virtual_index_type
from sqa_plus.sqaTensor import tensor, creOp, desOp, creDesTensor
from sqa_plus.sqaTerm import term
from sqa_plus.sqaSymmetry import symmetry
from sqa_plus.sqaCommutator import commutator
from sqa_plus.sqaIndexList import indexLists
from sqa_plus.sqaOptions import options


def Tamplitude(order, indices_lists):

    cor_alpha_inds, vir_alpha_inds, cor_beta_inds, vir_beta_inds = indices_lists
    # Define t amplitude according to their order
    if (order == 1):
        tname = 't1'
        tname_dex = 't1d'
    elif (order == 2):
        tname = 't2'
        tname_dex = 't2d'
    elif (order == 3):
        tname = 't3'

    T = []

    # Define tensors symmetries
    t1_ten_symm  = [symmetry((1,0), 1)]
#    t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]
    t2_ten_symm_ppqq = []

    if (order != 1):
         # Core-External: t_i^a
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()

        t1_ten = tensor(tname, [cor_1, vir_2])
        t1_ten_dex = tensor(tname_dex, [cor_1, vir_2])
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten_dex, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)
        
        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()

        t1_ten = tensor(tname, [cor_1, vir_2])
        t1_ten_dex = tensor(tname_dex, [cor_1, vir_2])
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten_dex, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)


    # Doubles
    ## Core-External-Core-External: t_{ij}^{ab}
    cor_1 = cor_alpha_inds.new_index()
    cor_2 = cor_alpha_inds.new_index()
    vir_3 = vir_alpha_inds.new_index()
    vir_4 = vir_alpha_inds.new_index()

    t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
    t2_ten_dex = tensor(tname_dex, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
    T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
    T2_dex = term(-0.25, [], [t2_ten_dex, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
    T.append(T2_ex)
    T.append(T2_dex)

    cor_1 = cor_beta_inds.new_index()
    cor_2 = cor_beta_inds.new_index()
    vir_3 = vir_beta_inds.new_index()
    vir_4 = vir_beta_inds.new_index()

    t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
    t2_ten_dex = tensor(tname_dex, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
    T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
    T2_dex = term(-0.25, [], [t2_ten_dex, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
    T.append(T2_ex)
    T.append(T2_dex)

    cor_1 = cor_alpha_inds.new_index()
    cor_2 = cor_beta_inds.new_index()
    vir_3 = vir_alpha_inds.new_index()
    vir_4 = vir_beta_inds.new_index()

    t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], [])
    t2_ten_dex = tensor(tname_dex, [cor_1, cor_2, vir_3, vir_4], [])
    T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
    T2_dex = term(-1.00, [], [t2_ten_dex, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
    T.append(T2_ex)
    T.append(T2_dex)

    return T



def dyallH(indices_lists):
    "Construct spin-integrated Dyall Hamiltonian operator."
    cor_alpha_inds, vir_alpha_inds, cor_beta_inds, vir_beta_inds = indices_lists
    

    DyallH = []

    # Core and vitual parts : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
    cor_1 = cor_alpha_inds.new_index()
    e_cor_ten = tensor('e', [cor_1], [])
    DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

    cor_1 = cor_beta_inds.new_index()
    e_cor_ten = tensor('e', [cor_1], [])
    DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

    vir_1 = vir_alpha_inds.new_index()
    e_vir_ten = tensor('e', [vir_1], [])
    DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))

    vir_1 = vir_beta_inds.new_index()
    e_vir_ten = tensor('e', [vir_1], [])
    DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))


    return DyallH

def Vperturbation(indices_lists):
    "Construct spin-integrated perturbation operator V."
    cor_alpha_inds, vir_alpha_inds, cor_beta_inds, vir_beta_inds = indices_lists
    
    v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
    v2e_sym_braket = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
    h1e_sym = [symmetry((1,0), 1)]
    
#    v2e_sym = [symmetry((1,0,3,2), 1)]
#    v2e_sym_braket = [symmetry((1,0,3,2), 1)]
#    h1e_sym = [symmetry((1,0), 1)]

    V = []

    # External-Core: <a|h|i> a^{\dag}_a a_i
    cor_1 = cor_alpha_inds.new_index()
    vir_2 = vir_alpha_inds.new_index()
    v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
    V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))
#
    cor_1 = cor_beta_inds.new_index()
    vir_2 = vir_beta_inds.new_index()
    v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
    V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))
#
    # Core-External: <i|h|a> a^{\dag}_i a_a
    vir_1 = vir_alpha_inds.new_index()
    cor_2 = cor_alpha_inds.new_index()
    v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
    V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))
#
    vir_1 = vir_beta_inds.new_index()
    cor_2 = cor_beta_inds.new_index()
    v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
    V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))


    # <rs||pq> a^{\dag}_p a^{\dag}_q a_s a_r
    tens_ind_cor = [options.core_type,    options.core_type,    options.core_type,    options.core_type]
    tens_ind_vir = [options.virtual_type, options.virtual_type, options.virtual_type, options.virtual_type]

    for ind_type_1 in (cor_alpha_inds, vir_alpha_inds):
        ind_1 = ind_type_1.new_index()

        for ind_type_2 in (cor_alpha_inds, vir_alpha_inds):
            ind_2 = ind_type_2.new_index()

            for ind_type_3 in (cor_alpha_inds, vir_alpha_inds):
                ind_3 = ind_type_3.new_index()

                for ind_type_4 in (cor_alpha_inds, vir_alpha_inds):
                    ind_4 = ind_type_4.new_index()

                    tens_ind_type = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                    if tens_ind_type == tens_ind_cor:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                        V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                    elif tens_ind_type == tens_ind_vir:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                        V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                    else:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                        V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

    for ind_type_1 in (cor_beta_inds, vir_beta_inds):
        ind_1 = ind_type_1.new_index()

        for ind_type_2 in (cor_beta_inds, vir_beta_inds):
            ind_2 = ind_type_2.new_index()

            for ind_type_3 in (cor_beta_inds, vir_beta_inds):
                ind_3 = ind_type_3.new_index()

                for ind_type_4 in (cor_beta_inds, vir_beta_inds):
                    ind_4 = ind_type_4.new_index()

                    tens_ind_type = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                    if tens_ind_type == tens_ind_cor:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                        V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                    elif tens_ind_type == tens_ind_vir:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                        V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                    else:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                        V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

    for ind_type_1 in (cor_alpha_inds, vir_alpha_inds):
        ind_1 = ind_type_1.new_index()

        for ind_type_2 in (cor_beta_inds, vir_beta_inds):
            ind_2 = ind_type_2.new_index()

            for ind_type_3 in (cor_alpha_inds, vir_alpha_inds):
                ind_3 = ind_type_3.new_index()

                for ind_type_4 in (cor_beta_inds, vir_beta_inds):
                    ind_4 = ind_type_4.new_index()

                    tens_ind_type = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                    if tens_ind_type == tens_ind_cor:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                        V.append(term(-1.00, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                    elif tens_ind_type == tens_ind_vir:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                        V.append(term(1.00, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                    else:
                        v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                        V.append(term(1.00, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))


    # External-Core-Core-External: <bj||ia> a^{\dag}_i a^{\dag}_a a_j a_b
    cor_1 = cor_alpha_inds.new_index()
    vir_2 = vir_alpha_inds.new_index()
    vir_3 = vir_alpha_inds.new_index()
    cor_4 = cor_alpha_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

    cor_1 = cor_beta_inds.new_index()
    vir_2 = vir_beta_inds.new_index()
    vir_3 = vir_beta_inds.new_index()
    cor_4 = cor_beta_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

    cor_1 = cor_alpha_inds.new_index()
    vir_2 = vir_beta_inds.new_index()
    vir_3 = vir_beta_inds.new_index()
    cor_4 = cor_alpha_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

    cor_1 = cor_beta_inds.new_index()
    vir_2 = vir_alpha_inds.new_index()
    vir_3 = vir_alpha_inds.new_index()
    cor_4 = cor_beta_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

    # External-Core-Core-External: <bj||ia> a^{\dag}_a a_b a_j a^{\dag}_i
    cor_1 = cor_alpha_inds.new_index()
    vir_2 = vir_alpha_inds.new_index()
    vir_3 = vir_alpha_inds.new_index()
    cor_4 = cor_alpha_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

    cor_1 = cor_beta_inds.new_index()
    vir_2 = vir_beta_inds.new_index()
    vir_3 = vir_beta_inds.new_index()
    cor_4 = cor_beta_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

    cor_1 = cor_alpha_inds.new_index()
    vir_2 = vir_beta_inds.new_index()
    vir_3 = vir_beta_inds.new_index()
    cor_4 = cor_alpha_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

    cor_1 = cor_beta_inds.new_index()
    vir_2 = vir_alpha_inds.new_index()
    vir_3 = vir_alpha_inds.new_index()
    cor_4 = cor_beta_inds.new_index()
    v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
    V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

    return V

def Heff(order=None, op = None):
    "Construct effective Hamiltonian(L)."
    indices_lists = indexLists.core_alpha, indexLists.virtual_alpha, indexLists.core_beta, indexLists.virtual_beta


#    print("--------------------------------- Hamiltonian({:}) ---------------------------------".format(order))
#    sys.stdout.flush()
#    #   order = 0 : L(0) = H(0)
#    #   order = 1 : L(1) = V + [H(0),T(1) - T'(1)]
#    #   order = 2 : L(2) = [H(0),(T(2) - T'(2))] + 1/2[(V + L(1)), (T(1) - T'(1))]
#
#
#    if (order == 0):
#        # L(0) = H(0)
#        L = dyallH(indices_lists)
#
#    elif (order == 1):
#        # L(1) = V + [H(0),T(1) - T'(1)]
#        L = []
#
#        effH = dyallH(indices_lists)
#
#        V = Vperturbation(indices_lists)
#
#        L.extend(V)
#
#        T1 = Tamplitude(1, indices_lists)
#
#        com1 = commutator(effH, T1)
#        print ("Commutation: Done ...")
#        sys.stdout.flush()
#
#        L.extend(com1)
#
#    elif (order == 2):
#        # L(2) = [H(0),T(2) - T'(2)]+ 1/2 [V + L(1),T(1) - T'(1)]
#        L = []
#
#        effH = dyallH(indices_lists)
#
#        T2 = Tamplitude(2, indices_lists)
#
#        com1 = commutator(effH, T2)
#        print ("First Commutation: Done ...")
#        sys.stdout.flush()
#
#        L.extend(com1)
#
#        V = Vperturbation(indices_lists)
#
#        T1 = Tamplitude(1, indices_lists)
#
#        com2 = commutator(V, T1)
#        print ("Second Commutation: Done ...")
#        sys.stdout.flush()
#
#        L.extend(com2)
#
#        effH = dyallH(indices_lists)
#
#        T1_1 = Tamplitude(1, indices_lists)
#
#        com3 = commutator(effH, T1_1)
#
#        T1_2 = Tamplitude(1, indices_lists)
#
#        com4 = commutator(com3, T1_2)
#        print ("Third Commutation: Done ...")
#        sys.stdout.flush()
#
#        for t in com4:
#            t.scale(0.5)
#        L.extend(com4)
#
#    elif (order == 3):    #   L(2) = [H(0),T(3) - T'(3)]+ [V ,T(2) - T'(2)] + 1/2 * [[H(0) ,T(1) - T'(1)],T(2) - T'(2)] 
#                          #          + 1/2 * [[H(0) ,T(2) - T'(2)],T(1) - T'(1)] + 1/2 * [[V ,T(1) - T'(1)],T(1) - T'(1)]
#                          #          + 1/6 * [[[H(0),T(1) - T'(1)] ,T(1) - T'(1)],T(1) - T'(1)]
#        L = []
#
#        effH = dyallH(indices_lists)
#
#        T3 = Tamplitude(3, indices_lists)
#
#        com1 = commutator(effH, T3)
#        print ("First Commutation: Done ...")
#        sys.stdout.flush()
#
#        L.extend(com1)
#
#        V = Vperturbation(indices_lists)
#
#        T2 = Tamplitude(2, indices_lists)
#
#        com2 = commutator(V, T2)
#        print ("Second Commutation: Done ...")
#        sys.stdout.flush()
#
#        L.extend(com2)
#
#        effH = dyallH(indices_lists)
#
#        T1_1 = Tamplitude(1, indices_lists)
#
#        com3 = commutator(effH, T1_1)
#
#        T2_2 = Tamplitude(2, indices_lists)
#
#        com4 = commutator(com3, T2_2)
#        print ("Third Commutation: Done ...")
#        sys.stdout.flush()
#
#        for t in com4:
#            t.scale(0.5)
#        L.extend(com4)
#
#        effH = dyallH(indices_lists)
#
#        T2_2_0 = Tamplitude(2, indices_lists)
#
#        com5 = commutator(effH, T2_2_0)
#
#        T1_1_0 = Tamplitude(1, indices_lists)
#
#        com6 = commutator(com5, T1_1_0)
#        print ("Fourth Commutation: Done ...")
#        sys.stdout.flush()
#
#        for t in com6:
#            t.scale(0.5)
#        L.extend(com6)
#
#        V = Vperturbation(indices_lists)
#
#        T1_1_1 = Tamplitude(1, indices_lists)
#
#        com7 = commutator(V, T1_1_1)
#
#        T1_1_2 = Tamplitude(1, indices_lists)
#
#        com8 = commutator(com7, T1_1_2)
#        print ("Fith Commutation: Done ...")
#        sys.stdout.flush()
#
#        for t in com8:
#            t.scale(0.5)
#        L.extend(com8)
#
#        effH = dyallH(indices_lists)
#
#        T1_1_3 = Tamplitude(1, indices_lists)
#
#        com9 = commutator(effH, T1_1_3)
#
#        T1_1_4 = Tamplitude(1, indices_lists)
#
#        com10 = commutator(com9, T1_1_4)
#
#        T1_1_5 = Tamplitude(1, indices_lists)
#
#        com11 = commutator(com10, T1_1_5)
#
#        print ("Sixth Commutation: Done ...")
#        sys.stdout.flush()
#
#        for t in com11:
#            t.scale(0.166666667)
#            #t.scale(1/6)
#        L.extend(com11)
################################################################################################################################

####Below is for bch expansions of general operators that are not Heff#####
    if (order == 0):
        # L(0) = H(0)
        L = op

    elif (order == 1):
        # L(1) = [H(0),T(1) - T'(1)]
        L = []

        T1 = Tamplitude(1, indices_lists)

        com1 = commutator(op, T1)
        print ("Commutation: Done ...")
        sys.stdout.flush()

        L.extend(com1)

    elif (order == 2):
        # L(2) = [H(0),T(2) - T'(2)]+ 1/2 [L(1),T(1) - T'(1)]
        L = []

        effH = op

        T2 = Tamplitude(2, indices_lists)

        com1 = commutator(effH, T2)
        print ("First Commutation: Done ...")
        sys.stdout.flush()

        L.extend(com1)

        T1 = Tamplitude(1, indices_lists)

        effH = op

        T1_1 = Tamplitude(1, indices_lists)

        com3 = commutator(effH, T1_1)

        T1_2 = Tamplitude(1, indices_lists)

        com4 = commutator(com3, T1_2)
        print ("Third Commutation: Done ...")
        sys.stdout.flush()

        for t in com4:
            t.scale(0.5)
        L.extend(com4)

    #if order is None:
    elif (order == 3):    #   L(2) = [H(0),T(3) - T'(3)]+ 1/2 * [[H(0) ,T(1) - T'(1)],T(2) - T'(2)] 
                          #          + 1/2 * [[H(0) ,T(2) - T'(2)],T(1) - T'(1)] + 1/2
                          #          + 1/6 * [[[H(0),T(1) - T'(1)] ,T(1) - T'(1)],T(1) - T'(1)]
        L = []

        T3 = Tamplitude(3, indices_lists)

        com1 = commutator(op, T3)
        print ("First Commutation: Done ...")
        sys.stdout.flush()

        L.extend(com1)
##################################yes to above

        T1_1 = Tamplitude(1, indices_lists)

        com3 = commutator(op, T1_1)

        T2_2 = Tamplitude(2, indices_lists)

        com4 = commutator(com3, T2_2)
        print ("second Commutation: Done ...")
        sys.stdout.flush()

        for t in com4:
            t.scale(0.5)
        L.extend(com4)
##################################yes to above

        T2_2_0 = Tamplitude(2, indices_lists)

        com5 = commutator(op, T2_2_0)

        T1_1_0 = Tamplitude(1, indices_lists)

        com6 = commutator(com5, T1_1_0)
        print ("third Commutation: Done ...")
        sys.stdout.flush()

        for t in com6:
            t.scale(0.5)
        L.extend(com6)
############################yes to the above
        T1_1_3 = Tamplitude(1, indices_lists)

        com9 = commutator(op, T1_1_3)

        T1_1_4 = Tamplitude(1, indices_lists)

        com10 = commutator(com9, T1_1_4)

        T1_1_5 = Tamplitude(1, indices_lists)

        com11 = commutator(com10, T1_1_5)

        print ("fourth Commutation: Done ...")
        sys.stdout.flush()

        for t in com11:
            t.scale(0.166666667)
            #t.scale(1/6)
        L.extend(com11)
####################yes to the above
    else:
        raise Exception ('Unknown type of effective Hamiltonian of order = %s' % (order))

    print ("Done ...")
    print ("""--------------------------------------------------------------""")
    sys.stdout.flush()



    return L
