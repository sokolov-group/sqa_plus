# Copyright 2009-2022 SecondQuantizationAlgebra Developers. All Rights Reserved.
#
# Licensed under the GNU General Public License v3.0;
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# In addition, any modification or use of this software should
# cite the following paper:
#
#   E. Neuscamman, T. Yanai, and G. K.-L. Chan.
#   J. Chem. Phys. 130, 124102 (2009)
#
# Authors: Koushik Chatterjee <koushikchatterjee7@gmail.com>
#          Carlos E. V. de Moura <carlosevmoura@gmail.com>
#

import sys

from sqaIndex import index, get_spatial_index_type, is_core_index_type, is_active_index_type, is_virtual_index_type
from sqaTensor import tensor, creOp, desOp, creDesTensor
from sqaTerm import term
from sqaOptions import options
from sqaSymmetry import symmetry
from sqaCommutator import commutator
from sqaIndexList import create_dummy_indices_list

def Heff(order = 0, spin_integrated = False, explicit_spin_cases = True, cvs_indices = False, internal_excit = True):
    "Construct effective Hamiltonian(L)."

    print("--------------------------------- Hamiltonian({:}) ---------------------------------".format(order))
    sys.stdout.flush()
    #   order = 0 : L(0) = H(0)
    #   order = 1 : L(1) = V + [H(0),T(1) - T'(1)]
    #   order = 2 : L(2) = [H(0),(T(2) - T'(2))] + 1/2[(V + L(1)), (T(1) - T'(1))]

    indices_lists = create_dummy_indices_list(spin_integrated, cvs_indices)

    if (order == 0):
        # L(0) = H(0)
        L = dyallH(indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

    elif (order == 1):
        # L(1) = V + [H(0),T(1) - T'(1)]
        L = []

        effH = dyallH(indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        V = Vperturbation(indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        L.extend(V)

        T1 = Tamplitude(1, indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        com1 = commutator(effH, T1)
        print("Commutation: Done ...")
        sys.stdout.flush()

        L.extend(com1)

    elif (order == 2):
        # L(2) = [H(0),T(2) - T'(2)]+ 1/2 [V + L(1),T(1) - T'(1)]
        L = []

        effH = dyallH(indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        T2 = Tamplitude(2, indices_lists, spin_integrated, explicit_spin_cases, internal_excit, cvs_indices)

        com1 = commutator(effH, T2)
        print("First Commutation: Done ...")
        sys.stdout.flush()

        L.extend(com1)

        V = Vperturbation(indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        T1 = Tamplitude(1, indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        com2 = commutator(V, T1)
        print("Second Commutation: Done ...")
        sys.stdout.flush()

        L.extend(com2)

        effH = dyallH(indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        T1_1 = Tamplitude(1, indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        com3 = commutator(effH, T1_1)

        T1_2 = Tamplitude(1, indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

        com4 = commutator(com3, T1_2)
        print("Third Commutation: Done ...")
        sys.stdout.flush()

        for t in com4:
            t.scale(0.5)
        L.extend(com4)

    else:
        raise Exception ('Unknown type of effective Hamiltonian of order = %s' % (order))

    print("Done ...")
    print("----------------------------------------------------------------------------------")
    sys.stdout.flush()
    return L

def dyallH(indices_lists, spin_integrated = False, explicit_spin_cases = False, cvs_indices = False):

    def dyallH_spin_orbital(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cor_inds, act_inds, vir_inds = indices_lists

        DyallH = []
        # E_fc
        c = index('Const.', [], True)
        Efc = tensor('E_fc',[c], [])
        DyallH.append(term(1.0, [], [Efc]))

        # Core and vitual parts : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
        cor_1 = cor_inds.new_index()
        e_cor_ten = tensor('e', [cor_1], [])
        DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        vir_1 = vir_inds.new_index()
        e_vir_ten = tensor('e', [vir_1], [])
        DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))

        # Active part : DyallH_act
        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_inds.new_index()
        cor_2 = cor_inds.new_index()
        act_3 = act_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        DyallH.extend(DyallH_act)

        return DyallH

    def dyallH_spin_integrated(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH = []
        # E_fc
        c = index('Const.', [], True)
        Efc = tensor('E_fc',[c], [])
        DyallH.append(term(1.0, [], [Efc]))

        # Core and vitual parts : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
        for cor_inds in [cor_alpha_inds, cor_beta_inds]:
            cor_1 = cor_inds.new_index()
            e_cor_ten = tensor('e', [cor_1], [])
            DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        for vir_inds in [vir_alpha_inds, vir_beta_inds]:
            vir_1 = vir_inds.new_index()
            e_vir_ten = tensor('e', [vir_1], [])
            DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))

        # Active part : DyallH_act
        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        for act_inds in [act_alpha_inds, act_beta_inds]:
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
            DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        for act_inds in [act_alpha_inds, act_beta_inds]:
            for cor_inds in [cor_alpha_inds, cor_beta_inds]:
                act_1 = act_inds.new_index()
                cor_2 = cor_inds.new_index()
                act_3 = act_inds.new_index()
                v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
                DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        for act_inds in [act_alpha_inds, act_beta_inds]:
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
            DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(1.00, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        DyallH.extend(DyallH_act)

        return DyallH

    def dyallH_spin_integrated_explicit_cases(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH = []
        # E_fc
        c = index('Const.', [], True)
        Efc = tensor('E_fc',[c], [])
        DyallH.append(term(1.0, [], [Efc]))

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

        # Active part : DyallH_act
        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        DyallH.extend(DyallH_act)

        return DyallH

    def dyallH_cvs_spin_orbital(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cvs_cor_inds, cvs_val_inds, act_inds, vir_inds = indices_lists

        DyallH = []
        # E_fc
        c = index('Const.', [], True)
        Efc = tensor('E_fc',[c], [])
        DyallH.append(term(1.0, [], [Efc]))

        # Core and vitual parts : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            e_cor_ten = tensor('e', [cor_1], [])
            DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        vir_1 = vir_inds.new_index()
        e_vir_ten = tensor('e', [vir_1], [])
        DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))

        # Active part : DyallH_act
        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            act_1 = act_inds.new_index()
            cor_2 = cvs_inds.new_index()
            act_3 = act_inds.new_index()
            v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
            DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        DyallH.extend(DyallH_act)

        return DyallH

    def dyallH_cvs_spin_integrated(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH = []
        # E_fc
        c = index('Const.', [], True)
        Efc = tensor('E_fc',[c], [])
        DyallH.append(term(1.0, [], [Efc]))

        # Core and vitual parts : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
        for cvs_inds in [cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds]:
            cor_1 = cvs_inds.new_index()
            e_cor_ten = tensor('e', [cor_1], [])
            DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        for vir_inds in [vir_alpha_inds, vir_beta_inds]:
            vir_1 = vir_inds.new_index()
            e_vir_ten = tensor('e', [vir_1], [])
            DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))

        # Active part : DyallH_act
        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        for act_inds in [act_alpha_inds, act_beta_inds]:
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
            DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        for cvs_inds in [cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds]:
            for act_inds in [act_alpha_inds, act_beta_inds]:
                act_1 = act_alpha_inds.new_index()
                cor_2 = cvs_inds.new_index()
                act_3 = act_alpha_inds.new_index()
                v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
                DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        for act_inds in [act_alpha_inds, act_beta_inds]:
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
            DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(1.00, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        DyallH.extend(DyallH_act)

        return DyallH

    def dyallH_cvs_spin_integrated_explicit_cases(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH = []
        # E_fc
        c = index('Const.', [], True)
        Efc = tensor('E_fc',[c], [])
        DyallH.append(term(1.0, [], [Efc]))

        # Core and vitual parts : SUM_i E_i {a_i a^+_i} + SUM_a E_a {a^+_a a_a}
        cor_1 = cvs_cor_alpha_inds.new_index()
        e_cor_ten = tensor('e', [cor_1], [])
        DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        e_cor_ten = tensor('e', [cor_1], [])
        DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        e_cor_ten = tensor('e', [cor_1], [])
        DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        e_cor_ten = tensor('e', [cor_1], [])
        DyallH.append(term(-1.0, [],[e_cor_ten, desOp(cor_1), creOp(cor_1)]))

        vir_1 = vir_alpha_inds.new_index()
        e_vir_ten = tensor('e', [vir_1], [])
        DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        e_vir_ten = tensor('e', [vir_1], [])
        DyallH.append(term( 1.0, [],[e_vir_ten, creOp(vir_1), desOp(vir_1)]))

        # Active part : DyallH_act
        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        DyallH.extend(DyallH_act)

        return DyallH

    if indices_lists is None:
        raise Exception("List of dummy indices should be specified.")

    if spin_integrated:
        if explicit_spin_cases:
            if cvs_indices:
                dyallH = dyallH_cvs_spin_integrated_explicit_cases(indices_lists)
            else:
                dyallH = dyallH_spin_integrated_explicit_cases(indices_lists)
        else:
            if cvs_indices:
                dyallH = dyallH_cvs_spin_integrated(indices_lists)
            else:
                dyallH = dyallH_spin_integrated(indices_lists)
    else:
        if cvs_indices:
            dyallH = dyallH_cvs_spin_orbital(indices_lists)
        else:
            dyallH = dyallH_spin_orbital(indices_lists)

    return dyallH

def dyallH_act(indices_lists, spin_integrated = False, explicit_spin_cases = False, cvs_indices = False):

    def dyallH_act_spin_orbital(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cor_inds, act_inds, vir_inds = indices_lists

        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_inds.new_index()
        cor_2 = cor_inds.new_index()
        act_3 = act_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        return DyallH_act

    def dyallH_act_spin_integrated(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        for act_inds in [act_alpha_inds, act_beta_inds]:
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
            DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        for act_inds in [act_alpha_inds, act_beta_inds]:
            for cor_inds in [cor_alpha_inds, cor_beta_inds]:
                act_1 = act_inds.new_index()
                cor_2 = cor_inds.new_index()
                act_3 = act_inds.new_index()
                v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
                DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        for act_inds in [act_alpha_inds, act_beta_inds]:
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
            DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(1.00, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        return DyallH_act

    def dyallH_act_spin_integrated_explicit_cases(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        return DyallH_act

    def dyallH_act_cvs_spin_orbital(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cvs_cor_inds, cvs_val_inds, act_inds, vir_inds = indices_lists

        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            act_1 = act_inds.new_index()
            cor_2 = cvs_inds.new_index()
            act_3 = act_inds.new_index()
            v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
            DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        return DyallH_act

    def dyallH_act_cvs_spin_integrated(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        for act_inds in [act_alpha_inds, act_beta_inds]:
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
            DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        for cvs_inds in [cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds]:
            for act_inds in [act_alpha_inds, act_beta_inds]:
                act_1 = act_alpha_inds.new_index()
                cor_2 = cvs_inds.new_index()
                act_3 = act_alpha_inds.new_index()
                v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
                DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(1.00, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        return DyallH_act

    def dyallH_act_cvs_spin_integrated_explicit_cases(indices_lists):
        "Construct spin-integrated Dyall Hamiltonian operator."

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        DyallH_act = []

        h1e_sym = [symmetry((1,0), 1)]
        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        h1e_ten = tensor('h', [act_2, act_1], h1e_sym)
        DyallH_act.append(term(1.0, [], [h1e_ten, creOp(act_1), desOp(act_2)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, cor_2, act_1, cor_2], v2e_sym)
        DyallH_act.append(term(1.0, [], [v2e_ten, creOp(act_3), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        v2e_ten = tensor('v', [act_3, act_4, act_1, act_2], v2e_sym)
        DyallH_act.append(term(0.25, [], [v2e_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(act_3)]))

        return DyallH_act

    if indices_lists is None:
        raise Exception("List of dummy indices should be specified.")

    if spin_integrated:
        if explicit_spin_cases:
            if cvs_indices:
                dyallH_act = dyallH_act_cvs_spin_integrated_explicit_cases(indices_lists)
            else:
                dyallH_act = dyallH_act_spin_integrated_explicit_cases(indices_lists)
        else:
            if cvs_indices:
                dyallH_act = dyallH_act_cvs_spin_integrated(indices_lists)
            else:
                dyallH_act = dyallH_act_spin_integrated(indices_lists)
    else:
        if cvs_indices:
            dyallH_act = dyallH_act_cvs_spin_orbital(indices_lists)
        else:
            dyallH_act = dyallH_act_spin_orbital(indices_lists)

    return dyallH_act

def Tamplitude(order = 1, indices_lists = None, spin_integrated = False, explicit_spin_cases = False, cvs_indices = False, internal_excit = True):
    # Cluster operator  : T - T^dag, Where T = T1 + T2
    # Single excitation : T1

    def Tamplitude_spin_orbital(order, indices_lists, internal_excit):

        cor_inds, act_inds, vir_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        if (order == 2 and internal_excit):
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            t1_ten = tensor(tname, [act_1, act_2], t1_ten_symm)
            T1_ex  = term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    def Tamplitude_spin_integrated(order, indices_lists, internal_excit):

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-1.00, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        if (order == 2 and internal_excit):
            t1_asym = [symmetry((1,0), -1)]
            act_1 = act_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

            act_1 = act_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    def Tamplitude_spin_integrated_explicit_cases(order, indices_lists, internal_excit):

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        if (order == 2 and internal_excit):
            t1_asym = [symmetry((1,0), -1)]
            act_1 = act_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

            act_1 = act_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    def Tamplitude_cvs_spin_orbital(order, indices_lists, internal_excit):

        cvs_cor_inds, cvs_val_inds, act_inds, vir_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            vir_2 = vir_inds.new_index()
            t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
            T.append(T1_ex)
            T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            vir_3 = vir_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        # Core-Active: t_i^x
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
            T.append(T1_ex)
            T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            vir_3 = vir_inds.new_index()
            vir_4 = vir_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            vir_3 = vir_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            vir_3 = vir_inds.new_index()
            vir_4 = vir_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        if (order == 2 and internal_excit):
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            t1_ten = tensor(tname, [act_1, act_2], t1_ten_symm)
            T1_ex  = term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    def Tamplitude_cvs_spin_integrated(order, indices_lists, internal_excit):

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm  = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            vir_2 = vir_alpha_inds.new_index()
            t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
            T.append(T1_ex)
            T.append(T1_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            vir_2 = vir_beta_inds.new_index()
            t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
            T.append(T1_ex)
            T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            vir_3 = vir_alpha_inds.new_index()
            act_4 = act_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_beta_inds.new_index()
            vir_3 = vir_alpha_inds.new_index()
            act_4 = act_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_beta_inds.new_index()
            vir_3 = vir_beta_inds.new_index()
            act_4 = act_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            vir_3 = vir_beta_inds.new_index()
            act_4 = act_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            vir_3 = vir_beta_inds.new_index()
            act_4 = act_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            vir_3 = vir_alpha_inds.new_index()
            act_4 = act_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        # Core-Active: t_i^x
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
            T.append(T1_ex)
            T.append(T1_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
            T.append(T1_ex)
            T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            act_3 = act_alpha_inds.new_index()
            act_4 = act_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_beta_inds.new_index()
            act_3 = act_alpha_inds.new_index()
            act_4 = act_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
            T2_ex  = term( 1.0, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            act_3 = act_beta_inds.new_index()
            act_4 = act_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            act_3 = act_beta_inds.new_index()
            act_4 = act_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
            T2_ex  = term( 1.0, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            cor_2 = cvs_alpha_inds.new_index()
            vir_3 = vir_alpha_inds.new_index()
            vir_4 = vir_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            cor_2 = cvs_beta_inds.new_index()
            vir_3 = vir_beta_inds.new_index()
            vir_4 = vir_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.50, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.50, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.50, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.50, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            cor_2 = cvs_alpha_inds.new_index()
            vir_3 = vir_alpha_inds.new_index()
            act_4 = act_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            cor_2 = cvs_beta_inds.new_index()
            vir_3 = vir_beta_inds.new_index()
            act_4 = act_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
                cor_1 = cvs_alpha_inds.new_index()
                cor_2 = cvs_beta_inds.new_index()
                vir_3 = vir_alpha_inds.new_index()
                act_4 = act_beta_inds.new_index()
                t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
                T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
                T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
                T.append(T2_ex)
                T.append(T2_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
                cor_1 = cvs_beta_inds.new_index()
                cor_2 = cvs_alpha_inds.new_index()
                vir_3 = vir_beta_inds.new_index()
                act_4 = act_alpha_inds.new_index()
                t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
                T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
                T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
                T.append(T2_ex)
                T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            vir_3 = vir_alpha_inds.new_index()
            vir_4 = vir_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_beta_inds.new_index()
            vir_3 = vir_alpha_inds.new_index()
            vir_4 = vir_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            vir_3 = vir_beta_inds.new_index()
            vir_4 = vir_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            vir_3 = vir_beta_inds.new_index()
            vir_4 = vir_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            cor_2 = cvs_alpha_inds.new_index()
            act_3 = act_alpha_inds.new_index()
            act_4 = act_alpha_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            cor_2 = cvs_beta_inds.new_index()
            act_3 = act_beta_inds.new_index()
            act_4 = act_beta_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
            T.append(T2_ex)
            T.append(T2_dex)

        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
                cor_1 = cvs_alpha_inds.new_index()
                cor_2 = cvs_beta_inds.new_index()
                act_3 = act_alpha_inds.new_index()
                act_4 = act_beta_inds.new_index()
                t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
                T2_ex  = term( 1.00, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
                T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
                T.append(T2_ex)
                T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-1.00, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        if (order == 2 and internal_excit):
            act_1 = act_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            t1_ten = tensor(tname, [act_1, act_2], t1_ten_symm)
            T1_ex  = term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

            act_1 = act_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            t1_ten = tensor(tname, [act_1, act_2], t1_ten_symm)
            T1_ex  = term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    def Tamplitude_cvs_spin_integrated_explicit_cases(order, indices_lists, internal_excit):

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_ex)
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        if (order == 2 and internal_excit):
            t1_asym = [symmetry((1,0), -1)]
            act_1 = act_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

            act_1 = act_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    if indices_lists is None:
        raise Exception("List of dummy indices should be specified.")

    if spin_integrated:
        if explicit_spin_cases:
            if cvs_indices:
                T = Tamplitude_cvs_spin_integrated_explicit_cases(order, indices_lists, internal_excit)
            else:
                T = Tamplitude_spin_integrated_explicit_cases(order, indices_lists, internal_excit)
        else:
            if cvs_indices:
                T = Tamplitude_cvs_spin_integrated(order, indices_lists, internal_excit)
            else:
                T = Tamplitude_spin_integrated(order, indices_lists, internal_excit)
    else:
        if cvs_indices:
            T = Tamplitude_cvs_spin_orbital(order, indices_lists, internal_excit)
        else:
            T = Tamplitude_spin_orbital(order, indices_lists, internal_excit)

    return T

def Tamplitude_excitation(order = 1, indices_lists = None, spin_integrated = False, explicit_spin_cases = False, cvs_indices = False):

    def Tamplitude_excitation_spin_orbital(order, indices_lists):

        cor_inds, act_inds, vir_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Core-Active: t_i^x
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Active-External: t_x^a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        # Active-External: t_{xy}^{az}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        return T

    def Tamplitude_excitation_spin_integrated(order, indices_lists):

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 4.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Core-Active: t_i^x
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        return T

    def Tamplitude_excitation_spin_integrated_explicit_cases(order, indices_lists):

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Core-Active: t_i^x
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        return T

    def Tamplitude_excitation_cvs_spin_orbital(order, indices_lists, internal_excit):

        cvs_cor_inds, cvs_val_inds, act_inds, vir_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            vir_2 = vir_inds.new_index()
            t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
            T.append(T1_ex)

        # Core-External: t_{ix}^{ay}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            vir_3 = vir_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T.append(T2_ex)

        # Core-Active: t_i^x
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
            T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
            T.append(T1_ex)

        # Core-Active: t_{ix}^{yz}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
            T.append(T2_ex)

        # Active-External: t_x^a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        # Active-External: t_{xy}^{az}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            vir_3 = vir_inds.new_index()
            vir_4 = vir_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
            T.append(T2_ex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            vir_3 = vir_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T.append(T2_ex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Active-External: t_{ix}^{ab}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            vir_3 = vir_inds.new_index()
            vir_4 = vir_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
            T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
            T.append(T2_ex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
            T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
            T.append(T2_ex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        if (order == 2 and internal_excit):
            act_1 = act_inds.new_index()
            act_2 = act_inds.new_index()
            t1_ten = tensor(tname, [act_1, act_2], t1_ten_symm)
            T1_ex  = term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    def Tamplitude_excitation_cvs_spin_integrated(order, indices_lists, internal_excit):

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm  = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Core-Active: t_i^x
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 4.0, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_ex)
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 4.0, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 2.0, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 1.00, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        if (order == 2 and internal_excit):
            act_1 = act_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            t1_ten = tensor(tname, [act_1, act_2], t1_ten_symm)
            T1_ex  = term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

            act_1 = act_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            t1_ten = tensor(tname, [act_1, act_2], t1_ten_symm)
            T1_ex  = term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    def Tamplitude_excitation_cvs_spin_integrated_explicit_cases(order, indices_lists, internal_excit):

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_ex  = term( 1.0, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Core-Active: t_i^x
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(act_2), desOp(cor_1)])
        T.append(T1_ex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_ex  = term( 1.0, [], [t1_ten, creOp(vir_2), desOp(act_1)])
        T.append(T1_ex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_ex  = term( 0.5, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(act_3), creOp(act_4), desOp(cor_2), desOp(cor_1)])
        T.append(T2_ex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_ex  = term( 0.25, [], [t2_ten, creOp(vir_3), creOp(vir_4), desOp(act_2), desOp(act_1)])
        T.append(T2_ex)

        if (order == 2 and internal_excit):
            t1_asym = [symmetry((1,0), -1)]
            act_1 = act_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

            act_1 = act_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            t1_ten =  tensor(tname, [act_1, act_2], t1_asym)
            T1_ex =  term(1.0, [], [t1_ten,  creOp(act_2), desOp(act_1)])
            T.append(T1_ex)

        return T

    if indices_lists is None:
        raise Exception("List of dummy indices should be specified.")

    if spin_integrated:
        if explicit_spin_cases:
            if cvs_indices:
                T = Tamplitude_excitation_cvs_spin_integrated_explicit_cases(order, indices_lists)
            else:
                T = Tamplitude_excitation_spin_integrated_explicit_cases(order, indices_lists)
        else:
            if cvs_indices:
                T = Tamplitude_excitation_cvs_spin_integrated(order, indices_lists)
            else:
                T = Tamplitude_excitation_spin_integrated(order, indices_lists)
    else:
        if cvs_indices:
            T = Tamplitude_excitation_cvs_spin_orbital(order, indices_lists)
        else:
            T = Tamplitude_excitation_spin_orbital(order, indices_lists)

    return T

def Tamplitude_deexcitation(order = 1, indices_lists = None, spin_integrated = False, explicit_spin_cases = False, cvs_indices = False):

    def Tamplitude_deexcitation_spin_orbital(order, indices_lists):

        cor_inds, act_inds, vir_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_inds.new_index()
        cor_2 = cor_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        return T

    #TODO: To be updated
    def Tamplitude_deexcitation_spin_integrated(order, indices_lists):

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(4.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(2.0, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(2.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(1.00, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        return T

    def Tamplitude_deexcitation_spin_integrated_explicit_cases(order, indices_lists):

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_alpha_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cor_beta_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        return T

    def Tamplitude_deexcitation_cvs_spin_orbital(order, indices_lists):

        cvs_cor_inds, cvs_val_inds, act_inds, vir_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            vir_2 = vir_inds.new_index()
            t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
            T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            vir_3 = vir_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
            T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_dex)

        # Core-Active: t_i^x
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
            T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
            T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
            T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            vir_3 = vir_inds.new_index()
            vir_4 = vir_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_dex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            vir_3 = vir_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
            T.append(T2_dex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            vir_3 = vir_inds.new_index()
            vir_4 = vir_inds.new_index()
            t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
            T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
            T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            cor_2 = cvs_inds.new_index()
            act_3 = act_inds.new_index()
            act_4 = act_inds.new_index()
            t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
            T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
            T.append(T2_dex)

        cor_1 = cvs_cor_inds.new_index()
        cor_2 = cvs_val_inds.new_index()
        act_3 = act_inds.new_index()
        act_4 = act_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        vir_4 = vir_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        return T

    def Tamplitude_deexcitation_cvs_spin_integrated(order, indices_lists):

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm  = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-4.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-4.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-2.0, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-2.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-1.00, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        return T

    def Tamplitude_deexcitation_cvs_spin_integrated_explicit_cases(order, indices_lists):

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        # Define t amplitude according to their order
        if (order == 1):
            tname = 't1'
        elif (order == 2):
            tname = 't2'

        T = []

        # Define tensors symmetries
        t1_ten_symm = [symmetry((1,0), 1)]
        t2_ten_symm_ppqq = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        t2_ten_symm_ppqr = [symmetry((1,0,2,3), -1)]
        t2_ten_symm_pqrr = [symmetry((0,1,3,2), -1)]

        # Core-External: t_i^a
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(vir_2)])
        T.append(T1_dex)

        # Core-External: t_{ix}^{ay}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, act_4], [])
        T2_dex = term(-1.0, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Core-Active: t_i^x
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        t1_ten = tensor(tname, [cor_1, act_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(cor_1), desOp(act_2)])
        T.append(T1_dex)

        # Core-Active: t_{ix}^{yz}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, act_3, act_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        # Active-External: t_x^a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        t1_ten = tensor(tname, [act_1, vir_2], t1_ten_symm)
        T1_dex = term(-1.0, [], [t1_ten, creOp(act_1), desOp(vir_2)])
        T.append(T1_dex)

        # Active-External: t_{xy}^{az}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        # Doubles
        ## Core-External-Core-External: t_{ij}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Core-Active: t_{ij}^{ax}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, vir_3, act_4], t2_ten_symm_ppqr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-External-Active-External: t_{ix}^{ab}
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, act_2, vir_3, vir_4], t2_ten_symm_pqrr)
        T2_dex = term(-0.5, [], [t2_ten, creOp(cor_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        ## Core-Active-Core-Active: t_{ij}^{xy}
        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_alpha_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_cor_beta_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_alpha_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        cor_1 = cvs_val_beta_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()
        t2_ten = tensor(tname, [cor_1, cor_2, act_3, act_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(cor_1), creOp(cor_2), desOp(act_4), desOp(act_3)])
        T.append(T2_dex)

        ## Active-External-Active-External: t_{xy}^{ab}
        act_1 = act_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        vir_4 = vir_alpha_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        act_1 = act_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        vir_4 = vir_beta_inds.new_index()
        t2_ten = tensor(tname, [act_1, act_2, vir_3, vir_4], t2_ten_symm_ppqq)
        T2_dex = term(-0.25, [], [t2_ten, creOp(act_1), creOp(act_2), desOp(vir_4), desOp(vir_3)])
        T.append(T2_dex)

        return T

    if indices_lists is None:
        raise Exception("List of dummy indices should be specified.")

    if spin_integrated:
        if explicit_spin_cases:
            if cvs_indices:
                T = Tamplitude_deexcitation_cvs_spin_integrated_explicit_cases(order, indices_lists)
            else:
                T = Tamplitude_deexcitation_spin_integrated_explicit_cases(order, indices_lists)
        else:
            if cvs_indices:
                T = Tamplitude_deexcitation_cvs_spin_integrated(order, indices_lists)
            else:
                T = Tamplitude_deexcitation_spin_integrated(order, indices_lists)
    else:
        if cvs_indices:
            T = Tamplitude_deexcitation_cvs_spin_orbital(order, indices_lists)
        else:
            T = Tamplitude_deexcitation_spin_orbital(order, indices_lists)

    return T

def Vperturbation(indices_lists, spin_integrated = False, explicit_spin_cases = False, cvs_indices = False):

    def Vperturbation_spin_orbital(indices_lists):
        "Construct spin-orbital perturbation operator V."

        cor_inds, act_inds, vir_inds = indices_lists

        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        v2e_sym_braket = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]
        h1e_sym = [symmetry((1,0), 1)]

        V = []

        # Active-Core: <x|h|i> a^{\dag}_x a_i
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        # Core-Active: <i|h|x> a^{\dag}_i a_x
        act_1 = act_inds.new_index()
        cor_2 = cor_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        # Active-External: <x|h|a> a^{\dag}_x a_a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        # External-Active: <a|h|x> a^{\dag}_a a_x
        vir_1 = vir_inds.new_index()
        act_2 = act_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        # External-Core: <a|h|i> a^{\dag}_a a_i
        cor_1 = cor_inds.new_index()
        vir_2 = vir_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        # Core-External: <i|h|a> a^{\dag}_i a_a
        vir_1 = vir_inds.new_index()
        cor_2 = cor_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        # Core-Active-Core-Active: <ix||jy> \gamma^x_y a^{\dag}_i a_j
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        cor_3 = cor_inds.new_index()
        act_4 = act_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        # External-Active-External-Active: <ax||by> \gamma^x_y a^{\dag}_a a_b
        vir_1 = vir_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        # <rs||pq> a^{\dag}_p a^{\dag}_q a_s a_r
        exceptions_list = [[options.active_type,  options.active_type,  options.active_type,  options.active_type]]

        for ind_type in [options.active_type, options.virtual_type]:
            exceptions_list.extend([[options.core_type, ind_type, options.core_type, ind_type],
                                    [options.core_type, ind_type, ind_type, options.core_type],
                                    [ind_type, options.core_type, options.core_type, ind_type],
                                    [ind_type, options.core_type, ind_type, options.core_type]])

        all_cor_types = [[options.core_type, options.core_type, options.core_type, options.core_type]]
        all_vir_types = [[options.virtual_type, options.virtual_type, options.virtual_type, options.virtual_type]]

        for ind_type_1 in (cor_inds, act_inds, vir_inds):
            for ind_type_2 in (cor_inds, act_inds, vir_inds):
                for ind_type_3 in (cor_inds, act_inds, vir_inds):
                    for ind_type_4 in (cor_inds, act_inds, vir_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        # External-Core-Core-External: <bj||ia> a^{\dag}_a a_b a_j a^{\dag}_i
        cor_1 = cor_inds.new_index()
        vir_2 = vir_inds.new_index()
        vir_3 = vir_inds.new_index()
        cor_4 = cor_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        # Active-Core-Core-Active: <yj||ix> a_j a^{\dag}_i a^{\dag}_x a_y
        cor_1 = cor_inds.new_index()
        act_2 = act_inds.new_index()
        act_3 = act_inds.new_index()
        cor_4 = cor_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        sys.stdout.flush()

        return V

    def Vperturbation_spin_integrated(indices_lists):
        "Construct spin-integrated perturbation operator V."

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        v2e_sym_braket = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]
        h1e_sym = [symmetry((1,0), 1)]

        V = []

        # Active-Core: <x|h|i> a^{\dag}_x a_i
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        # Core-Active: <i|h|x> a^{\dag}_i a_x
        act_1 = act_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        # Active-External: <x|h|a> a^{\dag}_x a_a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        # External-Active: <a|h|x> a^{\dag}_a a_x
        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        # External-Core: <a|h|i> a^{\dag}_a a_i
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        # Core-External: <i|h|a> a^{\dag}_i a_a
        vir_1 = vir_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        # Core-Active-Core-Active: <ix||jy> \gamma^x_y a^{\dag}_i a_j
        for cor_inds in [cor_alpha_inds, cor_beta_inds]:
            for act_inds in [act_alpha_inds, act_beta_inds]:
                cor_1 = cor_inds.new_index()
                act_2 = act_inds.new_index()
                cor_3 = cor_inds.new_index()
                act_4 = act_inds.new_index()

                v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
                rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
                V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        # External-Active-External-Active: <ax||by> \gamma^x_y a^{\dag}_a a_b
        for vir_inds in [vir_alpha_inds, vir_beta_inds]:
            for act_inds in [act_alpha_inds, act_beta_inds]:
                vir_1 = vir_inds.new_index()
                act_2 = act_inds.new_index()
                vir_3 = vir_inds.new_index()
                act_4 = act_inds.new_index()

                v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
                rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
                V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        # <rs||pq> a^{\dag}_p a^{\dag}_q a_s a_r
        exceptions_list = [[options.active_type, options.active_type, options.active_type, options.active_type]]

        for ind_type in [options.active_type, options.virtual_type]:
            exceptions_list.extend([[options.core_type, ind_type, options.core_type, ind_type],
                                    [options.core_type, ind_type, ind_type, options.core_type],
                                    [ind_type, options.core_type, options.core_type, ind_type],
                                    [ind_type, options.core_type, ind_type, options.core_type]])

        all_cor_types = [[options.core_type, options.core_type, options.core_type, options.core_type]]
        all_vir_types = [[options.virtual_type, options.virtual_type, options.virtual_type, options.virtual_type]]

        for ind_type_1 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
            for ind_type_2 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                for ind_type_3 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    for ind_type_4 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
            for ind_type_2 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                for ind_type_3 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                    for ind_type_4 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
            for ind_type_2 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                for ind_type_3 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    for ind_type_4 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-1.00, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(1.00, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(1.00, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

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
        vir_3 = vir_alpha_inds.new_index()
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

        # Active-Core-Core-Active: <yj||ix> a_j a^{\dag}_i a^{\dag}_x a_y
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        sys.stdout.flush()

        return V

    def Vperturbation_spin_integrated_explicit_cases(indices_lists):
        "Construct spin-integrated perturbation operator V."

        cor_alpha_inds, cor_beta_inds, act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        v2e_sym_braket = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]
        h1e_sym = [symmetry((1,0), 1)]

        V = []

        # Active-Core: <x|h|i> a^{\dag}_x a_i
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        # Core-Active: <i|h|x> a^{\dag}_i a_x
        act_1 = act_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        # Active-External: <x|h|a> a^{\dag}_x a_a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        # External-Active: <a|h|x> a^{\dag}_a a_x
        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        # External-Core: <a|h|i> a^{\dag}_a a_i
        cor_1 = cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        cor_1 = cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        # Core-External: <i|h|a> a^{\dag}_i a_a
        vir_1 = vir_alpha_inds.new_index()
        cor_2 = cor_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        cor_2 = cor_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        # Core-Active-Core-Active: <ix||jy> \gamma^x_y a^{\dag}_i a_j
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cor_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cor_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cor_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cor_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        # External-Active-External-Active: <ax||by> \gamma^x_y a^{\dag}_a a_b
        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        # <rs||pq> a^{\dag}_p a^{\dag}_q a_s a_r
        for ind_type_1 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cor_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cor_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

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
        vir_3 = vir_alpha_inds.new_index()
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

        # Active-Core-Core-Active: <yj||ix> a^{\dag}_i a^{\dag}_x a_j a_y
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

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
        vir_3 = vir_alpha_inds.new_index()
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

        # Active-Core-Core-Active: <yj||ix> a_j a^{\dag}_i a^{\dag}_x a_y
        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        sys.stdout.flush()

        return V

    def Vperturbation_cvs_spin_orbital(indices_lists):
        "Construct spin-orbital perturbation operator V."

        cvs_cor_inds, cvs_val_inds, act_inds, vir_inds = indices_lists

        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        v2e_sym_braket = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]
        h1e_sym = [symmetry((1,0), 1)]

        V = []

        # Active-Core: <x|h|i> a^{\dag}_x a_i
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            act_2 = act_inds.new_index()
            v_ten = tensor('h', [cor_1, act_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        # Core-Active: <i|h|x> a^{\dag}_i a_x
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            act_1 = act_inds.new_index()
            cor_2 = cvs_inds.new_index()
            v_ten = tensor('h', [act_1, cor_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        # Active-External: <x|h|a> a^{\dag}_x a_a
        act_1 = act_inds.new_index()
        vir_2 = vir_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        # External-Active: <a|h|x> a^{\dag}_a a_x
        vir_1 = vir_inds.new_index()
        act_2 = act_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        # External-Core: <a|h|i> a^{\dag}_a a_i
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            cor_1 = cvs_inds.new_index()
            vir_2 = vir_inds.new_index()
            v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        # Core-External: <i|h|a> a^{\dag}_i a_a
        for cvs_inds in [cvs_cor_inds, cvs_val_inds]:
            vir_1 = vir_inds.new_index()
            cor_2 = cvs_inds.new_index()
            v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        # Core-Active-Core-Active: <ix||jy> \gamma^x_y a^{\dag}_i a_j
        for cvs_inds_1 in [cvs_cor_inds, cvs_val_inds]:
            for cvs_inds_2 in [cvs_cor_inds, cvs_val_inds]:
                cor_1 = cvs_inds_1.new_index()
                act_2 = act_inds.new_index()
                cor_3 = cvs_inds_2.new_index()
                act_4 = act_inds.new_index()

                v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
                rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
                V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        # External-Active-External-Active: <ax||by> \gamma^x_y a^{\dag}_a a_b
        vir_1 = vir_inds.new_index()
        act_2 = act_inds.new_index()
        vir_3 = vir_inds.new_index()
        act_4 = act_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        # <rs||pq> a^{\dag}_p a^{\dag}_q a_s a_r
        exceptions_list = [[options.active_type, options.active_type, options.active_type, options.active_type]]

        for cvs_ind_type_1 in [options.cvs_core_type, options.cvs_valence_type]:
            for cvs_ind_type_2 in [options.cvs_core_type, options.cvs_valence_type]:
                for ind_type in [options.active_type, options.virtual_type]:
                    exceptions_list.extend([[cvs_ind_type_1, ind_type, cvs_ind_type_2, ind_type],
                                            [cvs_ind_type_1, ind_type, ind_type, cvs_ind_type_2],
                                            [ind_type, cvs_ind_type_1, cvs_ind_type_2, ind_type],
                                            [ind_type, cvs_ind_type_1, ind_type, cvs_ind_type_2]])

        all_cor_types = []
        for cvs_ind_type_1 in [options.cvs_core_type, options.cvs_valence_type]:
            for cvs_ind_type_2 in [options.cvs_core_type, options.cvs_valence_type]:
                for cvs_ind_type_3 in [options.cvs_core_type, options.cvs_valence_type]:
                    for cvs_ind_type_4 in [options.cvs_core_type, options.cvs_valence_type]:
                        all_cor_types.append([cvs_ind_type_1, cvs_ind_type_2, cvs_ind_type_3, cvs_ind_type_4])

        all_vir_types = [[options.virtual_type, options.virtual_type, options.virtual_type, options.virtual_type]]

        for ind_type_1 in (cvs_cor_inds, cvs_val_inds, act_inds, vir_inds):
            for ind_type_2 in (cvs_cor_inds, cvs_val_inds, act_inds, vir_inds):
                for ind_type_3 in (cvs_cor_inds, cvs_val_inds, act_inds, vir_inds):
                    for ind_type_4 in (cvs_cor_inds, cvs_val_inds, act_inds, vir_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        # External-Core-Core-External: <bj||ia> a^{\dag}_a a_b a_j a^{\dag}_i
        for cvs_inds_1 in [cvs_cor_inds, cvs_val_inds]:
            for cvs_inds_2 in [cvs_cor_inds, cvs_val_inds]:
                cor_1 = cvs_inds_1.new_index()
                vir_2 = vir_inds.new_index()
                vir_3 = vir_inds.new_index()
                cor_4 = cvs_inds_2.new_index()
                v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        # Active-Core-Core-Active: <yj||ix> a_j a^{\dag}_i a^{\dag}_x a_y
        for cvs_inds_1 in [cvs_cor_inds, cvs_val_inds]:
            for cvs_inds_2 in [cvs_cor_inds, cvs_val_inds]:
                cor_1 = cvs_inds_1.new_index()
                act_2 = act_inds.new_index()
                act_3 = act_inds.new_index()
                cor_4 = cvs_inds_2.new_index()
                v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        sys.stdout.flush()

        return V

    def Vperturbation_cvs_spin_integrated(indices_lists):
        "Construct spin-integrated perturbation operator V."

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        v2e_sym_braket = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]
        h1e_sym = [symmetry((1,0), 1)]

        V = []

        # Active-Core: <x|h|i> a^{\dag}_x a_i
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            act_2 = act_alpha_inds.new_index()
            v_ten = tensor('h', [cor_1, act_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            act_2 = act_beta_inds.new_index()
            v_ten = tensor('h', [cor_1, act_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        # Core-Active: <i|h|x> a^{\dag}_i a_x
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            act_1 = act_alpha_inds.new_index()
            cor_2 = cvs_alpha_inds.new_index()
            v_ten = tensor('h', [act_1, cor_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            act_1 = act_beta_inds.new_index()
            cor_2 = cvs_beta_inds.new_index()
            v_ten = tensor('h', [act_1, cor_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        # Active-External: <x|h|a> a^{\dag}_x a_a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        # External-Active: <a|h|x> a^{\dag}_a a_x
        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        # External-Core: <a|h|i> a^{\dag}_a a_i
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            cor_1 = cvs_alpha_inds.new_index()
            vir_2 = vir_alpha_inds.new_index()
            v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            cor_1 = cvs_beta_inds.new_index()
            vir_2 = vir_beta_inds.new_index()
            v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        # Core-External: <i|h|a> a^{\dag}_i a_a
        for cvs_alpha_inds in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            vir_1 = vir_alpha_inds.new_index()
            cor_2 = cvs_alpha_inds.new_index()
            v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        for cvs_beta_inds in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            vir_1 = vir_beta_inds.new_index()
            cor_2 = cvs_beta_inds.new_index()
            v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
            V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        # Core-Active-Core-Active: <ix||jy> \gamma^x_y a^{\dag}_i a_j
        for cvs_inds_1 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            for cvs_inds_2 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
                for act_inds in [act_alpha_inds, act_beta_inds]:
                    cor_1 = cvs_inds_1.new_index()
                    act_2 = act_inds.new_index()
                    cor_3 = cvs_inds_2.new_index()
                    act_4 = act_inds.new_index()

                    v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
                    rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
                    V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        for cvs_inds_1 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            for cvs_inds_2 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
                for act_inds in [act_alpha_inds, act_beta_inds]:
                    cor_1 = cvs_inds_1.new_index()
                    act_2 = act_inds.new_index()
                    cor_3 = cvs_inds_2.new_index()
                    act_4 = act_inds.new_index()

                    v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
                    rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
                    V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        # External-Active-External-Active: <ax||by> \gamma^x_y a^{\dag}_a a_b
        for vir_inds in [vir_alpha_inds, vir_beta_inds]:
            for act_inds in [act_alpha_inds, act_beta_inds]:
                vir_1 = vir_inds.new_index()
                act_2 = act_inds.new_index()
                vir_3 = vir_inds.new_index()
                act_4 = act_inds.new_index()

                v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
                rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
                V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        # <rs||pq> a^{\dag}_p a^{\dag}_q a_s a_r
        exceptions_list = [[options.active_type, options.active_type, options.active_type, options.active_type]]

        for cvs_ind_type_1 in [options.cvs_core_type, options.cvs_valence_type]:
            for cvs_ind_type_2 in [options.cvs_core_type, options.cvs_valence_type]:
                for ind_type in [options.active_type, options.virtual_type]:
                    exceptions_list.extend([[cvs_ind_type_1, ind_type, cvs_ind_type_2, ind_type],
                                            [cvs_ind_type_1, ind_type, ind_type, cvs_ind_type_2],
                                            [ind_type, cvs_ind_type_1, cvs_ind_type_2, ind_type],
                                            [ind_type, cvs_ind_type_1, ind_type, cvs_ind_type_2]])

        all_cor_types = []
        for cvs_ind_type_1 in [options.cvs_core_type, options.cvs_valence_type]:
            for cvs_ind_type_2 in [options.cvs_core_type, options.cvs_valence_type]:
                for cvs_ind_type_3 in [options.cvs_core_type, options.cvs_valence_type]:
                    for cvs_ind_type_4 in [options.cvs_core_type, options.cvs_valence_type]:
                        all_cor_types.append([cvs_ind_type_1, cvs_ind_type_2, cvs_ind_type_3, cvs_ind_type_4])

        all_vir_types = [[options.virtual_type, options.virtual_type, options.virtual_type, options.virtual_type]]

        for ind_type_1 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
            for ind_type_2 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                for ind_type_3 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    for ind_type_4 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
            for ind_type_2 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                for ind_type_3 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                    for ind_type_4 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
            for ind_type_2 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                for ind_type_3 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    for ind_type_4 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_1 = ind_type_1.new_index()
                        ind_2 = ind_type_2.new_index()
                        ind_3 = ind_type_3.new_index()
                        ind_4 = ind_type_4.new_index()

                        tens_spatial_types = [get_spatial_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)]

                        if tens_spatial_types in all_cor_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-1.00, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_spatial_types in all_vir_types:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(1.00, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif tens_spatial_types not in exceptions_list:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(1.00, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        # External-Core-Core-External: <bj||ia> a^{\dag}_a a_b a_j a^{\dag}_i
        for cvs_alpha_inds_1 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            for cvs_alpha_inds_2 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
                cor_1 = cvs_alpha_inds_1.new_index()
                vir_2 = vir_alpha_inds.new_index()
                vir_3 = vir_alpha_inds.new_index()
                cor_4 = cvs_alpha_inds_2.new_index()
                v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

                cor_1 = cvs_alpha_inds_1.new_index()
                vir_2 = vir_beta_inds.new_index()
                vir_3 = vir_beta_inds.new_index()
                cor_4 = cvs_alpha_inds_2.new_index()
                v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        for cvs_beta_inds_1 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            for cvs_beta_inds_2 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
                cor_1 = cvs_beta_inds_1.new_index()
                vir_2 = vir_beta_inds.new_index()
                vir_3 = vir_beta_inds.new_index()
                cor_4 = cvs_beta_inds_2.new_index()
                v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

                cor_1 = cvs_beta_inds_1.new_index()
                vir_2 = vir_alpha_inds.new_index()
                vir_3 = vir_alpha_inds.new_index()
                cor_4 = cvs_beta_inds_2.new_index()
                v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        for cvs_alpha_inds_1 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            for cvs_beta_inds_2 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
                cor_1 = cvs_alpha_inds_1.new_index()
                vir_2 = vir_beta_inds.new_index()
                vir_3 = vir_alpha_inds.new_index()
                cor_4 = cvs_beta_inds_2.new_index()
                v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        for cvs_beta_inds_1 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            for cvs_alpha_inds_2 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
                cor_1 = cvs_beta_inds_1.new_index()
                vir_2 = vir_alpha_inds.new_index()
                vir_3 = vir_beta_inds.new_index()
                cor_4 = cvs_alpha_inds_2.new_index()
                v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        # Active-Core-Core-Active: <yj||ix> a_j a^{\dag}_i a^{\dag}_x a_y
        for cvs_alpha_inds_1 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            for cvs_alpha_inds_2 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
                cor_1 = cvs_alpha_inds_1.new_index()
                act_2 = act_alpha_inds.new_index()
                act_3 = act_alpha_inds.new_index()
                cor_4 = cvs_alpha_inds_2.new_index()
                v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

                cor_1 = cvs_alpha_inds_1.new_index()
                act_2 = act_beta_inds.new_index()
                act_3 = act_beta_inds.new_index()
                cor_4 = cvs_alpha_inds_2.new_index()
                v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        for cvs_beta_inds_1 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            for cvs_beta_inds_2 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
                cor_1 = cvs_beta_inds_1.new_index()
                act_2 = act_beta_inds.new_index()
                act_3 = act_beta_inds.new_index()
                cor_4 = cvs_beta_inds_2.new_index()
                v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

                cor_1 = cvs_beta_inds_1.new_index()
                act_2 = act_alpha_inds.new_index()
                act_3 = act_alpha_inds.new_index()
                cor_4 = cvs_beta_inds_2.new_index()
                v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        for cvs_alpha_inds_1 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
            for cvs_beta_inds_2 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
                cor_1 = cvs_alpha_inds_1.new_index()
                act_2 = act_beta_inds.new_index()
                act_3 = act_alpha_inds.new_index()
                cor_4 = cvs_beta_inds_2.new_index()
                v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        for cvs_beta_inds_1 in [cvs_cor_beta_inds, cvs_val_beta_inds]:
            for cvs_alpha_inds_2 in [cvs_cor_alpha_inds, cvs_val_alpha_inds]:
                cor_1 = cvs_beta_inds_1.new_index()
                act_2 = act_alpha_inds.new_index()
                act_3 = act_beta_inds.new_index()
                cor_4 = cvs_alpha_inds_2.new_index()
                v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
                V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        sys.stdout.flush()

        return V

    def Vperturbation_cvs_spin_integrated_explicit_cases(indices_lists):
        "Construct spin-integrated perturbation operator V."

        cvs_cor_alpha_inds, cvs_cor_beta_inds, cvs_val_alpha_inds, cvs_val_beta_inds, \
        act_alpha_inds, act_beta_inds, vir_alpha_inds, vir_beta_inds = indices_lists

        v2e_sym = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1)]
        v2e_sym_braket = [symmetry((1,0,2,3), -1), symmetry((0,1,3,2), -1), symmetry((2,3,0,1), 1)]
        h1e_sym = [symmetry((1,0), 1)]

        V = []

        # Active-Core: <x|h|i> a^{\dag}_x a_i
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(cor_1)]))

        # Core-Active: <i|h|x> a^{\dag}_i a_x
        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        act_1 = act_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        v_ten = tensor('h', [act_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(act_1)]))

        # Active-External: <x|h|a> a^{\dag}_x a_a
        act_1 = act_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        act_1 = act_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [act_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(act_1)]))

        # External-Active: <a|h|x> a^{\dag}_a a_x
        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, act_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(act_2), desOp(vir_1)]))

        # External-Core: <a|h|i> a^{\dag}_a a_i
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        v_ten = tensor('h', [cor_1, vir_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(cor_1)]))

        # Core-External: <i|h|a> a^{\dag}_i a_a
        vir_1 = vir_alpha_inds.new_index()
        cor_2 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        cor_2 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        vir_1 = vir_alpha_inds.new_index()
        cor_2 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        cor_2 = cvs_val_beta_inds.new_index()
        v_ten = tensor('h', [vir_1, cor_2], h1e_sym)
        V.append(term(1.0, [], [v_ten, creOp(cor_2), desOp(vir_1)]))

        # Core-Active-Core-Active: <ix||jy> \gamma^x_y a^{\dag}_i a_j
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_cor_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_cor_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_cor_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_cor_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_val_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_val_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_val_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_val_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_val_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_val_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_val_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_val_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_cor_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_cor_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        cor_3 = cvs_cor_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        cor_3 = cvs_cor_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [cor_1, act_2, cor_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(1.0, [], [v_ten, rdm_ten, desOp(cor_1), creOp(cor_3)]))

        # External-Active-External-Active: <ax||by> \gamma^x_y a^{\dag}_a a_b
        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        vir_1 = vir_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        act_4 = act_alpha_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        vir_1 = vir_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        act_4 = act_beta_inds.new_index()

        v_ten = tensor('v', [vir_1, act_2, vir_3, act_4], v2e_sym)
        rdm_ten = creDesTensor([creOp(act_2), desOp(act_4)])
        V.append(term(-1.0, [], [v_ten, rdm_ten, creOp(vir_3), desOp(vir_1)]))

        # <rs||pq> a^{\dag}_p a^{\dag}_q a_s a_r
        for ind_type_1 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        for ind_type_1 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
            ind_1 = ind_type_1.new_index()

            for ind_type_2 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                ind_2 = ind_type_2.new_index()

                for ind_type_3 in (cvs_cor_alpha_inds, cvs_val_alpha_inds, act_alpha_inds, vir_alpha_inds):
                    ind_3 = ind_type_3.new_index()

                    for ind_type_4 in (cvs_cor_beta_inds, cvs_val_beta_inds, act_beta_inds, vir_beta_inds):
                        ind_4 = ind_type_4.new_index()

                        tens_cor_type = all([is_core_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_act_type = all([is_active_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])
                        tens_vir_type = all([is_virtual_index_type(ind.indType) for ind in (ind_1, ind_2, ind_4, ind_3)])

                        if tens_cor_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(-0.25, [], [v_ten, desOp(ind_3), desOp(ind_4), creOp(ind_1), creOp(ind_2)]))

                        elif tens_vir_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym_braket)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

                        elif not tens_act_type:
                            v_ten = tensor('v', [ind_3, ind_4, ind_1, ind_2], v2e_sym)
                            V.append(term(0.25, [], [v_ten, creOp(ind_1), creOp(ind_2), desOp(ind_4), desOp(ind_3)]))

        # External-Core-Core-External: <bj||ia> a^{\dag}_i a^{\dag}_a a_j a_b
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(vir_2), desOp(cor_4), desOp(vir_3)]))

        # Active-Core-Core-Active: <yj||ix> a^{\dag}_i a^{\dag}_x a_j a_y
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(-1.0, [], [v_ten, creOp(cor_1), creOp(act_2), desOp(cor_4), desOp(act_3)]))

        # External-Core-Core-External: <bj||ia> a^{\dag}_a a_b a_j a^{\dag}_i
        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        vir_2 = vir_beta_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        cor_1 = cvs_val_beta_inds.new_index()
        vir_2 = vir_alpha_inds.new_index()
        vir_3 = vir_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [vir_3, cor_4, cor_1, vir_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, creOp(vir_2), desOp(vir_3), desOp(cor_4), creOp(cor_1)]))

        # Active-Core-Core-Active: <yj||ix> a_j a^{\dag}_i a^{\dag}_x a_y
        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_val_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_cor_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_val_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_alpha_inds.new_index()
        act_2 = act_beta_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_beta_inds.new_index()
        cor_4 = cvs_cor_alpha_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        cor_1 = cvs_val_beta_inds.new_index()
        act_2 = act_alpha_inds.new_index()
        act_3 = act_alpha_inds.new_index()
        cor_4 = cvs_cor_beta_inds.new_index()
        v_ten = tensor('v', [act_3, cor_4, cor_1, act_2], v2e_sym)
        V.append(term(1.0, [], [v_ten, desOp(cor_4), creOp(cor_1), creOp(act_2), desOp(act_3)]))

        sys.stdout.flush()

        return V

    if spin_integrated:
        if explicit_spin_cases:
            if cvs_indices:
                V = Vperturbation_cvs_spin_integrated_explicit_cases(indices_lists)
            else:
                V = Vperturbation_spin_integrated_explicit_cases(indices_lists)
        else:
            if cvs_indices:
                V = Vperturbation_cvs_spin_integrated(indices_lists)
            else:
                V = Vperturbation_spin_integrated(indices_lists)
    else:
        if cvs_indices:
            V = Vperturbation_cvs_spin_orbital(indices_lists)
        else:
            V = Vperturbation_spin_orbital(indices_lists)

    return V

def getT(order = 1, spin_integrated = False, explicit_spin_cases = True, cvs_indices = False):
    "Get T amplitudes (order)."

    indices_lists = create_dummy_indices_list(spin_integrated, cvs_indices)

    T = Tamplitude(order, indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

    return T

def getV(spin_integrated = False, explicit_spin_cases = True, cvs_indices = False):
    "Get V pertubation terms."

    indices_lists = create_dummy_indices_list(spin_integrated, cvs_indices)

    V = Vperturbation(indices_lists, spin_integrated, explicit_spin_cases, cvs_indices)

    return V
