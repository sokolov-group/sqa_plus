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
# Author: Carlos E. V. de Moura <carlosevmoura@gmail.com>
#

import sys, time
import itertools

from sqaMatrixBlock import dummyLabel, reorder_tensor_indices
from sqaIndex import get_spatial_index_type, get_spin_index_type, is_alpha_index_type, is_beta_index_type
from sqaTerm import combineTerms, termChop
from sqaTensor import creOp, desOp, creDesTensor, kroneckerDelta
from sqaOptions import options
from sqaSymmetry import symmetry

def convertSpinIntegratedToAdapted(terms_si):
    "Convert Spin-Integrated Terms to Spin-Adapted Quantities."

    startTime = time.time()
    options.print_header("Converting Spin-Integrated Tensors to Spin-Adapted")
    sys.stdout.flush()

    # Convert Cre/Des Objects to RDM Objects
    options.print_divider()
    convert_credes_to_rdm(terms_si, trans_rdm = False)

    reorder_rdm_spin_indices(terms_si)
    dummyLabel(terms_si)
    len_terms_si = len(terms_si)
    
    # Convert Kronecker Delta to Spin-Adapted Formulation
    terms_sa = convert_kdelta_si_to_sa(terms_si)

    # Convert eigenvalues to Spin-Adapted Formulation
    terms_sa = convert_e_si_to_sa(terms_sa)

    # Convert 1e- integrals to Spin-Adapted Formulation
    terms_sa = convert_h1e_si_to_sa(terms_sa)

    # Convert 2e- integrals to Spin-Adapted Formulation
    terms_sa = convert_v2e_si_to_sa(terms_sa)

    # Convert T amplitudes to Spin-Adapted Formulation
    terms_sa = convert_t_amplitudes_si_to_sa(terms_sa)

    # Convert RDMs to Spin-Adapted Formulation
    terms_sa = convert_rdms_si_to_sa(terms_sa)

    # Combine Spin-Adapted Terms
    options.print_divider()
    for term_sa in terms_sa:
        term_sa.isInCanonicalForm = False

    num_terms_sa = len(terms_sa)
    print("\nCombining {:} spin-adapted terms...\n".format(num_terms_sa))
    combineTerms(terms_sa)

    # Reorder tensors to Chemist's Notation
    reorder_v2e_indices_notation(terms_sa)
    reorder_rdm_indices_notation(terms_sa)
    options.chemists_notation = True

    reorder_tensor_indices(terms_sa)
    dummyLabel(terms_sa)
    print("\n{:} spin-adapted terms combined.".format(num_terms_sa - len(terms_sa)))
    options.print_divider()

    # Set terms as spin-adapted in sqaOptions class
    options.spin_adapted = True

    # Print Spin-Adapted Equations
    options.print_header("Spin-adapted equations")
    terms_sa.sort()
    for term in terms_sa:
        print("{:}".format(term))

    print("\nTotal spin-integrated terms: {:}".format(len_terms_si))
    print("Total spin-adapted terms: {:}".format(len(terms_sa)))
    print("Spin-adaptation automation time :  {:.3f} seconds".format(time.time() - startTime))
    options.print_divider()
    sys.stdout.flush()

    return terms_sa

def convert_credes_to_rdm(_terms_credes, trans_rdm = False):
    "Convert Cre/Des Objects to RDM Objects"

    print("Convert Cre/Des objects to RDM objects...")
    sys.stdout.flush()

    for term_credes_ind, term_credes in enumerate(_terms_credes):

        ## List for storing cre/des operators
        credes_ops = []

        ## Append all cre/des operators to list
        for tens_credes in term_credes.tensors:
            if isinstance(tens_credes, creOp) or isinstance(tens_credes, desOp):
                credes_ops.append(tens_credes)

        ## Modify term in list to use creDesTensor object instead of cre/des objects
        if credes_ops:
            _terms_credes[term_credes_ind].tensors = [tens for tens in term_credes.tensors if tens not in credes_ops]
            ten_rdm = creDesTensor(credes_ops, trans_rdm)
            _terms_credes[term_credes_ind].tensors.append(ten_rdm)

    print("Done!")
    options.print_divider()
    return

def reorder_v2e_indices_notation(_terms_v2e):
    print("Reorder 2e- integrals in Chemists' notation...")
    sys.stdout.flush()

    if options.verbose:
        print("\nv2e notation in spin-adapted chemists' notation:")
        print("-> v2e[p,q,r,s] => v2e[p,r,q,s]")

    for term_v2e_ind, term_v2e in enumerate(_terms_v2e):

        ## List for storing v2e objects
        v2e_tensors = []
        v2e_tensors_ind = []

        ## Define indices symmetries
        v2e_symm = [symmetry((1,0,2,3), 1), symmetry((2,3,0,1), 1)]

        ## Append all v2e objects to list
        for ten_v2e_ind, ten_v2e in enumerate(term_v2e.tensors):
            if (ten_v2e.name == 'v') and (len(ten_v2e.indices) == 4):
                v2e_tensors.append(ten_v2e)
                v2e_tensors_ind.append(ten_v2e_ind)

        ## Convert from Physicists's to Chemists's Notation
        if v2e_tensors:
            for ten_v2e_ind, ten_v2e in zip(v2e_tensors_ind, v2e_tensors):

                # v2e[p,q,r,s] => v2e[p,r,q,s]
                ten_v2e.indices = [ten_v2e.indices[i] for i in [0, 2, 1, 3]]
                ten_v2e.symmetries = v2e_symm

                _terms_v2e[term_v2e_ind].tensors[ten_v2e_ind] = ten_v2e.copy()

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return

def reorder_rdm_indices_notation(_terms_rdm):
    print("Reorder RDM tensor indices in Chemists' notation...")
    sys.stdout.flush()

    if options.verbose:
        print("\nRDMs notation in spin-adapted chemists' notation:")
        print("-> rdm2[p,q,r,s] = < p^+ q^+ s r >")
        print("-> rdm3[p,q,r,s,t,u] = < p^+ q^+ r^+ u t s >")
        print("-> rdm4[p,q,r,s,t,u,v,w] = < p^+ q^+ r^+ s^+ w v u t >\n")

    for term_rdm_ind, term_rdm in enumerate(_terms_rdm):

        ## List for storing RDM objects
        rdm_tensors = []
        rdm_tensors_ind = []

        ## Define indices symmetries
        rdm2_symm = [symmetry((1,0,3,2), 1), symmetry((2,3,0,1), 1)]

        rdm3_symm = [symmetry((1,0,2,4,3,5), 1), symmetry((0,2,1,3,5,4), 1), symmetry((3,4,5,0,1,2), 1)]

        rdm4_symm = [symmetry((1,0,2,3,5,4,6,7), 1), symmetry((0,2,1,3,4,6,5,7), 1),
                     symmetry((0,1,3,2,4,5,7,6), 1), symmetry((4,5,6,7,0,1,2,3), 1)]

        ## Append all RDM objects to list
        for ten_rdm_ind, ten_rdm in enumerate(term_rdm.tensors):
            if isinstance(ten_rdm, creDesTensor):
                rdm_tensors.append(ten_rdm)
                rdm_tensors_ind.append(ten_rdm_ind)

        ## Convert from Physicists's to Chemists's Notation
        if rdm_tensors:
            for ten_rdm_ind, ten_rdm in zip(rdm_tensors_ind, rdm_tensors):

                # rdm2[p,q,r,s] = < p^+ q^+ s r >
                if len(ten_rdm.indices) == 4:
                    ten_rdm.indices = [ten_rdm.indices[i] for i in [0, 1, 3, 2]]
                    ten_rdm.symmetries = rdm2_symm

                # rdm3[p,q,r,s,t,u] = < p^+ q^+ r^+ u t s >
                if len(ten_rdm.indices) == 6:
                    ten_rdm.indices = [ten_rdm.indices[i] for i in [0, 1, 2, 5, 4, 3]]
                    ten_rdm.symmetries = rdm3_symm

                # rdm4[p,q,r,s,t,u,v,w] = < p^+ q^+ r^+ s^+ w v u t >
                if len(ten_rdm.indices) == 8:
                    ten_rdm.indices = [ten_rdm.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    ten_rdm.symmetries = rdm4_symm

                _terms_rdm[term_rdm_ind].tensors[ten_rdm_ind] = ten_rdm.copy()

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return

def reorder_rdm_spin_indices(_terms_rdm):
    "Reorder RDM tensor indices according to: alpha < beta"

    print("Reorder RDM tensor indices according to: alpha < beta")
    sys.stdout.flush()

    for term_rdm_ind, term_rdm in enumerate(_terms_rdm):

        ## List for storing RDM objects
        rdm_tensors = []
        rdm_tensors_ind = []

        ## Append all RDM objects to list
        for ten_rdm_ind, ten_rdm in enumerate(term_rdm.tensors):
            #TODO: Consider if this approach can be applied to other RDMs
            if isinstance(ten_rdm, creDesTensor) and len(ten_rdm.indices) == 8:
                rdm_tensors.append(ten_rdm)
                rdm_tensors_ind.append(ten_rdm_ind)

        if rdm_tensors:
            for ten_rdm_ind, ten_rdm in zip(rdm_tensors_ind, rdm_tensors):

                # Reordering indices according to alpha > beta
                original_rank = []
                original_name_rank = []
                for ind in ten_rdm.indices:
                    original_name_rank.append(ind.name)
                    if is_alpha_index_type(ind):
                        original_rank.append(1)
                    elif is_beta_index_type(ind):
                        original_rank.append(0)

                if sum(original_rank) > len(ten_rdm.indices)//2:
                    high_spin = True
                else:
                    high_spin = False

                full_permutes_indices, full_permutes_factors = ten_rdm.symPermutes()

                # Remove bra-ket permutes
                permutes_indices = []
                permutes_factors = []
                for permute_indices, permute_factors in zip(full_permutes_indices, full_permutes_factors):
                    if permute_indices[0] < len(ten_rdm.indices) // 2:
                        permutes_indices.append(permute_indices)
                        permutes_factors.append(permute_factors)

                permutes_rank = []
                permutes_name_rank = []

                for permute_indices in permutes_indices:
                    permute_rank = []
                    permute_name_rank = []

                    permuted_indices = [ten_rdm.indices[ind] for ind in permute_indices]
                    for ind in permuted_indices:
                        permute_name_rank.append(ind.name)
                        if is_alpha_index_type(ind):
                            permute_rank.append(1)
                        elif is_beta_index_type(ind):
                            permute_rank.append(0)

                    permutes_name_rank.append(permute_name_rank)
                    permutes_rank.append(permute_rank)

                permutes_rank, permutes_name_rank, permutes_indices, permutes_factors = zip(*sorted(zip(permutes_rank, permutes_name_rank, permutes_indices, permutes_factors), reverse=True))

                if high_spin:
                    permutes_rank_ind = [ind for ind, rank in enumerate(permutes_rank) if rank == permutes_rank[-1]][-1]
                else:
                    permutes_rank_ind = [ind for ind, rank in enumerate(permutes_rank) if rank == permutes_rank[0]][-1]

                permuted_indices = [ten_rdm.indices[ind] for ind in permutes_indices[permutes_rank_ind]]
                order_factor = permutes_factors[permutes_rank_ind]

                ordered_tensor = ten_rdm.copy()
                ordered_tensor.indices = permuted_indices

                if options.verbose:
                    print("{:}    --->    {:} (factor = {:})".format(ten_rdm, ordered_tensor, order_factor))

                _terms_rdm[term_rdm_ind].tensors[ten_rdm_ind] = ordered_tensor.copy()
                _terms_rdm[term_rdm_ind].scale(order_factor)

    print("Done!")
    options.print_divider()
    sys.stdout.flush()
    return

def remove_spin_index_type(_tensor):
    "Remove spin index type of a given spin-integrated tensor"

    for index in _tensor.indices:
        index_spatial_type = get_spatial_index_type(index)
        index.indType = (index_spatial_type,)

    return _tensor

def convert_h1e_si_to_sa(_terms_h1e_si):
    "Convert h1e Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting 1e- integrals to spin-adapted formulation...")

    # Define 1e- indices lists
    inds_aa = [options.alpha_type, options.alpha_type]
    inds_bb = [options.beta_type,  options.beta_type]

    # Convert 1e- objects in each term
    terms_h1e_sa = []
    for term_h1e_si in _terms_h1e_si:
        ## Select the 1e- object
        ten_h1e = False
        for ten_ind, ten in enumerate(term_h1e_si.tensors):
            if ten.name == 'h' and len(ten.indices) == 2:
                ten_h1e = ten
                ten_h1e_ind = ten_ind

        if ten_h1e:
            ten_h1e_spin_inds = [get_spin_index_type(ind) for ind in ten_h1e.indices]
            if ten_h1e_spin_inds in [inds_aa, inds_bb]:
                term_h1e_sa = term_h1e_si.copy()
                term_h1e_sa.tensors[ten_h1e_ind] = remove_spin_index_type(ten_h1e)
                terms_h1e_sa.append(term_h1e_sa)

            else:
                term_h1e_sa = term_h1e_si.copy()
                term_h1e_sa.scale(0.0)
                terms_h1e_sa.append(term_h1e_sa)
        else:
            terms_h1e_sa.append(term_h1e_si)

    termChop(terms_h1e_sa)

    # Print 1e- spin-adapted terms
    if options.verbose:
        print("")
        for term_h1e_sa in terms_h1e_sa:
            print(term_h1e_sa)
        print("")

    print("Done!")
    return terms_h1e_sa

def convert_v2e_si_to_sa(_terms_v2e_si):
    "Convert v2e Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting 2e- integrals to spin-adapted formulation...")

    # Define Spin-Adapted 2e- integrals Symmetries
    v2e_sa_symm = [symmetry((1,0,3,2), 1), symmetry((2,3,0,1), 1)]

    # Define 2e- indices lists
    inds_aaaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type]
    inds_bbbb = [options.beta_type,  options.beta_type,  options.beta_type,  options.beta_type]

    inds_abab = [options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type]
    inds_baba = [options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type]
    inds_abba = [options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type]
    inds_baab = [options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type]

    # Convert 2e- objects in each term
    terms_v2e_sa = []
    for term_v2e_si in _terms_v2e_si:
        # List for storing 2e- objects
        tens_v2e = []
        tens_v2e_ind = []

        for ten_ind, ten in enumerate(term_v2e_si.tensors):
            if ten.name == 'v' and len(ten.indices) == 4:
                tens_v2e.append(ten)
                tens_v2e_ind.append(ten_ind)

        ## Convert 2e- object
        if tens_v2e:
            tens_v2e_sa = []
            consts_v2e_sa = []

            for ten_v2e in tens_v2e:

                ten_v2e_inds = [get_spatial_index_type(ind) for ind in ten_v2e.indices]
                ten_v2e_spin_inds = [get_spin_index_type(ind) for ind in ten_v2e.indices]

                if ((ten_v2e_inds[0] == ten_v2e_inds[1] == ten_v2e_inds[2] == ten_v2e_inds[3]) or

                   (((ten_v2e_inds[0] == ten_v2e_inds[1]) and (ten_v2e_inds[2] == ten_v2e_inds[3]) and
                     (ten_v2e_inds[0] != ten_v2e_inds[2]) and (ten_v2e_inds[1] != ten_v2e_inds[3])) or

                    ((ten_v2e_inds[0] != ten_v2e_inds[1]) and (ten_v2e_inds[2] == ten_v2e_inds[3])) or
                    ((ten_v2e_inds[0] != ten_v2e_inds[1]) and (ten_v2e_inds[2] != ten_v2e_inds[3])))):

                    if ten_v2e_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten1_v2e = ten_v2e.copy()
                        const1_v2e = 1.0

                        ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                        ten2_v2e = ten_v2e.copy()
                        ten2_v2e.indices = [ten2_v2e.indices[i] for i in [0, 1, 3, 2]]
                        const2_v2e = - 1.0

                        tens_v2e_sa.append([ten1_v2e, ten2_v2e])
                        consts_v2e_sa.append([const1_v2e, const2_v2e])

                    elif ten_v2e_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten1_v2e = ten_v2e.copy()
                        const1_v2e = 1.0

                        tens_v2e_sa.append([ten1_v2e])
                        consts_v2e_sa.append([const1_v2e])

                    elif ten_v2e_spin_inds in [inds_abba, inds_baab]:
                        ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                        ten1_v2e = ten_v2e.copy()
                        ten1_v2e.indices = [ten1_v2e.indices[i] for i in [0, 1, 3, 2]]
                        const1_v2e = - 1.0

                        tens_v2e_sa.append([ten1_v2e])
                        consts_v2e_sa.append([const1_v2e])

                    else:
                        ten1_v2e = ten_v2e.copy()
                        const1_v2e = 0.0

                        tens_v2e_sa.append([ten1_v2e])
                        consts_v2e_sa.append([const1_v2e])

                elif ((ten_v2e_inds[0] == ten_v2e_inds[1]) and (ten_v2e_inds[2] != ten_v2e_inds[3])):
                    if ten_v2e_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten1_v2e = ten_v2e.copy()
                        const1_v2e = 1.0

                        ## Spin-Adapted 2e- term: v2e(q,p,r,s)
                        ten2_v2e = ten_v2e.copy()
                        ten2_v2e.indices = [ten2_v2e.indices[i] for i in [1, 0, 2, 3]]
                        const2_v2e = - 1.0

                        tens_v2e_sa.append([ten1_v2e, ten2_v2e])
                        consts_v2e_sa.append([const1_v2e, const2_v2e])

                    elif ten_v2e_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten1_v2e = ten_v2e.copy()
                        const1_v2e = 1.0

                        tens_v2e_sa.append([ten1_v2e])
                        consts_v2e_sa.append([const1_v2e])

                    elif ten_v2e_spin_inds in [inds_abba, inds_baab]:
                        ## Spin-Adapted 2e- term: v2e(q,p,r,s)
                        ten1_v2e = ten_v2e.copy()
                        ten1_v2e.indices = [ten1_v2e.indices[i] for i in [1, 0, 2, 3]]
                        const1_v2e = - 1.0

                        tens_v2e_sa.append([ten1_v2e])
                        consts_v2e_sa.append([const1_v2e])

                    else:
                        ten1_v2e = ten_v2e.copy()
                        const1_v2e = 0.0

                        tens_v2e_sa.append([ten1_v2e])
                        consts_v2e_sa.append([const1_v2e])

            tens_v2e_sa_permut = []
            for item in list(itertools.product(*tens_v2e_sa)):
                tens_v2e_sa_permut.append(list(item))

            consts_v2e_sa_permut = []
            for item in list(itertools.product(*consts_v2e_sa)):
                consts_v2e_sa_permut.append(list(item))

            consts_v2e_sa_prod = []
            for iter in consts_v2e_sa_permut:
                prod = 1.0
                for const in iter:
                    prod = prod * const
                consts_v2e_sa_prod.append(prod)

            for tens_v2e_sa_ind, tens_v2e_sa in enumerate(tens_v2e_sa_permut):
                term_v2e_sa = term_v2e_si.copy()
                term_v2e_sa.scale(consts_v2e_sa_prod[tens_v2e_sa_ind])

                for ten_v2e_sa_ind, ten_v2e_sa in zip(tens_v2e_ind, tens_v2e_sa):
                    term_v2e_sa.tensors[ten_v2e_sa_ind] = remove_spin_index_type(ten_v2e_sa)
                    term_v2e_sa.tensors[ten_v2e_sa_ind].symmetries = v2e_sa_symm
                terms_v2e_sa.append(term_v2e_sa)

        else:
            terms_v2e_sa.append(term_v2e_si)

    termChop(terms_v2e_sa)

    # Print 2e- spin-adapted terms
    if options.verbose:
        print("")
        for term_v2e_sa in terms_v2e_sa:
            print(term_v2e_sa)
        print("")

    print("Done!")
    return terms_v2e_sa

def convert_rdms_si_to_sa(_terms_rdm_si):
    "Convert RDM Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting RDMs to spin-adapted formulation...\n")

    # Convert One-Body RDMs
    print("Converting 1-RDMs to spin-adapted formulation...")

    # Define Spin-Adapted One-Body RDMs Symmetries
    rdm1_sa_symm = [symmetry((1,0), 1)]

    # Define 1e- indices lists
    inds_aa = [options.alpha_type, options.alpha_type]
    inds_bb = [options.beta_type,  options.beta_type]

    terms_rdm1_sa = []
    for term_rdm1_si in _terms_rdm_si:
        term_rdm1_sa = term_rdm1_si.copy()
        for ten_ind, ten in enumerate(term_rdm1_si.tensors):
            if isinstance(ten, creDesTensor) and len(ten.indices) == 2:
                ten_rdm1_spin_inds = [get_spin_index_type(ind) for ind in ten.indices]

                if ten_rdm1_spin_inds in [inds_aa, inds_bb]:
                    term_rdm1_sa.scale(1.0 / 2.0)
                    term_rdm1_sa.tensors[ten_ind] = remove_spin_index_type(ten)
                    term_rdm1_sa.tensors[ten_ind].symmetries = rdm1_sa_symm
                else:
                    term_rdm1_sa.scale(0.0)
                    term_rdm1_sa.tensors[ten_ind] = remove_spin_index_type(ten)

        terms_rdm1_sa.append(term_rdm1_sa)

    termChop(terms_rdm1_sa)

    # Print 1e- spin-adapted terms
    if options.verbose:
        print("")
        for term_rdm1_sa in terms_rdm1_sa:
            print(term_rdm1_sa)
        print("")

    # Convert Two-Body RDMs
    print("Converting 2-RDMs to spin-adapted formulation...")

    # Define Spin-Adapted Two-Body RDMs Symmetries
    rdm2_sa_symm = [symmetry((1,0,3,2), 1), symmetry((3,2,1,0), 1)]

    # Define 2e- indices lists
    inds_aaaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type]
    inds_bbbb = [options.beta_type,  options.beta_type,  options.beta_type,  options.beta_type]

    inds_abab = [options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type]
    inds_baba = [options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type]
    inds_abba = [options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type]
    inds_baab = [options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type]

    terms_rdm2_sa = []
    for term_rdm2_si in terms_rdm1_sa:
        # List for storing 2-RDMs
        tens_rdm2 = []
        tens_rdm2_ind = []

        # Append all 2-RDMs to list
        for ten_ind, ten in enumerate(term_rdm2_si.tensors):
            if isinstance(ten, creDesTensor) and len(ten.indices) == 4:
                tens_rdm2.append(ten)
                tens_rdm2_ind.append(ten_ind)

        if tens_rdm2:
            tens_rdm2_sa = []
            consts_rdm2_sa = []

            for ten_rdm2 in tens_rdm2:
                ten_rdm2_spin_inds = [get_spin_index_type(ind) for ind in ten_rdm2.indices]

                if ten_rdm2_spin_inds in [inds_aaaa, inds_bbbb]:
                    ## Spin-Adapted RDM term: rdm(u,v,y,x)
                    ten1_rdm2 = ten_rdm2.copy()
                    const1_rdm2 = 1.0 / 6.0

                    ## Spin-Adapted RDM term: rdm(u,v,x,y)
                    ten2_rdm2 = ten_rdm2.copy()
                    ten2_rdm2.indices = [ten2_rdm2.indices[i] for i in [0, 1, 3, 2]]
                    const2_rdm2 = - 1.0 / 6.0

                    tens_rdm2_sa.append([ten1_rdm2, ten2_rdm2])
                    consts_rdm2_sa.append([const1_rdm2, const2_rdm2])

                elif ten_rdm2_spin_inds in [inds_abab, inds_baba]:
                    ## Spin-Adapted RDM term: rdm(u,v,y,x)
                    ten1_rdm2 = ten_rdm2.copy()
                    const1_rdm2 = - 1.0 / 6.0

                    ## Spin-Adapted RDM term: rdm(u,v,x,y)
                    ten2_rdm2 = ten_rdm2.copy()
                    ten2_rdm2.indices = [ten2_rdm2.indices[i] for i in [0, 1, 3, 2]]
                    const2_rdm2 = - 1.0 / 3.0

                    tens_rdm2_sa.append([ten1_rdm2, ten2_rdm2])
                    consts_rdm2_sa.append([const1_rdm2, const2_rdm2])

                elif ten_rdm2_spin_inds in [inds_abba, inds_baab]:
                    ## Spin-Adapted RDM term: rdm(u,v,y,x)
                    ten1_rdm2 = ten_rdm2.copy()
                    const1_rdm2 = 1.0 / 3.0

                    ## Spin-Adapted RDM term: rdm(u,v,x,y)
                    ten2_rdm2 = ten_rdm2.copy()
                    ten2_rdm2.indices = [ten2_rdm2.indices[i] for i in [0, 1, 3, 2]]
                    const2_rdm2 = 1.0 / 6.0

                    tens_rdm2_sa.append([ten1_rdm2, ten2_rdm2])
                    consts_rdm2_sa.append([const1_rdm2, const2_rdm2])

                else:
                    ten1_rdm2 = ten_rdm2.copy()
                    const1_rdm2 = 0.0

                    tens_rdm2_sa.append([ten1_rdm2])
                    consts_rdm2_sa.append([const1_rdm2])

                tens_rdm2_sa_permut = []
                for item in list(itertools.product(*tens_rdm2_sa)):
                    tens_rdm2_sa_permut.append(list(item))

                consts_rdm2_sa_permut = []
                for item in list(itertools.product(*consts_rdm2_sa)):
                    consts_rdm2_sa_permut.append(list(item))

                consts_rdm2_sa_prod = []
                for iter in consts_rdm2_sa_permut:
                    prod = 1.0
                    for const in iter:
                        prod = prod * const
                    consts_rdm2_sa_prod.append(prod)

                for tens_rdm2_sa_ind, tens_rdm2_sa in enumerate(tens_rdm2_sa_permut):
                    term_rdm2_sa = term_rdm2_si.copy()
                    term_rdm2_sa.scale(consts_rdm2_sa_prod[tens_rdm2_sa_ind])

                    for ten_rdm2_sa_ind, ten_rdm2_sa in zip(tens_rdm2_ind, tens_rdm2_sa):
                        term_rdm2_sa.tensors[ten_rdm2_sa_ind] = remove_spin_index_type(ten_rdm2_sa)
                        term_rdm2_sa.tensors[ten_rdm2_sa_ind].symmetries = rdm2_sa_symm
                    terms_rdm2_sa.append(term_rdm2_sa)
        else:
            terms_rdm2_sa.append(term_rdm2_si)

    termChop(terms_rdm2_sa)

    # Print 2e- spin-adapted terms
    if options.verbose:
        print("")
        for term_rdm2_sa in terms_rdm2_sa:
            print(term_rdm2_sa)
        print("")

    # Convert Three-Body RDMs
    print("Converting 3-RDMs to spin-adapted formulation...")

    # Define Spin-Adapted Three-Body RDMs Symmetries
    rdm3_sa_symm = [symmetry((1,0,2,3,5,4), 1), symmetry((0,2,1,4,3,5), 1), symmetry((5,4,3,2,1,0), 1)]

    # Define 3e- indices lists
    inds_aaaaaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type]
    inds_bbbbbb = [options.beta_type,  options.beta_type,  options.beta_type,  options.beta_type,  options.beta_type,  options.beta_type]

    inds_aabaab = [options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type]
    inds_bbabba = [options.beta_type,  options.beta_type,  options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type]

    inds_aababa = [options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type]
    inds_bbabab = [options.beta_type,  options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type]

    inds_aabbaa = [options.alpha_type, options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type, options.alpha_type]
    inds_bbaabb = [options.beta_type,  options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type,  options.beta_type]

    inds_abaaab = [options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type]
    inds_babbba = [options.beta_type,  options.alpha_type, options.beta_type,  options.beta_type,  options.beta_type,  options.alpha_type]

    inds_abaaba = [options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type]
    inds_babbab = [options.beta_type,  options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type, options.beta_type]

    inds_ababaa = [options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type]
    inds_bababb = [options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type,  options.beta_type]

    inds_baaaab = [options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type]
    inds_abbbba = [options.alpha_type, options.beta_type,  options.beta_type,  options.beta_type,  options.beta_type,  options.alpha_type]

    inds_baaaba = [options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type]
    inds_abbbab = [options.alpha_type, options.beta_type,  options.beta_type,  options.beta_type,  options.alpha_type, options.beta_type]

    inds_baabaa = [options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type]
    inds_abbabb = [options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type, options.beta_type,  options.beta_type]

    terms_rdm3_sa = []
    for term_rdm3_si in terms_rdm2_sa:
        # List for storing 3-RDMs
        tens_rdm3 = []
        tens_rdm3_ind = []

        # Append all 3-RDMs to list
        for ten_ind, ten in enumerate(term_rdm3_si.tensors):
            if isinstance(ten, creDesTensor) and len(ten.indices) == 6:
                tens_rdm3.append(ten)
                tens_rdm3_ind.append(ten_ind)

        if tens_rdm3:
            tens_rdm3_sa = []
            consts_rdm3_sa = []

            if options.verbose:
                print("\n<<< {:}".format(term_rdm3_si))

            for ten_rdm3 in tens_rdm3:
                ten_rdm3_spin_inds = [get_spin_index_type(ind) for ind in ten_rdm3.indices]

                if ten_rdm3_spin_inds in [inds_aaaaaa, inds_bbbbbb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = 1.0 / 12.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3])

                elif ten_rdm3_spin_inds in [inds_aabbaa, inds_bbaabb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,z,x,y)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 3, 5, 4]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_aababa, inds_bbabab]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,z,x)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 4, 3, 5]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_aabaab, inds_bbabba]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,y,z)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 5, 4, 3]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_ababaa, inds_bababb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,z,x)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 4, 3, 5]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_abaaba, inds_babbab]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,y,z)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 5, 4, 3]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_abaaab, inds_babbba]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,z,x,y)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 3, 5, 4]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_baabaa, inds_abbabb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,y,z)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 5, 4, 3]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_baaaba, inds_abbbab]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,z,x,y)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 3, 5, 4]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                elif ten_rdm3_spin_inds in [inds_baaaab, inds_abbbba]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten2_rdm3 = ten_rdm3.copy()
                    ten2_rdm3.indices = [ten2_rdm3.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const2_rdm3 = - 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten3_rdm3 = ten_rdm3.copy()
                    ten3_rdm3.indices = [ten3_rdm3.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const3_rdm3 = 1.0 / 12.0

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,z,x)
                    ten4_rdm3 = ten_rdm3.copy()
                    ten4_rdm3.indices = [ten4_rdm3.indices[i] for i in [0, 1, 2, 4, 3, 5]]
                    const4_rdm4 = - 1.0 / 6.0

                    tens_rdm3_sa.append([ten1_rdm3, ten2_rdm3, ten3_rdm3, ten4_rdm3])
                    consts_rdm3_sa.append([const1_rdm3, const2_rdm3, const3_rdm3, const4_rdm4])

                else:
                    ten1_rdm3 = ten_rdm3.copy()
                    const1_rdm3 = 0.0

                    tens_rdm3_sa.append([ten1_rdm3])
                    consts_rdm3_sa.append([const1_rdm3])

                tens_rdm3_sa_permut = []
                for item in list(itertools.product(*tens_rdm3_sa)):
                    tens_rdm3_sa_permut.append(list(item))

                consts_rdm3_sa_permut = []
                for item in list(itertools.product(*consts_rdm3_sa)):
                    consts_rdm3_sa_permut.append(list(item))

                consts_rdm3_sa_prod = []
                for iter in consts_rdm3_sa_permut:
                    prod = 1.0
                    for const in iter:
                        prod = prod * const
                    consts_rdm3_sa_prod.append(prod)

                for tens_rdm3_sa_ind, tens_rdm3_sa in enumerate(tens_rdm3_sa_permut):
                    term_rdm3_sa = term_rdm3_si.copy()
                    term_rdm3_sa.scale(consts_rdm3_sa_prod[tens_rdm3_sa_ind])

                    for ten_rdm3_sa_ind, ten_rdm3_sa in zip(tens_rdm3_ind, tens_rdm3_sa):
                        term_rdm3_sa.tensors[ten_rdm3_sa_ind] = remove_spin_index_type(ten_rdm3_sa)
                        term_rdm3_sa.tensors[ten_rdm3_sa_ind].symmetries = rdm3_sa_symm

                    if options.verbose:
                        print("--> {:} (factor = {:.5f})".format(term_rdm3_sa, consts_rdm3_sa_prod[tens_rdm3_sa_ind]))

                    terms_rdm3_sa.append(term_rdm3_sa)


        else:
            terms_rdm3_sa.append(term_rdm3_si)

    termChop(terms_rdm3_sa)

    # Print 3e- spin-adapted terms
    if options.verbose:
        print("")
        for term_rdm3_sa in terms_rdm3_sa:
            print(term_rdm3_sa)
        print("")

    # Convert Four-Body RDMs
    print("Converting 4-RDMs to spin-adapted formulation...")

    #TODO: Update symm
    # Define Spin-Adapted Four-Body RDMs Symmetries
    rdm4_sa_symm = [symmetry((1,0,2,3,5,4,6,7), 1), symmetry((0,2,1,3,4,6,5,7), 1), symmetry((0,1,3,2,4,5,7,6), 1),
                    symmetry((2,1,0,3,6,5,4,7), 1), symmetry((3,1,2,0,7,5,6,4), 1), symmetry((0,1,2,3,4,5,6,7), 1),
                    symmetry((0,3,2,1,4,7,6,5), 1)]

    # Define 4e- indices lists
    inds_aaaaaaaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type]

    inds_bbbbbbbb = [options.beta_type, options.beta_type, options.beta_type, options.beta_type,
                     options.beta_type, options.beta_type, options.beta_type, options.beta_type]

    inds_abbbabbb = [options.alpha_type, options.beta_type, options.beta_type, options.beta_type,
                     options.alpha_type, options.beta_type, options.beta_type, options.beta_type]

    inds_baaabaaa = [options.beta_type, options.alpha_type, options.alpha_type, options.alpha_type,
                     options.beta_type, options.alpha_type, options.alpha_type, options.alpha_type]

    inds_aabbaabb = [options.alpha_type, options.alpha_type, options.beta_type, options.beta_type,
                     options.alpha_type, options.alpha_type, options.beta_type, options.beta_type]

    inds_bbaabbaa = [options.beta_type, options.beta_type, options.alpha_type, options.alpha_type,
                     options.beta_type, options.beta_type, options.alpha_type, options.alpha_type]

    terms_rdm4_sa = []
    for term_rdm4_si in terms_rdm3_sa:
        # List for storing 4-RDMs
        tens_rdm4 = []
        tens_rdm4_ind = []

        # Append all 4-RDMs to list
        for ten_ind, ten in enumerate(term_rdm4_si.tensors):
            if isinstance(ten, creDesTensor) and len(ten.indices) == 8:
                tens_rdm4.append(ten)
                tens_rdm4_ind.append(ten_ind)

        if tens_rdm4:
            tens_rdm4_sa = []
            consts_rdm4_sa = []

            for ten_rdm4 in tens_rdm4:
                ten_rdm4_spin_inds = [get_spin_index_type(ind) for ind in ten_rdm4.indices]

                if ten_rdm4_spin_inds in [inds_aaaaaaaa, inds_bbbbbbbb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten1_rdm4 = ten_rdm4.copy()
                    ten1_rdm4.indices = [ten1_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const1_rdm4 = 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten2_rdm4 = ten_rdm4.copy()
                    ten2_rdm4.indices = [ten2_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const2_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten3_rdm4 = ten_rdm4.copy()
                    ten3_rdm4.indices = [ten3_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const3_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten4_rdm4 = ten_rdm4.copy()
                    ten4_rdm4.indices = [ten4_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const4_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten5_rdm4 = ten_rdm4.copy()
                    ten5_rdm4.indices = [ten5_rdm4.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const5_rdm4 = 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten6_rdm4 = ten_rdm4.copy()
                    ten6_rdm4.indices = [ten6_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const6_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten7_rdm4 = ten_rdm4.copy()
                    ten7_rdm4.indices = [ten7_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const7_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten8_rdm4 = ten_rdm4.copy()
                    ten8_rdm4.indices = [ten8_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const8_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten9_rdm4 = ten_rdm4.copy()
                    ten9_rdm4.indices = [ten9_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const9_rdm4 = 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten10_rdm4 = ten_rdm4.copy()
                    ten10_rdm4.indices = [ten10_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const10_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten11_rdm4 = ten_rdm4.copy()
                    ten11_rdm4.indices = [ten11_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const11_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten12_rdm4 = ten_rdm4.copy()
                    ten12_rdm4.indices = [ten12_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const12_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten13_rdm4 = ten_rdm4.copy()
                    ten13_rdm4.indices = [ten13_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const13_rdm4 = - 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten14_rdm4 = ten_rdm4.copy()
                    ten14_rdm4.indices = [ten14_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const14_rdm4 = 1.0 / 30.0

                    tens_rdm4_sa.append([ten1_rdm4,  ten2_rdm4,  ten3_rdm4,  ten4_rdm4, ten5_rdm4,
                                         ten6_rdm4,  ten7_rdm4,  ten8_rdm4,  ten9_rdm4, ten10_rdm4,
                                         ten11_rdm4, ten12_rdm4, ten13_rdm4, ten14_rdm4])

                    consts_rdm4_sa.append([const1_rdm4, const2_rdm4, const3_rdm4, const4_rdm4, const5_rdm4,
                                           const6_rdm4, const7_rdm4, const8_rdm4, const9_rdm4, const10_rdm4,
                                           const11_rdm4, const12_rdm4, const13_rdm4, const14_rdm4])

                elif ten_rdm4_spin_inds in [inds_abbbabbb, inds_baaabaaa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten1_rdm4 = ten_rdm4.copy()
                    ten1_rdm4.indices = [ten1_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const1_rdm4 = - 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten2_rdm4 = ten_rdm4.copy()
                    ten2_rdm4.indices = [ten2_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const2_rdm4 = - 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten3_rdm4 = ten_rdm4.copy()
                    ten3_rdm4.indices = [ten3_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const3_rdm4 = - 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten4_rdm4 = ten_rdm4.copy()
                    ten4_rdm4.indices = [ten4_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const4_rdm4 = - 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten5_rdm4 = ten_rdm4.copy()
                    ten5_rdm4.indices = [ten5_rdm4.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const5_rdm4 = - 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten6_rdm4 = ten_rdm4.copy()
                    ten6_rdm4.indices = [ten6_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const6_rdm4 = - 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten7_rdm4 = ten_rdm4.copy()
                    ten7_rdm4.indices = [ten7_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const7_rdm4 = - 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten8_rdm4 = ten_rdm4.copy()
                    ten8_rdm4.indices = [ten8_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const8_rdm4 = - 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten9_rdm4 = ten_rdm4.copy()
                    ten9_rdm4.indices = [ten9_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const9_rdm4 = - 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten10_rdm4 = ten_rdm4.copy()
                    ten10_rdm4.indices = [ten10_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const10_rdm4 = - 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten11_rdm4 = ten_rdm4.copy()
                    ten11_rdm4.indices = [ten11_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const11_rdm4 = 1.0 / 15.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten12_rdm4 = ten_rdm4.copy()
                    ten12_rdm4.indices = [ten12_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const12_rdm4 = 1.0 / 15.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten13_rdm4 = ten_rdm4.copy()
                    ten13_rdm4.indices = [ten13_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const13_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten14_rdm4 = ten_rdm4.copy()
                    ten14_rdm4.indices = [ten14_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const14_rdm4 = 1.0 / 20.0

                    tens_rdm4_sa.append([ten1_rdm4,  ten2_rdm4,  ten3_rdm4,  ten4_rdm4, ten5_rdm4,
                                         ten6_rdm4,  ten7_rdm4,  ten8_rdm4,  ten9_rdm4, ten10_rdm4,
                                         ten11_rdm4, ten12_rdm4, ten13_rdm4, ten14_rdm4])

                    consts_rdm4_sa.append([const1_rdm4, const2_rdm4, const3_rdm4, const4_rdm4, const5_rdm4,
                                           const6_rdm4, const7_rdm4, const8_rdm4, const9_rdm4, const10_rdm4,
                                           const11_rdm4, const12_rdm4, const13_rdm4, const14_rdm4])

                elif ten_rdm4_spin_inds in [inds_aabbaabb, inds_bbaabbaa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten1_rdm4 = ten_rdm4.copy()
                    ten1_rdm4.indices = [ten1_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const1_rdm4 = - 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten2_rdm4 = ten_rdm4.copy()
                    ten2_rdm4.indices = [ten2_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const2_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten3_rdm4 = ten_rdm4.copy()
                    ten3_rdm4.indices = [ten3_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const3_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten4_rdm4 = ten_rdm4.copy()
                    ten4_rdm4.indices = [ten4_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const4_rdm4 = - 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten5_rdm4 = ten_rdm4.copy()
                    ten5_rdm4.indices = [ten5_rdm4.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const5_rdm4 = - 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten6_rdm4 = ten_rdm4.copy()
                    ten6_rdm4.indices = [ten6_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const6_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten7_rdm4 = ten_rdm4.copy()
                    ten7_rdm4.indices = [ten7_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const7_rdm4 = - 1.0 / 20.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten8_rdm4 = ten_rdm4.copy()
                    ten8_rdm4.indices = [ten8_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const8_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten9_rdm4 = ten_rdm4.copy()
                    ten9_rdm4.indices = [ten9_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const9_rdm4 = - 7.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten10_rdm4 = ten_rdm4.copy()
                    ten10_rdm4.indices = [ten10_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const10_rdm4 = 1.0 / 30.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten11_rdm4 = ten_rdm4.copy()
                    ten11_rdm4.indices = [ten11_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const11_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten12_rdm4 = ten_rdm4.copy()
                    ten12_rdm4.indices = [ten12_rdm4.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const12_rdm4 = 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten13_rdm4 = ten_rdm4.copy()
                    ten13_rdm4.indices = [ten13_rdm4.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const13_rdm4 = - 1.0 / 60.0

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten14_rdm4 = ten_rdm4.copy()
                    ten14_rdm4.indices = [ten14_rdm4.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const14_rdm4 = 1.0 / 30.0

                    tens_rdm4_sa.append([ten1_rdm4,  ten2_rdm4,  ten3_rdm4,  ten4_rdm4, ten5_rdm4,
                                         ten6_rdm4,  ten7_rdm4,  ten8_rdm4,  ten9_rdm4, ten10_rdm4,
                                         ten11_rdm4, ten12_rdm4, ten13_rdm4, ten14_rdm4])

                    consts_rdm4_sa.append([const1_rdm4, const2_rdm4, const3_rdm4, const4_rdm4, const5_rdm4,
                                           const6_rdm4, const7_rdm4, const8_rdm4, const9_rdm4, const10_rdm4,
                                           const11_rdm4, const12_rdm4, const13_rdm4, const14_rdm4])

                else:
                    ten1_rdm4 = ten_rdm4.copy()
                    const1_rdm4 = 0.0

                    tens_rdm4_sa.append([ten1_rdm4])
                    consts_rdm4_sa.append([const1_rdm4])

                tens_rdm4_sa_permut = []
                for item in list(itertools.product(*tens_rdm4_sa)):
                    tens_rdm4_sa_permut.append(list(item))

                consts_rdm4_sa_permut = []
                for item in list(itertools.product(*consts_rdm4_sa)):
                    consts_rdm4_sa_permut.append(list(item))

                consts_rdm4_sa_prod = []
                for iter in consts_rdm4_sa_permut:
                    prod = 1.0
                    for const in iter:
                        prod = prod * const
                    consts_rdm4_sa_prod.append(prod)

                for tens_rdm4_sa_ind, tens_rdm4_sa in enumerate(tens_rdm4_sa_permut):
                    term_rdm4_sa = term_rdm4_si.copy()
                    term_rdm4_sa.scale(consts_rdm4_sa_prod[tens_rdm4_sa_ind])

                    for ten_rdm4_sa_ind, ten_rdm4_sa in zip(tens_rdm4_ind, tens_rdm4_sa):
                        term_rdm4_sa.tensors[ten_rdm4_sa_ind] = remove_spin_index_type(ten_rdm4_sa)
                        term_rdm4_sa.tensors[ten_rdm4_sa_ind].symmetries = rdm4_sa_symm
                    terms_rdm4_sa.append(term_rdm4_sa)
        else:
            terms_rdm4_sa.append(term_rdm4_si)

    termChop(terms_rdm4_sa)

    # Print 4e- spin-adapted terms
    if options.verbose:
        print("")
        for term_rdm4_sa in terms_rdm4_sa:
            print(term_rdm4_sa)
        print("")

    print("Done!")
    return terms_rdm4_sa

def convert_kdelta_si_to_sa(_terms_kdelta_si):
    "Convert Kronecker Delta Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting Kronecker Deltas to spin-adapted formulation...")

    # Converting kdelta objects in each term
    terms_kdelta_sa = []
    for term_kdelta_si in _terms_kdelta_si:
        term_kdelta_sa = term_kdelta_si.copy()
        for ten_ind, ten in enumerate(term_kdelta_si.tensors):
            if isinstance(ten, kroneckerDelta):
                term_kdelta_sa.tensors[ten_ind] = remove_spin_index_type(ten)
        terms_kdelta_sa.append(term_kdelta_sa)

    # Print kdelta spin-adapted terms
    if options.verbose:
        print("")
        for term_kdelta_sa in terms_kdelta_sa:
            print(term_kdelta_sa)
        print("")

    print("Done!")
    return terms_kdelta_sa

def convert_e_si_to_sa(_terms_e_si):
    "Convert Eigenvalue Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting eigenvalues to spin-adapted formulation...")

    # Converting eigenvalue objects in each term
    terms_e_sa = []
    for term_e_si in _terms_e_si:
        term_e_sa = term_e_si.copy()
        for ten_ind, ten in enumerate(term_e_si.tensors):
            if ten.name == 'e' and len(ten.indices) == 1:
                term_e_sa.tensors[ten_ind] = remove_spin_index_type(ten)
        terms_e_sa.append(term_e_sa)

    # Print e spin-adapted terms
    if options.verbose:
        print("")
        for term_e_sa in terms_e_sa:
            print(term_e_sa)
        print("")

    print("Done!")
    return terms_e_sa

def convert_t_amplitudes_si_to_sa(_terms_t_si):
    "Convert T Amplitudes Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting T amplitudes to spin-adapted formulation...")

    # Define Spin-Adapted Amplitudes Symmetries
    t1_sa_symm = [symmetry((1,0), 1)]
    t2_sa_symm = [symmetry((1,0,3,2), 1), symmetry((2,3,0,1), 1)]

    # Define 1e- indices lists
    inds_aa = [options.alpha_type, options.alpha_type]
    inds_bb = [options.beta_type,  options.beta_type]

    # Convert One-Body Amplitudes
    terms_t1_sa = []
    for term_t1_si in _terms_t_si:
        term_t1_sa = term_t1_si.copy()

        for ten_ind, ten in enumerate(term_t1_sa.tensors):
            if ten.name[0] == 't' and len(ten.indices) == 2:
                ten_t1_spin_inds = [get_spin_index_type(ind) for ind in ten.indices]

                if ten_t1_spin_inds in [inds_aa, inds_bb]:
                    ten_t1 = ten.copy()
                    term_t1_sa.tensors[ten_ind] = remove_spin_index_type(ten_t1)
                    term_t1_sa.tensors[ten_ind].symmetries = t1_sa_symm
                else:
                    term_t1_sa.scale(0.0)
                    term_t1_sa.tensors[ten_ind] = remove_spin_index_type(ten)

        terms_t1_sa.append(term_t1_sa)

    termChop(terms_t1_sa)

    # Define 2e- indices lists
    inds_aaaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type]
    inds_bbbb = [options.beta_type,  options.beta_type,  options.beta_type,  options.beta_type]

    inds_abab = [options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type]
    inds_baba = [options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type]
    inds_abba = [options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type]
    inds_baab = [options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type]

    # Convert Two-Body Amplitudes
    terms_t2_sa = []
    for term_t2_si in terms_t1_sa:
        # List for storing Amplitudes
        tens_t2 = []
        tens_t2_ind = []

        for ten_ind, ten in enumerate(term_t2_si.tensors):
            if ten.name[0] == 't' and len(ten.indices) == 4:
                tens_t2.append(ten)
                tens_t2_ind.append(ten_ind)

        ## Convert 2e- object
        if tens_t2:
            tens_t2_sa = []
            consts_t2_sa = []

            for ten_t2 in tens_t2:

                ten_t2_inds = [get_spatial_index_type(ind) for ind in ten_t2.indices]
                ten_t2_spin_inds = [get_spin_index_type(ind) for ind in ten_t2.indices]

                if (ten_t2_inds[0] == ten_t2_inds[1]) and (ten_t2_inds[2] == ten_t2_inds[3]):

                    if ten_t2_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        ## Spin-Adapted 2e- term: t(p,q,s,r)
                        ten2_t2 = ten_t2.copy()
                        ten2_t2.indices = [ten2_t2.indices[i] for i in [0, 1, 3, 2]]
                        const2_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2, ten2_t2])
                        consts_t2_sa.append([const1_t2, const2_t2])

                    elif ten_t2_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    elif ten_t2_spin_inds in [inds_abba, inds_baab]:
                        ten1_t2 = ten_t2.copy()
                        ten1_t2.indices = [ten1_t2.indices[i] for i in [0, 1, 3, 2]]
                        const1_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    else:
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 0.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                elif (ten_t2_inds[0] != ten_t2_inds[1]) and (ten_t2_inds[2] == ten_t2_inds[3]):

                    if ten_t2_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        ## Spin-Adapted 2e- term: t(p,q,s,r)
                        ten2_t2 = ten_t2.copy()
                        ten2_t2.indices = [ten2_t2.indices[i] for i in [0, 1, 3, 2]]
                        const2_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2, ten2_t2])
                        consts_t2_sa.append([const1_t2, const2_t2])

                    elif ten_t2_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    elif ten_t2_spin_inds in [inds_abba, inds_baab]:
                        ten1_t2 = ten_t2.copy()
                        ten1_t2.indices = [ten1_t2.indices[i] for i in [0, 1, 3, 2]]
                        const1_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    else:
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 0.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                elif (ten_t2_inds[0] == ten_t2_inds[1]) and (ten_t2_inds[2] != ten_t2_inds[3]):

                    if ten_t2_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        ## Spin-Adapted 2e- term: t(p,q,s,r)
                        ten2_t2 = ten_t2.copy()
                        ten2_t2.indices = [ten2_t2.indices[i] for i in [1, 0, 2, 3]]
                        const2_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2, ten2_t2])
                        consts_t2_sa.append([const1_t2, const2_t2])

                    elif ten_t2_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    elif ten_t2_spin_inds in [inds_abba, inds_baab]:
                        ten1_t2 = ten_t2.copy()
                        ten1_t2.indices = [ten1_t2.indices[i] for i in [1, 0, 2, 3]]
                        const1_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    else:
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 0.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                elif ((ten_t2_inds[0] != ten_t2_inds[1]) and (ten_t2_inds[2] != ten_t2_inds[3])):

                    if ten_t2_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        ## Spin-Adapted 2e- term: t(p,q,s,r)
                        ten2_t2 = ten_t2.copy()
                        ten2_t2.indices = [ten2_t2.indices[i] for i in [0, 1, 3, 2]]
                        const2_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2, ten2_t2])
                        consts_t2_sa.append([const1_t2, const2_t2])

                    elif ten_t2_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    elif ten_t2_spin_inds in [inds_abba, inds_baab]:
                        ten1_t2 = ten_t2.copy()
                        ten1_t2.indices = [ten1_t2.indices[i] for i in [0, 1, 3, 2]]
                        const1_t2 = - 1.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

                    else:
                        ten1_t2 = ten_t2.copy()
                        const1_t2 = 0.0

                        tens_t2_sa.append([ten1_t2])
                        consts_t2_sa.append([const1_t2])

            tens_t2_sa_permut = []
            for item in list(itertools.product(*tens_t2_sa)):
                tens_t2_sa_permut.append(list(item))

            consts_t2_sa_permut = []
            for item in list(itertools.product(*consts_t2_sa)):
                consts_t2_sa_permut.append(list(item))

            consts_t2_sa_prod = []
            for iter in consts_t2_sa_permut:
                prod = 1.0
                for const in iter:
                    prod = prod * const
                consts_t2_sa_prod.append(prod)

            for tens_t2_sa_ind, tens_t2_sa in enumerate(tens_t2_sa_permut):
                term_t2_sa = term_t2_si.copy()
                term_t2_sa.scale(consts_t2_sa_prod[tens_t2_sa_ind])

                for ten_t2_sa_ind, ten_t2_sa in zip(tens_t2_ind, tens_t2_sa):
                    term_t2_sa.tensors[ten_t2_sa_ind] = remove_spin_index_type(ten_t2_sa)
                    term_t2_sa.tensors[ten_t2_sa_ind].symmetries = t2_sa_symm
                terms_t2_sa.append(term_t2_sa)

        else:
            terms_t2_sa.append(term_t2_si)

    termChop(terms_t2_sa)

    # Print kdelta spin-adapted terms
    if options.verbose:
        print("")
        for term_t_sa in terms_t2_sa:
            print(term_t_sa)
        print("")

    print("Done!")
    return terms_t2_sa
