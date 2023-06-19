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
from sqaIndex import get_spatial_index_type, get_spin_index_type, is_alpha_index_type, is_beta_index_type, is_cvs_index_type
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

    dummyLabel(terms_si)
    len_terms_si = len(terms_si)

    # Temporarily remove spin-integrated symmetries of non-RDM tensors
    remove_si_tensors_symmetries(terms_si)

    # Convert 1e- integrals to Spin-Adapted Formulation
    terms_sa = convert_h1e_si_to_sa(terms_si)

    # Convert 2e- integrals to Spin-Adapted Formulation
    terms_sa = convert_v2e_si_to_sa(terms_sa)

    # Convert T amplitudes to Spin-Adapted Formulation
    terms_sa = convert_t_amplitudes_si_to_sa(terms_sa)

   # Convert RDMs to Canonical Form before Spin-Adaptation
    for term_sa in terms_sa:
        has_high_rdms = False
        for tensor_sa in term_sa.tensors:
            if isinstance(tensor_sa, creDesTensor) and len(tensor_sa.indices) in [6, 8]:
                has_high_rdms = True
        
        if has_high_rdms:
            term_sa.isInCanonicalForm = False
            term_sa.makeCanonical()
    dummyLabel(terms_sa)

    # Convert RDMs to Spin-Adapted Formulation
    terms_sa = convert_rdms_si_to_sa(terms_sa)

    # Update Spin-Adapted Symmetries in tensors
    update_sa_tensors_symmetries(terms_sa)

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
        v2e_symm = [symmetry((1,0,3,2), 1), symmetry((2,3,0,1), 1)]

        ## Append all v2e objects to list
        for ten_v2e_ind, ten_v2e in enumerate(term_v2e.tensors):
            if (ten_v2e.name == 'v') and (len(ten_v2e.indices) == 4):
                v2e_tensors.append(ten_v2e)
                v2e_tensors_ind.append(ten_v2e_ind)

        ## Convert from Physicists's to Chemists's Notation
        if v2e_tensors:
            for ten_v2e_ind, ten_v2e in zip(v2e_tensors_ind, v2e_tensors):

                if options.verbose:
                    unordered_ten_v2e = ten_v2e.copy()

                # v2e[p,q,r,s] => v2e[p,r,q,s]
                ten_v2e.indices = [ten_v2e.indices[i] for i in [0, 2, 1, 3]]
                ten_v2e.symmetries = v2e_symm

                if options.verbose:
                    print("{:}    --->    {:}".format(unordered_ten_v2e, ten_v2e))

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
            if (isinstance(ten_rdm, creDesTensor)) and (len(ten_rdm.indices) > 2):
                rdm_tensors.append(ten_rdm)
                rdm_tensors_ind.append(ten_rdm_ind)

        ## Convert from Physicists's to Chemists's Notation
        if rdm_tensors:
            for ten_rdm_ind, ten_rdm in zip(rdm_tensors_ind, rdm_tensors):

                if options.verbose:
                    unordered_tensor = ten_rdm.copy()

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

                if options.verbose:
                    print("{:}    --->    {:}".format(unordered_tensor, ten_rdm))

                _terms_rdm[term_rdm_ind].tensors[ten_rdm_ind] = ten_rdm.copy()

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
        for ten in term_h1e_si.tensors:
            if ten.name == 'h' and len(ten.indices) == 2:
                ten_h1e = ten

        if ten_h1e:
            ten_h1e_spin_inds = [get_spin_index_type(ind) for ind in ten_h1e.indices]

            if ten_h1e_spin_inds in [inds_aa, inds_bb]:
                term_h1e_sa = term_h1e_si.copy()
                terms_h1e_sa.append(term_h1e_sa)

            else:
                term_h1e_sa = term_h1e_si.copy()
                term_h1e_sa.scale(0.0)
                terms_h1e_sa.append(term_h1e_sa)
        else:
            terms_h1e_sa.append(term_h1e_si)

    termChop(terms_h1e_sa)

    print("Done!")
    return terms_h1e_sa

def convert_v2e_si_to_sa(_terms_v2e_si):
    "Convert v2e Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting 2e- integrals to spin-adapted formulation...")

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

            if options.verbose:
                print("\n<<< {:}".format(term_v2e_si))

            for ten_v2e in tens_v2e:

                ten_v2e_inds = [get_spatial_index_type(ind) for ind in ten_v2e.indices]
                ten_v2e_spin_inds = [get_spin_index_type(ind) for ind in ten_v2e.indices]

                ten_v2e_tens_sa = []
                const_v2e_tens_sa = []

                if ((ten_v2e_inds[0] == ten_v2e_inds[1] == ten_v2e_inds[2] == ten_v2e_inds[3]) or

                   (((ten_v2e_inds[0] == ten_v2e_inds[1]) and (ten_v2e_inds[2] == ten_v2e_inds[3]) and
                     (ten_v2e_inds[0] != ten_v2e_inds[2]) and (ten_v2e_inds[1] != ten_v2e_inds[3])) or

                    ((ten_v2e_inds[0] != ten_v2e_inds[1]) and (ten_v2e_inds[2] == ten_v2e_inds[3])))):

                    if ten_v2e_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten_v2e_sa = ten_v2e.copy()
                        const_v2e_sa = 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)
                
                        ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                        ten_v2e_sa = ten_v2e.copy()
                        ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [0, 1, 3, 2]]
                        const_v2e_sa = - 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                    elif ten_v2e_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten_v2e_sa = ten_v2e.copy()
                        const_v2e_sa = 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                    elif ten_v2e_spin_inds in [inds_abba, inds_baab]:
                        ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                        ten_v2e_sa = ten_v2e.copy()
                        ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [0, 1, 3, 2]]
                        const_v2e_sa = - 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                    else:
                        ten_v2e_sa = ten_v2e.copy()
                        const_v2e_sa = 0.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                elif (ten_v2e_inds[0] != ten_v2e_inds[1]) and (ten_v2e_inds[2] != ten_v2e_inds[3]):

                    if options.cvs_approach:
                        ten_v2e_cvs_inds = [is_cvs_index_type(ind) for ind in ten_v2e.indices]

                    if options.cvs_approach and (ten_v2e_cvs_inds[0] == ten_v2e_cvs_inds[1]):
                        if ten_v2e_spin_inds in [inds_aaaa, inds_bbbb]:
                            ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                            ten_v2e_sa = ten_v2e.copy()
                            const_v2e_sa = 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                            ## Spin-Adapted 2e- term: v2e(q,p,r,s)
                            ten_v2e_sa = ten_v2e.copy()
                            ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [1, 0, 2, 3]]
                            const_v2e_sa = - 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                        elif ten_v2e_spin_inds in [inds_abab, inds_baba]:
                            ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                            ten_v2e_sa = ten_v2e.copy()
                            const_v2e_sa = 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                        elif ten_v2e_spin_inds in [inds_abba, inds_baab]:
                            ## Spin-Adapted 2e- term: v2e(q,p,r,s)
                            ten_v2e_sa = ten_v2e.copy()
                            ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [1, 0, 2, 3]]
                            const_v2e_sa = - 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                        else:
                            ten_v2e_sa = ten_v2e.copy()
                            const_v2e_sa = 0.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                    else:
                        if ten_v2e_spin_inds in [inds_aaaa, inds_bbbb]:
                            ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                            ten_v2e_sa = ten_v2e.copy()
                            const_v2e_sa = 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                            ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                            ten_v2e_sa = ten_v2e.copy()
                            ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [0, 1, 3, 2]]
                            const_v2e_sa = - 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                        elif ten_v2e_spin_inds in [inds_abab, inds_baba]:
                            ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                            ten_v2e_sa = ten_v2e.copy()
                            const_v2e_sa = 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                        elif ten_v2e_spin_inds in [inds_abba, inds_baab]:
                            ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                            ten_v2e_sa = ten_v2e.copy()
                            ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [0, 1, 3, 2]]
                            const_v2e_sa = - 1.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                        else:
                            ten_v2e_sa = ten_v2e.copy()
                            const_v2e_sa = 0.0

                            ten_v2e_tens_sa.append(ten_v2e_sa)
                            const_v2e_tens_sa.append(const_v2e_sa)

                elif ((ten_v2e_inds[0] == ten_v2e_inds[1]) and (ten_v2e_inds[2] != ten_v2e_inds[3])):
                    if ten_v2e_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten_v2e_sa = ten_v2e.copy()
                        const_v2e_sa = 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                        ## Spin-Adapted 2e- term: v2e(q,p,r,s)
                        ten_v2e_sa = ten_v2e.copy()
                        ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [1, 0, 2, 3]]
                        const_v2e_sa = - 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                    elif ten_v2e_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten_v2e_sa = ten_v2e.copy()
                        const_v2e_sa = 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                    elif ten_v2e_spin_inds in [inds_abba, inds_baab]:
                        ## Spin-Adapted 2e- term: v2e(q,p,r,s)
                        ten_v2e_sa = ten_v2e.copy()
                        ten_v2e_sa.indices = [ten_v2e_sa.indices[i] for i in [1, 0, 2, 3]]
                        const_v2e_sa = - 1.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                    else:
                        ten_v2e_sa = ten_v2e.copy()
                        const_v2e_sa = 0.0

                        ten_v2e_tens_sa.append(ten_v2e_sa)
                        const_v2e_tens_sa.append(const_v2e_sa)

                tens_v2e_sa.append(ten_v2e_tens_sa)
                consts_v2e_sa.append(const_v2e_tens_sa)

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
                    term_v2e_sa.tensors[ten_v2e_sa_ind] = ten_v2e_sa

                if options.verbose:
                    print("--> {:} (factor = {:.5f})".format(term_v2e_sa, consts_v2e_sa_prod[tens_v2e_sa_ind]))

                terms_v2e_sa.append(term_v2e_sa)

        else:
            terms_v2e_sa.append(term_v2e_si)

    termChop(terms_v2e_sa)

    print("Done!")
    return terms_v2e_sa

def convert_rdms_si_to_sa(_terms_rdm_si):
    "Convert RDM Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting RDMs to spin-adapted formulation...\n")

    # Convert One-Body RDMs
    print("Converting 1-RDMs to spin-adapted formulation...")

    # Define 1e- indices lists
    inds_aa = [options.alpha_type, options.alpha_type]
    inds_bb = [options.beta_type,  options.beta_type]

    terms_rdm1_sa = []
    for term_rdm1_si in _terms_rdm_si:
        term_rdm1_sa = term_rdm1_si.copy()
        for ten_ind, ten in enumerate(term_rdm1_si.tensors):
            if isinstance(ten, creDesTensor) and len(ten.indices) == 2:
                ten_rdm1_spin_inds = [get_spin_index_type(ind) for ind in ten.indices]

                if options.verbose:
                    print("\n<<< {:}".format(term_rdm1_si))

                if ten_rdm1_spin_inds in [inds_aa, inds_bb]:
                    consts_rdm1_sa_prod = 1.0 / 2.0 
                    term_rdm1_sa.scale(consts_rdm1_sa_prod)

                else:
                    consts_rdm1_sa_prod = 0.0
                    term_rdm1_sa.scale(consts_rdm1_sa_prod)

                if options.verbose:
                    print("--> {:} (factor = {:.5f})".format(term_rdm1_sa, consts_rdm1_sa_prod))

        terms_rdm1_sa.append(term_rdm1_sa)

    termChop(terms_rdm1_sa)

    # Convert Two-Body RDMs
    print("Converting 2-RDMs to spin-adapted formulation...")

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

            if options.verbose:
                print("\n<<< {:}".format(term_rdm2_si))

            for ten_rdm2 in tens_rdm2:
                ten_rdm2_spin_inds = [get_spin_index_type(ind) for ind in ten_rdm2.indices]

                ten_rdm2_tens_sa = []
                const_rdm2_tens_sa = []

                if ten_rdm2_spin_inds in [inds_aaaa, inds_bbbb]:
                    ## Spin-Adapted RDM term: rdm(u,v,y,x)
                    ten_rdm2_sa = ten_rdm2.copy()
                    const_rdm2_sa = 1.0 / 6.0

                    ten_rdm2_tens_sa.append(ten_rdm2_sa)
                    const_rdm2_tens_sa.append(const_rdm2_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,x,y)
                    ten_rdm2_sa = ten_rdm2.copy()
                    ten_rdm2_sa.indices = [ten_rdm2_sa.indices[i] for i in [0, 1, 3, 2]]
                    const_rdm2_sa = - 1.0 / 6.0

                    ten_rdm2_tens_sa.append(ten_rdm2_sa)
                    const_rdm2_tens_sa.append(const_rdm2_sa)

                elif ten_rdm2_spin_inds in [inds_abab, inds_baba]:
                    ## Spin-Adapted RDM term: rdm(u,v,y,x)
                    ten_rdm2_sa = ten_rdm2.copy()
                    const_rdm2_sa = - 1.0 / 6.0

                    ten_rdm2_tens_sa.append(ten_rdm2_sa)
                    const_rdm2_tens_sa.append(const_rdm2_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,x,y)
                    ten_rdm2_sa = ten_rdm2.copy()
                    ten_rdm2_sa.indices = [ten_rdm2_sa.indices[i] for i in [0, 1, 3, 2]]
                    const_rdm2_sa = - 1.0 / 3.0

                    ten_rdm2_tens_sa.append(ten_rdm2_sa)
                    const_rdm2_tens_sa.append(const_rdm2_sa)

                elif ten_rdm2_spin_inds in [inds_abba, inds_baab]:
                    ## Spin-Adapted RDM term: rdm(u,v,y,x)
                    ten_rdm2_sa = ten_rdm2.copy()
                    const_rdm2_sa = 1.0 / 3.0

                    ten_rdm2_tens_sa.append(ten_rdm2_sa)
                    const_rdm2_tens_sa.append(const_rdm2_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,x,y)
                    ten_rdm2_sa = ten_rdm2.copy()
                    ten_rdm2_sa.indices = [ten_rdm2_sa.indices[i] for i in [0, 1, 3, 2]]
                    const_rdm2_sa = 1.0 / 6.0

                    ten_rdm2_tens_sa.append(ten_rdm2_sa)
                    const_rdm2_tens_sa.append(const_rdm2_sa)

                else:
                    ten_rdm2_sa = ten_rdm2.copy()
                    const_rdm2_sa = 0.0

                    ten_rdm2_tens_sa.append(ten_rdm2_sa)
                    const_rdm2_tens_sa.append(const_rdm2_sa)

                tens_rdm2_sa.append(ten_rdm2_tens_sa)
                consts_rdm2_sa.append(const_rdm2_tens_sa)

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
                    term_rdm2_sa.tensors[ten_rdm2_sa_ind] = ten_rdm2_sa
                    term_rdm2_sa.tensors[ten_rdm2_sa_ind].symmetries = []

                if options.verbose:
                    print("--> {:} (factor = {:.5f})".format(term_rdm2_sa, consts_rdm2_sa_prod[tens_rdm2_sa_ind]))

                terms_rdm2_sa.append(term_rdm2_sa)

        else:
            terms_rdm2_sa.append(term_rdm2_si)

    termChop(terms_rdm2_sa)

    # Convert Three-Body RDMs
    print("Converting 3-RDMs to spin-adapted formulation...")

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

                ten_rdm3_tens_sa = []
                const_rdm3_tens_sa = []

                if ten_rdm3_spin_inds in [inds_aaaaaa, inds_bbbbbb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_aabbaa, inds_bbaabb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,z,x,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 3, 5, 4]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_aababa, inds_bbabab]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,z,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 3, 5]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_aabaab, inds_bbabba]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,y,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 4, 3]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_ababaa, inds_bababb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,z,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 3, 5]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_abaaba, inds_babbab]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,y,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 4, 3]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_abaaab, inds_babbba]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,z,x,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 3, 5, 4]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_baabaa, inds_abbabb]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,y,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 4, 3]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_baaaba, inds_abbbab]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,z,x,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 3, 5, 4]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                elif ten_rdm3_spin_inds in [inds_baaaab, inds_abbbba]:
                    ## Spin-Adapted RDM term: rdm(u,v,w,z,y,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,x,z,y)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 5, 3, 4]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,x,z)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 5, 3]]
                    const_rdm3_sa = - 1.0 / 12.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                    ## Spin-Adapted RDM term: rdm(u,v,w,y,z,x)
                    ten_rdm3_sa = ten_rdm3.copy()
                    ten_rdm3_sa.indices = [ten_rdm3_sa.indices[i] for i in [0, 1, 2, 4, 3, 5]]
                    const_rdm3_sa = - 1.0 / 6.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                else:
                    ten_rdm3_sa = ten_rdm3.copy()
                    const_rdm3_sa = 0.0

                    ten_rdm3_tens_sa.append(ten_rdm3_sa)
                    const_rdm3_tens_sa.append(const_rdm3_sa)

                tens_rdm3_sa.append(ten_rdm3_tens_sa)
                consts_rdm3_sa.append(const_rdm3_tens_sa)

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
                    term_rdm3_sa.tensors[ten_rdm3_sa_ind] = ten_rdm3_sa
                    term_rdm3_sa.tensors[ten_rdm3_sa_ind].symmetries = []

                if options.verbose:
                    print("--> {:} (factor = {:.5f})".format(term_rdm3_sa, consts_rdm3_sa_prod[tens_rdm3_sa_ind]))

                terms_rdm3_sa.append(term_rdm3_sa)

        else:
            terms_rdm3_sa.append(term_rdm3_si)

    termChop(terms_rdm3_sa)

    # Convert Four-Body RDMs
    print("Converting 4-RDMs to spin-adapted formulation...")

    # Define 4e- indices lists
    inds_aaaaaaaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.alpha_type, options.alpha_type]

    inds_bbbbbbbb = [options.beta_type, options.beta_type, options.beta_type, options.beta_type,
                     options.beta_type, options.beta_type, options.beta_type, options.beta_type]

    inds_bbbaabbb = [options.beta_type,  options.beta_type, options.beta_type, options.alpha_type,
                     options.alpha_type, options.beta_type, options.beta_type, options.beta_type]

    inds_bbbababb = [options.beta_type, options.beta_type,  options.beta_type, options.alpha_type,
                     options.beta_type, options.alpha_type, options.beta_type, options.beta_type]

    inds_bbbabbab = [options.beta_type, options.beta_type, options.beta_type,  options.alpha_type,
                     options.beta_type, options.beta_type, options.alpha_type, options.beta_type]

    inds_bbbabbba = [options.beta_type, options.beta_type, options.beta_type, options.alpha_type,
                     options.beta_type, options.beta_type, options.beta_type, options.alpha_type]

    inds_bbababbb = [options.beta_type,  options.beta_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.beta_type, options.beta_type,  options.beta_type]

    inds_bbabbabb = [options.beta_type, options.beta_type,  options.alpha_type, options.beta_type,
                     options.beta_type, options.alpha_type, options.beta_type,  options.beta_type]

    inds_bbabbbab = [options.beta_type, options.beta_type, options.alpha_type, options.beta_type,
                     options.beta_type, options.beta_type, options.alpha_type, options.beta_type]

    inds_bbabbbba = [options.beta_type, options.beta_type, options.alpha_type, options.beta_type,
                     options.beta_type, options.beta_type, options.beta_type,  options.alpha_type]

    inds_babbabbb = [options.beta_type,  options.alpha_type, options.beta_type, options.beta_type,
                     options.alpha_type, options.beta_type,  options.beta_type, options.beta_type]

    inds_babbbabb = [options.beta_type, options.alpha_type, options.beta_type, options.beta_type,
                     options.beta_type, options.alpha_type, options.beta_type, options.beta_type]

    inds_babbbbab = [options.beta_type, options.alpha_type, options.beta_type,  options.beta_type,
                     options.beta_type, options.beta_type,  options.alpha_type, options.beta_type]

    inds_babbbbba = [options.beta_type, options.alpha_type, options.beta_type, options.beta_type,
                     options.beta_type, options.beta_type,  options.beta_type, options.alpha_type]

    inds_abbbabbb = [options.alpha_type, options.beta_type, options.beta_type, options.beta_type,
                     options.alpha_type, options.beta_type, options.beta_type, options.beta_type]

    inds_abbbbabb = [options.alpha_type, options.beta_type,  options.beta_type, options.beta_type,
                     options.beta_type,  options.alpha_type, options.beta_type, options.beta_type]

    inds_abbbbbab = [options.alpha_type, options.beta_type, options.beta_type,  options.beta_type,
                     options.beta_type,  options.beta_type, options.alpha_type, options.beta_type]

    inds_abbbbbba = [options.alpha_type, options.beta_type, options.beta_type, options.beta_type,
                     options.beta_type,  options.beta_type, options.beta_type, options.alpha_type]

    inds_bbaaaabb = [options.beta_type,  options.beta_type,  options.alpha_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.beta_type,  options.beta_type]

    inds_bbaaabab = [options.beta_type,  options.beta_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.beta_type, options.alpha_type, options.beta_type]

    inds_bbaaabba = [options.beta_type,  options.beta_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.beta_type, options.beta_type,  options.alpha_type]

    inds_bbaabaab = [options.beta_type, options.beta_type,  options.alpha_type, options.alpha_type,
                     options.beta_type, options.alpha_type, options.alpha_type, options.beta_type]

    inds_bbaababa = [options.beta_type, options.beta_type,  options.alpha_type, options.alpha_type,
                     options.beta_type, options.alpha_type, options.beta_type,  options.alpha_type]

    inds_bbaabbaa = [options.beta_type, options.beta_type, options.alpha_type, options.alpha_type,
                     options.beta_type, options.beta_type, options.alpha_type, options.alpha_type]

    inds_babaaabb = [options.beta_type,  options.alpha_type, options.beta_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.beta_type, options.beta_type]

    inds_babaabab = [options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type,
                     options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type]

    inds_babaabba = [options.beta_type,  options.alpha_type, options.beta_type, options.alpha_type,
                     options.alpha_type, options.beta_type,  options.beta_type, options.alpha_type]

    inds_bababaab = [options.beta_type, options.alpha_type, options.beta_type,  options.alpha_type,
                     options.beta_type, options.alpha_type, options.alpha_type, options.beta_type]

    inds_babababa = [options.beta_type, options.alpha_type, options.beta_type, options.alpha_type,
                     options.beta_type, options.alpha_type, options.beta_type, options.alpha_type]

    inds_bababbaa = [options.beta_type, options.alpha_type, options.beta_type,  options.alpha_type,
                     options.beta_type, options.beta_type,  options.alpha_type, options.alpha_type]

    inds_abbaaabb = [options.alpha_type, options.beta_type,  options.beta_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.beta_type, options.beta_type]

    inds_abbaabab = [options.alpha_type, options.beta_type, options.beta_type,  options.alpha_type,
                     options.alpha_type, options.beta_type, options.alpha_type, options.beta_type]

    inds_abbaabba = [options.alpha_type, options.beta_type, options.beta_type, options.alpha_type,
                     options.alpha_type, options.beta_type, options.beta_type, options.alpha_type]

    inds_abbabaab = [options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type,
                     options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type]

    inds_abbababa = [options.alpha_type, options.beta_type,  options.beta_type, options.alpha_type,
                     options.beta_type,  options.alpha_type, options.beta_type, options.alpha_type]

    inds_abbabbaa = [options.alpha_type, options.beta_type, options.beta_type,  options.alpha_type,
                     options.beta_type,  options.beta_type, options.alpha_type, options.alpha_type]

    inds_baabaabb = [options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.alpha_type, options.beta_type,  options.beta_type]

    inds_baababab = [options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type]

    inds_baababba = [options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.beta_type,  options.beta_type,  options.alpha_type]

    inds_baabbaab = [options.beta_type, options.alpha_type, options.alpha_type, options.beta_type,
                     options.beta_type, options.alpha_type, options.alpha_type, options.beta_type]

    inds_baabbaba = [options.beta_type, options.alpha_type, options.alpha_type, options.beta_type,
                     options.beta_type, options.alpha_type, options.beta_type,  options.alpha_type]

    inds_baabbbaa = [options.beta_type, options.alpha_type, options.alpha_type, options.beta_type,
                     options.beta_type, options.beta_type,  options.alpha_type, options.alpha_type]

    inds_ababaabb = [options.alpha_type, options.beta_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.alpha_type, options.beta_type, options.beta_type]

    inds_abababab = [options.alpha_type, options.beta_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.beta_type, options.alpha_type, options.beta_type]

    inds_abababba = [options.alpha_type, options.beta_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.beta_type, options.beta_type,  options.alpha_type]

    inds_ababbaab = [options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type,
                     options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type]

    inds_ababbaba = [options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type,
                     options.beta_type,  options.alpha_type, options.beta_type,  options.alpha_type]

    inds_ababbbaa = [options.alpha_type, options.beta_type, options.alpha_type, options.beta_type,
                     options.beta_type,  options.beta_type, options.alpha_type, options.alpha_type]

    inds_aabbaabb = [options.alpha_type, options.alpha_type, options.beta_type, options.beta_type,
                     options.alpha_type, options.alpha_type, options.beta_type, options.beta_type]

    inds_aabbabab = [options.alpha_type, options.alpha_type, options.beta_type,  options.beta_type,
                     options.alpha_type, options.beta_type,  options.alpha_type, options.beta_type]

    inds_aabbabba = [options.alpha_type, options.alpha_type, options.beta_type, options.beta_type,
                     options.alpha_type, options.beta_type,  options.beta_type, options.alpha_type]

    inds_aabbbaab = [options.alpha_type, options.alpha_type, options.beta_type,  options.beta_type,
                     options.beta_type,  options.alpha_type, options.alpha_type, options.beta_type]

    inds_aabbbaba = [options.alpha_type, options.alpha_type, options.beta_type, options.beta_type,
                     options.beta_type,  options.alpha_type, options.beta_type, options.alpha_type]

    inds_aabbbbaa = [options.alpha_type, options.alpha_type, options.beta_type,  options.beta_type,
                     options.beta_type,  options.beta_type,  options.alpha_type, options.alpha_type]

    inds_aaabbaaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type,
                     options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type]

    inds_aaababaa = [options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type]

    inds_aaabaaba = [options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type]

    inds_aaabaaab = [options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type,
                     options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type]

    inds_aababaaa = [options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type,
                     options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type]

    inds_aabaabaa = [options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type,
                     options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type]

    inds_aabaaaba = [options.alpha_type, options.alpha_type, options.beta_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.beta_type, options.alpha_type]

    inds_aabaaaab = [options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type,
                     options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type]

    inds_abaabaaa = [options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type,
                     options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type]

    inds_abaaabaa = [options.alpha_type, options.beta_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.beta_type, options.alpha_type, options.alpha_type]

    inds_abaaaaba = [options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type]

    inds_abaaaaab = [options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type]

    inds_baaabaaa = [options.beta_type, options.alpha_type, options.alpha_type, options.alpha_type,
                     options.beta_type, options.alpha_type, options.alpha_type, options.alpha_type]

    inds_baaaabaa = [options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.beta_type,  options.alpha_type, options.alpha_type]

    inds_baaaaaba = [options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.beta_type,  options.alpha_type]

    inds_baaaaaab = [options.beta_type,  options.alpha_type, options.alpha_type, options.alpha_type,
                     options.alpha_type, options.alpha_type, options.alpha_type, options.beta_type]

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

            if options.verbose:
                print("\n<<< {:}".format(term_rdm4_si))

            for ten_rdm4 in tens_rdm4:
                ten_rdm4_spin_inds = [get_spin_index_type(ind) for ind in ten_rdm4.indices]

                ten_rdm4_tens_sa = []
                const_rdm4_tens_sa = []

                if ten_rdm4_spin_inds in [inds_aaaaaaaa, inds_bbbbbbbb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aaabbaaa, inds_bbbaabbb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 7.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aaababaa, inds_bbbababb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aaabaaba, inds_bbbabbab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aaabaaab, inds_bbbabbba]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aababaaa, inds_bbababbb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabaabaa, inds_bbabbabb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabaaaba, inds_bbabbbab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabaaaab, inds_bbabbbba]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_abaabaaa, inds_babbabbb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_abaaabaa, inds_babbbabb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 7.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_abaaaaba, inds_babbbbab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_abaaaaab, inds_babbbbba]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baaabaaa, inds_abbbabbb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baaaabaa, inds_abbbbabb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baaaaaba, inds_abbbbbab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baaaaaab, inds_abbbbbba]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 7.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabbbbaa, inds_bbaaaabb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabbbaba, inds_bbaaabab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabbbaab, inds_bbaaabba]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabbabba, inds_bbaabaab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabbabab, inds_bbaababa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_aabbaabb, inds_bbaabbaa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 7.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_ababbbaa, inds_babaaabb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_ababbaba, inds_babaabab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_ababbaab, inds_babaabba]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 7.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_abababba, inds_bababaab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_abababab, inds_babababa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_ababaabb, inds_bababbaa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baabbbaa, inds_abbaaabb]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baabbaba, inds_abbaabab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baabbaab, inds_abbaabba]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baababba, inds_abbabaab]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baababab, inds_abbababa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = - 7.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 60.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                elif ten_rdm4_spin_inds in [inds_baabaabb, inds_abbabbaa]:
                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,w,t)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 4, 7]]
                    const_rdm4_sa = 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,w,t,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 4, 7, 6]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 4, 6]]
                    const_rdm4_sa = 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,v,w,u)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 5, 4, 6]]
                    const_rdm4_sa = 1.0 / 5.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,w,u,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 4, 6, 7, 5]]
                    const_rdm4_sa = - 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,w,t,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 4, 7, 5]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,w,u,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 4, 6, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 4, 5]]
                    const_rdm4_sa = 1.0 / 30.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,w,v)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 4, 5]]
                    const_rdm4_sa = 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,u,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 6, 7, 4]]
                    const_rdm4_sa = - 2.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,v,t,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 5, 7, 4]]
                    const_rdm4_sa = - 3.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,v,t,u,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 5, 7, 6, 4]]
                    const_rdm4_sa = - 1.0 / 15.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,u,t,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 6, 7, 5, 4]]
                    const_rdm4_sa = - 1.0 / 10.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                    ## Spin-Adapted RDM term: rdm(p,q,r,s,t,u,v,w)
                    ten_rdm4_sa = ten_rdm4.copy()
                    ten_rdm4_sa.indices = [ten_rdm4_sa.indices[i] for i in [0, 1, 2, 3, 7, 6, 5, 4]]
                    const_rdm4_sa = - 1.0 / 20.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                else:
                    ten_rdm4_sa = ten_rdm4.copy()
                    const_rdm4_sa = 0.0

                    ten_rdm4_tens_sa.append(ten_rdm4_sa)
                    const_rdm4_tens_sa.append(const_rdm4_sa)

                tens_rdm4_sa.append(ten_rdm4_tens_sa)
                consts_rdm4_sa.append(const_rdm4_tens_sa)

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
                    term_rdm4_sa.tensors[ten_rdm4_sa_ind] = ten_rdm4_sa
                    term_rdm4_sa.tensors[ten_rdm4_sa_ind].symmetries = []

                if options.verbose:
                    print("--> {:} (factor = {:.5f})".format(term_rdm4_sa, consts_rdm4_sa_prod[tens_rdm4_sa_ind]))

                terms_rdm4_sa.append(term_rdm4_sa)

        else:
            terms_rdm4_sa.append(term_rdm4_si)

    termChop(terms_rdm4_sa)

    print("Done!")
    return terms_rdm4_sa

def convert_t_amplitudes_si_to_sa(_terms_t_si):
    "Convert T Amplitudes Objects from Spin-Integrated to Spin-Adapted"

    options.print_divider()
    print("Converting T amplitudes to spin-adapted formulation...")

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

                if options.verbose:
                    print("\n<<< {:}".format(term_t1_si))

                if ten_t1_spin_inds in [inds_aa, inds_bb]:
                    ten_t1 = ten.copy()
                    consts_t1_sa_prod = 1.0
                    term_t1_sa.tensors[ten_ind] = ten_t1
                else:
                    ten_t1 = ten.copy()
                    consts_t1_sa_prod = 0.0
                    term_t1_sa.scale(consts_t1_sa_prod)
                    term_t1_sa.tensors[ten_ind] = ten_t1

                if options.verbose:
                    print("--> {:} (factor = {:.5f})".format(term_t1_sa, consts_t1_sa_prod))

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

            if options.verbose:
                print("\n<<< {:}".format(term_t2_si))

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

                    if options.cvs_approach:
                        ten_t2_cvs_inds = [is_cvs_index_type(ind) for ind in ten_t2.indices]

                    if options.cvs_approach and (ten_t2_cvs_inds[0] == ten_t2_cvs_inds[1]):
                        if ten_t2_spin_inds in [inds_aaaa, inds_bbbb]:
                            ## Spin-Adapted 2e- term: t(p,q,r,s)
                            ten1_t2 = ten_t2.copy()
                            const1_t2 = 1.0

                            ## Spin-Adapted 2e- term: t(q,p,r,s)
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
                            ## Spin-Adapted 2e- term: t(q,p,r,s)
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

                    else:
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

                if options.verbose:
                    print("--> {:} (factor = {:.5f})".format(term_t2_sa, consts_t2_sa_prod[tens_t2_sa_ind]))

                terms_t2_sa.append(term_t2_sa)

        else:
            terms_t2_sa.append(term_t2_si)

    termChop(terms_t2_sa)

    print("Done!")
    return terms_t2_sa

def remove_si_tensors_symmetries(_terms_si):
    "Remove Symmetries of non-RDM Spin-Integrated Tensors"

    for _term_ind, _term_si in enumerate(_terms_si):
        for _tensor_ind, _tensor_si in enumerate(_term_si.tensors):

            if isinstance(_tensor_si, kroneckerDelta):
                _terms_si[_term_ind].tensors[_tensor_ind].symmetries = []

            elif _tensor_si.name == 'e' and len(_tensor_si.indices) == 1:
                _terms_si[_term_ind].tensors[_tensor_ind].symmetries = []

            elif _tensor_si.name == 'h' and len(_tensor_si.indices) == 2:
                _terms_si[_term_ind].tensors[_tensor_ind].symmetries = []

            elif _tensor_si.name == 'v' and len(_tensor_si.indices) == 4:
                _terms_si[_term_ind].tensors[_tensor_ind].symmetries = []

            elif _tensor_si.name in ['t1', 't2'] and len(_tensor_si.indices) in [2, 4]:
                _terms_si[_term_ind].tensors[_tensor_ind].symmetries = []

def update_sa_tensors_symmetries(_terms_sa):
    "Update Symmetries of Spin-Adapted Tensors"

    # Define Spin-Adapted Symmetries
    kdelta_sa_symm = [symmetry((1,0), 1)]

    h1e_sa_symm = [symmetry((1,0), 1)]
    v2e_sa_symm = [symmetry((1,0,3,2), 1), symmetry((2,3,0,1), 1)]

    t1_sa_symm = [symmetry((1,0), 1)]
    t2_sa_symm = [symmetry((1,0,3,2), 1), symmetry((2,3,0,1), 1)]

    rdm1_sa_symm = [symmetry((1,0), 1)]
    rdm2_sa_symm = [symmetry((1,0,3,2), 1), symmetry((3,2,1,0), 1)]
    rdm3_sa_symm = [symmetry((1,0,2,3,5,4), 1), symmetry((0,2,1,4,3,5), 1), symmetry((5,4,3,2,1,0), 1)]
    rdm4_sa_symm = [symmetry((1,0,2,3,4,5,7,6), 1), symmetry((0,2,1,3,4,6,5,7), 1), symmetry((0,1,3,2,5,4,6,7), 1),
                    symmetry((7,6,5,4,3,2,1,0), 1)]

    for _term_ind, _term_sa in enumerate(_terms_sa):
        for _tensor_ind, _tensor_sa in enumerate(_term_sa.tensors):
            _terms_sa[_term_ind].tensors[_tensor_ind] = remove_spin_index_type(_tensor_sa)

            if isinstance(_tensor_sa, kroneckerDelta):
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = kdelta_sa_symm

            elif _tensor_sa.name == 'h' and len(_tensor_sa.indices) == 2:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = h1e_sa_symm

            elif _tensor_sa.name == 'v' and len(_tensor_sa.indices) == 4:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = v2e_sa_symm

            elif _tensor_sa.name in ['t1', 't2'] and len(_tensor_sa.indices) == 2:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = t1_sa_symm

            elif _tensor_sa.name in ['t1', 't2'] and len(_tensor_sa.indices) == 4:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = t2_sa_symm

            elif isinstance(_tensor_sa, creDesTensor) and len(_tensor_sa.indices) == 2:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = rdm1_sa_symm

            elif isinstance(_tensor_sa, creDesTensor) and len(_tensor_sa.indices) == 4:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = rdm2_sa_symm

            elif isinstance(_tensor_sa, creDesTensor) and len(_tensor_sa.indices) == 6:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = rdm3_sa_symm

            elif isinstance(_tensor_sa, creDesTensor) and len(_tensor_sa.indices) == 8:
                _terms_sa[_term_ind].tensors[_tensor_ind].symmetries = rdm4_sa_symm
