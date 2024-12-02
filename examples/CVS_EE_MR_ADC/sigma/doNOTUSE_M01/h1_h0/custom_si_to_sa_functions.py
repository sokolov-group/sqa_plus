import sqa_plus

def convert_X_si_to_sa(_terms_x_si, options):
    print("Converting Sigma vector to spin-adapted formulation ...")

    from sqa_plus.sqaSpinAdapted import remove_spin_index_type
    from sqa_plus.sqaIndex import get_spatial_index_type, get_spin_index_type

    from sqa_plus.sqaTerm import termChop

    import itertools

    sqa_plus.options.print_divider()
    print("Converting custom vector to spin-adapted formulation...")

    # Define 1e- indices lists
    inds_aa = [sqa_plus.options.alpha_type, sqa_plus.options.alpha_type]
    inds_bb = [sqa_plus.options.beta_type,  sqa_plus.options.beta_type]

    # Convert One-Body Amplitudes
    terms_x_sa = []
    for term_x_si in _terms_x_si:
        term_x_sa = term_x_si.copy()

        for ten_ind, ten in enumerate(term_x_sa.tensors):
            if ten.name[0] == 'X':
                ten_x_spin_inds = [get_spin_index_type(ind) for ind in ten.indices]

                print("\n<<< {:}".format(term_x_si))

                if ten_x_spin_inds in [inds_aa, inds_bb]:
                    ten_x = ten.copy()
                    consts_x_sa_prod = 1.0
                    term_x_sa.tensors[ten_ind] = ten_x
                else:
                    ten_x = ten.copy()
                    consts_x_sa_prod = 0.0
                    term_x_sa.scale(consts_x_sa_prod)
                    term_x_sa.tensors[ten_ind] = ten_x

                print("--> {:} (factor = {:.5f})".format(term_x_sa, consts_x_sa_prod))

        terms_x_sa.append(term_x_sa)

    termChop(terms_x_sa)
    print("Done!")
    return terms_x_sa

#def convert_X_si_to_sa(_terms_x_si, options):
#    options.print_divider()
#    options.print_header("Converting Sigma vector to spin-adapted formulation")
#
#    from sqa_plus.sqaSpinAdapted import remove_spin_index_type
#
#    # Define Spin-Adapted 2e- integrals Symmetries
#    x_sa_symm = []
#
#    # Convert objects in each term
#    for term_ind, term_x_si in enumerate(_terms_x_si):
#        for ten_ind, ten in enumerate(term_x_si.tensors):
#            if ten.name[0] == 'X':
#                _terms_x_si[term_ind].tensors[ten_ind] = remove_spin_index_type(ten)
#                _terms_x_si[term_ind].tensors[ten_ind].symmetries = x_sa_symm
#
#    options.print_divider()
#
#    return _terms_x_si
