import sqa_plus

def convert_X_si_to_sa(_terms_x_si):
    print("Converting Sigma vector to spin-adapted formulation ...")

    from sqa_plus.sqaSpinAdapted import remove_spin_index_type
    from sqa_plus.sqaIndex import get_spatial_index_type, get_spin_index_type

    import itertools

    sqa_plus.options.print_divider()
    print("Converting custom vector to spin-adapted formulation...")

###
#    # Define 1e- indices lists
#    inds_aa = [options.alpha_type, options.alpha_type]
#    inds_bb = [options.beta_type,  options.beta_type]
#
#    # Convert One-Body Amplitudes
#    terms_t1_sa = []
#    for term_t1_si in _terms_t_si:
#        term_t1_sa = term_t1_si.copy()
#
#        for ten_ind, ten in enumerate(term_t1_sa.tensors):
#            if ten.name[0] == 't' and len(ten.indices) == 2:
#                ten_t1_spin_inds = [get_spin_index_type(ind) for ind in ten.indices]
#
#                if options.verbose:
#                    print("\n<<< {:}".format(term_t1_si))
#
#                if ten_t1_spin_inds in [inds_aa, inds_bb]:
#                    ten_t1 = ten.copy()
#                    consts_t1_sa_prod = 1.0
#                    term_t1_sa.tensors[ten_ind] = ten_t1
#                else:
#                    ten_t1 = ten.copy()
#                    consts_t1_sa_prod = 0.0
#                    term_t1_sa.scale(consts_t1_sa_prod)
#                    term_t1_sa.tensors[ten_ind] = ten_t1
#
#                if options.verbose:
#                    print("--> {:} (factor = {:.5f})".format(term_t1_sa, consts_t1_sa_prod))
#
#        terms_t1_sa.append(term_t1_sa)
#
#    termChop(terms_t1_sa)
###

    # Define 2e- indices lists
    inds_aaaa = [sqa_plus.options.alpha_type, sqa_plus.options.alpha_type, sqa_plus.options.alpha_type, sqa_plus.options.alpha_type]
    inds_bbbb = [sqa_plus.options.beta_type,  sqa_plus.options.beta_type,  sqa_plus.options.beta_type,  sqa_plus.options.beta_type]

    inds_abab = [sqa_plus.options.alpha_type, sqa_plus.options.beta_type,  sqa_plus.options.alpha_type, sqa_plus.options.beta_type]
    inds_baba = [sqa_plus.options.beta_type,  sqa_plus.options.alpha_type, sqa_plus.options.beta_type,  sqa_plus.options.alpha_type]
    inds_abba = [sqa_plus.options.alpha_type, sqa_plus.options.beta_type,  sqa_plus.options.beta_type,  sqa_plus.options.alpha_type]
    inds_baab = [sqa_plus.options.beta_type,  sqa_plus.options.alpha_type, sqa_plus.options.alpha_type, sqa_plus.options.beta_type]

    # Convert objects in each term
    terms_x_sa = []
    for term_x_si in _terms_x_si:
        # List for storing X objects
        tens_x = []
        tens_x_ind = []

        for ten_ind, ten in enumerate(term_x_si.tensors):
            if ten.name[0] == 'X' and len(ten.indices) == 4:
                tens_x.append(ten)
                tens_x_ind.append(ten_ind)

        ## Convert object
        if tens_x:
            tens_x_sa = []
            consts_x_sa = []

            print("\n<<< {:}".format(term_x_si))

            for ten_x in tens_x:

                ten_x_inds = [get_spatial_index_type(ind) for ind in ten_x.indices]
                ten_x_spin_inds = [get_spin_index_type(ind) for ind in ten_x.indices]

                ten_x_tens_sa = []
                const_x_tens_sa = []

                if (ten_x_inds[0] == ten_x_inds[1]) and (ten_x_inds[2] == ten_x_inds[3]):

                    if ten_x_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted term: x(p,q,r,s)
                        ten_x_sa = ten_x.copy()
                        const_x_sa = 1.0
    
                        ten_x_tens_sa.append(ten_x_sa)
                        const_x_tens_sa.append(const_x_sa)
                        
                        ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                        ten_x_sa = ten_x.copy()
                        ten_x_sa.indices = [ten_x_sa.indices[i] for i in [0, 1, 3, 2]]
                        const_x_sa = -1.0
    
                        ten_x_tens_sa.append(ten_x_sa)
                        const_x_tens_sa.append(const_x_sa)
    
                    elif ten_x_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: v2e(p,q,r,s)
                        ten_x_sa = ten_x.copy()
                        const_x_sa = 1.0
    
                        ten_x_tens_sa.append(ten_x_sa)
                        const_x_tens_sa.append(const_x_sa)
    
                    elif ten_x_spin_inds in [inds_abba, inds_baab]:
                        ## Spin-Adapted 2e- term: v2e(p,q,s,r)
                        ten_x_sa = ten_x.copy()
                        ten_x_sa.indices = [ten_x_sa.indices[i] for i in [0, 1, 3, 2]]
                        const_x_sa = - 1.0
    
                        ten_x_tens_sa.append(ten_x_sa)
                        const_x_tens_sa.append(const_x_sa)
    
                    else:
                        ten_x_sa = ten_x.copy()
                        const_x_sa = 0.0
    
                        ten_x_tens_sa.append(ten_x_sa)
                        const_x_tens_sa.append(const_x_sa)

                elif (ten_x_inds[0] != ten_x_inds[1]) and (ten_x_inds[2] == ten_x_inds[3]):

                    if ten_x_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten_x = ten_x.copy()
                        const_x = 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                        ## Spin-Adapted 2e- term: t(p,q,s,r)
                        ten_x = ten_x.copy()
                        ten_x.indices = [ten_x.indices[i] for i in [0, 1, 3, 2]]
                        const_x = - 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                    elif ten_x_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten_x = ten_x.copy()
                        const_x = 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                    elif ten_x_spin_inds in [inds_abba, inds_baab]:
                        ten_x = ten_x.copy()
                        ten_x.indices = [ten_x.indices[i] for i in [0, 1, 3, 2]]
                        const_x = - 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                    else:
                        ten_x = ten_x.copy()
                        const_x = 0.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                elif (ten_x_inds[0] == ten_x_inds[1]) and (ten_x_inds[2] != ten_x_inds[3]):

                    if ten_x_spin_inds in [inds_aaaa, inds_bbbb]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten_x = ten_x.copy()
                        const_x = 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                        ## Spin-Adapted 2e- term: t(p,q,s,r)
                        ten_x = ten_x.copy()
                        ten_x.indices = [ten_x.indices[i] for i in [1, 0, 2, 3]]
                        const_x = - 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                    elif ten_x_spin_inds in [inds_abab, inds_baba]:
                        ## Spin-Adapted 2e- term: t(p,q,r,s)
                        ten_x = ten_x.copy()
                        const_x = 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                    elif ten_x_spin_inds in [inds_abba, inds_baab]:
                        ten_x = ten_x.copy()
                        ten_x.indices = [ten_x.indices[i] for i in [1, 0, 2, 3]]
                        const_x = - 1.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                    else:
                        ten_x = ten_x.copy()
                        const_x = 0.0

                        ten_x_tens_sa.append(ten_x)
                        const_x_tens_sa.append(const_x)

                elif ((ten_x_inds[0] != ten_x_inds[1]) and (ten_x_inds[2] != ten_x_inds[3])):

                    if options.cvs_approach:
                        ten_x_cvs_inds = [is_cvs_index_type(ind) for ind in ten_x.indices]

                    if options.cvs_approach and (ten_x_cvs_inds[0] == ten_x_cvs_inds[1]):
                        if ten_x_spin_inds in [inds_aaaa, inds_bbbb]:
                            ## Spin-Adapted 2e- term: t(p,q,r,s)
                            ten_x = ten_x.copy()
                            const_x = 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                            ## Spin-Adapted 2e- term: t(q,p,r,s)
                            ten_x = ten_x.copy()
                            ten_x.indices = [ten_x.indices[i] for i in [1, 0, 2, 3]]
                            const_x = - 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                        elif ten_x_spin_inds in [inds_abab, inds_baba]:
                            ## Spin-Adapted 2e- term: t(p,q,r,s)
                            ten_x = ten_x.copy()
                            const_x = 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                        elif ten_x_spin_inds in [inds_abba, inds_baab]:
                            ## Spin-Adapted 2e- term: t(q,p,r,s)
                            ten_x = ten_x.copy()
                            ten_x.indices = [ten_x.indices[i] for i in [1, 0, 2, 3]]
                            const_x = - 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                        else:
                            ten_x = ten_x.copy()
                            const_x = 0.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                    else:
                        if ten_x_spin_inds in [inds_aaaa, inds_bbbb]:
                            ## Spin-Adapted 2e- term: t(p,q,r,s)
                            ten_x = ten_x.copy()
                            const_x = 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                            ## Spin-Adapted 2e- term: t(p,q,s,r)
                            ten_x = ten_x.copy()
                            ten_x.indices = [ten_x.indices[i] for i in [0, 1, 3, 2]]
                            const_x = - 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                        elif ten_x_spin_inds in [inds_abab, inds_baba]:
                            ## Spin-Adapted 2e- term: t(p,q,r,s)
                            ten_x = ten_x.copy()
                            const_x = 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                        elif ten_x_spin_inds in [inds_abba, inds_baab]:
                            ten_x = ten_x.copy()
                            ten_x.indices = [ten_x.indices[i] for i in [0, 1, 3, 2]]
                            const_x = - 1.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)

                        else:
                            ten_x = ten_x.copy()
                            const_x = 0.0

                            ten_x_tens_sa.append(ten_x)
                            const_x_tens_sa.append(const_x)


            tens_x_sa.append(ten_x_tens_sa)
            const_x_sa.append(const_x_tens_sa)


            tens_x_sa_permut = []
            for item in list(itertools.product(*tens_x_sa)):
                tens_x_sa_permut.append(list(item))

            consts_x_sa.append(const_x_tens_sa)

            consts_x_sa_permut = []
            for item in list(itertools.product(*consts_x_sa)):
                consts_x_sa_permut.append(list(item))

            consts_x_sa_prod = []
            for iter in consts_x_sa_permut:
                prod = 1.0
                for const in iter:
                    prod = prod * const
                consts_x_sa_prod.append(prod)

            for tens_x_sa_ind, tens_x_sa in enumerate(tens_x_sa_permut):
                term_x_sa = term_x_si.copy()
                term_x_sa.scale(consts_x_sa_prod[tens_x_sa_ind])

                for ten_x_sa_ind, ten_x_sa in zip(tens_x_ind, tens_x_sa):
                    term_x_sa.tensors[ten_x_sa_ind] = ten_x_sa
                    term_x_sa.tensors[ten_x_sa_ind].symmetries = [sqa_plus.symmetry((1,0,3,2), 1)]

                print("--> {:} (factor = {:.5f})".format(term_x_sa, consts_x_sa_prod[tens_x_sa_ind]))

                terms_x_sa.append(term_x_sa)

        else:
            terms_x_sa.append(term_x_si)

    print("Done!")
    return terms_x_sa

