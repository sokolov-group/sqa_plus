def convert_X_si_to_sa(_terms_x_si):
    print(">>> converting Sigma vector to spin-adapted formulation ...")

    x_sa_symm = []

    for term_ind, term_x_si in enumerate(_terms_x_si):
        for ten_ind, ten in enumerate(term_x_si.tensors):
            if ten.name[0] == 'X' and len(ten.indices) == 2:
                _terms_x_si[term_ind].tensors[ten_ind].symmetries = x_sa_symm
    return _terms_x_si
