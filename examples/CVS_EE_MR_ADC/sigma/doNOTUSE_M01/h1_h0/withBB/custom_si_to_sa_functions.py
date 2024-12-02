import sqa_plus

def convert_X_si_to_sa(_terms_x_si):
    sqa_plus.options.print_divider()

    from sqa_plus.sqaSpinAdapted import remove_spin_index_type

    # Define Spin-Adapted 2e- integrals Symmetries
    x_sa_symm = []

    # Convert objects in each term
    for term_ind, term_x_si in enumerate(_terms_x_si):
        for ten_ind, ten in enumerate(term_x_si.tensors):
            if ten.name[0] == 'X':
                _terms_x_si[term_ind].tensors[ten_ind] = remove_spin_index_type(ten)
                _terms_x_si[term_ind].tensors[ten_ind].symmetries = x_sa_symm

    sqa_plus.options.print_divider()

    return _terms_x_si
