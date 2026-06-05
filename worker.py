def process_chunk(terms_chunk):
    for t in terms_chunk:
        t.makeCanonical(rename_user_defined = False)
    return terms_chunk

