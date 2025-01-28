cdef bint match(
    const unsigned char* string, 
    const unsigned char* pattern, 
    const int allowed_mismatches
):
    cdef int n_mismatch = 0
    cdef int i = 0
    cdef int pattern_len = len(pattern)
    cdef int string_len = len(string)

    if string_len < pattern_len:
        return False

    for i in range(pattern_len):
        if not string[i] == pattern[i]:
            n_mismatch += 1
        
        if n_mismatch > allowed_mismatches:
            return False

    return True


cdef find_match(
    const unsigned char* seq, 
    dict bc_dict
):
    cdef char* match_name = b''
    cdef const unsigned char* bc
    cdef dict bcinfo
    cdef bint m
    for bc, bcinfo in bc_dict.items():
        m = match(seq, bc, bcinfo['mismatch'])
        if m:
            match_name = bcinfo['name']
            break
        
    return match_name


def match_with_errors(
    const unsigned char* seq, 
    dict bc_dict, 
    const int laxity
):
    match = None
    for match_start in range(laxity):
        match = find_match(seq[match_start:], bc_dict)
        if match:
            break

    return match, match_start
