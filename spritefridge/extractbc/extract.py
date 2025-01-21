import logging


def find_match(seq, bc_dict):
    match_name = b''
    for _, bcinfo in bc_dict.items():
        re = bcinfo['regex']
        m = re.match(seq)
        if m:
            match_name = bcinfo['name']
            break
        
    return match_name


def regex_match(seq, bc_dict, laxity):
    match = None
    for match_start in range(laxity):
        match = find_match(seq[match_start:], bc_dict)
        if match:
            break

    return match, match_start


def hash_match(seq, bc_dict, min_len, max_len):
    # if we only deal with one length we skip the loop
    if min_len == max_len:
        match = bc_dict.get(seq)
        return match['name'] if match else b'', max_len

    match_name = b''
    for bc_len in range(min_len, max_len + 1):
        match = bc_dict.get(seq[:bc_len])
        if match:
            match_name = match['name']
            break

    return match_name, bc_len


def extract_barcodes(read, bc_dicts, layout, laxity = 6):
    start = 0
    read_bcs = []
    readseq = memoryview(read['seq'])
    # print(readseq)
    for bc_cat, min_bc_len, max_bc_len, allowed_mismatches in layout:
        # this is a shortcut to avoid matching the full SPACER cat
        if bc_cat.startswith('S'):
            # print(bc_cat, start, start + max_bc_len)
            # print(' '* start + readseq[start: start + max_bc_len])
            start += max_bc_len
            continue

        if not allowed_mismatches:
            bc_match, match_len = hash_match(
                readseq[start: start + max_bc_len],
                bc_dicts[bc_cat],
                min_bc_len,
                max_bc_len
            )
            # print(bc_cat, start, start + match_len)
            # print(' '* start + readseq[start: start + match_len])
            read_bcs.append(bc_match)
            start += match_len
            continue
        
        bc_match, match_pos = regex_match(
            readseq[start: start + max_bc_len + laxity],
            bc_dicts[bc_cat],
            laxity
        )
        # print(bc_cat, start, start + match_pos + max_bc_len)
        # print(' '* (start + match_pos)  + readseq[start + match_pos : start + match_pos + max_bc_len])
        read_bcs.append(bc_match)
        start += (match_pos + max_bc_len)

    return read_bcs


def extract_parallel(
    bc_dicts,
    layout1,
    layout2,
    laxity,
    input_queue, 
    output_queue
):
    logging.info('start extractor process')
    while True:
        reads = input_queue.get()
        if not reads:
            break
        
        for readpair in reads:
            read1, read2 = readpair
            readpair.append(
                extract_barcodes(read1, bc_dicts, layout1, laxity) +
                extract_barcodes(read2, bc_dicts, layout2, laxity)
            )

        output_queue.put(reads)

    # termination signal
    logging.info('received empty readlist, shutting down extractor')
    output_queue.put([])
