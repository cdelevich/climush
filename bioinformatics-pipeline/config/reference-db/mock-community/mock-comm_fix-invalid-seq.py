# import the mock community sequences from the mock community .fasta file from the climush repository
with open(mockcomm_seqs_pathin, 'r') as seqfile_in:
    mockcomm_seqs = SeqIO.to_dict(SeqIO.parse(seqfile_in, 'fasta'))

# error-causing seqrecord info (value)
seqrecord_error = []
for rec in mockcomm_seqs.values():
    try:
        print(rec, end='\n\n')
    except UnicodeDecodeError:
        seqrecord_error.append(rec)

err_match_seqs = {}
for seqrec_err in seqrecord_error:
    seqrec_err_id = seqrec_err.id
    seqrec_err_regex = re.compile(seqrec_err_id)
    with open(mockcomm_seqs_pathin, 'r') as seqfile_in:
        seqrecord_lines = seqfile_in.readlines()
        for l,line in enumerate(seqrecord_lines):
            if line.startswith('>'):
                err_match = seqrec_err_regex.search(line)
                if err_match:
                    seqrec_err_seq = ''
                    for s,seq in enumerate(seqrecord_lines[(l+1):]):
                        if seq.startswith('>'):
                            break
                        else:
                            seqrec_err_seq += seq.strip()
                    err_match_seqs.update({seqrec_err_id:seqrec_err_seq})
                else:
                    continue
            else:
                continue

err_match_seqs_fixed = {}
for mockcomm_id, mockcomm_seq in err_match_seqs.items():

    mockcomm_seq_fixed = keep_only_nucleotides(
        nt_str=mockcomm_seq,
        accepted_nt='any',
    )

    # add the sequence record ID and the corrected sequence to the fixed dictionary
    err_match_seqs_fixed.update({mockcomm_id:mockcomm_seq_fixed})
