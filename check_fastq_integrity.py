#!/usr/bin/env python3
import sys, os, gzip

# Read fastq file and keep read ids on a list
def read_fq_read_ids(fastq_f):
    fh = open(fastq_f, 'r')
    if fastq_f.endswith('.gz'):
        fh = gzip.open(fastq_f, 'rt')
    print(f'Processing FASTQ file:\n    {fastq_f}', flush=True)
    records = 0
    headers = list()
    for i, line in enumerate(fh):
        line = line.strip('\n')
        if i%4==0:
            if line[0] != '@':
                sys.exit(f'Error: missing FASTQ header {i} {line}')
            records += 1
            line = line.split(' ')[0]
            headers.append(line)
    print(f'    Read {i+1:,} lines in file.', flush=True)
    print(f'    Composed of {records:,} FASTQ records.', flush=True)
    return headers

# Compare the list of FASTQ headers
def compare_fq_headers(for_headers, rev_headers):
    if len(for_headers) == 0:
        sys.exit('No R1 headers found')
    elif len(rev_headers) == 0:
        sys.exit('No R2 headers found')
    elif len(for_headers) != len(rev_headers):
        sys.exit(f'Lengths of R1 ({len(for_headers):,}) and R2 ({len(rev_headers):,}) do not match.')
    else:
        matching = 0
        for i in range(len(for_headers)):
            r1 = for_headers[i]
            r2 = rev_headers[i]
            if r1 != r2:
                sys.exit('Records do not match:\n    Record {i} in R1: {r1}\n    Record {i} in R2: {r2}')
            else:
                matching += 1
    print(f'Matched {matching:,} FASTQ records in R1 and R2', flush=True)

# CMD line Opts
R1 = sys.argv[1]
R2 = sys.argv[2]

# Check the input files files
if not os.path.exists(R1):
    sys.exit(f'Error: {R1} not found')
if not os.path.exists(R1):
    sys.exit(f'Error: {R2} not found')

# Run
for_headers = read_fq_read_ids(sys.argv[1])
rev_headers = read_fq_read_ids(sys.argv[2])
compare_fq_headers(for_headers, rev_headers)
