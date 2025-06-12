#!/usr/bin/env python3
import sys, os, argparse
from datetime import datetime

# Some constants
PROG = sys.argv[0].split('/')[-1]
FA_LINE_WIDTH = 60

def parse_args(prog=PROG):
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-s', '--sco-id-list', required=True, 
                   help='(str) Path to file containing the list of single-copy orthogroup IDs.')
    p.add_argument('-n', '--n-zero-tsv', required=True, 
                   help='(str) Path to the orthofinder N0 hierarchical orthogroup file.')
    p.add_argument('-c', '--cds-in-dir', required=True, 
                   help='(str) Path to the directory containing the input per-taxon CDS sequences.')
    p.add_argument('-o', '--out-dir', required=False, default='.',
                   help='(str) Path to output directory [default=./].')
    p.add_argument('-t', '--threads', required=False, type=int, default=1,
                   help='(int) Number of threads to run in parallel sections of code [default=1].')
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.sco_id_list)
    assert os.path.exists(args.n_zero_tsv)
    assert os.path.exists(args.cds_in_dir)
    args.out_dir = args.out_dir.rstrip('/')
    assert os.path.exists(args.out_dir)
    assert args.threads >= 1
    # Set some constants
    return args

def date() -> str:
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.now().strftime("%Y-%m-%d")

def time() -> str:
    '''Print the current time in HH:MM:SS format.'''
    return datetime.now().strftime("%H:%M:%S")

def load_sco_ids(sco_ids_f:str)->list:
    '''
    Load the SCO IDs from the input file.
    Args:
        sco_ids_f: (str) Path to the SCO IDs file
    Returns:
        sco_ids: (list) List of SCO IDs
    '''
    print('\nLoading the Single-Copy Orthogroup IDs...', flush=True)
    sco_ids = list()
    with open(sco_ids_f) as fh:
        for line in fh:
            line = line.strip('\n')
            if line.startswith('#') or len(line)==0:
                continue
            if not line.startswith('N0'):
                sys.exit(f'Error: {line} is invalid N0 orthogroup ID. They must start with N0.NOG[...].')
            sco_ids.append(line)
    print(f'    Loaded {len(sco_ids):,} SCO IDs from input file.')
    return sco_ids

def parse_n0_table(n0_tsv_f:str, sco_ids:list, out_dir:str)->dict:
    '''
    Parse the orthofinder N0 table and extract the per-taxon transcript IDs
    for each Single-Copy Ortholog.
    Args:
        n0_tsv_f: (str) Path to the N0.tsv file
        sco_ids: (list) List of single-copy ortholog IDs
        out_dir: (str) Path to output directory, for reports.
    Returns:
        transcript_ids: (dict) Dictionary of per-taxon transcript IDs.'''
    print('\nParsing the orthoginder N0 table...')
    records = 0
    kept = 0
    taxa = list()
    scos = set(sco_ids) # To make checking faster
    # Main output
    transcript_ids = dict()
    with open(n0_tsv_f) as fh:
        for line in fh:
            line = line.strip('\n')
            if line.startswith('#') or len(line)==0:
                continue
            # Parse the specific header and extract the taxon IDs
            if line.startswith('HOG\tOG'):
                fields = line.split('\t')
                for taxon in fields[3:]:
                    taxa.append(taxon)
                print(f'    Parsing data for {len(taxa):,} input taxa:')
                for t in taxa:
                    print(f'        {t}')
                continue
            # Parse all the other lines
            fields = line.split('\t')
            records += 1
            # Skip the orthogroups that are not single-copy
            hog_id = fields[0]
            node = fields[2]
            # if hog_id not in scos:
            #     continue
            # if not node != 'n0':
            #     continue
            # TODO: Check
            # Loop over the gene IDs per taxon
            keep = True
            for i, spp_genes in enumerate(fields[3:]):
                spp_genes_l = list()
                for gene in spp_genes.split(' '):
                    gene = gene.rstrip(',').rstrip('\'').lstrip('\'')
                    if len(gene) > 0:
                        spp_genes_l.append(gene)
                n_genes = len(spp_genes_l)
                # print(spp_genes_l, n_genes)
                if n_genes != 1:
                    keep = False
                    continue
            # Only process the orthogroups kept
            if keep:
                for i, gene in enumerate(fields[3:]):
                    taxon = taxa[i]
                    # Add to the dictionary
                    transcript_ids.setdefault(taxon, list())
                    transcript_ids[taxon].append((gene, hog_id))
                    # print(hog_id, taxon, gene)
                kept += 1

                # print(fields)

            # if records > 1000:
            #     break

    # print(transcript_ids)
    print(f'    Retained {kept:,} hierarchical orthogroups from the N0 table.', flush=True)

def main():
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    # Loading target IDs
    sco_ids = load_sco_ids(args.sco_id_list)
    # Parse the N0 table
    transcript_ids = parse_n0_table(args.n_zero_tsv, sco_ids, args.out_dir)

    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
