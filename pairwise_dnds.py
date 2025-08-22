#!/usr/bin/env python3
from datetime import datetime
from Bio import SeqIO
import sys
import os
import argparse
import math
import dnds

# Some constants
PROG = sys.argv[0].split('/')[-1]
MIN_ALN_LEN = 25
ALN_PREFIX = 'fa'

def parse_args():
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-s', '--sco-list',
                   required=True,
                   help='(str) Path to the single-copy orthgroup table \
                    (produced by `extract_orthogroups_cds.py`).')
    p.add_argument('-a', '--alignments',
                   required=True,
                   help='(str) Path to the directory with the trimmed \
                    multiple sequence alignments.')
    p.add_argument('-i', '--ingroup',
                   required=True,
                   help='(str) ID of the ingroup (focal) species in the alignment. \
                    Used to report gene/transcript IDs in the output.')
    p.add_argument('-o', '--out-dir',
                   required=False,
                   default='.',
                   help='(str) Path to output directory [default=./].')
    p.add_argument('-m', '--min-aln-len',
                   required=False,
                   type=int,
                   default=MIN_ALN_LEN,
                   help=f'(int) Minimum length required to keep an \
                    alignment [default={MIN_ALN_LEN}]')
    p.add_argument('-p', '--aln-prefix',
                   required=False,
                   default=ALN_PREFIX,
                   help=f'(str) Prefix for the alignment FASTA files \
                    [default={ALN_PREFIX}].')

    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.sco_list)
    assert os.path.exists(args.alignments)
    args.alignments = args.alignments.rstrip('/')
    args.out_dir = args.out_dir.rstrip('/')
    assert os.path.exists(args.out_dir)
    # Set some constants
    return args

def date() -> str:
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.now().strftime("%Y-%m-%d")

def time() -> str:
    '''Print the current time in HH:MM:SS format.'''
    return datetime.now().strftime("%H:%M:%S")

def parse_sco_list(sco_f:str)->dict:
    '''
    Parse the Single Copy Orthogroup list file
    Args
        sco_f (str): Path to the Single Copy Orthogroup list file.
    Returns:
        sco_dict (dict): Dictionary with orthogroup IDs as keys and a 
                         pair of taxon-gene IDs as values.
        sco_dict = { orthogroup_id: { taxon_id : gene_id, 
                                      taxon_id : gene_id },
                     orthogroup_id: { taxon_id : gene_id, 
                                      taxon_id : gene_id }, }
    '''
    print('\nParsing Single Copy Orthogroup input table...', flush=True)
    sco_dict = {}
    records = 0
    taxa = []
    with open(sco_f, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip('\n')
            if line.startswith('#') or len(line)==0:
                continue
            records += 1
            fields = line.split('\t')
            sco_id = fields[0]
            taxon = fields[1]
            gene_id = fields[2]
            # Add the orthogroup ID to the dictionary if not already present
            sco_dict.setdefault(sco_id, {})
            sco_dict[sco_id][taxon] = gene_id
            # Add the taxon to the taxa list if not already present
            if taxon not in taxa:
                taxa.append(taxon)
    # Report to log
    print(f'    Parsed {records:,} records.')
    print(f'    Found {len(taxa):,} taxa: {", ".join(taxa)}.')
    print(f'    And {len(sco_dict):,} Single Copy Orthogroups.', flush=True)
    return sco_dict

def extract_sequences(in_msa_f:str, ingroup:str)->tuple[str, str]:
    '''
    Extract the sequences from an input MSA FASTA. The
    alignment must contain a max of two sequences.
    Args:
        in_msa_f (str): Path to input MSA FASTA
        ingroup (str): ID of ingroup focal species.
    Returns:
        (ingroup_seq, outgroup_seq): tuple of sequences for 
            the ingroup and outgroup taxa.
    '''
    ingroup_seq = ''
    outgroup_seq = ''
    ingroup_found = False
    count = 0
    # Parse the input fasta
    for record in SeqIO.parse(in_msa_f, 'fasta'):
        count += 1
        # Process if ingroup
        if record.id == ingroup:
            ingroup_found = True
            ingroup_seq = record.seq
        else:
            outgroup_seq = record.seq
    # The fasta must contain only two sequences.
    if count != 2:
        sys.exit(f'Error: MSA FASTA {in_msa_f} does not have two sequences.\
Alignment must only be between a pair of ingroup and outgroup sequences.')
    # The FASTA must contain the ingroup
    if not ingroup_found:
        sys.exit(f'Error: Ingroup sequence ID ({ingroup}) not found in input MSA FASTA:\
{in_msa_f}')
    # And the sequences must be of the same length
    if len(ingroup_seq) != len(outgroup_seq):
        sys.exit(f'Error: The sequences in {in_msa_f} are of unequal length.')
    return ingroup_seq, outgroup_seq

def calculate_dnds(seq1:str, seq2:str)->list:
    '''
    Calculate synomynous vs non-synomymous rates between a
    pair of sequences.
    Args:
        seq1 (str): Sequence, for ingroup
        seq2 (str): Sequence, for outgroup
    Returns:
        dnds_out (list): List of the dN/dS output.

    Based on the code by S. Small
    '''
    # Tally the number of synynymous vs non-synonymous
    # sites
    syn_sites = dnds.syn_sum(seq1, seq2)
    non_sites = len(seq1) - syn_sites
    # Tally the number of synynymous vs non-synonymous
    # substitutions
    syn_subs, non_subs = dnds.substitutions(seq1, seq2)
    # Calculate the proportions to substitutions vs sites
    pN = non_subs / non_sites
    pS = syn_subs / syn_sites
    # Calculate using the Jukes-Cantor correction
    dN = -(3 / 4) * math.log(1 - (4 * pN / 3))
    dS = -(3 / 4) * math.log(1 - (4 * pS / 3))
    # Calculate the dN/dS ratio
    dnds_rat = 0
    if dS>0:
        dnds_rat = dN/dS
    # Prepare to return
    # These two are type Fraction and not Float, and
    # can cause some formatting issues downstream
    non_sites = float(non_sites)
    syn_sites = float(syn_sites)
    # The rest leave as is.
    dnds_out = [non_sites, # Number of non-syn sites
                syn_sites, # Number of syn sites
                non_subs,  # Number of non-syn subs
                syn_subs,  # Number of syn subs
                pN,        # Proportion of non-syn subs vs sites
                pS,        # Proportion of syn subs vs sites
                dN,        # Corrected proportion of non-syn subs vs sites
                dS,        # Corrected proportion of syn subs vs sites
                dnds_rat]  # dN/dS ratio
    return dnds_out

def process_orthogroups(sco_list:dict, ingroup:str, alignments:str,
                        outdir:str='.', min_aln_len:int=MIN_ALN_LEN,
                        aln_prefix:str=ALN_PREFIX)->None:
    '''
    Process the alignments for all orthogroups. Extract sequences,
    calculate dN/dS, and report.
    Args:
        sco_dict (dict): Dictionary with orthogroup IDs as keys and a 
                pair of taxon-gene IDs as values.
        ingroup (str): ID of ingroup focal species.
        alignments (str): Path to the alignment files.
        out_dir (str): Path to output directory [default=.].
        min_aln_len (int): Minimum length required to keep an 
                           alignment [default=MIN_ALN_LEN].
        aln_prefix (str): Prefix for all MSA fasta files.
    '''
    print('\nProcessing orthogroups...', flush=True)
    print(f'    Looking for multiple-sequence alignment FASTAs in the following format\n\
        {alignments}/OGNNNNNNN.{aln_prefix}', flush=True)
    # Prepare outputs
    out_tsv = f'{outdir}/pairwide_dNdS_{ingroup}.tsv'
    with open(out_tsv, 'w', encoding='utf-8') as fh:
        header = ['Orthogroup',
                  'GeneID',
                  'SeqLen',
                  'NonSynSites',
                  'SynSites',
                  'NonSynSubs',
                  'SynSubs',
                  'pN',
                  'pS',
                  'dN',
                  'dS',
                  'dNdS']
        header = '\t'.join(header)
        fh.write(f'{header}\n')

        # Loop over all orthogroups
        n_found = 0
        n_processed = 0
        for sco_id in sco_list:
            # Select the genes for the focal taxon
            sco_genes = sco_list[sco_id]
            gene = sco_genes.get(ingroup, None)
            if gene is None:
                sys.exit(f'Error: {ingroup} not found among the taxa for orthogroup {sco_id}')
            # Select the input multiple sequence alignment and
            # extract the input sequences.
            in_msa_f = f'{alignments}/{sco_id}.{aln_prefix}'
            # Some of these MSAs will not exists. This is expected.
            # Just skip them.
            if not os.path.exists(in_msa_f):
                continue
            n_found += 1
            ingroup_seq, outgroup_seq = extract_sequences(in_msa_f, ingroup)
            # Skip sequences if they are too small
            if len(ingroup_seq) < min_aln_len:
                continue
            # Now calculate dN/dS
            try:
                dnds_out = calculate_dnds(ingroup_seq, outgroup_seq)
            except Exception:
                continue
            # Save to output
            # Writing this in long form for record keeping purposes
            non_sites = dnds_out[0]
            syn_sites = dnds_out[1]
            non_subs  = dnds_out[2]
            syn_subs  = dnds_out[3]
            pN        = dnds_out[4]
            pS        = dnds_out[5]
            dN        = dnds_out[6]
            dS        = dnds_out[7]
            dnds_rat  = dnds_out[8]
            row = [f'{sco_id}',
                   f'{gene}',
                   f'{len(ingroup_seq)}',
                   f'{non_sites}',
                   f'{syn_sites}',
                   f'{non_subs}',
                   f'{syn_subs}',
                   f'{pN}',
                   f'{pS}',
                   f'{dN}',
                   f'{dS}',
                   f'{dnds_rat}']
            row = '\t'.join(row)
            fh.write(f'{row}\n')
            n_processed += 1

    # Report to log.
    print(f'\nFound alignments for a total of {n_found:,} orthogroups.', flush=True)
    print(f'    Calculated dN/dS for {n_processed} pairs of sequences.', flush=True)

def main():
    '''
    Main function
    '''
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    print(f'\nReporting dN/dS values for {args.ingroup}')

    # Parse the Single Copy Orthologues list
    sco_list = parse_sco_list(args.sco_list)

    # Process all the orthogroups
    process_orthogroups(sco_list, args.ingroup, args.alignments,
                        args.out_dir, args.min_aln_len,
                        args.aln_prefix)

    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
