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
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.sco_id_list)
    assert os.path.exists(args.n_zero_tsv)
    assert os.path.exists(args.cds_in_dir)
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
        transcript_ids: (dict) Dictionary of per-taxon transcript IDs.
            transcript_ids = { taxon_1 : [(transcript_1, ortholog_id_1), (transcript_2, ortholog_id_2)],
                               taxon_2 : [(transcript_1, ortholog_id_1), (transcript_2, ortholog_id_2)]}
    '''
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
    print(f'    Retained {kept:,} hierarchical orthogroups from the N0 table.', flush=True)
    return transcript_ids

def parse_taxon_transcripts(taxon_id:str, cds_fa_path:str, transcript_ids:set)->dict:
    '''
    Parse the transcript sequence from input FASTA and obtain the subset
    of target sequences of interest.
    Args:
        taxon_id: (str) ID of target taxon.
        cds_fa_path: (str) Path to the CDS FASTA of target taxon.
        transcript_ids: (set) Set of target sequences of interest.
    Returns:
        taxon_transcripts: (dict) Pairs of transcript IDs and sequences.
            taxon_transcripts = { transcript_1_id : transcript_1_seq,
                                  transcript_2_id : transcript_2_seq }
    '''
    taxon_transcripts = dict()
    records = 0
    kept = 0

    # Parse input fasta
    with open(cds_fa_path) as fh:
        seq_id = None
        seq = list()
        for line in fh:
            line = line.strip('\n')
            # Ignore empty lines
            if len(line) == 0:
                continue
            # Ignore comment lines
            if line[0] in {'#', '.'}:
                continue
            # Check if its ID or sequence
            if line.startswith('>'):
                records += 1
                if seq_id is not None:
                    sequence = ''.join(seq)
                    if seq_id in transcript_ids:
                        taxon_transcripts[seq_id] = sequence
                        kept += 1
                seq_id = line.lstrip('>').split()[0]
                seq = list()
            else:
                # if it is sequence
                seq.append(line.upper())
    sequence = ''.join(seq)
    if seq_id in transcript_ids:
        taxon_transcripts[seq_id] = sequence
        kept += 1

    print(f'    Parsed input CDS for {taxon_id}: Read {records:,} records, retained {kept:,} sequences.', flush=True)
    return taxon_transcripts

def load_taxon_cds(transcript_ids:dict, cds_dir:str)->dict:
    '''
    Load the CDSs sequences for each taxon.
    Args:
        transcript_ids: (dict) Dictionary of per-taxon transcript IDs.
        cds_dir: (str) Path to the directory containing the per-sample CDS FASTA.
    Returns:
        taxon_cds: (dict) Dictionary of per-taxon coding sequences.
            taxon_cds = { taxon_1 : { transcript_1_id : transcript_1_seq,
                                      transcript_2_id : transcript_2_seq },
                          taxon_2 : { transcript_1_id : transcript_1_seq,
                                      transcript_2_id : transcript_2_seq } }
    '''
    print('\nLoading CDSs...', flush=True)
    taxon_cds = dict()
    for taxon in transcript_ids:
        cds_fa_path = f'{cds_dir}/{taxon}.fa'
        if not os.path.exists(cds_fa_path):
            sys.exit(f'Error: {cds_fa_path} not found.')
        # Extract all the transcript of that sample into a single set for ease of search
        taxon_transcript_ids = set([ transcript[0] for transcript in transcript_ids[taxon] ])
        taxon_transcripts = parse_taxon_transcripts(taxon, cds_fa_path, taxon_transcript_ids)
        taxon_cds[taxon] = taxon_transcripts
    return taxon_cds

def regroup_transcripts(transcript_ids:dict)->dict:
    '''
    Regroup the per-taxon transcript IDs into per-ortogroup transcript IDs.
    Args:
        transcript_ids: (dict) Dictionary of per-taxon transcript IDs.
    Returns:
        orthogroup_transcripts: (dict) Dictionary of per-orthogroup transcript IDs.
            orthogroup_transcripts = { hog_id : [ (transcript_id, taxon_id),
                                                  (transcript_id, taxon_id) ] }
    '''
    orthogroup_transcripts = dict()
    for taxon in transcript_ids:
        for transcript in transcript_ids[taxon]:
            trans_id = transcript[0]
            hog_id = transcript[1]
            # Initialize the output directory
            orthogroup_transcripts.setdefault(hog_id, list())
            orthogroup_transcripts[hog_id].append((trans_id, taxon))
    return orthogroup_transcripts

def group_sco_sequences(transcript_ids:dict, taxon_cds:dict,
                        out_dir:str, fa_line_width:int=FA_LINE_WIDTH)->None:
    '''
    Group the taxon CDS sequences for each single-copy ortholog, and save as 
    individual FASTAs.
    Args:
        transcript_ids: (dict) Dictionary of per-taxon transcript IDs.
        taxon_cds: (dict) Dictionary of per-taxon coding sequences.
        out_dir: (str) Path to output directory to save the FASTA.
        fa_line_width: (int) Line width for output FASTA.
    Returns:
        None
    '''
    print('\nGrouping sequences and processing outputs...')
    # First, you have to re-organize the transcript IDs to be per orthogroup.
    orthogroup_transcripts = regroup_transcripts(transcript_ids)
    # Loop across each orthogroup and extract its sequences.
    for hog_id in orthogroup_transcripts:
        # Prepare a new output FASTA for that orthogroup
        out_fa = f'{out_dir}/{hog_id}.cds.fa'
        with open(out_fa, 'w') as fa_fh:
            # Now, process each taxon
            for taxon_transcript in orthogroup_transcripts[hog_id]:
                trans_id = taxon_transcript[0]
                taxon = taxon_transcript[1]
                # Get the sequences for that transcript
                sequence = taxon_cds[taxon].get(trans_id, None)
                if sequence is None:
                    sys.exit(f'Error: {trans_id} transcript not found for {taxon}.')
                # Prepare the new FASTA record
                fa_fh.write(f'>{trans_id}_{taxon}\n')
                # Wrap the sequence lines up to `fa_line_width` characters
                for start in range(0, len(sequence), fa_line_width):
                    seq_line = sequence[start:(start+fa_line_width)]
                    fa_fh.write(f'{seq_line}\n')

def main():
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    # Loading target IDs
    sco_ids = load_sco_ids(args.sco_id_list)
    # Parse the N0 table
    transcript_ids = parse_n0_table(args.n_zero_tsv, sco_ids, args.out_dir)
    # Load the per-taxon CDSs
    taxon_cds = load_taxon_cds(transcript_ids, args.cds_in_dir)
    # Prepare the outputs
    group_sco_sequences(transcript_ids, taxon_cds, args.out_dir)


    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
