#!/usr/bin/env python3
import sys
import os
import argparse
import json
from datetime import datetime

# Some constants
PROG = sys.argv[0].split('/')[-1]
VALID_ANALYSES = ['aBSREL', 'BUSTED']
ALPHA = 0.05

def parse_args():
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-s', '--sco-list',
                   required=True,
                   help='(str) Path to orthofinder Orthogroups/Orthogroups_SingleCopyOrthologues.txt file.')
    p.add_argument('-y', '--hyphy-outs',
                   required=True,
                   help='(str) Path to the HyPhy output files directory.')
    p.add_argument('-o', '--out-dir',
                   required=False,
                   default='.',
                   help='(str) Path to output directory [default=./].')
    p.add_argument('-a', '--analysis-type',
                   choices=VALID_ANALYSES,
                   nargs='+',
                   default=['aBSREL'],
                   help='(str) Type of HyPhy analysis to parse [choices: aBSREL, BUSTED] [default=aBSREL].')
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.sco_list)
    assert os.path.exists(args.hyphy_outs)
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
    with open(sco_f, 'r') as fh:
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

def parse_busted_outputs(sco_list:dict, hyphy_outs:str, out_dir:str)->None:
    '''
    Parse BUSTED output files.
    Args:
        sco_list (dict): Dictionary with orthogroup IDs as keys and a pair 
                         of taxon-gene IDs as values.
        hyphy_outs (str): Path to the HyPhy output files directory.
        out_dir (str): Path to output directory.
    '''
    print('\nParsing BUSTED outputs...', flush=True)
    n_found = 0
    n_signif = 0
    # Generate output file
    out_file = f'{out_dir}/parsed_BUSTED.tsv'
    with open(out_file, 'w') as fh:
        # Write header
        # TODO: Add needed columns
        fh.write('Orthogroup\tTaxon\tGeneID\tAlnLen\tSignifBranch\tCorrPVal\tLRT\n')

        # Loop through each orthogroup in the sco_list
        for sco_id in sco_list:
            sco_genes = sco_list[sco_id]
            # Select the target BUSTED output file
            busted_file = f'{hyphy_outs}/{sco_id}.BUSTED.json'
            # The BUSTED outputs are not going to be present for some
            # orthogroups. This is expected.
            if not os.path.exists(busted_file):
                continue
            # Open the BUSTED output file and extract the results per taxa.
            busted_model_outs, signif = read_busted_file(busted_file, sco_id, sco_genes)
            if signif:
                n_signif += 1
            n_found += 1
            # Save to the ouput file
            for row in busted_model_outs:
                fh.write(row)
    print(f'    Found and extracted BUSTED results for {n_found:,} orthogroups.', flush=True)
    print(f'    A total of {n_signif:,} showed evidence of selection on any branch based on BUSTED\'s model.', flush=True)

def read_busted_file(busted_file:str, sco_id:str, sco_genes:dict, alpha=ALPHA)->list:
    '''
    Read and parse BUSTED output json file.
    Args:
        busted_file (str): Path to the BUSTED output file.
        scp_id (str): ID for the single-copy orthogroup
        sco_genes (dict): Dictionary with taxon-gene IDs for the orthogroup.
        alpha (float): p-value threshold for significance [default=0.05]
    Returns:
        model_outs (list): list of the parsed BUSTED outputs per taxa
        signif (bool): Is there significant evidence of selection in *any* branch?
    '''
    model_outs = []
    signif = False
    with open(busted_file, 'r') as fh:
        busted_data = json.load(fh)

        # Check the general test results to see if significant
        test_results = busted_data['test results']
        print(test_results)
        # TODO: Continue parsing the BUSTED file.

    return model_outs, signif


def read_absrel_file(absrel_file:str, sco_id:str, sco_genes:dict)->list:
    '''
    Read and parse aBSREL output json file.
    Args:
        absrel_file (str): Path to the aBSREL output file.
        scp_id (str): ID for the single-copy orthogroup
        sco_genes (dict): Dictionary with taxon-gene IDs for the orthogroup.
    Returns:
        model_outs (list): list of the parsed aBSREL outputs per taxa
        signif (bool): Is there significant evidence of selection in *any* branch?
    '''
    model_outs = []
    signif = False
    with open(absrel_file, 'r') as fh:
        absrel_data = json.load(fh)

        # Check the general test results to see if significant
        # e.g., {'P-value threshold': 0.05, 'positive test results': 0, 'tested': 5}
        test_results = absrel_data['test results']
        n_signif = test_results['positive test results']
        p_val_cutoff = test_results['P-value threshold']
        # Is *any* branch significant?
        if n_signif > 0:
            signif = True

        # Check the inputs of the model to check for alignment length
        # e.g., {'file name': 'input_msa.fa', 'number of sequences': 5,
        #        'number of sites': 588, 'partition count': 1,
        #        'trees': {'0': '((A:0.0394154,B:0.0554306)Node3:0.0567855,
        #                         (C:0.103642,D:0.0897399)Node6:0.0363832,E:0.106333)'}}
        inputs = absrel_data['input']
        aln_len = inputs['number of sites']

        # Check the per-taxon branch attributes. These have the per-taxon results
        # of the model.
        taxon_branch_attributes = absrel_data['branch attributes']['0']
        # Loop over each taxon
        for taxon in taxon_branch_attributes:
            # First, choose the branches corresponding to the focal taxa
            if taxon not in sco_genes:
                continue
            # Then, select the gene for target taxon
            gene = sco_genes.get(taxon, None)
            if gene is None:
                sys.exit(f'Error: Taxon ({taxon}) mismatch in orthogroup {sco_id}.')
            signif_branch = 0
            taxon_attrs = taxon_branch_attributes[taxon]
            # e.g., {'Baseline MG94xREV': 0.4675796811578989,
            #        'Baseline MG94xREV omega ratio': 0.05156506382256223,
            #        'Corrected P-value': 1,
            #        'Full adaptive model': 0.43577628844732,
            #        'Full adaptive model (non-synonymous subs/site)': 0.0840848856637228,
            #        'Full adaptive model (synonymous subs/site)': 0.3516914027835968,
            #        'LRT': 0,
            #        'Nucleotide GTR': 0.2170857999678362,
            #        'Rate Distributions': [[0.05012255288081785, 1]],
            #        'Rate classes': 1,
            #        'Uncorrected P-value': 1,
            #        'original name': 'Name',
            #        'rate at which 2 nucleotides are changed instantly within a single codon': 0.9084505089461798,
            #        'rate at which 3 nucleotides are changed instantly within a single codon': 2.546331370694687}
            corr_p_val = taxon_attrs['Corrected P-value']
            lrt = taxon_attrs['LRT']
            nonsyn_site = taxon_attrs['Full adaptive model (non-synonymous subs/site)']
            syn_site = taxon_attrs['Full adaptive model (synonymous subs/site)']
            omega_ratio = taxon_attrs['Baseline MG94xREV omega ratio']
            n_rate_classes = taxon_attrs['Rate classes']
            rate_distributions = taxon_attrs['Rate Distributions']
            # Is gene significant?
            if corr_p_val < p_val_cutoff:
                signif_branch = 1
            # Format the rate distribution for the output row
            rate_str = []
            for rate in rate_distributions:
                rate_str.append(f'{rate[0]}:{rate[1]}')
            rate_str = ';'.join(rate_str)
            # Construct the output row
            row =  f'{sco_id}\t{taxon}\t{gene}\t{aln_len}\t{signif_branch}\t{corr_p_val}\t{lrt}\t'
            row += f'{omega_ratio}\t{nonsyn_site}\t{syn_site}\t{n_rate_classes}\t{rate_str}\n'
            model_outs.append(row)
    return model_outs, signif

def parse_absrel_outputs(sco_list:dict, hyphy_outs:str, out_dir:str)->None:
    '''
    Parse aBSREL output files.
    Args:
        sco_list (dict): Dictionary with orthogroup IDs as keys and a pair 
                         of taxon-gene IDs as values.
        hyphy_outs (str): Path to the HyPhy output files directory.
        out_dir (str): Path to output directory.
    '''
    print('\nParsing aBSREL outputs...', flush=True)
    n_found = 0
    n_signif = 0
    # Generate output file
    out_file = f'{out_dir}/parsed_aBSREL.tsv'
    with open(out_file, 'w') as fh:
        # Write header
        fh.write('Orthogroup\tTaxon\tGeneID\tAlnLen\tSignifBranch\tCorrPVal\tLRT\tOmegaRatio\tNonSynSite\tSynSite\tRateClasses\tRateDists\n')

        # Loop through each orthogroup in the sco_list
        for sco_id in sco_list:
            sco_genes = sco_list[sco_id]
            # Select the target aBSREL output file
            absrel_file = f'{hyphy_outs}/{sco_id}.aBSREL.json'
            # The absrel outputs are not going to be present for some
            # orthogroups. This is expected.
            if not os.path.exists(absrel_file):
                continue
            # Open the aBSREL output file and extract the results per taxa.
            absrel_model_outs, signif = read_absrel_file(absrel_file, sco_id, sco_genes)
            if signif:
                n_signif += 1
            n_found += 1
            # Save to the ouput file
            for row in absrel_model_outs:
                fh.write(row)
    print(f'    Found and extracted aBSREL results for {n_found:,} orthogroups.', flush=True)
    print(f'    A total of {n_signif:,} showed evidence of selection on any branch based on aBSREL\'s model.', flush=True)

def parse_hyphy_outs(sco_list:dict, hyphy_outs:str, out_dir:str, analysis_type:list)->None:
    '''
    Parse HyPhy output files for the specified analysis type.
    Args:
        sco_list (dict): Dictionary with orthogroup IDs as keys and a 
                         pair of taxon-gene IDs as values.
        hyphy_outs (str): Path to the HyPhy output files directory.
        out_dir (str): Path to output directory.
        analysis_type (list): List of HyPhy analysis types to parse.
    '''
    # Loop over the different analysis types and process accordingly
    for analysis in analysis_type:
        if analysis == 'aBSREL':
            # Parse aBSREL output files
            parse_absrel_outputs(sco_list, hyphy_outs, out_dir)
        elif analysis == 'BUSTED':
            # Parse BUSTED output files
            # TODO: Implement BUSTED output parsing
            parse_busted_outputs(sco_list, hyphy_outs, out_dir)
            continue

def main():
    '''
    Main function
    '''
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    # Report the analysis type outputs parsed.
    print('\nParsing HyPhy output files for analysis type(s):', flush=True)
    for a in args.analysis_type:
        print(f'    {a}', flush=True)

    # Parse the Single Copy Orthologues list
    sco_list = parse_sco_list(args.sco_list)

    # Parse HyPhy output files
    parse_hyphy_outs(sco_list, args.hyphy_outs, args.out_dir, args.analysis_type)

    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
