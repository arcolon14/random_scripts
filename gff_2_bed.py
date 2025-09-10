#!/usr/bin/env python3
import sys
import os
import argparse
from datetime import datetime

# Some constants
PROG = sys.argv[0].split('/')[-1]
MIN_ALN_LEN = 25
ALN_PREFIX = 'fa'

def parse_args():
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser()
    p.add_argument('-g', '--gff',
                   required=True,
                   help='(str) Path to input GFF.')
    p.add_argument('-f', '--feature',
                   required=False,
                   default=None,
                   help='(str) GFF feature to be extracted.')
    p.add_argument('-o', '--out-dir',
                   required=False,
                   default='.',
                   help='(str) Path to output directory [default=./].')

    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.gff)
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

def parse_gff(gff:str, feature:str|None=None, outdir:str='.')->None:
    '''
    Parse a GFF and write the selected feature elements
    as a BED file.
    '''
    print(f'\nParsing GFF file: {gff}', flush=True)
    if feature is not None:
        print(f'    Extracting feature: {feature}', flush=True)
    records = 0
    kept = 0
    seen_features = set()
    # Prepare output
    name = gff.split('/')[-1]
    name = name.rstrip('.gz').rstrip('.gff').rstrip('.gff3')
    bed = f'{outdir}/{name}.{feature}.bed'
    if feature is None:
        bed = f'{outdir}/{name}.ALL.bed'
    with open(bed, 'w', encoding='utf-8') as fh:
        # Parse the GFF
        for line in open(gff, encoding='utf-8'):
            line = line.strip('\n')
            if line.startswith('#') or len(line)==0:
                continue
            records += 1
            fields = line.split('\t')
            primary_tag = fields[2]
            seen_features.add(primary_tag)
            # If feature is selected, filter.
            if primary_tag != feature and feature is not None:
                continue
            # Select the components for the target elements
            chrom = fields[0]
            start = int(fields[3])-1
            end = int(fields[4])
            attributes = fields[8]
            fid = None
            parent = None
            for attribute in attributes.split(';'):
                if attribute.startswith('ID='):
                    fid = attribute.split('=')[1]
                    fid = fid.lstrip(' \'\"').rstrip(' \'\"')
                elif attribute.startswith('Parent='):
                    parent = attribute.split('=')[1]
                    parent = parent.lstrip(' \'\"').rstrip(' \'\"')
            if parent is None:
                parent = fid
            kept += 1
            # Prepare output
            row = f'{chrom}\t{start}\t{end}\t{primary_tag}\t{fid}\t{parent}'
            fh.write(f'{row}\n')
        # # Check if any elements were recovered.
        if kept < 1:
            seen = ", ".join(list(seen_features))
            print(f'\nError: {kept} records of the feature {feature} found.')
            print(f'    The input GFF contains the following features: {seen}', flush=True)
            sys.exit()

        # Print to log
        print(f'\nParsed {records:,} total records from input GFF.')
        print(f'    Kept {kept:,} records from the target feature(s).')

def main():
    '''
    Main function
    '''
    print(f'{PROG} started on {date()} {time()}.')
    # Parse args
    args = parse_args()
    # Parse GFF
    parse_gff(args.gff, args.feature, args.out_dir)

    # Done!
    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
