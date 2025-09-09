#!/usr/bin/env python3
import sys, os, argparse, statistics
from datetime import datetime

PROG = sys.argv[0].split('/')[-1]
MIN_LEN = 10_000
MIN_INTERVAL=1
TARGET_FEATURES = [] # TODO: remove this.
DESC = """Determine the proportion of each of the annotated genetic \
elements in a BED across the sites present in an phastCons conserved \
sites BED file. Provide other general stats for the phastCons BED."""

def parse_args():
    '''Set and verify command line options.'''
    p = argparse.ArgumentParser(prog=PROG, description=DESC)
    p.add_argument('-f', '--fai', required=True, 
                   help='(str) Path to genome index in FAI format.')
    p.add_argument('-a', '--annotation', required=True, 
                   help='(str) Path to the annotation in BED format.')
    p.add_argument('-p', '--phastcons', required=True,
                   help='(str) Path to the phastCons conserved sited BED.')
    p.add_argument('-o', '--out-dir', required=False, default='.',
                   help='(str) Path to output directory [default=.].')
    # TODO: Min len arg
    # Check inputs
    args = p.parse_args()
    assert os.path.exists(args.fai)
    assert os.path.exists(args.annotation)
    assert os.path.exists(args.phastcons)
    args.out_dir = args.out_dir.rstrip('/')
    return args

def date() -> str:
    '''Print the current date in YYYY-MM-DD format.'''
    return datetime.now().strftime("%Y-%m-%d")

def time() -> str:
    '''Print the current time in HH:MM:SS format.'''
    return datetime.now().strftime("%H:%M:%S")

def load_fai(fai_f:int, min_len:int=MIN_LEN)->dict:
    f'''
    Load the chromosome sizes from a genome FAI index.
    Args:
        fai_f: (str) Path to input FAI.
        min_len: (int/float) Minimum length to load a sequence [default={MIN_LEN:,}].
    Returns:
        chromosomes: (dict) chromosome ids/length pairs.
    '''
    print(f'\nLoading chromosomes from FAI (retaining sequences larger than {min_len:,} bp).', flush=True)
    chromosomes = dict()
    records = 0
    with open(fai_f) as fh:
        for line in fh:
            line = line.strip('\n')
            if len(line)==0 or line.startswith('#'):
                continue
            records += 1
            fields = line.split('\t')
            chrom_id = fields[0]
            chrom_len = int(fields[1])
            if chrom_len < min_len:
                continue
            chromosomes[chrom_id] = chrom_len
    # report to log
    n_chrs = len(chromosomes)
    print(f'    Read {records:,} records from the input FAI.', flush=True)
    print(f'    Retained {n_chrs:,} ({(n_chrs/records):0.2%}) records as chromosomes objects.', flush=True)
    return chromosomes

def load_phastcons_bed(phastcons_bed_f:str, chromosomes:dict,
                       min_len:int=MIN_INTERVAL)->dict:
    '''
    Load the conserved genomic intervals from a phastCons BED.
    Args:
        phastcons_bed_f: (str) Path to input phastCons conserved BED.
        chromosomes: (dict) chromosome ids/length pairs.
        min_len: (int) Minimum length required to retain an interval.
    Returns:
        phastcons: (dict) Dictionary of per-chromosome conserved intervals.
            { chrom_1_id : [ (start_1, end_1), (start_2, end_2), 
                             ..., (start_n, end_n) ], ... }
    '''
    phastcons = dict()
    print('\nLoading conserved genomic intervals from phastCons BED...')
    seen = 0
    kept = 0
    with open(phastcons_bed_f) as fh:
        for line in fh:
            line = line.strip('\n')
            if len(line)==0 or line.startswith('#'):
                continue
            seen += 1
            fields = line.split('\t')
            chrom = fields[0]
            start = int(fields[1]) # BED start are 0-based inclusive
            end = int(fields[2])   # BED end are 0-based exclusive
            # Check the formatting of the BED
            assert end>start
            # Skip small elements
            if (end-start) < min_len:
                continue
            # Only process elements in the target chromosome
            if chrom not in chromosomes:
                continue
            # Check the coordinates
            if start<0:
                sys.exit(f'Error: start coordinate must be > 0, Line: {line}')
            if end>chromosomes[chrom]:
                sys.exit(f'Error: end coordinate must be smaller than chromosome length ({chromosomes[chrom]}), Line: {line}')
            # Add to the intervals to the output dictionary
            phastcons.setdefault(chrom, [])
            phastcons[chrom].append((start, end))
            kept += 1
    # report to log
    print(f'    Read {seen:,} records from the input BED.', flush=True)
    print(f'    Kept {kept:,} ({(kept/seen):0.2%}) records as genomic intervals.', flush=True)
    return phastcons

def calculate_phastcons_stats(phastcons:dict, chromosomes:dict, outdir:str='.')->None:
    '''
    Calculate the per-chromosome stats of the phastCons conserved sites.
    Args:
        phastcons: (dict) Dictionary of per-chromosome conserved intervals.
        chromosomes: (dict) chromosome ids/length pairs.
        outdir: (str) Path to output directory [default='.']
    Returns:
        None
    '''
    out_f = f'{outdir}/phastCons_stats.tsv'
    print(f'\nCalculating stats for phastCons elements and saving to file:\n    {out_f}', flush=True)
    with open(out_f, 'w') as fh:
        header = ['chromID', 'chromLen', 'phastConsNum', 'phastConsLen',
                  'phastConsProp', 'meanLen', 'medianLen','sdLen', 'minLen', 'maxLen']
        header = '\t'.join(header)
        fh.write(f'{header}\n')
        # Loop over the chromosomes and calculate per-chrom stats
        for chromosome in chromosomes:
            chr_len = chromosomes[chromosome]
            chr_phastcons = phastcons.get(chromosome, [])
            size_dist = [ phastcon[1]-phastcon[0] for phastcon in chr_phastcons ]
            # Calculate stats from the size distribution
            phasc_n  = len(size_dist)
            total_l  = sum(size_dist)
            perc_l   = total_l/chr_len
            # Initialize these to 0 to prevent stat errors
            mean_l   = 0
            median_l = 0
            sd_l     = 0
            min_l    = 0
            max_l    = 0
            # Now, if values are present calculate stats
            if phasc_n>0:
                mean_l   = statistics.mean(size_dist)
                median_l = statistics.median(size_dist)
                min_l    = min(size_dist)
                max_l    = max(size_dist)
            if phasc_n>1:
                sd_l     = statistics.stdev(size_dist)
            # Save to output file
            row = [chromosome,              # chromID
                   f'{chr_len}',            # chromLen
                   f'{phasc_n}',            # numPhastCons
                   f'{total_l}',            # totalLen
                   f'{perc_l:0.6f}',        # totalProp
                   f'{mean_l:0.3f}',        # meanLen
                   f'{median_l:0.3g}',      # medianLen
                   f'{sd_l:0.3f}',          # sdLen
                   f'{min_l}',              # minLen
                   f'{max_l}' ]             # maxLen
            row = '\t'.join(row)
            fh.write(f'{row}\n')

class Annotation:
    '''Store the annotation intervals from a BED.'''
    def __init__(self, chrom:str, feature:str, start:int, end:int):
        self.chr = chrom
        self.fet = feature
        self.sta = start
        self.end = end
    def __str__(self):
        return f'{self.chr}\t{self.fet}\t{self.sta}\t{self.end}'

def load_annotation_bed(annotation_bed_f:str, 
                        chromosomes:dict)->dict:
    '''
    Parse a GFF file and load the annotated features.
    Args:
        annotation_bed_f: (str) Path to input annotation BED.
        chromosomes: (dict) chromosome ids/length pairs.
    Returns:
        annotations: (dict) per-chromosome Annotation objects
            { chr_id : [ Annotation_1, Annotation_2, ..., Annotation_n ], ... }
    '''
    print('\nLoading the annotatated features from the input BED...', flush=True)
    annotations = dict()
    feat_tally = dict()
    seen = 0
    kept = 0
    # Parse the file
    with open(annotation_bed_f) as fh:
        for line in fh:
            line = line.strip('\n')
            if len(line)==0 or line.startswith('#'):
                continue
            seen += 1
            fields = line.split('\t')
            # The input BED needs to have 4 fields
            assert len(fields) == 4, 'Input BED must have four columns: Chr<tab>Start<tab>End<tab>Feature/Annotation'
            chrom = fields[0]
            # Only process elements in the target chromosome
            if chrom not in chromosomes:
                continue
            # Set the coordinates
            start  = int(fields[1]) # BED start are 0-based inclusive
            end    = int(fields[2]) # BED end are 0-based exclusive
            # Check the coordinates
            if start<0:
                sys.exit(f'Error: start coordinate must be > 0, Line: {line}')
            if end>chromosomes[chrom]:
                sys.exit(f'Error: end coordinate must be smaller than chromosome length ({chromosomes[chrom]}), Line: {line}')
            assert start >=0, f'Error: starting coordinate must be '
            # Keep a track of the different features/annotations
            feature = fields[3]
            feat_tally.setdefault(feature, 0)
            feat_tally[feature] += 1
            # Create the Annotation object from the parsed elements
            annotation = Annotation(chrom, feature, start, end)
            # Add the annotation feature to the output dict
            annotations.setdefault(chrom, [])
            annotations[chrom].append(annotation)
            kept += 1
    # report to log
    print(f'    Read {seen:,} records from the input GFF.', flush=True)
    print(f'    Kept {kept:,} ({(kept/seen):0.2%}) records as annotations from the following features:', flush=True)
    for feat in sorted(feat_tally):
        n = feat_tally[feat]
        print(f'        {feat} : {n:,}',
              flush=True)
    return annotations

def tally_phastcons_annotations(phastcons:dict, annotations:dict,
                                chromosomes:dict, outdir:str='.',
                                target_features=TARGET_FEATURES)->None:
    '''
    Calculate the overlap between the phastcons and annotation elements
    and tally across different annotation elements.
    Args:
        phastcons: (dict) per-chromosome conserved intervals.
        annotations: (dict) per-chromosome Annotation objects
        chromosomes: (dict) chromosome ids/length pairs.
        outdir: (str) Path to output directory [default='.']
    Returns:
        None
    '''
    out_f = f'{outdir}/phastCons_annotations.tsv'
    print(f'\nCalculating the overlap between the phastcons and annotation elements and saving to file:\n    {out_f}', flush=True)
    with open(out_f, 'w') as fh:
        header = ['chromID', 'phastConsLen', 'featType',
                  'featLen', 'featProp']
        header = '\t'.join(header)
        fh.write(f'{header}\n')
        # Loop and process the chromosomes
        for i, chromosome in enumerate(chromosomes):
            # First, generate a object to hold the results for that chrom
            results = dict()
            # Get the annotations for that chromosome
            feature_sites = set_feature_sites(chromosome, annotations,
                                              target_features)
            # Process the phastCons for the target chromosome
            chr_phastcons = phastcons.get(chromosome, None)
            if chr_phastcons is None:
                continue
            total_phastcon_len = 0
            # Then, iterare over the phastCon intervals
            for interval in chr_phastcons:
                phastcon_set = set(range(interval[0], interval[1]))
                total_phastcon_len += len(phastcon_set)
                # Compare that phastCon set against the different feature
                # sites. The intersection between the two set is the 
                # overlap between the genomic interval of the two 
                # elements. As you do, keep a tally in the results object.
                for feature in feature_sites:
                    feat_s = feature_sites[feature]
                    overlap = len(phastcon_set.intersection(feat_s))
                    results.setdefault(feature, 0)
                    results[feature] += overlap
            # Determine the "other" category (non-target features)
            other = total_phastcon_len
            for feature in results:
                other -= results[feature]
            results['other'] = other
            # Prepare the output
            for feature in results:
                span = results[feature]
                prop = 0.0
                if total_phastcon_len>0:
                    prop = span/total_phastcon_len
                row = [chromosome,               # chromID
                       f'{total_phastcon_len}',  # phastConsLen
                       feature,                  # featType
                       f'{span}',                # featLen
                       f'{prop:0.8f}' ]         # featProp
                row = '\t'.join(row)
                fh.write(f'{row}\n')

def set_feature_sites(chromosome:int, annotations:dict,
                         target_features:list=TARGET_FEATURES)->dict:
    '''
    Process the annotations for a target chromsome, generating a set
    of site per annotation feature.
    Args:
        chromosome: (int) Target chromosome to process
        annotations: (dict) per-chromosome Annotation objects
        target_features: (list) list of features to retain from the GFF.
    Returns:
        feature_sites: (dict) feature feature sites (positions)
            { feature_1 : [ site1, site2, ..., siten ], ... }
    '''
    # Initialize outputs
    feature_sites = dict()
    for feature in target_features:
        feature_sites.setdefault(feature, set())
    # Get the annotations for the target chromosome
    chr_annotations = annotations.get(chromosome, [])
    if chr_annotations is None:
        return feature_sites
    # Loop over the annotations
    for annotation in chr_annotations:
        assert isinstance(annotation, Annotation)
        feature = annotation.fet
        start = annotation.sta
        end = annotation.end
        # Loop over each site in the annotation and add to 
        # the corresponding set
        for site in range(start, end):
            feature_sites[feature].add(site)
    return feature_sites

def find_overlapping_annotations(annotations:dict,
                                 chromosomes:dict,
                                 outdir:str='.')->dict:
    '''
    Identify and tag overlapping annotation windows.
    Args:
        annotations: (dict) per-chromosome Annotation objects
        chromosomes: (dict) chromosome ids/length pairs.
        outdir (str): Path to output directory
    Returns:
        annotations_no_olap: (dict) per-chromosome annotation
            containing the new "multiple" feature for
            overlapping annotations.
    '''
    out_bed = f'{outdir}/annotations_no_overlap.bed'
    print(f'\nProcessing annotations to identify and collapse overlapping sites. Saving to:\n    {out_bed}')
    annotations_no_olap = dict()

    # TODO: This merges adjacent (but non-overlapping) windows
    # with the same annotation. This doesn't matter for the final
    # result; however it would be good if this did not happen to
    # avoid skewing the length statistics later on. These should be
    # taken from the BED or GFF anyway, but making sure.

    with open(out_bed, 'w') as fh:

        # Loop over chromosome and process
        for chromosome in chromosomes:
            chr_len = chromosomes[chromosome]
            chr_annotations = annotations[chromosome]

            # Initialize the spans of the whole chromosome.
            # By default these will be annotated to "other"
            # and filled as annotations are read
            per_site_feats = { s : 'other' for s in range(chr_len) }

            # Loop over the annotations in the chromosome
            for annotation in chr_annotations:
                assert isinstance(annotation, Annotation)
                feature = annotation.fet
                start = annotation.sta
                end = annotation.end
                # Loop over the range of the annotation
                for s in range(start, end):
                    # "other" is the initialized state.
                    # If a feature has been added to this site already
                    if per_site_feats[s] != 'other':
                        # And if its not the same as the current feature,
                        # this means this site has multiple annotations,
                        # e.g., an intron+transposable element
                        if per_site_feats[s] != feature:
                            per_site_feats[s] = 'multiple'
                        # Otherwise, keep the current feature
                        else:
                            per_site_feats[s] = feature
                    # No feature has been added to the site, so just update
                    else:
                        per_site_feats[s] = feature

            # With info on all sites, collapse them down to new Annotation
            # objects by finding consecutive sites annotated with the same feature.
            curr_start = 0
            curr_feature = None

            for site in sorted(per_site_feats.keys()):
                site_feature = per_site_feats[site]

                # If this is the first position or we have a new feature
                if curr_feature is None or site_feature != curr_feature:
                    # If we had a previous interval, save it
                    if curr_feature is not None:
                        new_annotation = Annotation(chromosome, curr_feature, 
                                                        curr_start, site)
                        annotations_no_olap.setdefault(chromosome, [])
                        annotations_no_olap[chromosome].append(new_annotation)
                        fh.write(f'{new_annotation.chr}\t{new_annotation.sta}\t{new_annotation.end}\t{new_annotation.fet}\n')
                    # Start a new interval
                    curr_start = site
                    curr_feature = site_feature
            # Don't forget the last interval
            if curr_feature is not None:
                new_annotation = Annotation(chromosome, curr_feature, 
                                                curr_start, chr_len)
                annotations_no_olap.setdefault(chromosome, [])
                annotations_no_olap[chromosome].append(new_annotation)
                fh.write(f'{new_annotation.chr}\t{new_annotation.sta}\t{new_annotation.end}\t{new_annotation.fet}\n')

    return annotations_no_olap


def calculate_annotation_stats(annotations:dict, chromosomes:dict,
                               outdir:str='.')->None:
    '''
    Calculate the per-chromosome stats of the annotations BED.
    Args:
        annotations: (dict) per-chromosome Annotation objects
        chromosomes: (dict) chromosome ids/length pairs.
        outdir: (str) Path to output directory [default='.']
    Returns:
        None
    '''
    out_f = f'{outdir}/annotation_stats.tsv'
    print(f'\nCalculating stats for the annotation and saving to file:\n    {out_f}', flush=True)
    with open(out_f, 'w') as fh:
        header = ['chromID', 'chromLen', 'featType', 'featNum', 
                  'featLen', 'featProp', 'meanLen', 'medianLen',
                  'sdLen', 'minLen', 'maxLen']
        header = '\t'.join(header)
        fh.write(f'{header}\n')
        # Loop over the chromosomes and calculate per-chrom stats
        for chromosome in chromosomes:
            chr_len = chromosomes[chromosome]
            # Set an "other" categery for non-annotated regions
            other = chr_len
            chr_annotations = annotations.get(chromosome, [])
            # This will be the output for that chromosome
            feature_tally = dict()
            # Loop over the annotations in the chromosome
            for annotation in chr_annotations:
                assert isinstance(annotation, Annotation)
                feat  = annotation.fet
                start = annotation.sta
                end   = annotation.end
                span  = end-start
                # Add to the tally dict
                feature_tally.setdefault(feat, [])
                feature_tally[feat].append(span)
                other -= span
            # Prepare the output for each tallied feature
            for feature in feature_tally:
                size_dist = feature_tally[feature]
                feat_n  = len(size_dist)
                total_l = sum(size_dist)
                prop_l  = total_l/chr_len
                # Initialize these to 0 to prevent stat errors
                mean_l   = 0
                median_l = 0
                sd_l     = 0
                min_l    = 0
                max_l    = 0
                # Now, if values are present calculate stats
                if feat_n>0:
                    mean_l   = statistics.mean(size_dist)
                    median_l = statistics.median(size_dist)
                    min_l    = min(size_dist)
                    max_l    = max(size_dist)
                if feat_n>1:
                    sd_l     = statistics.stdev(size_dist)
                # Save to output file
                row = [
                    chromosome,              # chromID
                    f'{chr_len}',            # chromLen
                    feature,                 # featType
                    f'{feat_n}',             # featNum
                    f'{total_l}',            # featLen
                    f'{prop_l:0.8f}',        # featProp
                    f'{mean_l:0.3f}',        # meanLen
                    f'{median_l:0.3g}',      # medianLen
                    f'{sd_l:0.3f}',          # sdLen
                    f'{min_l}',              # minLen
                    f'{max_l}' ]             # maxLen
                row = '\t'.join(row)
                fh.write(f'{row}\n')
            # Prepare output for other class
            feature  = 'other'
            total_l  = other
            prop_l   = total_l/chr_len
            feat_n   = -1
            mean_l   = -1
            median_l = -1
            sd_l     = -1
            min_l    = -1
            max_l    = -1
            row = [
                chromosome,              # chromID
                f'{chr_len}',            # chromLen
                feature,                 # featType
                f'{feat_n}',             # featNum
                f'{total_l}',            # featLen
                f'{prop_l:0.8f}',        # featProp
                f'{mean_l:0.3f}',        # meanLen
                f'{median_l:0.3g}',      # medianLen
                f'{sd_l:0.3f}',          # sdLen
                f'{min_l}',              # minLen
                f'{max_l}' ]             # maxLen
            row = '\t'.join(row)
            fh.write(f'{row}\n')

def main():
    print(f'{PROG} started on {date()} {time()}.')
    args = parse_args()
    # First, load the genome lengths from the FAI
    chromosomes = load_fai(args.fai)
    # Load and tally the phastCons bed
    phastcons = load_phastcons_bed(args.phastcons, chromosomes)
    calculate_phastcons_stats(phastcons, chromosomes, args.out_dir)

    # Load the annotation BED
    annotations = load_annotation_bed(args.annotation, chromosomes)
    # Process to find any overlapping windows.
    annotations = find_overlapping_annotations(annotations,
                                               chromosomes,
                                               args.out_dir)
    # Get some stats on this new file.
    calculate_annotation_stats(annotations, chromosomes, args.out_dir)

    # Calculate overlaps between phastCons and GFF
    # tally_phastcons_annotations(phastcons, annotations, chromosomes, args.out_dir)

    print(f'\n{PROG} finished on {date()} {time()}.')

# Run Code
if __name__ == '__main__':
    main()
