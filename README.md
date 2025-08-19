# random_scripts

Just for random scripts that might have no other place... for now.

## Script to process phastCon results along with annotations

```sh
$ python tally_phastcon_elements.py -h
  tally_phastcon_elements.py started on 2025-05-14 01:56:33.
  usage: tally_phastcon_elements.py [-h] -f FAI -g GFF -b BED [-o OUT_DIR]

  Determine the proportion of each of the annotated genetic elements in a GFF across the
  sites present in an phastCons conserved sites BED file. Provide other general stats
  for the phastCons BED.

  options:
    -h, --help            show this help message and exit
    -f FAI, --fai FAI     (str) Path to genome index in FAI format.
    -g GFF, --gff GFF     (str) Path to the annotation in GFF format.
    -b BED, --bed BED     (str) Path to the phastCons conserved sited BED.
    -o OUT_DIR, --out-dir OUT_DIR
```

### Output example

```sh
chromID  phastConsLen  featType         featLen  featProp
chr01    5410028       CDS              447997   8.281%
chr01    5410028       intron           2205940  40.775%
chr01    5410028       five_prime_UTR   6832     0.126%
chr01    5410028       three_prime_UTR  45805    0.847%
chr01    5410028       other            2703454  49.971%
chr02    5008219       CDS              357766   7.144%
chr02    5008219       intron           1669904  33.343%
chr02    5008219       five_prime_UTR   4539     0.091%
chr02    5008219       three_prime_UTR  32199    0.643%
```

## Extract CDS sequences for single-copy orthologs based on an Orthofinder run

Take the results from Orthofinder (>3.1.0), identify the transcript/gene IDs for 
all single-copy orthogroups, and extract their corresponding CDS sequences.

```sh
$ python3 extract_orthogroups_cds.py -h
  extract_orthogroups_cds.py started on 2025-07-11 16:32:07.
  usage: extract_orthogroups_cds.py [-h] -s SCO_LIST -r ORTHOGROUPS_TSV -c CDS_IN_DIR
                                  [-o OUT_DIR]

  options:
    -h, --help            show this help message and exit
    -s SCO_LIST, --sco-list SCO_LIST
                          (str) Path to orthofinder
                          Orthogroups/Orthogroups_SingleCopyOrthologues.txt file.
    -r ORTHOGROUPS_TSV, --orthogroups-tsv ORTHOGROUPS_TSV
                          (str) Path to the orthofinder Orthogroups/Orthogroups.tsv
    -c CDS_IN_DIR, --cds-in-dir CDS_IN_DIR
                          (str) Path to the directory containing the input per-taxon CDS
                          sequences.
    -o OUT_DIR, --out-dir OUT_DIR
                          (str) Path to output directory [default=./].
```

### Output table

This is reformatting the Orthofinder `Orthogroups/Orthogroups.tsv` to include only 
single-copy orthogroups, one gene/transcript per line.

```sh
#OrthogroupID  Taxon     TranscriptID
OG0004756      Species1  mrna-16824
OG0004756      Species2  mrna-9033
OG0004756      Species3  mrna-5845
OG0004756      Species4  mrna-6124
OG0004756      Species5  mrna-33735
OG0004765      Species1  mrna-32745
OG0004765      Species2  mrna-15598
OG0004765      Species3  mrna-13968
OG0004765      Species4  mrna-15215
OG0004765      Species5  mrna-345
```

### Output FASTA

An output FASTA is created for each single-copy orthogroup. The sequence IDs of each CDS
sequence correspond to their respective taxon (done for compatibility reasons).

```sh
$ cat out/OG0004756.cds.fa
  >Species1
  ATGATCACTGTCCTGCTGCCGGAGGAGCTGACAAGGCAGCAGCAGGGCTC
  >Species2
  ATGTCTCACAGGTTCCTCGAGTCTGTTAACGACTGCTTTCTCACCCAACA
  >Species3
  ATGATCACTGTGCTGCTGCCGGAGGAGCTGGCCGGCCAGCAGCAGGGCTC
  >Species4
  ATGCTCCGACCCGTCCCGACACAAGGTCCGACCCGTCCCGACACACGGTC
  >Species5
  ATGCAGCGCCACCGCTCGGTCGTCTGGCTGAGCTGCGCCACGGCGCTCGT
```

## Parse outputs from `HyPhy`

```sh
$ python3 parse_hyphy_outs.py -h
  parse_hyphy_outs.py started on 2025-08-19 12:09:30.
  usage: parse_hyphy_outs.py [-h] -s SCO_LIST -y HYPHY_OUTS [-o OUT_DIR]
                           [-a {aBSREL,BUSTED} [{aBSREL,BUSTED} ...]] [-m MIN_ALN_LEN]

  options:
    -h, --help            show this help message and exit
    -s SCO_LIST, --sco-list SCO_LIST
                          (str) Path to the single-copy orthgroup table (produced by
                          `extract_orthogroups_cds.py`).
    -y HYPHY_OUTS, --hyphy-outs HYPHY_OUTS
                          (str) Path to the HyPhy output files directory.
    -o OUT_DIR, --out-dir OUT_DIR
                          (str) Path to output directory [default=./].
    -a {aBSREL,BUSTED} [{aBSREL,BUSTED} ...], --analysis-type {aBSREL,BUSTED} [{aBSREL,BUSTED} ...]
                          (str) Type of HyPhy analysis to parse [choices: aBSREL, BUSTED]
                          [default=aBSREL].
    -m MIN_ALN_LEN, --min-aln-len MIN_ALN_LEN
                          (int) Minimum length required to keep an alignment [default=24]
```

## Author

Angel G. Rivera-Colón  
Institute of Ecology and Evolution  
University of Oregon
