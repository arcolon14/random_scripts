# random_scripts

Just for random scripts that might have no other place... for now.

## Script to process phastCon results along with annotations

```sh
$ python3 tally_phastcon_elements.py -h
  tally_phastcon_elements.py started on 2025-09-09 14:22:42.
  usage: tally_phastcon_elements.py [-h] -f FAI -a ANNOTATION -p PHASTCONS
                                    [-o OUT_DIR] [-m MIN_SEQ_LEN]
                                    [-i MIN_INTERVAL_LEN]

  Determine the proportion of each of the annotated genetic elements in a
  BED across the sites present in an phastCons conserved sites BED file.
  Provide other general stats for the phastCons BED.

  options:
    -h, --help            show this help message and exit
    -f FAI, --fai FAI     (str) Path to genome index in FAI format.
    -a ANNOTATION, --annotation ANNOTATION
                          (str) Path to the annotation in BED format.
    -p PHASTCONS, --phastcons PHASTCONS
                          (str) Path to the phastCons conserved sited BED.
    -o OUT_DIR, --out-dir OUT_DIR
                          (str) Path to output directory [default=.].
    -m MIN_SEQ_LEN, --min-seq-len MIN_SEQ_LEN
                          (int|float) Min length of input sequences
                          [default=10,000].
    -i MIN_INTERVAL_LEN, --min-interval-len MIN_INTERVAL_LEN
                          (int) Min length of intervals in input BED files
                          [default=1].
```

### Example of input annotation BED

The input annotation BED follows the standard BED format:

1. sequence ID
2. start coordinate (0-based, inclusive)
3. end coordinate (0-based, exclusive)
4. feature ID (e.g., column 3 of a GFF/GTF)

```sh
chr01   93769   93936   TE
chr01   95152   95296   CDS
chr01   95152   95298   lncRNA
chr01   95296   95709   intron
chr01   95478   95524   TE
chr01   95630   95661   TE
chr01   95709   95866   CDS
chr01   95866   96140   intron
chr01   96140   96448   CDS
chr01   96378   96424   TE
chr01   96448   97002   intron
chr01   97002   97192   CDS
chr01   97192   97920   intron
```

### Output examples

First, the `phastCons_stats.tsv` file contains a breakdown of the distribution of
`phastCon` conserved intervals in the genome. The statistics are reported per-chromosome,
with a `AllGW` category reporting the distribution genome-wide.

```sh
chromID  chromLen    phastConsNum  phastConsLen  phastConsProp  meanLen  medianLen  sdLen     minLen  maxLen
chr01    80341233    25403         5410028       0.067338       212.968  38         908.569   1       21151
chr02    78161376    20391         5008219       0.064075       245.609  46         913.303   1       21350
chr03    63686915    11125         3408288       0.053516       306.363  44         1107.898  1       19916
...
AllGW    1033402680  196141        50849898      0.049206       259.252  42         990.404   1       29414
```

Secondly, the `annotation_stats.tsv` file reports the distribution of all annotation
features along the genome. This is done per-chromosome, with a `AllGW` category reporting
the genome-wide distributions.

```sh
chromID  chromLen    featType         featNum  featLen    featProp    meanLen   medianLen  sdLen     minLen  maxLen
chr01    80341233    CDS              6994     1819439    0.02264639  260.143   159        445.858   1       13941
chr01    80341233    TE               34842    33266476   0.41406479  954.781   387        1776.901  1       71192
chr01    80341233    five_prime_UTR   554      50504      0.00062862  91.162    54         110.332   1       832
chr01    80341233    intron           26380    12176014   0.15155374  461.562   202        907.093   1       33320
chr01    80341233    lncRNA           907      180908     0.00225175  199.458   67         318.948   1       2072
chr01    80341233    multiple         21941    14875054   0.18514844  677.957   296        1191.398  1       37401
chr01    80341233    other            36136    17707070   0.22039829  490.012   153        1071.050  1       48094
chr01    80341233    rRNA             2        1212       0.00001509  606.000   606        709.935   104     1108
chr01    80341233    tRNA             35       1040       0.00001294  29.714    8          31.578    3       84
chr01    80341233    three_prime_UTR  525      263516     0.00327996  501.935   380        443.673   1       2343
...
AllGW    1033402680  CDS              67111    17613714   0.01704439  262.456   155        418.576   1       15837
AllGW    1033402680  TE               476627   494889569  0.47889325  1038.316  405        1923.924  1       73062
AllGW    1033402680  five_prime_UTR   5695     564202     0.00054597  99.070    57         144.502   1       3439
AllGW    1033402680  intron           267933   120071863  0.11619078  448.141   193        919.243   1       68669
AllGW    1033402680  lncRNA           11887    2714087    0.00262636  228.324   102        367.418   1       10804
AllGW    1033402680  multiple         228641   156456043  0.15139891  684.287   297        1230.938  1       54289
AllGW    1033402680  other            490304   238571963  0.23086060  486.580   149        1122.196  1       106410
AllGW    1033402680  rRNA             13       2460       0.00000238  189.231   112        276.124   101     1108
AllGW    1033402680  tRNA             607      25166      0.00002435  41.460    37         32.147    1       88
AllGW    1033402680  three_prime_UTR  5349     2493613    0.00241301  466.183   339        435.273   1       3008
```

Note that this distribution includes two new features that might not be present in the input
annotations:
1. `other`: intervals of the genome without any annotations (this includes any gaps).
2. `multiple`: intervals of the genome with more than one annotation (e.g., overlap between
an intron and a transposable element).

Lastly, the `phastCons_annotations.tsv` file describes the overlap between `phastCon`
conserved regions and annotations. This is done per-chromosome, with a `AllGW` reporting
the genome-wide distributions.

```sh
chromID  phastConsLen  featType         featLen   featProp
chr01    5410028       CDS              321620    0.05944886
chr01    5410028       TE               1548507   0.28622902
chr01    5410028       five_prime_UTR   2865      0.00052957
chr01    5410028       intron           1059617   0.19586165
chr01    5410028       lncRNA           13186     0.00243733
chr01    5410028       multiple         1350787   0.24968207
chr01    5410028       other            1089707   0.20142354
chr01    5410028       rRNA             0         0.00000000
chr01    5410028       tRNA             394       0.00007283
chr01    5410028       three_prime_UTR  23345     0.00431513
...
AllGW    50849898      CDS              2225138   0.04375895
AllGW    50849898      TE               17505920  0.34426657
AllGW    50849898      five_prime_UTR   27679     0.00054433
AllGW    50849898      intron           8575960   0.16865245
AllGW    50849898      lncRNA           130744    0.00257118
AllGW    50849898      multiple         11537936  0.22690185
AllGW    50849898      other            10688761  0.21020221
AllGW    50849898      rRNA             0         0.00000000
AllGW    50849898      tRNA             5162      0.00010151
AllGW    50849898      three_prime_UTR  152598    0.00300095
```

This output files reports the total number of conserved sites, and the proportion
of those sites belonging to the input annotation classes. Using the example above,
of ~5.4 million `phastCon` conserved sites, 321 thousand (~6%) correspond to coding
sequences, 1.05 million (~20%) correspond to introns, etc. These values can be 
compared to the genome-wide distribution (in `annotation_stats.tsv`) to find
enrichment of any annotaion class among evolutionarily-conserved sites.

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
