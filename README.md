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

## Author

Angel G. Rivera-Colón  
Institute of Ecology and Evolution  
University of Oregon
