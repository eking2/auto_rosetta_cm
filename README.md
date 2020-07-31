# Auto RosettaCM

Automated protein homology modeling with Rosetta Comparative Modeling. 
Performs a BLAST search to identify template PDBs with high homology to the query sequence, and generates the input files (MSA, threaded template PDBs, RosettaScript XML and flags) to run RosettaCM.

## Instructions

#### Installation
```
git clone https://github.com/eking2/auto_rosetta_cm.git
```

#### BLAST parameters
```
usage: main.py [-h] -f FASTA [-b] {cm} ...

positional arguments:
  {cm}
    cm                  Setup RosettaCM for homology modeling

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Fasta file with template protein amino acid sequence
  -b, --blast           Run blast to identify homologs
```

#### CM subparser parameters
```
usage: main.py cm [-h] [-n NUM_TEMPLATES] [-d DECOYS]
                  [-t [TEMPLATES [TEMPLATES ...]]]

optional arguments:
  -h, --help            show this help message and exit
  -n NUM_TEMPLATES, --num_templates NUM_TEMPLATES
                        Number of templates to use for homology modeling
                        (default: 5)
  -d DECOYS, --decoys DECOYS
                        Number of decoys to generate (default: 500)
  -t [TEMPLATES [TEMPLATES ...]], --templates [TEMPLATES [TEMPLATES ...]]
                        Target PDB templates to use, format as XXXX_Chain (ie.
                        4e5n_A)
```


## Examples

Run a BLAST search to find homologous crystal structures sorted by sequence identity
```
python main.py -f chmo.fasta -b
```

Select 6 non-duplicate templates automatically and make RosettaCM inputs
```
python main.py -f chmo.fasta cm -n 6
```

Setup RosettaCM with manually input template structures
```
python main.py -f chmo.fasta cm 3GWD_A 4RG3_A
```



