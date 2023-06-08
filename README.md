# phyl_constr
construction of species tree based on large sequence data sets


## Phylogenetic tree constructor ##
Construct a phylogenetic tree based on all sequences across the species of interest.

```
python3 phyl_constr.py --baits <FILE> --in <FILE> --out <DIRECTORY>

Arguments:
--baits        STR    Bait sequence file
--in           STR    Config file
--out          STR    Output folder

optional:
--mode         STR    Analysis mode (cds|prot)[cds]

--blastp       STR    Path to blastp[blastp]
--makeblastdb  STR    Path to makeblastdb[makeblastdb]
--fasttree     STR    Path to FastTree2[FastTree]
--mafft        STR    Path to MAFFT[mafft]
				
--minsim       FLOAT  Minimal BLASTp hit similarity[80.0]
--minlen       INT    Minimal BLASTp hit length[50]
--repratio     FLOAT  Minimal ratio of represented species per gene tree[1.0]
--occupancy    FLOAT  Minimal alignment occupancy[0.1]
```

`--baits` specifies a FASTA input file.

`--in` specifies a config input file.

`--out` specifies an output folder. This folder will be created if it does not exist already.

`--mode` specifies the analysis mode. This needs to indicate whether the input comprizes coding sequences (cds) or polypeptide sequences (prot). Default: cds.

`--blastp` specifies the full path to the blastp binary. This variable can be used if blastp is not accessbile through your PATH variable. Default: blastp.

`--makeblastdb` specifies the full path to the makeblastdb binary. This variable can be used if blastp is not accessbile through your PATH variable. Default: makeblastdb.

`--fasttree` specifies the full path to the FastTree2 binary. This variable can be used if FastTree2 is not accessbile through your PATH variable. Default: fasttree.

`--mafft` specifies the full path to the MAFFT binary. This variable can be used if MAFFT is not accessbile through your PATH variable. Default: mafft.

`--minsim` specifies the minimal similarity of a BLAST hit to be considered in the following analysis. Default: 80.0.

`--minlen` specifies the minimal length of a BLAST hit to be considered in the following analysis. Default: 50.

`--repratio` specifies the number of species that need to contribute an ortholog to a set of sequences. Only sets of orthologous sequences with a sufficient ratio of representated species are considered for the phylogenetic tree construction. Default: 1.0 (all).

`--occupancy` specifies the minimal alignment occupancy. This determines which columns will be removed during the alignment trimming step. Default: 0.1.

## References
This repository.
