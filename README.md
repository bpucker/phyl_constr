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

`--baits` specifies a forward read input FASTQ file.
