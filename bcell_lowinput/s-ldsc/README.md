## Data:

### HapMap SNP list was downloaded from:

https://github.com/perslab/CELLECT/blob/master/data/ldsc/w_hm3.snplist

### LDSC files were obtained from:

https://zenodo.org/records/7768714

## Pipeline:

- Run setup.R to generate all input and auxiliary files;
- Run make_annot.slurm to generate annotations and compute LD score
- RUn ldsc.slurm to run S-LDSC and LDSC-SEG