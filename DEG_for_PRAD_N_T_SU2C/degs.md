## Prerequisite
`featureCounts` and `edgeR` will be used for gene expression quantification.

```
conda install bioconda::subread
conda install bioconda::bioconductor-edger
conda install bioconda::bioconductor-limma
```

## Gene expression quantification

Gene expression in all RNA-seq samples was quantified using `featureCounts`, and the resulting raw read count expression matrix will be saved in the file named `expression.txt`.

### DEG identification using `edgeR`.

```
Rscript degs.R
```
