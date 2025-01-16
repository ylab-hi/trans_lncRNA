## Prerequisite

In this manuscript, we used `bedtools`, `CPAT`, `txCdsPredict`, `CPC2`, `CNCI2`, `LncADeep` for coding putential prediction.

Bedtools, CPAT and CPC2 can be installed via conda.
```
conda install bioconda::bedtools
conda install bioconda::cpat
conda install bioconda::cpc2

```
CNCI2 is available at [GitHub](https://github.com/www-bioinfo-org/CNCI).
LncADeep is available at [GitHub](https://github.com/cyang235/LncADeep).

txCdsPredict is available at [UCSC](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v385/txCdsPredict).

## Coding potential prediction
In order to run `coding_potential.py`, you have to install the dependencies, including Biopython, HTSeq, and pyfaidx.

bash```
conda install conda-forge::biopython
conda install bioconda::htseq
conda install bioconda::pyfaidx
```
### Running the program

bash```
coding_potential.py -g novel_lncRNAs.gtf -o novel_lncRNAs.table.txt

```
