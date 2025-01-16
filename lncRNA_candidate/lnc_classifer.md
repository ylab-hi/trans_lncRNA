
## Prerequisite

We need `bedtools`, `gffcompare`, `stringtie` for lncRNA identification.

Bedtools, gffcompare and stringtie can be installed via conda.

```
conda install bioconda::bedtools
conda install bioconda::gffcompare
conda install bioconda::stringtie
```


## candiate novel lncRNA identification

### Step1: transcript assembly from RNA-seq data

transcript assembly for every sample, put them in the folder `transcripts`

```
stringtie --rf -G gencode/gencode.v37.gtf -f 0.01 sample.bam -o sample.gtf -A sample.tab -B -p 4
```

meta-assembly using all the samples, `gtf.list` includes all the sample.gtf in the previous step.

```
stringtie --merge --rf -G gencode/gencode.v37.gtf -c 1 -l G -p 16 -o merged.gtf gtf.list
```

### Step2: candidate novel lncRNA detection

```
lnc_classifer.py -i merged.gtf -f transcripts -o novel_lncRNAs
```
