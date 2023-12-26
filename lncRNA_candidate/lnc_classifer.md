
## Prerequisite

bedtools, gffcompare, stringtie


## candiate novel lncRNA identification

### Step1: transcript assembly from RNA-seq data

```
stringtie --rf -G gencode/gencode.v37.gtf -f 0.01 sample.bam -o sample.gtf -A sample.tab -B -p 4
```

```
stringtie --merge --rf -G gencode/gencode.v37.gtf -c 1 -l G -p 16 -o merged.gtf gtf.list
```


### Step2: candidate novel lncRNA detection
```
lnc_classifer.py -i ./meta-assembly/merge/merged.gtf -o novel_lncRNAs
```
