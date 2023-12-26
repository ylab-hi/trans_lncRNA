## LCR region from Heng Li
https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs37d5.bed.gz
https://figshare.com/articles/dataset/Low_complexity_regions_in_hs37d5/969685

### Liftover LCR region
```bash
CrossMap.py bed GRCh37ToHg19.over.chain.gz LCR-hs37d5.bed LCR-hs37d5.hg19.bed --chromid l
CrossMap.py bed hg19ToHg38.over.chain.gz LCR-hs37d5.hg19.bed LCR-hs37d5.hg38.bed
```

## blacklist from ENCODE
https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz


## Merge blacklist regions

```bash
cat LCR-hs37d5.hg38.bed hg38-blacklist.v2.bed |sort -V -k1,1 -k2,2|uniq |bedtools merge -i - > LCR.blacklist.hg38.bed
```
