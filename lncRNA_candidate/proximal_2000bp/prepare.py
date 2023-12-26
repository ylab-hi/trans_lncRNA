#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import re
import os
import glob
import HTSeq
import shutil

in_gtf = "gencode.v37.annotation.gff3"

gtf_file = HTSeq.GFF_Reader(in_gtf)
with open("proximal.b2000.hg38.bed", "w") as f:
    for feature in gtf_file:
        if feature.type == "gene":
            chrom = feature.iv.chrom
            start = feature.iv.start
            end = feature.iv.end
            gene_name = feature.attr["gene_name"]
            strand = feature.iv.strand
            new_start = start - 2000
            new_end = end + 2000
            if new_start < 0:
                new_start = 0
            f.write(f"{chrom}\t{new_start}\t{new_end}\t{gene_name}\t.\t{strand}\n")
