#!/usr/bin/env python
#-*- coding: utf-8 -*-
import glob
import re
import csv
import sys
import subprocess
import os
import argparse
import shutil
import random
import hashlib
from collections import Counter, defaultdict
import HTSeq
import shutil


CHROMS = {'chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
     'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
     'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'}

def status_message(msg):
    print(msg)
    sys.stdout.flush()

def remove(infile):
    if os.path.isfile(infile):
        os.remove(infile)
    elif os.path.isdir(infile):
        shutil.rmtree(infile)

def run_cmd(cmd, msg=None):
    '''
    '''
    status_message(cmd)
    if ',' in msg:
        begin, finish = msg.split(',')
        status_message(begin)
    else:
        finish = msg
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'Error happend!: {}\n{}'.format(err, err.output)
    else:
        error_msg = ''
    if not error_msg:
        status_message(finish)
        return True
    else:
        status_message(error_msg)
        return False


def external_tool_checking():
    """checking dependencies are installed"""
    software = ['gffcompare', 'bedtools']
    cmd = "which"
    for each in software:
        try:
            path = subprocess.check_output([cmd, each], stderr=subprocess.STDOUT)
            path = str(path, 'utf-8')
        except subprocess.CalledProcessError:
            sys.stderr.write("Checking for {0} : ERROR - could not find {0}".format(each))
            sys.stderr.write("Exiting.")
            sys.exit(0)
        print(("Checking for '" + each + "': found " + path))

def get_value_by_key(src, delimiter='='):
    out_dict = {}
    src = src.strip().strip(';')
    tmpList = re.split(r';\s{0,2}', src) 
    for i in tmpList:
        if i:
            m, n = i.split(delimiter)
            out_dict[m] = n
    return out_dict

def mapability_filter(in_gtf, cutoff=0.04, map_blacklist='./blacklist_regions/LCR.blacklist.hg38.bed'):
    '''
    Step 1: mapability filtering based LCR from Heng Li and blacklist region from ENCODE project
    :param in_gtf: Input GTF
    :type in_gtf: str
    :param cutoff: coverage cutoff
    :type cutoff: float
    :return: (out_gtf, removal_transcripts_list)
    :rtype: tuple
    '''
    pid = os.getpid()
    coverage_file = 'mapability.coverage.{}.bed'.format(pid)
    cmd = 'bedtools coverage -a {} -b {} > {}'.format(in_gtf, map_blacklist, coverage_file)
    ret = run_cmd(cmd, 'Mapability filtering begins!')

    retained_transcript = {}
    removal_transcript = {}
    # remove transcripts > %4 of which are overed by mapability regions
    if ret:
        with open(coverage_file, 'r') as f:
            for line in f:
                l = line.rstrip('\n').split('\t')
                chrm = l[0]
                id_type = l[2]
                coverage = float(l[-1])
                anno = l[8]
                assembly_dict = get_value_by_key(anno,delimiter=' ')
                if id_type == 'transcript':
                    txid = assembly_dict['transcript_id'].strip('"')
                    if chrm in CHROMS and coverage < cutoff:
                        retained_transcript[txid] = True
                    else:
                        removal_transcript[txid] = True
    #remove(coverage_file)
    # gtf file excluding mapbility regions
    gtf_exclu = '{}.reliable.gtf'.format(pid)
    output = open(gtf_exclu, 'w')
    with open(in_gtf, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                l = line.rstrip('\n').split('\t')
                anno = l[8]
                assembly_dict = get_value_by_key(anno,delimiter=' ')
                txid = assembly_dict['transcript_id'].strip('"')
                if txid in retained_transcript:
                    output.write(line)
    output.close()
    removal_transcripts_list = list(removal_transcript.keys())
    print('mapability filtering finished!\ngenerated {}\n'.format(gtf_exclu))
    print('Removed {} transcripts\n'.format(len(removal_transcripts_list)))
    print('{} transcripts retained in mapability step!\n'.format(len(list(retained_transcript.keys()))))

    return gtf_exclu, removal_transcripts_list


def length_filter(in_gtf):
    ''' Step 3
        transcripts shorter than 200 bp and without strand infromation were discarded
    '''
    trx_length = defaultdict(int)
    trx_strand = {}
    gtf_file = HTSeq.GFF_Reader(in_gtf)
    for feature in gtf_file:
        if feature.type == 'exon':
            tx_id = feature.attr['transcript_id']
            _length = int(feature.iv.length)
            _strand = feature.iv.strand
            trx_length[tx_id] += _length
            trx_strand[tx_id] = _strand

    retained_transcript = set()
    removed_transcript = set()
    for i in trx_length:
        if trx_length[i] > 200 and trx_strand[i] in {'+', '-'}:
            retained_transcript.add(i)
        else:
            removed_transcript.add(i)
    print('Transcript length filtering (<200bp) finished!\n')
    print('Removed {} transcripts that length shorter than 200bp or without strand information\n'.format(len(removed_transcript)))
    return removed_transcript

def gff_comparison(in_gtf, ref_gtf):
    ## gffcompare -R -r /home/tywang/database/lncipedia/protein_coding.lncipedia.hg19.gff3 -o cmp assembly.gtf
    '''
    compare with annotations including protein-coding genes, canonical ncRNAs (gencode) and long ncRNAs (lncipedia)
    '''
    uid = random.getrandbits(25)
    # tmap and refmap files are in the same folder as the query GTF file
    # loci, tracking, stats and annotated.gtf files are in the current directory 
    map_dir = os.path.dirname(in_gtf)
    gtf_name = os.path.basename(in_gtf)

    cwd = os.getcwd()

    tmap_file = os.path.join(map_dir, f'{uid}.{gtf_name}.tmap')
    refmap_file = os.path.join(map_dir, f'{uid}.{gtf_name}.refmap')

    cmd = f'gffcompare -R -r {ref_gtf} -o {uid} {in_gtf}'
    ret = run_cmd(cmd, 'gffcompare is running!\n')
    if ret:
        remove(os.path.join(cwd, f'{uid}.loci'))
        remove(os.path.join(cwd, f'{uid}.stats'))
        remove(os.path.join(cwd, f'{uid}.tracking'))
        remove(os.path.join(cwd, f'{uid}.annotated.gtf'))
        remove(refmap_file)
        return tmap_file
    else:
        return False


def trace_back_filter(in_gtf, folder):
    '''traceback the merged transcripts in individual assembled GTF files,
       if the number of supporting samples less than 3 OR FPKM < 1.0; then discard this transcript
    '''
    ref_tx_ids = set([])
    gtf_file = HTSeq.GFF_Reader(in_gtf)
    for feature in gtf_file:
        if feature.type == 'transcript':
            ref_tx_ids.add(feature.attr['transcript_id'])

    ref_tx_dict = dict((k, set([])) for k in ref_tx_ids)
    tx_traceback_dict = dict((k, set([])) for k in ref_tx_ids)
    matched_codes = {'=', 'c', 'k', 'j'}
    best_matched_codes = {'=', 'c', 'k'}
    path = os.path.abspath(folder)

    for g in glob.iglob(os.path.join(path, '*.gtf')):

        expressed_tx = set([])
        ind_gtf = HTSeq.GFF_Reader(g)
        for feature in ind_gtf:
            #if feature.type == 'transcript' and float(feature.attr['FPKM']) >=1.0 and float(feature.attr['TPM']) >= 1.0:
            if feature.type == 'transcript' and float(feature.attr['FPKM']) >=1.0:
                expressed_tx.add(feature.attr['transcript_id'])

        ind_name = os.path.splitext(os.path.basename(g))[0]
        tmap = gff_comparison(in_gtf=g, ref_gtf=in_gtf)
        with open(tmap) as f:
            f.readline()
            for line in f:
                l = line.rstrip('\n').split('\t')
                ref_id = l[1]
                class_code = l[2]
                qry_id = l[4]
                if class_code in best_matched_codes and ref_id != '-' and qry_id != '-' and qry_id in expressed_tx:
                    tx_traceback_dict[ref_id].add('{}:{}'.format(qry_id, ind_name))
                    ref_tx_dict[ref_id].add(ind_name)
        remove(tmap)
    removed_set = set([])
    uid = os.getpid()
    out = open(f'{uid}.traceback', 'w')
    for j in ref_tx_dict:
        if len(ref_tx_dict[j]) < 5:
            removed_set.add(j)
        else:
            out.write('{}\t{}\n'.format(j, ','.join(tx_traceback_dict[j])))
    out.close()
    print('Removed {} transcripts that are low expressed in individual samples and appears less than 3 samples\n'.format(len(removed_set)))
    return removed_set


def transcripts_remove_from_gtf(ids, in_gtf, out_gtf):
    ''' remove transcripts according to input "ids";
        only keep canonical chromosomes
    '''
    uniq_ids = set([])
    for i in ids:
        uniq_ids.add(i)

    input_trx = set([])
    output_trx = set([])

    output = open(out_gtf, 'w')
    gtf_file = HTSeq.GFF_Reader(in_gtf)
    for feature in gtf_file:
        if feature.type == 'transcript' or feature.type == 'exon':
            trx_id = feature.attr['transcript_id']
            input_trx.add(trx_id)
            if not trx_id in uniq_ids:
                output.write(feature.get_gff_line())
                output_trx.add(trx_id)
    output.close()
    print('{} generated!'.format(out_gtf))
    print('Input: {} transcripts, Output: {} transcripts.'.format(len(input_trx), len(output_trx)))
    print('{} transcripts removed!'.format(len(uniq_ids)))


def annotation_filter(in_gtf, ref_gtf):
    '''
      transcripts that meet the classification codes: =, c, k, j of gffcompare results
      on (protein-coding, long-noncoding and small-noncoding) in annotations.

      return : assembled transcripts have hit with annotations
    '''
    coding_genes = {'protein_coding'}
    lncRNAs = {"processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding",
        "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA",
        "bidirectional_promoter_lncrna"}
    can_ncRNAs = {"snRNA", "snoRNA", "rRNA", "Mt_tRNA", "Mt_rRNA", "misc_RNA", "miRNA", "ribozyme", "sRNA",
        "scaRNA", "vaultRNA"}
    gtf_file = HTSeq.GFF_Reader(ref_gtf)
    trx_type = {}
    for feature in gtf_file:
        if feature.type == 'exon':
            trx_id = feature.attr['transcript_id']
            gene_id = feature.attr['gene_name']
            gene_biotype = feature.attr['gene_type']
            if gene_biotype in coding_genes:
                trx_type[trx_id] = 'protein_coding'
            elif gene_biotype in can_ncRNAs:
                trx_type[trx_id] = 'canonical_ncRNA'
            elif gene_biotype in lncRNAs:
                trx_type[trx_id] = 'lncRNA'
            else:
                trx_type[trx_id] = 'other'
    print(trx_type)
    ## ref_gene_id  ref_id  class_code  qry_gene_id  qry_id
    tmap_file = gff_comparison(in_gtf, ref_gtf)
    gene_dict = {'protein_coding':set([]),'canonical_ncRNA':set([]), 'lncRNA':set([]), 'novel':set([]), 'other':set([])}
    trx_dict = {'protein_coding':set([]),'canonical_ncRNA':set([]), 'lncRNA':set([]), 'novel':set([]), 'other':set([])}
    removal_codes = {'=', 'c', 'k', 'j'}
    best_matched_codes = {'=', 'c', 'k'}

    with open(tmap_file, 'r') as f:
        f.readline()
        for line in f:
            l = line.rstrip('\n').split('\t')
            ref_gene = l[0]
            ref_trx = l[1]
            class_code = l[2]
            qry_gene = l[3]
            qry_trx = l[4]
            if class_code in best_matched_codes:
                if ref_trx in trx_type:
                    group_name = trx_type[ref_trx]
                    gene_dict[group_name].add(qry_gene)
                    trx_dict[group_name].add(qry_trx)
            else:
                gene_dict['novel'].add(qry_gene)
                trx_dict['novel'].add(qry_trx)
    print(gene_dict['protein_coding'])
    #remove(tmap_file)
    removed_transcripts = set([])
    removed_transcripts.update(trx_dict['protein_coding'])
    removed_transcripts.update(trx_dict['canonical_ncRNA'])
    removed_transcripts.update(trx_dict['lncRNA'])
    removed_transcripts.update(trx_dict['other'])

    print('Annotation filtering finished!\n')
    print('Removed {} genes, {} transcript overlapped with protein-coding genes'.format(len(gene_dict['protein_coding']), len(trx_dict['protein_coding'])))
    print('Removed {} genes, {} transcript overlapped with canonical ncRNA'.format(len(gene_dict['canonical_ncRNA']), len(trx_dict['canonical_ncRNA'])))
    print('Removed {} genes, {} transcript overlapped with lncRNAs'.format(len(gene_dict['lncRNA']), len(trx_dict['lncRNA'])))
    print('Removed {} genes, {} transcript overlapped with others'.format(len(gene_dict['other']), len(trx_dict['other'])))
    print('{} genes, {} transcript: canndidate novel lncRNAs retained!'.format(len(gene_dict['novel']), len(trx_dict['novel'])))

    return removed_transcripts


def proximal_monoexon_filter(in_gtf, ref_2000_bed="proximal_2000bp/proximal.b2000.hg38.bed"):
    '''
    Step 4
    signle-exon transcript proximal (2000bp) to protein-coding genes or ncRNAs on the same strand

    '''
    # get monoexonic transcripts
    trx_dict = defaultdict(int)

    gtf_file = HTSeq.GFF_Reader(in_gtf)
    for feature in gtf_file:
        if feature.type == 'exon':
            trx_id = feature.attr['transcript_id']
            trx_dict[trx_id] += 1

    # write monoexonic transcripts GTF file
    uid = os.getpid()
    monoexon_gtf = open(f'{uid}.monoexon.gtf', 'w')
    gtf_file = HTSeq.GFF_Reader(in_gtf)
    for feature in gtf_file:
        trx_id = feature.attr['transcript_id']
        if feature.type == 'exon' and trx_dict[trx_id] == 1:
            monoexon_gtf.write(feature.get_gff_line())
    monoexon_gtf.close()
    # overlap monoexonic transcripts with annotation 2000bp file
    overlap_file = f'{uid}.monoexonic.overlap.bed'
    cmd = 'bedtools intersect -s -wa -wb -a {} -b {} > {}'.format(ref_2000_bed, monoexon_gtf.name, overlap_file)
    ret = run_cmd(cmd, 'proximal monoexoinc transcripts begins!')

    # get the monexonic transcripts that overlapped with +- 2000bp annotation
    removed_transcripts = set([])
    if ret:
        with open(overlap_file, 'r') as f:
            for line in f:
                l = line.rstrip('\n').split('\t')
                assembly_anno = l[14]
                assembly_dict = get_value_by_key(assembly_anno,delimiter=' ')
                trx_id = assembly_dict['transcript_id'].strip('"')
                removed_transcripts.add(trx_id)

    remove(overlap_file)
    remove(monoexon_gtf.name)
    print('proximal monoexon filtering finished!\n')
    print('Removed {} monoexonic transcripts that are proximal to annotated genes\n'.format(len(removed_transcripts)))
    return monoexon_gtf.name, removed_transcripts


def print_list(inlist,name):
    out = open(name, 'w')
    for i in inlist:
        out.write(i + '\n')
    out.close()


def parse_args():
    parser = argparse.ArgumentParser(description="candidate lncRNA identification")
    parser.add_argument('-i', '--input', action='store', dest='input', help="RNA assembly GTF file", required=True)
    parser.add_argument('-o', '--output', action='store', dest='output', help="fitered GTF file for coding potential test", required=True)
    parser.add_argument('-c', '--cutoff', action='store', dest='cutoff', help="mapability filter cutoff (default: %(default)s)",type=float,default=0.04)
    parser.add_argument('-r', '--ref', action='store', dest='ref', help="hg19 or hg38 (default: %(default)s)", choices=['hg19','hg38'], default='hg38')
    parser.add_argument('-f', '--folder', action='store', dest='folder', help="folder of tanscript assemblies (default: %(default)s)", default='./transcripts/')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()
    return parser


def main():
    if sys.version_info < (3,4):
        sys.exit('Sorry, this code need Python 3.4 or higher. Please update. Aborting...')
    parser = parse_args()

    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        options = parser.parse_args()

    if options.ref == 'hg38':
        anno_gff = 'gencode/gencode.v37.annotation.gff3'
        ref_gene_2000_bed = 'proximal_2000bp/proximal.b2000.hg38.bed'
        map_blacklist='blacklist_regions/LCR.blacklist.hg38.bed'

    external_tool_checking()
    # Mapability region
    mapa_gtf, mapa_list = mapability_filter(options.input, cutoff=options.cutoff, map_blacklist=map_blacklist)
    # length > 200bp
    short_len_list = length_filter(mapa_gtf)
    # FPKM > 1 and Number of samples >= 3
    low_exp_sample_list = trace_back_filter(mapa_gtf, options.folder)
    exclu_ids = list(short_len_list) + list(low_exp_sample_list)

    uid = os.getpid()
    transcripts_remove_from_gtf(ids=exclu_ids, in_gtf=mapa_gtf, out_gtf=f'{uid}.expressed.gtf')
    # monoexon proximal to genes (2000bp)
    monoexon_gtf, monoexon_list = proximal_monoexon_filter(in_gtf=f'{uid}.expressed.gtf', ref_2000_bed=ref_gene_2000_bed)

    transcripts_remove_from_gtf(ids=monoexon_list, in_gtf=f'{uid}.expressed.gtf', out_gtf=f'{uid}.proxmono.removed.gtf')

    # annotated coding genes, canonical ncRNAs, lncRNAs and pseduogenes
    known_list = annotation_filter(in_gtf=f'{uid}.proxmono.removed.gtf', ref_gtf=anno_gff)
    transcripts_remove_from_gtf(ids=known_list, in_gtf=f'{uid}.proxmono.removed.gtf', out_gtf=f'{options.output}.CP.gtf')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
