#!/usr/bin/env python
#-*- coding: utf-8 -*-
import glob
import re
import csv
import subprocess
import os
import sys
import argparse
import shutil
from collections import Counter,defaultdict
from Bio import SeqIO
import HTSeq
from pyfaidx import Fasta


def get_value_by_key(src, delimiter='='):
    out_dict = {}
    src = src.strip(';')
    tmpList = re.split(r';\s{0,2}', src) 
    for i in tmpList:
        if i:
            m, n = i.split(delimiter)
            out_dict[m] = n
    return out_dict

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
    software = ['bedtools', 'cpat.py', 'txCdsPredict', 'CPC2.py', 'CNCI2.py', 'LncADeep.py']
    cmd = "which"
    for each in software:
        try:
            path = subprocess.check_output(
                [cmd, each], stderr=subprocess.STDOUT)
            path = str(path, 'utf-8')
        except subprocess.CalledProcessError:
            print("Checking for '" + each +
                  "': ERROR - could not find '" + each + "'", file=sys.stderr)
            print("Exiting.", file=sys.stderr)
            sys.exit(0)
        print("Checking for '" + each + "': found " + path)


def GTF2FASTA(in_gtf, ref_genome, out_fa):

    try:
        genome_fasta = Fasta(ref_genome, sequence_always_upper=True)
    except FastaNotFoundError as e:
        print('read reference genome ' + ref_genome + ' error!', e)
        sys.exit(1)

    trx_to_exons = defaultdict(list)
    gtf_file = HTSeq.GFF_Reader(in_gtf)
    for feature in gtf_file:
        trx_id = feature.attr['transcript_id']
        if feature.type == 'exon':
            trx_to_exons[trx_id].append(feature.iv)

    sorted_trx_to_exons = {}
    for trx_id in trx_to_exons:
        tmpList = trx_to_exons[trx_id]
        tmpList.sort(key=lambda x:x.start)
        sorted_trx_to_exons[trx_id] = tmpList
    trx_to_exons = None

    fasta_file = open(out_fa, 'w')
    for trx in sorted_trx_to_exons:
        exons = sorted_trx_to_exons[trx]
        strand = exons[0].strand
        seq = ''
        if strand == '+':
            for exon in exons:
                seq += genome_fasta[exon.chrom][exon.start:exon.end].seq
        elif strand == '-':
            for exon in exons[::-1]:
                seq += genome_fasta[exon.chrom][exon.start:exon.end].reverse.complement.seq
        fasta_file.write(f'>{trx}\n')
        fasta_file.write(f'{seq}\n')
    fasta_file.close()
    return out_fa


def genLncRNASeq(gtf, ref_fa="hg38.fa"):
    uid = os.getpid()

    out = open(f'{uid}.exons.gtf', 'w')
    with open(gtf) as f:
        for line in f:
            l = line.rstrip().split('\t')[2]
            if l[2] == 'exon':
                out.write(line)
    out.close()

    #ret = subprocess.check_call('grep -P "\texon\t" {} > {}.exons.gtf'.format(gtf, uid), shell=True)
    lncRNA_fa = '{}.tmp.lncRNA.fa'.format(uid)
    cmd = 'bedtools getfasta -s -fi {0} -bed {1}.exons.gtf -fo {2}'.format(ref_fa, uid, lncRNA_fa)
    ret = run_cmd(cmd, 'Generating FASTA!')
    if ret:
        fasta = {}
        for record in SeqIO.parse(lncRNA_fa, "fasta"):
            fasta[record.id] = str(record.seq)

        tx_dict = {}
        exon_dict = {}
        exon_gtf = f'{uid}.exons.gtf'
        with open(exon_gtf) as f:
            for line in f:
                l = line.rstrip('\n').split('\t')
                chrm = l[0]
                _type = l[2]
                strand = l[6]
                anno = l[8]
                start = int(l[3]) - 1
                end = int(l[4])
                f_dict = get_value_by_key(anno, delimiter=' ')
                tx_id = f_dict['transcript_id'].strip('"')
                gene_id = f_dict['gene_id'].strip('"')
                fasta_id = '{}:{}-{}({})'.format(chrm, start, end, strand)
                tx_dict[tx_id] = gene_id
                if tx_id in exon_dict:
                    exon_dict[tx_id].append(fasta_id)
                else:
                    exon_dict[tx_id] = [fasta_id]
        out_fa_name = '{}.lncRNA.fa'.format(os.getpid())
        out_fa = open(out_fa_name, 'w')
        for idx in exon_dict:
            # transcript_id - gene_id
            new_id = '{}-{}'.format(idx, tx_dict[idx])
            new_fasta_list = sorted(exon_dict[idx], key=lambda x: x.split(':')[1].split("-"))
            new_fasta = ''
            for j in new_fasta_list:
                if j in fasta:
                    new_fasta = new_fasta + fasta[j]
            out_fa.write('>'+new_id+'\n')
            out_fa.write(new_fasta+'\n')
        out_fa.close()

    remove(lncRNA_fa)
    remove(exon_gtf)

    return out_fa_name

def cpat_runner(in_fasta, cutoff=0.364):
    '''
    :param in_fasta: transcript DNA sequence
    :type in_fast: str
    :return: out_dict: transcript_id => coding probility
    :rtype: dict
    '''
    model = 'cpat/Human_logitModel.RData'
    hexamer =  'cpat/Human_Hexamer.tsv'
    result_id = '{}'.format(os.getpid())
    result = f'{result_id}.ORF_prob.best.tsv'
    result_all = f'{result_id}.ORF_prob.tsv'
    orf_seqs = f'{result_id}.ORF_seqs.fa'
    no_orf = f'{result_id}.no_ORF.txt'
    r_script = f'{result_id}.r'
    cmd = 'cpat.py -g {} -d {} -x {} -o {}'.format(in_fasta, model, hexamer, result_id)
    ret = run_cmd(cmd, 'CPAT done!')
    out_dict = {}
    if ret:
        with open(result) as f:
            f.readline()
            for line in f:
                l = line.rstrip('\n').split('\t')
                uid = l[0]
                prob = float(l[-1])
                out_dict[uid] = prob
                if prob >= cutoff:
                    judge = 'coding'
                else:
                    judge = 'noncoding'
                out_dict[uid] = {'prob':prob, 'label':judge}
    remove(result)
    remove(result_all)
    remove(orf_seqs)
    remove(no_orf)
    remove(r_script)
    remove('CPAT_run_info.log')
    return out_dict

def txCdsPredict_runner(in_fasta, cutoff=800):
    '''
    :param in_fasta: transcript DNA sequence
    :type in_fast: str
    :return: out_dict: transcript_id => coding score
    :rtype: dict
    '''
    result = '{}.txCdsPredict.result'.format(os.getpid())
    cmd = 'txCdsPredict {} {}'.format(in_fasta, result)
    ret = run_cmd(cmd, 'txCdsPredict done!')
    out_dict = {}
    if ret:
        with open(result) as f:
            for line in f:
                l = line.rstrip('\n').split('\t')
                uid = l[0]
                score = float(l[5])
                if score >= cutoff:
                    judge = 'coding'
                else:
                    judge = 'noncoding'
                out_dict[uid] = {'prob':score, 'label':judge}
    remove(result)
    return out_dict


def LncADeep_runner(in_fasta, cutoff=0.5):
    '''
    :param in_fasta: transcript DNA sequence
    :type in_fast: str
    :return: out_dict: transcript_id => coding probility
    :rtype: dict
    '''
    result_id = '{}'.format(os.getpid())
    result = f'{result_id}_LncADeep_lncRNA_results/{result_id}_LncADeep.results'
    cmd = 'LncADeep.py -MODE lncRNA -th 4 -m full -f {} -o {}'.format(in_fasta, result_id)
    ret = run_cmd(cmd, 'LncADeep done!')
    out_dict = {}
    if ret:
        with open(result) as f:
            f.readline()
            for line in f:
                l = line.rstrip('\n').split('\t')
                uid = l[0]
                prob = float(l[1])
                judge = l[2].lower()
                out_dict[uid] = {'prob':prob, 'label':judge}
    remove(f'{result_id}_LncADeep_lncRNA_results')
    return out_dict


def CPC2_runner(in_fasta, cutoff=0.5):
    '''
    :param in_fasta: transcript DNA sequence
    :type in_fast: str
    :return: out_dict: transcript_id => coding probility
    :rtype: dict
    '''
    result_id = '{}'.format(os.getpid())
    result = f'{result_id}.txt'
    cmd = 'CPC2.py -i {} -o {}'.format(in_fasta, result_id)
    ret = run_cmd(cmd, 'CPC2 done!')
    out_dict = {}
    if ret:
        with open(result) as f:
            f.readline()
            for line in f:
                l = line.rstrip('\n').split('\t')
                uid = l[0]
                prob = float(l[6])
                judge = l[7]
                out_dict[uid] = {'prob':prob, 'label':judge}
    remove(result)
    return out_dict


def CNCI2_runner(in_fasta, cutoff=0.5):
    '''
    :param in_fasta: transcript DNA sequence
    :type in_fast: str
    :return: out_dict: transcript_id => coding probility
    :rtype: dict
    '''
    result_id = '{}'.format(os.getpid())
    result = f'{result_id}/CNCI2.index'
    cmd = 'CNCI2.py -m "ve" -f {} -o {}'.format(in_fasta, result_id)
    ret = run_cmd(cmd, 'CNCI2 done!')
    out_dict = {}
    if ret:
        with open(result) as f:
            f.readline()
            for line in f:
                l = line.rstrip('\n').split('\t')
                uid = l[1].split()[0]
                judge = l[2]
                prob = float(l[3])
                out_dict[uid] = {'prob':prob, 'label':judge}
    remove(result_id)
    return out_dict


def removeIDs(ids, in_gtf, out_gtf):
    chrms = {'chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
         'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
         'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'}

    uniq_ids = {}
    for i in ids:
        uniq_ids[i.split('-')[0]] = True

    output = open(out_gtf, 'w')
    with open(in_gtf) as f:
        for line in f:
            l = line.rstrip('\n').split('\t')
            chrm = l[0]
            type = l[2]
            strand = l[6]
            anno = l[8]
            start = int(l[3])
            end = int(l[4])
            f_dict = get_value_by_key(anno,delimiter=' ')
            if ( type == "exon" or type =='transcript' ) and chrm in chrms:
                uid = f_dict['transcript_id'].strip('"')
                if not uid in uniq_ids:
                    output.write(line)
    output.close()
    print('{} generated!'.format(out_gtf))

def parse_args():
    parser = argparse.ArgumentParser(description="Coding potential prediction")
    #parser.add_argument('-i', '--input', action='store', dest='input', help="RNA assembly GTF file", required=True)
    parser.add_argument('-f', '--fasta', action='store', dest='fasta', help="RNA sequence seq (FASTA)")
    parser.add_argument('-g', '--gtf', action='store', dest='gtf', help="RNA sequence seq (GTF)")
    #parser.add_argument('-m', '--method', action='store', dest='method', help="hg19 or hg38 (default: %(default)s)", choices=['hg19','hg38'], default='hg38')
    #parser.add_argument('-o', '--output', action='store', dest='output', help="fitered GTF file for coding potential test")
    parser.add_argument('-o', '--output', action='store', dest='output', help="output table name for coding potential test")
    parser.add_argument('-r', '--ref', action='store', dest='ref', help="hg19 or hg38 (default: %(default)s)", choices=['hg19','hg38'], default='hg38')
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
        if not (options.gtf or options.fasta):
            parser.print_help()
            sys.exit('FASTA or GTF is needed!')

    if options.ref == 'hg19':
        ref_fa = 'hg19.fa'
    else:
        ref_fa = 'hg38.fa'

    external_tool_checking()

    if options.fasta and not options.gtf:
        in_fa = options.fasta
    elif options.fasta and options.gtf:
        in_fa = options.fasta
    elif options.gtf and not options.fasta:
        out_fa_name = os.path.splitext(options.gtf)[0]
        in_fa = GTF2FASTA(in_gtf=options.gtf, ref_genome=ref_fa, out_fa=f'{out_fa_name}.fasta')
    else:
        sys.exit('FASTA or GTF is needed!')
    out_table = options.output
    #fa = genLncRNASeq(args.input, ref_fa=ref_fa)
    cpat_dict = cpat_runner(in_fa)
    ucsc_dict = txCdsPredict_runner(in_fa)
    cpc2_dict = CPC2_runner(in_fa)
    lncadeep_dict = LncADeep_runner(in_fa)
    cnci2_dict = CNCI2_runner(in_fa)
    output = open(out_table, 'w')
    header = 'Transcript\tCPAT_prob\tCPAT_label\tCPC2_prob\tCPC2_label\tLncADeep_prob\tLncADeep_label\tCNCI2_prob\tCNCI2_label\tUCSC_score\tUCSC_label'
    output.write('{}\n'.format(header))
    for uid in cpat_dict:
        CPAT_prob = cpat_dict[uid]['prob']
        CPAT_label = cpat_dict[uid]['label']

        if uid in ucsc_dict:
            UCSC_score = ucsc_dict[uid]['prob']
            UCSC_label = ucsc_dict[uid]['label']
        else:
            UCSC_score = 'NA'
            UCSC_label = 'NA'

        if uid in cpc2_dict:
            CPC2_prob = cpc2_dict[uid]['prob']
            CPC2_label = cpc2_dict[uid]['label']
        else:
            CPC2_label = 'NA'
            CPC2_prob = 'NA'

        if uid in lncadeep_dict:
            LncADeep_prob = lncadeep_dict[uid]['prob']
            LncADeep_label = lncadeep_dict[uid]['label']
        else:
            LncADeep_label = 'NA'
            LncADeep_prob = 'NA'

        if uid in cnci2_dict:
            CNCI2_prob = cnci2_dict[uid]['prob']
            CNCI2_label = cnci2_dict[uid]['label']
        else:
            CNCI2_label = 'NA'
            CNCI2_prob = 'NA'
        line = f'{uid}\t{CPAT_prob}\t{CPAT_label}\t{CPC2_prob}\t{CPC2_label}\t{LncADeep_prob}\t{LncADeep_label}\t{CNCI2_prob}\t{CNCI2_label}\t{UCSC_score}\t{UCSC_label}'
        output.write(f'{line}\n')
    output.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
