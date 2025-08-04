# -*- coding: utf-8 -*-

## contributor: https://github.com/zhangrengang
## usage: python2 asssembly_statistic.py species.fa > statistic.txt

import sys
from Bio import SeqIO
from cutadapt import seqio
import re
import gzip
import __builtin__ as builtins  # Python 2

def open_file(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode)
    else:
        return builtins.open(filename, mode)

in_fa = sys.argv[1]
try: 
    format = sys.argv[2]
except IndexError: 
    format = 'fasta'

lst_len = []
ctg_len = []
gap_len = []
i = 0
nA, nT, nC, nG, nN = 0, 0, 0, 0, 0
try: 
    records = seqio.open(in_fa)
    tp = 'seqio'
except seqio.UnknownFileType: 
    records = SeqIO.parse(open_file(in_fa), format)
    tp = 'SeqIO'

try:
    for record in records:
        if tp == 'seqio':
            seq = str(record.sequence.upper())
        else:
            seq = str(record.seq.upper())
        lst_len.append(len(seq))
        i += 1
        for ctg in seq.split('N'):
            if ctg:
                ctg_len.append(len(ctg))

        if 'N' in seq:
            for gap in re.compile(r'N+').findall(seq):
                gap_len.append(len(gap))
                
        nA = nA + seq.count('A')
        nT = nT + seq.count('T')
        nC = nC + seq.count('C')
        nG = nG + seq.count('G')
        nN = nN + seq.count('N')
except EOFError as e: 
    print >> sys.stderr, 'WARN:', e

lst_len = sorted(lst_len, reverse=1)
gap_len = sorted(gap_len, reverse=1)
min_len = lst_len[-1]
max_len = lst_len[0]
median_len = lst_len[int(len(lst_len)/2)]
mean_len = sum(lst_len)/len(lst_len)
length = sum(lst_len)
others = length - (nG + nC + nT + nA + nN)
gc = round(100.0*(nG + nC)/length, 2)
pa, pt, pc, pg, pn, po = [round(100.0*v/length, 2) for v in [nA, nT, nC, nG, nN, others]]

sum_len = 0
l = 0
for scaf_len in lst_len:
    sum_len += scaf_len
    l += 1
    if sum_len > 0.5*length:
        N50 = scaf_len
        L50 = l
        break

sum_len = 0
l = 0
for scaf_len in lst_len:
    sum_len += scaf_len
    l += 1
    if sum_len > 0.1*length:
        N10 = scaf_len
        L10 = l
        break

sum_len = 0
l = 0
for scaf_len in lst_len:
    sum_len += scaf_len
    l += 1
    if sum_len > 0.9*length:
        N90 = scaf_len
        L90 = l
        break

def cal_NXX(len_lst, XX):
    sum_len = sum(len_lst)
    i = 0
    acc_len = 0
    for LEN in len_lst:
        acc_len += LEN
        i += 1
        if acc_len >= sum_len*XX/100.0:
            return LEN, i

print 'genome size\t%s bp\t\t\t' % (length, )
if nN > 0:
    print 'genome size without N\t%s bp\t\t\t' % (length-nN,)

print 'GC content\t%s%%\t\t\t\n  A\t%s\t%s%%\t\t\n  T\t%s\t%s%%\t\t\n  C\t%s\t%s%%\t\t\n  G\t%s\t%s%%\t\t' % (gc, nA,pa, nT,pt, nC,pc, nG,pg)
if nN > 0:
    print '  N\t%s\t%s%%\t\t' % (nN,pn)
if others > 0:
    print '  Others\t%s\t%s%%\t\t' % (others,po)

ctg_len = sorted(ctg_len, reverse=1)
cN50, cL50 = cal_NXX(ctg_len, 50)
cN10, cL10 = cal_NXX(ctg_len, 10)
cN90, cL90 = cal_NXX(ctg_len, 90)
print 'contig number\t%s\t\t\t' % len(ctg_len)
print '  Max.\t%s bp\t\tMin.\t%s bp' % (ctg_len[0], ctg_len[-1])
print '  Mean\t%d bp\t\tMedian\t%s bp' % (1.0*sum(ctg_len)/len(ctg_len), ctg_len[len(ctg_len)/2])
print '  N10\t%s bp\t\tL10\t%s\n  N50\t%s bp\t\tL50\t%s\n  N90\t%s bp\t\tL90\t%s' % (cN10,cL10,cN50,cL50,cN90,cL90)
if nN > 0 :
    print 'scaffold number\t%s\t\t\t' % i
    print '  Max.\t%s bp\t\tMin.\t%s bp' % (max_len, min_len)
    print '  Mean\t%d bp\t\tMedian\t%s bp' % (mean_len, median_len)
    print '  N10\t%s bp\t\tL10\t%s\n  N50\t%s bp\t\tL50\t%s\n  N90\t%s bp\t\tL90\t%s' % (N10,L10,N50,L50,N90,L90)
    print 'gap number\t%s\t\t\t' % len(gap_len)
    print '  Max.\t%s bp\t\tMin.\t%s bp' % (gap_len[0], gap_len[-1])
    print '  Mean\t%d bp\t\tMedian\t%s bp' % (1.0*sum(gap_len)/len(gap_len), gap_len[len(gap_len)/2] )
    