#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import random
import math
import collections
import fileinput

def get_weight(cigar,readlen):
    pe=[]
    for fi in re.finditer(r"(\d+)([HSIDMN])", cigar):
        num=int(fi.group(1))
        id=fi.group(2)
        pe.append((num,id))
    matchsum=0
    for num,id in pe:
        if id=="M":
            matchsum+=num
        elif id=="D" or id=="N" or id=="I" or id=="S" or id=="H":
            pass
        else:
            raise Exception("unknown cigar"+id)
    return float(matchsum)/float(readlen)

def readfai(fai):
     toret={}
     telist=[]
     scglist=[]
     scgset=set(["AGO1","MPK6","PHYB"])
    
     for l in open(fai):
          # LTR65_te	669
          a=l.rstrip("\n").split("\t")
          leng=int(a[1])
          seqid=a[0]
          if seqid in scgset:
               scglist.append(seqid)
          else:
               telist.append(seqid)
          toret[seqid]=leng
     return(scglist,telist,toret)
          



parser = argparse.ArgumentParser(description="""           
Description
-----------
Summary statistics
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""

Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--min-mq", type=int, required=False, dest="minmq", default=0, help="min mapping quality")
parser.add_argument("--fai", type=str, required=True, dest="fai", default=None, help="fai fasta-index (samtools faidx) of the TE database")
args = parser.parse_args()
minmq=args.minmq

scglist,telist,lengthash=readfai(args.fai)
scgset=set(scglist)
teset=set(telist)

hte=collections.defaultdict(lambda:0.0)
hscg=collections.defaultdict(lambda:0.0)


sumall=0
summapped=0
summapq=0
sumweight=0.0

sumte=0.0
sumscg=0.0


for line in args.sam:
     """
0     1         2    3    4    5   6    7    8            9                        10                  11
r1	16	M14653_te	172	70	23M	*	0	0	ATGTCGAGTTTCGTGCCGAATAA	FFFFFFFFFFFFFFFFFFBBBBB	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:23
r2	0	M14653_te	240	70	27M	*	0	0	AACAGCTGCGGAATCGCACCGAATGCT	BBBBBFFFFFBFFFFFFFFFFFFFFFF	PG:Z:novoalign	AS:i:0	UQ:i:0	NM:i:0	MD:Z:27
     """
     if line.startswith("@"):
          continue
     line=line.rstrip("\n")
     a=line.split("\t")
     flag=int(a[1])
     
     sumall+=1
     
     if flag & 0x004 > 0:   # remove unmapped
          continue
     summapped+=1 # all mapped reads

     mq=int(a[4])
     if mq< minmq:          # remove with minmapq
          continue
     summapq+=1 # all  reads with minmapq
     
     weight=get_weight(a[5],len(a[9])) # compute weight of the read; 
     sumweight+=weight
     
     ref=a[2]
     if ref in scgset:
          scgseq=ref
          hscg[scgseq]+=weight
          sumscg+=weight
    
     elif ref in teset:
          teseq=ref
          hte[teseq]+=weight
          sumte+=weight
     else:
          raise Warning("unknown reference; ignoring "+ref)


meanscgcov=0.0
covscg=[]
for scg,count in hscg.items():
     length=lengthash[scg]
     coverage=float(count)/float(length)
     covscg.append(coverage)
     meanscgcov+=coverage
meanscgcov=meanscgcov/float(len(covscg))
     
 
print("{0}\t{1}\t{2}".format("summary","all_reads",sumall))
print("{0}\t{1}\t{2}".format("summary","mapped_reads",summapped))
print("{0}\t{1}\t{2}".format("summary","reads_with_mapq",summapq))
print("{0}\t{1}\t{2}".format("summary","weighted_reads_with_mapq",sumweight))
print("{0}\t{1}\t{2}".format("summary","mapping_to_te_weighted",sumte))
print("{0}\t{1}\t{2}".format("summary","mapping_to_scg_weighted",sumscg))


covscg=[]
for scg in scglist:
     count=hscg[scg]
     length=lengthash[scg]
     coverage=float(count)/float(length)
     normcov=coverage/meanscgcov
     print("{0}\t{1}\t{2}\t{3}\t{4}".format("scg",scg,length,count,normcov))

for te in telist:
     length=lengthash[te]
     count=hte[te]
     coverage=float(count)/float(length)
     normcov=coverage/meanscgcov
     print("{0}\t{1}\t{2}\t{3}\t{4}".format("te",te,length,count,normcov))