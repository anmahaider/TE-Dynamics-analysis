#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse


parser = argparse.ArgumentParser(description="""           
Reads RepeatMaskers ori.out files into bioinicle objects
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
#parser.add_argument('--fastq', type=str, default=None,dest="fasta", required=True, help="A fasta file")

args,unknown = parser.parse_known_args()
fh=None
if(len(unknown)>0):
	fh=open(unknown[0])
else:
	fh=sys.stdin

while(True):
	h1=fh.readline()
	if h1 =="":
		break
	striped=h1.lstrip().rstrip("\n")
	tmp=re.split("\s+",striped)
	topr=["rm"]
	topr.extend(tmp)
	print("\t".join(topr))