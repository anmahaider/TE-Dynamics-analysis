#!/usr/bin/env python
import os
import sys
import inspect
import re
import argparse
import collections




parser = argparse.ArgumentParser(description=
"""           
Cleans repeatmaskers weird format into something more accessible. Notable unifies notation for reverse complements.
Output example
rm	1995	9.31	CM034928.1	149195	149570	+	DNAREP1_INE-1	1	424	0.62	D.mel.AKA017
rm	1512	13.85	CM034928.1	150460	150800	C	DNAREP1_INE-1	1	377	0.55	D.mel.AKA017
rm	483	21.21	CM034928.1	151847	152023	C	DNAREP1_INE-1	1	195	0.28	D.mel.AKA017
rm	2111	6.98	CM034928.1	154650	155010	+	DNAREP1_INE-1	1	401	0.59	D.mel.AKA017

Columns
col1: bionicle object ID
col2: score
col3: substitutions
col4: chromsome of repeat match
col5: start of repeat match
col6: end of repeat match
col7: strand of TE, + or C, of repeat match
col8: repeat
col9: start in repeat
col10: end in repeat
col11: aligned fraction of repeat (eg 1 is the entire repeat is aligned)
col12: optional, e.g. filename


""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
"""
Documentation from the RepeatMasker webpage
https://www.repeatmasker.org/webrepeatmaskerhelp.html

How to read the results
The annotation file contains the cross_match output lines. It lists all best matches (above a set minimum score) between the query sequence and any of the sequences in the repeat database or with low complexity DNA. The term "best matches" reflects that a match is not shown if its domain is over 80% contained within the domain of a higher scoring match, where the "domain" of a match is the region in the query sequence that is defined by the alignment start and stop. These domains have been masked in the returned masked sequence file. In the output, matches are ordered by query name, and for each query by position of the start of the alignment.

Example:
 1306 15.6  6.2  0.0 HSU08988  6563  6781  (22462) C  MER7A    DNA/MER2_type    (0)   336   103
12204 10.0  2.4  1.8 HSU08988  6782  7714  (21529) C  TIGGER1  DNA/MER2_type    (0)  2418  1493
  279  3.0  0.0  0.0 HSU08988  7719  7751  (21492) +  (TTTTA)n Simple_repeat      1    33   (0)
 1765 13.4  6.5  1.8 HSU08988  7752  8022  (21221) C  AluSx    SINE/Alu        (23)   289     1
12204 10.0  2.4  1.8 HSU08988  8023  8694  (20549) C  TIGGER1  DNA/MER2_type  (925)  1493   827
 1984 11.1  0.3  0.7 HSU08988  8695  9000  (20243) C  AluSg    SINE/Alu         (5)   305     1
12204 10.0  2.4  1.8 HSU08988  9001  9695  (19548) C  TIGGER1  DNA/MER2_type (1591)   827     2
  711 21.2  1.4  0.0 HSU08988  9696  9816  (19427) C  MER7A    DNA/MER2_type  (224)   122     2
This is a sequence in which a Tigger1 DNA transposon has integrated into a MER7 DNA transposon copy. Subsequently two Alus integrated in the Tigger1 sequence. The simple repeat is derived from the poly A of the Alu element. The first line is interpreted like this:

  1306    = Smith-Waterman score of the match, usually complexity adjusted
        The SW scores are not always directly comparable. Sometimes
        the complexity adjustment has been turned off, and a variety of
        scoring-matrices are used.
  15.6    = % substitutions in matching region compared to the consensus
  6.2     = % of bases opposite a gap in the query sequence (deleted bp)
  0.0     = % of bases opposite a gap in the repeat consensus (inserted bp)
  HSU08988 = name of query sequence
  6563    = starting position of match in query sequence
  7714    = ending position of match in query sequence
  (22462) = no. of bases in query sequence past the ending position of match
  C       = match is with the Complement of the consensus sequence in the database
  MER7A   = name of the matching interspersed repeat
  DNA/MER2_type = the class of the repeat, in this case a DNA transposon 
            fossil of the MER2 group (see below for list and references)
  (0)     = no. of bases in (complement of) the repeat consensus sequence 
            prior to beginning of the match (so 0 means that the match extended 
            all the way to the end of the repeat consensus sequence)
  2418    = starting position of match in database sequence (using top-strand numbering)
  1465    = ending position of match in database sequence

An asterisk (*) in the final column (no example shown) indicates that there is a higher-scoring match whose domain partly (<80%) includes the domain of this match.

Note that the SW score and divergence numbers for the three Tigger1 lines are identical. This is because the information is derived from a single alignment (the Alus were deleted from the query before the alignment with the Tigger element was performed). The program makes educated guesses about many fragments if they are derived from the same element (e.g. it knows that the MER7A fragments represent one insert). In a next version I can identify each element with a unique ID, if interest exists (this could help to represent repeats cleaner in graphic displays).
"""

"""
0	1		2		3		4		5			6		7		8			9	10				11				12		13		14		15	
ob	score	subst	del		ins		query		qstart	qend	qpost		rc	repeat			class			start	end		post	file
rm	822		7.48	12.29	10.44	CM034928.1	157688	157866	(26121625)	+	DNAREP1_INE-1	Unspecified		1		182		(501)	D.mel.AKA017
rm	309		10.58	5.68	3.33	CM034928.1	157865	157952	(26121539)	+	DNAREP1_INE-1	Unspecified		296		385		(298)	D.mel.AKA017
rm	1745	13.14	5.42	4.39	CM034928.1	174637	175042	(26104449)	C	DNAREP1_INE-1	Unspecified		(273)	410		1		D.mel.AKA017
rm	385		21.21	8.08	0.00	CM034928.1	175542	175640	(26103851)	C	DNAREP1_INE-1	Unspecified		(461)	222		116		D.mel.AKA017
"""


parser.add_argument("--someflag", type=str, required=False, dest="bed", default="" ,help="not yet implemented")
args = parser.parse_args()

for line in sys.stdin:
	
	a=line.rstrip("\n").split("\t")
	# object check
	if a[0] != "rm":
		raise Exception("wrong bionicle object; only accept 'rm', got "+a[0])
	rc=a[9]
	rms,rme,rmp=a[12],a[13],a[14]
	rstart,rend,post=0,0,0
	if(rc=="+"):
		tmp=int(rmp.strip("()"))
		rstart,rend,post=int(rms),int(rme),int(tmp)
	elif(rc=="C"):
		tmp=int(rms.strip("()"))
		rstart,rend,post=int(rmp),int(rme),int(tmp)
		pass
	else:
		raise Exception("Invalid rc; must be + or C")
	rlenght,aligned=0,0
	rlength=rend+post
	aligned="{:.2f}".format(float(rend-rstart)/float(rlength))
	topr=[a[0],a[1],a[2],a[5],a[6],a[7],a[9],a[10],str(rstart),str(rend),aligned]
	if(len(a)==16):
		topr.append(a[15])
	print("\t".join(topr))
	
