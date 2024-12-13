TE-copy number estimation pipeline
================

## Generating TE reference library

1.  Download the reference genome

2.  Run EarlGrey to get the library with option “-c yes” (Docker
    version)

<!-- -->

    earlGrey -g genome.fasta -s species_name -o output/ -t 40 -c yes

3.  Manually download 3 single copy genes FASTA sequences for the
    species of interest

4.  Concatenate the SCGs to the generated repeat library

<!-- -->

    cat repeat-library.fasta scgs.fasta > repeat-library_scgs.fasta

5.  Index the merged library

<!-- -->

    bwa index repeat-library_scgs.fasta
    samtools faidx repeat-library_scgs.fasta

## Copy number estimation in short reads datasets

1.  Download the files

<!-- -->

    ena-file-download-read_run-PRJNA273563-fastq_ftp-20240704-1345.sh

2.  Merge the fastq files (\_1 and \_2) for each accession

<!-- -->

    find . -name '*_1.fastq.gz' | parallel -j4 'n={}; n=${n%_1.fastq.gz}; gzip -cd {} ${n}_2.fastq.gz | gzip -c > path/to/folder/merged/${n}.fastq.gz'

    cd path/to/folder/merged

3.  Map the short-reads to the TE reference library

<!-- -->


    find . -name "*.fastq.gz" | parallel -j 4 'n={/.}; bwa bwasw -t 10 path/to/Repeatlibrary/library.fasta {} | samtools sort -@ 4 -m 3G - > path/to/folder/mapped/${n}.sort.bam'

    cd path/to/folder/mapped

    for i in *bam;do samtools index $i;done

Before running the ara-mapstat-weight.py script, it should be modified
to contain three single copy genes of the species. At line 32 of the
script the default single copy genes for A.thaliana (AGO1, PHYB and
MPK6) can be found. These are to be replaced with the sequence names as
they are found in the Repeatlibrary corresponding to the new genes.

4.  Copy number estimation

<!-- -->

    samtools faidx path/to/Repeatlibrary/library.fasta

    find . -name '*.bam'|parallel -j 30'n={/.}
    samtools view {}|python path/to/script/ara-mapstat-weight.py --sam - --fai /path/to/Repeatlibrary/library.fasta.fai --min-mq 0 > path/to/outputfolder/${n}.deviate'

    cat *.deviate > sum_out.csv

6.  Upload the data `sum_out.csv` to the folder “data/distributions”.
    Make clear which species it is. Add some metadata to the Google
    Sheet file.

### Technical details

- The threshold for a TE family to be present/absent is 0.8 estimated
  copies.
- A copy is identified as at least 0.5 sequence coverage and 0.8
  sequence identity in long read assembly analysis.
- Only TE sequences longer than 500 bp are considered during the
  analysis

## Copy number estimation in assemblies

1.  Download long-read assemblies in fasta format

2.  Run RepeatMasker on the assemblies with the TE reference library

<!-- -->

    cd /path/to/assemblies

    for i in *.fna ; do RepeatMasker -pa 20 -no_is -s -nolow -dir path/to/output/directory -lib /path/to/library.fasta $i;done 

3.  Clean up output data

<!-- -->

    cd /path/to/RepeatMasker/output/directory

    for i in *.ori.out; do cat $i|reader-rm.py|rm-cleanup.py > $i.clean; done

    for i in *ori.out.clean; do awk '{print $0,FILENAME}' $i |perl -pe 's/\.fna\.ori\.out\.clean//'; done > merged.clean.sum

4.  Upload the data `merged.clean.sum` to the folder
    “data/distributions” specifying the species. Add some metadata to
    the Google Sheet file.
