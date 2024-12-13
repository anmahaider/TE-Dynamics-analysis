TE-copy number estimation pipeline
================
2024-12-06

-\>Acquisition and preparation of short-read dataset

    ena-file-download-read_run-PRJNA273563-fastq_ftp-20240704-1345.sh

-\>merging the fastq files for each accession


    find . -name '*_1.fastq.gz' | parallel -j4 'n={}; n=${n%_1.fastq.gz}; gzip -cd {} ${n}_2.fastq.gz | gzip -c > path/to/folder/merged/${n}.fastq.gz'

    cd path/to/folder/merged

-\>mapping the short-reads to the merged Repeatlibrary


    find . -name "*.fastq.gz" | parallel -j 4 'n={/.}; bwa bwasw -t 10 path/to/Repeatlibrary/library.fasta {} | samtools sort -@ 4 -m 3G - > path/to/folder/mapped/${n}.sort.bam'

    cd path/to/folder/mapped

    for i in *bam;do samtools index $i;done

Before running the ara-mapstat-weight.py script, it should be modifyed
to contain three single copy genes of the species. At line 32 of the
script the default single copy genes for A.thaliana AGO1, PHYB and MPK6
can be found. These are to be replaced with the sequence names as they
are found in the Repeatlibrary corresponding to the new genes.

-\>copy number estimation

    samtools faidx path/to/Repeatlibrary/library.fasta

    find . -name '*.bam'|parallel -j 30'n={/.}
    samtools view {}|python path/to/script/ara-mapstat-weight.py --sam - --fai /path/to/Repeatlibrary/library.fasta.fai --min-mq 0 > path/to/outputfolder/${n}.deviate'

    cat *.deviate > sum_out.csv

If the Transposon library in use does not contain the consensus sequence
of each family but rather entries for each identified insertion, it is
possible to collapse the copy number estimation for each of those
insertions per family. For example the copy number of a TE, here name
Transposon1 in an accession would be calculated as followed:

Output: <Transposon1@Chr1>:500000-504000 norm_count=1
<Transposon1@Chr4>:700000-704000 norm_count=3

Collapse of entries for Transposon1: Transposon1 norm_count=4
