GD-hit filtering
================
2025-02-10

-Prerequisites: Gather Genome assemblies of your Species of interest in
one folder. Choose suitable strain for short read data. (reference
strain, recently sampled,…).

0)  Index the assemblies:

<!-- -->

    for i in assemblyfolder/*fa; do bwa index "$i"; done

1)  Run GenomeDelta:

<!-- -->

    conda activate GenomeDelta

    for fa_file in /assemblyfolder/*.fa; do
        base_name=$(basename "$fa_file" .fa)
        file="$base_name"
        GenomeDelta --fq short-reads.fastq.gz --fa "$fa_file" --of /path/to/output/"$file" --prefix outputname --t 20
    done

2)  Run the GenomeDelta iterator script to get rid of duplicate
    candidates and sequences with strong coverage bias (x \< -0.3 & x \>
    0.3).

First rename the sequences in the candidate file to ensure that each
name is unique between accessions.

    for i in ./*; do sed -i.raw 's/>/>$i_/' "$i"/*candidates.fasta ; done

Then move all the \*.candidates.fasta files into one folder and rename
the sequences

    cp */*.candidates.fasta iteration/folder

Then pick an arbitrary starting file and run the GD-iderator.sh. It is
important that the find_new_candidates.R is in the same directory as the
script.

    bash GD-iterator.sh iteration/folder starterfile.fasta output/folder

3)  BLASTx to ensure that the candidates indentified contain Transposon
    domains

<!-- -->

    blastx -query output/folder/GD-final.fa -subject transposon_domains.fasta -outfmt 6 -out TE_blast.txt

Then extract only sequences that had significant hits. The sequences of
interest are filtered in R and extracted with the python script
filter-fasta.py. It is worth mentioning that the filtering is based on
the Bitscore value produced by blastx and the cutoff at 200 can lead to
false-negatives. The database used is by no means complete but is being
regularly expanded upon. More information about the current state of the
database can be found at the end of this document.

R:

    library(tidyverse)

    hits<-read_tsv(file="/path/to/TE_blast.txt",col_names=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))%>%
    group_by(qseqid) %>%
    filter(bitscore > 200) %>%
    select(qseqid)%>%
    unique()

    write(hits,file="/path/to/working/directory/hits.txt")

Bash:

    python filter-fasta.py GD-final.fa hits.txt GD-final-TEs.fa

4)  If a reference TE library for the species of interest exists, ensure
    that the GD-candidates are not already represented there by
    conducting a blastn with said library.

<!-- -->

    blastn -query GD-final-TEs.fa -subject TE-ref-library.fasta -outfmt 6 -out ref_blast.txt

The filtering is conducted similar to step 3, however this time only
sequences WITHOUT hits are to be kept.

Bash:

    samtools faidx GD-final-TEs.fa

R:

    library(tidyverse)

    fai<-read_tsv(file="path/to/GD-final-TEs.fa.fai",col_names=c(te,length,X1,X2,X3))

    hits<-read_tsv(file="/path/to/TE_blast.txt",colnames=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"))%>%
    group_by(qseqid) %>%
    filter(bitscore > 200) %>%
    select(qseqid)%>%
    unique()

    hits_unique<-fai%>%filter(!te%in%hits)

    write(hits_unique,file="/path/to/working/directory/unique_hits.txt")

Bash:

    python filter-fasta.py GD-final-TEs.fa unique_hits.txt GD-final-TEs-unique.fa

Transposon-domain database:

-REXdb (Viridiplantae 4.0, Metazoa 3.1), Neumann et al. Systematic
survey of plant LTR-retrotransposons elucidates phylogenetic
relationships of their polyprotein domains and provides a reference for
element classification, Mobile DNA 2019
