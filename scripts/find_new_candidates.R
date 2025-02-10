suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
theme_set(theme_bw())

# Define the command-line arguments
parser <- ArgumentParser(description = "")
parser$add_argument("blast", help = "")
parser$add_argument("fai1", help = "")
parser$add_argument("fai2", help = "")
parser$add_argument("output", help = "")

# Parse the command-line arguments
args <- parser$parse_args()

# Read fai1 file
fai1 <- read_tsv(args$fai1, col_names = c("qseqid","tot_len","x","x","z")) %>% select(qseqid, tot_len) %>% separate(qseqid, into = c("acc","clus","clusn","bias","n"), sep = "_") %>%
  type_convert() %>%
  filter(bias<0.3, bias>-0.3) %>%
  mutate(qseqid = paste0(acc, "_", clus, "_", clusn, "_", bias, "_", n)) %>%
  select(qseqid,tot_len)

# Read fai2 file
fai2 <- read_tsv(args$fai2, col_names = c("sseqid","tot_len","x","x","z")) %>% select(sseqid, tot_len)

# Read the input file1
candidates <- read_tsv(args$blast, col_names = c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")) %>%
  group_by(qseqid) %>%
  filter(bitscore == max(bitscore)) %>%
  ungroup()

small_aligned <- candidates %>%
  filter(pident < 90) %>%
  inner_join(fai1, by="qseqid") %>%
  mutate(len_prop = length/tot_len) %>%
  filter(len_prop < 0.5) %>%
  select(sseqid)

unblasted <- fai2 %>% filter(!(sseqid %in% candidates$sseqid)) %>% select(sseqid)

new_candidates <- bind_rows(unblasted, small_aligned) %>% separate(sseqid, into = c("acc","clus","clusn","bias","n"), sep = "_") %>%
  type_convert() %>%
  filter(bias<0.3, bias>-0.3) %>%
  mutate(sseqid = paste0(acc, "_", clus, "_", clusn, "_", bias, "_", n)) %>%
  select(sseqid)
  
write_tsv(new_candidates, args$output, col_names = FALSE)