Dmel data processing
================

``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(svglite))
theme_set(theme_bw())
```

## File import

Read metadata and copynumber estimates

``` r
metadata <- read_tsv("/Volumes/Storage/mining/GenomeDelta2.0/utilities/dmel-metadata.tsv")
```

    ## Rows: 879 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): sample, strain, publication, study, study_id, location
    ## dbl (3): year, lat, lon
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
te_metadata <- read_tsv("/Volumes/Storage/ara-droso/te-metadata.tsv")
```

    ## Rows: 125 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): te, te_name, te_type, te_subtype, class
    ## dbl (1): len
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
(cn <- read_tsv("/Volumes/Storage/mining/GenomeDelta2.0/OUTPUT/dmel-validation/copynumber.tsv") %>%
    mutate(sample = gsub("/Volumes/Storage/mining/GenomeDelta2.0/OUTPUT/dmel-validation/copynumber/", "", sample),sample = gsub(".fastq.sort.tsv", "", sample)) %>%
    inner_join(metadata, by = "sample") %>%
    filter(study_id != "Shpak2023", name %in% te_metadata$te, name!=c("Souslik","Transib1_new","MLE","Spoink")))
```

    ## Rows: 223660 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): sample, type, name
    ## dbl (3): len, raw, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 101,875 × 14
    ##    sample  type  name    len    raw copynumber strain publication study study_id
    ##    <chr>   <chr> <chr> <dbl>  <dbl>      <dbl> <chr>  <chr>       <chr> <chr>   
    ##  1 ERR647… te    DME9…  7411 9.36e4       22.6 Orego… https://do… PRJN… Burny20…
    ##  2 ERR647… te    DMIS…  7439 9.80e4       23.6 Orego… https://do… PRJN… Burny20…
    ##  3 ERR647… te    DMTN…  4648 3.99e4       15.3 Orego… https://do… PRJN… Burny20…
    ##  4 ERR647… te    DMIS…  6995 1.00e5       25.6 Orego… https://do… PRJN… Burny20…
    ##  5 ERR647… te    DM23…  6126 6.56e4       19.1 Orego… https://do… PRJN… Burny20…
    ##  6 ERR647… te    412    7567 1.81e5       42.7 Orego… https://do… PRJN… Burny20…
    ##  7 ERR647… te    DMAU…  4263 2.42e4       10.2 Orego… https://do… PRJN… Burny20…
    ##  8 ERR647… te    DMBA…  1728 3.24e4       33.5 Orego… https://do… PRJN… Burny20…
    ##  9 ERR647… te    BS     5142 3.33e4       11.6 Orego… https://do… PRJN… Burny20…
    ## 10 ERR647… te    DMU8…  6411 4.87e4       13.6 Orego… https://do… PRJN… Burny20…
    ## # ℹ 101,865 more rows
    ## # ℹ 4 more variables: year <dbl>, location <chr>, lat <dbl>, lon <dbl>

``` r
(patchyness_te <- cn %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(name) %>%
    summarise(present = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(rarity = (present/815)) %>%
    arrange(rarity))
```

    ## # A tibble: 125 × 4
    ##    name         present absent rarity
    ##    <chr>          <int>  <int>  <dbl>
    ##  1 FB               622    193  0.763
    ##  2 HELITRON1_DM     704    111  0.864
    ##  3 DM14101          746     69  0.915
    ##  4 TARTC            756     59  0.928
    ##  5 AY561850         763     52  0.936
    ##  6 TRANSIB1         798     17  0.979
    ##  7 LOOPER1_DM       800     15  0.982
    ##  8 AF541951         801     14  0.983
    ##  9 TRANSIB4         803     12  0.985
    ## 10 PPI251           807      8  0.990
    ## # ℹ 115 more rows

``` r
(patchyness_sample <- cn %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(sample) %>%
    summarise(load = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(families_present = (load/129)) %>%
    inner_join(metadata, by="sample") %>%
    rename(ID=sample, Lat=lat, Long=lon) %>%
    select(ID, load, families_present, Lat, Long))
```

    ## # A tibble: 815 × 5
    ##    ID           load families_present   Lat  Long
    ##    <chr>       <int>            <dbl> <dbl> <dbl>
    ##  1 ERR6474638    124            0.961    NA    NA
    ##  2 SRR097730     123            0.953    -2    14
    ##  3 SRR098323     125            0.969     1    33
    ##  4 SRR098324     125            0.969     1    33
    ##  5 SRR098915     125            0.969     6    10
    ##  6 SRR098916     125            0.969     6    10
    ##  7 SRR105048     123            0.953     6    10
    ##  8 SRR11460801   125            0.969    NA    NA
    ##  9 SRR11846555   123            0.953    54    34
    ## 10 SRR11846560   125            0.969    52     1
    ## # ℹ 805 more rows

``` r
write_tsv(patchyness_te, "/Volumes/Storage/ara-droso/data-patchyness/D.melanogaster.TEs.tsv")
write_tsv(patchyness_sample, "/Volumes/Storage/ara-droso/data-patchyness/D.melanogaster.tsv")
```

    # Define the folder containing the copynumber files
    folder_path <- "/Volumes/Storage/ara-droso/data-patchyness/dyak"

    columns <- c("type","name","len","raw","copynumber")

    # List all *.ori.out files in the folder
    file_list <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)

    # Iterate over each file in oriout_folder
    all_cn_data <- tibble(sample = character(), type = character(), name = character(), len = integer(), raw = integer(), copynumber = double())

    for (cn_file in file_list) {
      cn <- read_tsv(cn_file, col_names = columns, skip = 6) %>% mutate(sample = cn_file)
      print(cn_file)
      all_cn_data <- bind_rows(all_cn_data, cn)
      }
      
    write_tsv(all_cn_data, "/Volumes/Storage/ara-droso/data-patchyness/dyak/copynumber.tsv")

``` r
(dsim <- read_tsv("/Volumes/Storage/ara-droso/data-patchyness/dsim/copynumber.tsv") %>%
    mutate(sample = gsub("/Volumes/Storage/ara-droso/data-patchyness/dsim/", "", sample),sample = gsub(".fastq.sort.tsv", "", sample)) %>% filter(type=="te", !str_detect(name, "Unknown|Satellite|Simple_repeat|tRNA"), !(sample%in%c("ERR694696","ERR694699","ERR694698","ERR694696","ERR694700"))))
```

    ## Rows: 19158 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): sample, type, name
    ## dbl (3): len, raw, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 13,706 × 6
    ##    sample    type  name                             len     raw copynumber
    ##    <chr>     <chr> <chr>                          <dbl>   <dbl>      <dbl>
    ##  1 ERR668306 te    rnd-1_family-13#LINE/CR1        4494 232643.      69.1 
    ##  2 ERR668306 te    rnd-1_family-133#LTR/Gypsy      6334  40737.       8.59
    ##  3 ERR668306 te    rnd-1_family-143#LINE/I-Jockey  3575  35241.      13.2 
    ##  4 ERR668306 te    rnd-1_family-159#LTR/Gypsy      8025  16246.       2.70
    ##  5 ERR668306 te    rnd-1_family-84#LINE/I          5382  75620.      18.8 
    ##  6 ERR668306 te    rnd-1_family-1#LINE/I-Jockey    3369  19543.       7.75
    ##  7 ERR668306 te    rnd-1_family-104#RC/Helitron     135    277.       2.74
    ##  8 ERR668306 te    rnd-1_family-113#LINE/I-Jockey  4275  41647.      13.0 
    ##  9 ERR668306 te    rnd-1_family-16#LINE/I-Jockey   4724  48682.      13.8 
    ## 10 ERR668306 te    rnd-1_family-162#LTR/Gypsy        75    272        4.84
    ## # ℹ 13,696 more rows

``` r
(patchyness_te_sim <- dsim %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(name) %>%
    summarise(present = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(rarity = (present/89)) %>%
    arrange(rarity))
```

    ## # A tibble: 154 × 4
    ##    name                           present absent rarity
    ##    <chr>                            <int>  <int>  <dbl>
    ##  1 rnd-4_family-1044#LTR/Gypsy          7     82 0.0787
    ##  2 rnd-3_family-41#LINE                32     57 0.360 
    ##  3 rnd-4_family-4712#DNA/P             39     50 0.438 
    ##  4 rnd-4_family-1031#LTR/Pao           84      5 0.944 
    ##  5 rnd-1_family-154#LTR/Gypsy          86      3 0.966 
    ##  6 rnd-4_family-1931#DNA/P             86      3 0.966 
    ##  7 rnd-4_family-805#LINE/I-Jockey      86      3 0.966 
    ##  8 rnd-1_family-148#DNA/TcMar-Tc1      87      2 0.978 
    ##  9 rnd-1_family-23#LINE/I-Jockey       87      2 0.978 
    ## 10 rnd-1_family-52#RC/Helitron         87      2 0.978 
    ## # ℹ 144 more rows

``` r
(patchyness_sample <- dsim %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(sample) %>%
    summarise(load = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(families_present = (load/154)) %>%
    rename(ID=sample) %>%
    select(ID, load, families_present))
```

    ## # A tibble: 89 × 3
    ##    ID         load families_present
    ##    <chr>     <int>            <dbl>
    ##  1 ERR668306   151            0.981
    ##  2 ERR668307   151            0.981
    ##  3 ERR668308   151            0.981
    ##  4 ERR668309   151            0.981
    ##  5 ERR668310   151            0.981
    ##  6 ERR668311   151            0.981
    ##  7 ERR668312   151            0.981
    ##  8 ERR668313   151            0.981
    ##  9 ERR668314   151            0.981
    ## 10 ERR668315   150            0.974
    ## # ℹ 79 more rows

``` r
write_tsv(patchyness_sample, "/Volumes/Storage/ara-droso/data-patchyness/D.simulans.tsv")
write_tsv(patchyness_te_sim, "/Volumes/Storage/ara-droso/data-patchyness/D.simulans.TEs.tsv")
```

``` r
(dsec <- read_tsv("/Volumes/Storage/ara-droso/data-patchyness/dsec/copynumber.tsv") %>%
    mutate(sample = gsub("/Volumes/Storage/ara-droso/data-patchyness/dsec/", "", sample),sample = gsub(".fastq.sort.tsv", "", sample)) %>% filter(type=="te", !str_detect(name, "Unknown|Satellite|Simple_repeat|tRNA"), name!="rnd-4_family-4968#DNA/P"))
```

    ## Rows: 9630 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): sample, type, name
    ## dbl (3): len, raw, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 7,200 × 6
    ##    sample      type  name                             len     raw copynumber
    ##    <chr>       <chr> <chr>                          <dbl>   <dbl>      <dbl>
    ##  1 SRR14138506 te    rnd-1_family-50#LINE/I-Jockey   4740  44416.      29.5 
    ##  2 SRR14138506 te    rnd-1_family-124#LINE/R1        5946 262667.     139.  
    ##  3 SRR14138506 te    rnd-1_family-136#LINE/I-Jockey  4424  23635.      16.8 
    ##  4 SRR14138506 te    rnd-1_family-17#LINE/I-Jockey   4712  87624.      58.6 
    ##  5 SRR14138506 te    rnd-1_family-25#LINE/I-Jockey   4729  13743.       9.16
    ##  6 SRR14138506 te    rnd-1_family-30#LINE/I-Jockey   4748  59144.      39.3 
    ##  7 SRR14138506 te    rnd-1_family-77#LINE/R1-LOA     5452  82805.      47.9 
    ##  8 SRR14138506 te    rnd-1_family-120#LTR/Gypsy        86    785.      28.8 
    ##  9 SRR14138506 te    rnd-1_family-154#LINE/R1-LOA     643   3086.      15.1 
    ## 10 SRR14138506 te    rnd-1_family-21#DNA/P           1159  16579.      45.1 
    ## # ℹ 7,190 more rows

``` r
(patchyness_te_sec <- dsec %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(name) %>%
    summarise(present = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(rarity = (present/45)) %>%
    arrange(rarity))
```

    ## # A tibble: 160 × 4
    ##    name                           present absent rarity
    ##    <chr>                            <int>  <int>  <dbl>
    ##  1 rnd-4_family-2363#LTR/Gypsy         27     18  0.6  
    ##  2 rnd-4_family-1504#DNA/P             44      1  0.978
    ##  3 rnd-1_family-102#LTR/Gypsy          45      0  1    
    ##  4 rnd-1_family-103#LTR/Pao            45      0  1    
    ##  5 rnd-1_family-110#LINE/I-Jockey      45      0  1    
    ##  6 rnd-1_family-113#LTR/Gypsy          45      0  1    
    ##  7 rnd-1_family-117#LTR/Gypsy          45      0  1    
    ##  8 rnd-1_family-118#LINE/I-Jockey      45      0  1    
    ##  9 rnd-1_family-120#LTR/Gypsy          45      0  1    
    ## 10 rnd-1_family-124#LINE/R1            45      0  1    
    ## # ℹ 150 more rows

``` r
(patchyness_sample <- dsec %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(sample) %>%
    summarise(load = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(families_present = (load/160)) %>%
    rename(ID=sample) %>%
    select(ID, load, families_present))
```

    ## # A tibble: 45 × 3
    ##    ID           load families_present
    ##    <chr>       <int>            <dbl>
    ##  1 SRR14138506   159            0.994
    ##  2 SRR14138507   159            0.994
    ##  3 SRR5860570    159            0.994
    ##  4 SRR5860573    159            0.994
    ##  5 SRR5860582    159            0.994
    ##  6 SRR5860583    160            1    
    ##  7 SRR5860584    159            0.994
    ##  8 SRR5860625    160            1    
    ##  9 SRR5860626    160            1    
    ## 10 SRR5860627    160            1    
    ## # ℹ 35 more rows

``` r
write_tsv(patchyness_sample, "/Volumes/Storage/ara-droso/data-patchyness/D.sechellia.tsv")
write_tsv(patchyness_te_sec, "/Volumes/Storage/ara-droso/data-patchyness/D.sechellia.TEs.tsv")
```

``` r
(dyak <- read_tsv("/Volumes/Storage/ara-droso/data-patchyness/dyak/copynumber.tsv") %>%
    mutate(sample = gsub("/Volumes/Storage/ara-droso/data-patchyness/dyak/", "", sample),sample = gsub(".fastq.sort.tsv", "", sample)) %>% filter(type=="te", !str_detect(name, "Unknown|Satellite|Simple_repeat|tRNA"), !(name %in% c("rnd-1_family-314#LTR/Gypsy","rnd-1_family-316#LTR/Gypsy","rnd-1_family-96#LTR/Gypsy","rnd-1_family-226#LTR/Gypsy"))))
```

    ## Rows: 16131 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): sample, type, name
    ## dbl (3): len, raw, copynumber
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 11,799 × 6
    ##    sample       type  name                                 len    raw copynumber
    ##    <chr>        <chr> <chr>                              <dbl>  <dbl>      <dbl>
    ##  1 SAMN04044077 te    rnd-1_family-79#DNA/TcMar-Tc1       1691  8516.       8.68
    ##  2 SAMN04044077 te    rnd-1_family-17#LINE/CR1            4419 47422.      18.5 
    ##  3 SAMN04044077 te    rnd-4_family-1317#LINE/I-Jockey     3471 15289.       7.59
    ##  4 SAMN04044077 te    rnd-1_family-150#LINE/I-Jockey      4695 31536.      11.6 
    ##  5 SAMN04044077 te    rnd-1_family-284#LINE/I-Jockey      4169  2878.       1.19
    ##  6 SAMN04044077 te    rnd-1_family-327#LINE/I-Jockey      4855  3769.       1.34
    ##  7 SAMN04044077 te    rnd-1_family-95#LINE/R1             5387 26307.       8.42
    ##  8 SAMN04044077 te    rnd-1_family-11#RC/Helitron          264  5967.      39.0 
    ##  9 SAMN04044077 te    rnd-1_family-135#DNA/TcMar-Mariner  1308 15185.      20.0 
    ## 10 SAMN04044077 te    rnd-1_family-181#DNA/CMC-Transib    1148  1527.       2.29
    ## # ℹ 11,789 more rows

``` r
(patchyness_te_yak <- dyak %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(name) %>%
    summarise(present = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(rarity = (present/57)) %>%
    arrange(rarity))
```

    ## # A tibble: 207 × 4
    ##    name                          present absent rarity
    ##    <chr>                           <int>  <int>  <dbl>
    ##  1 rnd-1_family-60#RC/Helitron        42     15  0.737
    ##  2 rnd-4_family-1072#LTR/Gypsy        42     15  0.737
    ##  3 rnd-4_family-4385#LTR/Gypsy        48      9  0.842
    ##  4 rnd-1_family-252#LTR/Gypsy         53      4  0.930
    ##  5 rnd-4_family-1380#RC/Helitron      55      2  0.965
    ##  6 rnd-4_family-76#LTR/Copia          55      2  0.965
    ##  7 rnd-4_family-74#LTR/Copia          56      1  0.982
    ##  8 rnd-1_family-102#LTR/Gypsy         57      0  1    
    ##  9 rnd-1_family-106#LTR/Gypsy         57      0  1    
    ## 10 rnd-1_family-109#LTR/Pao           57      0  1    
    ## # ℹ 197 more rows

``` r
(patchyness_sample <- dyak %>%
    mutate(presence = ifelse(copynumber > 0.3, "present", "absent")) %>%
    group_by(sample) %>%
    summarise(load = sum(presence == "present"), absent = sum(presence == "absent")) %>%
    mutate(families_present = (load/207)) %>%
    rename(ID=sample) %>%
    select(ID, load, families_present))
```

    ## # A tibble: 57 × 3
    ##    ID            load families_present
    ##    <chr>        <int>            <dbl>
    ##  1 SAMN04044077   203            0.981
    ##  2 SAMN04044078   202            0.976
    ##  3 SRR3739672     205            0.990
    ##  4 SRR3739673     204            0.986
    ##  5 SRR3739674     205            0.990
    ##  6 SRR3739675     205            0.990
    ##  7 SRR3739676     205            0.990
    ##  8 SRR3739677     206            0.995
    ##  9 SRR3739678     206            0.995
    ## 10 SRR3739679     205            0.990
    ## # ℹ 47 more rows

``` r
write_tsv(patchyness_sample, "/Volumes/Storage/ara-droso/data-patchyness/D.yakuba.tsv")
write_tsv(patchyness_te_yak, "/Volumes/Storage/ara-droso/data-patchyness/D.yakuba.TEs.tsv")
```

``` r
# Define the folder containing the *.ori.out files
folder_path <- "/Volumes/Storage/ara-droso/data-longreads/dmel"

columns <- c("SWscore", "pdiv", "pdel", "pins", "contig", "qstart", "qend", "qleft",  "strand", "te", "class", "position_in_repeat_begin", "position_in_repeat_end", "position_in_repeat_left", "ID")

# List all *.ori.out files in the folder
file_list <- list.files(path = folder_path, pattern = "\\.ori\\.out$", full.names = TRUE)

# Iterate over each file in oriout_folder
max_data <- tibble(Species = character(), te = character(), contig = character(), qstart = double(), qend = double(), strand = character(), SWscore = double(), pdiv = double())

for (oriout_file in file_list) {
  (oo <- read_table(oriout_file, col_names = columns) %>% mutate(Species = gsub(".ori.out", "", oriout_file), Species = gsub("/Volumes/Storage/ara-droso/data-longreads/dmel//", "", Species)) %>% select(Species, te, contig, qstart, qend, strand, SWscore, pdiv))
  
  max_data <- bind_rows(max_data, oo)
  }
```

    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )
    ## 
    ## 
    ## ── Column specification ────────────────────────────────────────────────────────
    ## cols(
    ##   SWscore = col_double(),
    ##   pdiv = col_double(),
    ##   pdel = col_double(),
    ##   pins = col_double(),
    ##   contig = col_character(),
    ##   qstart = col_double(),
    ##   qend = col_double(),
    ##   qleft = col_character(),
    ##   strand = col_character(),
    ##   te = col_character(),
    ##   class = col_character(),
    ##   position_in_repeat_begin = col_character(),
    ##   position_in_repeat_end = col_double(),
    ##   position_in_repeat_left = col_character(),
    ##   ID = col_logical()
    ## )

``` r
write_tsv(max_data, "/Volumes/Storage/ara-droso/data-longreads/dmel/DmelRM.tsv")
```

``` r
(dmel_tes <- read_tsv("/Volumes/Storage/ara-droso/te-metadata.tsv"))
```

    ## Rows: 125 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): te, te_name, te_type, te_subtype, class
    ## dbl (1): len
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 125 × 6
    ##    te       te_name    te_type                          te_subtype  class    len
    ##    <chr>    <chr>      <chr>                            <chr>       <chr>  <dbl>
    ##  1 TC3      Tc3        terminal_inverted_repeat_element Tc1-Mariner DNA     1743
    ##  2 DMTN1731 1731       LTR_retrotransposon              Copia       LTR     4648
    ##  3 DMMDG3   mdg3       LTR_retrotransposon              Gypsy       LTR     5519
    ##  4 DMRTMGD1 mdg1       LTR_retrotransposon              Gypsy       LTR     7480
    ##  5 DMGYPF1A gypsy      LTR_retrotransposon              Gypsy       LTR     7469
    ##  6 DMCOPIA  copia      LTR_retrotransposon              Copia       LTR     5143
    ##  7 DMRER2DM R2-element non_LTR_retrotransposon          R2          nonLTR  3607
    ##  8 DMBARI1  Bari1      terminal_inverted_repeat_element Tc1-Mariner DNA     1728
    ##  9 DMLINEJA jockey     non_LTR_retrotransposon          Jockey      nonLTR  5020
    ## 10 DMIS176  17.6       LTR_retrotransposon              Gypsy       LTR     7439
    ## # ℹ 115 more rows

``` r
(dmel_LR <- read_tsv("/Volumes/Storage/ara-droso/data-longreads/dmel/DmelRM.tsv") %>% mutate(Species = gsub("/Volumes/Storage/ara-droso/data-longreads/dmel/", "", Species)) %>% group_by(te, Species) %>% filter(SWscore==max(SWscore)) %>% ungroup() %>% group_by(te) %>% mutate(score = SWscore/max(SWscore)) %>% filter(te %in% dmel_tes$te) %>% mutate(presence_absence = ifelse(score>0.5, "present", "absent")) %>% select(Species, te, presence_absence) %>% distinct() %>% group_by(te) %>% mutate(present_in = sum(presence_absence=="present")) %>% select(Species, te, present_in))
```

    ## Rows: 2203886 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): Species, te, contig, strand
    ## dbl (4): qstart, qend, SWscore, pdiv
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 6,112 × 3
    ## # Groups:   te [125]
    ##    Species     te         present_in
    ##    <chr>       <chr>           <int>
    ##  1 D.mel.A2.fa DMTNFB             49
    ##  2 D.mel.A2.fa 1360               49
    ##  3 D.mel.A2.fa 412                49
    ##  4 D.mel.A2.fa DMMDG3             49
    ##  5 D.mel.A2.fa FB                 13
    ##  6 D.mel.A2.fa DMTRDNA            49
    ##  7 D.mel.A2.fa STALKER3           49
    ##  8 D.mel.A2.fa McCLINTOCK         49
    ##  9 D.mel.A2.fa INVADER5           48
    ## 10 D.mel.A2.fa TRANSIB3           49
    ## # ℹ 6,102 more rows

``` r
write_tsv(dmel_LR, "/Volumes/Storage/ara-droso/data-longreads/DmelRM.tsv")
```
