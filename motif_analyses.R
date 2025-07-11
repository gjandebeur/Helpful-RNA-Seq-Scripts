### investigation of motif for each sample/species

### dependencies ###
library(stringr)
library(dplyr)
library(ggplot2)

### prevents scientific notation in graph
options(scipen=999)

### producing updated bedfiles 
### change code from 17802 (pseU) to a (m6A) if preferred ###

native_filtered %>%
  filter(code == "17802") %>%
  select(chrom, start, end, code, n_valid_cov, strand) %>%
  distinct(chrom, start, end, .keep_all = TRUE) %>%
  arrange(chrom, start) %>%
  mutate(
    start = start - 3,
    end = end + 3
  ) %>%
  write.table(
    "./for_motif/at_pseU_3pad.bed",
    col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
  )

#### this is IVT motif samples that are NOT found in the native replicates

ivt_antifiltered %>%
  filter(code == "17802") %>%
  select(chrom, start, end, code, n_valid_cov, strand) %>%
  distinct(chrom, start, end, .keep_all = TRUE) %>%
  arrange(chrom, start) %>%
  mutate(
    start = start - 3,
    end = end + 3
  ) %>%
  write.table(
    "./for_motif/at_pseU_IVT_3pad.bed",
    col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
  )

### FOLLOWING THIS IS 1 SMALL IMPORTANT STEP. Can use MEME for motif analyses OR 
### Use a UNIX SSH to run bedtools intersect with the updated fasta reference file to generate bed2fasta files (.fasta) 


# after producing fasta files in a seperate SSH that is UNIX rather than R.
#read the files created, change pseU for m6A if needed

lines <- readLines("./for_motif/at_pseU_3pad.fa")
IVT_lines <- readLines("./for_motif/at_pseU_IVT_3pad.fa")

### small filtering for the files above

stopifnot(length(lines) %% 2 == 0)
stopifnot(length(IVT_lines) %% 2 == 0)

headers <- lines[seq(1, length(lines), by = 2)]
motifs  <- lines[seq(2, length(lines), by = 2)]

IVT_headers <- IVT_lines[seq(1, length(IVT_lines), by = 2)]
IVT_motifs  <- IVT_lines[seq(2, length(IVT_lines), by = 2)]

parsed <- str_match(headers, "^>([^:]+):(\\d+)-(\\d+)\\(([+-])\\)")
IVT_parsed <- str_match(IVT_headers, "^>([^:]+):(\\d+)-(\\d+)\\(([+-])\\)")

#### now transform into df for graphing.

at_pseU_3pad <- data.frame(
  chrom  = paste0("chr", parsed[, 2]),
  start  = as.integer(parsed[, 3]),
  end    = as.integer(parsed[, 4]),
  strand = parsed[, 5],
  motif  = motifs,
  stringsAsFactors = FALSE
)

at_pseU_IVT_3pad <- data.frame(
  chrom  = paste0("chr", IVT_parsed[, 2]),
  start  = as.integer(IVT_parsed[, 3]),
  end    = as.integer(IVT_parsed[, 4]),
  strand = IVT_parsed[, 5],
  motif  = IVT_motifs,
  stringsAsFactors = FALSE
)

#### collapse down to 5mer (uses 7mer at first for safety, but ONT pores are 5mers)
at_pseU_2pad <- at_pseU_3pad %>%
  mutate(
    start = start + 1,
    end = end - 1,
    motif = substr(motif, 2, nchar(motif) - 1)
  )

at_pseU_IVT_2pad <- at_pseU_IVT_3pad %>%
  mutate(
    start = start + 1,
    end = end - 1,
    motif = substr(motif, 2, nchar(motif) - 1)
  )

### view which motifs are shared between the two (quality control step)

shared_motifs <- intersect(
  unique(at_pseU_2pad$motif),
  unique(at_pseU_IVT_2pad$motif)
)

at_pseU_2pad_shared <- at_pseU_2pad %>%
  filter(motif %in% shared_motifs)

at_pseU_IVT_2pad_shared <- at_pseU_IVT_2pad %>%
  filter(motif %in% shared_motifs)

### count unique motifs and their percent of total motif

motif_ranked_top25 <- at_pseU_2pad %>%
  group_by(motif) %>%
  summarise(count = n()) %>%
  mutate(
    dataset = "Arabidopsis",
    percent = count / sum(count) * 100
  )

IVT_motif_ranked_top25 <- at_pseU_IVT_2pad %>%
  group_by(motif) %>%
  summarise(count = n()) %>%
  mutate(
    dataset = "IVT",
    percent = count / sum(count) * 100
  )

### combine into 1 dataset

combined_motifs <- bind_rows(
  motif_ranked_top25 %>% select(motif, percent, dataset),
  IVT_motif_ranked_top25 %>% select(motif, percent, dataset)
)

### limits to top 25 motifs, can change n= to any number ###

top_motifs <- combined_motifs %>%
  group_by(motif) %>%
  summarise(avg_percent = mean(percent)) %>%
  arrange(desc(avg_percent)) %>%
  slice_head(n = 25) %>%
  pull(motif)


##### example graph to produce to show these motifs 

combined_motifs %>%
  filter(motif %in% top_motifs) %>%
  mutate(motif = factor(motif, levels = rev(top_motifs))) %>%
  ggplot(aes(x = percent, y = motif, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  labs(
    title = "Top 25 pseU 5-mer Motif Frequencies (Normalized)",
    x = "Percent of Total Sites",
    y = "Motif"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "top"
  )
