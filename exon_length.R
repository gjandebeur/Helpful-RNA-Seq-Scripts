### determining exon lengths in methylated vs unmethylated samples

#### dependencies ####

library(GenomicFeatures)
library(dplyr)
library(stringr)
library(ggplot2)

### load in txdb (change to your gtf for species studied)

txdb <- makeTxDbFromGFF("../data/at_ensembl_plants.gtf")

### extract exon length from txdb file

exons_df <- exons(txdb, columns = c("EXONID", "TXNAME", "GENEID")) %>% as.data.frame()
exons_df_clean <- exons_df %>%
  rename(
    chrom = seqnames,
    exon_start = start,
    exon_end = end
  ) %>%
  mutate(
    exon_length = exon_end - exon_start + 1,
    GENEID = str_remove(GENEID, "^gene:")
  )

### isolate m6A sites from native_filtered (produced in libraryprep_analyses.R)

m6A_sites <- native_filtered %>%
  filter(mod == "m6A") %>%
  mutate(GENEID = str_remove(GENEID, "^gene:")) %>%
  select(GENEID, m6a_pos = start)

#####################################
# determine which transcripts are methylated or unmethylated

methylated_exons <- exons_df_clean %>%
  inner_join(m6A_sites, by = "GENEID") %>%
  filter(m6a_pos >= exon_start & m6a_pos <= exon_end) %>%
  distinct(chrom, exon_start, exon_end, GENEID, .keep_all = TRUE)

# Identify unmethylated exons (exons without any m6A sites)
unmethylated_exons <- anti_join(
  exons_df_clean, methylated_exons,
  by = c("chrom", "exon_start", "exon_end", "GENEID")
)

###############################################
#produce one combined length file with both meth/unmeth lengths

length_df <- bind_rows(
  methylated_exons %>% select(exon_length) %>% mutate(type = "methylated"),
  unmethylated_exons %>% select(exon_length) %>% mutate(type = "unmethylated")
)

#### example plot to view the difference in length (violin & boxplot)

ggplot(length_df, aes(x = type, y = exon_length, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  labs(
    title = "Exon Lengths: m6A Methylated vs Unmethylated",
    x = NULL,
    y = "Exon Length (nt)"
  ) +
  scale_fill_manual(values = c(methylated = "darkgreen", unmethylated = "maroon")) +
  coord_cartesian(ylim = c(0, 2000)) +
  theme_minimal()
