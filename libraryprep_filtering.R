# R script for RNA modification analyses (motif, enrichment, )

####################################
#load dependencies
####################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(AnnotationDbi)
library(GenomicFeatures)
library(ggplot2)
library(ggsci)
library(cowplot)
library(scales)
library(svglite)
library(readr)
library(tximport)
library(edgeR)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(dbplyr)

#####################################
#set colnames and read files

names <- c("chrom", "start", "end", "code", "score", "strand", "start_1", "end_1", "color", 
           "n_valid_cov", "percent_modified", "n_mod_reads", "n_canonical_reads", "n_other_mod_reads", 
           "n_deletion_readds", "n_fail_reads", "n_diff_reads", "n_nocall_reads", 
           "tx_chrom", "tx_start", "tx_end", "tx", "tx_score", "tx_strand", "tx_overlap")

native1 <- fread("../data/native1.bed", header = FALSE) %>% setNames(names)
native2 <- fread("../data/native2.bed", header = FALSE) %>% setNames(names)
IVT1 <- fread("../data/IVT1.bed", header = FALSE) %>% setNames(names)
IVT2 <- fread("../data/IVT2.bed", header = FALSE) %>% setNames(names)

###############################
#annotate for reproducibility. changing code to mod, loading in txdb 
###############################
txdb <- makeTxDbFromGFF("../data/annotation.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene_at <- AnnotationDbi::select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

process_data <- function(df, label) {
  df %>% 
    inner_join(tx2gene_at, by = c("tx" = "TXNAME")) %>%
    unite(position, chrom:end, remove = FALSE) %>%
    separate(code, c("code", "motif", "offset"), sep = ",") %>%
    filter((code == "a" & motif == "A") | (code == "m" & motif == "C") | 
           (code == "17596" & motif == "A") | (code == "69426" & motif == "A") | 
           (code == "19228" & motif == "C") | (code == "19229" & motif == "G") | 
           (code == "17802" & motif == "T") | (code == "19227" & motif == "T")) %>%
    mutate(experiment = label,
           mod = case_when(
             code == "a" ~ "m6A",
             code == "m" ~ "m5C",
             code == "17802" ~ "pseU",
             code == "17596" ~ "inosine",
             code == "69426" ~ "2'O-Me A",
             code == "19228" ~ "2'O-Me C",
             code == "19229" ~ "2'O-Me G",
             code == "19227" ~ "2'O-Me U"
           ))
}

native1_filtered <- process_data(native1, "native1")
native2_filtered <- process_data(native2, "native2")
IVT1_filtered <- process_data(IVT1, "IVT1")
IVT2_filtered <- process_data(IVT2, "IVT2")

#############################
#combine library prep for downstream analyses
############################

native_positions <- bind_rows(native1_filtered, native2_filtered) %>%
  group_by(position) %>%
  summarise(n_native = n_distinct(experiment), .groups = "drop") %>%
  filter(n_native == 2) %>% # keep only in both
  distinct(position)

# Identify all IVT positions
ivt_positions <- bind_rows(IVT1_filtered, IVT2_filtered) %>%
  distinct(position, .keep_all = TRUE)

# Find IVT positions found in both replicates
ivt_positions_both <- bind_rows(IVT1_filtered, IVT2_filtered) %>%
  group_by(position) %>%
  summarise(n_ivt = n_distinct(experiment), .groups = "drop") %>%
  filter(n_ivt == 2)

ivt_2filtered <- bind_rows(IVT1_filtered, IVT2_filtered) %>%
  inner_join(ivt_positions_both, by = "position")

# Filter for native-only reproducible sites
native_filtered <- bind_rows(native1_filtered, native2_filtered) %>%
  inner_join(native_positions, by = "position") %>%
  anti_join(ivt_positions, by = "position")
#######################################

# Save combined filtered data
write.csv(native_filtered, "../data/native_filtered_clean.csv", row.names = FALSE)

#downstream analyses to investigate further in the repository.
# - expression_stoich_polyA_analysis.R
# - motif_sequence_processing.R
# - enrichment_GO_analysis.R
