
### script to visualize poly-A tail lengths 

## dependencies ## 
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(tidytable)
  library(ggplot2)
  library(stringr)
  library(dplyr)
})


### load in pA tail lengths for each sample 

polyA_data <- readr::read_tsv("data/polyA_tail_lengths.tsv") %>%
  rename(
    tail_length = <tail_length_column>,  # <-- Replace with actual column name
    gene_annotation = <annotation_column> # <-- if gene info is embedded here
  ) %>%
  mutate(
    gene_id = str_extract(gene_annotation, 'gene_id "gene:[^"]+') %>% 
              str_replace('gene_id "gene:', "")
  ) %>%
  filter(!is.na(tail_length), !is.na(gene_id))

#### load the final filtered data from the library preprocessing script. ####
final_filtered_data <- readr::read_csv("data/final_filtered_data.csv") 

## seperating pA length into 3 categories, seperated by <Q1 , Q1-Q3, and >Q3 as "short", "medium", and "long"

gene_level_tails <- polyA_data %>%
  group_by(gene_id) %>%
  summarise(
    mean_tail_length = mean(tail_length, na.rm = TRUE),
    median_tail_length = median(tail_length, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    tail_length_group = case_when(
      mean_tail_length <= quantile(mean_tail_length, 0.25, na.rm = TRUE) ~ "short",
      mean_tail_length >= quantile(mean_tail_length, 0.75, na.rm = TRUE) ~ "long",
      TRUE ~ "medium"
    )
  )

## combine the two files for graphing and join by gene_id (important)

mod_with_tail <- filtered_mods %>%
  left_join(gene_level_tails %>% select(gene_id, tail_length_group), by = "gene_id") %>%
  filter(!is.na(tail_length_group))

### basic distribution of pA lengths 

ggplot(gene_level_tails, aes(x = mean_tail_length)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Mean PolyA Tail Lengths (Per Gene)",
       x = "Mean PolyA Tail Length (nt)",
       y = "Gene Count")


### and graph distributions of % mod to pA length, can change y axis to other columns as well

ggplot(mod_with_tail, aes(x = tail_length_group, y = avg_percent_modified, fill = tail_length_group)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", color = "gray20", alpha = 0.5) +
  facet_wrap(~ mod, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    title = "% Modified by PolyA Tail Length Group",
    x = "PolyA Tail Length Group",
    y = "Average % Modified"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")
