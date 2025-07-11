### Script to view stoichiometry (useful for downstream comparisons)

### dependencies ### 

library(tidyverse)
library(ggpubr)  
library(viridis) 

### load in the final filtered data file from preprocessing script.
### need the percent modified or n_mod and n_cov to determine stoich

filtered_mod_file <- "path/to/filtered_mod_data.csv"

## fixing for reproducibility 
filtered_mods <- readr::read_csv(filtered_mod_file) %>%
  mutate(gene_id = str_remove(gene_id, "^gene:"))

### percent modified should be present but if not here is calc for site-level stoich 

native_stoich <- native_filtered %>%
  select(position, n_valid_cov, n_mod_reads, GENEID, mod, experiment) %>%
  mutate(
    stoich = n_mod_reads / n_valid_cov,
    GENEID = str_remove(GENEID, "^gene:")
  )

######################
## basic graph 
## histogram of site level stoich, combines all samples 

ggplot(native_stoich, aes(x = stoich)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  labs(
    title = "Stoichiometry Distribution",
    x = "Stoichiometry",
    y = "Count"
  ) +
  theme_minimal()

