######### load dependencies

######
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(dplyr)
#######################

# Load in past library prep file
native_filtered <- fread("../data/native_filtered_clean.csv", row.names = F)

####
mod_gene_data <- native_filtered %>%
  mutate(
    gene_id = str_remove(GENEID, "^gene:"),
    mod = str_trim(mod)
  ) %>%
  select(gene_id, mod, experiment) %>%
  distinct() %>%
  as.data.frame()

mod_gene_list <- mod_gene_data %>%
  group_by(mod) %>%
  summarise(genes = list(unique(gene_id))) %>%
  as.data.frame()

mod_gene_list_named <- setNames(mod_gene_list$genes, mod_gene_list$mod)
############################################

#begin enrichment code

ego <- compareCluster(
  geneCluster = mod_gene_list_named,
  fun = "enricher",
  universe = background$gene,
  pvalueCutoff = 0.05,
  TERM2GENE = goterms
)

hyper_ego <- mutate(ego@compareClusterResult,
                    FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

#############################################
#some potential graphs to help visualize
#bubble plot

library(ggplot2)
library(cowplot)

ggplot(ego, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "forestgreen", high = "darkviolet", name = "Adj p-value") +
  theme_cowplot() +
  labs(
    title = "Ratio of GO Enrichment by Modification",
    x = "Modification",
    y = "GO Term",
    size = "Gene Count"
  )
