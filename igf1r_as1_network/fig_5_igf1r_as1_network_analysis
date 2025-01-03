library(tidyverse)
source('~/theme_manuscript.R')
# *****************************************************************************
# load in data --------------------
# *****************************************************************************


gene_peak_links <- read_tsv('data/analysis_data/gene_peak_links_smarca4_fdr_05_l2f_025.tsv')
peak_peak_links <- read_tsv('data/analysis_data/peak_peak_links_fdr_05_l2f_025.tsv')
diff_exp_df <- read_tsv('/projects/b1171/twp7981/backup/hpc/projects/lncRNAs/IGF1R-AS1_DEGs/meta_analysis/si025.meta.txt')
gene_peak_links <- gene_peak_links %>%
  mutate(significant = rna_log2fc < -0.25 & 
           rna_padj < 0.05)

# *****************************************************************************
#  Enhancer-gene Network --------------------
# *****************************************************************************

# *****************************************************************************
# build graph --------------------
# *****************************************************************************
library(tidygraph)
gpl_graph <- gene_peak_links %>%
  filter(significant) %>%
  # filter(!is.na(gene_name)) %>%
  mutate(log2fc_strength = peak_log2fc * rna_log2fc) %>%
  mutate(
    log2fc_strength = rank(log2fc_strength),
    cpm_strength = avg_replicate_cpm) %>%
  mutate(cpm_strength = rank(cpm_strength)) %>%
  mutate(rank_score = scales::rescale(rank(log2fc_strength*cpm_strength))) %>%
  mutate(connection_type = ifelse(gene_peak_link == 'promoter', 'promoter', 'loop')) %>%
  select(pos_id, gene_name, log2fc_strength, cpm_strength, connection_type, rank_score)  %>%
  unique %>%
  rename(to = pos_id,
         from = gene_name) %>%
  # some links are both 'promoter' and 'loop' 
  # (because loop is contained within promoter regions)
  # in this case, just change to promoter
  group_by(to, from) %>%
  mutate(dual = n() > 1) %>%
  mutate(connection_type = ifelse(dual, 'promoter', connection_type)) %>%
  select(!dual) %>%
  unique

gpl_graph_plus_background <- bind_rows(
  gpl_graph,
  tibble(to = 'IGF1R-AS1',
         from = c(gpl_graph$to, gpl_graph$from) %>% unique,
         log2fc_strength = 1e-2,
         cpm_strength = 1e-2,
         rank_score = 1e-2,
         connection_type = 'background')
)


pp_graph <-
  peak_peak_links %>%
  left_join(gene_peak_links %>%
              select(pos_id, peak_log2fc, avg_replicate_cpm) %>%
              unique %>%
              rename(pos1_log2fc = peak_log2fc,
                     pos1_cpm = avg_replicate_cpm), by = 'pos_id') %>%
  left_join(gene_peak_links %>%
              select(pos_id, peak_log2fc, avg_replicate_cpm) %>%
              unique %>%
              rename(pos2_log2fc = peak_log2fc,
                     pos2_cpm = avg_replicate_cpm), by = c('pos_id_2'= 'pos_id')) %>%
  mutate(log2fc_strength = pos1_log2fc * pos2_log2fc,
         log2fc_strength = rank(log2fc_strength)
  ) %>%
  mutate(cpm_strength = map2_dbl(pos1_cpm, pos2_cpm, ~ mean(c(.x, .y))))  %>%
  mutate(cpm_strength = rank(cpm_strength)) %>%
  mutate(rank_score = scales::rescale(rank(log2fc_strength*cpm_strength))) %>%
  select(pos_id, pos_id_2, log2fc_strength, cpm_strength, rank_score) %>%
  unique %>% 
  rename(to = pos_id,
         from = pos_id_2) %>%
  mutate(connection_type = 'peak_peak')


# combine ------
graph <- bind_rows(gpl_graph,
                   pp_graph) %>%
  ungroup %>%
  as_tbl_graph(directed = F) %E>%
  distinct

# add some metadata
graph <- graph %>%
  activate(nodes) %>%
  mutate(node_type = case_when(
    str_detect(name, 'chr') ~ 'PEAK',
    TRUE ~ 'GENE'
  ))

# combine with background links ------
graph_plus_background <- bind_rows(gpl_graph_plus_background,
                                   pp_graph) %>%
  ungroup %>%
  as_tbl_graph(directed = F)

# add some metadata
graph_plus_background <- graph_plus_background %>%
  activate(nodes) %>%
  mutate(node_type = case_when(
    str_detect(name, 'chr') ~ 'PEAK',
    TRUE ~ 'GENE'
  ))

# filter out peak nodes that are not connected to genes
graph_plus_background <- graph_plus_background %>%
  activate(nodes) %>%
  mutate(group = group_components()) %>%
  # group 1 is the subgraph with all genes-peaks
  filter(group == 1)

graph_plus_background %>%
  activate(edges) %>%
  distinct() %>%
  activate(nodes) %>%
  mutate(hub_score = centrality_eigen(
    weights = rank_score,
    scale = F
  )) %>%
  filter(node_type == 'GENE') %>%
  as_tibble() -> c_scores


graph_plus_background %>%
  activate(edges) %>%
  distinct() %>%
  activate(nodes) %>%
  mutate(hub_score = centrality_eigen(
    weights = rank_score,
    scale = F
  )) %>%
  filter(node_type == 'GENE') %>%
  as_tibble() -> c_scores

#### plot top 20 ----------------------------------------
c_scores %>%
  arrange(-hub_score) %>%
  filter(name != 'IGF1R-AS1') %>%
  head(n = 20) %>%
  ggplot(aes(x = reorder(name, -hub_score), y = log10(hub_score * 10/min(hub_score)))) + 
  geom_col(fill = '#ffb875', width = 0.75) +
  labs(x = NULL, y = 'Normalized\nEigencentrality') + 
  theme_manuscript_small() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.border = element_blank()) 

ggsave('plots/top_20_centrality_scores_gene.pdf', device = 'pdf',
       height = 25, width = 75, units = 'mm')  

graph_plus_background %>%
  activate(edges) %>%
  distinct() %>%
  activate(nodes) %>%
  mutate(hub_score = centrality_eigen(
    weights = rank_score,
    scale = F
  )) %>%
  filter(node_type == 'PEAK') %>%
  as_tibble() %>%
  arrange(-hub_score) %>%
  # filter(name != 'IGF1R-AS1') %>%
  head(n = 20) %>%
  ggplot(aes(x = reorder(name, -hub_score), y = log10(hub_score * 10/min(hub_score)))) + 
  geom_col(fill = '#ffb875', width = 0.75) +
  labs(x = NULL, y = '-log10(Centrality Score)') + 
  theme_manuscript_small() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        panel.border = element_blank()) 

ggsave('plots/top_20_centrality_scores_peak.pdf', device = 'pdf',
       height = 45, width = 75, units = 'mm')  


# *****************************************************************************
# identify subgraphs --------------------
# *****************************************************************************

graph %>%
  activate(nodes) %>%  
  mutate(subgraph_group = group_components()) %>%
  as_tibble() -> subgraph_groups

graph %>%
  activate(nodes) %>%  
  mutate(group = group_components()) %>%
  group_by(group) %>%
  mutate(n = n()) %>%
  ungroup %>%
  left_join(c_scores) %>%
  mutate(c_score_rank = rank(-hub_score)) %>%
  group_by(group) %>%
  mutate(sig_group = any(c_score_rank <= 20)) %>%
  ungroup %>%
  filter(sig_group) %>%
  as_tibble() %>% pull(name) -> nodes_to_show


# *****************************************************************************
# plot for figure 5c --------------------
# final plot is cleaned up for publication in illustrator
# *****************************************************************************


graph %>%
  activate(nodes) %>%
  mutate(group = group_components()) %>%
  mutate(gene = !str_detect(name, 'chr')) %>%
  mutate(label = ifelse(gene, name, '')) %>%
  # left_join(subgraph_groups) %>%
  # group 1 is a large subgraph with many low scoring nodes
  # it will be almost impossible to fit in main figure
  # thus, we'll only include it in supplemental graph
  filter(group != 1) %>%
  left_join(c_scores %>% select(!group)) %>%
  mutate(c_score_rank = rank(-hub_score)) %>%
  filter(name != 'IGF1R-AS1') %>%
  mutate(size = ifelse(c_score_rank <= 20, 21 - c_score_rank, -1)) %>%
  mutate(label = ifelse(size > 0, label, '')) %>%
  group_by(group) %>%
  mutate(sig_group = any(c_score_rank <= 20)) %>%
  ungroup %>%
  filter(sig_group) %>%
  mutate(color = ifelse(gene,
                        ifelse(name == 'MYC',
                               hub_score/5000,
                               hub_score),
                        NA)) %>%
  ggraph(layout = 'fr', ) +
  geom_edge_fan(aes(width = rank_score,),
                alpha = 0.8, color = "grey55", 
                show.legend = T,
                n = ) +
  geom_node_point(aes(color = color, size = size), alpha = 1, show.legend = T) +
  geom_node_text(aes(label = label),  colour = 'black', size=1.75,
                 show.legend = FALSE, repel = F) +
  theme_graph(background = "white") +
  scale_edge_width(range = c(0.1, 0.3)) +
  scale_size(range = c(0.2, 4.5)) +
  scale_color_gradient(low = 'grey85',high = '#ff5500',
                       na.value = '#6ece58')


ggsave('plots/igf1r_as1_network_v2.pdf', device = 'pdf',
       height = 100, width = 100, units = 'mm')
