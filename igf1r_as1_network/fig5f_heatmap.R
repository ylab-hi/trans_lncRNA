library(tidyverse)

# data was generated using standard code as described in methods section
# it's a lot of code and too distracting to upload here. 
data <- readRDS('data/heatmap_data.rds')

# *****************************************************************************
# build heatmap --------------------
# *****************************************************************************
library(tidygraph)
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
  filter(!str_starts(name, 'chr11')) %>%
  head(n = 10) %>%
  pull(name) %>%
  GRanges() -> top_100_peaks

overlaps <- findOverlaps(peak_ranges_hg38, top_100_peaks)
overlaps <- tibble(their_peaks = peak_ranges_hg38[overlaps@from] %>% as.character(),
                   my_peaks = top_100_peaks[overlaps@to] %>% as.character()) 
tmp <- peak_ranges_hg38[findOverlaps(peak_ranges_hg38, top_100_peaks) %>% .@from] 
overlaps <- overlaps %>%
  mutate(peak_id = tmp@elementMetadata$peak_id)
overlaps <- overlaps %>%
  bind_cols(peak_counts[tmp@elementMetadata$peak_id,] %>%
              select(!(1 | pos_id)))

overlaps <- overlaps %>%
  pivot_longer(4:43, names_to = 'sample', values_to = 'counts') %>%
  group_by(my_peaks, their_peaks) %>%
  mutate(avg_counts = mean(counts)) %>%
  group_by(my_peaks, sample) %>%
  slice_max(order_by = avg_counts, n = 1, with_ties = F) %>%
  select(my_peaks, sample, counts) %>%
  pivot_wider(names_from = 'sample', values_from = 'counts') %>%
  ungroup

mat <- overlaps %>%
  select(!1) %>%
  as.matrix
rownames(mat) <- overlaps$my_peaks


#### set annotations ----------------------------------------
myc_counts <- log_norm %>% filter(gene_name == 'MYC') %>%
  select(!(1:2)) %>%
  pivot_longer(everything(), names_to = 'sample', values_to = 'myc_exp') %>%
  with(setNames(log10(myc_exp), sample))
smarca4_counts <- log_norm %>% filter(gene_name == 'SMARCA4') %>%
  select(!(1:2)) %>%
  pivot_longer(everything(), names_to = 'sample', values_to = 'smarca4_exp') %>%
  with(setNames(log10(smarca4_exp), sample))
lnc_spm <- spm %>%
  with(setNames(spm, sample))

#### creat heatmap ----------------------------------------
colors <- list('Group' = setNames(viridis::viridis(4), 
                                  c('AR', 'NEPC', 'SCL', 'WNT')),
               'MYC' = circlize::colorRamp2(c(-2, max(myc_counts)), c("white", "indianred")),
               'IGF1R-AS1' = circlize::colorRamp2(c(0, min(lnc_spm[lnc_spm[colnames(mat)] > 0]), 4), c("white", "#f3f0f7", "#8968cd")),
               'SMARCA4' = circlize::colorRamp2(c(-1, max(smarca4_counts)), c("white", "lightsalmon2")),
               'AR activity' = circlize::colorRamp2(c(0, 1), c("white", "seagreen3"))
               
)

col_ann <- HeatmapAnnotation(
  Group = group_assignments_vector[colnames(mat)],
  `AR activity` = scales::rescale(ar_activity[colnames(mat)]),
  SMARCA4 = smarca4_counts[colnames(mat)],
  MYC = myc_counts[colnames(mat)],
  `IGF1R-AS1` = lnc_spm[colnames(mat)],
  col = colors,
  border = F,
  simple_anno_size = unit(0.35, 'cm'),
  gp = gpar(fontsize = 9),
  annotation_legend_param = list(nrow = 2)
)


set.seed(222)
hm <- ComplexHeatmap::Heatmap(matrix = mat,
                              column_km = 2,
                              top_annotation = col_ann,
                              column_names_gp = gpar(fontsize = 9),
                              row_names_gp = gpar(fontsize = 9),
                              border = T)
pdf('plots/main/science_cohort_heatmap.pdf', height = 3.5, width = 8.25)
ComplexHeatmap::draw(hm, annotation_legend_side = 'left')
dev.off()
