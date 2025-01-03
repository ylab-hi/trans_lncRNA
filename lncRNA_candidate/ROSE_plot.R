library(tidyverse)
source('~/theme_manuscript.R')

# *****************************************************************************
# plot --------------------
# *****************************************************************************
library(ggrepel)


enhancers_with_genes <- read_tsv('rose/ROSE/rose_out_tx_exclude/enhancers_with_genes.tsv')

lncRNAs_to_label <- c('CASC19', 'NEAT1', 'MALAT1', 'PCAT1', 'PVT1', 'IGF1R-AS1')

enhancers_with_genes_plots <- enhancers_with_genes %>%
  mutate(label = ifelse(gene_name %in% lncRNAs_to_label, 
                        gene_name,
                        '')) %>%
  mutate(color = case_when(`h3k27ac_merged.bam` < 16317.36 ~ 'enhancer',
                           str_starts(CLOSEST_GENE, 'G.') ~ 'novel_lnc',
                           gene_name == 'IGF1R-AS1' ~ 'igf1r_as1',
                           TRUE ~ 'superenhancer'))

enhancers_with_genes_plots %>%
  ggplot(aes(x = (max(stitchedPeakRank) - stitchedPeakRank), 
             y = `h3k27ac_merged.bam`)) +
  geom_point(aes(color = color), size = 0.25,
             alpha = 0.35, shape = 15) +
  geom_hline(yintercept = 16317.36,
             linetype="dashed",
             linewidth=0.5,
             alpha=0.2) +
  geom_vline(xintercept = max(enhancers$stitchedPeakRank) - enhancers[enhancers$h3k27ac_merged.bam == 16317.36,]$stitchedPeakRank, 
             linetype="dashed",
             size=0.5,
             alpha=0.2)+
  geom_text_repel(aes(label = label),
                  nudge_x=-5000,
                  direction="y",force=2,max.iter=4000,
                  hjust=1.5,segment.size=0.4,
                  segment.curvature = -1e-20,
                  segment.ncp = 2,
                  segment.angle = 45,
                  segment.shape = 0,
                  segment.square = T,
                  size = 2,
                  color = ifelse(enhancers_with_genes_plots$label=="IGF1R-AS1","red","black"),
                  fontface=ifelse(enhancers_with_genes_plots$label=="IGF1R-AS1","bold","plain")) +
  labs(x = 'Enhancers ranked by signal (x10^4)',
       y = 'Signal at enhancers (x10)') +
  xlim(0, nrow(enhancers_with_genes)) +
  scale_color_manual(values = c('grey', 'red', 'red', '#9bf095')) +
  theme_manuscript_small() +
  theme(legend.position = 'none')

ggsave('plots/rose.pdf', device = 'pdf',
       height = 40, width = 50, units = 'mm')



#### label genes ----------------------------------------

genes_to_label <- c('AR', 'TMPRSS2', 'KLK2', 'FOXA1', 'ACP2', 'NDRG1','MYC', 'ERG', 'CCND1')

lncRNAs_to_label <- c('CASC19', 'NEAT1', 'MALAT1', 'PCAT1', 'PVT1', 'IGF1R-AS1')

enhancers_with_genes_plots <- enhancers_with_genes %>%
  mutate(label = ifelse(gene_name %in% c(lncRNAs_to_label, genes_to_label), 
                        gene_name,
                        '')) %>%
  mutate(color = case_when(`h3k27ac_merged.bam` < 16317.36 ~ 'enhancer',
                           str_starts(CLOSEST_GENE, 'G.') ~ 'novel_lnc',
                           gene_name == 'IGF1R-AS1' ~ 'igf1r_as1',
                           TRUE ~ 'superenhancer'))

enhancers_with_genes_plots %>%
  # only label first occurence
  mutate(duplicated = duplicated(gene_name)) %>%
  mutate(label = ifelse(duplicated, '', label)) %>%
  ggplot(aes(x = (max(stitchedPeakRank) - stitchedPeakRank), 
             y = `h3k27ac_merged.bam`)) +
  geom_point(aes(color = color), size = 0.25,
             alpha = 0.35, shape = 15) +
  geom_hline(yintercept = 16317.36,
             linetype="dashed",
             linewidth=0.5,
             alpha=0.2) +
  geom_vline(xintercept = max(enhancers$stitchedPeakRank) - enhancers[enhancers$h3k27ac_merged.bam == 16317.36,]$stitchedPeakRank, 
             linetype="dashed",
             size=0.5,
             alpha=0.2)+
  geom_text_repel(aes(label = label),
                  nudge_x=-5000,
                  direction="y",force=2,max.iter=4000,
                  hjust=1.5,segment.size=0.2,
                  segment.curvature = -1e-20,
                  segment.ncp = 2,
                  segment.angle = 45,
                  segment.shape = 0,
                  segment.square = T,
                  size = 1.5,
                  color = ifelse(enhancers_with_genes_plots$label=="IGF1R-AS1","red","black"),
                  fontface=ifelse(enhancers_with_genes_plots$label=="IGF1R-AS1","bold","plain")) +
  labs(x = 'Enhancers ranked by signal (x10^4)',
       y = 'Signal at enhancers (x10)') +
  xlim(0, nrow(enhancers_with_genes)) +
  scale_color_manual(values = c('grey', 'red', '#9bf095')) +
  theme_manuscript_small() +
  theme(legend.position = 'none')


  ggsave('plots/rose_lnc_genes.pdf', device = 'pdf',
       height = 40, width = 50, units = 'mm')
