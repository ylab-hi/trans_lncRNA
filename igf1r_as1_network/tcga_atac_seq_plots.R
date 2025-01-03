#### This script uses TCGA ATAC data that is too large to host on github.

library(tidyverse)
library(GenomicRanges)

theme_manuscript_small <- theme_manuscript_small() %+replace% 
  theme(panel.border = element_rect(linewidth = 0.5, fill = NA))

# *****************************************************************************
# load data --------------------
# *****************************************************************************


gtf <- rtracklayer::import('data/annotations/gencode.v37.igf1r_as1.gtf')

prad_fpkm <- read_tsv('/projects/b1171/twp7981/backup/hpc/projects/lncRNAs/lncRNAs_data/SU2C_PolyA/IGF1R-AS1_count_PacBio_forms_PRAD_N_T_SU2C/PRAD_T/PRAD.tumors.counts.FPKM.txt')

tcga_atac <- read_tsv('/home/jfp0082/ylab/jfp0082/projects/igf1r_atac/data/tcga_atac/TCGA_ATAC_peak_Log2Counts_dedup_sample.gz')
tcga_atac_regions <- read_tsv('/home/jfp0082/ylab/jfp0082/projects/igf1r_atac/data/tcga_atac/peak_regions_to_id.tsv') %>%
  GRanges()

tcga_fpkm <- data.table::fread('data/tcga_atac/GDC-PANCAN.htseq_fpkm-uq.tsv.gz')

cols <- c('xena_sample', tcga_atac %>% select(!1) %>% colnames)
cols <- cols[cols %in% colnames(tcga_fpkm)]
tcga_fpkm[, ..cols]

tcga_fpkm_atac_samples <- tcga_fpkm[, ..cols] %>% as_tibble()
tcga_fpkm_atac_samples <- tcga_fpkm_atac_samples %>%
  dplyr::rename(gene_id = xena_sample) %>%
  mutate(gene_id = str_split_i(gene_id, '\\.', 1)) %>%
  left_join(gtf %>% 
              as_tibble() %>%
              select(gene_id, gene_name) %>%
              unique())
  

tcga_canonical_splice <- read_tsv('data/tcga_atac/tcga_prad_samples_canonical_splice.tsv')

prad_fpkm %>% select(!1) %>%
  colnames() %>%
  str_sub(1, 16) -> prad_samples

prad_atac <- tcga_atac %>%
  select(sample, any_of(prad_samples))


prad_igf1r_as1 <- prad_fpkm %>%
  filter(Geneid == 'PB.69') %>%
  pivot_longer(!1, names_to = 'tcga_long_id', values_to = 'lncRNA_exp')

get_prad_exp <- function (gene_id) {
  prad_fpkm %>%
    filter(Geneid == gene_id) %>%
    pivot_longer(!1, names_to = 'tcga_long_id', values_to = 'exp')
}

# *****************************************************************************
# CCAT1 --------------------
# *****************************************************************************

CCAT1_region <- tcga_atac_regions %>%
  subsetByOverlaps(GRanges('chr8:127215335-127215578')) %>%
  .$name

prad_atac %>%
  filter(sample == CCAT1_region) %>%
  pivot_longer(!1, names_to = 'tcga_id', values_to = 'CCAT1_enh_exp') %>%
  left_join(get_prad_exp('PCAT1') %>%
              mutate(tcga_id = str_sub(tcga_long_id, 1, 16))) %>%
  ggplot(aes(x = exp, y = CCAT1_enh_exp)) +
  geom_point() +
  geom_smooth(method = 'lm')

# CCAT1 enhancer region vs IGF1R-AS1 exp
prad_atac %>%
  filter(sample == CCAT1_region) %>%
  pivot_longer(!1, names_to = 'tcga_id', values_to = 'CCAT1_enh_exp') %>%
  left_join(get_prad_exp('PB.69') %>%
              mutate(tcga_id = str_sub(tcga_long_id, 1, 16))) %>%
  left_join(tcga_canonical_splice, by = c('tcga_id' = 'sample')) %>%
  mutate(exp = ifelse(lnc_splces > 0, exp, 0)) %>%
  mutate(status = ifelse(lnc_splces > 0, '+', '-')) %>%
  ggplot(aes(x = status, y = CCAT1_enh_exp, color = status)) +
  geom_boxplot(aes(fill = status), alpha = 0.2, outlier.shape = NA, 
               width = 0.2, linewidth = 0.5) +
  ggbeeswarm::geom_quasirandom(method = 'quasirandom', size = 0.25) + 
  ggpubr::stat_compare_means(method = 't.test', comparisons = list(c('+',
                                                                     '-')),
                             size = 2.5,
                             label = 'p.signif') +
  scale_color_manual(values = c('#ff5500', '#4b96fa')) +
  scale_fill_manual(values = c('#ff5500', '#4b96fa')) +
  theme_manuscript_small() +
  ylim(c(-.1, 4.3)) +
  theme(panel.border = element_rect(linewidth = 0.5, fill = NA),
        legend.position = 'none') +
  labs(x = 'IGF1R-AS1 status', y = 'PRAD CCAT1 Enhancer Signal')

ggsave('plots/supplemental/prad_ccat1_igf1r_as1.pdf', device = 'pdf',
       width = 26, height = 34, units = 'mm')

# CCAT1 enhancer region is correlated with MYC
data <- tcga_atac %>%
  filter(sample == CCAT1_region) %>%
  pivot_longer(!1, names_to = 'tcga_id', values_to = 'CCAT1_enh_exp') %>%
  left_join(tcga_fpkm_atac_samples %>%
              filter(gene_id == 'ENSG00000136997') %>%
              pivot_longer(!(gene_id | gene_name), names_to = 'tcga_id', values_to = 'MYC_exp')
  ) %>%
  mutate(PRAD = tcga_id %in% prad_samples)

data %>%
  ggplot(aes(x = CCAT1_enh_exp, y = MYC_exp, color = PRAD)) +
  geom_point(size = 0.25) + 
  geom_smooth(data = data %>% filter(PRAD), method = MASS::rlm) +
  ggpubr::stat_cor( 
                   label.x.npc = 'left', 
                   label.y.npc = 'bottom',
                   size = 2.5) +
  scale_color_manual(values = c('grey', '#ff5500')) +
  theme_manuscript_small() +
  theme(legend.position = 'none') +
  labs(x = 'CCAT1 Enhancer Signal', y = 'MYC Expression (FPKM-UQ > 0)')

ggsave('plots/supplemental/prad_ccat1_enhancer_myc.pdf', device = 'pdf',
       width = 45, height = 34, units = 'mm')
  

# CCAT1 enhancer region is correlated with CCAT1
data <- tcga_atac %>%
  filter(sample == CCAT1_region) %>%
  pivot_longer(!1, names_to = 'tcga_id', values_to = 'CCAT1_enh_exp') %>%
  left_join(tcga_fpkm_atac_samples %>%
              filter(gene_id == 'ENSG00000247844') %>%
              pivot_longer(!(gene_id | gene_name), 
                           names_to = 'tcga_id', values_to = 'CCAT1_exp')
  ) %>%
  mutate(PRAD = tcga_id %in% prad_samples)

data %>%
  filter(CCAT1_exp > 0) %>%
  ggplot(aes(x = CCAT1_enh_exp, y = CCAT1_exp, color = PRAD)) +
  geom_point(size = 0.5) + 
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(label.x.npc = 'left', 
                   label.y.npc = 'top',
                   size = 2.5) +
  scale_color_manual(values = c('grey', '#ff5500')) +
  theme_manuscript_small +
  theme(legend.position = 'none') +
  labs(x = 'CCAT1 Enhancer Signal', y = 'CCAT1 Expression (FPKM-UQ > 0)')

ggsave('plots/supplemental/prad_ccat1_enhancer_ccat1.pdf', device = 'pdf',
       width = 50, height = 40, units = 'mm')


# CCAT1 enhancer region is correlated with CCAT1
data <- tcga_fpkm_atac_samples %>%
  filter(gene_id == 'ENSG00000247844' | gene_name == 'MYC') %>%
  pivot_longer(!(gene_id | gene_name), 
               names_to = 'tcga_id', values_to = 'exp') %>%
  pivot_wider(names_from = 'gene_name', values_from = 'exp',
              id_cols = tcga_id) %>%
  mutate(PRAD = tcga_id %in% prad_samples)


data %>%
  filter(`NA` > 0 & MYC > 0) %>%
  ggplot(aes(x = MYC , y = `NA`, color = PRAD)) +
  geom_point(size = 0.5) + 
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(label.x.npc = 'left', 
                   label.y.npc = 'top',
                   size = 2.5) +
  scale_color_manual(values = c('grey', '#ff5500')) +
  theme_manuscript_small +
  theme(legend.position = 'none') +
  labs(x = 'MYC Expression (FPKM-UQ > 0)', y = 'CCAT1 Expression (FPKM-UQ > 0)')

ggsave('plots/supplemental/prad_myc_ccat1.pdf', device = 'pdf',
       width = 50, height = 40, units = 'mm')
