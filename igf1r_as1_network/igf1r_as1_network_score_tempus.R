# RAW TEMPUS DATA NOT AVAILABLE TO PUBLIC

library(tidyverse)


rna_wide <- read_tsv('local_data/g_rna_gene_expression_log2_corrected_wide.tsv')

gene_id_to_name <- read_tsv('local_data/gene_id_to_name.tsv')

# convert

rna <- rna_wide %>%
  pivot_longer(!(1:2), names_to = 'ensembl_gene', values_to = 'tpm') %>%
  left_join(gene_id_to_name)


# cscores
c_scores <- read_tsv('local_data/centrality_scores.tsv')

os_data <- read_tsv('local_data/clinical/os_times.tsv')

# *****************************************************************************
# calculate igf1r_as1 act scores --------------------
# *****************************************************************************

# get genes with the top 100 centrality scores
top_c_genes <- c_scores %>%
  filter(node_type == 'GENE') %>%
  #slice_max(hub_score, n = 200) %>%
  pull(name)


igf1r_as1_zscores <- rna %>%
  # only include genes with top centrality scores
  filter(human_gene_symbol %in% top_c_genes) %>%
  # pivot so that we have a sample_name column and expression column
  # pivot_longer(!(1), names_to = 'sample_name', values_to = 'exp') %>%
  # calculate z-score for each gene
  group_by(human_gene_symbol) %>%
  mutate(z = scale(tpm)[,1]) %>%
  # sum z-scores across samples
  group_by(patient_id) %>%
  summarize(igf1r_as1_activity = sum(z)) %>%
  # finally, take z-score of the sum of z-scores
  mutate(igf1r_as1_activity = scale(igf1r_as1_activity)[,1])


# *****************************************************************************
# section --------------------
# *****************************************************************************
rna_wide <- read_tsv('local_data/g_rna_gene_expression_log2_corrected_wide.tsv')
gene_id_to_name <- read_tsv('local_data/gene_id_to_name.tsv')

rna <- rna_wide %>%
  pivot_longer(!(1:2), names_to = 'ensembl_gene', values_to = 'tpm') %>%
  left_join(gene_id_to_name)


rna <- rna %>% select(patient_id, tpm, ensembl_gene) %>%
  pivot_wider(names_from = patient_id, values_from = tpm)

rna_mat <- rna %>% column_to_rownames('ensembl_gene')

ar_myc_net <- msigdbr::msigdbr(category = 'H') %>%
  filter(gs_name %in% c('HALLMARK_ANDROGEN_RESPONSE', 'HALLMARK_MYC_TARGETS_V1')) %>%
  select(gs_name, ensembl_gene) %>%
  mutate(mor = 1) %>%
  rename(source = gs_name, target = ensembl_gene) %>%
  distinct()

viper_ar_myc <- decoupleR::run_wmean(rna_mat, ar_myc_net)

igf1r_myc_ar <- viper_ar_myc %>%
  filter(statistic == 'norm_wmean') %>%
  as_tibble() %>%
  pivot_wider(names_from = source, values_from = score, id_cols = condition) %>%
  rename(patient_id = condition,
         AR_response = HALLMARK_ANDROGEN_RESPONSE,
         MYC_targets = HALLMARK_MYC_TARGETS_V1) %>%
  #left_join(gene_id_to_name) %>%
  left_join(igf1r_as1_zscores)


# *****************************************************************************
# get clinical plots, igf1r_as1 --------------------
# *****************************************************************************
library(survival)
library(survminer)

igf1r_myc_ar_clin <- igf1r_myc_ar %>%
  inner_join(os_data)

# pspline analysis shows that large values of igf1r_as1_activity
#   are linked to decreased survival, justifying optimal cutpoint visualization
coxph(Surv(times, status) ~ pspline(igf1r_as1_activity, df=4),
      data = igf1r_myc_ar_clin) %>% summary

surv_cut <- surv_cutpoint(data = igf1r_myc_ar_clin,
                          time = 'times',
                          event = 'status',
                          variables = c( 'igf1r_as1_activity'))

surv_cat <- surv_categorize(surv_cut)

fit <- survfit(Surv(times/ 30.44, status) ~ igf1r_as1_activity,
               data = surv_cat)

ggsurvplot(fit, data = surv_cat, pval = T,
           risk.table = T)

plot <-
  ggsurvplot(fit, data = surv_cat,
             pval =T,
             pval.coord = c(0.05, 0.1),
             risk.table = T,
             palette = c('red3', 'blue3'),
             tables.theme = theme_cleantable(),
             censor.size = 4.25,
             # lwd = 1.5,
             legend.title="",
             legend.text = element_text(size = 2, color = "black", face = "bold"),
             xlab = 'Months',
             ylab = 'Percent Survival',
             # title = 'FBL-high',
             #legend = 'bottom',
             legend = c(0.6, 0.95),
             legend.labs=c("IGF1R-AS1 high (N = 58)", "IGF1R-AS1 low (N = 441)"),
             ggtheme = theme_survminer()
  )

plot$plot + guides(color = guide_legend(nrow = 1))

ggsave('local_data/plots/igf1r_as1.pdf', plot = plot$plot,
       width = 85, height = 55, units = 'mm')

# *****************************************************************************
# get clinical plots, igf1r_as1, myc --------------------
# *****************************************************************************
library(survival)
library(survminer)

igf1r_myc_ar_clin <- igf1r_myc_ar %>%
  inner_join(os_data)

# pspline analysis shows that large values of igf1r_as1_activity
#   are linked to decreased survival, justifying optimal cutpoint visualization
coxph(Surv(times, status) ~ pspline(igf1r_as1_activity, df=4) + MYC_targets,
      data = igf1r_myc_ar_clin) %>% summary

surv_cut <- surv_cutpoint(data = igf1r_myc_ar_clin,
                          time = 'times',
                          event = 'status',
                          variables = c( 'igf1r_as1_activity',
                                         'MYC_targets'))

surv_cat <- surv_categorize(surv_cut) %>% as_tibble() %>%
  filter(MYC_targets == 'high')


fit <- survfit(Surv(times/ 30.44, status) ~ igf1r_as1_activity + MYC_targets,
               data = surv_cat)

ggsurvplot(fit, data = surv_cat, pval = T,
           risk.table = T)

plot <-
  ggsurvplot(fit, data = surv_cat,
             pval =T,
             pval.coord = c(0.05, 0.1),
             risk.table = T,
             palette = c('red3', 'blue3'),
             tables.theme = theme_cleantable(),
             censor.size = 4.25,
             # lwd = 1.5,
             legend.title="",
             legend.text = element_text(size = 2, color = "black", face = "bold"),
             xlab = 'Months',
             ylab = 'Percent Survival',
             # title = 'FBL-high',
             #legend = 'bottom',
             legend = c(0.6, 0.95),
             legend.labs=c("IGF1R-AS1 high, MYC high (N = 49)", "IGF1R-AS1 low, MYC high (N = 441)"),
             ggtheme = theme_survminer()
  )

plot$plot + guides(color = guide_legend(nrow = 1))

ggsave('local_data/plots/igf1r_as1_myc.pdf', plot = plot$plot,
       width = 85, height = 55, units = 'mm')

# *****************************************************************************
# get clinical plots --------------------
# *****************************************************************************
library(survival)
library(survminer)

igf1r_myc_ar_clin <- igf1r_myc_ar %>%
  filter(patient_id %in% adt_patients) %>%
  inner_join(adt_progression)

coxph(Surv(times, status) ~ AR_response + MYC_targets + igf1r_as1_activity,
      data = igf1r_myc_ar_clin) %>% summary

surv_cut <- surv_cutpoint(data = igf1r_myc_ar_clin,
                          time = 'times',
                          event = 'status',
                          variables = c( 'AR_response',
                                         'igf1r_as1_activity'))

surv_cat <- surv_categorize(surv_cut) %>% as_tibble() %>%
  filter((AR_response == 'low' & igf1r_as1_activity == 'high') | (AR_response == 'high' & igf1r_as1_activity == 'low'))

fit <- survfit(Surv(times, status) ~ AR_response + igf1r_as1_activity,
               data = surv_cat)

ggsurvplot(fit, data = surv_cat, pval = T)
