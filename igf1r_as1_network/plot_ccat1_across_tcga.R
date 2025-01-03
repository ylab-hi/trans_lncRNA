source('~/theme_manuscript.R')
library(tidyverse)

# *****************************************************************************
# across TCGA --------------------
# *****************************************************************************

# data downloaded from xena browser
ccat1 <- read_tsv('data/ccat1_pancancer/ccat1_pan_cancer_xena_browser.tsv',
                  col_names = c('sample', 'samples',
                                'ccat1', 'cancer_cohort'), skip = 1)

ccat1 <- ccat1 %>%
  filter(!is.na(ccat1)) %>%
  filter(!is.na(cancer_cohort))

cohort_order <- reorder(ccat1$cancer_cohort, ccat1$ccat1, FUN = median) %>% 
  attr('scores') %>% sort

ccat1 %>%
  ggplot(aes(x = reorder(cancer_cohort, ccat1, FUN = median), 
             y = ccat1,
             color = cancer_cohort)) +
  ggbeeswarm::geom_quasirandom(size = 0.25) +
  geom_boxplot(aes(fill = cancer_cohort), color = 'black', outlier.shape = NA,
               alpha = 0.25, linewidth = 0.25) +
  theme_manuscript_small() +
  scale_color_manual(values = viridisLite::turbo(33)[sapply(sort(names(cohort_order)), \(x) which(names(cohort_order) == x))]) +
  scale_fill_manual(values = viridisLite::turbo(33)[sapply(sort(names(cohort_order)), \(x) which(names(cohort_order) == x))]) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position = 'none') +
  labs(x = 'TCGA Cohort', y = 'CCAT1 Expression log2(RPKM+1)')

ggsave('plots/ccat1_tcga2.pdf',
       height = 50, width = 80, units = 'mm')

