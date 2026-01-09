#### Jinxin Meng, 20240125, 20251220 ####
setwd('F:/project/20240731_PR_IBD_IPA_zhangyn/Code-Data-available/git/Figure8/')
pacman::p_load(tidyverse, ggpubr, openxlsx)
source('F:/code/R_func/calcu_difference.R')

#### Figure 8a ####
group <- read.delim('HeQ_2017.PRJEB15371.sample_group')

data <- read.delim('HeQ_2017.PRJEB15371.tsv') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('Control', 'CD')))

comparisons <- calcu_diff(data, rela_ab ~ group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00'),
         ylab = 'Relative abundance(%)', outlier.shape = NA,  title = 'HeQ_2017') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .2)) +
  stat_compare_means(comparisons = comparisons, label = 'p.label', 
                     label.y = .17, tip.length = .01, step.increase = .07, 
                     vjust = .6, size = 3) +
  theme(aspect.ratio = 1.6) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))
ggsave('8a.pdf', width = 4, height = 4)

rstatix::wilcox_test(data, rela_ab ~ group, p.adjust.method = 'fdr')

#### Figure 8b ####
group <- read.delim('SchirmerM_2018.PRJNA389280.sample_group') %>% 
  select(sample = run, group = group2)

data <- read.delim('SchirmerM_2018.PRJNA389280.tsv') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('CT', 'CD', 'UC')))

comparisons <- calcu_diff(data, rela_ab ~ group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00','#f08178'),
         ylab = 'Relative abundance(%)', outlier.shape = NA,
         title = 'SchirmerM_2018') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  stat_compare_means(comparisons = comparisons, label = 'p.label', 
                     label.y = .135, tip.length = .03, step.increase = .07, 
                     vjust = .6, size = 3) +
  theme(aspect.ratio = 1.2) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))
ggsave('8b.pdf', width = 4, height = 4)

rstatix::wilcox_test(data, rela_ab ~ group, p.adjust.method = 'BH') # FDR replace

#### Figure 8c ####
group <- read.delim('LloydPriceJ_2019.PRJNA398089.sample_group') %>% 
  select(sample = run, group = group2)

data <- read.delim('LloydPriceJ_2019.PRJNA398089.tsv') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('CT', 'CD', 'UC')))

comparisons <- calcu_diff(data, rela_ab ~ group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00','#f08178'),
         ylab = 'Relative abundance(%)', outlier.shape = NA,
         title = 'LloydPriceJ_2019') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .15)) +
  stat_compare_means(comparisons = comparisons, label = 'p.label', 
                     label.y = .14, tip.length = .012, step.increase = .07, 
                     vjust = .6, size = 3) +
  theme(aspect.ratio = 1.2) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))
ggsave('8c', width = 4, height = 4)

rstatix::wilcox_test(data, rela_ab ~ group, p.adjust.method = 'BH')

#### Figure 8d ####
group <- read.delim('YanQ_2023c.PRJEB67456.sample_group') %>% 
  select(sample = run, group)

data <- read.delim('YanQ_2023c.PRJEB67456.tsv') %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = 'sample') %>% 
  mutate(group = factor(group, c('healthy', 'CD', 'UC')))

comparisons <- calcu_diff(data, rela_ab ~ group) %>% 
  filter(pval < 0.05) %>% 
  pull(comparison) %>% 
  strsplit(split = '_vs_')

ggviolin(data, 'group', 'rela_ab', color = 'group', legend = 'none', xlab = '',
         size = 1, width = .7, palette = c('#69a7bc','#E69F00','#f08178'),
         ylab = 'Relative abundance(%)', outlier.shape = NA, title = 'YanQ_2023c') +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .5)) +
  stat_compare_means(comparisons = comparisons, label = 'p.label', method = 'wilcox',
                     label.y = .46, tip.length = .01, step.increase = .07, 
                     vjust = .6, size = 3) +
  theme(aspect.ratio = 1.2) +
  theme(plot.title = element_text(hjust = .5, face = 'bold'))
ggsave('8d.pdf', width = 4, height = 4)

rstatix::wilcox_test(data, rela_ab ~ group, p.adjust.method = 'fdr')

#### Figure 8e ####
data <- read.delim('data_PR_vs_IPA.txt')

ggscatter(data, 'abso_ab', y = 'content', color = '#ec6644', 
          add = 'reg.line', ylab = 'concentration of IPA (umol/g)', 
          xlab = 'log10 (genomic DNA copies/mg feces) of P. russellii',
          cor.coef = T, cor.method = 'spearman', cor.coef.coord = c(1, 0.015)) +
  ggpmisc::stat_poly_eq(aes(label = paste(after_stat(eq.label), 
                                          after_stat(adj.rr.label), 
                                          after_stat(p.value.label), sep = '~~~~')), 
                        formula = y ~ x) +
  theme(aspect.ratio = 1)

fit <- lm(abso_ab ~ content, data)
summary(fit)
