#### Jinxin Meng, 20240801, 20250109 ####
setwd('F:/project/20240731_PR_IBD_IPA_ZhangYN/Code-Data-available/git/Figure3/')
pacman::p_load(tidyverse, openxlsx, ggpubr)
source('/Code/R_func/profile_process.R')

#### Figures 3b, 3d and S3b ####
group <- read.xlsx('RNASeq_data.xlsx', sheet = 'group')

rc <- read.xlsx('RNASeq_data.xlsx', sheet = 'rc') %>% 
  column_to_rownames('name')

fpkm <- read.xlsx("RNASeq_data.xlsx", sheet = 'fpkm') %>% 
  column_to_rownames('name')

markers <- list(
  aISC = c("Lgr5","Ascl2",'Slc12a2'),
  rISC = c("Hopx","Bmi1",'Lrig1'),
  PC = c("Ang4","Lyz1",'Agr2'),
  GC = c("Muc2","Ccl9"),
  TC = c("Dclk1","Cd24a"),
  EEC = c("Chga",'Neurod1')
)

marker_levels <- c("aISC", "rISC", "PC", "GC", "TC", "EEC")

# D0
group_x <- filter(group, sample %in% samples[[1]]) %>% 
  mutate(group = factor(group, c('PBS_D0', 'IPA_D0')))

fpkm_x <- dplyr::select(fpkm, all_of(group_x$sample)) %>%
  filter(rownames(.) %in% unlist(markers)) %>% 
  profile_transSqrt()

row_data <- markers %>% 
  map2_dfr(names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  arrange(match(name, rownames(fpkm_x))) %>% 
  mutate(class = factor(class, marker_levels))

col_data <-  group_x

row_split <- row_data$class
col_split <- col_data$group

pdf("fpkm_D0.pdf", width = 5, height = 5)
ComplexHeatmap::pheatmap(
  fpkm_x, scale = "row",
  color = colorRampPalette(c("#307cc0", "white", "#e43589"))(100),
  split = row_split, column_split = col_split,
  cluster_row_slices = F, cluster_column_slices = F, 
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 8),
  row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
  treeheight_col = 10, treeheight_row = 10,
  cellheight = 12, cellwidth = 12, fontsize = 10, 
  border_color = "white", border_gp = grid::gpar(col = "black"), 
  heatmap_legend_param = list(title = 'Scale FPKM') )
dev.off()

# D3
group_x <- filter(group, sample %in% samples[[2]]) %>% 
  mutate(group = factor(group, c('PBS_D3', 'IPA_D3')))

fpkm_x <- dplyr::select(fpkm, all_of(group_x$sample)) %>%
  filter(rownames(.) %in% unlist(markers))

row_data <- markers %>% 
  map2_dfr(names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  arrange(match(name, rownames(fpkm_x))) %>% 
  mutate(class = factor(class, marker_levels))

col_data <-  group_x

row_split <- row_data$class
col_split <- col_data$group

pdf("fpkm_D3.pdf", width = 5, height = 5)
ComplexHeatmap::pheatmap(
  fpkm_x, scale = "row",
  color = colorRampPalette(c("#307cc0", "white", "#e43589"))(100),
  split = row_split, column_split = col_split,
  cluster_row_slices = F, cluster_column_slices = F, 
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 8),
  row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
  treeheight_col = 10, treeheight_row = 10,
  cellheight = 12, cellwidth = 12, fontsize = 10, 
  border_color = "white", border_gp = grid::gpar(col = "black"), 
  heatmap_legend_param = list(title = 'Scale FPKM') )
dev.off()

# D7
group_x <- filter(group, sample %in% samples[[3]]) %>% 
  mutate(group = factor(group, c('PBS_D7', 'IPA_D7')))

fpkm_x <- dplyr::select(fpkm, all_of(group_x$sample)) %>%
  filter(rownames(.) %in% unlist(markers))

row_data <- markers %>% 
  map2_dfr(names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  arrange(match(name, rownames(fpkm_x))) %>% 
  mutate(class = factor(class, marker_levels))

col_data <-  group_x

row_split <- row_data$class
col_split <- col_data$group

pdf("fpkm_D7.pdf", width = 5, height = 5)
ComplexHeatmap::pheatmap(
  fpkm_x, scale = "row",
  color = colorRampPalette(c("#307cc0", "white", "#e43589"))(100),
  split = row_split, column_split = col_split,
  cluster_row_slices = F, cluster_column_slices = F, 
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 8), column_title_gp = grid::gpar(fontsize = 8),
  row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
  treeheight_col = 10, treeheight_row = 10,
  cellheight = 12, cellwidth = 12, fontsize = 10, 
  border_color = "white", border_gp = grid::gpar(col = "black"), 
  heatmap_legend_param = list(title = 'Scale FPKM') )
dev.off()

#### Figure 3a, 3c ####
gene_info <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, rownames(fpkm), 
                                   keytype = 'SYMBOL', columns = c('ENTREZID'))

fpkm_x <- profile_transLOG2(fpkm)

metadata <- data.frame(row.names = colnames(fpkm_x)) %>%
  mutate(group = group$group[match(rownames(.), group$sample)],
         group = factor(group))

grps <- metadata$group
design <- model.matrix(~ 0 + grps)
colnames(design) <- gsub('grps', '', colnames(design))

comparisons <- list(c('IPA_D0', 'PBS_D0'), 
                    c('IPA_D3', 'PBS_D3'), 
                    c('IPA_D7', 'PBS_D7'))

contrasts <- limma::makeContrasts(
  contrasts = map_vec(comparisons, ~ paste(.x, collapse = '-')), levels = design)

fit <- limma::lmFit(fpkm_x, design)
fit_contrast <- limma::contrasts.fit(fit, contrasts)
fit_contrast <- limma::eBayes(fit_contrast)

diffs <- map(comparisons, \(x) 
             limma::topTable(fit_contrast, adjust = "BH", number = Inf, 
                      coef = paste0(x, collapse = '-') ) %>% 
               rownames_to_column('name') %>% 
               mutate(enriched = ifelse(logFC > .5 & P.Value < 0.05, x[1], 
                                        ifelse(logFC < -.5 & P.Value < 0.05, 
                                               x[2], "none"))) %>% 
                 left_join(gene_info, by = c('name' = 'SYMBOL'))) %>% 
  set_names(map_vec(comparisons, ~ paste(.x, collapse = '_vs_')) )

eKEGGs <- map2(diffs, map_vec(comparisons, ~ .x[1]), \(x, y)
              filter(x, enriched == y) %>% 
                pull(ENTREZID) %>% 
                na.omit %>% 
                clusterProfiler::enrichKEGG(
                  organism = 'mmu', pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                data.frame())
write.xlsx(eKEGGs, 'eKEGGs_up.xlsx')

eGOs <- map2(diffs, map_vec(comparisons, ~ .x[1]), \(x, y)
            filter(x, enriched == y) %>% 
              pull(name) %>% 
              clusterProfiler::enrichGO(
                org.Mm.eg.db::org.Mm.eg.db, keyType = "SYMBOL", ont = 'ALL',
                pvalueCutoff = 1, qvalueCutoff = 1) %>% 
              data.frame())
write.xlsx(eGOs, 'eGOs_up.xlsx')

#### Figure 3j, S3d ####
# read
path <- list.files('scRNASeq_data/', '1', full.names = T)
name <- list.files('scRNASeq_data/', '1')

seurat_list <- map2(path, name, \(x, y) {
  message(paste0('Processing: ', y)) 
  count <- Read10X(x)
  CreateSeuratObject(counts = count, project = y, min.cells = 5, min.features = 100)
})

seurat <- merge(seurat_list[[1]], seurat_list[-1], add.cell.ids = name) %>% 
  JoinLayers()

# table(seurat@meta.data$orig.ident)
# IPA_1 Veh_1 
# 10049  8865

# qc 
mit_genes <- rownames(seurat)[grep('^mt-', rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat, features = mit_genes, col.name = 'percent_mit')
rib_genes <- rownames(seurat)[grep('^Rp[sl]', rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat, features = rib_genes, col.name = 'percent_rib')
hb_genes <- rownames(seurat)[grep('^Hb[^(p)]', rownames(seurat), ignore.case = T)]
seurat <- PercentageFeatureSet(seurat,  features = hb_genes, col.name = 'percent_hb')

p1 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('nFeature_RNA', 'nCount_RNA'), 
              pt.size = 0, ncol = 2) + 
  NoLegend()

p2 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('percent_mit', 'percent_rib', 'percent_hb'), 
              pt.size = 0, ncol = 3) + 
  scale_y_continuous(breaks = seq(0, 100, 5)) + 
  NoLegend()

cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
ggsave('qc.before.pdf', width = 6, height = 10)

# 过滤, 基因表达太多了，可能就是异常点，存在双细胞混一起的可能
cells <- purrr::reduce(list(WhichCells(seurat, expression = percent_mit < 20), 
                            WhichCells(seurat, expression = nFeature_RNA > 100),
                            WhichCells(seurat, expression = nFeature_RNA < 5000),
                            WhichCells(seurat, expression = percent_hb < 1)),
                       \(x, y) base::intersect(x, y))
seurat <- subset(seurat, cells = cells)

p1 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('nFeature_RNA', 'nCount_RNA'), 
              pt.size = 0, ncol = 2) + 
  NoLegend()

p2 <- VlnPlot(seurat, group.by = 'orig.ident', features = c('percent_mit', 'percent_rib', 'percent_hb'), 
              pt.size = 0, ncol = 3) + 
  NoLegend()

cowplot::plot_grid(p1, p2, ncol = 1, align = 'v')
ggsave('qc.after.pdf', width = 6, height = 10)

# dim(seurat)
# [1] 17396  8206

# PCA and scale
seurat <- NormalizeData(seurat, normalization.method = 'LogNormalize', scale.factor = 1e4, verbose = T) # 标准化
seurat <- FindVariableFeatures(seurat, selection.method = 'vst', nfeatures = 2000) # 筛选高变基因
VariableFeaturePlot(seurat, cols = c('grey','red'))
seurat <- ScaleData(seurat) # 数据归一化，这步骤有一步回归的分析，去除噪声
seurat <- RunPCA(seurat, features = VariableFeatures(seurat)) # PCA线性降维分析
ElbowPlot(seurat, ndims = 30) # 流石图协助选择PC维度 10
seurat <- RunHarmony(seurat, 'orig.ident') # Harmony去批次
seurat <- RunUMAP(seurat, dims = 1:10, reduction = 'harmony')
seurat <- RunTSNE(seurat, dims = 1:10, reduction = 'harmony')
seurat <- FindNeighbors(seurat, reduction = 'harmony', dims = 1:10)
seurat <- FindClusters(seurat, resolution = c(seq(0, 1.9, .2)))

map2(paste0('RNA_snn_res.', seq(0, 1.9, .2)), paste0('SNN: ', seq(0, 1.9, .2)),
     \(x, y) DimPlot(seurat, reduction = 'umap', group.by = x, label = T) +
       ggtitle(y) + 
       theme(aspect.ratio = 1) + 
       guides(color = 'none')) %>% 
  plot_grid(plotlist = ., nrow = 2)
ggsave('umap.SNN.dimplot.jpg', width = 20, height = 8)

# annotation
markers <- list(
  Stem_cells = c('Lgr5','Hopx','Bmi1','Ascl2','Olfm4','Smoc2'),
  Paneth_cells = c('Ang4','Lyz1','Defa17','Defa21','Defa22','Defa24','Gm14851','Defa30'),
  Goblet_cells = c('Muc2','Ccl9','Clca1','Tff3','Agr2','Fcgbp','Zg16'),
  Enteroendocrine_cells = c('Chga','Neurog3','Chgb','Tac1','Tph1'),
  Tuft_cells = c('Dclk1','Trpm5'),
  Enterocyte = c('Aldob', 'Apoa1', 'Apoa4', 'Gsta1','Fabp1', 'Prap1'),
  T_cells = c('Cd3g','Cd4','Cd8a','Tcrg-C1','Tcrg-C2','Tcrg-C4','Gzma','Gzmb'),
  B_cells = c('Cd79a','Jchain','Igha','Igkc'))

cell_col <- c(Stem_cells = '#8dd3c7', Paneth_cells = '#ffed6f', Goblet_cells = '#bebada',
              Enteroendocrine_cells = '#fb8072', Tuft_cells = '#80b1d3',
              Enterocyte = '#fdb462', T_cells = '#b3de69', B_cells = '#bc80bd')

seurat <- SetIdent(seurat, value = 'RNA_snn_res.1.4')
DimPlot(seurat, reduction = 'tsne', label = T, alpha = 1, raster = F) + 
  NoLegend() + 
  theme(aspect.ratio = 1)
ggsave('tsne.SNN.1.4.dimplot.pdf', width = 6, height = 6)

# 鉴定每个cluster中的marker
allmarkers <- FindAllMarkers(seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(allmarkers, 'cluster.allmarkers.rds')
write.xlsx(allmarkers, 'cluster.allmarkers.xlsx')

seurat <- RenameIdents(seurat, 
                       '0' = 'Goblet_cells','1' = 'Paneth_cells','2' = 'Enterocyte',
                       '3' = 'T_cells','4' = 'Enterocyte','5' = 'Enterocyte',
                       '6' = 'Goblet_cells','7' = 'Enterocyte', '8' = 'Goblet_cells',
                       '9' = 'Stem_cells','10' = 'Enterocyte','11' = 'Goblet_cells', 
                       '12' = 'Enterocyte','13' = 'Enterocyte','14' = 'Paneth_cells',
                       '15' = 'Enterocyte', '16' = 'Enterocyte','17' = 'Enterocyte',
                       '18' = 'Paneth_cells','19' = 'Enterocyte','20' = 'Goblet_cells',
                       '21' = 'T_cells','22' = 'Goblet_cells','23' = 'Tuft_cells',
                       '24' = 'Enteroendocrine_cells','25' = 'Goblet_cells','26' = 'T_cells',
                       '27' = 'B_cells','28' = 'Enterocyte')

DimPlot(seurat, reduction = 'tsne', label = T, alpha = 1, cols = cell_col) + 
  theme(aspect.ratio = 1) + NoLegend()
ggsave('tsne.marker.dimplot.pdf', width = 6, height = 6)

# find marker
allmarkers <- FindAllMarkers(seurat, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(allmarkers, 'cell.allmarkers.rds')
write.xlsx(allmarkers, 'cell.allmarkers.xlsx')

# marker violin plot
markers <- list(
  Stem_cells = c('Lgr5','Bmi1','Ascl2','Smoc2'),
  Paneth_cells = c('Ang4','Lyz1'),
  Goblet_cells = c('Muc2','Clca1','Tff3','Agr2'),
  Enteroendocrine_cells = c('Chga','Tac1','Tph1'),
  Tuft_cells = c('Dclk1','Trpm5'),
  Enterocyte = c('Aldob', 'Gsta1', 'Prap1'),
  T_cells = c('Cd3g','Cd8a','Gzma'),
  B_cells = c('Cd79a','Igkc'))

top_markers <- map2_df(markers, names(markers), \(x, y) 
                       rbind(filter(allmarkers, cluster == y & gene %in% x),
                             filter(allmarkers, cluster == y & !gene %in% x) %>% 
                               head(n = 4 - sum(allmarkers$cluster == y & allmarkers$gene %in% x)) ) %>% 
                         arrange(desc(avg_log2FC)) )

plot_data <- FetchData(seurat, vars = unique(top_markers$gene)) %>% 
  add_column(cell = seurat@active.ident) %>% 
  gather(key = 'gene', value = 'value', -cell) %>% 
  mutate(gene = factor(gene, unique(top_markers$gene)),
         cell = factor(cell, names(markers)))

ggplot(plot_data, aes(cell, value), color = factor(cell)) +
  geom_violin(aes(fill = cell), scale ='width') +
  facet_grid(gene ~ ., scales = 'free_y') +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cell_col) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1,
                                          vjust = NULL, color = 'black', size = 14),
        axis.text.y.left = element_blank(),
        legend.position = 'none',
        panel.spacing.y = unit(0, 'cm'),
        strip.text.y = element_text(angle = 0, size = 14, hjust = 0, face = 'italic'),
        strip.background.y = element_blank())
ggsave('Figure S3d.pdf', width = 5, height = 9)

#### Figure 3k ####
data <- seurat@meta.data %>% 
  dplyr::select(group = orig.ident) %>% 
  add_column(cells = Idents(seurat)) %>% 
  group_by(group, cells) %>% 
  summarise(value = n()) %>% 
  group_modify(~.x %>% mutate(prec = value / sum(value) * 100)) %>% 
  ungroup() %>% 
  mutate(group = factor(group, c('Veh_1', 'IPA_1')))

ggbarplot(data, 'group', 'prec', fill = 'cells', xlab = '', 
          ylab = 'Cell precentage (%)', legend = 'right', palette = cell_col)
ggsave('Figure 3k.pdf', width = 4, height = 4)

#### Figure 3l ####
cell_name <- levels(seurat@active.ident)
cell_col <- c(Stem_cells = '#8dd3c7', Paneth_cells = '#ffed6f', Goblet_cells = '#bebada',
              Enteroendocrine_cells = '#fb8072', Tuft_cells = '#80b1d3',
              Enterocyte = '#fdb462', T_cells = '#b3de69', B_cells = '#bc80bd')

gene_info <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, 
                                   AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db), 
                                   keytype = "ENTREZID", columns = "SYMBOL")

diffs <- map(cell_name, \(x) {
  seurat_x <- subset(seurat, idents = x)
  FindMarkers(seurat_x, ident.1 = 'IPA_1', ident.2 = 'Veh_1', 
              group.by = 'orig.ident') %>% 
    rownames_to_column('name') %>% 
    mutate(enriched = ifelse(avg_log2FC > .1 & p_val < .05, 'up', 
                             ifelse(avg_log2FC < -.1 & p_val < .05, 
                                    'down', 'none'))) }) %>% set_names(cell_name)

eKEGGs <- map(diffs, \(x)
              filter(x, p_val < .05 & abs(avg_log2FC) > .25) %>% 
                left_join(gene_info, by = c('name' = 'SYMBOL')) %>% 
                pull(ENTREZID) %>% 
                enrichKEGG(organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1) %>% 
                data.frame() )
write.xlsx(eKEGGs, 'cells.diff.Gene.eKEGG.xlsx')

eGOs <- map(diffs, \(x)
            filter(x, p_val < .05 & abs(avg_log2FC) > .25) %>% 
              pull(name) %>% 
              enrichGO(org.Mm.eg.db, keyType = 'SYMBOL', pvalueCutoff = 1, 
                       qvalueCutoff = 1, ont = 'ALL') ) 
write.xlsx(eGOs, 'cells.diff.Gene.eGO.xlsx')

terms <- c('GO:0042180','GO:0006631','GO:0042181','GO:0010565','GO:0050678',
           'GO:0007219','GO:0019827','GO:0045927')

eGOs_data <- map_dfr(cell_name, ~ 
                       read.xlsx('cells.diff.Gene.eGO.xlsx', sheet = .x) %>% 
                       filter(ID %in% terms) %>% 
                       dplyr::select(Description, FoldEnrichment, qvalue) %>%
                       add_column(cell = .x))

terms <- c('mmu04110','mmu03030','mmu00650','mmu04152','mmu03030')

eKEGGs_data <- map_dfr(cell_name, ~ 
                         read.xlsx('cells.diff.Gene.eKEGG.xlsx', sheet = .x) %>% 
                         filter(ID %in% terms) %>% 
                         dplyr::select(Description, FoldEnrichment, qvalue) %>%
                         add_column(cell = .x)) 

plot_data <- rbind(eGOs_data, eKEGGs_data) %>% 
  filter(qvalue < 0.2)

ggscatter(plot_data, 'cell', 'Description', fill = 'qvalue', shape = 21,
          size = 'FoldEnrichment', xlab = '', ylab = '', legend = 'right',
          x.text.angle = 30) +
  scale_fill_viridis_c(begin = .5, end = .9) +
  scale_size_continuous(range = c(3, 6)) +
  theme(aspect.ratio = 2,
        panel.grid.major = element_line(linewidth = .5, color = 'grey88'))

ggsave('Figure 3l.pdf', width = 7, height = 5)


#### Figure 3m ####
# fold
filter(data, !cells %in% c('T_cells','B_cells')) %>% 
  select(-value) %>% 
  spread('group', 'prec') %>% 
  mutate(prop = IPA_1 / Veh_1,
         prop = round(prop, 2))

# P values
filter(data, !cells %in% c('T_cells','B_cells')) %>% 
  select(-prec) %>% 
  spread('group', 'value') %>% 
  mutate(IPA_other = 4825 - IPA_1,
         Veh_total = 3381 - Veh_1) %>% 
  group_by(cells) %>% 
  group_modify(~ fisher.test(matrix(unlist(.x), 2))$p.value %>% 
                 data.frame(pval = .)) %>% 
  ungroup %>% 
  mutate(padj = p.adjust(pval, method = 'BH') %>% format(scientific = T, digit = 3))

# plot
filter(data, !cells %in% c('T_cells','B_cells')) %>%
  mutate(group = factor(group, c('IPA_1', 'Veh_1')),
         cells = factor(cells, c('Enterocyte','Goblet_cells','Paneth_cells',
                                 'Stem_cells','Tuft_cells','Enteroendocrine_cells'))) %>% 
  ggbarplot('cells', 'prec', fill = 'group', palette = c('#87d6f5', '#a09fa4'),
            xlab = '', ylab = 'Percentage (%)') +
  annotate('text', x = 1, y = 85, label = '1.34-Fold ***', size = 4) +
  annotate('text', x = 2, y = 60, label = '0.76-Fold ***', size = 4) +
  annotate('text', x = 3, y = 30, label = '0.92-Fold ns', size = 4) +
  annotate('text', x = 4, y = 9, label = '1.61-Fold ***', size = 4) +
  annotate('text', x = 5, y = 6, label = '0.96-Fold ns', size = 4) +
  annotate('text', x = 6, y = 3, label = '0.91-Fold ns', size = 4) +
  theme(aspect.ratio = 3/3.5)
ggsave('Figure 3m.pdf', width = 5, height = 4)

#### Figure 3n, 3o ####
library(ggunchull)
library(ggforce)

cell_col <- c(Stem_cells = '#8dd3c7', Paneth_cells = '#ffed6f', Goblet_cells = '#bebada',
              Enteroendocrine_cells = '#fb8072', Tuft_cells = '#80b1d3',
              Enterocyte = '#fdb462', T_cells = '#b3de69', B_cells = '#bc80bd')

data <- FetchData(seurat, vars = c('Hmgcs2')) %>% 
  add_column(cells = seurat@active.ident) %>% 
  rownames_to_column('name') %>% 
  left_join(seurat@reductions$tsne@cell.embeddings %>% 
              data.frame %>% 
              rownames_to_column('name'),
            by = 'name') %>% 
  mutate(.Hmgcs2 = ifelse(Hmgcs2 > 0, '1', '0'))

neg_data <- filter(data, .Hmgcs2 == 0)
pos_data <- filter(data, .Hmgcs2 == 1)

ggplot() +
  stat_unchull(aes(x = tSNE_1, y = tSNE_2, color = cells), data = data,
               alpha = 0, linetype = 'longdash', linewidth = .6, 
               delta = 2, th = 2) +
  scale_color_manual(values = cell_col) +
  geom_point(aes(tSNE_1, tSNE_2), data = neg_data, size = .3, color = 'grey80') +
  geom_point(aes(tSNE_1, tSNE_2, color = cells), data = pos_data, size = .3) +
  theme_pubr() +
  theme(aspect.ratio = 1)
ggsave('Figure 3n.pdf', width = 5.2, height = 5)

# 整体情况
data <- FetchData(seurat, vars = c('Hmgcs2')) %>% 
  add_column(cells = seurat@active.ident,
             group = seurat@meta.data$orig.ident) %>% 
  rownames_to_column('name') %>% 
  mutate(value = ifelse(Hmgcs2 > 0, 'pos', 'neg')) %>% 
  count(group, value) %>% 
  group_by(group) %>% 
  group_modify(~ .x %>% mutate(total = sum(n))) %>% 
  mutate(perc = n / total * 100) %>% 
  filter(value == 'pos')

44.6 / 30.5
fisher.test(matrix(c(2154, 2671, 1032, 2349), nrow = 2))

ggbarplot(data, 'value', 'perc', fill = 'group', xlab = '', 
          palette = c('#87d6f5','#a09fa4')) +
  scale_y_continuous(expand = c(.01, .01)) +
  annotate('text', x = 1, y = 60, label = '44.6%', size = 4) +
  annotate('text', x = 1, y = 20, label = '30.5%', size = 4) +
  annotate('text', x = 1, y = 70, label = '1.46-Fold ***', size = 4) +
  theme(aspect.ratio = 5)
ggsave('Figure 3n.pdf', width = 2, height = 6)

data <- FetchData(seurat, vars = c('Hmgcs2')) %>% 
  add_column(cells = seurat@active.ident) %>% 
  mutate(value = ifelse(Hmgcs2 > 0, '+', '-')) %>% 
  count(cells, value) %>% 
  group_by(cells) %>% 
  group_modify(~ mutate(.x, pct = n / sum(n) * 100)) %>% 
  filter(value == '+')

ggbarplot(data, 'cells', 'pct', fill = 'cells', legend = 'none', xlab = '', width = .6,
          ylab = 'Precentage (%)', sort.val = 'asc', sort.by.groups = F, rotate = T,
          title = 'Expression pattern of Hmgcs2 in CBCs and other cells',
          palette = cell_col) +
  theme(aspect.ratio = 4/3)
ggsave('Figure 3o.pdf', width = 6, height = 3.6)

#### Figure 3p ####
plot_data <- FetchData(seurat, vars = 'Hmgcs2') %>% 
  add_column(cell = Idents(seurat), group = seurat$orig.ident) %>% 
  mutate(value = ifelse(Hmgcs2 > 0, '+', '-')) %>% 
  count(cell, group, value) %>% 
  group_by(cell, group) %>% 
  group_modify(~ mutate(.x, pct = n / sum(n))) %>% 
  ungroup %>% 
  filter(value == '+') %>% 
  select(-value, -n) %>% 
  spread('group', 'pct') %>% 
  filter(!grepl('T|B', cell))

# fold
plot_data %>% 
  mutate(perc = IPA_1 / Veh_1,
         perc = round(perc, 2))

# P-value
FetchData(seurat, vars = 'Hmgcs2') %>% 
  add_column(cell = Idents(seurat), group = seurat$orig.ident) %>% 
  group_by(cell) %>% 
  rstatix::wilcox_test(Hmgcs2 ~ group, detailed = T) %>% 
  mutate(padj = p.adjust(p, method = 'BH'),
         padj = format(padj, scientific = T, digits = 3))

gather(plot_data, 'group', 'value', -cell) %>% 
  mutate(group = factor(group, c('Veh_1', 'IPA_1')),
         cell = factor(cell, c('Stem_cells','Enterocyte','Paneth_cells',
                               'Enteroendocrine_cells', 'Goblet_cells'))) %>% 
  ggbarplot('cell', 'value', fill = 'group', xlab = '', 
            palette = c('#a09fa4','#87d6f5'),
            position = position_dodge(width = .85), 
            ylab = 'Percentage of Hmgcs2+ cell') +
  scale_y_continuous(expand = c(.01, .01)) +
  annotate('text', x = 1, y = .9, label = '0.97-Fold\nFDR=1.70e-05', size = 4) +
  annotate('text', x = 2, y = .6, label = '1.15-Fold\nFDR=1.16e-04', size = 4) +
  annotate('text', x = 3, y = .25, label = '1.67-Fold\nFDR=6.20e-05', size = 4) +
  annotate('text', x = 4, y = .3, label = '1.41-Fold\nFDR=0.341', size = 4) +
  annotate('text', x = 5, y = .25, label = '1.31-Fold\nFDR=8.82e-03', size = 4) +
  theme(aspect.ratio = 3/5)
ggsave('Figure 3p.pdf', width = 6, height = 3.6)
