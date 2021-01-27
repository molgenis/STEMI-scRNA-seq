
# load lib
library('variancePartition')
library('edgeR')
library('BiocParallel')
library(RColorBrewer)
library(VennDiagram)

venn_gene_overlap <- function(genes_set1, genes_set2, genes_name1, genes_name2, output_loc){
  myCol <- brewer.pal(3, "Pastel2")[1:2]
  venn.diagram(x = list(genes_set1, genes_set2),
               main = paste(genes_name1, ' vs ', genes_name2, ' DE genes overlap', sep = ''),
               category.names = c(genes_name1, genes_name2),
               filename = paste(output_loc, genes_name1, 'vs', genes_name2, '.png', sep = ''),
               imagetype="png" ,
               height = 1000 , 
               width = 1000 , 
               resolution = 300,
               compression = "lzw",
               lwd = 2,
               lty = 'blank',
               fill = myCol,
               cex = .6,
               fontface = "bold",
               fontfamily = "sans",
               cat.cex = 0.6,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.pos = c(-27, 27),
               cat.dist = c(0.055, 0.055),
               cat.fontfamily = "sans")
}



# read result against Baseline
v2_mono_limma_full <- readRDS('/data/cardiology/differential_expression/limma/output/limmadream_monocyte_v2.rds')
v3_mono_limma_full <- readRDS('/data/cardiology/differential_expression/limma/output/limmadream_monocyte_v3.rds')

# grab the relevant tables from limma
v2_mono_limma_Baseline_t24h <- topTable(v2_mono_limma_full, coef = c('timepoint.finalt24h'))
v2_mono_limma_Baseline_t8w <- topTable(v2_mono_limma_full, coef = c('timepoint.finalt8w'))
v3_mono_limma_Baseline_t24h <- topTable(v3_mono_limma_full, coef = c('timepoint.finalt24h'))
v3_mono_limma_Baseline_t8w <- topTable(v3_mono_limma_full, coef = c('timepoint.finalt8w'))

# grab the relevant tables from MAST
v2_mono_mast_Baseline_t24h <- read.table('/data/cardiology/differential_expression/MAST/results/stemi_v2_paired_lores_lfc01minpct01ncountrna_20200707/rna/monocyteBaselinet24h.tsv', sep = '\t', header = T, row.names = 1)
v2_mono_mast_Baseline_t8w <- read.table('/data/cardiology/differential_expression/MAST/results/stemi_v2_paired_lores_lfc01minpct01ncountrna_20200707/rna/monocyteBaselinet8w.tsv', sep = '\t', header = T, row.names = 1)
v3_mono_mast_Baseline_t24h <- read.table('/data/cardiology/differential_expression/MAST/results/stemi_v3_paired_lores_lfc01minpct01ncountrna_20200707/rna/monocyteBaselinet24h.tsv', sep = '\t', header = T, row.names = 1)
v3_mono_mast_Baseline_t8w <- read.table('/data/cardiology/differential_expression/MAST/results/stemi_v3_paired_lores_lfc01minpct01ncountrna_20200707/rna/monocyteBaselinet8w.tsv', sep = '\t', header = T, row.names = 1)

# get the significant genes from limma
v2_mono_limma_Baseline_t24h_sigs <- rownames(v2_mono_limma_Baseline_t24h[v2_mono_limma_Baseline_t24h$adj.P.Val <= 0.05, ])
v2_mono_limma_Baseline_t8w_sigs <- rownames(v2_mono_limma_Baseline_t8w[v2_mono_limma_Baseline_t8w$adj.P.Val <= 0.05, ])
v3_mono_limma_Baseline_t24h_sigs <- rownames(v3_mono_limma_Baseline_t24h[v3_mono_limma_Baseline_t24h$adj.P.Val <= 0.05, ])
v3_mono_limma_Baseline_t8w_sigs <- rownames(v3_mono_limma_Baseline_t8w[v3_mono_limma_Baseline_t8w$adj.P.Val <= 0.05, ])

# get the significant genes from MAST
v2_mono_mast_Baseline_t24h_sigs <- rownames(v2_mono_mast_Baseline_t24h[v2_mono_mast_Baseline_t24h$p_val_adj <= 0.05, ])
v2_mono_mast_Baseline_t8w_sigs <- rownames(v2_mono_mast_Baseline_t8w[v2_mono_mast_Baseline_t8w$p_val_adj <= 0.05, ])
v3_mono_mast_Baseline_t24h_sigs <- rownames(v3_mono_mast_Baseline_t24h[v3_mono_mast_Baseline_t24h$p_val_adj <= 0.05, ])
v3_mono_mast_Baseline_t8w_sigs <- rownames(v3_mono_mast_Baseline_t8w[v3_mono_mast_Baseline_t8w$p_val_adj <= 0.05, ])

# output overlaps
venn_gene_overlap(v2_mono_limma_Baseline_t24h_sigs, v2_mono_mast_Baseline_t24h_sigs, 'v2 limma', 'v2 MAST', '/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/Baselinet24h_')
venn_gene_overlap(v3_mono_limma_Baseline_t24h_sigs, v3_mono_mast_Baseline_t24h_sigs, 'v3 limma', 'v3 MAST', '/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/Baselinet24h_')
venn_gene_overlap(v2_mono_limma_Baseline_t8w_sigs, v2_mono_mast_Baseline_t8w_sigs, 'v2 limma', 'v2 MAST', '/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/Baselinet8w_')
venn_gene_overlap(v3_mono_limma_Baseline_t8w_sigs, v3_mono_mast_Baseline_t8w_sigs, 'v3 limma', 'v3 MAST', '/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/Baselinet8w_')

# get overlaps
v2_mono_Baseline_t24h_sigs_overlaps <- intersect(v2_mono_limma_Baseline_t24h_sigs, v2_mono_mast_Baseline_t24h_sigs)
v2_mono_Baseline_t8w_sigs_overlaps <- intersect(v2_mono_limma_Baseline_t8w_sigs, v2_mono_mast_Baseline_t8w_sigs)
v3_mono_Baseline_t24h_sigs_overlaps <- intersect(v3_mono_limma_Baseline_t24h_sigs, v3_mono_mast_Baseline_t24h_sigs)
v3_mono_Baseline_t8w_sigs_overlaps <- intersect(v3_mono_limma_Baseline_t8w_sigs, v3_mono_mast_Baseline_t8w_sigs)

# plot lfc concordance
v2_mono_Baseline_t24h <- data.frame(limma=v2_mono_limma_Baseline_t24h[v2_mono_Baseline_t24h_sigs_overlaps, 'logFC'], mast=v2_mono_mast_Baseline_t24h[v2_mono_Baseline_t24h_sigs_overlaps, 'avg_logFC'])
ggplot(v2_mono_Baseline_t24h, aes(x=limma, y=mast)) + geom_point()
ggsave('/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/logFC_concordance_v2_mono_Baselinet24h.png', width = 5, height = 5)
v2_mono_Baseline_t8w <- data.frame(limma=v2_mono_limma_Baseline_t8w[v2_mono_Baseline_t8w_sigs_overlaps, 'logFC'], mast=v2_mono_mast_Baseline_t8w[v2_mono_Baseline_t8w_sigs_overlaps, 'avg_logFC'])
ggplot(v2_mono_Baseline_t8w, aes(x=limma, y=mast)) + geom_point()
ggsave('/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/logFC_concordance_v3_mono_Baselinet24h.png', width = 5, height = 5)
v3_mono_Baseline_t24h <- data.frame(limma=v3_mono_limma_Baseline_t24h[v3_mono_Baseline_t24h_sigs_overlaps, 'logFC'], mast=v3_mono_mast_Baseline_t24h[v3_mono_Baseline_t24h_sigs_overlaps, 'avg_logFC'])
ggplot(v3_mono_Baseline_t24h, aes(x=limma, y=mast)) + geom_point()
ggsave('/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/logFC_concordance_v2_mono_Baselinet8w.png', width = 5, height = 5)
v3_mono_Baseline_t8w <- data.frame(limma=v3_mono_limma_Baseline_t8w[v3_mono_Baseline_t8w_sigs_overlaps, 'logFC'], mast=v3_mono_mast_Baseline_t8w[v3_mono_Baseline_t8w_sigs_overlaps, 'avg_logFC'])
ggplot(v3_mono_Baseline_t8w, aes(x=limma, y=mast)) + geom_point()
ggsave('/data/cardiology/differential_expression/overlap/plots/limma_vs_mast/logFC_concordance_v3_mono_Baselinet8w.png', width = 5, height = 5)

