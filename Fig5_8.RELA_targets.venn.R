# Figures 5H, S5O, 8E
library(chipenrich)
library(eulerr)
library(ggplot2)

dir = "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/RNAseq/RSEM/"
setwd(dir)

################## 786-O parental ##################
# Link peaks to genes to find RELA target genes
# RELA ChIP
RELA_WT_bed <- "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/ChIPseq/PBRM1/q0.01/RCC523.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC501.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed"
RELA_WT_results <- chipenrich(peaks = RELA_WT_bed,
                              genome = 'hg19',
                              genesets = c("hallmark"),
                              locusdef = "nearest_tss",
                              qc_plots = FALSE,
                              out_name = NULL,
                              n_cores = 1)
results <- RELA_WT_results$results
print(results[1:5, 1:5])

# shRELA RNAseq
RELA_WT_expression <- read.delim("Deseq2_O786_WT_shRELA.txt")
RELA_WT_expression <- RELA_WT_expression$GeneSymbol[RELA_WT_expression$pvalue < 0.05 & 
                                                      RELA_WT_expression$log2FoldChange < -0.5]

# Combine RELA ChIPseq and shRELA RNAseq
RELA_WT_expression_chip <-  intersect(RELA_WT_expression, RELA_WT_results$peaks$gene_symbol)


# Plot Venn diagrams
pdf("PBRM1_WT_venn.pdf")
RELA_WT_euler <- eulerr::euler(combinations = list(
  RELA_WT = RELA_WT_expression_chip,
  PBRM1_KO_down = PBRM1_KO_expression_down,
  PBRM1_KO_up = PBRM1_KO_expression_up))

plot(RELA_WT_euler, quantities = list(type = c("percent", "counts")))
dev.off()

PBRM1_RELA_WT_intersect <- intersect(RELA_WT_expression_chip, PBRM1_KO_expression_up)
length(PBRM1_RELA_WT_intersect)
length(RELA_WT_expression_chip)
#41/259=0.16

################## 786-O PBRM1 KO ##################
# RELA ChIP
RELA_KO_bed <- "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/ChIPseq/PBRM1/q0.01/RCC524.merged.bam.sort.bam.mapQ.bam.bed_vs_RCC502.merged.bam.sort.bam.mapQ.bam.bed_peaks.bed"
RELA_KO_results <- chipenrich(peaks = RELA_KO_bed,
                              genome = 'hg19',
                              genesets = c("hallmark"),
                              locusdef = "nearest_tss",
                              qc_plots = FALSE,
                              out_name = NULL,
                              n_cores = 1)

results <- RELA_KO_results$results
print(results[1:5, 1:5])

# shRELA RNAseq
RELA_KO_expression <- read.delim("Deseq2_O786_KO_shRELA.txt")
RELA_KO_expression <- RELA_KO_expression$GeneSymbol[RELA_KO_expression$pvalue < 0.05 & 
                                                      RELA_KO_expression$log2FoldChange < -0.5]
RELA_KO_expression_chip <- intersect(RELA_KO_expression, RELA_KO_results$peaks$gene_symbol)

# Combine RELA ChIPseq and shRELA RNAseq
PBRM1_KO_expression <- read.delim("Deseq2_786O_PB_KO_noIFNg.txt")
PBRM1_KO_expression_up <- PBRM1_KO_expression$GeneSymbol[PBRM1_KO_expression$pvalue < 0.05 & 
                                                           PBRM1_KO_expression$log2FoldChange > 0.5]
PBRM1_KO_expression_down <- PBRM1_KO_expression$GeneSymbol[PBRM1_KO_expression$pvalue < 0.05 & 
                                                             PBRM1_KO_expression$log2FoldChange < -0.5]



PBRM1_RELA_KO_intersect <- intersect(RELA_KO_expression_chip, PBRM1_KO_expression_up)

length(PBRM1_RELA_KO_intersect)
length(RELA_KO_expression_chip)
#113/371=0.30


# Plot Venn diagrams
pdf("PBRM1_KO_venn.pdf")
RELA_KO_euler = eulerr::euler(combinations = list(
  RELA_KO = RELA_KO_expression_chip,
  PBRM1_KO_down = PBRM1_KO_expression_down,
  PBRM1_KO_up = PBRM1_KO_expression_up
  )
)
plot(RELA_KO_euler, quantities = list(type = c("percent", "counts")))
dev.off()

################### PB restore ##################
PBRM1_restore_expression <- read.delim("Deseq2_786OKO_PBrestore.DMSO.2.txt")
PBRM1_restore_expression_down <- PBRM1_restore_expression$GeneSymbol[PBRM1_restore_expression$pvalue < 0.05 &
                                                                       PBRM1_restore_expression$log2FoldChange < -0.3]
RELA_KO_restore_euler <- eulerr::euler(combinations = list(RELA_KO_targets = RELA_KO_expression_chip, 
                                                           PBRM1_restore = PBRM1_restore_expression_down)
)

KO_PB_restore_intersect <- intersect(RELA_KO_expression_chip, PBRM1_restore_expression_down)
pdf("PBRM1_KO_restore_venn.pdf")
plot(RELA_KO_restore_euler, 
     quantities = list(type = c("percent", "counts")))

dev.off()  

################### BTZ treatment ##################
BTZ_expression <- read.delim("Deseq2_786O_KO_DMSOvsBTZ.txt")
BTZ_expression_up <- BTZ_expression$GeneSymbol[BTZ_expression$pvalue < 0.05 &
                                                BTZ_expression$log2FoldChange > 0.5]
BTZ_expression_down <- BTZ_expression$GeneSymbol[BTZ_expression$pvalue <
                                                  0.05 & BTZ_expression$log2FoldChange < -0.5]

pdf("BTZ_venn.pdf")
RELA_KO_BTZ_euler <- eulerr::euler(combinations = list(
  RELA_KO_targets = RELA_KO_expression_chip,
  BTZ_down = unique(BTZ_expression_down),
  PBRM1_KO_up = PBRM1_KO_expression_up
  )
)
plot(RELA_KO_BTZ_euler, quantities = list(type = c("percent", "counts")))
dev.off()

gene_127=intersect(PBRM1_KO_expression_up,unique(BTZ_expression_down))[!intersect(PBRM1_KO_expression_up,unique(BTZ_expression_down)) %in% RELA_KO_expression_chip]




PBRM1_KO_BTZ_RELA <- intersect(PBRM1_RELA_KO_intersect, BTZ_expression_down)


### hypergeometric calculation
# overlap between RELA targets and upregulated targets after PBRM1 KO
q <- length(intersect(RELA_KO_expression_chip,PBRM1_KO_expression_up))-1
m <- length(which(PBRM1_KO_expression_up %in%  RELA_KO_expression$GeneSymbol))
n <- length(RELA_KO_expression$GeneSymbol) - m
k <- length(RELA_KO_expression_chip)

pval <- phyper(q,m,n,k,lower.tail = FALSE)
#8.405412e-101

# overlap between RELA targets and downregulated targets after BTZ treatment
q <- length(intersect(RELA_KO_expression_chip,BTZ_expression_down))-1
m <- length(which(BTZ_expression_down %in%  RELA_KO_expression$GeneSymbol))
n <- length(RELA_KO_expression$GeneSymbol) - m
k <- length(RELA_KO_expression_chip)

pval <- phyper(q,m,n,k,lower.tail = FALSE)
#3.74127e-55


