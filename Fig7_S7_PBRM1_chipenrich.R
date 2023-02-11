library(chipenrich)
library(eulerr)
library(ggplot2)

dir = "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/ChIPseq/diffbind/"
setwd(dir)

#pathway analysis on BRG1
peakdot = list()

BRG1_gain_bed = "diffbind.786O.BRG1.gain.bed"
BRG1_gain_results = chipenrich(
  peaks = BRG1_gain_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = BRG1_gain_results$results
print(results[1:5, 1:5])

results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[1]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue")) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() + ggtitle("SMARCA4 gained")


BRG1_lost_bed = "diffbind.786O.BRG1.lost.bed"
BRG1_lost_results = chipenrich(
  peaks = BRG1_lost_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = BRG1_lost_results$results
print(results[1:5, 1:5])

results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[2]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue")) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() +
  ggtitle("SMARCA4 lost")

#pathway analysis on ARID2

ARID2_gain_bed = "diffbind.786O.ARID2.gain.bed"
ARID2_gain_results = chipenrich(
  peaks = ARID2_gain_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = ARID2_gain_results$results
print(results[1:5, 1:5])
results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[3]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue"), drop = F) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() +
  ggtitle("ARID2 gained")


ARID2_lost_bed = "diffbind.786O.ARID2.lost.bed"
ARID2_lost_results = chipenrich(
  peaks = ARID2_lost_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = ARID2_lost_results$results
print(results[1:8, ])
results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[4]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue"), drop = F) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() +
  ggtitle("ARID2 lost")




#pathway analysis on BRD7

#early

BRD7_gain_bed = "diffbind.786O.BRD7.early.gain.bed"
BRD7_gain_results = chipenrich(
  peaks = BRD7_gain_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = BRD7_gain_results$results
print(results[1:5, 1:5])
results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[5]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue"), drop = F) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() +
  ggtitle("BRD7 early gained")


BRD7_lost_bed = "diffbind.786O.BRD7.early.lost.bed"
BRD7_lost_results = chipenrich(
  peaks = BRD7_lost_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = BRD7_lost_results$results
print(results[1:8, ])
results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[6]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue"), drop = F) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() +
  ggtitle("BRD7 early lost")


#late

BRD7_gain_bed = "diffbind.786O.BRD7.late.gain.bed"
BRD7_gain_results = chipenrich(
  peaks = BRD7_gain_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = BRD7_gain_results$results
print(results[1:5, 1:5])
results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[7]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue"), drop = F) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() +
  ggtitle("BRD7 late gained")


BRD7_lost_bed = "diffbind.786O.BRD7.late.lost.bed"
BRD7_lost_results = chipenrich(
  peaks = BRD7_lost_bed,
  genome = 'hg19',
  genesets = c("hallmark"),
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)
results = BRD7_lost_results$results
print(results[1:8, ])
results$Description = factor(as.character(results$Description), levels = rev(unique(results$Description)))
peakdot[[8]] = ggplot(head(results, 8), aes(y = Effect, x = Description, color =
                                              Status)) + scale_color_manual(values = c("red", "blue"), drop = F) +
  geom_point(stat = "identity", aes(size = -log10(P.value))) +
  coord_flip() +
  theme_bw() +
  ggtitle("BRD7 latelost")

pdf("peakdot.pdf", width=13, height=12)
gridExtra::grid.arrange(grobs=peakdot, nrow=4, ncol=2)
dev.off()

