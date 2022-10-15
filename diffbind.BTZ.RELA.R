library(DiffBind)
library(GenomicRanges)
library(rtracklayer)

cwd="/scratch/users/astar/gis/yaoxs/PBRM1/chip/diffbind"
setwd(cwd)

samples <- read.csv("project213.RELA.csv",stringsAsFactors=F)
chipseq <- dba(sampleSheet=samples)
chipseq <- dba.count(chipseq, filter=2)

# normalize background regions using csaw gives best differential
# 
chipseq <- dba.normalize(chipseq, library=DBA_LIBSIZE_BACKGROUND, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,background=TRUE)
saveRDS(chipseq, "project213.RELA.chipseq.rds")

#### WT ####

contrast.WT <- dba.contrast(chipseq,  design="~Replicate + Condition", contrast=c("Condition","WT_BTZ","WT_DMSO"))
contrast.WT <- dba.analyze(contrast.WT, bBlacklist=F, bGreylist =F)

saveRDS(contrast.WT, "project213.RELA.WT.BTZ.contrast.rds")
contrast.WT=readRDS("project213.RELA.WT.BTZ.contrast.rds")

#plot MA plot
pdf("project213.RELA.WT.BTZ.MAplot.pdf")
dba.plotMA(contrast.WT, th=0.1, fold=log2(1.5), yrange=c(-10,10))
dev.off()

#generate report
report=dba.report(contrast.WT, contrast=1, fold=log2(1.5), th=0.1)
report.gain=report[report$Fold>0]
report.lost=report[report$Fold<0]

report.all=dba.report(contrast.WT, contrast=1, th=1)
report.lost$name=paste0("peak_",1:length(report.lost))
report.gain$name=paste0("peak_",1:length(report.gain))
report.all$name=paste0("peak_",1:length(report.all))

export.bed(report.gain, "project213.diffbind.786O.RELA.WT.BTZ.gain.bed")
export.bed(report.lost, "project213.diffbind.786O.RELA.WT.BTZ.lost.bed")
export.bed(report.all, "project213.diffbind.786O.RELA.WT.BTZ.all.bed")



#### KO ####
contrast.KO <- dba.contrast(chipseq,  design="~Replicate + Condition", contrast=c("Condition","KO_BTZ","KO_DMSO"))
contrast.KO <- dba.analyze(contrast.KO, bBlacklist=F, bGreylist =F)

saveRDS(contrast.KO, "project213.RELA.KO.BTZ.contrast.rds")
contrast.KO=readRDS("project213.RELA.KO.BTZ.contrast.rds")

#plot MA plot
pdf("project213.RELA.KO.BTZ.MAplot.pdf")
dba.plotMA(contrast.KO, th=0.1, fold=log2(1.5), yrange=c(-10,10))
dev.off()

#generate report
report=dba.report(contrast.KO, contrast=1, fold=log2(1.5), th=0.1)
report.gain=report[report$Fold>0]
report.lost=report[report$Fold<0]

report.all=dba.report(contrast.KO, contrast=1, th=1)
report.lost$name=paste0("peak_",1:length(report.lost))
report.gain$name=paste0("peak_",1:length(report.gain))
report.all$name=paste0("peak_",1:length(report.all))

export.bed(report.gain, "project213.diffbind.786O.RELA.KO.BTZ.gain.bed")
export.bed(report.lost, "project213.diffbind.786O.RELA.KO.BTZ.lost.bed")
export.bed(report.all, "project213.diffbind.786O.RELA.KO.BTZ.all.bed")

