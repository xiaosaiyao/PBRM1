library(DiffBind)
library(GenomicRanges)
library(rtracklayer)

cwd="/scratch/users/astar/gis/yaoxs/PBRM1/chip/diffbind"
setwd(cwd)

#samples <- read.csv("project214.RELA.csv",stringsAsFactors=F)
#chipseq <- dba(sampleSheet=samples)
#chipseq <- dba.count(chipseq, filter=2)

# normalize background regions using csaw gives best differential
# 
#chipseq <- dba.normalize(chipseq, library=DBA_LIBSIZE_BACKGROUND, method=DBA_ALL_METHODS,normalize=DBA_NORM_NATIVE,background=TRUE)
#saveRDS(chipseq, "project214.RELA.chipseq.rds")
chipseq <- readRDS("project214.RELA.chipseq.rds")
#### WT ####

contrast.WT <- dba.contrast(chipseq,  design="~Replicate + Condition", contrast=c("Condition","WT_BRG1i","WT_DMSO"))
contrast.WT <- dba.analyze(contrast.WT, bBlacklist=F, bGreylist =F)

saveRDS(contrast.WT, "project214.RELA.WT.BRG1i.contrast.rds")
contrast.WT=readRDS("project214.RELA.WT.BRG1i.contrast.rds")

#plot MA plot
pdf("project214.RELA.WT.BRG1i.MAplot.pdf")
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

export.bed(report.gain, "project214.diffbind.786O.RELA.WT.BRG1i.gain.bed")
export.bed(report.lost, "project214.diffbind.786O.RELA.WT.BRG1i.lost.bed")
export.bed(report.all, "project214.diffbind.786O.RELA.WT.BRG1i.all.bed")



#### KO ####
contrast.KO <- dba.contrast(chipseq,  design="~Replicate + Condition", contrast=c("Condition","KO_BRG1i","KO_DMSO"))
contrast.KO <- dba.analyze(contrast.KO, bBlacklist=F, bGreylist =F)

saveRDS(contrast.KO, "project214.RELA.KO.BRG1i.contrast.rds")
contrast.KO=readRDS("project214.RELA.KO.BRG1i.contrast.rds")

#plot MA plot
pdf("project214.RELA.KO.BRG1i.MAplot.pdf")
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

export.bed(report.gain, "project214.diffbind.786O.RELA.KO.BRG1i.gain.bed")
export.bed(report.lost, "project214.diffbind.786O.RELA.KO.BRG1i.lost.bed")
export.bed(report.all, "project214.diffbind.786O.RELA.KO.BRG1i.all.bed")

