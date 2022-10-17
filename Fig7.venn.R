library(ChIPpeakAnno)
library(rtracklayer)
library(GenomicRanges)


dir <- "/scratch/users/astar/gis/yaoxs/PBRM1/chip/diffbind/"
setwd(dir)

# Load peaks

KO_BRG1i_RELA_loss_bed <- "project214.diffbind.786O.RELA.KO.BRG1i.lost.DBA_DESEQ2.bed"
KO_BRG1i_RELA_loss <-  import(KO_BRG1i_RELA_loss_bed) 
KO_BRG1i_BRG1_loss_bed <- "project214.diffbind.786O.BRG1.KO.BRG1i.lost.DBA_DESEQ2.bed"
KO_BRG1i_BRG1_loss <-  import(KO_BRG1i_BRG1_loss_bed)

WT_BRG1i_RELA_loss_bed <- "project214.diffbind.786O.RELA.WT.BRG1i.lost.DBA_DESEQ2.bed"
WT_BRG1i_RELA_loss <-  import(WT_BRG1i_RELA_loss_bed) 
WT_BRG1i_BRG1_loss_bed <- "project214.diffbind.786O.BRG1.WT.BRG1i.lost.DBA_DESEQ2.bed"
WT_BRG1i_BRG1_loss <-  import(WT_BRG1i_BRG1_loss_bed)


# Plot venn diagrams
pdf("venn.BRG1.loss.pdf")
makeVennDiagram(list(KO_BRG1i_RELA_loss, KO_BRG1i_BRG1_loss), 
		NameOfPeaks=c("KO_BRG1i_RELA_loss", "KO_BRG1i_BRG1_loss")) 

makeVennDiagram(list(WT_BRG1i_RELA_loss, WT_BRG1i_BRG1_loss), 
		NameOfPeaks=c("WT_BRG1i_RELA_loss", "WT_BRG1i_BRG1_loss")) 

makeVennDiagram(list(KO_BRG1i_RELA_loss, WT_BRG1i_RELA_loss), 
		NameOfPeaks=c("KO_BRG1i_RELA_loss", "WT_BRG1i_RELA_loss"))
dev.off()


