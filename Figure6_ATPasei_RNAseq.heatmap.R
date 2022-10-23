#Figure 6C
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

dir <- "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/RNAseq/RSEM/"
setwd(dir)

datmatrix <- "BRG1i.header.txt"
sampletable <- "samplet_table_O786_BRG1i.csv"
outfile <- "Deseq2_786O_WT_KO_BRG1i.txt"

#-----load data-------------------

colData <- read.delim2(sampletable, row.names=1, header = TRUE, sep = ",", quote = "\"", dec = ",", fill = TRUE, comment.char = "#", check.names=FALSE)
countData <- read.delim2(datmatrix, row.names=1, header = TRUE, sep = "\t", quote = "\"", dec = ",", fill = TRUE, comment.char = "#", check.names=FALSE)

countData <- as.matrix(countData)
storage.mode(countData) <- "double"
countData <- round(countData)

#-----select samples-------------------
samplesToUse <- rownames(colData) #colnames(countData)
designToUse <- c("design")
countDataSel <- countData[,samplesToUse]
colDataSel <- data.frame(as.matrix(colData)[samplesToUse, designToUse])
colnames(colDataSel) <- designToUse

#-----DESEQ2---------------------------
dds <- DESeqDataSetFromMatrix(countData = countDataSel, colData = colDataSel, design = ~ design)
dds <- dds[ rowMeans(counts(dds)) > 1, ] #filter those not expressed. use your own cutoff..

dds <- DESeq(dds)
res <- results(dds, contrast=c("design","KO_DMSO","KO_BRG1i"))
nrow(dds)

rld <- rlog(dds, blind=TRUE)
rldMat<-assay(rld)

#-----Annotation-------------------
E2G<-read.delim2("GRCh37.p13.E2G.txt", row.names=1, header = FALSE, sep = "\t", quote = "\"", dec = ",", fill = TRUE, comment.char = "#", check.names=FALSE)

#-----Output Exp Levels-------------------
resOrdered <- res[order(res$pvalue),]
resOrderedDF <- as.data.frame(resOrdered)
expressedIds<-rownames(resOrderedDF)
expLog2<-assay(rld)[expressedIds,] #rlog normalized...
exp<-countData[expressedIds,rownames(colData)] #countDataSel[expressedIds,]
gn<-as.matrix(E2G[expressedIds,])
ID<-rownames(exp)
colnames(gn)<-c("GeneSymbol")


write.table(cbind(gn, resOrderedDF,expLog2, exp), file=outfile, sep="\t", col.names = NA)

#-----Plot Heatmap --------------------

BRG1i <- read.delim(file = "Deseq2_786O_WT_KO_BRG1i.txt")
BRG1i_matrix <- BRG1i[,c( "WT_DMSO_1",	"WT_DMSO_2",	"WT_BRG1i_1",	"WT_BRG1i_2",
                       "KO_DMSO_1",	"KO_DMSO_2",	"KO_BRG1i_1",	"KO_BRG1i_2")]

index <- match(direct_NFKB, BRG1i$GeneSymbol)
TNF_matrix <- BRG1i_matrix[index,]
direct_NFKB <- unlist(t(read.delim("direct_tnf_targets.txt")))
rownames(TNF_matrix) <- direct_NFKB

# genes downregulated more than 2 fold
BRG1i_downgenes <- BRG1i$GeneSymbol[which(BRG1i$log2FoldChange> log2(2))]
index_down_TNF <- intersect(index, match(BRG1i_downgenes, BRG1i$GeneSymbol))
BRG1i_downgenes_matrix <- BRG1i_matrix[index_down_TNF,]
rownames(BRG1i_downgenes_matrix) <- BRG1i$GeneSymbol[index_down_TNF]

pdf("BRG1i.pdf", height = 5)
pheatmap(TNF_matrix,
         scale = "row",
         show_colnames = TRUE,
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE,  
         breaks = seq(-2, 2, by = 4/31), 
         color = viridis::plasma(30))

pheatmap(BRG1i_downgenes_matrix,
         scale = "row",
         show_colnames = TRUE,
         show_rownames = FALSE, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE,  
         breaks = seq(-2, 2, by = 4/31), 
         color = viridis::plasma(30))


dev.off()

