# 786O WT/KO BTZ Fig8

library(DESeq2)
library(RColorBrewer)
library(pheatmap)

dir <- "/Users/abc/Desktop/ASTAR/Projects/PBRM1/Analysis/RNAseq/RSEM/"
setwd(dir)

### DESEQ2 to determine differential gene expression

datmatrix <- 'O786_DMSOBTZCFZ9hr_WTvsKO.genes.results.header.txt'
sampletable <- 'sample_table_786O_DMSOCTZBTZ9hr.csv'
outfile <- "Deseq2_786O_PB_KO_noIFNg.txt"


colData <- read.delim2(sampletable, row.names=1, header = TRUE, 
                       sep = ",", quote = "\"", dec = ",", fill = TRUE, 
                       comment.char = "#", check.names=FALSE)
countData <- read.delim2(datmatrix, row.names=1, header = TRUE, 
                         sep = "\t", quote = "\"", dec = ",", fill = TRUE, 
                         comment.char = "#", check.names=FALSE)

countData <- as.matrix(countData)
storage.mode(countData) <- "double"
countData <- round(countData)

### select samples 
samplesToUse <- rownames(colData) #colnames(countData)
designToUse <- c("design")
countDataSel <- countData[,samplesToUse]
colDataSel <- data.frame(as.matrix(colData)[samplesToUse, designToUse])
colnames(colDataSel) <- designToUse


### Differential analysis
dds <- DESeqDataSetFromMatrix(countData = countDataSel, colData = colDataSel, design = ~ design)
dds <- dds[ rowMeans(counts(dds)) > 1, ] #filter those not expressed. use your own cutoff..

dds <- DESeq(dds)
res <- results(dds, contrast=c("design","KO","WT"))


rld <- rlog(dds, blind=TRUE)
rldMat<-assay(rld)

### Annotate to gene symbol 
E2G <- read.delim2("GRCh37.p13.E2G.txt", row.names=1, header = FALSE, sep = "\t", 
                   quote = "\"", dec = ",", fill = TRUE, comment.char = "#", 
                   check.names=FALSE, stringsAsFactors = F)

### Output Exp Level
resOrdered <- res[order(res$pvalue),]
resOrderedDF <- as.data.frame(resOrdered)
expressedIds <- rownames(resOrderedDF)
expLog2 <- assay(rld)[expressedIds,] #rlog normalized...
exp <- countData[expressedIds,rownames(colData)] #countDataSel[expressedIds,]
gn <- as.matrix(E2G[expressedIds,])
ID <- rownames(exp)
colnames(gn) <- c("GeneSymbol")

### plot heatmap 
TNFgenes <- c("CCL2","TNF","LAMC2","LAMB3","BCL2A1","PDLIM4","APOL3","TNFAIP3","GCH1","TNFRSF9","JAG1",
              "CXCL2", "SOX9","TMEM132A","ADAMTS9","DNMBP","ABLIM3","IL11","NCEH1")


TNFsymbols <- rownames(E2G)[match(TNFgenes,E2G[,1])]
TNFexpLog2 <- expLog2[match(TNFsymbols, rownames(expLog2)),1:12]
rownames(TNFexpLog2) <- TNFgenes

pdf("heatmap.BTZ.pdf")
pheatmap( TNFexpLog2[,c(1:3,7:9,4:6,10:12)],
         scale = "row",
         show_colnames = TRUE,
         show_rownames = TRUE, 
         cluster_cols = FALSE, 
         breaks = seq(-2, 2, by = 4/31), 
         color = viridis::plasma(30))

dev.off()



