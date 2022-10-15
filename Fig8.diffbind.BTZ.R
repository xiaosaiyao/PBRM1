library(DiffBind)
library(GenomicRanges)
library(rtracklayer)

cwd <- "/scratch/users/astar/gis/yaoxs/PBRM1/chip/diffbind"
setwd(cwd)

factors <- c("BRG1","RELA")

for (factor in factors) {

    # Load samples
    #samples <- read.csv("project213.BRG1.csv",stringsAsFactors=F)
    #chipseq <- dba(sampleSheet=samples)
    #chipseq <- dba.count(chipseq, filter=2)

    # Normalize background regions using csaw gives best differential
    #chipseq <- dba.normalize(chipseq, 
#			     library=DBA_LIBSIZE_BACKGROUND, 
#			     method=DBA_ALL_METHODS,
#			     normalize=DBA_NORM_NATIVE,
#			     background=TRUE)

    filename <- paste0("project213.", factor, ".chipseq.rds")
#    saveRDS(chipseq, filename)
    chipseq <- readRDS(filename)
    
    for (PBRM1_status in c("WT","KO")){
        contrast <- dba.contrast(chipseq,
				 design="~Replicate + Condition", 
				 contrast=c("Condition",
					    paste0(PBRM1_status,"_BTZ"),
					    paste0(PBRM1_status,"_DMSO")))

	contrast <- dba.analyze(contrast, 
				bBlacklist=F, 
				bGreylist=F, 
				method=DBA_ALL_METHODS)
		
	filename <- paste0("project213.", factor,".", PBRM1_status, ".BTZ.contrast.rds")
	saveRDS(contrast, filename)
	contrast <- readRDS(filename)
	
	for (diff_method in c("DBA_EDGER","DBA_DESEQ2")){
		# Plot MA plot
		pdf(paste0("project213.", factor, ".", PBRM1_status, ".BTZ.MAplot.", 
			   diff_method,".pdf"))

		if (diff_method == "DBA_EDGER"){
			dba.plotMA(contrast, th=0.1, yrange=c(-10,10), 
				   fold=log2(1.5), method=DBA_EDGER)
		
		}else if (diff_method == "DBA_DESEQ2"){
			dba.plotMA(contrast, th=0.1, yrange=c(-10,10),
                                   fold=log2(1.5), method=DBA_DESEQ2)

		}
		dev.off()
		
		# Generate report
		
		if (diff_method == "DBA_EDGER"){
			report <- dba.report(contrast, contrast=1, fold=log2(1.5), 
					     th=0.1, method=DBA_EDGER)
		} else if (diff_method =="DBA_DESEQ2"){
			report <- dba.report(contrast, contrast=1, fold=log2(1.5),
					     th=0.1, method=DBA_DESEQ2)
		}
		
		if (!is.null(report)){
			report.gain <- report[report$Fold>0]
			report.lost <- report[report$Fold<0]
		}
		

		report.all <- dba.report(contrast, contrast=1, th=1)
		report.all$name <- paste0("peak_",1:length(report.all))
		
		if (exists("report.gain")){
			if (!is.null(report.gain)){
				message("exporting gained peaks")
				report.lost$name <- 
					paste0("peak_",1:length(report.lost))
				export.bed(report.gain,
					   paste0("project213.diffbind.786O.", factor, ".", 
						  PBRM1_status, ".BTZ.gain.", 
						  diff_method,".bed")
					   )
			}
		}
		
		if (exists("report.lost")) {
			if (!is.null(report.lost)){
				message("exporting lost peaks")
				report.lost$name <- 
					paste0("peak_",1:length(report.lost))
				export.bed(report.lost,
					   paste0("project213.diffbind.786O.", factor, ".", 
						  PBRM1_status, ".BTZ.lost.", 
						  diff_method, ".bed")
					   )
			}
		}
		
		if (exists("report.all")) {
			message("exporting all peaks")
			head(report.all)
			export.bed(report.all,
				   paste0("project213.diffbind.786O.", factor, ".", 
					  PBRM1_status,".BTZ.all.", diff_method,".bed")
				   )
		}
	}
    }

}


