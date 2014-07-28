rm(list=ls(all=TRUE))
#Libraries
library(oligo)
library("pd.110413.mglo.cag.exp")
library(genefilter)
library(limma)

#++++++++++++++

#Device
pdf("~/Documents/Yulien_data/Scripts/Results_Differential_Expression.pdf")
#+++++++++++

#Input Data
#xys_files
xys_files_path <- "~/Documents/Yulien_data/Raw_Data_Files/Pair_Files"
xys.files <- list.xysfiles( xys_files_path, full.names=TRUE ) #list xys_files
#phenoData
pheno_file <- "~/Documents/Yulien_data/Design_information/chip_data.tab"
theData = read.table(pheno_file, row.names=1, header=TRUE,sep="\t")
theData = subset(theData, select= c("Sample_Name"))

lvls <- c( "channel1", "channel2", "_ALL_" )
metadata <- data.frame( channel=factor("_ALL_", levels=lvls), labelDescription="Channels" )

#AnnotatedDataFrame building
pheno_data <- new("AnnotatedDataFrame", data=theData, varMetadata=metadata)
#ExpressionFeatureSet building using xys_files
raw_exprs_set <- read.xysfiles( xys.files, pkgname="pd.110413.mglo.cag.exp", phenoData=pheno_data )


ComparingTreatments = function(treatment="4DS", control="4DST80") {

	index_1 <- which(raw_exprs_set[["Sample_Name"]] == treatment)
	index_2 <- which(raw_exprs_set[["Sample_Name"]] == control)
	all_indexes <- c(index_1,index_2)


	exprs_set_subset <- raw_exprs_set[,all_indexes]
	boxplot(exprs_set_subset, main=paste("Raw Data", treatment,"vs",control))
	hist(exprs_set_subset, main=paste("Raw Data", treatment, "vs", control))
	
	 
	exprs_set_rma <-  rma(exprs_set_subset)
	
	boxplot(exprs_set_rma, main=paste("After RMA", treatment,"vs",control))
	hist(exprs_set_rma, main=paste("After RMA", treatment, "vs", control))
	

	e_matrix <- exprs(exprs_set_rma)

	d_means <- rowMeans(e_matrix[, 4:6])-rowMeans(e_matrix[, 1:3])
	a_means <- rowMeans(e_matrix)
	sum(abs(d_means)>1)

	tt <- rowttests(e_matrix, factor(exprs_set_rma[["Sample_Name"]]))
	lod <- -log10(tt[["p.value"]])
	smoothScatter(a_means, d_means, xlab="Average Intensity", ylab="Log-ratio", main=paste(treatment,"vs",control))
	abline(h=c(-2, 2), col=2)
	o1 <- order(abs(d_means), decreasing=TRUE)[1:25]
	o2 <- order(abs(tt), decreasing=TRUE)[1:25]
	o <- union(o1, o2)
	

	smoothScatter(d_means, lod, main="t-test", xlab="Log-Ratio", ylab="LOD")
	points(d_means[o1], lod[o1], pch=18,col="blue")
	points(d_means[o2], lod[o2], pch=1,col="red")
	abline(h=2, v=c(-1, 1))

	design <- model.matrix(~factor(exprs_set_rma[["Sample_Name"]]))
	fit <- lmFit(exprs_set_rma, design)
	ebayes <- eBayes(fit)
	lod <- -log10(ebayes[["p.value"]][,2])
	mtstat <- ebayes[["t"]][,2]

	o1 <- order(abs(d_means), decreasing=TRUE)[1:25]
	o2 <- order(abs(mtstat), decreasing=TRUE)[1:25]
	o <- union(o1, o2)
	smoothScatter(d_means, lod, main="Moderated t", xlab="Log-ratio", ylab="LOD")
	points(d_means[o1], lod[o1], pch=18,col="blue")
	points(d_means[o2], lod[o2], pch=1,col="red")
	abline(h=2, v=c(-1, 1))

	tab <- topTable(ebayes, coef=2, adjust="fdr", n=10)
	tab

}

ComparingTreatments(treatment="4DST80", control="4DS")
ComparingTreatments(treatment="4DST80", control="M.furfur")
ComparingTreatments(treatment="4DST80", control="M.globosa")
ComparingTreatments(treatment="4DS", control="M.furfur")
ComparingTreatments(treatment="4DS", control="M.globosa")
ComparingTreatments(treatment="M.furfur", control="M.globosa")


dev.off()

