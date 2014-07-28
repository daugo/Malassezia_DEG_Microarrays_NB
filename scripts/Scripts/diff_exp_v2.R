rm(list=ls(all=TRUE))
#Libraries
library(oligo)
library("pd.110413.mglo.cag.exp")
library(genefilter)
library(limma)
#++++++++++++++

#Device
pdf("~/Documents/Yulien_data/Scripts/4DS_vs_4DST80.pdf")
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
pd <- new("AnnotatedDataFrame", data=theData, varMetadata=metadata)
#ExpressionFeatureSet building using xys_files
maqc <- read.xysfiles( xys.files, pkgname="pd.110413.mglo.cag.exp", phenoData=pd )


#exprs(maqc)[10001:10010, 1:2]
#Figures
boxplot(maqc, main="Raw Data")
hist(maqc, main="Raw Data")

#RMA algorithm
#Expresion set
eset <- rma(maqc)

#boxplot(eset, transfo=identity, main="After RMA")
#hist(eset, transfo=identity, main="After RMA")

#Assesing differential expression
e_matrix <- exprs(eset)


ComparingTreatments = function(treatment="4DS", control="4DST80", treatments_order=c("4DS", "4DST80", "4DST80", "4DS", "4DST80", "4DS") ) {
	

	index_1 <- which(eset[["Sample_Name"]] == treatment)
	index_2 <- which(eset[["Sample_Name"]] == control)
	all <- c(index_1,index_2)

	d_means <- rowMeans(e_matrix[, index_1])-rowMeans(e_matrix[, index_2])
	a_means <- rowMeans(e_matrix[,all])

	relevant_treatments <- treatments_order
	treatment_levels <- c(treatment, control)

	tt <- rowttests(e_matrix[,all], factor( relevant_treatments, levels= treatment_levels ))
	lod <- -log10(tt[["p.value"]])
	smoothScatter(a_means, d_means, xlab="Average Intensity", ylab="Log-ratio", main=paste(treatment,"vs",control))
	abline(h=c(-1, 1), col=2)

	design <- model.matrix(~factor ( relevant_treatments, levels= treatment_levels ))
	fit <- lmFit(eset[,all], design)
	ebayes <- eBayes(fit)
	lod <- -log10(ebayes[["p.value"]][,2])
	mtstat<- ebayes[["t"]][,2]

	o1 <- order(abs(d_means), decreasing=TRUE)[1:25]
	o2 <- order(abs(mtstat), decreasing=TRUE)[1:25]
	o <- union(o1, o2)

	#Plot
	smoothScatter(d_means, lod, main= paste("Moderated t",treatment,"vs",control), xlab="Log-ratio", ylab="LOD")
	points(d_means[o1], lod[o1], pch=18,col="blue")
	points(d_means[o2], lod[o2], pch=1,col="red")
}


ComparingTreatments("4DS", "4DST80", c("4DS", "4DST80", "4DST80", "4DS", "4DST80", "4DS"))
ComparingTreatments("4DS", "M.furfur", c("M.furfur", "M.furfur", "4DS", "M.furfur", "4DS", "4DS"))
ComparingTreatments("4DST80", "M.furfur", c("M.furfur", "M.furfur", "4DST80", "M.furfur", "4DST80", "4DST80"))
ComparingTreatments("M.furfur", "M.globosa", c("M.furfur", "M.globosa", "M.globosa", "M.globosa", "M.furfur","M.furfur"))


dev.off()