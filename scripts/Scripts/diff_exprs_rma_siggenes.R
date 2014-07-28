rm(list=ls(all=TRUE))
#Libraries
library(oligo)
library("pd.110413.mglo.cag.exp")
library(genefilter)
library(siggenes)

#++++++++++++++

#Device
pdf("~/Documents/Yulien_data/Scripts/Results_Differential_Expression_siggenes.pdf")
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

treatment <- "4DST80"
control <- "4DS"

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


	exprs.cl <- c(rep(0,3),rep(1,3))

	sam.out <- sam(exprs_set_rma, exprs.cl, method= d.stat,  rand= 123)
	summary(sam.out)
	plot(sam.out, 4.2)
	sum.sam.out <- summary(sam.out, 4.2)
	sum.sam.out
}

ComparingTreatments(treatment="4DST80", control="4DS")

dev.off()

