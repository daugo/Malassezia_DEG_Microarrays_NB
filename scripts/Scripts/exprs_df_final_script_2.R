rm(list=ls(all=TRUE))
#Libraries
library(oligo)
library("pd.110413.mglo.cag.exp")

#Functions

#build data.frame with signals information to be used for ggplot graphics
ggplot_dataframe = function(x) {

	exprs_df <- data.frame(signal= numeric(), Dataset=factor(),  Treatment= factor())
	if (class(x) == "ExpressionFeatureSet") {
		for (i in 1:length(pm(x)[1,])) {
			
			Dataset <- colnames(pm(x))[i]
			signals <- pm(x)[,i]
			Treatment <- as.character(x[["Sample_Name"]][i])
			exprs_df <- rbind(
							exprs_df,
							data.frame(
								signal=signals,
								Dataset=rep(Dataset,length(signals)),
								Treatment=rep(Treatment,length(signals)))
							)
		}
	}
	else {
		for (i in 1:length(exprs(x)[1,])) {
		
			Dataset <- colnames(exprs(x))[i]
			signals <- exprs(x)[,i]
			Treatment <- as.character(x[["Sample_Name"]][i])
			exprs_df <- rbind(
							exprs_df,
							data.frame(
								signal=signals,
								Dataset=rep(Dataset,length(signals)),
								Treatment=rep(Treatment,length(signals)))
							)
		}
	}
	return(exprs_df)
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
out_path = "~/Dropbox/Malassezia_2013_01/Yulien_data/Scripts/"
out_base_name = "A_VSN_Final_Results_Differential_Expression"
#Treatment and control to analyse 
treatment <- "M.furfur"
control <- "M.globosa"
out_name = paste(out_base_name,treatment,"vs",control,sep="")
#Device
pdf(paste(out_path, out_name, ".pdf", sep=""))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Input Data
#xys_files
xys_files_path <- "~/Dropbox/Malassezia_2013_01/Yulien_data/Raw_Data_Files/Pair_Files"
xys.files <- list.xysfiles( xys_files_path, full.names=TRUE ) #list xys_files
#phenoData
pheno_file <- "~/Dropbox/Malassezia_2013_01/Yulien_data/Design_information/chip_data.tab"
theData = read.table(pheno_file, row.names=1, header=TRUE,sep="\t")
theData = subset(theData, select= c("Sample_Name"))

lvls <- c( "channel1", "channel2", "_ALL_" )
metadata <- data.frame( channel=factor("_ALL_", levels=lvls), labelDescription="Channels" )

#AnnotatedDataFrame building
pheno_data <- new("AnnotatedDataFrame", data=theData, varMetadata=metadata)
#ExpressionFeatureSet building using xys_files
raw_exprs_set <- read.xysfiles( xys.files,
								pkgname="pd.110413.mglo.cag.exp",
								phenoData=pheno_data,
								sampleNames=c("M.furfur-1",
											"M.globosa-1",
											"M.globosa-2",
											"M.globosa-3",
											"M.furfur-2",
											"4DS-1",
											"4DST80-1",
											"M.furfur-3",
											"4DST80-2",
											"4DS-2",
											"4DST80-3",
											"4DS-3")
								)


#Idexes of the datasets related to the treatment and control specified
index_1 <- which(raw_exprs_set[["Sample_Name"]] == treatment)
index_2 <- which(raw_exprs_set[["Sample_Name"]] == control)
all_indexes <- c(index_1,index_2)
#Create a ExpressioFeatureSet subset with the datasets related to the specified treatment and control
exprs_set_subset <- raw_exprs_set[,all_indexes]
#Oligo default distribution graphics for raw data
#boxplot(exprs_set_subset, main=paste("Raw Data", treatment,"vs",control))
#hist(exprs_set_subset, main=paste("Raw Data", treatment, "vs", control))

#ggplot distribution graphics for raw data
exprs_df <- ggplot_dataframe(exprs_set_subset)
library(ggplot2)
p <- ggplot(exprs_df,aes(y=log2(signal),x=Dataset))
p + geom_boxplot(aes(fill=Treatment)) #+ geom_jitter(alpha = I(1/300))

p <- ggplot(exprs_df, aes(x=log2(signal)))
p + geom_density(aes(fill=Dataset), alpha= 0.3)

#ExpressionFeatureSet to matrix

eset <- rma(exprs_set_subset, normalize=FALSE, background=FALSE)

e_matrix <- 2^exprs(eset)

#Basic MAplot applying log2 to raw data
d_means <- log2(rowMeans(e_matrix[, 4:6])) - log2(rowMeans(e_matrix[, 1:3]))
a_means <- (log2(rowMeans(e_matrix[,4:6])) + log2(rowMeans(e_matrix[,1:3]))) / 2
smoothScatter(a_means, d_means, xlab="Average Intensity", ylab="Log2-ratio", main=paste(treatment,"vs",control, "(Applying log2 to raw data)"))

#Signal Data Processing 
#bacgroundCorrect (oligo method)
raw <- exprs_set_subset
#raw <- backgroundCorrect(exprs_set_subset, "rma") # Add noise to the dataset, in this case
raw_begin <- raw
pms = pm(raw) #ExpressionFeatureSet to matrix, ignoring spots with NA signal
pmsVSN = vsn::vsnMatrix(pms) #vsn "normalization" of signal data
require("vsn")
#meanvsSd plots 
meanSdPlot(pmsVSN, ranks=TRUE)
meanSdPlot(pmsVSN, ranks=FALSE)
#vsn normalized data, store in a vsn object to ExpressionFeatureSet
pm(raw) <- exprs(pmsVSN)

#ggplot distribution graphics for normalized data
exprs_df <- ggplot_dataframe(raw)

p <- ggplot(exprs_df, aes(y=signal, x=Dataset))
p + geom_boxplot(aes(fill=Treatment)) #+ geom_jitter(alpha = I(1/300))

p <- ggplot(exprs_df, aes(x=signal))
p + geom_density(aes(fill=Dataset), alpha= 0.3)

#oligo defaults graphics
#boxplot(eset, main=paste("Raw Data", treatment,"vs",control), transfo=identity)
#hist(eset, main=paste("Raw Data", treatment,"vs",control), transfo=identity)

#summarization of normalize data using rma oligo method = median.polish
eset <- rma(raw, normalize = FALSE, background = FALSE)

#ExpressionSet to matrix
e_matrix <- 2^exprs(eset)

eset_treatments_info = as.character(theData[all_indexes,1])
eset_theData = data.frame(Sample_Name= factor(eset_treatments_info, levels=unique(eset_treatments_info)))
rownames(eset_theData) = rownames(theData)[all_indexes]
pd <- new("AnnotatedDataFrame", data=eset_theData, varMetadata=metadata)
eset <- new("ExpressionSet", exprs = e_matrix, phenoData = pd, annotation = "pd.110413.mglo.cag.exp")

exprs_df <- ggplot_dataframe(eset)

p <- ggplot(exprs_df, aes(y=signal, x=Dataset))
p + geom_boxplot(aes(fill=Treatment)) #+ geom_jitter(alpha = I(1/300))

p <- ggplot(exprs_df, aes(x=signal))
p + geom_density(aes(fill=Dataset), alpha= 0.3)

#MAplot of normalize data 
d_means <- rowMeans(e_matrix[, 4:6]) - rowMeans(e_matrix[, 1:3])
a_means <- rowMeans(e_matrix)
smoothScatter(a_means, d_means, xlab="Average Intensity", ylab="Ratio", main=paste(treatment,"vs",control, "(After VSN normalization)"))

#SAM analysis
library(siggenes)
exprs.cl <- c(rep(0,3),rep(1,3))

sam.out <- sam(eset, exprs.cl, method= d.stat,  rand= 123)
summary(sam.out)
plot(sam.out)
thres <- 6
plot(sam.out, thres)
sum.sam.out <- summary(sam.out, thres)
sam_csv_file <- paste(out_path, out_name,"sam.csv",sep="")
sam2excel(sam.out,thres,sam_csv_file)
sam_results <- read.csv(sam_csv_file, skip=19, header=TRUE)
num_sig_genes  <- dim(sam_results)[1]

#eBayes analysis
library(limma)

d_means <- rowMeans(eset_matrix[, 4:6])-rowMeans(eset_matrix[, 1:3])

design <- model.matrix(~factor(eset[["Sample_Name"]]))
fit <- lmFit(eset, design)
ebayes <- eBayes(fit)
lod <- -log10(ebayes[["p.value"]][,2])
mtstat <- ebayes[["t"]][,2]

o1 <- order(abs(d_means), decreasing=TRUE)[1:num_sig_genes]
o2 <- order(abs(mtstat), decreasing=TRUE)[1:num_sig_genes]
o <- union(o1, o2)

smoothScatter(d_means, lod, main="Moderated t", xlab="Ratio", ylab="LOD")
points(d_means[o1], lod[o1], pch=18,col="blue")
points(d_means[o2], lod[o2], pch=8,col="red")
#abline(h=2, v=c(-1, 1))

tab <- topTable(ebayes, coef=2, adjust="fdr", n=num_sig_genes)
ebayes_csv_file <- paste(out_path, out_name,"ebayes.csv",sep="")
write.csv(tab, file=ebayes_csv_file, row.names=FALSE)

imp_genes <-  head((tab$ID),num_sig_genes)
imp_genes_exprs <- matrix(nrow=num_sig_genes,ncol=6)
rownames(imp_genes_exprs) <- imp_genes
colnames(imp_genes_exprs) <- colnames(e_matrix)
for (i in imp_genes) {
	print(i)
	print(e_matrix[i,])
	imp_genes_exprs[i,] <-  e_matrix[i,]
}

e_matrix <- exprs(eset)
#MAplot of normalize data 
d_means <- rowMeans(e_matrix[, 4:6]) - rowMeans(e_matrix[, 1:3])
a_means <- rowMeans(e_matrix)
smoothScatter(a_means, d_means, xlab="Average Intensity", ylab="Ratio", main=paste(treatment,"vs",control, "(After VSN normalization)"))
points(a_means[rownames(imp_genes_exprs)], d_means[rownames(imp_genes_exprs)], pch=8,col="red")
points(a_means[as.character(sam_results[["Name"]])], d_means[as.character(sam_results[["Name"]])], pch=20,col="green")

save(e_matrix, eset, file=paste(out_path, out_name,".RData", sep=""))

dev.off()