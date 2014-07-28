#install.packages("~/Documents/Yulien_data/oligo_input_files/pd.110413.mglo.cag.exp",repos=NULL, type = "source")

library(oligo)
library("pd.110413.mglo.cag.exp")
library(genefilter)
library(limma)

#extdata <- system.file("extdata", package="pd.110413.mglo.cag.exp")
xys.files <- list.xysfiles("~/Documents/Yulien_data/Raw_Data_Files/Pair_Files",full.names=TRUE)
theData <- data.frame(Key=rep(c("brain", "universal reference"), each=6))
rownames(theData) <- basename(xys.files)
lvls <- c("channel1", "channel2", "_ALL_")
vMtData <- data.frame(channel=factor("_ALL_", levels=lvls),labelDescription="Sample type")
pd <- new("AnnotatedDataFrame", data=theData, varMetadata=vMtData)
maqc <- read.xysfiles(xys.files, phenoData=pd)
class(maqc)
exprs(maqc)[10001:10010, 1:2]

boxplot(maqc, main="MAQC Sample Data")
hist(maqc, main="MAQC Sample Data")

#RMA algorithm
eset <- rma(maqc)
class(eset)
show(eset)
exprs(eset)[1:10, 1:2]
boxplot(eset, transfo=identity, main="After RMA")
hist(eset, transfo=identity, main="After RMA")


#Assesing differential expression
e <- exprs(eset)
index <- which(eset[["Key"]] == "brain")
d <- rowMeans(e[, index])-rowMeans(e[, -index])
a <- rowMeans(e)
sum(abs(d)>1)
tt <- rowttests(e, factor(eset[["Key"]]))
lod <- -log10(tt[["p.value"]])
smoothScatter(a, d, xlab="Average Intensity", ylab="Log-ratio", main="MAQC Sample")


design <- model.matrix(~factor(eset[["Key"]]))
fit <- lmFit(eset, design)
ebayes <- eBayes(fit)
lod <- -log10(ebayes[["p.value"]][,2])
mtstat<- ebayes[["t"]][,2]
o1 <- order(abs(d), decreasing=TRUE)[1:25]
o2 <- order(abs(mtstat), decreasing=TRUE)[1:25]
o <- union(o1, o2)
smoothScatter(d, lod, main="Moderated t", xlab="Log-ratio", ylab="LOD")
points(d[o1], lod[o1], pch=18,col="blue")
points(d[o2], lod[o2], pch=1,col="red")
abline(h=2, v=c(-1, 1))