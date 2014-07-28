#Create annotation package
library(pdInfoBuilder)
ndf <- "~/Desktop/Documents_Mac/Yulien_data/Design_information/110413_Mglo_CAG_exp.ndf"
xys <- "~/Desktop/Documents_Mac/Yulien_data/Raw_Data_Files/Pair_Files/559446A01_532.xys"

seed <- new("NgsExpressionPDInfoPkgSeed",
            ndfFile = ndf, xysFile = xys,
            author = "David Urbina",
            email = "daugo182@gmail.com",
            biocViews = "AnnotationData",
            genomebuild = "NZ_AAYY00000000.1",
            organism = "Malassezia", species = "Malassezia globosa")

makePdInfoPackage(seed, destDir = "~/Desktop/Documents_Mac/Yulien_data/oligo_input_files/")