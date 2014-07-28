#Create annotation package
library(pdInfoBuilder)
#++++++ Parameters ++++++++++
ndf_file <- "~/Desktop/Documents_Mac/Yulien_data/Design_information/110413_Mglo_CAG_exp.ndf"
xys_file <- "~/Desktop/Documents_Mac/Yulien_data/Raw_Data_Files/Pair_Files/559446A01_532.xys"
dest_path <- "~/Desktop/davidurbina/Yulien_data/oligo_input_files/"

#++++++ Annotation Package Builder ++++++++++++
seed <- new("NgsExpressionPDInfoPkgSeed",
            ndfFile = ndf_file, xysFile = xys_file,
            author = "David Urbina",
            email = "daugo182@gmail.com",
            biocViews = "AnnotationData",
            genomebuild = "NZ_AAYY00000000.1",
            organism = "Malassezia", species = "Malassezia globosa")

makePdInfoPackage(seed, destDir = dest_path )