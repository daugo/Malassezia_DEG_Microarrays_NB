rm(list=ls())
library(topGO)
library(RMySQL)

con <- dbConnect(MySQL(),
                 user='malassezia',
                 password='malassezia',
                 dbname='Malassezia_JGI',
                 host='localhost')

#Parameters
treatment <- "M.furfur"
control <- "M.globosa"

file_geneID2GO <- "~/Dropbox/Malassezia_2013_01/Yulien_data/Scripts/GOannTopGO_biological_process.txt"
geneID2GO <- readMappings(file=file_geneID2GO)
all_genesScores_SAM_df <- dbGetQuery(con,"SELECT name, p_value  FROM SAM_Results WHERE treatment = 'M.furfur' AND control= 'M.globosa'")
all_genesScores_Ebayes_df <- dbGetQuery(con,"SELECT name, adj_p_value  FROM Ebayes_Results WHERE treatment = 'M.furfur' AND control= 'M.globosa'") 

all_genesScores_SAM <- all_genesScores_SAM_df$p_value
all_genesScores_Ebayes <- all_genesScores_Ebayes_df$adj_p_value

names(all_genesScores_SAM) <- all_genesScores_SAM_df$name
names(all_genesScores_Ebayes) <- all_genesScores_Ebayes_df$name

topDiffGenes_SAM <- function(allScore) {
  return (allScore < 0.005)
}

topDiffGenes_Ebayes <- function(allScore) {
  return (allScore < 0.005)
}

GOdata_SAM <- new("topGOdata",
                  description = paste('GO Analysis for SAM Results:',treatment ,'vs' ,control,'.'),
                  ontology = "BP",
                  allGenes = all_genesScores_SAM,
                  geneSel = topDiffGenes_SAM,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
                  nodeSize = 5)

GOdata_Ebayes <- new("topGOdata",
                  description = paste('GO Analysis for Ebayes Results:',treatment ,'vs' ,control,'.'),
                  ontology = "BP",
                  allGenes = all_genesScores_Ebayes,
                  geneSel = topDiffGenes_Ebayes,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
                  nodeSize = 5)

description(GOdata_SAM)
description(GOdata_Ebayes)

test.stat_classicFisher <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher tests")
test.stat_elimFisher <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher tests")
test.stat_weightFisher <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher tests")
test.stat_weight01Fisher <- new("weight01Count", testStatistic = GOFisherTest, name = "Fisher tests")
test.stat_leaFisher <- new("leaCount", testStatistic = GOFisherTest, name = "Fisher tests")
test.stat_parentchildFisher <- new("parentChild", testStatistic = GOFisherTest, name = "Fisher tests")

test.stat_classicKS <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
test.stat_elimKS  <- new("elimScore", testStatistic = GOKSTest, name = "KS tests")
test.stat_weight01KS <- new("weight01Score", testStatistic = GOKSTest, name = "KS tests")
test.stat_leaKS <- new("leaScore", testStatistic = GOKSTest, name = "KS tests")

test.stat_classict <- new("classicScore", testStatistic = GOtTest, name = "t tests")
test.stat_elimt  <- new("elimScore", testStatistic = GOtTest, name = "t tests")
test.stat_weight01t <- new("weight01Score", testStatistic = GOtTest, name = "t tests")
test.stat_leat <- new("leaScore", testStatistic = GOtTest, name = "t tests")


#Fisher
result_classicFisher_SAM <- getSigGroups(GOdata_SAM, test.stat_classicFisher)
result_elimFisher_SAM <- getSigGroups(GOdata_SAM, test.stat_elimFisher)
result_weightFisher_SAM <- getSigGroups(GOdata_SAM, test.stat_weightFisher)
result_weight01Fisher_SAM <- getSigGroups(GOdata_SAM, test.stat_weight01Fisher)
result_leaFisher_SAM <- getSigGroups(GOdata_SAM, test.stat_leaFisher)
result_parentchildFisher_SAM <- getSigGroups(GOdata_SAM, test.stat_parentchildFisher)

result_classicFisher_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_classicFisher)
result_elimFisher_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_elimFisher)
result_weightFisher_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_weightFisher)
result_weight01Fisher_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_weight01Fisher)
result_leaFisher_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_leaFisher)
result_parentchildFisher_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_parentchildFisher)

#KS
result_classicKS_SAM <- getSigGroups(GOdata_SAM, test.stat_classicKS)
result_elimKS_SAM <- getSigGroups(GOdata_SAM, test.stat_elimKS)
result_weight01KS_SAM <- getSigGroups(GOdata_SAM, test.stat_weight01KS)
result_leaKS_SAM <- getSigGroups(GOdata_SAM, test.stat_leaKS)

result_classicKS_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_classicKS)
result_elimKS_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_elimKS)
result_weight01KS_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_weight01KS)
result_leaKS_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_leaKS)

#t
result_classict_SAM <- getSigGroups(GOdata_SAM, test.stat_classict)
result_elimt_SAM <- getSigGroups(GOdata_SAM, test.stat_elimt)
result_weight01t_SAM <- getSigGroups(GOdata_SAM, test.stat_weight01t)
result_leat_SAM <- getSigGroups(GOdata_SAM, test.stat_leat)

result_classict_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_classict)
result_elimt_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_elimt)
result_weight01t_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_weight01t)
result_leat_Ebayes <- getSigGroups(GOdata_Ebayes, test.stat_leat)



