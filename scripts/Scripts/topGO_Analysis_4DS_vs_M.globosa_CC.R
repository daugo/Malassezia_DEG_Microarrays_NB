
rm(list=ls())
library(topGO)
library(RMySQL)

con <- dbConnect(MySQL(),
                 user='malassezia',
                 password='malassezia',
                 dbname='Malassezia_JGI',
                 host='localhost')

#Parameters
treatment <- "4DS"
control <- "M.globosa"
annotation_file <- '~/Dropbox/Malassezia_2013_01/Yulien_data/Scripts/GOannTopGO_cellular_component.txt'
p_value_SAM <- 0.005
p_value_Ebayes <- 0.01
ontology <- 'CC'
nodeSize <- 5
nodes_table = 30
nodes_graph = 5
table_outfile_SAM <- 'GeneTable_4DS_vs_M.globosa_CC_SAM.txt'
table_outfile_Ebayes <- 'GeneTable_4DS_vs_M.globosa_CC_Ebayes.txt'
graph_outfile_SAM <- 'Graph_4DS_vs_M.globosa_CC_SAM.pdf'
graph_outfile_Ebayes <- 'Graph_4DS_vs_M.globosa_CC_Ebayes.pdf'


file_geneID2GO <- annotation_file
geneID2GO <- readMappings(file=file_geneID2GO)
all_genesScores_SAM_df <- dbGetQuery(con,"SELECT name, p_value  FROM SAM_Results WHERE treatment = '4DS' AND control= 'M.globosa'")
all_genesScores_Ebayes_df <- dbGetQuery(con,"SELECT name, adj_p_value  FROM Ebayes_Results WHERE treatment = '4DS' AND control= 'M.globosa'") 

all_genesScores_SAM <- all_genesScores_SAM_df$p_value
all_genesScores_Ebayes <- all_genesScores_Ebayes_df$adj_p_value

names(all_genesScores_SAM) <- all_genesScores_SAM_df$name
names(all_genesScores_Ebayes) <- all_genesScores_Ebayes_df$name

topDiffGenes_SAM <- function(allScore) {
  return (allScore < p_value_SAM)
}

topDiffGenes_Ebayes <- function(allScore) {
  return (allScore < p_value_Ebayes)
}

GOdata_SAM <- new("topGOdata",
                  description = paste('GO Analysis for SAM Results:',treatment ,'vs' ,control,'.'),
                  ontology = ontology,
                  allGenes = all_genesScores_SAM,
                  geneSel = topDiffGenes_SAM,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
                  nodeSize = nodeSize)

GOdata_Ebayes <- new("topGOdata",
                  description = paste('GO Analysis for Ebayes Results:',treatment ,'vs' ,control,'.'),
                  ontology = ontology,
                  allGenes = all_genesScores_Ebayes,
                  geneSel = topDiffGenes_Ebayes,
                  annot = annFUN.gene2GO,
                  gene2GO = geneID2GO,
                  nodeSize = nodeSize)

description(GOdata_SAM)
description(GOdata_Ebayes)

result_ks.classic_SAM <- runTest(GOdata_SAM, algorithm='classic', statistic='ks')
result_ks.classic_Ebayes <- runTest(GOdata_Ebayes, algorithm='classic', statistic='ks')
result_ks.elim_SAM <- runTest(GOdata_SAM, algorithm='elim', statistic='ks')
result_ks.elim_Ebayes <- runTest(GOdata_Ebayes, algorithm='elim', statistic='ks')
result_ks.weight01_SAM <- runTest(GOdata_SAM, algorithm='weight01', statistic='ks')
result_ks.weight01_Ebayes <- runTest(GOdata_Ebayes, algorithm='weight01', statistic='ks')
result_ks.lea_SAM <- runTest(GOdata_SAM, algorithm='lea', statistic='ks')
result_ks.lea_Ebayes <- runTest(GOdata_Ebayes, algorithm='lea', statistic='ks')

result_t.classic_SAM <- runTest(GOdata_SAM, algorithm='classic', statistic='t')
result_t.classic_Ebayes <- runTest(GOdata_Ebayes, algorithm='classic', statistic='t')
result_t.elim_SAM <- runTest(GOdata_SAM, algorithm='elim', statistic='t')
result_t.elim_Ebayes <- runTest(GOdata_Ebayes, algorithm='elim', statistic='t')
result_t.weight01_SAM <- runTest(GOdata_SAM, algorithm='weight01', statistic='t')
result_t.weight01_Ebayes <- runTest(GOdata_Ebayes, algorithm='weight01', statistic='t')
result_t.lea_SAM <- runTest(GOdata_SAM, algorithm='lea', statistic='t')
result_t.lea_Ebayes <- runTest(GOdata_Ebayes, algorithm='lea', statistic='t')

pdf(graph_outfile_SAM)
GenTable_ks.classic_SAM <- GenTable(GOdata_SAM, ks.classic = result_ks.classic_SAM,  topNodes=nodes_table)
write.table(GenTable_ks.classic_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_ks.classic_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_ks.elim_SAM <- GenTable(GOdata_SAM, ks.elim = result_ks.elim_SAM,  topNodes=nodes_table)
write.table(GenTable_ks.elim_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_ks.elim_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_ks.weight01_SAM <- GenTable(GOdata_SAM, ks.weight01 = result_ks.weight01_SAM,  topNodes=nodes_table)
write.table(GenTable_ks.weight01_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_ks.weight01_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_ks.lea_SAM <- GenTable(GOdata_SAM, ks.lea = result_ks.lea_SAM,  topNodes=nodes_table)
write.table(GenTable_ks.lea_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_ks.lea_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.classic_SAM <- GenTable(GOdata_SAM, t.classic = result_t.classic_SAM,  topNodes=nodes_table)
write.table(GenTable_t.classic_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_t.classic_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.elim_SAM <- GenTable(GOdata_SAM, t.elim = result_t.elim_SAM,  topNodes=nodes_table)
write.table(GenTable_t.elim_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_t.elim_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.weight01_SAM <- GenTable(GOdata_SAM, t.weight01 = result_t.weight01_SAM,  topNodes=nodes_table)
write.table(GenTable_t.weight01_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_t.weight01_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.lea_SAM <- GenTable(GOdata_SAM, t.lea = result_t.lea_SAM,  topNodes=nodes_table)
write.table(GenTable_t.lea_SAM, file=table_outfile_SAM, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_SAM, score(result_t.lea_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))
dev.off()
pdf(graph_outfile_Ebayes)
GenTable_ks.classic_Ebayes <- GenTable(GOdata_Ebayes, ks.classic = result_ks.classic_Ebayes, topNodes= nodes_table)
write.table(GenTable_ks.classic_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_ks.classic_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_ks.elim_Ebayes <- GenTable(GOdata_Ebayes, ks.elim = result_ks.elim_Ebayes, topNodes= nodes_table)
write.table(GenTable_ks.elim_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_ks.elim_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_ks.weight01_Ebayes <- GenTable(GOdata_Ebayes, ks.weight01 = result_ks.weight01_Ebayes, topNodes= nodes_table)
write.table(GenTable_ks.weight01_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_ks.weight01_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_ks.lea_Ebayes <- GenTable(GOdata_Ebayes, ks.lea = result_ks.lea_Ebayes, topNodes= nodes_table)
write.table(GenTable_ks.lea_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_ks.lea_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.classic_Ebayes <- GenTable(GOdata_Ebayes, t.classic = result_t.classic_Ebayes, topNodes= nodes_table)
write.table(GenTable_t.classic_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_t.classic_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.elim_Ebayes <- GenTable(GOdata_Ebayes, t.elim = result_t.elim_Ebayes, topNodes= nodes_table)
write.table(GenTable_t.elim_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_t.elim_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.weight01_Ebayes <- GenTable(GOdata_Ebayes, t.weight01 = result_t.weight01_Ebayes, topNodes= nodes_table)
write.table(GenTable_t.weight01_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_t.weight01_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
GenTable_t.lea_Ebayes <- GenTable(GOdata_Ebayes, t.lea = result_t.lea_Ebayes, topNodes= nodes_table)
write.table(GenTable_t.lea_Ebayes, file=table_outfile_Ebayes, sep='\t',append=TRUE)
try(showSigOfNodes(GOdata_Ebayes, score(result_t.lea_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))
dev.off()
