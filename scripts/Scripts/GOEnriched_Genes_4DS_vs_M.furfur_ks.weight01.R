rm(list=ls())
library(RMySQL)
library(GO.db)

con <- dbConnect(MySQL(),
                 user='malassezia',
                 password='malassezia',
                 dbname='Malassezia_JGI',
                 host='localhost')

treatment = '4DS'
control = 'M.furfur'
test = 'ks'
algorithm = 'weight01'
DEG_Type_Analysis = 'SAM'

p_value = 0.005
outfile_table = 'GOEnriched_Genes_4DS_vs_M.furfur_ks.weight01.txt'

GO_sig_Terms <- dbGetQuery(con, paste("SELECT TopGO_Analysis.GO_ID FROM TopGO_Analysis WHERE treatment ='",
                                      treatment,
                                      "' AND control ='",
                                      control,
                                      "' AND algorithm ='",
                                      algorithm,
                                      "' AND test ='",
                                      test,
                                      "' AND DEG_Type_Analysis ='",
                                      DEG_Type_Analysis,
                                      "' ORDER BY  p_value",sep=''))


get_offspring = function(x) {
  if(Ontology(get(x,GOTERM)) == "MF") { 
    get(x,GOMFOFFSPRING) 
  }
  else if(Ontology(get(x,GOTERM)) == "BP") {
    get(x,GOBPOFFSPRING)
  }
  else {
    get(x,GOCCOFFSPRING)
  }
}



related_genes = function(x) dbGetQuery(con, paste("SELECT DISTINCT(Genes.name) FROM Go_Annotation JOIN Genes ON Go_Annotation.transcriptId = Genes.transcriptId JOIN SAM_Results ON SAM_Results.name = Genes.name WHERE Go_Annotation.goAcc ='",x, "' AND SAM_Results.treatment ='", treatment ,"' AND SAM_Results.control ='", control,"' AND p_value <", p_value, sep=''))

get_offspring_genes <- function(x) {
    l <- lapply(get_offspring(x), related_genes)
    names(l) <- get_offspring(x)
    return(Filter(function(x) length(x)>0, l))
}

GO_terms <- NA
GO_offspring <- NA
GO_genes <- NA
index <- 0
for (i in GO_sig_Terms[,1]) {
  if (length(related_genes(i)) > 0) {
    for (m in related_genes(i)[,1]) {
      GO_terms[index + 1] = i
      GO_offspring[index + 1] = i
      GO_genes[index+1] = m
      index = index + 1
    }
  }
  offspring_genes <- get_offspring_genes(i)
  for (j in get_offspring(i)) {
    if (length(related_genes(j)) > 0) {
      for (m in get(j,offspring_genes)[,1]){
      	GO_terms[index + 1] = i
      	GO_offspring[index + 1] = j
      	GO_genes[index+1] = m
      	index = index + 1
      }
    }
  }
}

df_alldata <- data.frame(GO_term=GO_terms, GO_offspring=GO_offspring, Gene=GO_genes)
write.table(df_alldata,file=outfile_table,sep='	')



