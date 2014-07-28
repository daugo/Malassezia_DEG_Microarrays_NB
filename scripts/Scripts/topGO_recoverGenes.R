rm(list=ls())
library(RMySQL)
library(GO.db)

con <- dbConnect(MySQL(),
                 user='malassezia',
                 password='malassezia',
                 dbname='Malassezia_JGI',
                 host='localhost')

treatment = '4DST80'
control = '4DS'
algorithm = 'weight01'
test = 'ks'
DEG_Type_Analysis = 'SAM'

p_value = '0.005'

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
                                      "'",sep=''))


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
#offspring <- apply(GO_sig_Terms,1:2,get_offspring)
#names(offspring) <- GO_sig_Terms$GO_ID

#Annotated_number = function(x) dbGetQuery(con, paste("SELECT COUNT(goAcc) FROM Go_Annotation WHERE goAcc ='", x,"'",sep=''))

related_genes = function(x) dbGetQuery(con, paste("SELECT DISTINCT(Genes.name) FROM Go_Annotation JOIN Genes ON Go_Annotation.transcriptId = Genes.transcriptId JOIN SAM_Results ON SAM_Results.name = Genes.name WHERE Go_Annotation.goAcc ='",x, "' AND SAM_Results.treatment ='", treatment ,"' AND SAM_Results.control ='", control,"' AND p_value <", p_value, sep=''))

get_offspring_genes <- function(x) {
    l <- lapply(get_offspring(x), related_genes)
    names(l) <- get_offspring(x)
    return(Filter(function(x) length(x)>0, l))
}
#GO_Genes <- lapply(GO_sig_Terms[,1], get_offspring_present)
#names(GO_Genes) <- GO_sig_Terms[,1]


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
df_GO_Gene <- data.frame(GO_term=GO_terms, Gene=GO_genes)
df_GO_Gene <- df_GO_Gene[!duplicated(df_GO_Gene), ]
write.table(df,file='~/Dropbox/Malassezia_2013_01/Yulien_data/Scripts/GO_Gene.txt',sep='\t')

