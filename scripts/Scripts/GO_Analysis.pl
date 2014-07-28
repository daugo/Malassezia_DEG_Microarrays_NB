#!/usr/bin/perl
use warnings;
use strict;
use DBD::mysql;

&topGo_analysis('molecular_function', 20, 10, '4DST80', 'M.furfur');

sub topGo_analysis { 
	my $ontology_group = shift
	my $table_show = shift
	my $graph_show_nodes = shift
	my $treatment  = shift
	my $control = shift

	%ontology_group_ID = (
		'molecular_function' => 'MF',
		'biological_process' => 'BP',
		'cellular_component' => 'CC',
		);

	my $topGO_script = "
rm(list=ls())
library(topGO)
library(RMySQL)
con <- dbConnect(MySQL(), user='malassezia', password='malassezia',dbname='Malassezia_JGI', host='localhost')
geneID2GO <- readMappings(file='GOannTopGO_$ontology_group.txt')
geneNames <- names(geneID2GO)
interesting_Genes <- dbGetQuery(con,\"SELECT name, d_value FROM SAM_Results WHERE treatment = '$treatment' AND control = '$control'\")
myInterestingGenes <- interesting_Genes\$name_ID
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
GOdata <- new('topGOdata', ontology = '$ontology_group_ID{$ontology_group}', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#=====
test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'Fisher test')
resultFis <- getSigGroups(GOdata, test.stat)

test.stat <- new('weightCount', testStatistic = GOFisherTest, name = 'Fisher test', sigRatio = 'ratio')
resultWeight <- getSigGroups(GOdata, test.stat)

pvalFis <- score(resultFis)
pvalWeight <- score(resultWeight, whichGO = names(pvalFis))

cat ('Correlation fisher vs weigths')
cor(pvalFis, pvalWeight)

allRes <- GenTable(GOdata, classic = resultFis, weight = resultWeight,orderBy = 'weight', ranksOf = 'classic', topNodes = topNo)
write.table(allRes,file='GO_Analysis$treatment\_vs\_$control\_$ontology_group\_results.tab',sep='\t',quote=F)";

	open my $fh,'>',"GO_Analysis$treatment\_vs\_$control\_$ontology_group\_results.R" or die "Cannot create R script";
	print $fh "$top_GO_script";
	close $fh;
	!system "Rscript GO_Analysis$treatment\_vs\_$control\_$ontology_group\_results.R 2> GO_Analysis$treatment\_vs\_$control\_$ontology_group\_results.report" or die "Cannot execute Rscript";
	
return 1;

}