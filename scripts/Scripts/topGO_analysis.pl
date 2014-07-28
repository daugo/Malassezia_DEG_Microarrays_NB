#!/usr/bin/perl
use warnings;
use strict;


#treatments = ("4DST80","4DST80", "4DST80", "4DS", "4DS", "M.furfur")
#controls =  ("4DS", "M.furfur", "M.globosa", "M.furfur", "M.globosa", "M.globosa")
my $treatments = shift @ARGV;
my $controls = shift @ARGV;
my $p_value_cutoff_SAM = shift @ARGV;
my $p_value_cutoff_Ebayes = shift @ARGV;
my $nodeSize = shift @ARGV;
my $nodes_table = shift @ARGV;
my $nodes_graph = shift @ARGV;

&Manage_Jobs($treatments, $controls);
#&Manage_Jobs($treatments, $controls, $p_value_cutoff_SAM, $p_value_cutoff_Ebayes, $nodeSize, $nodes_table, $nodes_graph);

sub Manage_Jobs { 
	my @treatments = split(',',shift);
	my @controls =   split(',',shift);
	my %annnotation_file_ontology = (
			'~/Dropbox/Malassezia_2013_01/Yulien_data/Scripts/GOannTopGO_molecular_function.txt' => 'MF',
			'~/Dropbox/Malassezia_2013_01/Yulien_data/Scripts/GOannTopGO_biological_process.txt' => 'BP',
			'~/Dropbox/Malassezia_2013_01/Yulien_data/Scripts/GOannTopGO_cellular_component.txt' => 'CC',
		);
	#my $p_value_cutoff_SAM = shift;
	#my $p_value_cutoff_Ebayes = shift;
	#my $nodeSize = shift;
	#my $nodes_table = shift;
	#my $nodes_graph = shift;

	for (my $i = 0; $i < @treatments; $i++) {
		foreach my $anno_file (keys %annnotation_file_ontology) { 
			&TopGO_job($treatments[$i],
					$controls[$i],
					$anno_file,
					$annnotation_file_ontology{$anno_file},
					);
		}
		
	}
}


sub TopGO_job {

	my $db_user = 'malassezia';
	my $db_psw = 'malassezia';
	my $db_name = 'Malassezia_JGI';
	my $db_host = 'localhost';

	my $treatment = shift;
	my $control = shift;

	my $annnotation_file = shift;
	my $ontology = shift;

	#my $p_value_cutoff_SAM = shift;
	#my $p_value_cutoff_Ebayes = shift;
	#my $nodeSize = shift;

	my %test_algorithms = (
			#"fisher" => ["classic", "elim", "weight", "weight01", "lea", "parentchild"],
			"ks" => ["classic", "elim", "weight01", "lea"],
			"t" => ["classic", "elim", "weight01", "lea"],
		);
	
	#my $nodes_table = shift;
	#my $nodes_graph = shift;

	my $table_outfile_base = "GeneTable\_$treatment\_vs\_$control\_$ontology\_";
	my $table_outfile_SAM = $table_outfile_base.'SAM.txt';
	my $table_outfile_Ebayes = $table_outfile_base.'Ebayes.txt';
	unlink $table_outfile_SAM if (-e $table_outfile_SAM);
	unlink $table_outfile_Ebayes if (-e $table_outfile_Ebayes);

	my $graph_base_name = "Graph\_$treatment\_vs\_$control\_$ontology\_";
	my $graph_outfile_SAM = $graph_base_name.'SAM.pdf';
	my $graph_outfile_Ebayes  = $graph_base_name.'Ebayes.pdf';

	my $Rscript_outfile = "topGO_Analysis_$treatment\_vs\_$control\_$ontology.R";
	my $Rscript_report1 = "topGO_Analysis_$treatment\_vs\_$control\_$ontology.report1";
	my $Rscript_report2 = "topGO_Analysis_$treatment\_vs\_$control\_$ontology.report2";

	my $RScript_start = <<EOF

rm(list=ls())
library(topGO)
library(RMySQL)

con <- dbConnect(MySQL(),
                 user='$db_user',
                 password='$db_psw',
                 dbname='$db_name',
                 host='$db_host')

#Parameters
treatment <- "$treatment"
control <- "$control"
annotation_file <- '$annnotation_file'
p_value_SAM <- $p_value_cutoff_SAM
p_value_Ebayes <- $p_value_cutoff_Ebayes
ontology <- '$ontology'
nodeSize <- $nodeSize
nodes_table = $nodes_table
nodes_graph = $nodes_graph
table_outfile_SAM <- '$table_outfile_SAM'
table_outfile_Ebayes <- '$table_outfile_Ebayes'
graph_outfile_SAM <- '$graph_outfile_SAM'
graph_outfile_Ebayes <- '$graph_outfile_Ebayes'


file_geneID2GO <- annotation_file
geneID2GO <- readMappings(file=file_geneID2GO)
all_genesScores_SAM_df <- dbGetQuery(con,"SELECT name, p_value  FROM SAM_Results WHERE treatment = '$treatment' AND control= '$control'")
all_genesScores_Ebayes_df <- dbGetQuery(con,"SELECT name, adj_p_value  FROM Ebayes_Results WHERE treatment = '$treatment' AND control= '$control'") 

all_genesScores_SAM <- all_genesScores_SAM_df\$p_value
all_genesScores_Ebayes <- all_genesScores_Ebayes_df\$adj_p_value

names(all_genesScores_SAM) <- all_genesScores_SAM_df\$name
names(all_genesScores_Ebayes) <- all_genesScores_Ebayes_df\$name

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

EOF
;
	#Result tests
	my $Rscript_tests;
	foreach my $test (keys %test_algorithms) {
		foreach my $algorithm (@{$test_algorithms{$test}}) {
			$Rscript_tests .= "result_$test.$algorithm\_SAM <- runTest(GOdata_SAM, algorithm='$algorithm', statistic='$test')\n";
			$Rscript_tests .= "result_$test.$algorithm\_Ebayes <- runTest(GOdata_Ebayes, algorithm='$algorithm', statistic='$test')\n";
		}
		$Rscript_tests .= "\n";
	}


	my $Rscript_GeneTables_SAM = '';
	my @GeneTables_SAM_args; 
	my $Rscript_GeneTables_Ebayes = '';
	my @GeneTables_Ebayes_args; 
	foreach my $test (keys %test_algorithms) {
		foreach my $algorithm (@{$test_algorithms{$test}}) {
			$Rscript_GeneTables_SAM .= "GenTable_$test.$algorithm\_SAM <- GenTable(GOdata_SAM, $test.$algorithm = result_$test.$algorithm\_SAM,  topNodes=nodes_table)\n";
			$Rscript_GeneTables_SAM .= "write.table(GenTable_$test.$algorithm\_SAM, file=table_outfile_SAM, sep='\\t',append=TRUE)\n";
			$Rscript_GeneTables_SAM .= "try(showSigOfNodes(GOdata_SAM, score(result_$test.$algorithm\_SAM), firstSigNodes = nodes_graph, useInfo = 'all'))\n";
			$Rscript_GeneTables_Ebayes .= "GenTable_$test.$algorithm\_Ebayes <- GenTable(GOdata_Ebayes, $test.$algorithm = result_$test.$algorithm\_Ebayes, topNodes= nodes_table)\n";
			$Rscript_GeneTables_Ebayes .= "write.table(GenTable_$test.$algorithm\_Ebayes, file=table_outfile_Ebayes, sep='\\t',append=TRUE)\n";
			$Rscript_GeneTables_Ebayes .= "try(showSigOfNodes(GOdata_Ebayes, score(result_$test.$algorithm\_Ebayes), firstSigNodes = nodes_graph, useInfo = 'all'))\n";
		}
	}

	my $Rscript_pdfdev_SAM = "pdf(graph_outfile_SAM)\n";
	my $Rscript_pdfdev_Ebayes = "pdf(graph_outfile_Ebayes)\n";
	my $Rscript_devoff = "dev.off()\n";

	my $Rscript = $RScript_start.$Rscript_tests.$Rscript_pdfdev_SAM.$Rscript_GeneTables_SAM.$Rscript_devoff.$Rscript_pdfdev_Ebayes.$Rscript_GeneTables_Ebayes.$Rscript_devoff;
	open OUT, '>', $Rscript_outfile or die "Cannot write in $Rscript_outfile";
	print OUT $Rscript;
	!system("Rscript $Rscript_outfile > $Rscript_report1 2> $Rscript_report2");
}


