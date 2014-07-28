#!/usr/bin/perl
use warnings;
use strict;
use DBD::mysql;

#=========== Database login info =============
my $dsn = "DBI:mysql:Malassezia_JGI:localhost";
my $user = "malassezia";
my $passwd = "malassezia";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#===================================



my @treatments = ("4DST80","4DST80", "4DST80", "4DS", "4DS", "M.furfur");
my @controls =  ("4DS", "M.furfur", "M.globosa", "M.furfur", "M.globosa", "M.globosa");
my $test = 'ks';
my $algorithm = 'weight01';
my $DEG_Type_Analysis = 'SAM';
my $p_value = "0.005"; 
my $base_output = 'GOEnriched_Genes';

my @out_files;
for(my $i=0; $i < scalar(@treatments); $i++) {
	my $outfile_fn = $base_output.'_'.$treatments[$i].'_vs_'.$controls[$i].'_'.$test.'.'.$algorithm;
	my $outfile = get_involve_genes($treatments[$i], $controls[$i], $outfile_fn);
	push(@out_files, $outfile);
}


for (my $i=0; $i < scalar(@treatments); $i++) {
	my $report_outfile = 'Report_GO_'.$treatments[$i].'_vs_'.$controls[$i].'_'.$test.'.'.$algorithm.'.txt';
	&GO_Genes_Report($out_files[$i],$treatments[$i],$controls[$i],$report_outfile);
}


sub get_involve_genes {
	my $treatment = shift;
	my $control = shift;
	my $outfile = shift;

	my $outfile_table = $outfile.'.txt';
	my $outfile_Rscript = $outfile.'.R';
	my $Rscript = 
<<EOF
rm(list=ls())
library(RMySQL)
library(GO.db)

con <- dbConnect(MySQL(),
                 user='malassezia',
                 password='malassezia',
                 dbname='Malassezia_JGI',
                 host='localhost')

treatment = '$treatment'
control = '$control'
test = '$test'
algorithm = '$algorithm'
DEG_Type_Analysis = '$DEG_Type_Analysis'

p_value = $p_value
outfile_table = '$outfile_table'

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
write.table(df_alldata,file=outfile_table,sep='\t')



EOF
;
	open OUT, '>',$outfile_Rscript or die "Cannot write into $outfile_Rscript";
	print OUT $Rscript;
	!system("Rscript $outfile_Rscript") or die "Cannot execute $outfile_Rscript";
	return $outfile_table;
	close OUT;
}

sub GO_Genes_Report {

	my $file = shift;
	my $treatment = shift;
	my $control = shift;
	my $out_file = shift;
	my %GO_genes;
	open IN,'<',$file or die "Cannot open $file";
	while (<IN>) {
		chomp;
		next if($_ =~ /"GO_term"/);
		$_ =~ s/"//g;
		my (undef,$go_term, $go_offspring, $gene) = split (/\t/,$_);
		unless (grep($gene eq $_,@{$GO_genes{$go_term}})) {
			push( @{ $GO_genes { $go_term } }, $gene); 	
		}
	}
	close IN;

	open OUT, '>', $out_file or die "Cannot write in $out_file";
	my @GO_keys_quote = map("'".$_."'", keys %GO_genes);
	print "@GO_keys_quote\n";
	my $go_ids = join(', ', @GO_keys_quote);

	my $get_TopGO_Result_sql = "SELECT treatment, control, GO_ID, Term, Term_type, Annotated, Significant, Expected, test, algorithm, p_value, DEG_Type_Analysis FROM
												TopGO_Analysis WHERE
												treatment = '$treatment' AND
												control = '$control' AND
												GO_ID IN ($go_ids) AND
												test = '$test' AND
												algorithm = '$algorithm' AND
												DEG_Type_Analysis = '$DEG_Type_Analysis'
												GROUP BY GO_ID ORDER BY Term_type, p_value;";
	my $get_TopGO_Result_sth = $dbh->prepare($get_TopGO_Result_sql);
	$get_TopGO_Result_sth->execute();

	while (my $ref = $get_TopGO_Result_sth->fetchrow_arrayref()) {
		my @fields = @{$ref};
		my $go_id = $fields[2];
		my %gene_rank;
		foreach  my $gene (@{$GO_genes{$go_id}}) {
			my $get_gene_rank_sql = "SELECT Rank FROM SAM_Results WHERE
																name = '$gene' AND
																treatment = '$treatment' AND
																control = '$control'";
			$gene_rank{$gene} = $dbh->selectrow_array($get_gene_rank_sql);
		}

		my @gene_list;
		foreach my $gene (keys %gene_rank) {
			push (@gene_list,"$gene($gene_rank{$gene})");
		}
		my $gene_str = join(', ', @gene_list);

		print OUT join("\t",(@fields,$gene_str))."\n";
	}
	close OUT;
	return $out_file;
}



