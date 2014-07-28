#!/usr/bin/perl
use warnings;
use strict;
use DBD::mysql;

#REMIENDER: I was in a rush but this with the hash global definition is Super-BAD coding

#=========== Database login info =============
my $dsn = "DBI:mysql:Malassezia_JGI:localhost";
my $user = "malassezia";
my $passwd = "malassezia";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#===================================

#=======Hashes that map Genes with their GO IDS supported by the differente annotations used=======
my %geneName_goId_M;
my %geneName_goId_P;
my %geneName_goId_C;

#=======SQL Queries that retrive the mapping for the available annotations that I have for M.globosa===
#JGI annotation (Performed 2009)
my $get_GO_Annotation_sql = "SELECT Genes.name, GeneOntologyInfo.goAcc 
						FROM Genes JOIN GeneOntologyInfo ON 
							Genes.transcriptId = GeneOntologyInfo.transcriptId 
							WHERE gotermType = 'GOTERMTYPE'";

#B2GO automatic annotation (I performed this annotation It's up to dated for 2013)
my $get_B2GO_annot_sql = "SELECT Genes.name, Blast2GO_Annotation.goAcc
						FROM Genes JOIN Blast2GO_Annotation ON
							Genes.transcriptId = Blast2GO_Annotation.transcriptId
							WHERE go_group = 'GOTERMTYPE'";

#===================================================================

&get_annotation_topGO($get_GO_Annotation_sql, 'molecular_function');
&get_annotation_topGO($get_GO_Annotation_sql, 'biological_process');
&get_annotation_topGO($get_GO_Annotation_sql, 'cellular_component');

&get_annotation_topGO($get_B2GO_annot_sql, 'molecular_function');
&get_annotation_topGO($get_B2GO_annot_sql, 'biological_process');
&get_annotation_topGO($get_B2GO_annot_sql, 'cellular_component');

&Print_Annotation('molecular_function', \%geneName_goId_M);
&Print_Annotation('biological_process', \%geneName_goId_P);
&Print_Annotation('cellular_component', \%geneName_goId_C);

sub Print_Annotation {
	my ($gotermType, $geneName_goId_ref) = @_;

	open OUT, '>', 'GOannTopGO_'.$gotermType.'.txt';

	foreach my $gene_name (keys %{$geneName_goId_ref}) {
		my $go_list = join(', ',@{$geneName_goId_ref->{$gene_name}});
		print OUT $gene_name.'	'.$go_list."\n";
	}
	close OUT;
}

sub get_annotation_topGO {
	my $get_Ann_sql = shift;
	my $gotermType = shift;
	
	$get_Ann_sql =~ s/GOTERMTYPE/$gotermType/;
	my $get_Ann_sth = $dbh->prepare($get_Ann_sql);
	$get_Ann_sth->execute();


	while (my $ref = $get_Ann_sth->fetchrow_arrayref()) {
		my ($gene_name, $go_id) = ($ref->[0], $ref->[1]);
		if ($gotermType eq 'molecular_function') {
			push( @{$geneName_goId_M{$gene_name}}, $go_id) unless (grep {$_ eq $go_id} @{$geneName_goId_M{$gene_name}});
		}
		elsif ($gotermType eq 'biological_process') {
			push( @{$geneName_goId_P{$gene_name}}, $go_id) unless (grep {$_ eq $go_id} @{$geneName_goId_P{$gene_name}});
		}
		elsif ($gotermType eq 'cellular_component') {
			push( @{$geneName_goId_C{$gene_name}}, $go_id) unless (grep {$_ eq $go_id} @{$geneName_goId_C{$gene_name}});
		}
	}

}


$dbh->disconnect;