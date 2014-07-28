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

get_GO_annotation_topGO('molecular_function');
get_GO_annotation_topGO('biological_process');
get_GO_annotation_topGO('cellular_component');

sub get_GO_annotation_topGO {

	my $gotermType = shift;
	
	my $get_Ann_sql = "SELECT Genes.name, GeneOntologyInfo.goAcc FROM Genes JOIN GeneOntologyInfo ON Genes.transcriptId = GeneOntologyInfo.transcriptId WHERE gotermType = '$gotermType'";
	my $get_Ann_sth = $dbh->prepare($get_Ann_sql);
	$get_Ann_sth->execute();

	my %geneName_goId = {};
	while (my $ref = $get_Ann_sth->fetchrow_arrayref()) {
		my ($gene_name, $go_id) = ($ref->[0], $ref->[1]);
		push( @{$geneName_goId{$gene_name}}, $go_id);
	}

	open OUT, '>', 'GOannTopGO_'.$gotermType.'.txt';

	foreach my $gene_name (keys %geneName_goId) {
		my $go_list = join(', ',@{$geneName_goId{$gene_name}});
		print OUT $gene_name.'	'.$go_list."\n";
	}

$dbh->disconnect;
