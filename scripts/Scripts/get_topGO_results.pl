#!/usr/bin/perl
use warnings;
use strict;
use DBD::mysql;
use File::Basename;

#=========== Database login info =============
my $dsn = "DBI:mysql:Malassezia_JGI:localhost";
my $user = "malassezia";
my $passwd = "malassezia";

my $dbh = DBI->connect($dsn,$user,$passwd) or die "Cannot connect to DB!\n"; #conection to the DB
#===================================

my @treatments = ("4DST80", "4DST80",  "4DS", "4DS", "M.furfur");
my @controls =  ("4DS", "4DS", "M.furfur", "M.globosa", "M.globosa");
my $test = 'ks';
my $algorithm = 'weight01';
my @term_types = ('biological_process', 'molecular_function', 'cellular_component');

for (my $i = 0; $i < @treatments; $i++) {
	foreach my $term_type (@term_types) { 
		my $get_terms_sql = "SELECT TopGO_Analysis.treatment,
						 		TopGO_Analysis.control,
						 		Go_Annotation.goAcc,
						 		Go_Annotation.term,
						 		Go_Annotation.goTermType,
						 		TopGO_Analysis.test,
						 		TopGO_Analysis.algorithm,
						 		TopGO_Analysis.p_value 
						 	FROM Go_Annotation JOIN TopGO_Analysis ON Go_Annotation.goAcc = TopGO_Analysis.GO_ID 
						 		WHERE  TopGO_Analysis.treatment = '$treatments[$i]' AND 
						 		TopGO_Analysis.control = '$controls[$i]'  AND 
						 		TopGO_Analysis.test = '$test' AND 
						 		TopGO_Analysis.algorithm = '$algorithm' AND 
						 		Go_Annotation.goTermType = '$term_type' AND 
						 		TopGO_Analysis.DEG_Type_Analysis = 'SAM' 

						 		GROUP BY  Go_Annotation.goAcc  
						 		ORDER BY  TopGO_Analysis.p_value";
		my $get_terms_sth = $dbh->prepare($get_terms_sql);
		$get_terms_sth->execute();
		
		print "Treatment\tControl\tgoAcc\tTerm\tTermType\tTest\tAlgorithm\tp_value\n";
		while(my $ref= $get_terms_sth->fetchrow_arrayref()) {
			my $line = join("\t", @{$ref});
			#my ($treatment, $control, $goAcc, $term, $goTermType, $test, $algorithm, $p_value) = @{$ref};
			print $line."\n";
		}
		print "\n";
	}
	print "\n";
}
	