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

#Insert Genes table
my $questions = '?'x12;
my $ques_sep = join(', ', split('', $questions));
my $insert_TopGO_Analysis_sql = "INSERT INTO TopGO_Analysis(treatment, control, GO_ID, Term, Term_type , Annotated, Significant, Expected, test, algorithm, p_value, DEG_Type_Analysis) VALUES ($ques_sep)";
my $insert_TopGO_Analysis_sth = $dbh->prepare($insert_TopGO_Analysis_sql);


my @files = @ARGV;

foreach (@files) {
	&parse_topGO_tab_file($_,$insert_TopGO_Analysis_sth);
}


sub parse_topGO_tab_file {



	my $file = shift;
	my $insert_sth = shift;
	my ($fn,$dir,$suf) = fileparse($file,qr/\.[^.]*/);
	my ($treatment, $control, $term_type_av, $DEG_type) = $fn =~ /^GeneTable_([^_]+)_vs_([^_]+)_([^_]+)_([^_]+)$/;
	my %typeAv_type = (
                        'BP' => 'biological_process',
                        'CC' => 'cellular_component',
                        'MF' => 'molecular_function',
                        );
 
	open IN,'<',"$file" or die "Cannot open $file";
	
	my $test;
	my $algorithm;
	while(<IN>) {
		chomp;
		$_ =~ s/"//g;
		if ($_ =~ /^GO.ID/) {
			($test, $algorithm) = $_ =~ /\t([A-Z0-9]+).([A-Z0-9]+)$/i;
		}
		else {
			my (@fields) = split('\t',$_);
			#print "$treatment, $control=> @fields\n";
			my (undef, $GO_ID, $term, $Annotated, $Significant, $Expected, $p_value) = @fields;
			my $check_GoAnalysisEntry_sql = "SELECT COUNT(*) FROM TopGO_Analysis WHERE treatment = '$treatment' AND control = '$control' AND GO_ID = '$GO_ID' AND 'test' = '$test' AND algorithm = '$algorithm' AND DEG_Type_Analysis = '$DEG_type'";
			unless (check_existence($check_GoAnalysisEntry_sql)) { 
				$insert_sth->execute($treatment, $control, $GO_ID, $term, $typeAv_type{$term_type_av}, $Annotated, $Significant, $Expected, $test, $algorithm, $p_value, $DEG_type);
			}
		}
	}

	
}

sub check_existence { #send COUNT(*) to check if register is new
        my $check_if_new = $dbh->selectrow_array($_[0]); #if it's already in the table value = 1
        if ($check_if_new) {
                return 1;
        }
        else {
                return 0;
        }
}

sub recover_info { #recover scalar info from the DB
        my ($recover_info) = $dbh->selectrow_array($_[0]);
        return $recover_info;
}
