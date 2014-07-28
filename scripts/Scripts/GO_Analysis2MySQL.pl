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
my $questions = '?'x10;
my $ques_sep = join(', ', split('', $questions));
my $insert_goAnalysis_sql = "INSERT INTO GO_Analysis(treatment, control, Rank_weight, GO_ID, Term, Term_type, Annotated, Significant, Expected, Rank_classic, p_classic, p_weigth) VALUES ($ques_sep)";
my $insert_goAnalysis_sth = $dbh->prepare($insert_goAnalysis_sql);


my @files = @ARGV;

&parse_topGO_tab_file($_,$insert_goAnalysis_sth) foreach (@files);

sub parse_topGO_tab_file {



	my $file = shift;
	my $insert_sth = shift;
	my ($fn,$dir,$suf) = fileparse($file,qr/\.[^.]*/);
	my ($treatment, $control, $term_type_av) = $fn =~ /^GO_Analysis([^_]+)_vs_([^_]+)_([^_]+)/;
	my %typeAv_type = (
			'BP' => 'biological_process',
			'CC' => 'cellular_component',
			'MF' => 'molecular_function',
			);
	
	open IN,'<',"$file" or die "Cannot open $file";

	while(<IN>) {
		chomp;
		next if $_ =~ /^GO.ID/;
		my (@fields) = split('\t',$_);
		#print "$treatment, $control=> @fields\n";
		my ($Rank_weight, $GO_ID, $term, $Annotated, $Significant, $Expected, $Rank_classic, $p_classic, $p_weight) = @fields;
		my $check_GoId_sql = "SELECT goAcc FROM GeneOntologyInfo WHERE goAcc = '$GO_ID'";
		if (recover_info($check_GoId_sql)) {
			my $check_GoAnalysisEntry_sql = "SELECT COUNT(*) FROM GO_Analysis WHERE treatment = '$treatment' AND control = '$control' AND GO_ID = '$GO_ID'";
			unless (check_existence($check_GoAnalysisEntry_sql)) { 
				$insert_sth->execute($treatment, $control, $Rank_weight, $GO_ID, $term, $typeAv_type{$term_type_av}, $Annotated, $Significant, $Expected, $Rank_classic, $p_classic, $p_weight);
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
