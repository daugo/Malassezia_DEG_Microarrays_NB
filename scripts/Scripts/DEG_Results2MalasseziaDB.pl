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

&main;

sub main {

	#============ User Input ======================

	my @csv_files = @ARGV;

	#=========== Insert Queries ==================



	#Insert SAM_Results query

	my $ques_sep_sam = &InsertSQL_question_str(10);
	my $insert_samResults_sql = "INSERT INTO SAM_Results(
								name,
								d_value,
							 	stdev,
							 	p_value,
							 	q_value,
							 	R_fold,
							 	Name_long,
							 	treatment,
							 	control,
							 	Rank)
								VALUES ($ques_sep_sam)";
	my $insert_samResults_sth = $dbh->prepare($insert_samResults_sql);


	#Insert Ebayes_Results query

	my $ques_sep_ebayes = &InsertSQL_question_str(11);
	my $insert_ebayesResults_sql = "INSERT INTO Ebayes_Results(
								name,
								logFC,
							 	AveExpr,
							 	t,
							 	p_value,
							 	adj_p_value,
							 	B,
							 	Name_long,
							 	treatment,
							 	control,
							 	Rank)
								VALUES ($ques_sep_ebayes)";
	my $insert_ebayesResults_sth = $dbh->prepare($insert_ebayesResults_sql);

	foreach(@csv_files) { 
		if ($_ =~ /sam.csv$/) {
			#&fill_SAM_Results_Table($_,$insert_samResults_sth,$dbh);
		}
		elsif ($_ =~ /ebayes.csv$/) {
			&fill_Ebayes_Results_Table($_,$insert_ebayesResults_sth,$dbh);
		}
	}

}



sub fill_SAM_Results_Table {
        my ($sam_results_file, $insert_samexpr_sth, $dbh) = @_;
        print "Filling SAM_Results table...";
        open IN, '<', $sam_results_file or die "Cannot open $sam_results_file";
        my ($fn,$dir,$suf) = fileparse($sam_results_file, qr/\.[^.]*/);
        my ($treatment, $control) = $fn =~ /_([^_]+)vs([^_]+)_sam$/; 
        my $rank = 1;
        while (<IN>) {
                chomp;
                next if $_ =~ /^"Name",/;
                $_ =~ s/"//g;
                my ($name_full, $d_value, $stdev, $p_value, $q_value, $r_fold) = split(',');
                my $name = substr($name_full,0,8);
                my $check_genename_sql = "SELECT name FROM Genes WHERE name = '$name'";
                if (recover_info($check_genename_sql)) {
                        my $check_exprentry_sql = "SELECT COUNT(*) FROM  SAM_Results WHERE name = '$name' AND treatment = '$treatment' AND control = '$control'";
                        unless (check_existence($check_exprentry_sql)) {
                                $insert_samexpr_sth->execute($name, $d_value, $stdev, $p_value, $q_value, $r_fold, $name_full, $treatment, $control, $rank);
                        }
                }
                $rank += 1;
        }
        close IN;
        print "Done!\n";
        1;
}


sub fill_Ebayes_Results_Table {
        my ($ebayes_results_file, $insert_ebayesexpr_sth, $dbh) = @_;
        print "Filling Ebayes_Results table...";
        open IN, '<', $ebayes_results_file or die "Cannot open $ebayes_results_file";
        my ($fn,$dir,$suf) = fileparse($ebayes_results_file, qr/\.[^.]*/);
        my ($treatment, $control) = $fn =~ /_([^_]+)vs([^_]+)_ebayes$/; 
        my $rank = 1;
        while (<IN>) {
                chomp;
                next if $_ =~ /^"ID",/;
                $_ =~ s/"//g;
                my ($name_full, $logFC, $AveExpr, $t, $p_value, $adj_p_value, $b) = split(',');
                my $name = substr($name_full,0,8);
                my $check_genename_sql = "SELECT name FROM Genes WHERE name = '$name'";
                if (recover_info($check_genename_sql)) {
                        my $check_exprentry_sql = "SELECT COUNT(*) FROM  Ebayes_Results WHERE name = '$name' AND treatment = '$treatment' AND control = '$control'";
                        unless (check_existence($check_exprentry_sql)) {
                                $insert_ebayesexpr_sth->execute($name, $logFC, $AveExpr, $t, $p_value, $adj_p_value, $b, $name_full, $treatment, $control, $rank);
                        }
                }
                $rank += 1;
        }
        close IN;
        print "Done!\n";
        1;
}


sub InsertSQL_question_str {
	my $ques_num = shift;
	my $questions = '?'x$ques_num;
	return join(', ',split('', $questions));
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

sub sth { #prepare,execute a sql
        my $insert_sth = $dbh->prepare($_[0]);#prepare , look for basic sql errors
        $insert_sth->execute(); # action takes part here!!
        return 1;
}
