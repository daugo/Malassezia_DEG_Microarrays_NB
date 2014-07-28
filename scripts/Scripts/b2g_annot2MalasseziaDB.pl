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

my $Seq_anno_b2g_file = $ARGV[0];

&Fill_Blast2GO_Annotation($Seq_anno_b2g_file);


sub Fill_Blast2GO_Annotation{
	my $annotation_file = shift;

	#=========== Insert Queries ==================

	#Insert Blast2GO_Annotation query

	my $ques_sep_sam = &InsertSQL_question_str(5);
	my $insert_B2G_annot_sql = "INSERT INTO Blast2GO_Annotation(
								transcriptId,
								hit_desc,
							 	go_group,
							 	goAcc,
							 	term)
								VALUES ($ques_sep_sam)";
	my $insert_B2G_annot_sth = $dbh->prepare($insert_B2G_annot_sql);

	open IN, '<', $annotation_file or die "Cannot read $annotation_file";
	while (<IN>) {
		chomp;
		next if $_ =~ /^SeqName/;
		my ($seq_name, $hit_desc, $go_group, $goAcc, $term)= split(/\t/, $_);
		my ($transcriptId) = $seq_name =~ /^[^|]+\|[^|]+\|(\d+)\|/;
		
		if ($go_group eq 'F') {
			$go_group = 'molecular_function';
		}
		elsif ($go_group eq 'P') {
			$go_group = 'biological_process'; 
		}
		elsif ($go_group eq 'C') {
			$go_group = 'cellular_component';
		}

		my $check_transcriptId_sql = "SELECT transcriptId FROM Genes WHERE transcriptId = '$transcriptId'";
		if (recover_info($check_transcriptId_sql)) {
			my $check_Blast2GO_Annotation_entry_sql = "SELECT COUNT(*) FROM Blast2GO_Annotation WHERE
																				transcriptId = '$transcriptId' AND
																				goAcc = '$goAcc'";
			unless (check_existence($check_Blast2GO_Annotation_entry_sql)) {
				$insert_B2G_annot_sth->execute($transcriptId, $hit_desc, $go_group, $goAcc, $term);
			}
		}

	}
	close IN;

}

sub InsertSQL_question_str {
	my $ques_num = shift;
	my $questions = '?'x$ques_num;
	return join(', ',split('', $questions));
}


sub recover_info { #recover scalar info from the DB
        my ($recover_info) = $dbh->selectrow_array($_[0]);
        return $recover_info;
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
