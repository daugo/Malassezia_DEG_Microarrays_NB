#!/usr/bin/perl
use warnings;
use strict;
use Bio::SeqIO;
use File::Basename;

my $fasta_file = shift @ARGV;
my @ids_files = @ARGV;

foreach my $ids_file (@ids_files){
	&Parse_SAMcsv2completeSeqIDs($fasta_file,$ids_file);
}


sub Parse_SAMcsv2completeSeqIDs {
	my $fasta_file = shift;
	my $ids_file = shift;

	my $id_completeId_ref = &build_ID_completeID($fasta_file);

	open IN, '<', $ids_file or die "Cannot read $ids_file";

	my ($fn,$dir,$suf) = fileparse($ids_file,qr/\.[^.]*/);
	my $out_file  = $dir.$fn.'_seqIDs_complete'.'.txt';
	open OUT,'>', $out_file or die "Cannot write $out_file";
	
	while (<IN>) {
		chomp;
		next if ($_ !~ /,/);
		next if ($_ =~ /Row,|"ID",/);
		my @fields = split(',', $_);
		my $id = '';
		($id) = $fields[-1] =~ /(^MGL_[0-9]{4})/ if ($ids_file =~ /_sam\.csv$/);
		($id) = $fields[0] =~ /^"(MGL_[0-9]{4})/ if ($ids_file =~ /_ebayes\.csv$/);
		if (!$id) { 
			print "HERE: $_ $ids_file\n";
		}

		if ($id and ${$id_completeId_ref}{$id}) { 
			print OUT ${$id_completeId_ref}{$id}."\n";
		}
		else {
			warn "problem with $id in $ids_file";
			print OUT "jgi|Malgl1|2191|MGL_2189B\n";
		}
	}
	close OUT;
	close IN;
}

sub build_ID_completeID {
	my $fasta_file = shift;

	my $fasta_fh = SeqIO_fasta_fh($fasta_file,'<');

	my %id_completeId;
	while (my $transcript = $fasta_fh->next_seq) {
		my ($id) = $transcript->id =~ /([^|]+$)/;
		$id_completeId{$id} = $transcript->id;
	}
	return \%id_completeId;
}

sub SeqIO_fasta_fh {
		my $fasta_file = shift;
		my $io = shift;
		die "Not valid input/ouput sign for Bio::SeqIO->new" if $io !~ /^[><]$/ ;
		my $fh = Bio::SeqIO->new(
						-format => 'fasta',
						-file => $io.$fasta_file,
						);
		return $fh;
}
