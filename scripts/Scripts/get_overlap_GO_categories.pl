#!/usr/bin/perl
use warnings;
use strict;

my $go_gene_imp_file = shift @ARGV;
open IN, '<', $go_gene_imp_file or die "Cannot open $go_gene_imp_file";

my %go_linecount;

while (<IN>) {
	chomp;
	next unless $_ =~ /GO:/;
	my ($GO_ID) = $_ =~ /\t(GO:\d+)\t/;
	if ($go_linecount{$GO_ID}) {
		$go_linecount{$GO_ID} += 1;
	}
	else {
		$go_linecount{$GO_ID} = 1;
	}
}

close IN;

open IN, '<', $go_gene_imp_file or die "Cannot open $go_gene_imp_file";
while (<IN>) {
	chomp;
	next unless $_ =~ /GO:/;
	my ($GO_ID) = $_ =~ /\t(GO:\d+)\t/;
	if ($go_linecount{$GO_ID} == 1) {
		print "$_\n";
	}
}
close IN;
