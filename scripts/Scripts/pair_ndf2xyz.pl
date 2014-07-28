#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use File::Basename;

#pair_ndf2xyz.pl
#=====Procedure===============
my $ndf_file = shift @ARGV;
my $pair_file = shift @ARGV;

&write_xys_file($ndf_file, $pair_file);
#=============================

sub write_xys_file {
	my $ndf_file = shift;
	my $pair_file = shift;

	open IN,'<',$pair_file or die "Cannot open pair file";
	my $header_line;
	while (<IN>) {
		chomp;
		last unless $_ =~ /^#/;
		$header_line = $_;
	}
 	close IN;

	my $ndf_info_ref = &get_ndf_info($ndf_file);
	my $pair_info_ref = &get_pair_info($pair_file);
	

	my ($fn_pair,$dir_pair,$suf_pair) = fileparse($pair_file,qr/\.[^.]*/);
	open OUT, '>',"$dir_pair/$fn_pair.xys" or die "Cannot write .xys file"; 
	print OUT $header_line."\n";
	print OUT "X	Y	SIGNAL	COUNT\n";

	foreach my $y (sort { $a <=> $b } keys %{$ndf_info_ref}) {
		foreach my $x (sort { $a <=> $b } keys %{ $ndf_info_ref->{$y} }) {
			my $xys_line = '';
			my $count = 1;
			if ($pair_info_ref->{$y}->{$x} && $ndf_info_ref->{$y}->{$x} =~ /experimental|random/i ) {
				my $signal = $pair_info_ref->{$y}->{$x};
				$xys_line = "$x	$y	$signal	$count\n";
			}
			elsif ( $pair_info_ref->{$y}->{$x} && $ndf_info_ref->{$y}->{$x} !~ /experimental|random/i ) {
				warn "coordinate reported as experimental according to pair file";
			} 
			elsif ( !$pair_info_ref->{$y}->{$x} && $ndf_info_ref->{$y}->{$x} =~ /experimental|random/i ) {
				warn "coordinate reported as experimental according to ndf file";
			}
			else {
				$xys_line = "$x	$y	NA	NA\n";
			}
			print OUT $xys_line;
		}
	}
	close OUT;
}	


#=====Process ndf file ============================
sub get_ndf_info {
	$ndf_file = shift @_;
	open IN,'<',$ndf_file or die "Cannot open $ndf_file";

	my %ndf_info;

	while (<IN>) {
		chomp;
		my @fields = split('\t');
		my ($treatment, $x, $y) = ($fields[11], $fields[15], $fields[16]);
		$treatment = 'random' unless $treatment;
		next if ($x =~ /^X$/ && $y =~ /^Y$/);
		$ndf_info{$y}{$x} = $treatment;  
		#my @coordinates = ($x, $y);
		#push @{ $ndf_info{ $treatment } }, \@coordinates ;           
	}

	close IN;
	return \%ndf_info
}

#=====Process pair files============================
sub get_pair_info {
	my $pair_file = shift;
	open IN,'<',$pair_file or die "Cannot open $pair_file";
	
	my %pair_info;
	while (<IN>) {
		chomp;
		next if $_ =~ /^(#|IMAGE_ID\t)/;
		my @fields = split('\t');
		my ($x, $y, $signal) = ($fields[5], $fields[6], $fields[9]);
		#my @coordinates = ($x,$y);
		$pair_info{$y}{$x} = $signal;
	}
	return \%pair_info;
	close IN;
}

=c
#====Display %ndf_info==================
sub display_ndf_info {
	my $ndf_info_ref  = shift;
	foreach ( keys %{$ndf_info_ref} ) {
		foreach (@{$ndf_info_ref->{$_}}) {
			print "@$_\n";
		}
	}
}
=cut

=c
sub display_xy_object {
	my $pair_info_ref = shift;
	foreach my $x (sort { $a <=> $b } keys %{$pair_info_ref}) {
		foreach my $y (sort { $a <=> $b } keys %{ $pair_info_ref->{$x} }) {
			my $signal = $pair_info_ref->{$x}->{$y};
			print "$x	$y	$signal\n";
		}
	}
}
=cut

	
		
