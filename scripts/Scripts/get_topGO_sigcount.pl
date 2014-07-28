#!/usr/bin/perl
use warnings;
use strict;



my @files = glob('*.tab');

&get_topGO($_) foreach (@files);

sub get_topGO {
	my $file = shift;
	my $R_script = "

library(ggplot2)

pdf(paste('$file','topGO_sigcount.pdf',sep=''))
topGO_df <- read.table('$file', sep='\t',header=TRUE)
topGO_df <- topGO_df[topGO_df\$weight <= 0.05,]
topGO_df\$Term <- reorder(topGO_df\$Term, -topGO_df\$weight)
ggplot(data=topGO_df, aes(x= Term, y=Significant)) + geom_bar(stat='identity') + coord_flip()
dev.off()
";

	open OUT, '>', $file.'topGO_sigcount.R' or die "Cannot write into file";
	print OUT $R_script;
	close OUT;

	!system('Rscript '.$file.'topGO_sigcount.R') or die "Cannot execute R script";
	return 1;
}