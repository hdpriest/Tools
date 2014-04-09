#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my ($geneList,$ANN,$col)=@ARGV;
die "usage perl $0 <gene list> <annotation file (tab del)> <column of tab delim. file to take as annotation>\n\n" unless $#ARGV==2;

my @List=@{hdpTools->LoadFile($geneList)};
my @ANN =@{hdpTools->LoadFile($ANN)};
my %ANN;
foreach my $line (@ANN){
	my @line=split(/\t/,$line);
	my $id=shift @line;
	next unless defined $line[$col];
	if(defined($ANN{$id})){
		my $L=length($line[$col]);
		$ANN{$id}=$line[$col] if $L>length($ANN{$id});
	}else{
		$ANN{$id}=$line[$col];
	}
}

foreach my $gene (@List){
	if(defined($ANN{$gene})){
		print $gene."\t".$ANN{$gene}."\n";
	}elsif(defined($ANN{$gene.".1"})){
		print $gene."\t".$ANN{$gene.".1"}."\n";
	}else{
		print $gene."\tnone\n";
	}
}
