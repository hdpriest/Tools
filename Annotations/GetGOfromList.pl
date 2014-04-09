#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my ($ANN,$list)=@ARGV;

die "usage: perl $0 <file of go terms> <list of genes of interest>\n" unless $#ARGV==1;
my $output=$list;
$output=~s/\.txt//;
$output=~s/\.list//;
$output.=".AgriGoIn.txt";

my @ANN=@{hdpTools->LoadFile($ANN)};
my %ANN;
foreach my $line (@ANN){
	my ($ID,$GO)=split(/\t/,$line);
	if(defined($ANN{$ID})){
		push @{$ANN{$ID}}, $GO;
	}else{
		$ANN{$ID}=[];
		push @{$ANN{$ID}}, $GO;
	}
}

my @list=@{hdpTools->LoadFile($list)};
my %U;
my @output;
my $unAnnotated=0;
my $Total=0;
foreach my $gene (@list){
	$gene=~s/[FR]\_[as]t//;
	next if defined($U{$gene});
	$Total++;
	$U{$gene}=1;
	if(defined($ANN{$gene})){
		foreach my $GO (@{$ANN{$gene}}){
#			print $gene."\t".$GO."\n";
			push @output, $gene."\t".$GO;
		}
	}else{
#		warn $gene." not annotated\n";
		$unAnnotated++;
	}
}
warn "$unAnnotated of $Total unannotated\n";
hdpTools->printToFile($output,\@output);
