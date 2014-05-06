#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use Alignment::BWAtools;

my $config=$ARGV[0];


die "usage: perl ./runBWA.pl <BWAtools configuration file>\n\nThis script runs BWA, converts the output, and aggregates the results into alignment ratios and counts\n\n" unless $#ARGV==0;

die "$config does not exist!\n" unless -e $config;

my $BWA=BWAtools->new($config);
$BWA->getFilesFromConfig();
my $alg=$BWA->runBWAalign();
$BWA->runBWAconvert() unless $alg eq "bwasw";
my $ref=$BWA->getAlignmentStats();
my %hash=%$ref;
foreach my $file (sort {$a cmp $b} keys %hash){
	my $Num=$hash{$file};
	my $outLine=$file;
	$outLine=~s/.+\///g;
	$outLine.=",$Num";
	print $outLine."\n";
}
$BWA->findUnaligned();
$BWA->cleanUp();
exit(0);
