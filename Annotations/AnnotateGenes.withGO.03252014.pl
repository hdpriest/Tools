#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use Annotation::Annotation;

my $config	 =$ARGV[0];
my $inputFile =$ARGV[1];
my $outputFile=$ARGV[2];
die "Usage: perl $0 <config file> <input file (gene list)> <output file>\n\n" unless $#ARGV==2;

my $Annotation=Annotation->new($config);
warn "Loading GFF\n";
$Annotation->loadAnnotation();
$Annotation->IndexByName();
warn "Loading Genome\n";
$Annotation->loadGenome();
warn "Loading Annotations\n";
$Annotation->addDescriptions();
warn "Loading GO\n";
$Annotation->addGOannotations();

my @file=@{hdpTools->LoadFile($inputFile)};
my @output;
foreach my $line (@file){
	if(($line=~m/\t/) || ($line=~m/,/)){
		die "line: $line has delimiters\n";
	}
	my %H=%{$Annotation->getAnnotationOfEntireTree($line)};
	my $newLine=$line;
	my %uniq;
	foreach my $key (sort {$a cmp $b} keys %{$H{'go'}}){
		my  @keys=split(/\,/,$key);
		map{$uniq{$_}=1} @keys;
	}
	my @go=keys %uniq;
	$newLine.="\t".join(",",@go);
	push @output, $newLine if $newLine=~m/GO/;
}
hdpTools->printToFile($outputFile,\@output);
