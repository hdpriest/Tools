#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use Annotation::Annotation;

my $config	 =$ARGV[0];
my $inputDir =$ARGV[1];
my $outputDir=$ARGV[2];
die "Usage: perl $0 <config file> <input directory (module*.txt)> <output directory>\n\n" unless $#ARGV==2;

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

opendir(DIR,$inputDir) || die "cannot open directory!!\n$!\nexiting...\n";
my @files=grep {m/module.+\.txt/i} readdir DIR;
closedir DIR;

foreach my $file (@files){
	my $root=$file;
	$root=~s/\.txt//;
	my $inPath =$inputDir."/".$file;
	my $outPath=$outputDir."/".$root.".annotated.txt";
	my @file=@{hdpTools->LoadFile($inPath)};
	my @output;
	foreach my $line (@file){
		if(($line=~m/\t/) || ($line=~m/,/)){
			die "line: $line has delimiters\n";
		}
		my %H=%{$Annotation->getAnnotationOfEntireTree($line)};
		my $newLine=$line;
		foreach my $key (sort {$a cmp $b} keys %{$H{'go'}}){
			$line.="\t".$key;
		}
#		foreach my $key (sort {$a cmp $b} keys %{$H{'func'}}){
#			$line.="\t".$key;
#		}
		push @output, $line;
	}
	hdpTools->printToFile($outPath,\@output);
}
