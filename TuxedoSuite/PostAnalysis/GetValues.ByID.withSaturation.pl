#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;

die "usage $0 <input file of gene_exp.diff> <output file> <List of IDs within file>\n\n" unless $#ARGV==2;

my $infile	=$ARGV[0];
my $outfile	=$ARGV[1];
my @List	=@{hdpTools->LoadFile($ARGV[2])};

my %D=%{parse_diffexp($infile)};

my @header=("Locus");

my @keys=keys %D;

my @Output;

for(my$i=0;$i<=$#List;$i++){
	my $key=$List[$i];
	push @header, $D{$key}{"s1n"} if $i==0;
	push @header, $D{$key}{"s2n"} if $i==0;
	next unless defined $D{$key};
	my @output;
	push @output, $key;
	push @output, $D{$key}{"s1v"};
	push @output, $D{$key}{"s2v"};
	push @Output, join("\t",@output);
}
unshift @Output, join("\t",@header);
hdpTools->printToFile($outfile,\@Output);

exit(0);

sub parse_diffexp {
	my $file=shift @_;
	my @infile=@{hdpTools->LoadFile($file)};
	my %D;
	my $header=shift @infile;
	foreach my $line (@infile){
		my @line	=split(/\t/,$line);
		my $gene	=$line[1];
		$gene=~s/\.g//;
		my $s1name	=$line[4];
		my $s2name	=$line[5];
		my $value1	=$line[7];
		my $value2	=$line[8];
		my $log2fc	=$line[9];
		my $sigans	=$line[$#line];
		if(defined($D{$gene})){
			die "$gene found twice!\n";
		}else{
			$D{$gene}={};
			$D{$gene}{"s1n"}=$s1name;
			$D{$gene}{"s2n"}=$s2name;
			$D{$gene}{"s1v"}=$value1;
			$D{$gene}{"s2v"}=$value2;
			$D{$gene}{"l2f"}=$log2fc;
			$D{$gene}{"sig"}=$sigans;
		}
	}
	return \%D;
}
