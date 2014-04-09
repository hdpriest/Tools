#!/usr/bin/perl
use warnings;
use strict;

my $fq1=$ARGV[0];
my $fq2=$ARGV[1];
my $out=$ARGV[2];

die "usage: perl $0 <fastq 1> <fastq 2> <output fastq file>\n\n" unless $#ARGV==2;

my %data;

open(SEQ1,"<",$fq1) || die "cannot open $fq1!\n$!\nexiting...\n";
until(eof(SEQ1)){
	my $header=<SEQ1>;
	my $sequence=<SEQ1>;
	my $had2=<SEQ1>;
	my $qual=<SEQ1>;
	chomp $header;
	chomp $sequence;
	chomp $qual;
	$header=~s/\@//;
	$header=~s/\_\d$//;
	$header=~s/\s.+//;
	if($data{$header}){
		die "$header found already!\n";
	}else{
		$data{$header}={};
		$data{$header}{'1s'}=$sequence;
		$data{$header}{'1q'}=$qual;
	}
}
close SEQ1;

open(SEQ2,"<",$fq2) || die "cannot open $fq2!\n$!\nexiting...\n";
until(eof(SEQ2)){
	my $header=<SEQ2>;
	my $sequence=<SEQ2>;
	my $had2=<SEQ2>;
	my $qual=<SEQ2>;
	chomp $header;
	chomp $sequence;
	chomp $qual;
	$header=~s/\@//;
	$header=~s/\_\d$//;
	$header=~s/\s.+//;
	if($data{$header}){
		$data{$header}{'2s'}=$sequence;
		$data{$header}{'2q'}=$qual;
	}else{
		die "$header does not have a pair\n";
	}
}
close SEQ2;

open(OUT,">",$out) || die "Cannot open $out!\n$!\nexiting...\n";
foreach my $key (keys %data){
	print OUT "\@".$key."\\1\n".$data{$key}{'1s'}."\n";
	print OUT "\+".$key."\\1\n".$data{$key}{'1q'}."\n";
	print OUT "\@".$key."\\2\n".$data{$key}{'2s'}."\n";
	print OUT "\+".$key."\\2\n".$data{$key}{'2q'}."\n";
}
close OUT;

