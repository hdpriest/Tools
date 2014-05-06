#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use hdpTools;

die "usage: perl $0 <input FPKM tracking file> <output table> <min FPKM to report> <Ordered list of BAM files\n\n" unless $#ARGV>=1;

my $inFile	= shift @ARGV;
my $outFile	= shift @ARGV;
my $minCutoff=shift @ARGV;
my @BAMs	= @ARGV;
for(my$x=0;$x<=$#BAMs;$x++){
	$BAMs[$x]=~s/_accepted_hits\.bam//;
}
my @output;
my @FPKM=@{hdpTools->LoadFile($inFile)};
my $header=shift @FPKM;
my %header=%{ParseHeader($header)};
foreach my $line (@FPKM){
	my @line=split(/\t/,$line);
	my $gene=$line[$header{"gene_id"}];
	my @Values;
	for(my$i=0;$i<=$#BAMs;$i++){
		my $b=$i+1;
		my $phrase="q".$b."_FPKM";
		my $BAM=$BAMs[$i]."_FPKM";
		my $value=0;
		if(defined($header{$phrase})){
			$value=$line[$header{$phrase}];
		}elsif(defined($header{$BAM})){
			$value=$line[$header{$BAM}];
		}else{
			die "more BAMs provided than is in the file!\n";
		}
		push @Values, $value;
	}
	my $output=$gene."\t".join("\t",@Values);
	push @output, $output if check(\@Values,$minCutoff);
}
my $H="gene_id\t".join("\t",@BAMs);
unshift @output, $H;

hdpTools->printToFile($outFile,\@output);

sub check {
	my @values=@{$_[0]};
	my $cutoff=$_[1];
	foreach my $value (@values){
		return 0 if $value<$cutoff;
	}
	return 1;
}


sub ParseHeader {
	my $header=shift;
	my @header=split(/\t/,$header);
	my %Hash;
	for(my$i=0;$i<=$#header;$i++){
		my $field=$header[$i];
		$Hash{$field}=$i;
	}
	return \%Hash;

}
