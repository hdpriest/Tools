#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use DataFrame;

my ($file1,$outfile)=@ARGV;

my $usage="usage: perl $0
Arguments:
<data file 1>
<outfile>
";

die $usage unless $#ARGV==1;
my $del="\t";
my $DF1=new DataFrame;
my $n=$DF1->loadFile($file1,$del);
my @histogram=@{$DF1->getHistogramOfValues()};
my @output;
for(my$i=0;$i<=$#histogram;$i++){
	my $line;
	if(defined($histogram[$i])){
		$line=$i.",".$histogram[$i];
	}else{	
		$line=$i.",". 0;
	}
	push @output, $line;
}
hdpTools->printToFile($outfile,\@output);
exit;
