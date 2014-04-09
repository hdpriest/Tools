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
my $n;
$n=$DF1->loadFile($file1,$del);
my @names=@{$DF1->logNTransform()};
$DF1->printPerturbedFrame($outfile,$del,\@names);
exit;
