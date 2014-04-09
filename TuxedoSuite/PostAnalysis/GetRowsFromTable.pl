#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use DataFrame;

my ($file1,$rows,$outfile)=@ARGV;

my $usage="usage: perl $0
Arguments:
<data file 1>
<files of row IDs to get>
<outfile>
";

die $usage unless $#ARGV==2;
my $del="\t";
my $DF1=new DataFrame;
my @Rows=@{hdpTools->LoadFile($rows)};
my $n;
$n=$DF1->loadFile($file1,$del);
my @D1=@{$DF1->getDataByRowID(\@Rows)};
my @H1=@{$DF1->getHeader()};

my @Output;
push @Output, "rowID".$del.join($del,@H1);
for(my$i=0;$i<=$#Rows;$i++){
	my $line=$Rows[$i].$del.join($del,@{$D1[$i]});
	push @Output, $line;
}
hdpTools->printToFile($outfile,\@Output);
exit;
