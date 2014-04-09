#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use DataFrame;

my ($file1,$file2,$outfile,$minValue,$ratioAbove,$minFC)=@ARGV;

my $usage="usage: perl $0
Arguments:
<data file 1>
<data file 2>
<output file>
<minumum value>
<ratio of data points above min value for inclusion>
<min FC cutoff for inclusion>
";

die $usage unless $#ARGV==5;

my $DF1=new DataFrame;
my $DF2=new DataFrame;
my $del="\t";


my $n;
$n=$DF1->loadFile($file1,$del);
my @R1=@{$DF1->getIDsByMinValueRatio($minValue,$ratioAbove)};
$n=$DF2->loadFile($file2,$del);
my @R2=@{$DF2->getIDsByMinValueRatio($minValue,$ratioAbove)};


my @Rows=@{GetCommonRows(\@R1,\@R2)};
my @D1=@{$DF1->getDataByRowID(\@Rows)};
my @D2=@{$DF2->getDataByRowID(\@Rows)};
my @H1=@{$DF1->getHeader()};
my @H2=@{$DF2->getHeader()};

my @Output;
push @Output, "rowID".$del.join($del,@H1).$del.join($del,@H2);
for(my$i=0;$i<=$#Rows;$i++){
	my $line=$Rows[$i].$del.join($del,@{$D1[$i]}).$del.join($del,@{$D2[$i]});
	push @Output, $line;
}
hdpTools->printToFile($outfile,\@Output);
exit;

sub GetCommonRows {
	my @array1=@{$_[0]};
	my @array2=@{$_[1]};
	my %H;
	map {$H{$_}=1} @array1;
	my @good;
	foreach my $ent (@array2){
		push @good, $ent if defined $H{$ent};
	}
	return \@good;
}


exit;
