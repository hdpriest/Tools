#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;

my $file=$ARGV[0];
die "usage: perl $0 <file in>\n" unless $#ARGV==0;

my @file=@{hdpTools->LoadFile($file)};
foreach my $line (@file){
	my ($entry,$gene,$transcript,$protein,$PF,$PTHR,$KOG,$mapmap,$K,$GO,$ATN,$ATLN,$ATDESC,$OSID,$OSN,$OSDEC)=split(/\t/,$line);
	unless($OSID eq ""){
#		print $gene."\nGO:".$GO."\nKOG:$KOG\nATN:$ATN\nPTHR:$PTHR\nPF:$PF\n";
		print $gene."\t".$OSID."\n";
	}
}
