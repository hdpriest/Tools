#!/usr/bin/perl
use warnings;
use strict;

use lib '/home/hpriest/Scripts/Library';
use hdpTools;

my ($genelist,$iprscan)=@ARGV;

die "usage: perl $0 <gene list> <ipr scan file>\n\n" unless $#ARGV==1;

my @List=@{hdpTools->LoadFile($genelist)};
my @IPS= @{hdpTools->LoadFile($iprscan)};
my %A;
foreach my $gene (@List){
	$A{$gene}={};
}
foreach my $line (@IPS){
	my @line=split(/\t/,$line);
	my $id=$line[0];
	$id=~s/\.\d+//;
	next if $id eq "Bradi1g02575";
	my $x=0;
	next unless defined ($line[13]);
	while($line[13]=~m/(GO:\d\d\d\d\d\d\d)/g){
		my $GO=$1;
		$x++;
		if(defined($A{$id})){
			$A{$id}{$GO}=1;
		}else{
			die "could not find $id in hash!\n";
		}
	}
}

foreach my $gene (keys %A){
	foreach my $GO (keys %{$A{$gene}}){
		print $gene."\t".$GO."\n";
	}
}
