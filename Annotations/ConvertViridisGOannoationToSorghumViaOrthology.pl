#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;

die "usage: perl $0 <file to convert IDs within> <Map of Ids to gene IDs>\n" unless $#ARGV==1;

my @file=@{hdpTools->LoadFile($ARGV[0])};
my @Map=@{hdpTools->LoadFile($ARGV[1])};

my %Map;
foreach my $line (@Map){
	my ($sorghum,$viridis,$type,$bit,$e)=split(/\t/,$line);
	$Map{$viridis}=$sorghum;
}
my %File;
my $missed=0;
foreach my $line (@file){
	my ($id,$go)=split(/\t/,$line);
	if(defined($Map{$id})){
		print $Map{$id}."\t".$go."\n";
	}else{
		$missed+=1;
	}
}

warn "$missed genes of ". scalar(@file) ." did not have a hit\n";
