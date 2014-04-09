#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;

my $gff=$ARGV[0];

die "usage: perl $0 <gff> \n\n" unless $#ARGV==0;

my @GFF=@{hdpTools->LoadFile($gff)};
my %G;
foreach my $line (@GFF){
	next if $line=~m/^\#/;
	my ($chr,$src,$type,$start,$stop,$dot1,$strand,$dot2,$comment)=split(/\t/,$line);
	next unless $type eq "CDS";
	my @cmt=split(/\;/,$comment);
	my $unqid=$cmt[$#cmt];
	if(defined($G{$unqid})){
		push @{$G{$unqid}{'c'}}, $start;
		push @{$G{$unqid}{'c'}}, $stop;
	}else{
		$G{$unqid}={};
		$G{$unqid}{'c'}=[];
		$G{$unqid}{'chr'}=$chr;
		push @{$G{$unqid}{'c'}}, $start;
		push @{$G{$unqid}{'c'}}, $stop;
	}
}

foreach my $gene (keys %G){
	my $c=$G{$gene}{'c'};
	my $chr=$G{$gene}{'chr'};
	my $i=0;
	while(my $r=getFeatureIntronCoordByNum($i,$c)){
		my $line=$chr."\tJGI\tIntron\t".$$r[0]."\t".$$r[1]."\t.\t.\t.\tID=Intron.$gene\n";
		print $line;
		$i++;
	}
	
}


sub getFeatureIntronCoordByNum {
	my $i=shift;
	my @c=sort {$a <=> $b} @{$_[0]};
	my $i1=(2*($i+1))-1;
	my $i2=(2*($i+1));
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1]+1;
		my $c2=$c[$i2]-1;
		my $r=[$c1,$c2]; ######## Propagate this change to all other subs!
		return $r;
	}else{
		return 0;
	}
}
