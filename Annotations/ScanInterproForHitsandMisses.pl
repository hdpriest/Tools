#!/usr/bin/perl
use warnings;
use strict;
use hdpTools;

die "usage: perl $0 <space-del list of annotation files>\n\n" unless $#ARGV>=0;

my %IDs;

foreach my $file (@ARGV){
	my @file=@{hdpTools->LoadFile($file)};
	foreach my $line (@file){
		my @annotations=split(/\t/,$line);
		my $id=shift @annotations;
		$id=~s/\.\d+//;
		my @Annotations;
		for(my$i=2;$i<=$#annotations;$i++){
			next if $annotations[$i] eq "Seq";
			next if $annotations[$i] eq "seq";
			next if $annotations[$i] eq "?";
			next if $annotations[$i] eq "NA";
			next if $i==5;
			next if $i==6;
			next if $annotations[$i]=~m/\d\d\-\w+\-\d+/;
			push @Annotations, $annotations[$i];
		}
		if(defined($IDs{$id})){
			foreach my $ann (@Annotations){
				$IDs{$id}{$ann}=1;
			}
		}else{
			$IDs{$id}={};
			foreach my $ann (@Annotations){
				$IDs{$id}{$ann}=1;
			}
		}
	}
}

## now, figure out which are annotated, and which are not
my $unAn="Genes.not.annotated.tab";
my $An="Genes.annotated.tab";
my @UnAn;
my @An;
foreach my $id (keys %IDs){
	my @Annotations=keys %{$IDs{$id}};
	my $annotated=0;
	my @passed;
	foreach my $ann (@Annotations){
		$ann=~s/^\s//;
		if($ann=~m/DUF\d+ domain containing protein/){
		}elsif($ann eq "NULL"){
		}elsif($ann eq "seg"){
		}elsif($ann =~m/AT\dG\d+/){
		}elsif($ann eq "Seg"){
		}elsif($ann =~ m/Protein of unknown function \(DUF\d+\)/){
		}elsif($ann =~ m/Plant protein of unknown function \(DUF\d+\)/){
		}elsif($ann eq "coiled-coil"){
		}elsif($ann eq "coil"){
		}elsif($ann eq "Coil"){
		}elsif($ann =~m/DUF\d+/){
		}else{
			push @passed, $ann;
			$annotated=1;
		}
	}
	if($annotated==1){
		my $line=$id."\t".join(",",@passed);
		push @An, $line;
	}else{
		my $line=$id."\t".join(",",@Annotations);
		push @UnAn, $line;
	}
}

hdpTools->printToFile($unAn,\@UnAn);
hdpTools->printToFile($An,\@An);
