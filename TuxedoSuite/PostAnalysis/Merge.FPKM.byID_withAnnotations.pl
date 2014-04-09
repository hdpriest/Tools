#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;

my ($tlist,$fpkm1,$fpkm2,$annf)=@ARGV;

die "usage: perl $0 <Tab-del list of paired genes (1\\t2)> <fpkm file 1> <fpkm file2> <annotation file>\n\n" unless $#ARGV==3;

my @Tlist=@{hdpTools->LoadFile($tlist)};
my @ANN  =@{hdpTools->LoadFile($annf)};
my %Ann;
map {my @a=split(/\t/,$_);$Ann{$a[0]}=$a[1];} @ANN;
my %FPKM1=%{LoadFPKM($fpkm1)};
my %FPKM2=%{LoadFPKM($fpkm2)};
my $Header="Locus1\tLocus2\tWC1\tsFR1\tWC2\tsFR2";
print $Header."\n";
my $undef=0;
foreach my $line (@Tlist){
	my ($Id1,$Id2)=split(/\t/,$line);
	if((defined($FPKM1{$Id1})) && (defined($FPKM2{$Id2}))){
		my @V1=@{$FPKM1{$Id1}};
		my @V2=@{$FPKM2{$Id2}};
		print $Id1."\t".$Id2."\t".join("\t",@V1)."\t".join("\t",@V2)."\t".$Ann{$Id1}."\n";
	}else{
#		warn "1: $Id1 and $Id2 not both defined!\n" unless defined $FPKM1{$Id1};
#		warn "2: $Id1 and $Id2 not both defined!\n" unless defined $FPKM2{$Id2};
		$undef++;
	}
}

warn $undef." pairs are not both defined\n";

sub LoadFPKM {
	my $file=shift @_;
	my @file=@{hdpTools->LoadFile($file)};
	my %D;
	my $head=shift @file;
	foreach my $line (@file){
		my ($ID,@d)=split(/\t/,$line);
		$D{$ID}=\@d;
	}
	return \%D;
}
