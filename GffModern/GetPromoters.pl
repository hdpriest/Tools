#!/usr/bin/perl
use warnings;
use strict;

use FindBin;
use lib "$FindBin::Bin/../Library";
use GFF3::GffTools;
use hdpTools;
my ($gff,$fastaRef,$phrase,$out,$length)=@ARGV;

die "usage: perl $0 <gff> <genome> <gene phrase> <output> <length>\n\n" unless $gff && $fastaRef && $phrase && $out && $length;

my %Fasta=%{hdpTools->LoadFasta($fastaRef)};
my $Gff=GffTools->new();
$Gff->loadGFF($gff);

my $OutFile=$out.".promoters.$length.fasta";

my @objects=@{$Gff->getObjectsOfType($phrase)};
my @output;
foreach my $gene (@objects){
	my $object=$Gff->getObjectByID($gene);
	my $contig=$object->getContig();
	if(defined($Fasta{$contig})){
	}else{
		die "Could not find fasta contig: $contig\n";
	}
	my $promoter=$object->getUpstreamRegion($Fasta{$contig},$length,0);
	push @output, ">".$gene;
	push @output, $promoter;
}

hdpTools->printToFile($OutFile,\@output);
