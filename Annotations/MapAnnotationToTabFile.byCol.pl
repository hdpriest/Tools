#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use Annotation::Annotation;

my $config	=$ARGV[0];
my $tabfile	=$ARGV[1];
my $col	=$ARGV[2];
die "Usage: perl $0 <config file> <tab del file> <column (0 index) containing ID to map with>\n\n" unless $#ARGV==2;

my $Annotation=Annotation->new($config);
warn "Loading GFF\n";
$Annotation->loadAnnotation();
$Annotation->IndexByName();
warn "Loading Genome\n";
$Annotation->loadGenome();
warn "Loading Annotations\n";
$Annotation->addDescriptions();
warn "Loading GO\n";
$Annotation->addGOannotations();

my @file=@{hdpTools->LoadFile($tabfile)};
my @output;
my @head;
my $c=0;
my @sources=@{$Annotation->getFunctionalKeys()};
foreach my $line (@file){
	my @line;
	if($line=~m/,/){
		@line=split(/,/,$line);
	}elsif($line=~m/\t/){
		@line=split(/\t/,$line);
	}else{
		warn "Could not determine delimiter:\n$line\n";
		next;
	}
	if(($line=~m/Locus/) || ($line =~ m/gene_id/)){
		@head=@line;
		next;
	}
	my $ID=$line[$col];
	my %H=%{$Annotation->getAnnotationOfEntireTree($ID)};
	my @go=@{$H{'go'}};
	push @line, join(";",@go);
	push @head, "GO" if $c==0;;
	foreach my $source (@sources){
		push @head, $source if $c==0;
		my @F=keys %{$H{'func'}{$source}};
		if(defined($F[0])){
			push @line, join(";",@F);
		}else{
			push @line, "N/A";
		}
	}
	push @output, join("\t",@line);
	$c++;
}
unshift @output, join("\t",@head);
map {print $_."\n"} @output;
