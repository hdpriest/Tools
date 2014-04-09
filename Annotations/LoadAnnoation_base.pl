#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use Annotation::Annotation;

my $config=$ARGV[0];

die "Usage: perl $0 <config file>\n\n" unless defined $config;

my $Annotation=Annotation->new($config);
warn "Loading GFF\n";
$Annotation->loadAnnotation();
warn "Loading Genome\n";
$Annotation->loadGenome();
warn "Loading Homology\n";
$Annotation->addHomologies();
my @targets=@{$Annotation->getHomologyTargets()};

my $allKids=$Annotation->getAllParentsOfID("PAC:28400253.CDS.3");
foreach my $child (@$allKids){
	print $child."\n";
}
