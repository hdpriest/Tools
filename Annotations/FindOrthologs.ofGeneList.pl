#!/usr/bin/perl
use warnings;
use strict;

use hdpTools;
use Homology;

die "usage: perl $0 <map file> <list of source IDs>\n\n" unless $#ARGV==1;

my $Hm=Homology->new();
$Hm->loadFile($ARGV[0]);

my @IDs=@{hdpTools->LoadFile($ARGV[1])};
my $missed=0;
foreach my $id (@IDs){
	if($Hm->getTargetForSource($id)){
		my $target=$Hm->getTargetForSource($id);
		my $type  =$Hm->getTypeForSource($id);
		print $target."\n";
	}else{
		$missed++;
	}
}

warn "$missed IDs didn't have an entry (of ".scalar(@IDs)." IDs)\n\n";


