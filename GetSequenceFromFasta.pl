#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/Lib";
use Tools;

die "usage: perl $0 <fasta file> <id>\n\n";

my $file=$ARGV[0];
my $id=$ARGV[1];

my %Fasta=%{Tools->LoadFile($file)};
if(defined($Fasta{$id})){
	print ">".$id."\n".$Fasta{$id}."\n";
}else{
	warn "$id not found in file\n";
}
exit(0);
