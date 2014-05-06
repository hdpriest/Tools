#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use Tuxedo::TopHatTools;

my $config=$ARGV[0];


die "usage: perl $0 <TopHatTools configuration file>\n\nThis script runs TopHat and collects the output for you.\n\n" unless $#ARGV==0;

die "$config does not exist!\n" unless -e $config;

my $TH=TopHatTools->new($config);
$TH->getFilesFromConfig();
my $alg=$TH->runTopHat();
$TH->collectOutput();
$TH->cleanUp();
exit(0);
