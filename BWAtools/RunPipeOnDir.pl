#!/usr/bin/perl

use warnings;
use strict;

my $dir=$ARGV[0];
my $tmp=$ARGV[1];
my $out=$ARGV[2];
my $ref=$ARGV[3];

die "usage: perl $0 <input directory of fastq files> <temp directory for sai files> <output directory of sorted bam files> <bwa ref>\n\n" unless $#ARGV==3;

opendir(DIR,$dir) || die "cannot open directory $dir!\n$!\nexiting...\n";
my @files=grep {m/fastq/} readdir(DIR);
closedir DIR;
my @cmds;
foreach my $file (@files){
	my $root=$file;
	$root=~s/\.fastq//;
	my $ipath=$dir."/".$file;
	my $opath=$out."/".$root.".bam";
	my $tpath=$tmp."/".$root.".sai";
#	bwa aln -t 4 ./hg19.fasta ./s1_1.fastq > ./s1_1.sai
#	bwa aln -t 4 ./hg19.fasta ./s1_2.fastq > ./s1_2.sai
#	bwa sampe ./hg19.fasta ./s1_1.sai ./s1_2.sai ./s1_1.fastq ./s1_2.fastq | \
#	samtools view -Shu - | \
#	samtools sort - - | \
#	samtools rmdup -s - - | \
#	tee s1_sorted_nodup.bam | \
#	bamToBed > s1_sorted_nodup.bed
	my $line ="bwa aln -t 16 $ref $ipath > $tpath";
	unshift @cmds,$line;
	my $pipe ="bwa samse $ref $tpath $ipath | samtools view -Sbhu - | samtools sort -o - - | samtools rmdup -s - - > $opath";
	push @cmds, $pipe;
}

print join("\n",@cmds)."\n";
