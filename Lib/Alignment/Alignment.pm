#!/usr/bin/perl

package Alignment;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/";
use Alignment::AlignmentFeature;
use Alignment::AlignmentReference;
use Alignment::Read;

sub new {
	my $class=shift;
	my $self = {
		References	=>{},
		ReadFiles	=>{},
		FIDI		=> undef,
		FIDS		=> [],
	};
	bless $self, $class;
	return $self;
}

sub addReadFile {
	my $self=shift;
	my $rFile=shift;
	open(READS,"<",$rFile) || die "cannot open $rFile!\n$!\nexiting...\n";
	my $firstline=<READS>;
	my $format;
	if($firstline=~m/^\@/){
		$format="fastq";
		warn "format of $rFile is fastq\n";
	}elsif($firstline=~m/^\>/){
		$format="fasta";
		warn "format of $rFile is fasta\n";
	}else{
		die "$rFile is neither a fasta or a fastq file - what is it ?\n";
	}
	seek(READS,0,0);
	$self->{ReadFiles}{$rFile}=ReadFile->new($rFile,$format);
	my $numAdded=0;
	until(eof(READS)){
		if($format eq "fastq"){
			my $h1=<READS>;
			my $s1=<READS>;
			my $h2=<READS>;
			my $s2=<READS>;
			chomp $h1;
			chomp $h2;
			chomp $s1;
			chomp $s2;
			$h1=~s/\@//;
			$h1=~s/\s.+//;
			my $res=$self->{ReadFiles}{$rFile}->addRead($h1,$s1,$h2,$s2);
			$numAdded++ if $res==1;
		}elsif($format eq "fasta"){
			my $h1=<READS>;
			my $s1=<READS>;
			chomp $h1;
			chomp $s1;
			$h1=~s/\>//;
			$h1=~s/\s.+//;
			my $res=$self->{ReadFiles}{$rFile}->addRead($h1,$s1,"NULL","NULL");
			$numAdded++ if $res==1;
		}else{
			die "unrecognized format\n";
		}
	}
	warn "Added $numAdded reads from file $rFile\n";
	close READS;
}

sub getUnalignedReads {
	my $self=shift;
	my $rfile=shift;
	my $ref;
	if(defined($self->{ReadFiles}{$rfile})){
		$ref=$self->{ReadFiles}{$rfile}->getUnalignedReads();
	}else{
		die "$rfile is an undefined value for a ReadFile object\n";
	}
	#### This returns a hashref where each key is a header, and 's' is sequence, 'q' is quality string
	return $ref;
}

sub getRatioOfAligningReads {
	my $self=shift;
	my $rfile=shift;
	my $aligned=$self->{ReadFiles}{$rfile}->getNumAlignedReads();
	my $total=$self->{ReadFiles}{$rfile}->getNumOfReads();
	my $ratio=$aligned/$total;
	return $ratio;
}

sub addFastaReferences {
	my $self=shift;
	my $refFile=shift;
	open(REF,"<",$refFile) || die "cannt open $refFile!\n$!\nexiting...\n";
	my %refs;
	my $head;
	until(eof(REF)){
		my $line=<REF>;
		chomp $line;
		if($line=~m/\>/){
			$head=$line;
			$head=~s/\>//;
			$head=ucfirst($head);
			$refs{$head}="";
		}else{
			$refs{$head}.=$line;
		}
	}
	foreach my $key (keys %refs) {
		$self->{Reference}{$key}=AlignmentReference->new($key,$refs{$key});
		push @{$self->{FIDS}}, $key;
	}
	$self->{FIDI}=0;

}

sub getNumTargets {
	my $self=shift;
	my $ref=shift;
	my $numTargs=$self->{Reference}{$ref}->getNumTargets;
	return $numTargs;
}

sub parseSamToWiggle {
	my $self=shift;
	my $samFile=shift;
	open(SAM,"<",$samFile) || die "Cannot open $samFile!\n$!\nexiting...\n";
	until(eof(SAM)){
		my $line=<SAM>;
		next if $line=~m/^\@/;
		chomp $line;
		my @line=split(/\t/,$line);
		my $len=length($line[9]);
		next if $line[2] eq "*";
		my $chrome=ucfirst($line[2]);
		my $start=$line[3];
		my $stop=$start+$len-1;
		die "start not defined: $line\n" unless $start;
		die "stop not defined: $line\n" unless $stop;
		die "space 2 on $line not defined\n" unless $chrome;
		if(defined($self->{Reference}{$chrome})){
			$self->{Reference}{$chrome}->addCoverage($start,$stop,1);
		}else{
			die "Never heard of $chrome before....\n";
		}
	}
	close SAM;
}

sub getAlignmentsAll {
	my $self=shift;
	my %align;
	foreach my $key (keys %{$self->{Reference}}){
		$align{$key}=$self->{Reference}{$key}->getNumAlignments();
	}
	return \%align;
}

sub getReadsByRegion {
	my $self=shift;
	my $chrome=shift;
	my $start=shift;
	my $stop=shift;
	if(defined($self->{Reference}{$chrome})){
		my @reads=@{$self->{Reference}{$chrome}->getReadList($start,$stop)};
		return \@reads;
	}else{
		die "Never heard of $chrome before....\n";
	}
}

sub parseSamToHitsComplex {
	my $self=shift;
	my $samFile=shift;
	open(SAM,"<",$samFile) || die "Cannot open $samFile!\n$!\nexiting...\n";
	until(eof(SAM)){
		my $line=<SAM>;
		next if $line=~m/^\@/;
		chomp $line;
		my @line=split(/\t/,$line);
		my $len=length($line[9]);
		next if $line[2] eq "*";
		my $chrome=ucfirst($line[2]);
		my $start=$line[3];
		my $stop=$start+$len-1;
		die "start not defined: $line\n" unless $start;
		die "stop not defined: $line\n" unless $stop;
		die "space 2 on $line not defined\n" unless $chrome;
		if(defined($self->{Reference}{$chrome})){
			$self->{Reference}{$chrome}->addCovByRead($start,$stop,$line[0]);
		}else{
			die "Never heard of $chrome before....\n";
		}
	}
	close SAM;

}

sub parseSamToHitsSimple {
	my $self=shift;
	my $samFile=shift;
	open(SAM,"<",$samFile) || die "Cannot open $samFile!\n$!\nexiting...\n";
	until(eof(SAM)){
		my $line=<SAM>;
		next if $line=~m/^\@/;
		chomp $line;
		my @line=split(/\t/,$line);
		my $len=length($line[9]);
		next if $line[2] eq "*";
		my $chrome=ucfirst($line[2]);
		my $start=$line[3];
		my $stop=$start+$len-1;
		die "start not defined: $line\n" unless(defined($start));
		die "stop not defined: $line\n" unless $stop;
		die "space 2 on $line not defined\n" unless $chrome;
		if(defined($self->{Reference}{$chrome})){
			$self->{Reference}{$chrome}->addAlignment();
		}else{
			die "Never heard of $chrome before....\n";
		}
	}
	close SAM;

}

sub countHitsInSam {
	my $self=shift;
	my $samFile=shift;
	my $hits=0;
	open(SAM,"<",$samFile) || die "Cannot open $samFile!\n$!\nexiting...\n";
	until(eof(SAM)){
		my $line=<SAM>;
		next if $line=~m/^\@/;
		chomp $line;
		my @line=split(/\t/,$line);
		next if $line[2] eq "*";
		$hits++;
	}
	close SAM;
	return $hits;
}

sub printWiggleToFile {
	my $self=shift;
	my $targetFile=shift;
	open(OUT,">",$targetFile) || die "cannot open $targetFile!\n$!\nexiting...\n";
	foreach my $ref (@{$self->{FIDS}}){
		my @coverage=@{$self->{Reference}{$ref}->getWiggleCoverage};
		my $initLine="fixedStep chrom=$ref start=1 step=1 span=1";
		print OUT $initLine."\n";
		print OUT join("\n",@coverage)."\n";
#		fixedStep chrom=scaffold_23 start=1 step=1 span=1
	}
	close OUT;
	return 1;
}

sub parseWiggleToWiggle {
	my $self=shift;
	my $wFile=shift;
	open(WIG,"<",$wFile) || die "cannot open $wFile!\n$!\nexiting...\n";
	my $chrome;
	my $index=0;
	until(eof(WIG)){
		my $line=<WIG>;
		chomp $line;
		if($line=~m/fixedStep\schrom=(.+)\sstart=1\sstep=1\sspan=1/){
			$chrome=ucfirst($1);
			$index=0;
		}else{
			if(defined($self->{Reference}{$chrome})){
				$self->{Reference}{$chrome}->addCoverageByNT($index,$line);
			}else{
				die "$wFile contains chromosome: $chrome which is not known from reference\n";
			}
			$index++;
		}
	}
	close WIG;
	
}

#sub printWiggleToFile {
#	my $self=shift;
#	my $outFile=shift;
	#	getWiggleCoverage
#	foreach my $chr (keys %{$self->{Reference}}){
#		my @wiggle=
#	}
#}

sub parseSamToAlignments {
	my $self=shift;
	my $samFile=shift;
	my $readFile=shift;
	open(SAM,"<",$samFile) || die "Cannot open $samFile!\n$!\nexiting...\n";
	until(eof(SAM)){
		my $line=<SAM>;
		next if $line=~m/^\@/;
		chomp $line;
		my @line=split(/\t/,$line);
		my $len=length($line[9]);
		next if $line[2] eq "*";
		my $chrome=ucfirst($line[2]);
		my $start=$line[3];
		my $stop=$start+$len-1;
		if(defined($readFile)){
			if(defined($self->{ReadFiles}{$readFile})){
				$self->{ReadFiles}{$readFile}->flagReadAligned($line[0]);
			}else{
				die "Readfiles doesn't know about $readFile\n";
			}
		}else{
			die "Cannot call ". (caller(0))[3]." without passing a source readFile (Arg 3)\n";
		}
	}
	close SAM;
}

sub parseSoapToWiggle{
	my $self=shift;
	my $File=shift;
	open(FILE,"<",$File) || die "Cannot open $File!\n$!\nexiting...\n";
	until(eof(FILE)){
		my $line=<FILE>;
		chomp $line;
		my @line=split(/\t/,$line);
		my $start=$line[3];
		my $len=length($line[1]);
		my $stop=$start+$len-1;
		my $chrome=ucfirst($line[2]);
		$self->{Reference}{$chrome}->addCoverage($start,$stop,1);
	}
	close FILE;
}

sub printUnalignedReadsFasta {
	my $self=shift;
	my $rfile=shift;
	my $outFile=shift;
	my $ref;
	if(defined($self->{ReadFiles}{$rfile})){
		my $result=$self->{ReadFiles}{$rfile}->printUnalignedReadsToFasta($outFile);
	}else{
		die "$rfile is an undefined value for a ReadFile object\n";
	}
	return 1;
	
}

sub parseHmToWiggle {
	my $self=shift;
	my $hmFile=shift;
	open(HM,"<",$hmFile) || die "cannot open $hmFile!\n$!\nexiting...\n";
	until(eof(HM)){
		my $line=<HM>;
		chomp $line;
		my @line=split(/\t/,$line);
		$self->{Reference}{$line[0]}->addCoverage($line[1],$line[2],$line[5]);
	}
	close HM;
}


sub getFeatureObjById {
	my $self=shift;
	my $id=shift;
	die "Can't call ". (caller(0))[3]." without data loaded into the HM Object\n\n" unless defined($self->{Reference});
	if(defined($self->{Reference}{$id})){
		return $self->{Reference}{$id};
	}else{
		return 0;
	}	
}

sub restartFeatureIter {
	my $self=shift;
	$self->{FIDI}=0;
	return 1;
}

sub getNextFeature {
	my $self=shift;
	die "Can't call ". (caller(0))[3]." without data loaded into the GffTools3 Object\n\n" unless defined($self->{Reference});
	if(defined(${$self->{FIDS}}[$self->{FIDI}])){
		my $nextId=${$self->{FIDS}}[$self->{FIDI}];
		$self->{FIDI}++;	
		return $self->{Reference}{$nextId}; 
	}else{
		return 0;
	}
}

sub getCurrentFeatureId {
	my $self=shift;
	die "Can't call ". (caller(0))[3]." without data loaded into the GffTools3 Object\n\n" unless $self->dataIsLoaded();
	my $thisIDI=$self->{FIDI}-1;
	if(defined(${$self->{FIDS}}[$thisIDI])){
		my $thisID=${$self->{FIDS}}[$thisIDI];
		return $thisID; 
	}else{
		return 0;
	}
		
}
1;
