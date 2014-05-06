#!/usr/bin/perl

package AlignmentReference;
use FindBin;
use lib "$FindBin::Bin/../Library";
use strict;
use warnings;
our $FIDI;
our @FIDS;

sub new {
	my $class=shift;
	my $self = {
      	Name		=> shift,
		Sequence	=> shift,
		Coverage	=> [],
		Alignments	=> 0,
		CovByRead	=> [],
	};
	die "Can't call ". (caller(0))[3]." without passing a reference name\n\n" unless defined($self->{Name});
	die "Can't call ". (caller(0))[3]." without passing a reference length\n\n" unless defined($self->{Sequence});
      bless $self, $class;
	$self->{Length}=length($self->{Sequence});
	for(my$i=0;$i<=$self->{Length};$i++){
		${$self->{Coverage}}[$i]=0;
	}
	return $self;
}

sub addAlignment {
	my $self=shift;
	$self->{Alignments}++;
	return 1;
}

sub getReadList {
	my $self=shift;
	my $start=shift;
	my $stop=shift;
	my @reads;
	for(my$i=$start;$i<=$stop;$i++){
		if(defined(${$self->{CovByRead}}[$i])){
			push @reads, @{$self->{CovByRead}[$i]};
		}else{
		}
	}
	return \@reads;
}

sub addCovByRead {
	my $self=shift;
	my $start=shift;
	my $stop=shift;
	my $head=shift;
	if(defined($self->{CovByRead}[$start])){
		push @{$self->{CovByRead}[$start]}, $head;
	}else{
		$self->{CovByRead}[$start]=[];
		push @{$self->{CovByRead}[$start]}, $head;
	}
	if(defined($self->{CovByRead}[$stop])){
		push @{$self->{CovByRead}[$stop]}, $head;
	}else{
		$self->{CovByRead}[$stop]=[];
		push @{$self->{CovByRead}[$stop]}, $head;
	}
}

sub getNumAlignments {
	my $self=shift;
	return $self->{Alignments};
}


sub getSeq {
	my $self=shift;
	return $self->{Sequence};
}

sub getName {
	my $self=shift;
	return $self->{Name};
}

sub getNumTargets {
	my $self=shift;
	my $numTargets=scalar(keys(%{$self->{Name}}));
	return $numTargets;
}

sub sum_data {
    my $data = shift;
    my $sum = 0;
    foreach (@$data) {
        $sum += $_;
    }
    $sum;
}

sub median {
	my $self=shift;
	my $arrayref=shift;
	my @a=@$arrayref;
	die "Can't take the median of an empty array\n" if $#a==0;
	my @array=sort {$a <=> $b} @a;
	my $midInt=int(scalar(@array)/2);
	return $array[$midInt];
}

sub mean {
	my $self=shift;
	my $arrayref=shift;
	my $result;
	foreach (@$arrayref) { $result += $_ }
	return $result/scalar(@$arrayref);
}


sub stdev {
    my $arrayref = shift;
    my $mean = mean($arrayref);
    return sqrt(mean( [map(($_-$mean)**2,@$arrayref)]));
}

sub addCoverageByNT {
	my $self=shift;
	my $coord=shift;
	my $copy=shift;
	if(defined(${$self->{Coverage}}[$coord])){
		${$self->{Coverage}}[$coord]+=$copy;
	}else{
		die "$coord undefined on ".$self->{Name}."!\n";
	}
}

sub addCoverage {
	my $self=shift;
	my $start=shift;
	my $stop=shift;
	my $copy=shift;
	$self->{Alignments}++;
	for(my$i=$start;$i<=$stop;$i++){
		${$self->{Coverage}}[$i]+=$copy;
	}
	
}

sub parseWiggleToCoverage {
	my $self=shift;
	my $wigFile=shift;
	my $index=0;
	open(WIG,"<",$wigFile) || die "cannot open $wigFile!\n$!\nexiting...\n";
	until(eof(WIG)){
		my $line=<WIG>;
		chomp $line;
		if($line=~m/fixedStep\schrom=(.+)\sstart=1\sstep=1\sspan=1/){
			my $header=ucfirst($1);
			if($header ne $self->{Name}){
				die "Found $header in wiggle file, but this is being parsed as part of ". $self->{Name} ."\nexiting...\n\n";
			}
			$index=0;
		}else{
			${$self->{Coverage}}[$index]+=$line;
			$index++;
		}
	}
	close WIG;
}

sub getWiggleCoverage {
	my $self=shift;
	my @wiggle=@{$self->{Coverage}};
	return \@wiggle;
}

sub getNumOfBasesInRegionCovered {
	my $self=shift;
	my $start=shift;
	my $stop=shift;
	my $total=0;
	for(my$i=$start;$i<=$stop;$i++){
		$total++ if(${$self->{Coverage}}[$i]>0);
	}
	return $total;
}

sub getTotalCovOverContigRegion {
	my $self=shift;
	my $start=shift;
	my $stop=shift;
	my $total=0;
	my $length=$stop-$start;
	for(my$i=$start;$i<=$stop;$i++){
		$total+=${$self->{Coverage}}[$i];
	}
	return $total;
}

sub findPutativeExons {
	my $self=shift;
	my %exons;
	my $e=1;
	my $EF=0;
	my $start;
	my $stop;
	my @cov=@{$self->{Coverage}};
	for(my$i=0;$i<=$#cov;$i++){
		if($cov[$i]==0){
			if($EF==0){
			}else{
				$stop=$i-1;
				$exons{$e}=$start."-".$stop;
				$start=undef;
				$stop=undef;
				$EF=0;
				$e++;
			}
		}else{
			if($EF==0){
				$EF=1;
				$start=$i;
			}else{
			}
		}
	}
	return \%exons;
}
1;
