#!/usr/bin/perl

package AlignmentFeature;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../Library";
our $FIDI;
our @FIDS;

sub new {
	my $class=shift;
	my $self = {
      	target	=> shift,
		totalAlignments=>0,
		hitlocations=>{},
	};
      bless $self, $class;
      die "Cannot call ". (caller(0))[3] ." without passing a target name.\n\n" unless $self->{target};
	return $self;
}

sub addAlignment {
	my ($self,$start,$stop,$strand,$copy)=@_;
	if($copy=~m/Repeat\|/){
		$copy=~s/Repeat\|//; #### this just ignores repeats for now
	}
	my $ident="$start-$stop";
	if($self->{hitlocations}{$ident}){
		$self->{hitlocations}{$ident}{copy}+=$copy;
	}else{
		$self->{hitlocations}{$ident}={};
		$self->{hitlocations}{$ident}{copy}=0;
		$self->{hitlocations}{$ident}{copy}+=$copy;
	}
}

sub getNumAlignmentsAboveXcopy {
	my $self=shift;
	my $min=shift;
	my $numAlign=0;
	foreach my $hit (keys %{$self->{hitlocations}}){
		$numAlign++ if ($self->{hitlocations}{$hit}{copy}>$min);
	}
	return $numAlign;
}

sub getTotalReads {
	my $self=shift;
	my $numReads=0;
	foreach my $hit (keys %{$self->{hitlocations}}){
		$numReads+=$self->{hitlocations}{$hit}{copy};
	}
	return $numReads;
}

sub getTargName {
	my $self=shift;
	return $self->{target};
}

sub sum_data {
    my $data = shift;
    my $sum = 0;
    foreach (@$data) {
        $sum += $_;
    }
    $sum;
}


sub mean {
	my $self=shift;
	my ($arrayref)=shift;
	my $result;
	foreach (@$arrayref) { $result += $_ }
	return $result/scalar(@$arrayref);
}


sub stdev {
    my $arrayref = shift;
    my $mean = mean($arrayref);
    return sqrt(mean( [map(($_-$mean)**2,@$arrayref)]));
}
1;
