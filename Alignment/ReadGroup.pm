#!/usr/bin/perl

package ReadGroup;
use strict;
use warnings;
use lib '/home/hpriest/Scripts/Library/';
use Alignment::Read;
use hdpTools;
use FindBin;
use lib "$FindBin::Bin/";
our $FIDI;
our @FIDS;

sub new {
	my $class=shift;
	my $self = {
		Name		=>shift,
		Reads		=>[],
		ReadNumber	=>0,
	};
	die "Cannot call ". (caller(0))[3] ." without passing a group Name!\n" unless $self->{Name};
	bless $self, $class;
	return $self;
}

sub addRead {
	my $self=shift;
	my $readID=shift;
	push @{$self->{Reads}},Read->new($readID);
	$self->{ReadNumber}+=1;
	return $self->{ReadNumber}-1;
}

sub getNumReads {
	my $self=shift;
	return $self->{ReadNumber};
}

sub addAlignmentByIndex {
	my $self=shift;
	my $index=shift;
	my $position=shift;
	my $chr=shift;
	my $str=shift;
	my $length=shift;
	$self->{Reads}[$index]->setAligned($chr,$position,$str,$length);
	return 1;
}

sub getAligned {
	my $self=shift;
	my @aligned;
	for(my$i=0;$i<$self->{ReadNumber};$i++){
		push @aligned, $i if ${$self->{Reads}}[$i]->isAligned();
	}
	return \@aligned;
}

sub getAlignmentByIndex {
	my $self=shift;
	my $index=shift;
	return ${$self->{Reads}}[$index]->getAlignment();
}


1;
