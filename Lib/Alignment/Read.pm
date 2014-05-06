#!/usr/bin/perl

package Read;
use strict;
use warnings;
use lib '/home/hpriest/Scripts/Library/';
use hdpTools;
use FindBin;
use lib "$FindBin::Bin/";


sub new {
	my $class=shift;
	my $self = {
		ReadID	=>shift,
		Sequence	=>undef,
		Quality	=>undef,
		AlignRef	=>undef,
		AlignPos	=>undef,
		AlignStr	=>undef,
		AlignLen	=>undef,
	};
	die "Cannot call ". (caller(0))[3] ." without passing a Read Identifier!\n" unless $self->{ReadID};
	bless $self, $class;
	return $self;
}

sub setAligned {
	my $self=shift;
	my ($Ref,$Pos,$Str,$Len)=@_;
	die "Cannot call ". (caller(0))[3] ." without passing a Reference ID\n" unless defined($Ref);
	die "Cannot call ". (caller(0))[3] ." without passing an Alignment Position\n" unless defined($Pos);
	die "Cannot call ". (caller(0))[3] ." without passing a Alignment String\n" unless defined($Str);
	$self->{AlignRef}=$Ref;
	$self->{AlignPos}=$Pos;
	$self->{AlignStr}=$Str;
	$self->{AlignLen}=$Len;
	return 1;
}

sub setSequence {
	my $self=shift;
	my $Seq=shift;
	die "Cannot call ". (caller(0))[3] ." without passing a sequence\n" unless defined($Seq);
	$self->{Sequence}=$Seq;
	return 1;
}

sub setQuality {
	my $self=shift;
	my $Qual=shift;
	die "Cannot call ". (caller(0))[3] ." without passing a Quality String\n" unless defined($Qual);
	$self->{Quality}=$Qual;
	return 1;
}

sub isAligned {
	my $self=shift;
	return 0 unless defined $self->{AlignRef};
	return 1;
}

sub getAlignment {
	my $self=shift;
	my $string=$self->{AlignRef}.",".$self->{AlignPos}.",".$self->{AlignStr}.",".$self->{AlignLen};
	return $string;
}


1;
