#!/usr/bin/perl

package GffObject;
use FindBin;
use lib "$FindBin::Bin/../Library";
use strict;
use warnings;

sub new {
        my $class=shift;
        my $self = {
            ID		=> shift,
		IAM		=> shift,
		contig	=> shift,
            start		=> undef,
		stop		=> undef,
		strand	=> undef,
		primaryGO	=> {},
		homology	=> {},
		children	=> [],
		function	=> {},
		parents	=> [],
        };
        bless $self, $class;
        die "Cannot call ". (caller(0))[3] ." without passing an ID.\n\n" unless $self->{"ID"};
        die "Cannot call ". (caller(0))[3] ." without passing a type.\n\n" unless $self->{"IAM"};
        die "Cannot call ". (caller(0))[3] ." without passing a contig.\n\n" unless $self->{"contig"};
        return $self;
}

sub getChildren {
	my $self=shift;
	return $self->{children};
}
sub getParents {
	my $self=shift;
	return $self->{parents};
}

sub getAllInformation {
	my $self=shift;
	return $self;
}

sub addGO {
	my $self=shift;
	my $go=shift;
	$self->{primaryGO}{$go}=1;
	return 1;
}

sub getGOAnnotation {
	my $self=shift;
	my @GO=keys %{$self->{primaryGO}};
	return \@GO;
}

sub getFunctionalAnnotation {
	my $self=shift;
	return $self->{function};
}

sub addFunctionalAnnotation {
	my $self=shift;
	my $phrase=shift;
	my $stanza=shift;
	die "Need all of : annotation phrase (arg0) and stanza (arg1) to define functional annotation\n" unless((defined $phrase)&&(defined $stanza));
	$self->{function}{$phrase}=$stanza;
	return 1;
}

sub getHomologyByPhrase {
	my $self=shift;
	my $phrase=shift;
	if(defined($self->{homology}{$phrase})){
		return $self->{homology}{$phrase};
	}else{
		return undef;
	}
	return undef;
}

sub addHomology {
	my $self=shift;
	my $phrase=shift;
	my $target=shift;
	my $relation=shift;
	die "Need all of : species phrase (arg0), target ID (arg1), and relationship type (arg2) to defined a homology relation\n" unless((defined $phrase)&&(defined $target)&&(defined $relation));
	my %h;
	$h{target}=$target;
	$h{relation}=$relation;
	$self->{homology}{$phrase}=\%h;
	return 1;
}

sub getContig {
	my $self=shift;
	return $self->{contig};
	
}

sub getUpstreamRegion {
	## enables users to get upstream region... user responsible for not using it on things like exons...
	my $self=shift;
	my $reference=shift;
	my $length=shift;
	my $overlap=shift;   #### overlap into gene locus
	my @coordinates=sort {$a <=> $b} ($self->{start},$self->{stop});
	my $promoter;
	if($self->{strand} eq "-"){
		my $start=$coordinates[1]-$overlap+1;
		$start=0 if $start<0;
		if(($start+$length)>length($reference)){
			$length=length($reference)-$start-1;
		}
		$promoter=substr($reference,$start,$length);
		$promoter=_reverse_transcribe($promoter);
	}else{
		my $start=$coordinates[0]-$length-1+$overlap;
		$promoter=substr($reference,$start,$length);
	}
	die "no promoter with $length and ...\n" unless defined $promoter;
	return $promoter;
	
}

sub getType {
	my  $self=shift;
	return $self->{IAM};
}

sub addChild {
	my $self=shift;
	die "No child sent!\n" unless defined ($_[0]);
	push @{$self->{children}}, $_[0];
	return $_[0];
}

sub addParent {
	my $self=shift;
	die "No parent sent!\n" unless defined ($_[0]);
	push @{$self->{parents}}, $_[0];
	return $_[0];
}

sub addStrand {
	my $self=shift;
	die "No strand sent!\n" unless defined ($_[0]);
	$self->{strand}=$_[0];
	return $_[0];
}

sub addStop {
	my $self=shift;
	die "No stop coordinate sent!\n" unless defined ($_[0]);
	$self->{stop}=$_[0];
	return $_[0];
}

sub addStart {
	my $self=shift;
	die "No start coordinate sent!\n" unless defined ($_[0]);
	$self->{start}=$_[0];
	return $_[0];
}

sub addName {
	my $self=shift;
	die "No name sent!\n" unless defined ($_[0]);
	$self->{name}=$_[0];
	return $_[0];
}

sub _reverse_transcribe {
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return $seq;
}
return 1;
