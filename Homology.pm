#!/usr/bin/perl
use warnings;
use strict;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin/../Library";
use Configuration;
use hdpTools;
package Homology;




sub new {
	my $class=shift;
	my $self = {
		data		=> {},
		currentRow	=> undef,
      };
      bless $self, $class;
      return $self;
}

sub getTypeForSource {
	my $self=shift;
	my $source=shift;
	if(defined($self->{data}{$source})){
		return $self->{data}{$source}{'type'};
	}else{
		return undef;
	}
	die "arrived here is seemingly impossible...\n";
}

sub getTargetForSource {
	my $self=shift;
	my $source=shift;
	if(defined($self->{data}{$source})){
		return $self->{data}{$source}{'target'};
	}else{
		return undef;
	}
	die "arrived here is seemingly impossible...\n";
}

sub loadFile {
	my $self=shift;
	my $file=shift;
	my @File=@{hdpTools->LoadFile($file)};
	foreach my $line (@File) {
		my ($source,$target,$type,$bit,$eval)=split(/\t/,$line);
		my %h;
		$h{'target'}=$target;
		$h{'type'}=$type;
		$h{'bit'}=$bit;
		$h{'eval'}=$eval;
		die "Cannot load two relations for a source ID:\n$source\n$line\n" if defined ($self->{data}{$source});
		$self->{data}{$source}=\%h;
	}
	return scalar(keys(%{$self->{data}}));
}

sub _launchJob {
	my $command=shift;
	`$command`;
	return 1;
}

1;

