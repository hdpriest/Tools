#!/usr/bin/perl

package Coverage;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/";
our $FIDI;
our @FIDS;

sub new {
	my $class=shift;
	my $self = {
		refID		=>shift,
		refLength	=>shift,
	};
	die "Cannot call ". (caller(0))[3] ." without passing a reference Identifier!\n" unless $self->{refID};
	die "Cannot call ". (caller(0))[3] ." without passing a reference length!\n" unless $self->{refLength};
	bless $self, $class;
	initObj($self);
	return $self;
}

sub initObj {
	my $self=shift;
	$self->{Reference}=[];
	for(my$i=0;$i<=$self->{refLength}-1;$i++){
		${$self->{Reference}}[$i]=0;
	}
}

sub getCov {
	my $self=shift;
	return $self->{Reference};
}

sub getWig {
	my $self=shift;
	my $output="fixedStep\tchrom=".$self->{refID}."\tstart=1\tstep=1\n";
	for(my$i=0;$i<=$self->{refLength}-1;$i++){
		$output.=${$self->{Reference}}[$i]."\n";
	}
	return $output;
}

sub addCoverage {
	my $self=shift;
	my $start=shift;
	my $stop=shift;
	my $value=shift;
	for(my$i=$start-1;$i<$stop;$i++){
		${$self->{Reference}}[$i]+=$value;
	}
	return 1;
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
		next unless $line[2] eq $self->{refID};
		my $start=$line[3];
		my $stop=$start+$len-1;
		for(my$i=$start;$i<=$stop;$i++){
			if(defined(${$self->{Reference}}[$i])){
				${$self->{Reference}}[$i]++;
			}else{
				warn "Coordinate outside of index: ". $self->{refID} ."\t$i\n";
			}
		}
	}
	close SAM;
}
1;
