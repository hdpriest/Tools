#!/usr/bin/perl

package GffFeature;
use lib '/deepspace/mockler/priesth/Scripts';
use strict;
use warnings;

sub new {
        my $class=shift;
        my $self = {
            name		=> shift,
		strand	=> shift,
		coordinates	=>[],
        };
        bless $self, $class;
        die "Cannot call ". (caller(0))[3] ." without passing a name.\n\n" unless $self->{name};
        die "Cannot call ". (caller(0))[3] ." without passing a strand.\n\n" unless $self->{strand};
        return $self;
}

sub getFeatureName {
	my $self=shift;
	return $self->{name};
}

sub addStanza {
	my $self=shift;
	my $c1=shift;
	my $c2=shift;
	push @{$self->{coordinates}}, $c1;
	push @{$self->{coordinates}}, $c2;
	return 1;
}

sub generateGeneGff {
	my $self=shift;
	my @c=@{$self->{coordinates}};
	my %conv=%{getGenomeToGeneConverter($self)};
	my $name=$self->{name};
	my $gene=$name;
	$gene=~s/\.\d+$//;
	my $output;
	@c=sort {$a <=> $b} @c;
	$output="$gene\tSubGene\tgene\t".$conv{$c[0]}."\t".$conv{$c[$#c]}."\t.\t+\t.\tID=$gene;Name=$gene\n";
	$output.="$gene\tSubGene\tmRNA\t".$conv{$c[0]}."\t".$conv{$c[$#c]}."\t.\t+\t.\tID=$name;Parent=$gene\n";
	for(my$E=0;$E<=$#c-1;$E+=2){
		$output.="$gene\tSubGene\texon\t".$conv{$c[$E]}."\t".$conv{$c[$E+1]}."\t.\t+\t.\tParent=$name\n";
	}
	return $output;
}

sub getGenomeToGeneConverter {
	my $self=shift;
	my @c=@{$self->{coordinates}};
	my %conv;
	if($self->{strand} eq "+"){
		@c=sort {$a <=> $b} @c;
		my $T=1;
		for(my$E=0;$E<=$#c-1;$E+=1){
			for(my$i=$c[$E];$i<=$c[$E+1];$i++){
				$conv{$i}=$T;
				$T++ unless $i==$c[$E+1];
			}
		}
	}else{
		@c=sort {$b <=> $a} @c;
		my $T=1;
		for(my$E=0;$E<=$#c-1;$E+=1){
			for(my$i=$c[$E];$i>=$c[$E+1];$i--){
				$conv{$i}=$T;
				$T++ unless $i==$c[$E+1];
			}
		}

	}
	return \%conv;
}

sub getGenomeToCdnaConverter {
	my $self=shift;
	my @c=@{$self->{coordinates}};
	my %conv;
	if($self->{strand} eq "+"){
		@c=sort {$a <=> $b} @c;
		my $T=0;
		warn "strand = +\n";
		for(my$E=0;$E<=$#c-1;$E+=2){
			warn $c[$E]."-".$c[$E+1]."\n";
			for(my$i=$c[$E];$i<=$c[$E+1];$i++){
				$conv{$i}=$T;
				$T++;
			}
		}
	}else{
		@c=sort {$b <=> $a} @c;
		my $T=0;
		warn "strand = -\n";
		for(my$E=0;$E<=$#c-1;$E+=2){
			warn $c[$E]."-".$c[$E+1]."\n";
			for(my$i=$c[$E];$i>=$c[$E+1];$i--){
				$conv{$i}=$T;
				$T++;
			}
		}

	}
	return \%conv;
}

sub getCdnaToGenomeConverter {
	my $self=shift;
	my @c=@{$self->{coordinates}};
	my @conv;
	if($self->{strand} eq "+"){
		@c=sort {$a <=> $b} @c;
		my $T=1;
		for(my$E=0;$E<=$#c-1;$E+=2){
			for(my$i=$c[$E];$i<=$c[$E+1];$i++){
				$conv[$T]=$i;
				$T++;
			}
		}
	}else{
		@c=sort {$b <=> $a} @c;
		my $T=1;
		for(my$E=0;$E<=$#c-1;$E++){
			for(my$i=$c[$E];$i>=$c[$E+1];$i--){
				$conv[$T]=$i;
				$T++;
			}
		}

	}
	return \@conv;
}

sub getCoordinates {
	my $self=shift;
	my @c=@{$self->{coordinates}};
	return \@c;
}

sub getFeatureWholeSequence {
	my $self=shift;
	my $refSeq=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $f="";
	for(my$j=0;$j<=$#c;$j+=2){
		my $k=$j+1;
		my $b=$c[$j]-1;
		my $l=$c[$k]-$c[$j]+1;
		my $s=substr($refSeq,$b,$l);
		$f.=$s;
	}
	if($self->{strand} eq "-"){
		$f=reverse_transcribe($f);
	}
	return $f;
}

sub getFeatureIntronSeqByNum {
	my $self=shift;
	my $refSeq=shift;
	my $i=shift;
	die "malformed input to ". (caller(0))[3] ." : dying\n" unless(defined($refSeq) && defined($i));
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $i1=(2*($i+1))-1;
	my $i2=(2*($i+1));
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1]+1;
		my $c2=$c[$i2]-1;
		my $l=$c2-$c1+1;
		my $b=$c1-1;
		my $s=substr($refSeq,$b,$l);
		if($self->{strand} eq "-"){
			$s=reverse_transcribe($s);
		}
		return $s;
	}else{
		return 0;
	}
}

sub getFeatureCdnaSeq {
	my $self=shift;
	my $refSeq=shift;
	my $i=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $cdna="";
	for(my$i=0;$i<=$#c-1;$i+=2){
		my $j=$i+1;
		if((defined($c[$i])) && (defined($c[$j]))){
			my $c1=$c[$i];
			my $c2=$c[$j];
			my $b=$c1-1;
			my $l=$c2-$c1+1;
			my $s=substr($refSeq,$b,$l);
			$cdna.=$s;
		}else{
			return undef;
		}
	}
	if($self->{strand} eq "-"){
		$cdna=reverse_transcribe($cdna);
	}
	return $cdna;
}

sub getFeatureExonSeqByNum {
	my $self=shift;
	my $refSeq=shift;
	my $i=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $i1=2*$i;
	my $i2=(2*$i)+1;
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1];
		my $c2=$c[$i2];
		my $b=$c1-1;
		my $l=$c2-$c1+1;
		my $s=substr($refSeq,$b,$l);
		if($self->{strand} eq "-"){
			$s=reverse_transcribe($s);
		}
		return $s;
	}else{
		return 0;
	}
}

sub getSjITDNSeqByNum {
	my $self=shift;
	my $refSeq=shift;
	my $i=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $i1=(2*($i+1))-1;
	my $i2=(2*($i+1));
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1];
		my $c2=$c[$i2];
		my $l=2;
		my $b1=$c1+1-1; ### add one to move from Eend to ITDN start, subtract one for 1->0 conversion
		my $b2=$c2-1-2; ### subtract 2 to go from Estart to ITDN start, subtract one for coord conversion
		my $s1=substr($refSeq,$b1,$l);
		my $s2=substr($refSeq,$b2,$l);
		my $s=$s1."-".$s2;
		if($self->{strand} eq "-"){
			$s=reverse_transcribe($s);
		}
		return $s;
	}else{
		return 0;
	}
}

sub getFeatureSJSeqByNum {
	my $self=shift;
	my $refSeq=shift;
	my $kmer=shift;
	my $i=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $i1=(2*($i+1))-1;
	my $i2=(2*($i+1));
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1];
		my $c2=$c[$i2];
		my $l=$kmer-1;
		my $b1=$c1-$kmer+2-1;
		my $b2=$c2-1;
		my $s1=substr($refSeq,$b1,$l);
		my $s2=substr($refSeq,$b2,$l);
		my $s=$s1.$s2;
		if($self->{strand} eq "-"){
			$s=reverse_transcribe($s);
		}
		return $s;
	}else{
		return 0;
	}
}

sub getFeatureExonCoordByNum {
	my $self=shift;
	my $i=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $i1=2*$i;
	my $i2=(2*$i)+1;
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1];
		my $c2=$c[$i2];
		my $r=[$c1,$c2]; ######## Propagate this change to all other subs!
		return $r;
	}else{
		return 0;
	}
}

sub getFeatureIntronCoordByNum {
	my $self=shift;
	my $i=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $i1=(2*($i+1))-1;
	my $i2=(2*($i+1));
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1]+1;
		my $c2=$c[$i2]-1;
		my $r=[$c1,$c2]; ######## Propagate this change to all other subs!
		return $r;
	}else{
		return 0;
	}
}

sub getFeatureSJCoordByNum {
	my $self=shift;
	my $i=shift;
	my @c=sort {$a <=> $b} @{$self->{coordinates}};
	my $i1=(2*($i+1))-1;
	my $i2=(2*($i+1));
	if((defined($c[$i1])) && (defined($c[$i2]))){
		my $c1=$c[$i1];
		my $c2=$c[$i2];
		my $r=[$c1,$c2]; ######## Propagate this change to all other subs!
		return $r;
	}else{
		return 0;
	}
}

sub reverse_transcribe {
        my $sequence=shift;
        $sequence=reverse($sequence);
        $sequence=~tr/ATCGatcg/TAGCtagc/;
        return $sequence;
}

1;

