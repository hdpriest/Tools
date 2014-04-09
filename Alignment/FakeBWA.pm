#!/usr/bin/perl
use warnings;
use strict;
use FindBin;
use lib "$FindBin::Bin/../Library";
package FakeBWA;


sub new {
        my $class=shift;
        my $self = {
		    numReplicates => shift,
		    noiseRatio	=> shift,
		    dexRatio	=> shift,
		    librarySize	=> shift,
		    key		=> {},
        };
        bless $self, $class;
        die "Cannot call ". (caller(0))[3] ." without passing a number of Replicates!\n" unless $self->{numReplicates};
        die "Cannot call ". (caller(0))[3] ." without passing a number of noise Ratio!\n" unless $self->{noiseRatio};
        die "Cannot call ". (caller(0))[3] ." without passing a number of ratio of differential expression!\n" unless $self->{dexRatio};
        die "Cannot call ". (caller(0))[3] ." without passing a Library Size!\n" unless $self->{librarySize};
	  die "Noise ratio must be between 0 and 1!\n\n" unless(($self->{noiseRatio}>0)&&($self->{noiseRatio}<1));
	  my $Success=InitObj($self);
        return $self;
}

sub InitObj {
	my $self=shift;
	$self->{Probabilities}={};
	$self->{SimulatedData}={};
	for(my $g=1;$g<=2;$g++){
		$self->{Probabilities}{$g}={};
		$self->{SimulatedData}{$g}={};
		for(my$n=1;$n<=$self->{numReplicates};$n++){
			$self->{Probabilities}{$g}{$n}={};
			$self->{SimulatedData}{$g}{$n}={};
			my $LibVar=rand(.2)-.1; #### rand var -.1 -- +.1
			my $size=$self->{librarySize}+int($self->{librarySize}*$LibVar); ### i.e., libraries of 30m reads vary in size between 27 and 33 mReads
			$self->{Probabilities}{$g}{$n}{'LS'}=$size;
			$self->{Probabilities}{$g}{$n}{'EXP'}=0;
		}
	}
	return 1;
}

sub makeData {
	my $self=shift;
	my $numTargets=shift;
	my $keyFile=shift;
	warn "Setting probabilities...\n";
	open(KEY,">",$keyFile) || die "cannot open $keyFile!\n$!\nexiting...\n";
	print KEY "Target\tDEX BOOLEAN\tFold-Change\tExpression level\tdirection\n";
	for(my$i=0;$i<=$numTargets;$i++){
		my $ref=MakeProbs($self,$i);
		my @info=@$ref;
		my ($targetName,$dex,$FC,$exp,$direction)=@info;
		if(defined($self->{key}{$targetName})){
			die "$targetName exists twice! oh crap!\n";
		}else{
			$self->{key}{$targetName}={};
			$self->{key}{$targetName}{'dex'}=1;
			$self->{key}{$targetName}{'FC'}=$FC;
			$self->{key}{$targetName}{'exp'}=$exp;
			$self->{key}{$targetName}{'direction'}=$direction;
			print KEY $targetName."\t".$dex."\t".$FC."\t".$exp."\t".$direction."\n";
		}
	}
	close KEY;
	warn "Done.\n";
	warn "Normalizing....\n";
	NormProbs($self);
	warn "Done.\n";
#	warn "Assigning reads...\n";
#	AssignReads($self);
#	warn "Done.\n";
	return 1;
}

sub AssignReads {
	my $self=shift;
	for(my$g=1;$g<=2;$g++){
		for(my$r=1;$r<=$self->{numReplicates};$r++){
			my $total=0;
			while($total<=$self->{Probabilities}{$g}{$r}{'LS'}){
				foreach my $target (keys %{$self->{Probabilities}{$g}{$r}}){
					next if (($target eq "EXP")||($target eq "LS"));
					my $R=rand(1);
					if($R<$self->{Probabilities}{$g}{$r}{$target}){
						$self->{SimulatedData}{$g}{$r}{$target}++;
						$total++;
					}else{
					}
				}
			}
		}
	}
	return 1;
}

sub NormProbs {
	my $self=shift;
	for(my$g=1;$g<=2;$g++){
		for(my$r=1;$r<=$self->{numReplicates};$r++){
			foreach my $target (keys %{$self->{Probabilities}{$g}{$r}}){
				next if (($target eq "EXP")||($target eq "LS"));
				$self->{Probabilities}{$g}{$r}{$target}=$self->{Probabilities}{$g}{$r}{$target}/$self->{Probabilities}{$g}{$r}{'EXP'};
				$self->{SimulatedData}{$g}{$r}{$target}=int($self->{Probabilities}{$g}{$r}{$target}*$self->{Probabilities}{$g}{$r}{'LS'});
			}
		}
	}
	return 1;
}

sub getProbs {
	my $self=shift;
	my $g=shift;
	my $r=shift;
	my @values;
	if(defined($self->{Probabilities}{$g}{$r})){
		foreach my $target (keys %{$self->{Probabilities}{$g}{$r}}){
			next if (($target eq "EXP")||($target eq "LS"));
			push @values, $self->{Probabilities}{$g}{$r}{$target};
		}
	}else{
		die "sample $g and replicate $r not found!\n";
	}
	return \@values;
}

sub MakeProbs {
	my $self=shift;
	my $num=shift;
	my $dex=0;
	my $targetName="Target$num";
	$dex=1 if(rand(1)<$self->{dexRatio});
	my $exp=int(rand(1000))+1; ### random expression level between 1 and 10000
	my $FC=undef;
	my $UD=int(rand(2));
	if($dex==1){
		if($UD==0){
			### DOWN
			$FC=rand(.49) + 0.01; ### FC between 0.001 and .5, or 2 and 1000 fold down reg
		}elsif($UD==1){
			#### UP
			$FC=int(rand(98) + 2); ### FC between 2 and 1000 fold up reg, integer
		}else{
			die "Randomization of direction failed: $UD\n";
		}
	}else{
		$FC=1;
	}

	##### The above generates the fold changes, directions, and boolean DEX values
	##### The below needs to set the expression of each target, normalized against
	##### The expression of the total set. This is so that 
	
	for(my$g=1;$g<=2;$g++){
		for(my$r=1;$r<=$self->{numReplicates};$r++){
			my $thisExp=$exp;;
			if($g==1){
				## WT
				$thisExp=$exp*1;
			}else{
				## MUT
				$thisExp=$exp*$FC;
			}
			my $noise=(rand(2*$self->{noiseRatio}) - $self->{noiseRatio});
			$thisExp=int($thisExp+($thisExp*$noise));
			$self->{Probabilities}{$g}{$r}{$targetName}=$thisExp;
			$self->{Probabilities}{$g}{$r}{'EXP'}+=$thisExp;
			$self->{SimulatedData}{$g}{$r}{$targetName}=0;
		}
	}
	#die "unfinished!!!\n";
	#### now have expression values, and total expression vlaues
	###3 need to normalize expression values to probability based on total expressoion
	##### and randomly assign reads until library size reached to each target

#				$self->{Probabilities}{$g}{$r}{$targetName}=$thisExp/$self->{Probabilities}{$g}{$r}{'LS'};
	return [$targetName,$dex,$FC,$exp,$UD];
}

sub printDataDirectory {
	my $self=shift;
	my $outDir=shift;
	for(my$g=1;$g<=2;$g++){
		for(my$r=1;$r<=$self->{numReplicates};$r++){
			if($g==1){
				my $outfile=$outDir."/WildType.fakeData.$r.sam";
				open(OUT,">",$outfile) || die "cannot open $outfile!\n$!\nexiting...\n";
				foreach my $target (keys %{$self->{SimulatedData}{$g}{$r}}){
					for(my$i=0;$i<=$self->{SimulatedData}{$g}{$r}{$target};$i++){
						my $output="";
						my $readId="Read.$r.$i";
						$output.=$readId;
						$output.="\t0";
						$output.="\t$target";
						$output.="\t1";
						$output.="\t1";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\n";
						############ OK! this is a line of SAM output - might need to be altered.
						print OUT $output;
					}
				}
				close OUT;
			}else{
				my $outfile=$outDir."/Mutant.fakeData.$r.sam";
				open(OUT,">",$outfile) || die "cannot open $outfile!\n$!\nexiting...\n";
				foreach my $target (keys %{$self->{SimulatedData}{$g}{$r}}){
					for(my$i=0;$i<=$self->{SimulatedData}{$g}{$r}{$target};$i++){
						my $output="";
						my $readId="Read.$r.$i";
						$output.=$readId;
						$output.="\t0";
						$output.="\t$target";
						$output.="\t1";
						$output.="\t1";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\tnull";
						$output.="\n";
						############ OK! this is a line of SAM output - might need to be altered.
						print OUT $output;
					}
				}
				close OUT;
			}
		}
	}

}

return 1;
