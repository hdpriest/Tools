#!/usr/bin/perl
use warnings;
use strict;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin/../Lib";
use Configuration;
use Alignment::Alignment;
use hdpTools;
package TopHatTools;




sub new {
	my $class=shift;
	my $self = {
      	configFile	=> shift,
      	config	=> undef,
		files		=> {},
      };
      bless $self, $class;
	$self->{config} = Configuration->new($self->{configFile});
      die "Cannot call ". (caller(0))[3] ." without passing a configuration file.\n\n" unless $self->{configFile};
	my $logFile=$self->{config}->get('PATHS','RootDir')."/TopHatLog$$";
	open(LOG,">",$logFile) || die "cannot open $logFile!\n$!\nexiting...\n";
      return $self;
}
sub GetAlignmentStatsFromFile {
	my $File=shift;
	my @File=@{hdpTools->LoadFile($File)};
	$File[1]=~m/\s+Input\s+\:\s+(\d+)/i;
	my $IR=$1;
	$File[2]=~m/\s+Mapped\s+\:\s+(\d+)\s.+/i;
	my $MR=$1;
	$File[3]=~m/\s+of\sthese\:\s+(\d+)\s.+/i;
	my $RR=$1;
	my @a=($IR,$MR,$RR);
	return \@a;
}
sub collectOutput {
	my $self	=shift;
	my $Threads	=$self->{config}->get('OPTIONS','Threads');
	my $rootDir =$self->{config}->get('PATHS','RootDir');
	my $collect =$rootDir."/Collected_Output";
	if(-e $collect){
		die "$collect exists! cannot re-collect to an existing directory. You'll have to do this manually.\n\n";
	}else{
		mkdir $collect;
	}
	my %Align;
	foreach my $file (keys %{$self->{files}}){
		my $alias   =$self->{files}{$file}{'alias'};
		unless(defined($alias)){
			print LOG "Dying inside alignment. Alias value not defined for $file\n\n";
			die "$file did not have a corresponding alias value associated with it under 'ALIAS' in the configuration file!\n". (caller(0))[3] ."\n";
		}
		my $oDir =$rootDir."/".$alias."_temp";
		my $abam =$oDir."/accepted_hits.bam";
		my $sFile=$oDir."/align_summary.txt";
		my $cFile=$collect."/".$alias."_accepted_hits.bam";
		$Align{$alias}=GetAlignmentStatsFromFile($sFile);
		my $command="mv $abam $cFile";
		my $thr=threads->create(\&launchJob,$command);
		while(threads->list()>=$Threads){
			my @thr=threads->list();
			$thr[0]->join();
		}
	}
	while(threads->list()>=1){
		my @thr=threads->list();
		$thr[0]->join();
	}
	MakeSummary($collect,\%Align);
	return 1;
}

sub MakeSummary {
	my $Dir=shift;
	my %Stats=%{$_[0]};
	my @Output;
	push @Output, "Library,Input,Mapped,Repeat";
	foreach my $alias (sort {$a cmp $b} keys %Stats){
		push @Output, $alias.",".join(",",@{$Stats{$alias}});
	}
	my $Path=$Dir."/OverallSummary.csv";
	hdpTools->printToFile($Path,\@Output);
	return 1;
}

sub runTopHat {
	my $self	=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $bThreads=$self->{config}->get('OPTIONS','Threads');
	my $binary 	=$self->{config}->get('PATHS','TopHat');
	my $dataDir	=$self->{config}->get('PATHS','dataDir');
	my $rootDir	=$self->{config}->get('PATHS','RootDir');
	my $index	=$self->{config}->get('PATHS','Index');
	foreach my $file (keys %{$self->{files}}){
		my $mMatch	=$self->{files}{$file}{'mM'};
		my $alias	=$self->{files}{$file}{'alias'};
		unless(defined($mMatch)){
			print LOG "Dying inside alignment. Mis-Match value not defined for $file\n\n";
			die "$file did not have a corresponding misMatch value associated with it under 'FILES' in the configuration file!\n". (caller(0))[3] ."\n";
		}
		unless(defined($alias)){
			print LOG "Dying inside alignment. Alias value not defined for $file\n\n";
			die "$file did not have a corresponding alias value associated with it under 'ALIAS' in the configuration file!\n". (caller(0))[3] ."\n";
		}


		my $dFile=$dataDir."/".$file;
		my $oDir =$rootDir."/".$alias."_temp";
# tophat2 -o ./output_temp --no-novel-juncs --segment-length 17 -p 16 /home/hpriest/TBrutnell/References/Indexes/seteria_viridis_db sFR_R1_Sv_A10_L002_R1.final.fastq
		my $command=$binary." -o $oDir --no-novel-juncs -p $bThreads $index $dFile";
		print LOG "$command\n";
		warn $command."\n";
		my $thr=threads->create(\&launchJob,$command);
		while(threads->list()>=$mThreads){
			my @thr=threads->list();
			$thr[0]->join();
		}
	}
	print LOG "Jobs launched.... waiting....\n";
	while(threads->list()>=1){
		my @thr=threads->list();
		$thr[0]->join();
	}
	print LOG "Alignment completed\n";
	return 1;
}

sub getFilesFromConfig {
	my $self=shift;
	foreach my $file ($self->{config}->getAll('FILES')){
		my $misMatches	=$self->{config}->get('FILES',$file);
		my $alias		=$self->{config}->get('ALIAS',$file);
		die "$file did not have a corresponding misMatch value associated with it under 'FILES' in the configuration file!\n". (caller(0))[3] ."\n" unless defined $misMatches;
		die "$file did not have a corresponding misMatch value associated with it under 'ALIAS' in the configuration file!\n". (caller(0))[3] ."\n" unless defined $alias;
		$self->{files}{$file}={};
		$self->{files}{$file}{'mM'}=$misMatches;
		$self->{files}{$file}{'alias'}=$alias;
	}
}

sub launchJob {
	my $command=shift;
	`$command`;
	return 1;
}

sub cleanUp {
	my $self=shift;
	print LOG "Cleanup routine called.\n";
	print LOG "All Jobs complete. Closing.\n";
	close LOG;
}

1;

