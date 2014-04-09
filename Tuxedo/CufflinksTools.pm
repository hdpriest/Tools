#!/usr/bin/perl
use warnings;
use strict;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin/../Library";
use Configuration;
use Alignment;
use hdpTools;
package CufflinksTools;




sub new {
	my $class=shift;
	my $self = {
      	configFile	=> shift,
      	config	=> undef,
		files		=> {},
      };
	die "THIS PACKAGE IS NOT FINISHED\n";
      bless $self, $class;
	$self->{config} = Configuration->new($self->{configFile});
      die "Cannot call ". (caller(0))[3] ." without passing a configuration file.\n\n" unless $self->{configFile};
	my $logFile=$self->{config}->get('PATHS','RootDir')."/TopHatLog$$";
	open(LOG,">",$logFile) || die "cannot open $logFile!\n$!\nexiting...\n";
      return $self;
}

sub runCuffLinks {
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

