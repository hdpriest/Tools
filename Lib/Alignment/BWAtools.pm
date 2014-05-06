#!/usr/bin/perl
use warnings;
use strict;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin/../Library";
use Configuration;
use Alignment;
package BWAtools;

our %files;
my %return:shared;



sub new {
	my $class=shift;
	my $self = {
      	configFile	=> shift,
      	config	=> undef,
      };
      bless $self, $class;
	$self->{config} = Configuration->new($self->{configFile});
      die "Cannot call ". (caller(0))[3] ." without passing a configuration file.\n\n" unless $self->{configFile};

	my $logFile=$self->{config}->get('PATHS','LogFile').$$;
	open(LOG,">",$logFile) || die "cannot open $logFile!\n$!\nexiting...\n";
	warn "finishing...\n";
      return $self;
}

sub getFilesFromDir {
	my $self=shift;
	my $dir=$self->{config}->get('PATHS','dataDir');
	my $misMatches;
	if($self->{config}->get('OPTIONS','misMatches')){
		$misMatches=$self->{config}->get('OPTIONS','misMatches');
	}else{
		die "If you want to read input files from a data directory, a global mismatch value must exist (OPTIONS ==> misMatches)\n" . (caller(0))[3] ."\n";
	}
	opendir(DIR,$dir) || die "Cannot open $dir!\n$!\nexiting...\n";
	my @files=grep {m/fas/} readdir(DIR);
	closedir DIR;
	foreach my $file (@files){
		$files{$file}=$misMatches;
		print LOG "$misMatches permissable with $file\n";
	}
	return 1;
}

sub runBWAalign {
	my $self	=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $bThreads=$self->{config}->get('OPTIONS','bwaThreads');
	my $iDir	=$self->{config}->get('PATHS','dataDir');
	my $oDir	=$self->{config}->get('PATHS','outDir');
	my $index	=$self->{config}->get('PATHS','bwaIndex');
	my $alg	=$self->{config}->get('OPTIONS','algorithm');
	print LOG "Running BWA alignment....\n";
	foreach my $file (keys %files){
		my $iPath=$iDir."/$file";
		my $oPath=$oDir."/$file";
		my $misMatches=$files{$file};
		unless(defined($misMatches)){
			print LOG "Dying inside alignment. Mis-Match value not defined for $file\n\n";
			die "$file did not have a corresponding misMatch value associated with it under 'OPTIONS' in the configuration file!\n". (caller(0))[3] ."\n";
		}
		my $command;
		if($alg eq "bwasw"){
			$oPath=~s/fa.+$/sam/;
			$command="bwa bwasw -t $bThreads -f $oPath $index $iPath";
		}else{
			$oPath=~s/fa.+$/bin/;
			$command="bwa aln -o 0 -t $bThreads -R 5 -n $misMatches -f $oPath $index $iPath";
		}
		print LOG "$command\n";
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
	return $alg;
}

sub indexBam {
	my $self	=shift;
	my $subDir	=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $bThreads=$self->{config}->get('OPTIONS','bwaThreads');
	my $iDir	=$self->{config}->get('PATHS','dataDir');
	my $oDir	=$self->{config}->get('PATHS','outDir');
	my $samTools=$self->{config}->get('PATHS','samTools');
	die "Samtools executable not specified!\n" unless defined $samTools;
	my $outDir=$oDir;
	$outDir.="/".$subDir if defined($subDir);
	mkdir $outDir unless -e $outDir;
	foreach my $file (keys %files){
		my $cPath=$outDir."/$file";
		$cPath=~s/fas.+/sort/;
#		samtools index FC187_1.sort.bam
		my $command=$samTools." index $cPath";
		print LOG "$command\n";
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
	print LOG "Conversion completed\n";
}

sub sortBam {
	my $self	=shift;
	my $subDir	=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $bThreads=$self->{config}->get('OPTIONS','bwaThreads');
	my $iDir	=$self->{config}->get('PATHS','dataDir');
	my $oDir	=$self->{config}->get('PATHS','outDir');
	my $samTools=$self->{config}->get('PATHS','samTools');
	die "Samtools executable not specified!\n" unless defined $samTools;
	my $outDir=$oDir;
	$outDir.="/".$subDir if defined($subDir);
	mkdir $outDir unless -e $outDir;
	foreach my $file (keys %files){
		my $oPath=$outDir."/$file";
		my $cPath=$outDir."/$file";
		$oPath=~s/fas.+/bam/;
		$cPath=~s/fas.+/sort/;
		my $command=$samTools." sort $oPath $cPath";
		print LOG "$command\n";
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
	print LOG "Conversion completed\n";
}

sub convertSamToBam {
	my $self	=shift;
	my $subDir	=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $bThreads=$self->{config}->get('OPTIONS','bwaThreads');
	my $iDir	=$self->{config}->get('PATHS','dataDir');
	my $oDir	=$self->{config}->get('PATHS','outDir');
	my $samTools=$self->{config}->get('PATHS','samTools');
	die "Samtools executable not specified!\n" unless $samTools;
	my $outDir=$oDir;
	$outDir.="/".$subDir if defined($subDir);
	mkdir $outDir unless -e $outDir;
	foreach my $file (keys %files){
		my $oPath=$outDir."/$file";
		my $cPath=$oDir."/$file";
		$oPath=~s/fas.+/bam/;
		$cPath=~s/fas.+/sam/;
		my $command=$samTools." view -bS $cPath > $oPath";
		print LOG "$command\n";
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
	print LOG "Conversion completed\n";
}

sub runBWAconvert {
	my $self	=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $bThreads=$self->{config}->get('OPTIONS','bwaThreads');
	my $iDir	=$self->{config}->get('PATHS','dataDir');
	my $oDir	=$self->{config}->get('PATHS','outDir');
	my $index	=$self->{config}->get('PATHS','bwaIndex');
	my $convert =$self->{config}->get('OPTIONS','converter');
	print LOG "Running BWA conversion....\n";
	foreach my $file (keys %files){
		my $iPath=$iDir."/$file";
		my $oPath=$oDir."/$file";
		my $cPath=$oDir."/$file";
		$oPath=~s/fas.+/bin/;
		$cPath=~s/fas.+/sam/;
		my $command="bwa $convert $index $oPath $iPath > $cPath";
		print LOG "$command\n";
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
	print LOG "Conversion completed\n";
}

sub getAlignmentStats {
	my $self	=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $oDir	=$self->{config}->get('PATHS','outDir');
	print LOG "Collecting alignment numbers....\n";
	foreach my $file (sort {$a cmp $b} keys %files){
		my $oPath=$oDir."/$file";
		$oPath=~s/fas.+/sam/;
		my $thr=threads->create(\&Count,$oPath);
		while(threads->list()>=$mThreads){
			my @thr=threads->list();
			$thr[0]->join();
		}
	}
	print LOG "Counting jobs launched.... waiting....\n";
	while(threads->list()>=1){
		my @thr=threads->list();
		$thr[0]->join();
	}
	print LOG "Alignment counts collected\n";
	return \%return;
}

sub Count {
	my $file=shift;
	my $AOBJ=Alignment->new();
	my $AN=$AOBJ->countHitsInSam($file);
	$return{$file}=$AN;
	return 1;
}



sub getFilesFromConfig {
	my $self=shift;
	foreach my $file ($self->{config}->getAll('FILES')){
		my $misMatches=$self->{config}->get('FILES',$file);
		die "$file did not have a corresponding misMatch value associated with it under 'OPTIONS' in the configuration file!\n". (caller(0))[3] ."\n" unless defined $misMatches;
		$files{$file}=$misMatches;
	}
}

sub findUnaligned {
	my $self=shift;
	my $mThreads=$self->{config}->get('OPTIONS','managerThreads');
	my $bThreads=$self->{config}->get('OPTIONS','bwaThreads');
	my $iDir	=$self->{config}->get('PATHS','dataDir');
	my $oDir	=$self->{config}->get('PATHS','outDir');
	my $index	=$self->{config}->get('PATHS','bwaIndex');
	my $convert =$self->{config}->get('OPTIONS','converter');
	return 1 unless $self->{config}->get('OPTIONS','extractUnaligned');
	print LOG "extracting unaligned reads.... \n";
	foreach my $file (keys %files){
		my $iPath=$iDir."/$file";
		my $cPath=$oDir."/$file";
		my $uPath=$oDir."/$file";
		$cPath=~s/fas.+/sam/;
		$uPath.=".unaligned";
		my $thr=threads->create(\&Unaligned,$iPath,$cPath,$uPath);
		while(threads->list()>=$mThreads){
			my @thr=threads->list();
			$thr[0]->join();
		}
	}
	print LOG "Jobs launched.\n";
	while(threads->list()>=1){
		my @thr=threads->list();
		$thr[0]->join();
	}
	print LOG "Unaligned reads collected.\n";
}

sub Unaligned {
	my $readFile=shift;
	my $alignFile=shift;
	my $output=shift;
	my $AOBJ=Alignment->new();
	$AOBJ->addReadFile($readFile);
	$AOBJ->parseSamToAlignments($alignFile,$readFile);
	my $ratio=$AOBJ->printUnalignedReadsFasta($readFile,$output);
}


#sub loadConfig {
#	my $self=shift;
#my $dir         = $conf->get('PATHS','data_dir');
#my $dbase       = $conf->get('SQL','database');
#my $table       = $conf->get('SQL','table');
#my $convert     = $conf->get('PATHS',"convert");
#	my $confObj=$self->{config};
#	foreach my $id (split(/\,/,$confObj->get('OPTIONS','idList'))){
#	foreach my $id ($self->{config}->getAll('FILES')){
#		print $id."\n";
#	}	
#}

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

