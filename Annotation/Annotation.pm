#!/usr/bin/perl

package Annotation;

use warnings;
use strict;
use lib '/home/hpriest/Scripts/Library';
use GFF3::GffTools;
use hdpTools;
use Configuration;

sub new {
        my $class=shift;
	  my $config=shift;
	  die "Can't init an Annotation object without a configuration file!\n" unless defined $config;
	  my $C=Configuration->new($config);
	  _checkConfig($C);
        my $self = {
		config	=> $C,
        };
        bless $self, $class;
	  
        return $self;
}

sub getAllRelatives {
	my $self=shift;
	my $ID=shift;
	if($self->{annotation}->checkForID($ID)){
	}else{
		if($self->{annotation}->getIDbyName($ID)){
			$ID=$self->{annotation}->getIDbyName($ID);
		}else{
			die "Cannot find this ID in annotation:\n$ID\n";
		}
	}
	my $children=$self->{annotation}->getAllChildrenOfID($ID);
	my $parents =$self->{annotation}->getAllParentsOfID($ID);
	my @object_ids=(@$children,@$parents);
	unshift @object_ids, $ID;
	return \@object_ids;
}

sub getAnnotationOfEntireTree {
	my $self=shift;
	my $ID=shift;
	if($self->{annotation}->checkForID($ID)){
	}else{
		if($self->{annotation}->getIDbyName($ID)){
			$ID=$self->{annotation}->getIDbyName($ID);
		}else{
			die "Cannot find this ID in annotation:\n$ID\n";
		}
	}
	my $children=$self->{annotation}->getAllChildrenOfID($ID);
	my $parents =$self->{annotation}->getAllParentsOfID($ID);
	my @object_ids=(@$children,@$parents);
	unshift @object_ids, $ID;
	my $finalID;
	my %func;
	my %go;
	foreach my $obj_id (@object_ids){
		next	unless(defined($self->{annotation}->getObjectByID($obj_id)));
		my $obj=$self->{annotation}->getObjectByID($obj_id);
		my $IAM=$obj->getType();
		if(($IAM =~ m/mrna/i) || ($IAM=~m/gene/i)){
			$finalID=$obj_id if $IAM=~m/gene/i;
			my %Func=%{$obj->getFunctionalAnnotation()};
			my @GO  =@{$obj->getGOAnnotation()};
			map {$go{$_}=1} @GO;
			my @z=keys %go;
			if(scalar(@z)==0){
				$go{"N/A"}=1;
			}else{
			}
			foreach my $k (keys %Func){
				if(defined($func{$k})){
					$func{$k}{$Func{$k}}=1;
				}else{
					$func{$k}={};
					$func{$k}{$Func{$k}}=1;
				}
			}
		}else{
		}
	}
	my %h;
	$go{"N/A"}=1 unless ((scalar(keys(%go)))>0);
	my @k=keys %go;
	$h{'go'}=\@k;
	$h{'func'}=\%func;
	return \%h;
}

sub getObjByID {
	my $self=shift;
	my $id=shift;
	return $self->{annotation}->getObjectByID($id);
}

sub getAllIDs {
	my $self=shift;
	return $self->{annotation}->getListOfAllObjects();
}

sub getAllInformationByID {
	my $self=shift;
	my $id=shift;
	if(my $obj=$self->{annotation}->getObjectByID($id)){
		return $obj->getAllInformation();
	}else{
		die "Cannot find information for ID : $id\n";
	}
	return -1;
}

sub getAllParentsOfID {
	my $self=shift;
	my $id=shift;
	return $self->{annotation}->getAllParentsOfID($id);
}

sub getAllChildrenOfID {
	my $self=shift;
	my $id=shift;
	return $self->{annotation}->getAllChildrenOfID($id);
}

sub loadGenome {
	my $self=shift;
	my $path=$self->{config}->get('PATHS','AnnotationDir')."/".$self->{config}->get('SEQUENCEFILES','Genome');
	$self->{genome}=hdpTools->LoadFasta($path);
	return 1;
}

sub IndexByName {
	my $self=shift;
	my $path=$self->{config}->get('PATHS','AnnotationDir')."/".$self->{config}->get('ANNOTATION','GFF');
	$self->{annotation}->nameByID($path);
	return 1;
}

sub checkForID {
	my $self=shift;
	my $id=shift;
	return 1 if $self->{annotation}->checkForID($id);
	return 0;
}

sub getIDbyName {
	my $self=shift;
	my $name=shift;
	return $self->{annotation}->getIDbyName($name);
}

sub loadAnnotation {
	my $self=shift;
	my $path=$self->{config}->get('PATHS','AnnotationDir')."/".$self->{config}->get('ANNOTATION','GFF');
	$self->{annotation}=GffTools->new();
	$self->{annotation}->loadGFF($path);
	return 1;
}

sub getHomologyTargets {
	my $self=shift;
	my @HomologyPhrases=$self->{config}->getAll('HOMOLOGY');
	return \@HomologyPhrases;
}

sub addGOannotations {
	my $self=shift;
	my @HomologyPhrases=$self->{config}->getAll('FUNCTIONAL');
	foreach my $phrase (@HomologyPhrases){
		my $path=$self->{config}->get('PATHS','AnnotationDir')."/".$self->{config}->get('FUNCTIONAL',$phrase);
	warn $path."\n";
		die "cannot find file for functional descriptions from $phrase!\n" unless -e $path;
		$self->{annotation}->addGOannotation($phrase,$path);
	}
	return 1;
}
sub getFunctionalKeys {
	my $self=shift;
	my @HomologyPhrases=$self->{config}->getAll('DESCRIPTIONS');
	return \@HomologyPhrases;
}

sub addDescriptions {
	my $self=shift;
	my @HomologyPhrases=$self->{config}->getAll('DESCRIPTIONS');
	foreach my $phrase (@HomologyPhrases){
		my $path=$self->{config}->get('PATHS','AnnotationDir')."/".$self->{config}->get('DESCRIPTIONS',$phrase);
		die "cannot find file for functional descriptions from $phrase!\n" unless -e $path;
		$self->{annotation}->addFunctionalAnnotation($phrase,$path);
	}
	return 1;
}

sub getHomologySources {
	my $self=shift;
	my @HomologyPhrases=$self->{config}->getAll('HOMOLOGY');
	return \@HomologyPhrases;
}

sub addHomologies {
	my $self=shift;
	my @HomologyPhrases=$self->{config}->getAll('HOMOLOGY');
	foreach my $phrase (@HomologyPhrases){
		my $path=$self->{config}->get('PATHS','AnnotationDir')."/".$self->{config}->get('HOMOLOGY',$phrase);
		die "cannot find file for homology to $phrase!\n" unless -e $path;
		$self->{annotation}->addHomology($phrase,$path);
	}
	return 1;
}

sub _checkConfig {
	my $Config=shift;
	die "Annotation directory undefined in config file!\n" unless(-e $Config->get('PATHS','AnnotationDir'));
	my $AnnoDir= $Config->get('PATHS','AnnotationDir');
	die "No GFF defined.\n" 					unless($Config->get('ANNOTATION','GFF'));
	die "Cannot find GFF file.\n" 				unless(-e $AnnoDir."/".$Config->get('ANNOTATION','GFF'));
	die "No Genome defined.\n" 					unless($Config->get('SEQUENCEFILES','Genome'));
	die "Cannot find Genome Fasta file.\n" 			unless(-e $AnnoDir."/".$Config->get('SEQUENCEFILES','Genome'));
	warn "No Protein fasta defined.\n" 				unless($Config->get('SEQUENCEFILES','Proteins'));
	warn "No Promoter fasta defined.\n" 			unless($Config->get('SEQUENCEFILES','Promoters'));
	warn "No Gene fasta defined.\n" 				unless($Config->get('SEQUENCEFILES','Genes'));
	warn "No CDNA fasta defined.\n" 				unless($Config->get('SEQUENCEFILES','CDNA'));
	warn "No Primary GO annotation file defined.\n" 	unless($Config->get('FUNCTIONAL','PrimaryGO'));
	warn "No Module Lists Defined\n" 				unless($Config->get('PATHS','ModuleDir'));
	return 1;
}


1;
