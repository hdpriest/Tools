#!/usr/bin/perl

package GffTools;
use strict;
use warnings;
use lib '/home/hpriest/Scripts/Library';
use GFF3::GffObject;
our $GIDI;
our @GIDS;

sub new {
	my $class=shift;
	my $self = {
		GIDI	=> undef,
		GIDS	=> [],
		Annotation	=> undef,
		NameByID => {},
	};
	bless $self, $class;
	return $self;
}

#scaffold100|size9316    maker   gene    2881    3935    .       -       .       ID=scaffold100|size9316G1;Name=genemark-scaffold100%257Csize9316-abinit-gene-0.3;
#scaffold100|size9316    maker   mRNA    2881    3935    .       -       .       ID=scaffold100|size9316G1-mRNA-1;Parent=scaffold100|size9316G1;Name=genemark-scaffold100%257Csize9316-abinit-gene-0.3-mRNA-1;_AED=1.00;_eAED=1.00;_QI=0|0|0|0|1|1|3|0|250;
#scaffold100|size9316    maker   exon    2881    3345    .       -       .       ID=scaffold100|size9316G1-mRNA-1:exon-1;Parent=scaffold100|size9316G1-mRNA-1;
#scaffold100|size9316    maker   exon    3440    3578    .       -       .       ID=scaffold100|size9316G1-mRNA-1:exon-2;Parent=scaffold100|size9316G1-mRNA-1;
#scaffold100|size9316    maker   exon    3787    3935    .       -       .       ID=scaffold100|size9316G1-mRNA-1:exon-3;Parent=scaffold100|size9316G1-mRNA-1;

sub nameByID { ### Because people can't follow standards
	my $self=shift;
	my $file=shift;
	open(GFF,"<",$file) || die "cannot open GFF file: $file!\n$!\nExiting...\n";
	my $lid=0;
	until(eof(GFF)){
		my $line=<GFF>;
		$lid++;
		next if $line=~m/^\#/;
		next unless $line=~m/\t/;
		chomp $line;
		my ($contig,$source,$type,$start,$stop,$score,$strand,$frame,$comment)=split(/\t/,$line);
		my %comment=%{_parseComment($comment,$type,$lid)};
		my $id=$comment{"ID"};
		my $name;
		if(defined($comment{"Name"})){
			$name=$comment{"Name"};
			$self->{NameByID}{$name}=$id;
		}else{
		}
		$lid++;
	}
	return 1;
}

sub loadGFF {
	my $self=shift;
	if(defined($self->{Annotation})){
		die "Cannot load a GFF if a gff is already loaded!\n";
	}
	$self->{Annotation}={};
	my $file=shift;
	open(GFF,"<",$file) || die "cannot open GFF file: $file!\n$!\nExiting...\n";
	my $lid=0;
	until(eof(GFF)){
		my $line=<GFF>;
		$lid++;
		next if $line=~m/^\#/;
		next unless $line=~m/\t/;
		chomp $line;
		my ($contig,$source,$type,$start,$stop,$score,$strand,$frame,$comment)=split(/\t/,$line);
		next if $type=~m/UTR/;
		my %comment=%{_parseComment($comment,$type,$lid)};
		my $id=$comment{"ID"};
		my $name;
		my $Object=GffObject->new($id,$type,$contig);
		$Object->addStop($stop);
		$Object->addStart($start);
		$Object->addStrand($strand);
		if(defined($comment{"Name"})){
			$name=$comment{"Name"};
		}else{
			$name="NoName";
		}
		$Object->addName($name);
		if($type =~m/gene/){
		}elsif($type eq "chromosome"){
		}elsif($type eq "contig"){
		}else{
			my $p=$comment{"Parent"};
			die "cannot find parent $p for $type line \"$comment{ID}\"\n" unless defined $self->{Annotation}{$p};
			$self->{Annotation}{$p}->addChild($id);
			$Object->addParent($p);
		}
		$self->{Annotation}{$id}=$Object;
	}
	close GFF;
}

sub checkForID {
	my $self=shift;
	my $id=shift;
	return 1 if defined $self->{Annotation}{$id};
	return 0;
}

sub getIDbyName {
	my $self=shift;
	my $name=shift;
	if(defined($self->{NameByID}{$name})){
		return $self->{NameByID}{$name};
	}else{
		return -1;
	}
	return -1;
}

sub addGOannotation {
	my $self=shift;
	my $phrase=shift;
	my $file=shift;
	my @description=@{hdpTools->LoadFile($file)};
	my $missed=0;
	foreach my $line (@description) {
		my ($source,$stanza)=split(/\t/,$line);
		my @GO=split(/\;/,$stanza);
		if(defined($self->{Annotation}{$source})){
			map {$self->{Annotation}{$source}->addGO($_)} @GO;
		}else{
			if(defined($self->{NameByID}{$source})){
				my $ID=getIDbyName($self,$source);
				map {$self->{Annotation}{$ID}->addGO($_)} @GO;
			}else{
				#die "Could not find the object or name for $source\n";
				$missed++;
			}
		}
	}
	warn "$missed entries in $file are not in the definition files\n";
	return 1;
}

sub addFunctionalAnnotation {
	my $self=shift;
	my $phrase=shift;
	my $file=shift;
	my @description=@{hdpTools->LoadFile($file)};
	my $missed=0;
	foreach my $line (@description) {
		my ($source,$stanza)=split(/\t/,$line);
		if(defined($self->{Annotation}{$source})){
			$self->{Annotation}{$source}->addFunctionalAnnotation($phrase,$stanza);
		}else{
			if(defined($self->{NameByID}{$source})){
				my $ID=getIDbyName($self,$source);
				$self->{Annotation}{$ID}->addFunctionalAnnotation($phrase,$stanza);
			}else{
				$missed++;
				#die "Could not find the object or name for:\n$source\n$stanza\n";
			}
		}
	}
	warn "$missed entries in $file are not in the definition files\n";
	return 1;
}

sub addHomology {
	my $self=shift;
	my $phrase=shift;
	my $file=shift;
	my @homology=@{hdpTools->LoadFile($file)};
	foreach my $line (@homology) {
		my ($source,$target,$type,$bit,$e)=split(/\t/,$line);
		if(defined($self->{Annotation}{$source})){
			$self->{Annotation}{$source}->addHomology($phrase,$target,$type);
		}else{
			die "Could not find the object for $source\n";
		}
	}
	return 1;
}

sub _getParentsByID {
	my $self=shift;
	my $id=shift;
	return [] unless defined $self->{Annotation}{$id};
	die "$id is undefind in annotation!\n" unless defined $self->{Annotation}{$id};
	return $self->{Annotation}{$id}->getParents();
}

sub getParentsOfID {
	my $self=shift;
	my $id=shift;
	my $parents=_getParentsByID($self,$id);
	return $parents;
}

sub getAllParentsOfID {
	my $self=shift;
	my $id=shift;
	my @parents;
	my @thisSet=@{_getParentsByID($self,$id)};
	while($#thisSet>=0){
		my $parent_id=shift @thisSet;
		my @newParents=@{_getParentsByID($self,$parent_id)};
		map {push @thisSet, $_} @newParents;
		push @parents, $parent_id;
	}
	return \@parents;
}

sub getAllChildrenOfID {
	my $self=shift;
	my $id=shift;
	my @children;
	my @thisSet=@{_getChildrenByID($self,$id)};
	while($#thisSet>=0){
		my $child_id=shift @thisSet;
		my @newKids=@{_getChildrenByID($self,$child_id)};
		map {push @thisSet, $_} @newKids;
		push @children, $child_id;
	}
	return \@children;
}

sub _getChildrenByID {
	my $self=shift;
	my $id=shift;
	if(defined($self->{Annotation}{$id})){
		return $self->{Annotation}{$id}->getChildren();
	}else{
		return [];
	}
	return [];
}

sub getChildrenOfID {
	my $self=shift;
	my $id=shift;
	my $children=_getChildrenByID($self,$id);
	return $children;
}

sub getObjectByID {
	my $self=shift;
	my $id=shift;
	if(defined($self->{Annotation}{$id})){
		return $self->{Annotation}{$id};		
	}else{
		return undef;
	}
	die "Could not find object with id $id\n";
}

sub getObjectsOfType {	
	my $self=shift;
	my $type=shift;
	my @IDs;
	foreach my $obj (keys %{$self->{Annotation}}){
		if($type eq $self->{Annotation}{$obj}->getType()){
			push @IDs, $obj;
		}
	}
	return \@IDs;
}

sub getListOfAllObjects {
	my $self=shift;
	my @Obj=keys %{$self->{Annotation}};
	return \@Obj;
}

sub _checkCommentFormat {
	my %comment=%{$_[0]};
	my $type=$_[1];
	my $lid=$_[2];
	if($type =~ m/gene/){
		die "$type: No ID defined at line $lid!\n" unless defined $comment{'ID'};
	}elsif($type eq "mRNA"){
		die "$type: No parent defined at line $lid!\n" unless defined $comment{'Parent'};
		die "$type: defined at line $lid!\n" unless defined $comment{'ID'};
	}elsif($type eq "transcript"){
		die "$type: No parent defined at line $lid!\n" unless defined $comment{'Parent'};
		die "$type: No ID defined at line $lid!\n" unless defined $comment{'ID'};
	}elsif($type eq "chromosome"){
	}elsif($type eq "contig"){
	}elsif($type eq "five_prime_UTR"){
	}elsif($type eq "three_prime_UTR"){
	}else{
		die "$type: Not handled\n No parent defined at line $lid!\n" unless defined $comment{'Parent'};
		die "$type: Not handled\n No ID defined at line $lid!\n" unless defined $comment{'ID'};
	}
	return 1;
}

sub _parseComment {
	my $comment=shift;
	my $type=shift;
	my $lid=shift;
	my @comment=split(/\;/,$comment);
	my %comment;
	foreach my $stanza (@comment){
		my ($key,$value)=split(/\=/,$stanza);
		die "Comment was not correctly formatted!\n$comment\nLine: $lid\n" unless((defined($key))&&(defined($value)));
		$comment{$key}=$value;
	}
	_checkCommentFormat(\%comment,$type,$lid);
	return \%comment;
}

return 1;
