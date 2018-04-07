#!/usr/bin/env perl

use strict;
use warnings;
use XML::LibXML::Reader;
use File::Basename qw/basename/;
use Getopt::Long;
use Data::Dumper;
use Getopt::Long;

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help level=i)) or die $!;
  $$settings{level}||=5;

  my($xmlfile,$termsfile)=@ARGV;

  die usage() if($$settings{help});

  my $ontology=readLangual($xmlfile,$settings);
  
  printSpreadsheet($ontology,$termsfile,$settings);

  return 0;
}

sub readLangual{
  my($xmlfile,$settings)=@_;
  my(%parentNode,%childNode,%term,%ontologyID);

  my $xml=XML::LibXML::Reader->new(
    location=>$xmlfile,
  ) or die "ERROR: cannot read $xmlfile";

  # Skip to the descriptors
  $xml->nextElement("DESCRIPTORS");

  # Map child to parent
  # Map term to termID
  # Map synonyms to termID
  while($xml->nextElement("DESCRIPTOR")){
    # <FTC> is the current identifier
    my $ftcXml=XML::LibXML::Reader->new(string=>$xml->readOuterXml());
    $ftcXml->read;
    $ftcXml->nextElement("FTC");
    my $ftc=uc($ftcXml->readInnerXml());
    
    # <BT> is the parent identifier
    my $btXml=XML::LibXML::Reader->new(string=>$xml->readOuterXml());
    $btXml->read;
    $btXml->nextElement("BT");
    my $bt=uc($btXml->readInnerXml());

    my $termXml=XML::LibXML::Reader->new(string=>$xml->readOuterXml());
    $termXml->read;
    $termXml->nextElement("TERM");
    my $term=lc($termXml->readInnerXml());

    $parentNode{$ftc}=$bt;
    $term{$ftc}=$term;
    $ontologyID{$term}=$ftc;
    
    # Other terms: grab the SYNONYMS tag
    my $synXml=XML::LibXML::Reader->new(string=>$xml->readOuterXml());
    $synXml->read;
    while($synXml->nextElement("SYNONYM")){
      $ontologyID{lc($synXml->readInnerXml())}=$ftc;
    }
    #die Dumper \%parentNode,\%term,\%ontologyID if(values(%parentNode) > 6);
    #last if(values(%parentNode) > 10);
  }

  # Now get the child nodes
  while(my($parent,$child)=each(%parentNode)){
    if($childNode{$parent}){
      die "ERROR: $childNode{$parent} is already a child of $parent, but you wanted to assign $child!";
    }
    $childNode{$parent}=$child;
  }

  return({
    parentNode  =>\%parentNode,
    childNode   =>\%childNode,
    term        =>\%term,
    ontologyID  =>\%ontologyID,
  });
}
sub printSpreadsheet{
  my($ontology,$termsfile,$settings)=@_;

  open(my $termsFh,"<", $termsfile) or die "ERROR: cannot read $termsfile: $!";
  while(<$termsFh>){
    chomp;
    my @F=split /\t/;
    my $term=$F[1] || $F[0]; # the term is in the second field if it exists

    $term=lc($term);
    $term=~s/_/ /g;
    $term=~s/^\s+|\s+$//g;
    my $ID=getOntologyIdAtLevel($term,$$settings{level},$ontology,$settings);
    $ID//="";

    # Print the label and the ID but skip the input file's second field
    print join("\t",$F[0],$ID)."\n";
  }

}


######
## helper subroutines

sub getOntologyIdAtLevel{
  my($term,$level,$ontology,$settings)=@_;

  my $ID=$$ontology{ontologyID}{$term};
  if(!$ID && $term=~/^[a-z]\d{4,}/i){
    $ID=uc($term);
  }
  return undef if(!$ID);
  
  my @node=($ID);
  my $parentID;
  do{
    $parentID=$$ontology{parentNode}{$ID};
    $ID=$parentID;
    unshift(@node,$ID);
  } while($parentID && $parentID ne '00000');

  # Return the most specific node possible if the level is too high
  if(@node - 1 < $level){
    return $node[@node - 1];
  } 
  # If possible, return exactly the level requested
  else {
    return $node[$level-1];
  }
}

sub nodeInfo{
  my($xml)=@_;
  printf "%d %d %s %d\n  %s\n", (
    $xml->depth,
    $xml->nodeType, 
    $xml->name, 
    $xml->isEmptyElement, 
    $xml->readInnerXml()
  );
}

sub usage{
  "$0: read a LanguaL XML file and a file with one term per line
  The input file with terms can optionally have a second column that matches LinguaL better

  Usage: $0 langual.xml terms.txt > mapping.tsv
  --level   2   What level to descend to in the ontology
  "
}
