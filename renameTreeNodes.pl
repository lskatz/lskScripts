#!/usr/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Bio::DB::EUtilities;
use XML::LibXML;
use POSIX qw(strftime ceil floor);
use Date::Parse qw(str2time);
use Time::ParseDate qw/parsedate/;

$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";};
exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help map=s delimiter=s suffix=s tailSuffix=s daysRecent=i markup=s)) or die;
  die usage() if($$settings{help});
  my $mapFile=$$settings{map} || die "ERROR: need the map file\n".usage();
  $$settings{delimiter}||='|';
  $$settings{suffix}||="";
  $$settings{daysRecent}||=0;
  $$settings{markup}||='!!';
  my $inTree=$ARGV[0] || die "ERROR: need tree file\n".usage();
  
  my $map=readMap($mapFile,$settings);
  printNewTree($inTree,$map,$settings);
  return 0;
}

sub readMap{
  my($map,$settings)=@_;
  my $d=$$settings{delimiter} || die "ERROR: there is no delimiter\n".usage();
  my %map;
  open(MAP,$map) or die "ERROR: could not read map file $map: $!";
  logmsg "Reading $map";
  while(<MAP>){
    s/,//g; # don't let commas in the outbreak codes mess with the code for the tree
    chomp;
    my @F=split /\t/;
    my $numF=@F;
    if($numF<2){
      logmsg "Nothing found on this line. Will not convert for it.\n  $_";
      next;
    }
    $map{$F[0]}=join($d,@F[1..$numF-1]);
    #$map{$F[0]}=$F[1]; # DEBUG
  }
  close MAP;
  return \%map;
}

sub printNewTree{
  my($inTree,$map,$settings)=@_;
  my $d=$$settings{delimiter} || die "ERROR: there is no delimiter\n".usage();
  my $suffix=$$settings{suffix};
  my $in=Bio::TreeIO->new(-file=>$inTree);
  my $out=Bio::TreeIO->new(); # stdout
  my %seen;
  my $i=0;
  while(my $tree=$in->next_tree){
    for my $node($tree->get_nodes){
      #if ($i++>100){
      #  logmsg "DEBUG";
      #  last;
      #}
      my $oldid=$node->id || "";
      next if(!$oldid);
      $oldid=~s/^'+|'+$//g; # remove quotes from the ID
      $oldid=~s/^"+|"+$//g; # remove quotes from the ID
      $oldid=~s/^\Q$d\E//;  # remove any leading delimiters

      $oldid=(split /\Q$d\E/, $oldid)[0];
      $oldid=~s/\Q$suffix\E$// if($suffix);
      if(defined(my $tail=$$settings{tailSuffix})){
        $oldid=~s/\Q$tail\E.*$//;
      }
      my $newid; # what the node will be turned into
      if(!defined($oldid) || !defined($$map{$oldid})){
        # if it is defined but just not in the mappings, then produce a warning
        if($oldid && !defined($$map{$oldid})){
          # see if this identifier can be found on NCBI
          $newid=ncbiInfo($oldid,$settings);
          logmsg "Warning: I do not understand $oldid and it was left unconverted" if(!$newid);
        }
      } 
      # if the new node name can be found in the mappings then grab it.
      else {
        $newid=$$map{$oldid};
      }
      next if(!$newid);
      next if($seen{$oldid}++); # I think that somehow nodes were being parsed twice

      # make some kind of marking on this taxon if it's new-ish
      $newid=highlightIfNew($newid,$$settings{daysRecent},$settings);

      $node->id($newid);
    }
    $out->write_tree($tree);
    print "\n";
  }
  logmsg scalar(keys(%seen))." entries were converted";
  return 1;
}

# go onto NCBI and grab information pertaining to this identifier.
sub ncbiInfo{
  my($oldid,$settings)=@_;
  my $newid;
  if($oldid =~/^(SAMN|CFSAN)/){
    $newid=biosampleInfo($oldid,$settings);
  } elsif($oldid=~/^GC\w/){
    $newid=assemblyInfo($oldid,$settings);
  } else{
    logmsg "Warning: I do not know how to query $oldid.";
  }
  return $newid;
}

sub biosampleInfo{
  my($oldid,$settings)=@_;
  # UID is the identifier without SAMN
  my $uid=$oldid; $uid=~s/^SAMN//i;

  # make the query
  # e.g. http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=02389791
  my $eutil=Bio::DB::EUtilities->new(
           -eutil => 'efetch',
           -db => 'biosample',
           -id=>$uid,
           -retmax => 10, # just in case it returns a huge amount when in reality it should be 1
           -email => "nobody\@cdc.gov",
  );
  my $xml=$eutil->get_Response->content;
  my $p=XML::LibXML->new;
  my $doc=$p->parse_string($xml);

  # set up key/value pairs from the XML attributes
  my %attr;
  foreach my $node($doc->findnodes("BioSampleSet/BioSample/Attributes/Attribute")){
    my $text=$node->textContent;
    #$text=~s/[,\s\|;\.:]+/_/g; # remove any offensive characters
    $text=~s/[^\w\d\-\/]+/_/g; # remove any offensive characters
    for my $attr($node->attributes){
      $attr{$attr->value}=$text;
    }
  }

  # not sure if all these headers are even in the attributes
  my @headers=('strain', 'strain', 'SAMN', 'collection date', 'geographic location (country and/or sea,region)', 'strain', 'isolation source', 'serovar', 'PFGE_PrimaryEnzyme_pattern', 'PFGE_SecondaryEnzyme_pattern', 'Outbreak'); 
  my $newid;
  #die Dumper [$xml,\%attr] if($uid eq 1125831); # see what's happening with GCA_000183845.1
  if(keys(%attr)<1){
    return $newid;
  }
  
  $attr{SAMN}=$oldid;
  $newid.=($attr{$_} || "").'|' for(@headers);
  $newid="" if($newid!~/[^\|]/); # if it's just pipe characters, then set it back to blank
  $newid=~s/\|$//; # remove trailing pipe
  logmsg "New ID found for $oldid through the BioSample API =>\n  $newid";

  return $newid;
}

sub assemblyInfo{
  my($oldid,$settings)=@_;
  #my $uid=$oldid; $uid=~s/^SAMN//i;

  # get the UID
  my $eutil=Bio::DB::EUtilities->new(
           -eutil => 'esearch',
           -db => 'assembly',
           -term=>$oldid,
           -retmax => 1, # just in case it returns a huge amount when in reality it should be 1
           -email => "nobody\@cdc.gov",
  );
  my $xml=$eutil->get_Response->content;
  my $p=XML::LibXML->new;
  my $doc=$p->parse_string($xml);
  my $uid;
  foreach my $node ($doc->findnodes("eSearchResult/IdList/Id")){
    $uid=$node->textContent;
  }
  die "ERROR: no UID for $oldid" if(!$uid);

  # get the ID to go to BioSample
  my $search=Bio::DB::EUtilities->new(
            -eutil=>'elink',
            -db=>"biosample",
            -dbfrom=>"assembly",
            -id=>$uid,
            -email=>"nobody\@cdc.gov",
  );
  my $xml2=$search->get_Response->content;
  my $p2=XML::LibXML->new;
  my $doc2=$p2->parse_string($xml2);
  my @link=($doc2->findnodes("eLinkResult/LinkSet/LinkSetDb/DbTo"));
  return "" if(!@link || !grep(/biosample/,@link));
  my $linkid=($doc2->findnodes("eLinkResult/LinkSet/LinkSetDb/Link/Id"))[0]->textContent;
  my $samn="SAMN$linkid";
  
  #die Dumper [$oldid,$uid,$linkid,$samn];
  return biosampleInfo($samn,$settings);

}

sub highlightIfNew{
  my($newid,$daysRecentThreshold,$settings)=@_;
  my @F=split(/\|/,$newid);

  # if any of these fields look like dates, judge how old they are
  my $newestTimestamp=1; # initialized to look "old"
  #logmsg $newid;
  for(@F){
    next if($_ !~/[\/\-]/); # all dates need to have a slash or dash for this to work
    my $timestamp=parsedate($_,UK=>0) || 0;
    next if($timestamp==0 || !$_); # zero indicates that it didn't parse as a date
    $newestTimestamp=$timestamp if($timestamp > $newestTimestamp);

    ## Debuging to figure out if the timestamps are generated correctly
    #my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($newestTimestamp);
    #$year+=1900;
    #$mon++;
    #my $date="$year-$mon-$mday";
    #logmsg "$_ => $timestamp => $date";
  }

  # if it is X days old, then make a mark on $newid
  my $daysOld=ceil((time - $newestTimestamp) / 3600 / 24);
  if($daysOld < $daysRecentThreshold){
    logmsg "Found a new sample and will mark it\n  $newid";
    my $markup=$$settings{markup};
    $newid="$markup$daysOld$markup $newid";
  }

  return $newid;
}

sub usage{
  "Transforms a dnd tree such that the nodes are renamed according to a spreadsheet.
  Usage: $0 tree.dnd -m map.tsv > converted.dnd
  -map map.tsv
    The map file is tab-delimited with two columns: FROM  TO [TO_2...]
    Where any node in the tree with the FROM name will be converted to the TO name
    Multiple TO fields will be concatenated with the delimiter (given by -d)
  --delimiter delimiter Default:|
    A delimiter for the FROM tree. The first field with be used to convert the identifier.
  -days 0 The number of days since today for a genome to be 'recent'
  --markup '!!' the markup given to new isolates in the format of '!!%d!!' with the number of days shown where %d is
  -s suffix Default: nothing
    If set, this suffix will be trimmed off the IDs in the source phylogeny
  --tailSuffix character. Everything after this character will be removed
  "
}
