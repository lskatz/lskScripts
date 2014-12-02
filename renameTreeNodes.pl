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
  GetOptions($settings,qw(help map=s delimiter=s suffix=s tailSuffix=s daysRecent=i markup=s discard-old newmap=s)) or die;
  die usage() if($$settings{help});
  my $mapFile=$$settings{map} || die "ERROR: need the map file\n".usage();
  $$settings{delimiter}||='|';
  $$settings{suffix}||="";
  $$settings{daysRecent}||=0;
  $$settings{markup}||='!!';
  $$settings{'discard-old'}||=0;
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
  #my $out=Bio::TreeIO->new(); # stdout
  my %seen;
  my $i=0; # node counter
  my $deleteCounter=0;
  my $NEWMAP;
  if($$settings{newmap}){ open($NEWMAP,">",$$settings{newmap}) or die "ERROR: could not open $$settings{newmap} for writing:$!";}
  while(my $tree=$in->next_tree){
    for my $node($tree->get_leaf_nodes){
      #if ($i++>100){
      #  logmsg "DEBUG";
      #  last;
      #}
      my $oldid=$node->id || "";
      $oldid=~s/^'+|'+$//g; # remove quotes from the ID
      $oldid=~s/^"+|"+$//g; # remove quotes from the ID
      $oldid=~s/^\Q$d\E//g;  # remove any leading delimiters
      next if(!$oldid);

      # split/grep:
      # Step 1: Split by the delimiter
      # Step 2: Only retrieve elements that have at least one non-whitespace character
      # Step 3: Return the zeroth element: that is now $oldid
      # Step 4: Go to the next iteration in the loop if that doesn't return anything
      $oldid=(grep{/\S/} split (/\Q$d\E/, $oldid))[0];
      next if(!defined($oldid) || $oldid eq "");

      $oldid=~s/\Q$suffix\E$//g if($suffix);
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

      # print to the TSV if this is a new identifier that should be saved locally
      if($newid ne $oldid && $$settings{newmap}){
        printNewmapId($oldid,$newid,$NEWMAP,$settings);
      }

      if($seen{$oldid}++){
        logmsg "Warning: I have already seen this id in the tree and so it appears more than once\n  $newid";
        $newid="$newid$$settings{delimiter}DuplicateNumber$seen{$oldid}";
      }

      # make some kind of marking on this taxon if it's new-ish
      if(isDocNew($newid,$$settings{daysRecent},$settings)){
        logmsg "Found a new sample and will mark it\n  $newid";
        $newid="$$settings{markup} $newid";
      } 
      # If the DOC is not new and the user only wants new nodes
      elsif($$settings{'discard-old'}) {
        logmsg "Discarding because it is older than $$settings{daysRecent} days old\n  $newid";
        #$tree->remove_Node($node);
        #$tree->splice(-remove_id=>[$oldid],-preserve_lengths=>0);
        $newid="-----".++$deleteCounter;
      }

      die "ERROR: $oldid => $newid has illegal characters: $1" if($newid=~/(,|\(|\)|:|;)+/);

      $node->id($newid);
    }
    print $tree->as_text('newick')."\n";next;
    #$out->write_tree($tree);
    #print "\n";
  }
  logmsg scalar(keys(%seen))." entries were converted";
  close $NEWMAP if($$settings{newmap});
  return 1;
}

sub printNewmapId{
  my($oldid,$newid,$NEWMAP,$settings)=@_;
  # TODO consider if there is an escaped delimiter like \|
  my $d=$$settings{delimiter};
  my @F=split(/\Q$d\E/,$newid);
  my $line=join("\t",$oldid,@F)."\n";
  print $NEWMAP $line;
  return length($line);
}

# go onto NCBI and grab information pertaining to this identifier.
sub ncbiInfo{
  my($oldid,$settings)=@_;
  my $newid;

  # Query for a new identifier using various methods
  if($oldid=~/^GC\w/){
    $newid=assemblyInfo($oldid,$settings);
  } elsif($oldid=~/^SAM\w/){
    $newid=biosampleInfo($oldid,$settings);
  } elsif($oldid=~/^PNUSA\w/){
    $newid=biosampleInfoFromName($oldid,$settings);
  } elsif($oldid!~/^\d+$/){ # if it's not just an integer
    $newid=biosampleInfo($oldid,$settings);
  } else{
    logmsg "Warning: I do not know how to query $oldid. (maybe it's just a bootstrap value)";
  }

  # Remove any characters that might mess up the output Newick format: ,;:()
  $newid=~s/[\,;\(\)\:]|\s+/_/g;

  # if the identifier is just a number then it didn't get converted correctly. Revert!
  my $d=$$settings{delimiter};
  my $tmpid=$newid;
  $tmpid=~s/\s+//g;    # remove whitespace
  $tmpid=~s/\Q$d\E//g; # remove delimiters
  if($tmpid=~/^\d+$/){ # if it's only numbers
    logmsg "Reverted! $newid back to $oldid";
    $newid=$oldid;
  }

  return $newid;
}

sub biosampleInfoFromName{
  my($oldid,$settings)=@_;
  return "" if(!defined($oldid) || !$oldid);
  my $query=$oldid;
  my $retryConnection=0;  # how many times to try if a connection times out

  my $eutil=Bio::DB::EUtilities->new(
      -eutil => 'esearch',
      -db    => 'BioSample',
      -term  => $oldid,
      -email => "nobody\@cdc.gov",
      -retmax => 1, # just in case it returns a huge amount when in reality it should be 1
  );
  my $xml=$eutil->get_Response->content;
  my $p=XML::LibXML->new;
  my $doc=$p->parse_string($xml);
  my $uid;
  foreach my $node ($doc->findnodes("eSearchResult/IdList/Id")){
    $uid=$node->textContent;
  }
  # $uid is an integer which is a stand-in for SAMN. Find it through biosampleinfo
  die "Internal error: uid cannot be a non-integer but it is: $uid" if($uid=~/\D/);
  logmsg "$oldid => $uid";
  return biosampleInfo($uid,$settings);
}

sub biosampleInfo{
  my($oldid,$settings)=@_;
  return "" if(!defined($oldid) || !$oldid);
  # UID is the identifier without SAMN
  my $uid=$oldid; $uid=~s/^SAMN//i;
  my $retryConnection=0;  # how many times to try if a connection times out

  EFETCHBIOSAMPLE:
  # make the query
  # e.g. http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=02389791
  my $eutil=Bio::DB::EUtilities->new(
           -eutil => 'efetch',
           -db => 'biosample',
           -id=>$uid,
           -retmax => 10, # just in case it returns a huge amount when in reality it should be 1
           -email => "nobody\@cdc.gov",
  );
  my $xml;
  eval{
    $xml=$eutil->get_Response->content;
  };
  if($@) {
    die "Server error querying $oldid\n $@." if(++$retryConnection > 5);
    logmsg "Server error. Retry number $retryConnection\n  $@";
    #sleep $retryConnection; # sleep as many seconds as is the retry, so that you can give the server time to get back to normal
    goto EFETCHBIOSAMPLE;
  }
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

  # Nothing found; see if anything else matches.
  # Don't do this blanket-search on integer identifiers. They are probably SAMN identifiers and need to be treated differently.
  if(keys(%attr)<1 && $oldid !~/^\d+$/){
    $eutil=Bio::DB::EUtilities->new(
      -eutil => 'esearch',
      -db    => 'BioSample',
      -term  => $oldid,
      -email => "nobody\@cdc.gov",
      -retmax => 1, # just in case it returns a huge amount when in reality it should be 1
    );
    my $xml=$eutil->get_Response->content;
    my $p=XML::LibXML->new;
    my $doc=$p->parse_string($xml);
    my $uid;
    foreach my $node ($doc->findnodes("eSearchResult/IdList/Id")){
      $uid=$node->textContent;
    }
    # $uid is an integer which is a stand-in for SAMN. Find it through biosampleinfo
    die "Internal error: uid cannot be a non-integer but it is: $uid" if($uid=~/\D/);
    return biosampleInfo($uid,$settings);
  }

  # Nothing found
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
  my $retryConnection=0;

  # get the UID
  ASSEMBLYESEARCH:
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
  #return "" if(!@link || !grep(/biosample/,@link));

  # Ok now hopefully we have a SAMN biosample identifier.
  if(grep(/biosample/,@link)){
    my $linkid=($doc2->findnodes("eLinkResult/LinkSet/LinkSetDb/Link/Id"))[0]->textContent;
    my $samn="SAMN$linkid";

    my $biosampleNewid=biosampleInfo($samn,$settings);
    return $biosampleNewid if($biosampleNewid && $biosampleNewid ne $oldid);
  }

  # If there is no biosample identifier then just get what is in the assembly NCBI database
  $eutil=Bio::DB::EUtilities->new(
           -eutil => 'esearch',
           -db => 'genome',
           -term=>$oldid,
           -retmax => 10, # just in case it returns a huge amount when in reality it should be 1
           -email => "nobody\@cdc.gov",
  );
  my $xml3=$eutil->get_Response->content;
  my $p3=XML::LibXML->new;
  my $doc3=$p->parse_string($xml3);
  @link=($doc3->findnodes("eSearchResult/IdList/Id"));
  my $genomeIdentifier=$link[0]->textContent;

  # Get information from the genome db and return an identifier with that information
  $eutil=Bio::DB::EUtilities->new(
          -eutil => 'esummary',
          -db    => 'genome',
          -id    => $genomeIdentifier,
          -email => "nobody\@cdc.gov",
  );
  my $ds=$eutil->next_DocSum;
  my %genomeInfo;
  while(my $item=$ds->next_Item('flattened')){
    $genomeInfo{$item->get_name}=$item->get_content || "";
  }
  my $assemblyId=join($$settings{delimiter},$genomeInfo{Organism_Name},$genomeInfo{DefLine},$genomeInfo{Create_Date},$oldid);
  logmsg "New ID found through the genome database API: $assemblyId";
  return $assemblyId;
}

# See if the date of collection (DOC) is newer than XX days
# Returns 0 if it's old; the number of days if it is new.
sub isDocNew{
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
  }
  # if it is X days old, then make a mark on $newid
  my $daysOld=ceil((time - $newestTimestamp) / 3600 / 24);
  if($daysOld < $daysRecentThreshold){
    return $daysOld;
  }
  return 0;
}

sub usage{
  "Transforms a dnd tree such that the nodes are renamed according to a spreadsheet.
  Usage: $0 tree.dnd -m map.tsv > converted.dnd
  -map map.tsv
    The map file is tab-delimited with two columns: FROM  TO [TO_2...]
    Where any node in the tree with the FROM name will be converted to the TO name
    Multiple TO fields will be concatenated with the delimiter (given by -d)
  --newmap new.tsv  Create a new TSV with the information you've downloaded, so that the process goes more quickly next time.
  --delimiter delimiter  Default:|
    A delimiter for the FROM tree. The first field with be used to convert the identifier.
  --daysRecent 0         The number of days since today for a genome to be 'recent'
  --discard-old          To discard taxa that are older than 'days'
  --markup '!!'          The markup given to new isolates in the format of '!!%d!!' with the number of days shown where %d is
  -s suffix              Default: nothing
    If set, this suffix will be trimmed off the IDs in the source phylogeny
  --tailSuffix character Everything after this character will be removed
  "
}
