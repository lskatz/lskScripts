#!/usr/bin/env perl
# Downloads new reads from a given bioproject
# Author: Lee Katz <lkatz@cdc.gov>

use FindBin qw($Bin);
use lib "$Bin/lib";

use strict;
use warnings;
use Fatal;
use Data::Dumper;
use Getopt::Long;
#use XML::Parser;
use LWP::Simple;
use File::Basename;
use Schedule::SGELK;

$0=fileparse $0;
local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDOUT "$0: ".(caller(1))[3].": @_\n";}

my $sge=Schedule::SGELK->new(warn_on_error=>1);
exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(project=s inbox=s force help maxslots=i));
  die usage() if($$settings{help});
  $$settings{project} || die "ERROR: need project\n".usage();
  $$settings{inbox}   || die "ERROR: need inbox\n".usage();
  mkdir $$settings{inbox} if(!-d $$settings{inbox});
  $$settings{"maxslots"}||=1;
  $sge->set("maxslots",$$settings{"maxslots"});

  logmsg "Getting SRA Experiment IDs for this project";
  my @sraId=getSraIds($$settings{project},$settings);
  logmsg "prefetching ".scalar(@sraId)." read sets. Fortunately prefetch does know whether something has already been prefetched.";
  prefetch(\@sraId,$settings);
  logmsg "Getting run IDs for the SRA Experiments";
  my @srrId=getSrrIds(\@sraId,$settings);
  logmsg "Found ".scalar(@srrId)." read sets associated with $$settings{project}";
  logmsg "Downloading...";
  downloadReads(\@srrId,$$settings{inbox},$settings);
  return 0;
}

sub getSraIds{
  my($project,$settings)=@_;
  $project=~s/^PRJNA//i; # strip the "PRJNA" prefix if present
  logmsg "Downloading";
  my $xml;
  while(!defined $xml){
    $xml=get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=bioproject&db=SRA&id=".$project);
  }
  my @xml=split(/\n/,$xml);
  my @id;
  my $is_linkSetDb=0;
  for my $line(@xml){
    $is_linkSetDb=1 if($line=~/<LinkSetDb>/);
    $is_linkSetDb=0 if($line=~/<\/LinkSet>/);
    next if(!$is_linkSetDb);
    next if($line !~/<Id>(\d+)<\/Id>/);
    push(@id,$1);
  }
  return @id if wantarray;
  return \@id;
}
sub getSrrIds{
  my($sraId,$settings)=@_;
  my @id;
  for (@$sraId){
    #logmsg "Downloading information on $_";
    my $xml=get("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?Db=SRA&id=".$_);
    my @xml=split(/\n/,$xml);
    for my $line(@xml){
      while($line=~/<PRIMARY_ID>(SRR\d+)<\/PRIMARY_ID>/g){
        push(@id,$1);
      }
    }
  }
  return @id if wantarray;
  return \@id;
}

sub prefetch{
  my($id,$settings)=@_;
  for(my $i=0;$i<@$id;$i++){
    system("prefetch $$id[$i]");
    die "ERROR with prefetch command from sra-tools" if $?;
  }
  return scalar(@$id);
}

sub downloadReads{
  my($srrId,$inbox,$settings)=@_;
  for my $id (@$srrId){
    my($read1,$read2,$shuffled)=("$inbox/$id"."_1.fastq","$inbox/$id"."_2.fastq","$inbox/$id.shuffled.fastq.gz");
    next if(!$$settings{force} && -e $shuffled);
    logmsg "Downloading and formatting $id";
    $sge->pleaseExecute("fastq-dump -O $inbox/ -I --split-files $id && run_assembly_shuffleReads.pl $read1 $read2 | gzip -c > $shuffled.tmp && rm -v $read1 $read2 && mv $shuffled.tmp $shuffled -v",{jobname=>"sra$id",numcpus=>1});
    die if $?;
  }
  $sge->wrapItUp();
}

sub usage{
  "Downloads new reads from a given bioproject
  usage: $0 -p bioprojectId -i inbox/
  -p The bioproject ID can be the numerical identifier or the PRJN number
  -i inbox for where reads will go
  -f for force. Otherwise, it will not re-download reads
  -m 1 maximum slots to use (on a cluster)
  "
}
