#!/usr/bin/env perl
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Bio::Perl;
use Getopt::Long;

exit(main());

sub main{
  my $settings={};
  die usage() if(!@ARGV);
  GetOptions($settings,qw(flanking=i name=s));
  $$settings{flanking}||=0;
  $$settings{revcom}  ||=0;
  my ($db,$query)=@ARGV[0,1];
  die "ERROR: need db and query\n".usage() if(!$query || !$db);
  
  my($contig,$start,$stop)=blastAgainstDb($db,$query,$settings);
  my $seq=extractSeq($contig,$start,$stop,$db,$settings);
  print "$seq\n";
  return 0;
}

sub extractSeq{
  my($contig,$start,$stop,$db,$settings)=@_;
  my $fastaStr;
  if($start>$stop){
    my $tmp=$start;
    $start=$stop;
    $stop=$tmp;
    $$settings{revcom}=1;
  }
  my $in=Bio::SeqIO->new(-file=>$db);
  while(my $seq=$in->next_seq){
    next if($seq->id ne $contig);
    $start=$start-$$settings{flanking};
    $stop=$stop+$$settings{flanking};
    $start=1 if($start<1);
    $stop=$seq->length if($stop>$seq->length);
    my $sequence=$seq->subseq($start,$stop);
    $sequence=~tr/ATCGatcg/TAGCtagc/ if($$settings{revcom});
    $sequence=reverse($sequence) if($$settings{revcom});
    $sequence=~s/(.{60})/$1\n/g;
    my $id=join("_",">$$settings{name}".$seq->id,$start,$stop);
    $fastaStr="$id\n$sequence";
  }
  die "Could not create a fasta in contig $contig with start/stop=$start/$stop\n" if(!$fastaStr);
  return $fastaStr;
}

sub blastAgainstDb{
  my($db,$query,$settings)=@_;
  my $command="legacy_blast.pl blastall -i '$query' -d '$db' -a 12 -e 0.05 -m 8 -p blastn -v 5 -b 5";
  my @result=split(/\n/,`$command`);
  die if(!@result);
  @result=sort{
    (split(/\t/,$b))[2]<=>(split(/\t/,$a))[2]
  } @result;
  
  my($contig,$start,$stop)=(split(/\t/,$result[0]))[1,8,9];
  return ($contig,$start,$stop);
}

sub usage{
  "Blasts a nucleotide sequence against a database and extracts the hit
  Usage: $0 database.fasta query.fasta > hit.fasta
    -f 100 to extract 100bp upstream/downstream
    -n custom genome name to put into the defline
  "
}
