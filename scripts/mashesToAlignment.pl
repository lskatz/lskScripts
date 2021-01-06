#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use Bio::Sketch::Mash;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help presence=s absence=s)) or die $!;

  $$settings{presence} //= "1";
  $$settings{absence}  //= "0";

  die usage($settings) if($$settings{help} || !@ARGV);

  # Get the sorted arguments so that they are uniform
  # and not random in each of the methods below,
  # resulting in a determininstic order for the pseudo
  # alignment
  my @infile = sort(@ARGV);

  # Find presence/absence of hashes
  logmsg "Finding presence/absence of ".scalar(@infile)." files";
  my $presence = readSketches(\@infile, $settings);
  # Make sequences out of the hashes
  logmsg "Determining the pseudosequence for each input file";
  logmsg "Present nucleotides will be $$settings{presence} and absent nucleotides will be $$settings{absence}";
  my $seqs = determinePseudoSequences(\@infile, $presence, $settings);
  # make an actual alignment string
  logmsg "Making the alignment from sequence";
  my $aln  = makeAlignment(\@infile, $seqs, $settings);

  print "$aln";

  return 0;
}

sub readSketches{
  my($sketches, $settings) = @_;

  my %p; # presence/absence
  for my $file(@$sketches){
    print STDERR ".";
    my $msh = Bio::Sketch::Mash->new($file);
    my $sketches=$$msh{sketches}[0]{hashes};
    for my $s(@$sketches){
      $p{$s}{$file}=1;
    }
  }
  print STDERR "\n";

  return \%p;
}

sub determinePseudoSequences{
  my($infiles, $p, $settings) = @_;

  my %sampleSeq;

  # Sort to help keep output stable
  my @hashInt = sort{$a<=>$b} keys(%$p);
  for my $h(@hashInt){
    for my $file(@$infiles){
      # If the hash is present, then give the "present" nucleotide
      if($$p{$h}{$file}){
        $sampleSeq{$file} .= $$settings{presence};
      }
      # If the hash is not present, then give the "absent" nucleotide
      else{
        $sampleSeq{$file} .= $$settings{absence};
      }
    }
  }

  return \%sampleSeq;
}

sub makeAlignment{
  my($infiles, $seqs, $settings) = @_;

  my $alnStr = "";
  for my $file(@$infiles){
    $alnStr .= ">$file\n";
    $alnStr .= $$seqs{$file}."\n";
  }
  return $alnStr;
}

sub usage{
  my($settings) = @_;
  print "$0: transforms a set of mash sketches to an alignment
  Usage: $0 [options] *.msh > aln.fasta
  --presence  The nucleotide to use for a present hash integer
              default: $$settings{presence}
  --absence   The nucleotide to use for an absent hash integer
              default: $$settings{absence}
  --help      This useful help menu

  suggested workflow:
  $0 ... 
  goalign reformat phylip -i binary.fasta > binary.fasta.phylip
  raxmlHPC -f a -s binary.fasta.phylip -n \$prefix -T \$numcpus -p \$RANDOM -x \$RANDOM -N 100 -m BINGAMMA
  \n";
}
