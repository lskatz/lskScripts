#!/usr/bin/env perl

# run-assembly-filterContigs: remove contigs that fit certain criteria
# Author: Lee Katz (lkatz@cdc.gov)
# Taken from CG-Pipeline and modified as a standalone script

use strict;
use warnings;

use Getopt::Long qw/GetOptions/;
use File::Temp ('tempdir');
use File::Basename qw/basename/;
use Data::Dumper;

$0 = basename($0);
#local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
sub logmsg {print STDERR "$0: @_\n";}

exit(main());

sub main {
	my $settings = {};

	my @cmd_options = qw(verbose cov|min_cov=i assembly_min_contig_length|length=i help numcpus=i rename);
	GetOptions($settings, @cmd_options) or die;
	$$settings{numcpus}||=1;
  $$settings{assembly_min_contig_length}||=0;
  $$settings{cov} ||= 0;
	usage($settings) if(@ARGV < 1 || $$settings{help});

  my %seq;
  for(@ARGV){
    my $seq=readMfa($_,$settings);
    $seq=filterSeqs($seq,$settings);
    %seq=(%seq,%$seq);
  }

  while(my($id,$sequence)=each(%seq)){
    $sequence=~s/(.{80})/$1\n/g;
    $sequence=~s/^\s+|\s+//g;     # trim the sequence's whitespace
    print ">$id\n$sequence\n";    # but put a newline after the sequence
  }

  return 0;
}

sub filterSeqs{
  my($seq,$settings)=@_;
  my %seq;
  my $i=0;
  while(my($id,$sequence)=each(%$seq)){
    # Remove short contigs
    if(length($sequence) < $$settings{assembly_min_contig_length}){
      logmsg "Removing short contig $id" if($$settings{verbose});
      next;
    }

    # Parse the ID for coverage: usually in the format of cov_\d+\.\d+
    # Remove the contig if the coverage is too low.
    my $cov;
    if($id=~/cov_(\d+(\.\d+)?)/){
      $cov=$1;
    }
    if(defined $cov && $cov < $$settings{cov}){
      logmsg "Removing low coverage contig $id with coverage $cov" if($$settings{verbose});
      next;
    }
    
    # Get rid of uncomplex contigs:
    #   Contigs of a K-mer repeat
    my $uncomplexkmer="";
    for(my $k=1;$k<20;$k++){
      my $kmer=substr($sequence,0,$k);
      # If we see a repeat of the kmer from start to finish, then it's uncomplex.
      # There is a {0,$k} to account for a truncated kmer at the end of the sequence.
      if($sequence=~/^($kmer)+(.{0,$k})$/i){
        $uncomplexkmer=$kmer;
        last;
      }
    }
    if($uncomplexkmer){
      logmsg "Removing uncomplex contig $id with kmer $uncomplexkmer";
      next;
    }

    if($$settings{rename}){
      my $newid = sprintf("contig%06d",++$i);
      $id=$newid;
      logmsg "Renaming $id to $newid" if($$settings{verbose});
    }

    $seq{$id}=$sequence;
  }
  return \%seq;
}

# Modified from AKUtils
sub readMfa {
	my ($mfa_file) = @_;
	open(FH, '<', $mfa_file) or die("Could not open file $mfa_file for reading: ".$!);
	my %seqs;
	my ($cur_seq_hdr, $seq);
	while (<FH>) {
		if (/^\>\s*(.+)/) {
			$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
			undef $seq;
			$cur_seq_hdr = (split /\s/, $1)[0];
		} else {
			chomp; 
      s/\s//g;
      $seq .= $_;
		}
	}
	close FH;
	$seqs{$cur_seq_hdr} = $seq if $cur_seq_hdr;
	die("Error: No sequences found in $mfa_file") unless %seqs;
	return \%seqs;
}

sub usage{
  my ($settings)=@_;
  print "Filters sequences from a fasta file that do not meet certain criteria and prints them to stdout
  Usage: $0 file.fasta [file2.fasta...] > out.fasta
  --length  $$settings{assembly_min_contig_length}   Minimum size of a contig
  --cov     $$settings{cov}     Minimum coverage of a contig - detected by searching deflines for cov_f where f is a floating point number
  --rename        rename contigs to 'contigXXXXXX' where XXXXXX is a number
  --verbose       
  \n";
  exit(0);
}

