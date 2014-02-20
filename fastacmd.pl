#!/usr/bin/env perl

# extracts a sequence from a fasta file
# Good especially for when you don't feel like making a database

use strict;
use warnings;
use Getopt::Long;

my $settings={
  errorcode_good=>0,
  errorcode_notFound=>1,
};

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
exit(main());

sub main{
  die usage() if(@ARGV<1);
  GetOptions($settings,qw(search=s database=s insensitive!));
  my $search=$$settings{search} || <>;
    chomp($search);
  my $db=$$settings{database} or die "Error: need database\n". usage();
  my @db=split(/\s*,\s*/,$db);

  my $fasta;
  for my $d(@db){
    $fasta=fastacmdWithGrep($search,$d,$settings);
    last if($fasta);
  }
  print $fasta;

  return $$settings{errorcode_notFound} if($fasta=~/^\s*$/);
  return $$settings{errorcode_good};
}

sub fastacmdWithGrep{
  my($search,$db,$settings)=@_;
  my $fasta;

  my $grepParam="-m 1 -n";
  $grepParam.=" -i" if($$settings{insensitive});
  my $command="grep $grepParam '$search' '$db' 2>&1|cut -f 1 -d \":\"";
  my $lineNumber=`$command`;
  return fastacmd($search,$db,$settings) if($? || $!);
  return "" if(!$lineNumber);
  
  my $i=0;
  open(DB,$db) or die "Could not open database $db because $!\n".usage();
  ENTRY:
  while(<DB>){
    $i++;
    next if($i<$lineNumber);
    $fasta.=$_;
    while(<DB>){
      last ENTRY if(/>/);
      $fasta.=$_;
    }
  }
  return $fasta;
}

sub fastacmd{
  my($search,$db,$settings)=@_;
  print "trying line by line";
  $search=lc($search) if($$settings{insensitive});
  my $fasta="";
  open(DB,$db) or die "Could not open database $db because $!\n".usage();
  ENTRY:
  while(<DB>){
    $_=lc($_) if($$settings{insensitive});
    if(/$search/ || ($$settings{insensitive} && /$search/i)){
      $fasta.=$_;
      while(<DB>){
        last ENTRY if(/>/);
        $fasta.=$_;
      }
    }
  }
  close DB;
  return $fasta;
}

sub usage{
  "
Finds a sequence in a fasta file and prints it  
perl $0 -s search -d database
  -s search term in the defline.
    search can also be from stdin
  -d fasta file to search
    comma-separated databases if multiples
  -i case-insensitive
";
}
