#!/usr/bin/env perl

# extracts a sequence from a fasta file
# Good especially for when you don't feel like making a database

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/fileparse/;

my @fastaExt=qw(.fasta .ffn .fas .fna .faa .mfa);
my @gbkExt=qw(.gbk .gb .gbf .embl);

my $settings={
  errorcode_good=>0,
  errorcode_notFound=>1,
};

sub logmsg { local $0=fileparse $0; print STDERR "$0: @_\n";}
local $SIG{'__DIE__'} = sub {local $0=fileparse $0; my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };
exit(main());

sub main{
  die usage() if(@ARGV<1);
  GetOptions($settings,qw(search=s database=s insensitive! help));
  die usage() if($$settings{help});
  my $search=$$settings{search} || <>;
    chomp($search);
  my $db=$$settings{database} or die "Error: need database\n". usage();
  my @db=split(/\s*,\s*/,$db);

  my $fasta="";
  for my $d(@db){
    my($name,$path,$ext)=fileparse($d,@fastaExt,@gbkExt);
    if(grep(/$ext/,@fastaExt)){
      $fasta=fastacmdWithGrep($search,$d,$settings);
    }elsif(grep(/$ext/,@gbkExt)){
      $fasta=fastacmdWithGrepGbk($search,$d,$settings);
    }else{
      die "ERROR: I do not understand the extension on $d";
    }
    last if($fasta);
  }
  print $fasta;

  return $$settings{errorcode_notFound} if($fasta=~/^\s*$/);
  return $$settings{errorcode_good};
}

sub fastacmdWithGrep{
  my($search,$db,$settings)=@_;
  my $fasta="";

  my $grepParam="-m 1 -n";
  $grepParam.=" -i" if($$settings{insensitive});
  my $command="grep $grepParam '$search' '$db' 2>&1|cut -f 1 -d \":\"";
  my $lineNumber=`$command`;
  return fastacmd($search,$db,$settings) if($? || $!);
  if(!$lineNumber){
    logmsg "WARNING: could not find search $search in $db";
    return $fasta;
  }
  
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

sub fastacmdWithGrepGbk{
  my($search,$db,$settings)=@_;
  my $gbk="";

  my $grepParam="-m 1";
  $grepParam.=" -i" if($$settings{insensitive});
  my $command="grep -n LOCUS '$db' | grep $grepParam '$search' 2>&1|cut -f 1 -d \":\"";
  my $lineNumber=`$command`;

  if(!$lineNumber){
    logmsg "WARNING: Could not find $search in $db";
    return $gbk;
  }
  
  my $i=1;
  open(DB,$db) or die "Could not open database $db because $!\n".usage();
  while(<DB>){
    $i++;
    next if($i<$lineNumber);
    last if($i-1 > $lineNumber && /^\s*LOCUS/);
    $gbk.=$_;
  }
  close DB;
  $gbk=~s/^\s*|\s*$//g;
  return $gbk;
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
  local $0=fileparse $0;
  "\n\n  Finds a sequence in a fasta file and prints it. Searches
  contig names if gbk; searches deflines if fasta.
  
  Usage: $0 -s search -d database
  -s  '' search term in the defline.
         search can also be from stdin
  -d  '' fasta or gbk file to search
         comma-separated databases if multiples
  -i     case-insensitive
";
}
