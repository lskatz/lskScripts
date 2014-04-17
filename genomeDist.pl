#!/usr/bin/env perl
# Finds the distance between any two assemblies or raw reads
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Spec;

my $start=time;

sub logmsg{ $|++; print STDERR ((time - $start)."\t"); print STDERR "@_\n"; $|--;}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help averages quiet method=s coverage=i kmerlength=i tempdir=s numcpus=i downsample=s)) or die $!;
  die usage() if(!@ARGV || $$settings{help});
  my @asm=@ARGV;
  $$settings{method}||="mummer";
  $$settings{method}=lc($$settings{method});
  $$settings{coverage}||=2; $$settings{coverage}=1 if($$settings{coverage}<1);
  $$settings{kmerlength}||=18;
  $$settings{tempdir}||="tmp";
  $$settings{numcpus}||=1;
  $$settings{downsample}||=0;

  if(!-d $$settings{tempdir}){
    logmsg "WARNING: I could not find $$settings{tempdir}; creating it now.";
    mkdir $$settings{tempdir};
    die $! if $?;
  }

  # GENOME DISTANCE METHODS
  if($$settings{method} eq 'mummer'){
    my $pdist=mummer(\@asm,$settings);
    
    print join("\t",".",@asm)."\n";
    for(my $i=0;$i<@asm;$i++){
      my $asm1=$asm[$i];
      print "$asm1\t";
      for(my $j=0;$j<@asm;$j++){
        my $asm2=$asm[$j];
        my $num=$$pdist{$asm1}{$asm2};
        $num=($$pdist{$asm1}{$asm2} + $$pdist{$asm2}{$asm1})/2 if($$settings{averages});
        print "$num\t";
      }
      print "\n";
    }
  }
  elsif($$settings{method} eq 'jaccard'){
    jaccardDistance(\@asm,$settings);
  }
  else {
    die "I do not know how to perform method $$settings{method}";
  }

  return 0;
}

sub jaccardDistance{
  my($genome,$settings)=@_;
  my %jDist;

  downsample($genome,$settings) if($$settings{downsample});

  for(my $i=0;$i<@$genome-1;$i++){
    logmsg "Comparing $$genome[$i]";
    my %k1=kmerCountJellyfish($$genome[$i],$settings);
    for(my $j=$i+1;$j<@$genome;$j++){
      my %k2=kmerCountJellyfish($$genome[$j],$settings);
      my $jDist=jDist(\%k1,\%k2,$settings);
      print join("\t",@$genome[$i,$j],$jDist)."\n";
    }
  }
}

# downsample each genome reads file into a temp directory
sub downsample{
  my($genome,$settings)=@_;
  logmsg "Downsampling the reads by $$settings{downsample}";
  for my $g(@$genome){
    my $newG="$$settings{tempdir}/".fileparse($g,qw(.fastq .fastq.gz)).".fastq";
    system("run_assembly_removeDuplicateReads.pl --downsample $$settings{downsample} $g > $newG");
    die "ERROR: problem with CGP run_assembly_removeDuplicateReads.pl" if $?;
    $g=$newG;
  }
}

sub jDist{
  my($k1,$k2,$settings)=@_;
  my $minKCoverage=$$settings{coverage};

  #logmsg "Finding intersection and union of kmers";
  my %kmerSet=kmerSets($k1,$k2,$settings);

  my $jDist=1-($kmerSet{intersection} / $kmerSet{union});

  #logmsg "$jDist=1-($kmerSet{intersection} / $kmerSet{union})";
  return $jDist;
}

sub kmerSets{
  my($k1,$k2,$settings)=@_;
  
  my($intersectionCount,%union);
  $intersectionCount=0;

  # Find uniq kmers in the first set of kmers.
  # Also find the union.
  for my $kmer(keys(%$k1)){
    $intersectionCount++ if(!$$k2{$kmer});
    $union{$kmer}=1;
  }

  # Find uniq kmers in the second set of kmers.
  # Also find the union.
  for my $kmer(keys(%$k2)){
    $intersectionCount++ if(!$$k1{$kmer});
    $union{$kmer}=1;
  }

  my $unionCount=scalar(keys(%union));

  return (intersection=>$intersectionCount,union=>$unionCount);
}

sub kmerCountJellyfish{
  my($genome,$settings)=@_;
  my $minKCoverage=$$settings{coverage};
  my($name,$path,$suffix)=fileparse($genome,qw(.fastq.gz .fastq .gz));
  my $tmpName="$$settings{tempdir}/$name.fastq";

  if($suffix =~/\.gz$/){
    my $newName="$$settings{tempdir}/$name.fastq";
    logmsg "Uncompressing $genome to $newName";
    system("gunzip -c $genome > $newName"); die if $?;
    $genome=$newName;
  }
  
  # use jellyfish to count kmers
  my $rand=int(rand(999999999)); # avoid clashing
  my $outprefix="$$settings{tempdir}/$rand.mer_counts";
  my $jfDb="$$settings{tempdir}/merged.$rand.jf";
  my $kmerTsv="$$settings{tempdir}/jf.$rand.tsv";
  system("jellyfish count -s 10000000 -m $$settings{kmerlength} -o $outprefix -t $$settings{numcpus} $genome");
  die "Error: problem with jellyfish" if $?;
  my @mer=glob("${outprefix}_*");
  my $merStr=join(" ",@mer);
  if(@mer >1){
    system("jellyfish merge $merStr -o $jfDb");
  } else {
    system("cp -v $merStr $jfDb 1>&2");
  }
  die if $?;
  system("rm ${outprefix}_*"); die if $?;
  system("jellyfish dump -L $minKCoverage --column --tab -o $kmerTsv $jfDb");
  die if $?;
  
  # load kmers to memory
  my %kmer;
  open(TSV,$kmerTsv) or die "ERROR: Could not open $kmerTsv: $!";
  while(<TSV>){
    chomp;
    my @F=split /\t/;
  
    #next if($F[1] < $minKCoverage); # remove kmers with low depth
    $kmer{$F[0]}=$F[1];
  }
  close TSV;

  unlink $tmpName if($tmpName && -e $tmpName);
  return %kmer;
}


sub mummer{
  my($asm,$settings)=@_;

  # Make symbolic links to all the genomes in the temp directory so that random files aren't spewed everywhere.
  my (@asm,@asmName);
  for(my $i=0;$i<@$asm;$i++){
    my $asmName=fileparse($$asm[$i]);
    my $target="$$settings{tempdir}/$asmName";
    symlink File::Spec->rel2abs($$asm[$i]), $target or warn "Warning while trying to make a symbolic link ($target):$!";
    push(@asm,$target);
    push(@asmName,$asmName);
  }

  # calculate distances with mummer/nucmer
  my %pdist=();
  for(my $i=0;$i<@asm;$i++){
    my $asm1=$asm[$i];
    for(my $j=0;$j<@asm;$j++){
      my $asm2=$asm[$j];
      my $prefix="$$settings{tempdir}/".join("_",$asmName[$i],$asmName[$j]);
      if(!-e "$prefix.snps"){
        my $nucErr="$$settings{tempdir}/nucmer.err";
        logmsg "Running on $prefix" unless($$settings{quiet});
        my $command="nucmer --prefix $prefix $asm1 $asm2 2>$nucErr";
        system($command); die "Error running nucmer with command:\n  $command\n".`cat $nucErr` if $?;
        system("show-snps -Clr $prefix.delta > $prefix.snps 2>/dev/null"); die if $?;
      }
      my $numSnps=countSnps("$prefix.snps",$settings);
      #$pdist{$asm1}{$asm2}=$numSnps;
      $pdist{$asmName[$i]}{$asmName[$j]}=$numSnps;
      #$pdist{$asm2}{$asm1}=$numSnps;
    }
  }
  return \%pdist;
}

sub countSnps{
  my($snpsFile,$settings)=@_;
  my $num=0;
  open(SNP,$snpsFile) or die "Could not open $snpsFile: $!";
  while(<SNP>){
    $num++;
  }
  close SNP;
  $num-=5; # header
  die "internal parsing error" if ($num<0);
  return $num;
}

sub usage{
  "Finds the p-distance between two assemblies using mummer. With more genomes, it creates a table.
  Usage: $0 assembly.fasta assembly2.fasta [assembly3.fasta ...]
  -a to note averages. Switching the subject and query can reveal some artifacts in the algorithm.
  -q for minimal stdout
  -m method.  Can be mummer (default) or jaccard
    Mummer: uses mummer to discover SNPs and counts the total number
    Jaccard: (kmer method) counts 18-mers and calculates 1 - (intersection/union)
  -c minimum kmer coverage. Default: 2
  -k kmer length. Default: 18
  -t tempdir Default: tmp
  -n numcpus Default: 1
  --downsample 0 Downsample the reads to a frequency, using CG-Pipeline. Useful for removing extraneous data and saving on computation time. Values are between 0 and 1; 0 means no downsampling. 1 will remove duplicate reads only.
  "
}

