#!/usr/bin/env perl
# Finds the distance between any two assemblies or raw reads
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Spec;
use File::Temp qw/tempdir/;

my $start=time;

sub logmsg{ $|++; print STDERR ((time - $start)."\t"); print STDERR "@_\n"; $|--;}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help averages quiet method=s coverage=i kmerlength=i tempdir=s numcpus=i downsample=s)) or die $!;
  die usage() if(!@ARGV || $$settings{help});
  my($ref,@query)=@ARGV;
  $$settings{method}||="mummer";
  $$settings{method}=lc($$settings{method});
  $$settings{coverage}||=2; $$settings{coverage}=1 if($$settings{coverage}<1);
  $$settings{kmerlength}||=18;
  $$settings{tempdir}||=tempdir(TEMPLATE=>"genomeDistXXXXXX",CLEANUP=>1,TMPDIR=>1);
  $$settings{numcpus}||=1;
  $$settings{downsample}||=0;

  if(!-d $$settings{tempdir}){
    logmsg "WARNING: I could not find $$settings{tempdir}; creating it now.";
    mkdir $$settings{tempdir};
    die $! if $?;
  }
  logmsg "Temp dir is $$settings{tempdir}";

  # GENOME DISTANCE METHODS
  if($$settings{method} eq 'mummer'){
    my $pdist=mummer($ref,\@query,$settings);
    ...;
    
#    # TODO move this to the mummer sub
#    print join("\t",".",$ref,@query)."\n";
#    for(my $i=0;$i<@asm;$i++){
#      my $asm1=$asm[$i];
#      print "$asm1\t";
#      for(my $j=0;$j<@asm;$j++){
#        my $asm2=$asm[$j];
#        my $num=$$pdist{$asm1}{$asm2};
#        $num=($$pdist{$asm1}{$asm2} + $$pdist{$asm2}{$asm1})/2 if($$settings{averages});
#        print "$num\t";
#      }
#      print "\n";
#    }
  }
  elsif($$settings{method} eq 'jaccard'){
    jaccardDistance($ref,\@query,$settings);
  }
  elsif($$settings{method} eq 'jaccardkanalyze'){
    jaccardKAnalyze($ref,\@query,$settings);
  }
  else {
    die "I do not know how to perform method $$settings{method}";
  }

  return 0;
}

sub jaccardDistance{
  my($ref,$query,$settings)=@_;
  my %jDist;

  if($$settings{downsample}){
    my @genome=downsample([$ref,@$query],$settings);
    $ref=shift(@genome);
    $query=\@genome;
  }
  
  my %refK=kmerCountJellyfish($ref,$settings);
  for(my $j=0;$j<@$query;$j++){
    my %queryK=kmerCountJellyfish($$query[$j],$settings);
    my $jDist=jDist(\%refK,\%queryK,$settings);
    print join("\t",$ref,$$query[$j],$jDist)."\n";
  }
}

sub jaccardKAnalyze{
  my($ref,$query,$settings)=@_;
  my %jDist;

  if($$settings{downsample}){
    my @genome=downsample([$ref,@$query],$settings);
    $ref=shift(@genome);
    $query=\@genome;
  }
  
  ...;
#  for(my $i=0;$i<@$genome-1;$i++){
#    logmsg "Comparing $$genome[$i]";
#    my $kmerFile1=kmerCountKAnalyze($$genome[$i],$settings);
#    for(my $j=$i+1;$j<@$genome;$j++){
#      my $kmerFile2=kmerCountKAnalyze($$genome[$j],$settings);
#      my $jDist=jDistSortedFile($kmerFile1,$kmerFile2,$settings);
#      print join("\t",@$genome[$i,$j],$jDist)."\n";
#    }
#  }
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
  return @$genome if(wantarray);
  return $genome;
}

sub jDistSortedFile{
  my($kmerFile1,$kmerFile2,$settings)=@_;
  my($matchCount,$uniqueCount)=(0,0);
  open(KC1,$kmerFile1) or die "ERROR: could not open $kmerFile1:$!";
  open(KC2,$kmerFile2) or die "ERROR: could not open $kmerFile2:$!";
  while (my $kc1=<KC1>){
    chomp($kc1);
    my($kmer1,$kmercount1)=split(/\t/,$kc1);
    next if($kmercount1 < $$settings{coverage});

    # Get the second file's kmer 
    my($kmer2,$kmercount2)=("",0);
    do{
      my $kc2=<KC2>; chomp($kc2);
      last if(!$kc2); # Quit out if this is the end of the file
      ($kmer2,$kmercount2)=split(/\t/,$kc2);
    } while($kmercount2 < $$settings{coverage});
    last if(!$kmer2);   # Stop reading KC1 if there is nothing to compare against

    # Do the test one time for the whole iteration
    my $cmp=$kmer1 cmp $kmer2;

    if($cmp==0){
      $matchCount++;
      #print "$kmer1\t$kmer2\t$matchCount\t$uniqueCount\n";
      undef($kmer2); # marks that a new one has to be obtained
      next;
    }
    ############
    # MISMATCH #
    ############

    # If the left comes before alphabetically, it has to catch up. 
    # Advance to the next kmer1.
    while($cmp==-1){
      $uniqueCount++;   # Count the current mismatch
      #print "$kmer1\t$kmer2\t$matchCount\t$uniqueCount\n";
      do{               # Move onto the next kmer1
        $kc1=<KC1>;
        chomp($kc1);
        last if(!$kc1); # quit out if this is the end of file
        ($kmer1,$kmercount1)=split(/\t/,$kc1);
      } while($kmercount1 < $$settings{coverage});  # Keep moving if the coverage isn't good yet
      $cmp=$kmer1 cmp $kmer2;
    }

    # If the left kmer comes after alphabetically, the right has to catch up.
    # Advance to the next kmer2.
    while($cmp==1){
      $uniqueCount++;
      #print "$kmer1\t$kmer2\t$matchCount\t$uniqueCount\n";
      do{
        my $kc2=<KC2>;
        chomp($kc2);
        last if(!$kc2); # quit out if this is the end of file
        ($kmer2,$kmercount2)=split(/\t/,$kc2);
      } while($kmercount2 < $$settings{coverage});
      $cmp=$kmer1 cmp $kmer2;
    }

    # now that the lines are supposedly equal, seek back to the previous line
    seek(KC1,-length("$kmer1\t$kmercount1\n"),1);
    seek(KC2,-length("$kmer2\t$kmercount2\n"),1);
  }
  # Close out the files and count mismatches
  while(<KC1>){
    $uniqueCount++;
    #print "KC1\tKC1\t$matchCount\t$uniqueCount\n";
  }
  while(<KC2>){
    $uniqueCount++;
    #print "KC2\tKC2\t$matchCount\t$uniqueCount\n";
  }
  #logmsg "1-$uniqueCount/($matchCount+$uniqueCount)";
  my $jDist=1-$uniqueCount/($matchCount+$uniqueCount);

  return $jDist;
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

sub kmerCountKAnalyze{
  my($genome,$settings)=@_;
  my $minKCoverage=$$settings{coverage};
  my($name,$path,$suffix)=fileparse($genome,qw(.fastq.gz .fastq .gz));
  my $rand=int(rand(999999999)); # avoid clashing
  my $outfile="$$settings{tempdir}/$name.$rand.kc";
  my $informat="";
  if($suffix eq '.fastq.gz'){
    $informat="fastqgz";
  } elsif ($suffix eq '.fastq'){
    $informat="fastq";
  } else {
    die "I have no idea what to do with the extension on $genome ($suffix)";
  }

  system("count -d $$settings{numcpus} -l $$settings{numcpus} -f fastqgz -k 18 -g 6553600 -o $outfile $genome");
  die if $?;
  return $outfile;
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

  for($tmpName, $jfDb, $kmerTsv){
    unlink $_ if($_ && -e $_);
  }
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
  Usage: $0 reference.fasta query1.fasta [query2.fasta ...]
  -q for minimal stdout
  -m method.  Can be the following choices
    Mummer: (default): uses mummer to discover SNPs and counts the total number
    Jaccard: (kmer method) counts 18-mers and calculates 1 - (intersection/union)
    JaccardKAnalyze: (Same method as 'Jaccard' but uses KAnalyze)
  -c minimum kmer coverage. Default: 2
  -k kmer length. Default: 18
  -t tempdir If one is not given, then one will be made under /tmp/
  -n numcpus Default: 1
  --downsample 0 Downsample the reads to a frequency, using CG-Pipeline. Useful for removing extraneous data and saving on computation time. Values are between 0 and 1; 0 means no downsampling. 1 will remove duplicate reads only.
  "
}

