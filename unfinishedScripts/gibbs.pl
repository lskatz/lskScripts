#!/usr/bin/env perl
#performs gibbs sampling
#usage: perl lab5.pl [seq]
#100% done

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw/basename/;
use Data::Dumper;

local $0=basename $0;
sub logmsg{print STDERR "$0: @_\n";}
my $file_name = $ARGV[0] || "seq";

#------------------------------
my @nt = ("A","T","G","C");
my $J = @nt; #total number of nucleotides
my $W = 6; #motif length
my $i=0;
#my $L = 24; # sequence length - 1
my $sequence_length=25;
my $total_reps = 10000; #total number of times to go through list of sequences
my $check_freq=1000; #how often to compare Fs to see if it should quit
#------------------------------

print "Total iterations: $total_reps, checking every $check_freq iterations.\n\n";

# Initialize the matrix
my @sequence;
my $i;
my %nucleotide_counter;
my @motif_start = ();
open (my $fh, "$file_name") || die "Couldn't open input file, $!\n";
while (<$fh>) {
    chomp;
    $sequence[$i] = $_;
   
    # randomly assign positions of motif start
    push (@motif_start, int(rand($L - $W)));

    $pos = 0;

    #motif positional counts and background counts

    for $j ( 0 .. $sequence_length-1 ) {

	$nucleotide = substr($sequence[$i], $j, 1);
        $nucleotide_counter{$nucleotide}++;

	if($j >= $motif_start[$i] && $j < ($motif_start[$i] + 6)){
	    $cMatrix{$nucleotide}{$pos}++;
	    $pos++;
	}

	else{
	    $BG_count{$nucleotide}++;
	    #$num_nt++;
	}
    }

    $i++;
}
close $fh;

#total number of sequences
$num_sequences=$i;

#MY CODE HERE

#pseudocounts
$B=0;
$coeff=1.1;
foreach $key (keys %c){
	$PsC{$key}=($coeff-1)*$c{$key}; #pseudocount of that nucleotide
	$B+=$PsC{$key};
}
$B+=1;

#calculate q matrix from motif hash matrix
for($i=0;$i<$W;$i++){
    for($j=0;$j<$J;$j++){
	$qMatrix{$nt[$j]}{$i} = ($PsC{$nt[$j]}+$cMatrix{$nt[$j]}{$i}) / ($B+$num_sequences);
    }
}
$hashSum=hashSum(\%BG_count);
for($j=0;$j<$J;$j++){
	$BG_q{$nt[$j]}=($BG_count{$nt[$j]}+$PsC{$nt[$j]})/($hashSum+$B);
}

#calculate F
$prev_score=-1;
$F=0;
for($i=0;$i<$W;$i++){
	for($j=0;$j<$J;$j++){
		$F+=$cMatrix{$nt[$j]}{$i}*log($qMatrix{$nt[$j]}{$i}/$BG_q{$nt[$j]});
	}
}

print "Iteration number :: F score\n";
#go through iterations
for($l=0;$l<$total_reps;$l++){
    #go through each sequence
    for($m=0;$m<$num_sequences;$m++){
        #take out the current sequence's motif from the cMatrix
	$num_sequences--;
        for($j=0;$j<$W;$j++){
            $nucleotide=substr($sequence[$m],$a[$m]+$j,1);
            $cMatrix{$nucleotide}{$j}--;
            $BG_count{$nucleotide}++;
        }
        #recalculate qMatrix
	for($i=0;$i<$W;$i++){
	    for($j=0;$j<$J;$j++){
		$qMatrix{$nt[$j]}{$i} = ($PsC{$nt[$j]}+$cMatrix{$nt[$j]}{$i}) / ($B+$num_sequences);
	    }
	}
        #calculate A values-----------
	#loop through each possible motif through the sequence
	$num=$L-$W+1;
	$total=0; #keep track of total so that A can be normalized
	for($k=0;$k<$num;$k++){
	    #initilize this A value
	    $A[$k]=1;
	    #make up the motif and add it into the cMatrix
	    for($i=0;$i<$W;$i++){
		$nucleotide=substr($sequence[$m],$i+$k,1);
		$cMatrix{$nucleotide}{$i}++;
		$BG_count{$nucleotide}--;
	    }
	    #recalculate qMatrix
	    for($i=0;$i<$W;$i++){
		for($j=0;$j<$J;$j++){
		    $qMatrix{$nt[$j]}{$i} = ($PsC{$nt[$j]}+$cMatrix{$nt[$j]}{$i}) / ($B+$num_sequences);
		}
	    }
	    #loop through residues and through positions
	    for($i=0;$i<$W;$i++){
		for($j=0;$j<$J;$j++){
		    $A[$k]*=($qMatrix{$nt[$j]}{$i}/$BG_count{$nt[$j]});
		}
	    }
	    #keep track of the total in A so that it can be divided out later
	    $total+=$A[$k];
	    #remove motif from CMatrix from previous iteration
	    for($i=0;$i<$W;$i++){
		$nucleotide=substr($sequence[$m],$i+$k,1);
		$cMatrix{$nucleotide}{$i}--;
		$BG_count{$nucleotide}++;
	    }
	} #end k loop for adding and subtracting motifs
	#normalize A
	$num_A=@A;
	for($k=0;$k<$num_A;$k++){
	    $A[$k] = ($A[$k]) / $total;
	}
	#make an actual scoring array with probabilities from zero to one
	$Ax[0]=$A[0];
	for($k=1;$k<$num_A;$k++){
	    $Ax[$k]=$Ax[$k-1]+$A[$k];
	}
	$Ax[$num_A]=1;
	#find out which numbers the random number lies in between
	$old_a[$m]=$a[$m];
	$rand=rand(0 . 1);
	for($k=0;$k<$num_A;$k++){
	    if($rand>=$Ax[$k] && $rand<$Ax[$k+1]){
		$a[$m]=$k;
		last;
	    }
	}
	$new_motif=substr($sequence[$m],$a[$m],$W);
	#insert the new motif
	$num_sequences++;
	for($i=0;$i<$W;$i++){
	    $nucleotide=substr($new_motif,$i,1);
	    $cMatrix{$nucleotide}{$i}++;
	    $BG_count{$nucleotide}--;
	}
        #recalculate qMatrix, F
        for($i=0;$i<$W;$i++){
	    for($j=0;$j<$J;$j++){
		$qMatrix{$nt[$j]}{$i} = ($PsC{$nt[$j]}+$cMatrix{$nt[$j]}{$i}) / ($B+$num_sequences);
	    }
	}
	$hashSum=hashSum(\%BG_count);
	for($j=0;$j<$J;$j++){
	  $BG_q{$nt[$j]}=($BG_count{$nt[$j]}+$PsC{$nt[$j]})/($hashSum+$B);
	}
	#F
	$prev_F=$F;
	$F=0;
	for($i=0;$i<$W;$i++){
	  for($j=0;$j<$J;$j++){
	    $F+=$cMatrix{$nt[$j]}{$i}*log($qMatrix{$nt[$j]}{$i}/$BG_q{$nt[$j]});
	  }
	}
        #compare with previous F
	#if the previous F is bigger, then keep that F and keep the old motifs
	if($prev_F>$F){
	  $F=$prev_F;
	  $a[$m]=$old_a[$m];
	  #subtract off current motif and add back on the old one
	  for($i=0;$i<$W;$i++){
	    $nucleotide=substr($new_motif,$i,1);
	    $cMatrix{$nucleotide}{$i}--;
	    $BG_count{$nucleotide}++;
	  }
	  for($i=0;$i<$W;$i++){
	    $nucleotide=substr($sequence[$m],$a[$m]+$i,1);
	    $cMatrix{$nucleotide}{$i}++;
	    $BG_count{$nucleotide}--;
	  }
	}
	#if the current F is bigger, then keep the current F and the current motifs (do nothing)
    }
    #check F every now and then to see if it is the same F

    if($l%$check_freq==0){
	print "$l :: $F\n";
	if($F==$prev_score){
	  last;
	}
	$prev_score=$F;
    }
}
print "\n";

#---- print the results ----
for ($i = 0; $i <= $#nt; $i++){
    print"$nt[$i] \t";
    for($l = 0; $l < $W; $l++){
    	if(!$number[$l]){
		$number[$l]=0;
	}
        #Your script HERE
        printf("%.3f   ",$qMatrix{$nt[$i]}{$l});
	if($number[$l]<$qMatrix{$nt[$i]}{$l}){
		$number[$l]=$qMatrix{$nt[$i]}{$l};
		$residue[$l]=$nt[$i];
	}
	$ttl[$l]+=$qMatrix{$nt[$i]}{$l};
    }print"\n";
}
#print totals row
print "ttl \t";
for ($i=0;$i<$W;$i++){
	$total=0;
	for($j=0;$j<$J;$j++){
		$total+=$number[$j];
	}
	printf("%.3f   ",$ttl[$i]);
}
print "\n";
for($l = 0; $l < $W; $l++){
	print "$residue[$l] ";
}
print "\n\n";

#print time that it takes
$time=time()-$time;
$minutes = $time/60;
$speed=$total_reps/$minutes;
print "Total Score: $F\n";
printf ("%8.2f iterations per minute\n",$speed);

sub hashSum{
  my($hash)=@_;
  $total=0;
  foreach $key (keys %$hash){
    $total+=$$hash{$key};
  }
  return $total;
}
