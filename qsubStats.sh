#!/bin/bash

# Author: Lee Katz
# Figures out some quick metrics from SGE

# Take a snapshot of qstat.
QSTAT=$(qstat -u '*')
if [ $? -gt 0 ]; then
  echo "ERROR with qstat" >&2
  exit 1;
fi

QSTAT=$(echo "$QSTAT"|tail -n +3) # qstat, minus the header

# How many of the cluster's slots I'm taking
# echo "$QSTAT" |tail -n +3| perl -MData::Dumper -e 'while(<>){s/^\s+|\s+$//g; @F=split /\s+/; next if($F[4] ne 'r'); $slots=$F[8]; if($F[3] eq $ENV{USER}){$mine+=$slots;} $total+=$slots; } print "$mine out of $total\n";'

# who is the current hog
echo "$QSTAT" | perl -lane 'BEGIN{print "USER\tSLOTS";} next if(!$F[3] || !$F[8]); next if($F[4] ne "r"); $slot{$F[3]}+=$F[8]; END{@user=sort{$slot{$b}<=>$slot{$a}} keys(%slot); print "$_\t$slot{$_}" for @user;}' | column -t
