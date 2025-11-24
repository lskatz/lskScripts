#!/usr/bin/env perl
################################################################################
# lincodes_LK.pl - Define LINcodes from cgMLST profiles (Flat File Version)
#
# OVERVIEW:
# LINcodes (Life Identification Numbers) are hierarchical taxonomic identifiers
# based on genetic distance thresholds. This script assigns LINcodes to cgMLST
# (core genome Multi-Locus Sequence Typing) profiles using flat files instead
# of a database.
#
# HOW IT WORKS:
# 1. Reads allelic profiles from tab-separated files
# 2. Calculates pairwise genetic distances between profiles
# 3. Uses Prim's algorithm to order profiles for optimal assignment
# 4. Assigns hierarchical LINcodes based on configurable distance thresholds
# 5. Each position in a LINcode represents a different similarity threshold
#
# EXAMPLE:
# LINcode "0_1_2_0" means:
#   - Position 0 (value=0): Same major lineage
#   - Position 1 (value=1): Different at second threshold
#   - Position 2 (value=2): Different at third threshold  
#   - Position 3 (value=0): Same at finest resolution
#
# Written by Keith Jolley
# Based on code by Melanie Hennart (https://gitlab.pasteur.fr/BEBP/LINcoding).
# Copyright (c) 2022-2025, University of Oxford
# E-mail: keith.jolley@biology.ox.ac.uk
# Modified for standalone flat file use: 2025
#
# This file is part of Bacterial Isolate Genome Sequence Database (BIGSdb).
#
# BIGSdb is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# BIGSdb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BIGSdb.  If not, see <http://www.gnu.org/licenses/>.
#
# Version: 20250220
################################################################################
use strict;
use warnings;
use 5.010;

# Core modules
use Getopt::Long qw(:config no_ignore_case);
use Digest::MD5;
use File::Temp;
use File::Path qw(make_path);
use File::Spec;

# Numerical computing
use PDL;
use PDL::IO::FastRaw;
use File::Map;

# Global variables
my %opts;
my $script;

# Enable PDL big data mode
# PDL 2.067 2022-01-13 changelog: "remove $PDL::BIGPDL check for 1GB data" 
$PDL::BIGPDL = 1 if !$PDL::BIGPDL;

# Cleanup handler for temp files
END {
	if (defined $script && $script->{'config'}->{'secure_tmp_dir'}) {
		my $tmp_dir = $script->{'config'}->{'secure_tmp_dir'};
		if (-d $tmp_dir) {
			opendir(my $dh, $tmp_dir);
			my @files = grep { !/^\./ } readdir($dh);
			closedir $dh;
			unlink File::Spec->catfile($tmp_dir, $_) foreach @files;
		}
	}
}

# Parse command line options
GetOptions(
	'dir=s'          => \$opts{'dir'},
	'debug'          => \$opts{'debug'},
	'batch_size=i'   => \$opts{'batch_size'},
	'create'         => \$opts{'create'},
	'input_profiles=s' => \$opts{'input_profiles'},
	'log=s'          => \$opts{'log'},
	'missing=i'      => \$opts{'missing'},
	'mmap'           => \$opts{'mmap'},
	'q|quiet'        => \$opts{'quiet'},
	's|scheme_id=i'  => \$opts{'scheme_id'},
	'x|min=s'        => \$opts{'x'},
	'y|max=s'        => \$opts{'y'},
	'help'           => \$opts{'help'},
);

if ( $opts{'help'} ) {
	show_help();
	exit;
}

if ( $opts{'create'} ) {
	if ( !$opts{'dir'} || !$opts{'scheme_id'} ) {
		say 'Usage: lincodes_LK.pl --create --dir [DIRECTORY] --scheme [SCHEME_ID]';
		say 'Use --help for more information.';
		exit;
	}
	create_example_schema();
	say "Example schema created in $opts{'dir'}/";
	exit;
}

if ( !$opts{'dir'} || !$opts{'scheme_id'} ) {
	say 'Usage: lincodes_LK.pl --dir [DIRECTORY] --scheme [SCHEME_ID]';
	say 'Use --help for more information.';
	exit;
}

# Validate and create directory if needed
if ( !-d $opts{'dir'} ) {
	make_path($opts{'dir'}) or die "Cannot create directory $opts{'dir'}: $!\n";
}

# Set defaults for optional parameters
$opts{'batch_size'} //= 10_000;  # Maximum profiles to process per batch
$opts{'missing'}    //= 0;        # Maximum missing loci allowed in profiles

# Initialize script object with file paths and data structures
$script = initialize_script();

# Concatenate input profiles if provided
if ( $opts{'input_profiles'} ) {
	concatenate_profiles();
}

# Verify that required files exist and are properly formatted
check_db();

# Prevent multiple instances from running simultaneously on same data
check_if_script_already_running();

# Initialize debug log file if requested
if ( $opts{'log'} ) {
	initiate_log_file( $opts{'log'} );
}

# Run main LINcode assignment loop
main();

# Clean up script object
undef $script;

# Remove lock file to allow future runs
remove_lock_file();

sub main {
	# Enable autoflush for real-time output
	local $| = 1;
	
	# Load any previously assigned LINcodes from file
	my $lincodes = get_lincode_definitions();
	if ( !@{ $lincodes->{'profile_ids'} } ) {
		say 'No LINcodes yet defined.' if !$opts{'quiet'};
	}
	
	# Main processing loop: continue until all profiles have LINcodes
	my $profiles_to_assign;
	do {
		# Get batch of profiles that don't have LINcodes yet
		$profiles_to_assign = get_profiles_without_lincodes();
		if ( !@$profiles_to_assign ) {
			say 'All profiles assigned.' if !$opts{'quiet'};
			return;
		}
		
		# Determine optimal order for assignment using Prim's algorithm
		# (skip if only one profile to assign)
		my $profiles = @$profiles_to_assign == 1 ? $profiles_to_assign : get_prim_order($profiles_to_assign);
		
		# If LINcodes already exist, adjust order to start with profile closest to existing ones
		if ( @{ $lincodes->{'profile_ids'} } ) {
			$profiles = adjust_prim_order( $lincodes->{'profile_ids'}, $lincodes->{'profiles'}, $profiles );
		}
		
		# Assign LINcodes to profiles in optimized order
		assign_lincodes($profiles);
		
		# Reload LINcode definitions to include newly assigned ones
		$lincodes = get_lincode_definitions();
	} while @$profiles_to_assign;
	return;
}

sub adjust_prim_order {
	my ( $assigned_profile_ids, $assigned_profiles, $new_profiles ) = @_;
	
	# Get locus information for distance calculations
	my $loci             = get_scheme_loci();
	my $locus_count      = @$loci;
	my $closest_distance = 100;  # Start with maximum distance (100%)
	my $closest_profile_index;    # Index of closest new profile
	my $index   = 0;
	my %missing = ( N => 0 );     # Map missing alleles (N) to 0 for PDL
	
	print 'Adjusting PRIM order ...' if !$opts{'quiet'};
	print "\n"                       if @$new_profiles >= 500;
	my $start_time = time;

	# Find which new profile is closest to any already-assigned profile
	foreach my $profile_id (@$new_profiles) {
		# Load profile data from file
		my $profile_array = get_profile($profile_id);
		next unless $profile_array;
		
		# Convert missing alleles (N) to 0 for comparison
		$_ = $missing{$_} // $_ foreach @$profile_array;
		my $profile = pdl($profile_array);
		
		# Compare this new profile against all assigned profiles
		for my $i ( 0 .. @$assigned_profile_ids - 1 ) {
			my $assigned_profile = $assigned_profiles->slice(",($i)");
			
			# Count allelic differences (excluding missing data)
			my $diffs = sum( ( $assigned_profile != $profile ) & ( $assigned_profile != 0 ) & ( $profile != 0 ) );
			# Count positions where either profile has missing data
			my $missing_in_either = sum( ( $assigned_profile == 0 ) | ( $profile == 0 ) );
			# Calculate distance as percentage difference
			my $distance          = 100.0 * $diffs / ( $locus_count - $missing_in_either ); #100.0 - force float not int
			
			# Track the closest match found so far
			if ( $distance < $closest_distance ) {
				$closest_distance      = $distance;
				$closest_profile_index = $index;
			}
		}
		$index++;
		if ( $opts{'debug'} ) {
			say "Profile $index ordered.";
		} elsif ( $index % 500 == 0 ) {
			if ( !$opts{'quiet'} ) {
				say "Order adjusted for $index profiles.";
			}
		}
	}
	my $reordered_profiles = [ @$new_profiles[ $closest_profile_index .. @$new_profiles - 1 ] ];
	if ( $closest_profile_index > 0 ) {
		push @$reordered_profiles, reverse @$new_profiles[ 0 .. $closest_profile_index - 1 ];
	}
	say 'Done.' if !$opts{'quiet'};
	my $stop_time = time;
	my $duration  = $stop_time - $start_time;
	say "Time taken (adjusting PRIM order): $duration second(s)." if !$opts{'quiet'};
	return $reordered_profiles;
}

sub initiate_log_file {
	my ($filename) = @_;
	open( my $fh, '>', $filename ) || die "Cannot write to log file $filename.\n";
	say $fh qq(profile_id\tclosest profile_id\tcommon alleles\tmissing alleles\tmissing in either\tidentity\tdistance\t)
	  . qq(chosen prefix\tnew LINcode);
	close $fh;
	return;
}

sub assign_lincodes {
	my ($profiles_to_assign) = @_;
	
	my $count                = @$profiles_to_assign;
	my $plural               = $count == 1 ? q() : q(s);
	say "Assigning LINcodes for $count profile$plural." if !$opts{'quiet'};
	
	# Load existing LINcode definitions and thresholds
	my $definitions = get_lincode_definitions();
	my $thresholds  = get_thresholds();
	my %missing     = ( N => 0 );  # Map missing alleles to 0

	foreach my $profile_id (@$profiles_to_assign) {
		my $lincode;
		
		# Load profile data from file
		my $profile = get_profile($profile_id);
		next unless $profile;
		
		# Convert missing alleles (N) to 0 for calculations
		$_ = $missing{$_} // $_ foreach @$profile;
		
		if ( !@{ $definitions->{'profile_ids'} } ) {
			# First profile gets all zeros as LINcode
			$lincode = [ (0) x @{ $thresholds->{'diffs'} } ];
			$definitions->{'profiles'} = pdl($profile);
		} else {
			# Calculate LINcode based on similarity to existing profiles
			$lincode = get_new_lincode( $definitions, $profile_id, $profile );
		}
		
		# Display assigned LINcode
		local $" = q(_);
		my $identifier = "profile-$profile_id";
		my $spaces     = q( ) x abs( 20 - length($identifier) );
		say "$identifier:$spaces@$lincode." if !$opts{'quiet'};
		
		# Write LINcode to file
		assign_lincode( $profile_id, $lincode );
		
		# Add to in-memory definitions for next iteration
		push @{ $definitions->{'profile_ids'} }, $profile_id;
		push @{ $definitions->{'lincodes'} },    $lincode;
	}
	return;
}

sub get_lincode_definitions {
	my $lincode_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_lincodes.tsv");
	my $profile_ids = [];  # List of profile IDs
	my $profiles    = [];  # List of profile allele arrays
	my $lincodes    = [];  # List of assigned LINcodes
	my %missing     = ( N => 0 );  # Map missing alleles to 0
	
	if (-e $lincode_file) {
		open(my $fh, '<', $lincode_file) or die "Cannot read $lincode_file: $!\n";
		my $header = <$fh>; # Skip header
		while (my $line = <$fh>) {
			chomp $line;
			next if $line =~ /^\s*$/;
			my @fields = split(/\t/, $line);
			my $profile_id = $fields[0];
			my $lincode_str = $fields[1];
			my @profile = split(/,/, $fields[2]);
			
			push @$profile_ids, $profile_id;
			$_ = $missing{$_} // $_ foreach @profile;
			push @$profiles, \@profile;
			push @$lincodes, [split(/_/, $lincode_str)];
		}
		close $fh;
	}
	
	return {
		profile_ids => $profile_ids,
		profiles    => @$profiles ? pdl($profiles) : pdl([]),
		lincodes    => $lincodes
	};
}

sub get_new_lincode {
	my ( $definitions, $profile_id, $profile ) = @_;
	
	# Get locus information
	my $loci        = get_scheme_loci();
	my $locus_count = @$loci;
	
	# Add new profile to the matrix of profiles
	$definitions->{'profiles'} = $definitions->{'profiles'}->glue( 1, pdl($profile) );
	my $j            = @{ $definitions->{'profile_ids'} };  # Index of newly added profile (last row)
	my $prof2        = $definitions->{'profiles'}->slice(",($j)");  # The new profile
	
	# Find the closest existing profile
	my $min_distance = 100;  # Start at maximum distance
	my $closest_index;
	my $closest = {};  # Store details for logging

	# Compare new profile against all existing profiles to find closest match
	for my $i ( 0 .. @{ $definitions->{'profile_ids'} } - 1 ) {
		my $prof1             = $definitions->{'profiles'}->slice(",($i)");
		
		# Count differences, excluding positions where either profile has missing data
		my $diffs             = sum( ( $prof1 != $prof2 ) & ( $prof1 != 0 ) & ( $prof2 != 0 ) );
		my $missing_in_either = sum( ( $prof1 == 0 ) | ( $prof2 == 0 ) );
		# Calculate distance as percentage of differing alleles
		my $distance          = 100.0 * $diffs / ( $locus_count - $missing_in_either );
		
		if ( $distance < $min_distance ) {
			$min_distance  = $distance;
			$closest_index = $i;
			# Store additional info for debug logging
			if ( $opts{'log'} ) {
				( $closest->{'common_alleles'} ) = sum( $prof1 == $prof2 );
				( $closest->{'missing'} )        = sum( $prof2 == 0 );
				$closest->{'missing_in_either'} = $missing_in_either;
			}
		}
		
		# If profiles are identical, reuse the same LINcode
		if ( !$diffs ) {
			return $definitions->{'lincodes'}->[$closest_index];
		}
	}
	my $identity        = 100 - $min_distance;
	my $thresholds      = get_thresholds();
	my $threshold_index = 0;
	foreach my $threshold_identity ( @{ $thresholds->{'identity'} } ) {
		if ( $identity >= $threshold_identity ) {
			$threshold_index++;
			next;
		}
		last;
	}
	my $new_lincode = increment_lincode( $definitions->{'lincodes'}, $closest_index, $threshold_index );
	if ( $opts{'log'} ) {
		open( my $fh, '>>', $opts{'log'} ) || die "Cannot append to $opts{'log'}.\n";
		my @chosen_prefix =
		  $threshold_index == 0 ? () : @{ $definitions->{'lincodes'}->[$closest_index] }[ 0 .. $threshold_index - 1 ];
		local $" = q(_);
		say $fh qq($profile_id\t$definitions->{'profile_ids'}->[$closest_index]\t$closest->{'common_alleles'}\t)
		  . qq($closest->{'missing'}\t$closest->{'missing_in_either'}\t$identity\t$min_distance\t@chosen_prefix\t)
		  . qq(@$new_lincode);
		close $fh;
	}
	return $new_lincode;
}

sub increment_lincode {
	my ( $lincodes, $closest_index, $threshold_index ) = @_;
	my $thresholds = get_thresholds();
	
	if ( $threshold_index == 0 ) {
		# Profile is very distant from all others - create new top-level LINcode
		my $max_first = 0;
		# Find highest value in first position across all LINcodes
		foreach my $lincode (@$lincodes) {
			if ( $lincode->[0] > $max_first ) {
				$max_first = $lincode->[0];
			}
		}
		# Increment and pad with zeros
		my @new_lincode = ( ++$max_first, (0) x ( @{ $thresholds->{'diffs'} } - 1 ) );
		return [@new_lincode];
	}
	$closest_index //= 0;
	my $closest_lincode     = $lincodes->[$closest_index];
	my @lincode_prefix      = @$closest_lincode[ 0 .. $threshold_index - 1 ];
	my $max_threshold_index = 0;
	foreach my $lincode (@$lincodes) {
		local $" = q(_);
		next if qq(@lincode_prefix) ne qq(@$lincode[ 0 .. $threshold_index - 1 ]);
		if ( $lincode->[$threshold_index] > $max_threshold_index ) {
			$max_threshold_index = $lincode->[$threshold_index];
		}
	}
	my @new_lincode = @lincode_prefix;
	push @new_lincode, ++$max_threshold_index;
	push @new_lincode, 0 while @new_lincode < @{ $thresholds->{'diffs'} };
	return [@new_lincode];
}

sub assign_lincode {
	my ( $profile_id, $lincode ) = @_;
	
	# Define output file paths
	my $lincode_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_lincodes.tsv");
	my $profile_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_profiles.tsv");
	
	# Retrieve the profile data to include in output
	my $profile = get_profile($profile_id);
	return unless $profile;
	
	# Create file with header if it doesn't exist
	if (!-e $lincode_file) {
		open(my $fh, '>', $lincode_file) or die "Cannot create $lincode_file: $!\n";
		say $fh "profile_id\tlincode\tprofile";
		close $fh;
	}
	
	# Append the lincode
	open(my $fh, '>>', $lincode_file) or die "Cannot append to $lincode_file: $!\n";
	local $" = '_';
	my $lincode_str = "@$lincode";
	local $" = ',';
	my $profile_str = "@$profile";
	say $fh "$profile_id\t$lincode_str\t$profile_str";
	close $fh;
	return;
}

sub get_profiles_without_lincodes {
	print "Retrieving up to $opts{'batch_size'} profiles without LINcodes ..." if !$opts{'quiet'};
	
	# Check that profiles file exists
	my $profile_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_profiles.tsv");
	if (!-e $profile_file) {
		die "Error: Profile file not found: $profile_file\n" .
		    "Please create this file with format: profile_id\tallele1,allele2,allele3,...\n";
	}
	
	# Read assigned lincodes
	my %assigned;
	my $lincode_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_lincodes.tsv");
	if (-e $lincode_file) {
		open(my $fh, '<', $lincode_file) or die "Cannot read $lincode_file: $!\n";
		my $header = <$fh>;
		while (my $line = <$fh>) {
			chomp $line;
			my ($profile_id) = split(/\t/, $line);
			$assigned{$profile_id} = 1;
		}
		close $fh;
	}
	
	# Get profiles without lincodes
	my @profiles;
	open(my $fh, '<', $profile_file) or die "Cannot read $profile_file: $!\n";
	my $header = <$fh>;
	while (my $line = <$fh>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		my ($profile_id, $profile_str) = split(/\t/, $line, 2);
		next if $assigned{$profile_id};
		
		# Check filters
		if ($opts{'x'} && $profile_id < $opts{'x'}) { next; }
		if ($opts{'y'} && $profile_id > $opts{'y'}) { next; }
		
		# Check missing count
		my @alleles = split(/,/, $profile_str);
		my $missing_count = grep { $_ eq 'N' || $_ eq '0' } @alleles;
		next if $missing_count > $opts{'missing'};
		
		push @profiles, $profile_id;
		last if @profiles >= $opts{'batch_size'};
	}
	close $fh;
	
	my $count = @profiles;
	say "$count retrieved." if !$opts{'quiet'};
	return \@profiles;
}

sub get_profile_order_term {
	# For flat files, we just use numeric ordering
	return 'profile_id';
}

sub get_prim_order {
	my ($profiles) = @_;
	
	# Calculate pairwise distance matrix between all profiles
	my ( $filename, $index, $dismat ) = get_distance_matrix($profiles);
	return $index if @$index == 1;  # Skip if only one profile
	
	print 'Calculating PRIM order ...' if !$opts{'quiet'};
	print "\n"                         if @$index >= 500;
	my $start_time = time;
	for my $i ( 0 .. @$index - 1 ) {
		$dismat->set( $i, $i, 999 );
	}
	my $ind = $dismat->flat->minimum_ind;
	my ( $x, $y ) = ( int( $ind / @$index ), $ind - int( $ind / @$index ) * @$index );
	my $index_order   = [ $x, $y ];
	my $profile_order = [ $index->[$x], $index->[$y] ];
	$dismat->set( $x, $y, 999 );
	$dismat->set( $y, $x, 999 );
	while ( @$profile_order != @$index ) {
		my $min = 101;
		my $v_min;
		foreach my $x (@$index_order) {
			my $this_min = $dismat->slice($x)->min;
			if ( $this_min < $min ) {
				$min   = $this_min;
				$v_min = $x;
			}
		}
		my $k = $dismat->slice($v_min)->flat->minimum_ind;
		for my $i (@$index_order) {
			$dismat->set( $i, $k, 999 );
			$dismat->set( $k, $i, 999 );
		}
		push @$index_order,   $k;
		push @$profile_order, $index->[$k];
		my $count = @$profile_order;
		if ( $opts{'debug'} ) {
			say "Profile $count ordered.";
		} elsif ( $count % 500 == 0 ) {
			if ( !$opts{'quiet'} ) {
				say "Order calculated for $count profiles.";
			}
		}
	}
	say 'Done.' if !$opts{'quiet'};
	unlink $filename;
	unlink "$filename.hdr";
	my $stop_time = time;
	my $duration  = $stop_time - $start_time;
	say "Time taken (calculating PRIM order): $duration second(s)." if !$opts{'quiet'};
	# reconnect is a no-op for flat files
	return $profile_order;
}

sub get_distance_matrix {
	my ($profile_ids) = @_;
	my $loci          = get_scheme_loci();
	my $locus_count   = @$loci;
	die "Scheme has no loci.\n" if !$locus_count;
	
	my $matrix      = [];
	my $index       = [];
	my $max_missing = $opts{'missing'};

	foreach my $profile_id (@$profile_ids) {
		my $profile_data = get_profile($profile_id);
		next unless $profile_data;
		
		my $Ns = 0;
		foreach my $allele ( @$profile_data ) {
			if ( $allele eq 'N' || $allele eq '0' ) {
				$Ns++;
				$allele = 0;
			}
		}
		next if $Ns > $max_missing;
		push @$index,  $profile_id;
		push @$matrix, $profile_data;
	}
	my $profile_matrix = pdl($matrix);
	my $count          = @$index;
	die "No profiles to assign.\n" if ( !$count );
	return $index                  if @$index == 1;
	my $msg = $count > 2000 ? ' (this will take a while)' : q();
	print "Calculating distance matrix$msg ..." if !$opts{'quiet'};
	print "\n"                                  if $count >= 500;
	my $start_time = time;
	my ( $fh, $filename ) = File::Temp::tempfile(
		'dismatXXXX',
		DIR    => $script->{'config'}->{'secure_tmp_dir'},
		SUFFIX => '.dismat',
		UNLINK => 0
	);
	close $fh;
	my $dismat =
	  $opts{'mmap'}
	  ? mapfraw( $filename, { Creat => 1, Dims => [ $count, $count ], Datatype => float } )
	  : zeroes( float, $count, $count );

	for my $i ( 0 .. $count - 1 ) {
		if ( $opts{'debug'} ) {
			say "Profile $i.";
		} elsif ( $i && $i % 500 == 0 ) {
			if ( !$opts{'quiet'} ) {
				say "Calculated for $i profiles.";
				if ( $i == 500 && $count > 2000 ) {
					say 'Note that it does speed up (matrix calculations are for upper triangle)!';
				}
			}
		}
		for my $j ( $i + 1 .. $count - 1 ) {
			my $prof1             = $profile_matrix->slice(",($i)");
			my $prof2             = $profile_matrix->slice(",($j)");
			my $diffs             = sum( ( $prof1 != $prof2 ) & ( $prof1 != 0 ) & ( $prof2 != 0 ) );
			my $missing_in_either = sum( ( $prof1 == 0 ) | ( $prof2 == 0 ) );
			my $distance = 100.0 * $diffs / ( $locus_count - $missing_in_either );    #100.0 - force float not int.
			$dismat->set( $i, $j, $distance );
			$dismat->set( $j, $i, $distance );
		}
	}
	say 'Done.' if !$opts{'quiet'};
	my $stop_time = time;
	my $duration  = $stop_time - $start_time;
	say "Time taken (distance matrix): $duration second(s)." if !$opts{'quiet'};
	# reconnect is a no-op for flat files
	return ( $filename, $index, $dismat );
}

sub get_thresholds {
	if ( $script->{'cache'}->{'thresholds'} ) {
		return $script->{'cache'}->{'thresholds'};
	}
	
	# Read thresholds from semicolon-separated file (e.g., "2;4;7;14;21;35;70;140")
	my $thresholds_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_thresholds.txt");
	open(my $fh, '<', $thresholds_file) or die "Cannot read $thresholds_file: $!\n";
	my $thresholds = <$fh>;
	chomp $thresholds;
	close $fh;
	
	my $diffs    = [ split /\s*;\s*/x, $thresholds ];
	my $identity = [];
	my $loci     = get_scheme_loci();
	foreach my $diff (@$diffs) {
		push @$identity, 100 * ( @$loci - $diff ) / @$loci;
	}
	$script->{'cache'}->{'thresholds'} = {
		diffs    => $diffs,
		identity => $identity
	};
	return $script->{'cache'}->{'thresholds'};
}

sub get_profile {
	my ($profile_id) = @_;
	
	# Read from tab-separated file: profile_id<TAB>allele1,allele2,allele3,...
	my $profile_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_profiles.tsv");
	open(my $fh, '<', $profile_file) or die "Cannot read $profile_file: $!\n";
	my $header = <$fh>;  # Skip header line
	
	while (my $line = <$fh>) {
		chomp $line;
		my ($id, $profile_str) = split(/\t/, $line, 2);
		
		# Return allele array when matching profile is found
		if ($id eq $profile_id) {
			close $fh;
			# Handle both comma-separated and tab-separated formats
			my @alleles;
			if ($profile_str =~ /,/) {
				@alleles = split(/,/, $profile_str);
			} else {
				# Tab-separated alleles
				@alleles = split(/\t/, $profile_str);
			}
			# Clean up alleles: handle empty strings and semicolon-separated values
			foreach my $allele (@alleles) {
				if ($allele eq '' || $allele eq '-') {
					$allele = 'N';  # Missing data
				} elsif ($allele =~ /;/) {
					# Multiple alleles - take the first one
					($allele) = split(/;/, $allele);
					$allele = 'N' if $allele eq '';  # Handle edge case
				}
			}
			return \@alleles;
		}
	}
	close $fh;
	return undef;  # Profile not found
}

sub concatenate_profiles {
	my $input_file = $opts{'input_profiles'};
	my $target_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_profiles.tsv");
	
	# Verify input file exists
	if ( !-e $input_file ) {
		die "Error: Input profiles file not found: $input_file\n";
	}
	
	say "Concatenating profiles from $input_file ..." if !$opts{'quiet'};
	
	# Read existing profile IDs to avoid duplicates
	my %existing_profiles;
	if ( -e $target_file ) {
		open(my $fh, '<', $target_file) or die "Cannot read $target_file: $!\n";
		my $header = <$fh>;  # Skip header
		while (my $line = <$fh>) {
			chomp $line;
			next if $line =~ /^\s*$/;
			my ($profile_id) = split(/\t/, $line);
			$existing_profiles{$profile_id} = 1;
		}
		close $fh;
	}
	
	# Open target file for appending
	open(my $out_fh, '>>', $target_file) or die "Cannot append to $target_file: $!\n";
	
	# Read and append new profiles
	open(my $in_fh, '<', $input_file) or die "Cannot read $input_file: $!\n";
	my $header = <$in_fh>;  # Skip header line
	my $added_count = 0;
	my $skipped_count = 0;
	
	while (my $line = <$in_fh>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		
		my @fields = split(/\t/, $line);
		my $profile_id = shift @fields;
		
		# Skip if profile already exists
		if ( $existing_profiles{$profile_id} ) {
			$skipped_count++;
			next;
		}
		
		# Convert tab-separated alleles to comma-separated
		# Clean up alleles: handle empty strings and semicolon-separated values
		foreach my $allele (@fields) {
			if ($allele eq '' || $allele eq '-') {
				$allele = 'N';  # Missing data
			} elsif ($allele =~ /;/) {
				# Multiple alleles - take the first one
				($allele) = split(/;/, $allele);
				$allele = 'N' if $allele eq '';  # Handle edge case
			}
		}
		my $profile_str = join(',', @fields);
		
		# Append the profile in comma-separated format
		say $out_fh "$profile_id\t$profile_str";
		$existing_profiles{$profile_id} = 1;
		$added_count++;
	}
	
	close $in_fh;
	close $out_fh;
	
	say "Added $added_count profile(s), skipped $skipped_count duplicate(s)." if !$opts{'quiet'};
	return;
}

sub check_db {
	my $dir = $opts{'dir'};
	my $scheme_id = $opts{'scheme_id'};
	
	# Check for thresholds file (required for LINcode calculation)
	if ( !are_lincodes_defined() ) {
		die "Error: LINcode thresholds are not defined for scheme $scheme_id.\n" .
		    "Please create file: $dir/scheme_${scheme_id}_thresholds.txt\n" .
		    "Format: semicolon-separated list of allelic difference thresholds (e.g., 2;4;7;14;21;35;70;140)\n";
	}
	
	# Check for loci file
	my $loci_file = "$dir/scheme_${scheme_id}_loci.txt";
	if ( !-e $loci_file ) {
		die "Error: Loci file not found: $loci_file\n" .
		    "Please create this file with one locus name per line.\n";
	}
}

sub check_if_script_already_running {
	my $lock_file = get_lock_file();
	if ( -e $lock_file ) {
		open( my $fh, '<', $lock_file ) || die "Cannot open lock file $lock_file for reading: $!\n";
		my $pid = <$fh>;
		close $fh;
		my $pid_exists = kill( 0, $pid );
		if ( !$pid_exists ) {
			say 'Lock file exists but process is no longer running - deleting lock.'
			  if !$opts{'quiet'};
			unlink $lock_file;
		} else {
			say 'Script already running with these parameters - terminating.' if !$opts{'quiet'};
			exit(1);
		}
	}
	open( my $fh, '>', $lock_file ) || die "Cannot open lock file $lock_file for writing: $!\n";
	say $fh $$;
	close $fh;
	return;
}

sub get_lock_file {
	my $hash      = Digest::MD5::md5_hex("$0||$opts{'dir'}||$opts{'scheme_id'}");
	my $lock_dir  = $script->{'config'}->{'lock_dir'};
	# Fallback if script not yet initialized
	unless ($lock_dir) {
		$lock_dir = File::Spec->catdir($opts{'dir'}, '.locks');
		make_path($lock_dir) unless -d $lock_dir;
	}
	my $lock_file = "$lock_dir/lincodes_$hash";
	return $lock_file;
}

sub remove_lock_file {
	my $lock_file = get_lock_file();
	unlink $lock_file;
	return;
}

sub show_help {
	print <<'HELP';
NAME
    lincodes_LK.pl - Define LINcodes from cgMLST profiles

USAGE
    lincodes_LK.pl --dir DIRECTORY --scheme SCHEME_ID [options]

REQUIRED
    --dir DIRECTORY       Directory with scheme files
    --scheme SCHEME_ID    Scheme ID number

OPTIONS
    --batch_size NUMBER   Max profiles per batch (default: 10000)
    --create             Create example schema directory structure
    --input_profiles FILE External profiles to add
    --log FILE           Debug log file
    --missing NUMBER     Max missing loci allowed (default: 0)
    --mmap               Use disk for distance matrix (slower, less memory)
    --quiet, -q          Only show errors
    --debug              Verbose output
    --min, -x ID         Minimum profile ID
    --max, -y ID         Maximum profile ID

FILES
    scheme_<ID>_profiles.tsv   - Profile data (id<TAB>allele1,allele2,...)
    scheme_<ID>_loci.txt       - Locus names (one per line)
    scheme_<ID>_thresholds.txt - Thresholds (e.g., 50;100;200;400;800;1600;3200)
    scheme_<ID>_lincodes.tsv   - Output file

EXAMPLES
    lincodes_LK.pl --create --dir example --scheme 1
    lincodes_LK.pl --dir example --scheme 1
    lincodes_LK.pl --dir data --scheme 1 --input_profiles new.tsv --batch_size 5000
HELP
	return;
}

################################################################################
# Standalone helper functions to replace BIGSdb dependencies
################################################################################

sub create_example_schema {
	my $dir = $opts{'dir'};
	my $scheme_id = $opts{'scheme_id'};
	
	# Create directory if it doesn't exist
	make_path($dir) unless -d $dir;
	
	# Create empty profiles file with header comment
	my $profiles_file = File::Spec->catfile($dir, "scheme_${scheme_id}_profiles.tsv");
	open(my $pf, '>', $profiles_file) or die "Cannot create $profiles_file: $!\n";
	print $pf "# profile_id<TAB>allele1,allele2,allele3,...\n";
	close $pf;
	
	# Create loci file with example loci
	my $loci_file = File::Spec->catfile($dir, "scheme_${scheme_id}_loci.txt");
	open(my $lf, '>', $loci_file) or die "Cannot create $loci_file: $!\n";
	print $lf "# One locus name per line\n";
	for my $i (1..10) {
		print $lf "locus_$i\n";
	}
	close $lf;
	
	# Create thresholds file with example thresholds
	my $thresholds_file = File::Spec->catfile($dir, "scheme_${scheme_id}_thresholds.txt");
	open(my $tf, '>', $thresholds_file) or die "Cannot create $thresholds_file: $!\n";
	print $tf "50;100;200;400;800;1600;3200\n";
	close $tf;
	
	# Create empty lincodes file
	my $lincodes_file = File::Spec->catfile($dir, "scheme_${scheme_id}_lincodes.tsv");
	open(my $lcf, '>', $lincodes_file) or die "Cannot create $lincodes_file: $!\n";
	close $lcf;
	
	# Create subdirectories
	my $locks_dir = File::Spec->catdir($dir, '.locks');
	my $tmp_dir = File::Spec->catdir($dir, '.tmp');
	make_path($locks_dir, $tmp_dir);
	
	return;
}

sub initialize_script {
	my $script = {};
	my $tmp_dir = File::Spec->catdir($opts{'dir'}, '.tmp');
	make_path($tmp_dir) unless -d $tmp_dir;
	$script->{'config'}->{'secure_tmp_dir'} = $tmp_dir;
	my $lock_dir = File::Spec->catdir($opts{'dir'}, '.locks');
	make_path($lock_dir) unless -d $lock_dir;
	$script->{'config'}->{'lock_dir'} = $lock_dir;
	$script->{'cache'} = {};
	return $script;
}

sub get_scheme_loci {
	my $loci_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_loci.txt");
	my $loci = [];
	if (-e $loci_file) {
		open(my $fh, '<', $loci_file) or die "Cannot read $loci_file: $!\n";
		while (my $line = <$fh>) {
			chomp $line;
			next if $line =~ /^\s*#/ || $line =~ /^\s*$/;
			push @$loci, $line;
		}
		close $fh;
	} else {
		die "Loci file not found: $loci_file\nPlease create this file with one locus name per line.\n";
	}
	return $loci;
}

sub are_lincodes_defined {
	my $thresholds_file = File::Spec->catfile($opts{'dir'}, "scheme_$opts{'scheme_id'}_thresholds.txt");
	return -e $thresholds_file ? 1 : 0;
}


