#!/bin/bash

# Get the Linux path from a windows path
WINDOWS=$1

# Change backslashes to forward slashes
LINUX=$(sed 's|\\|/|g' <<< $WINDOWS)
# Remove extra leading slashes
LINUX=$(sed 's|^/\+|/|' <<< $LINUX)

# Change the domain name
LINUX=$(sed 's|^/data.biotech.cdc.gov/|/scicomp/|' <<< $LINUX);

# Print the final linux path
echo $LINUX;
