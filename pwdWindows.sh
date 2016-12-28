#!/bin/bash

# Get the Windows network path

# Get the actual path with resolved symlinks
pwd=$(pwd -P);
if [ "$1" != "" ]; then
  pwd=$(realpath $1)
fi
# remove scicomp and leading slash
pwd=$(sed 's|^/scicomp||' <<< $pwd);
# is this in the home directory?
pwd=$(sed "s|^/home/$USER|/home|" <<< $pwd);
# / to \
pwd=$(sed 's|/|\\|g' <<< $pwd);
# tack on the domain name
pwd="\\\\data.biotech.cdc.gov$pwd"


echo $pwd;
