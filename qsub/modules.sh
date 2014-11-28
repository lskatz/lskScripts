
module() { eval `/usr/bin/modulecmd bash $*`; }
#export -f module

MODULESHOME=/usr/share/Modules
export MODULESHOME

if [ "${LOADEDMODULES:-}" = "" ]; then
  LOADEDMODULES=
  export LOADEDMODULES
fi

## july 2013 wgx9 chris dagdigian
## master node had incorrect modulepath for centos nodes and the correct
## path was not being set because the env was already present and defined
## commenting out the if-then test so we force the recreation of the
## proper centos-specific module file paths ...

#if [ "${MODULEPATH:-}" = "" ]; then
#  MODULEPATH=`sed -n 's/[ 	#].*$//; /./H; $ { x; s/^\n//; s/\n/:/g; p; }' ${MODULESHOME}/init/.modulespath`
#  export MODULEPATH
#fi

MODULEPATH=`sed -n 's/[ 	#].*$//; /./H; $ { x; s/^\n//; s/\n/:/g; p; }' ${MODULESHOME}/init/.modulespath`
export MODULEPATH


if [ ${BASH_VERSINFO:-0} -ge 3 ] && [ -r ${MODULESHOME}/init/bash_completion ]; then
 . ${MODULESHOME}/init/bash_completion
fi
