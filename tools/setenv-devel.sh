#!/bin/bash
# NEED To chage AGV_HOME and PGV_HOME to something more neutral... MALGEN_HOME ?
# For now AGV_HOME is assumed to point to the top directory of this software

stat ./tools/setenv-devel.sh 2> /dev/null 1> /dev/null

if [ "$?" = "0" ]; then
  export AGV_HOME=`pwd`/build/dist/layout-devel
  export PGV_HOME=`pwd`/build/dist/layout-devel
  export PERL5LIB=`pwd`/src/perl/lib:${PGV_HOME}/lib/perl
  export PATH=${PGV_HOME}/bin:${PATH}

  `pgv-dev query samtrak-env sh`
else
  echo "ERROR!"
  echo ""
  echo " You need to run or source this script from the project checkout base directory"
  echo ""
  echo "Sorry."
fi
