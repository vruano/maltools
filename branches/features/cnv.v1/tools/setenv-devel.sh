#!/bin/bash

stat ./tools/setenv-devel.sh 2> /dev/null 1> /dev/null

if [ "$?" = "0" ]; then
  export MALTOOLS_HOME=`pwd`/build/dist/layout-devel
  export PERL5LIB=`pwd`/src/perl/lib:${MALTOOLS_HOME}/lib/perl
  export PATH=${MALTOOLS_HOME}/bin:${PATH}

  `maltools-dev query samtrak-env sh`
else
  echo "ERROR!"
  echo ""
  echo " You need to run or source this script from the project checkout base directory"
  echo ""
  echo "Sorry."
fi
