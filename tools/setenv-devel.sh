#!/bin/bash
# NEED To chage AGV_HOME and PGV_HOME to something more neutral... MALGEN_HOME ?
# For now AGV_HOME is assumed to point to the top directory of this software

export AGV_HOME=`pwd`/build/dist/layout-devel
export PGV_HOME=`pwd`/build/dist/layout-devel
export PERL5LIB=${PGV_HOME}/lib/perl
export PATH=${PGV_HOME}/bin:${PATH}

`pgv-dev query samtrak-env sh`
