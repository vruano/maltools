#!/bin/tcsh
# NEED To chage AGV_HOME and PGV_HOME to something more neutral... MALGEN_HOME ?
# For now AGV_HOME is assumed to point to the top directory of this software

setenv AGV_HOME @install.dir@
setenv PGV_HOME @install.dir@
setenv PERL5LIB ${PGV_HOME}/lib/perl
setenv PATH ${PGV_HOME}/bin:${PATH}

`${PGV_HOME}/bin/@install.name@ query samtrak-env csh`

