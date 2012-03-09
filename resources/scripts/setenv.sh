#!/bin/bash
# NEED To chage AGV_HOME and PGV_HOME to something more neutral... MALGEN_HOME ?
# For now AGV_HOME is assumed to point to the top directory of this software

export AGV_HOME=@install.dir@
export PGV_HOME=@install.dir@
export PERL5LIB=${PGV_HOME}/lib/perl
export PATH=${PGV_HOME}/bin:${PATH}

`${PGV_HOME}/bin/@install.name@ query samtrak-env sh`
