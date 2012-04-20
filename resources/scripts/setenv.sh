#!/bin/bash

export MALTOOLS_HOME=@install.dir@
export PERL5LIB=${MALTOOLS_HOME}/lib/perl
export PATH=${MALTOOLS_HOME}/bin:${PATH}

`${MALTOOLS_HOME}/bin/@install.name@ query samtrak-env sh`
