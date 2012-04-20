#!/bin/tcsh

setenv MALTOOLS_HOME @install.dir@
setenv PERL5LIB ${MALTOOLS_HOME}/lib/perl
setenv PATH ${MALTOOLS_HOME}/bin:${PATH}

`${MALTOOLS_HOME}/bin/@install.name@ query samtrak-env csh`

