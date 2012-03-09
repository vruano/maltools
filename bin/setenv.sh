#!/bin/tcsh
# NEED To chage AGV_HOME and PGV_HOME to something more neutral... MALGEN_HOME ?
# For now AGV_HOME is assumed to point to the top directory of this software
export AGV_HOME=/nfs/team112/PGV_RD
export PGV_HOME=/nfs/team112/PGV_RD
export PERLVER=5.10.1
export PERLARCH=x86_64-linux-thread-multi
export PERLHOME=/software/perl-${PERLVER}
export IRODS_HOME=/software/irods/icommands
export SAMTRAK_HOME=/nfs/team112/samtrak
export JAVA_HOME=/software/jdk

export LD_LIBRARY_PATH=${AGV_HOME}/lib
export LIBRARY_PATH=${AGV_HOME}/lib
export C_INCLUDE_PATH=${AGV_HOME}/include

export PERL5LIB=${AGV_HOME}/perl/modules
export PERL5LIB=${PERL5LIB}:${AGV_HOME}/perl/lib/perl5/${PERLVER}/${PERLARCH}:${SAMTRAK_HOME}/lib/perl5/${PERLVER}/${PERLARCH}
export PERL5LIB=${PERL5LIB}:${SAMTRAK_HOME}/share/samtrak/modules
export PERL5LIB=${PERL5LIB}:${AGV_HOME}/perl/lib/perl5/${PERLVER}:${SAMTRAK_HOME}/lib/perl5/${PERLVER}
export PERL5LIB=${PERL5LIB}:${PERLHOME}/lib/site_perl/${PERLVER}/${PERLARCH}
export PERL5LIB=${PERL5LIB}:${PERLHOME}/lib/site_perl/${PERLVER}
export PERL5LIB=${PERL5LIB}:${PERLHOME}/lib/${PERLVER}/${PERLARCH}
export PERL5LIB=${PERL5LIB}:${PERLHOME}/lib/${PERLVER}

export PATH=${JAVA_HOME}/bin:${PATH}
export PATH=${PERLHOME}/bin:${PATH}
export PATH=${IRODS_HOME}/bin:${PATH}
export PATH=${AGV_HOME}/bin:${PATH}

export ANT_OPTS="-Xmx500m -Dhttp.proxyHost=webcache -Dhttp.proxyPort=3128"
