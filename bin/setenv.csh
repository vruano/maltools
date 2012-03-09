#!/bin/tcsh
# NEED To chage AGV_HOME and PGV_HOME to something more neutral... MALGEN_HOME ?
# For now AGV_HOME is assumed to point to the top directory of this software
setenv AGV_HOME /nfs/team112/PGV_RD
setenv PGV_HOME /nfs/team112/PGV_RD
setenv PERLVER  5.10.1
setenv PERLARCH x86_64-linux-thread-multi
setenv PERLHOME /software/perl-${PERLVER}
setenv IRODS_HOME /software/irods/icommands
setenv SAMTRAK_HOME /nfs/team112/samtrak
setenv JAVA_HOME /software/jdk

setenv LD_LIBRARY_PATH ${AGV_HOME}/lib
setenv LIBRARY_PATH ${AGV_HOME}/lib
setenv C_INCLUDE_PATH ${AGV_HOME}/include

setenv PERL5LIB ${AGV_HOME}/perl/modules
setenv PERL5LIB ${PERL5LIB}:${AGV_HOME}/perl/lib/perl5/${PERLVER}/${PERLARCH}:${SAMTRAK_HOME}/lib/perl5/${PERLVER}/${PERLARCH}
setenv PERL5LIB ${PERL5LIB}:${SAMTRAK_HOME}/share/samtrak/modules
setenv PERL5LIB ${PERL5LIB}:${AGV_HOME}/perl/lib/perl5/${PERLVER}:${SAMTRAK_HOME}/lib/perl5/${PERLVER}
setenv PERL5LIB ${PERL5LIB}:${PERLHOME}/lib/site_perl/${PERLVER}/${PERLARCH}
setenv PERL5LIB ${PERL5LIB}:${PERLHOME}/lib/site_perl/${PERLVER}
setenv PERL5LIB ${PERL5LIB}:${PERLHOME}/lib/${PERLVER}/${PERLARCH}
setenv PERL5LIB ${PERL5LIB}:${PERLHOME}/lib/${PERLVER}

setenv PATH ${JAVA_HOME}/bin:${PATH}
setenv PATH ${PERLHOME}/bin:${PATH}
setenv PATH ${IRODS_HOME}/bin:${PATH}
setenv PATH ${AGV_HOME}/bin:${PATH}

setenv ANT_OPTS "-Xmx500m -Dhttp.proxyHost=webcache -Dhttp.proxyPort=3128"
