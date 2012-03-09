use Test::More;

use strict;
use warnings;

my $ManifestPipeline = 'MalariaGEN::AGV::Manifest::Pipeline';

goto FINISH unless use_ok $ManifestPipeline;
goto FINISH unless my $man = new_ok $ManifestPipeline;

is($man->class,'pipeline',"checking class");
eval { $man->load(file => 'resources/test/pgv/Hb3xDd2.json') };

unless (ok(!$@,"loading manifest file")) {
  diag($@);
  goto FINISH;
} 

is($man->get_variable('basedir'),'/data/malariagen/pipelines/PGV-production-2.0',"checking get_variable");
eval { $man->set_variable(rampant => '///'); };
if (!ok(!$@,"checking set_variable")) {
    diag($@);
    goto FINISH;
}

$man = $ManifestPipeline->new();
eval { $man->set_variable(basedir => 'build/perl/test',datadir => 'resources/test/data') };
if (!ok(!$@,"checking multi variable set")) {
    diag($@);
    goto FINISH;
}
my ($basedir,$datadir) = (undef,undef);
eval { ($basedir,$datadir) = $man->get_variable(qw(basedir datadir)) };
if (!ok(!$@,"checking multi variable get")) {
    diag($@);
    goto FINISH;
}

is($basedir,'build/perl/test',"checking mutli variable get 1");
is($datadir,'resources/test/data',"checking multi variable get 2");


$man->load(file => 'resources/test/pgv/Hb3xDd2.json');

is($man->get_variable('basedir'),'build/perl/test',"testing load do not override variables \$basedir");
is($man->get_variable('datadir'),'resources/test/data',"testing load do not override variables \$datadir");

is($man->get_path('snpsdata'),$basedir . '/snpsdata',"checking interpolation has worked");


my $man_bak = $man->clone();

use Samtrak::Site;

$man->auto_complete(samtrak => Samtrak::Site->new());





FINISH:

done_testing();
