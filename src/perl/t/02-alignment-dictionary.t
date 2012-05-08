use warnings;
use strict;
use File::Spec::Functions qw(catfile);

use Test::More;

use Maltools::Config;

my $config = Maltools::Config->new(file => "$ENV{PGV_HOME}/conf/pgv-dev.json", program => 'pgv-dev');

my $data_dir = $config->get_datadir();

unless ($data_dir && -d $data_dir) {
  diag("Cannot test due to lack of data");
  goto FINISH;
}

my $AlignmentClass = 'Maltools::Alignment';

use_ok $AlignmentClass or goto FINISH;

my $aln = new_ok ($AlignmentClass => [ file => catfile($data_dir,'pgv','testdata','truncated-pgv.bam')]) or goto FINISH;

my $dictionary = $aln->dictionary();
ok($dictionary, "checking non-null dictonary");

is($#$dictionary+1,17,"checking dictionary length");
my $mal14 = $dictionary->[5];
is ($mal14->{name},'MAL14',"checking some chromosome name");
is ($mal14->{length},3291871,"checking some chromosome name");

FINISH:
done_testing();
