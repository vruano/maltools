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

my $pgv_aln = new_ok ($AlignmentClass => [ file => catfile($data_dir,'pgv','testdata','truncated-pgv.bam')]) or goto FINISH;
my $aln = new_ok ($AlignmentClass => [ file => catfile($data_dir,'pgv','testdata','truncated-fhdr.bam')]) or goto FINISH;

ok(!$pgv_aln->is_sorted, "check unsorted alignment");
ok($aln->is_sorted,"check sorted alignment");

FINISH:
done_testing();
