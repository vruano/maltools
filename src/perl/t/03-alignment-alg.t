use warnings;
use strict;
use File::Spec::Functions qw(catfile);

use Test::More;

use MalariaGEN::AGV::Config;

my $config = MalariaGEN::AGV::Config->new(file => "$ENV{PGV_HOME}/conf/pgv-dev.json", program => 'pgv-dev');

my $data_dir = $config->get_datadir();

unless ($data_dir && -d $data_dir) {
  diag("Cannot test due to lack of data");
  goto FINISH;
}

my $AlignmentClass = 'MalariaGEN::AGV::Alignment';


use_ok $AlignmentClass or goto FINISH;

my $pgv_aln = new_ok ($AlignmentClass => [ file => catfile($data_dir,'pgv','testdata','truncated-pgv.bam')]) or goto FINISH;
my $aln = new_ok ($AlignmentClass => [ file => catfile($data_dir,'pgv','testdata','truncated-fhdr.bam')]) or goto FINISH;

ok(!$pgv_aln->valid_algorithm('anything'), "check unsorted alignment with random algorithm");
ok(!$aln->valid_algorithm('anyhing'),"check sorted alignment with random algorithm");
ok(!$aln->valid_algorithm('anyhing'),"check sorted alignment with random algorithm");
ok($aln->valid_algorithm('bwa aln'),"check sorted alignment with bwa_aln");
ok($aln->valid_algorithm(name => 'bwa aln'),"check sorted alignment with bwa_aln");
ok($aln->valid_algorithm(name => "bwa aln"),"check sorted alignment with bwa_aln");
ok(!$aln->valid_algorithm(name => 'bwa aln',version=>"wrong"),"check sorted alignment with bwa_aln");
ok($aln->valid_algorithm(name => 'bwa aln',version=>"0.5.9-r16"),"check sorted alignment with bwa_aln");

FINISH:
done_testing();
