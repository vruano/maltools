#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use File::Spec::Functions qw(catfile);

use_ok 'Maltools::Reference::RODSet';
use_ok 'Maltools::Reference::ROD';
use_ok 'Maltools::Reference';
use_ok 'Maltools::Config' or BAIL_OUT('cannot instanciate the config object');

my $config = Maltools::Config->new() or BAIL_OUT('cannot instancate the config object');
my $ref_dir = $config->get_resource(qw(test mockups reference Ref1));
$ref_dir and -d $ref_dir or BAIL_OUT('Missing testing reference ' . $config->get_resource());
diag ("reference at '$ref_dir'");
my $ref1 = new_ok 'Maltools::Reference' , [ directory => $ref_dir ];

is($ref1->sequence_file,catfile($ref_dir,'sequences.fa'));
my $rods = $ref1->rods;
ok($rods->isa('Maltools::Reference::RODSet'));

my $uniqueness_rod = $rods->get('uniqueness');
ok($uniqueness_rod);
is($uniqueness_rod->type, 'UQN');
is($uniqueness_rod->file, catfile($ref_dir,'uniqueness.uq'));

done_testing();

1;
