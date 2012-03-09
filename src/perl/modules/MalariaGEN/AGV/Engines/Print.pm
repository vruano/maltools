package MalariaGEN::AGV::Engines::Print;

use strict;
use warnings;
use base 'MalariaGEN::AGV::Engine';
use MalariaGEN::AGV::Tool;

sub run_job {
  my ($self,%args) = @_;
  my $job = $args{job};
  my $tool = $args{tool};
  my $script => $job->script(E => \$self, J => \$job, T => \$tool , S => "/dev/null", SH => "/dev/null", %args);
  print STDOUT '-' x 10 . "\n";
  print STDOUT $script,"\n";
  return $? == 0; 
}

1;
