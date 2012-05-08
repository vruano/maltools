package Maltools::Engines::Print;

use strict;
use warnings;
use base 'Maltools::Engine';
use Maltools::Tool;

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
