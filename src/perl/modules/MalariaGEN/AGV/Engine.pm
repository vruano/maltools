package MalariaGEN::AGV::Engine;

use strict;
use warnings;

use Moose;

has 'name' => (is => 'ro', isa => 'Str', required => 1);
has 'running_options' => (is => 'rw', isa => 'HashRef', default => sub { {} });




sub run_job {
   my ($self,%args) = @_;
   die "yet not implemented\n";
}

sub engine_by_name {
   my ($class,$name,%args)  = @_;
   $name ||= "auto";
   return $class->guess_engine(%args) if ($name eq "auto");
   my $cls = "MalariaGEN::AGV::Engines::" . uc(substr($name,0,1)) . substr($name,1);
   eval "require $cls" or die "not such a engine '$name'";
   return $cls->new(%args);
}

sub guess_engine {
   my ($class,%args) = @_;
   my $domain = `hostname -d`;
   chomp $domain;
   my $name;
   if ($domain =~ /sanger.ac.uk$/) {
     `which lsid`;
     $name = $? ? 'local':'sanger'; 
   }
   elsif ($domain =~ /well.ox.ac.uk$/) {
     `which qstat`;
     $name = $? ? 'local':'oxford';
   }
   else {
     $name = 'local';
   }
   return $class->engine_by_name($name,%args);
}

sub import_resource {
   return $_[1];
}

1;
