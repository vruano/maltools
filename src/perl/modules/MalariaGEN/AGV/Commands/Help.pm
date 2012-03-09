package MalariaGEN::AGV::Commands::Help;

use warnings;
use strict;

use File::Spec::Functions qw(splitpath catpath);
use base 'MalariaGEN::AGV::Command'; 

sub help_summary {
  return 'shows extended information on available commands';
}

sub help_text {
  my ($self) = @_;
  my $prog = $self->cl_name;
  my $summary = $self->help_summary;
  my $command_table = join("\n  ", sort grep { defined $_ } map { command_summary_line($_) } command_list());
  return <<EOM
Program:
  
  $prog - pipeline and analysis work command line interface base script.

Commands:
  
  $command_table

EOM

}

sub command_summary_line {
  my ($command) = @_;  
  eval "require $command";
  warn $@ if $@;
  $@ and return undef;
  $command->isa('MalariaGEN::AGV::Command') or return undef;
  return $command->name .  (" " x (22 - length($command->name))) . "\t" . $command->help_summary;
}

sub command_list {
  my $self_file = $INC{'MalariaGEN/AGV/Commands/Help.pm'};
  my ($vol,$dir) = splitpath($self_file);
  my $command_dir = catpath($vol,$dir,'');
  my @result = ();
  opendir(my $dh, $command_dir) or die "cannot open command package directory $command_dir\n";
  
  while (my $file = readdir($dh)) {
    next unless $file =~ /(\w+)\.pm$/;
    push @result, 'MalariaGEN::AGV::Commands::' . $1;
  }
  closedir($dh);
  return @result;
}

sub execute {
  my ($self,$site,$params) = @_;
  my @args = @{$params->{arguments}};
  if (scalar(@args) == 0) {
    print $self->help_text();
    return $self->ok_return();
  }
  while (my $cmd_name = shift(@args)) {
    next unless $cmd_name =~ /\S/;
    my $cmd = MalariaGEN::AGV::Command->instance($cmd_name);
    if (defined $cmd) {
      print STDOUT  $cmd->help_text; 
      return $self->ok_return();
    }
    else {
      return $self->error_return('unknown command ' . $cmd_name);
    }
  }
}

1;
