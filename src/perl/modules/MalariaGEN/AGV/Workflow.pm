package MalariaGEN::AGV::Workflow;
use strict;
use warnings;
use Moose;
use Scalar::Util;
use MalariaGEN::AGV::Workflow::Step;
use MalariaGEN::AGV::Engines::Local;

has 'steps' => ( is => 'ro', isa => 'ArrayRef', required => 1 );
has 'locals' => ( is => 'ro', isa => 'ArrayRef', default => sub { [] });
has '_steps_by_id' => ( is => 'ro', isa => 'HashRef', lazy => 1, builder => '_steps_by_id_builder');
has '_deps' => ( is => 'ro', isa => 'ArrayRef', builder => '_deps_builder', lazy => 1 );
has 'vars' => ( is => 'ro', isa => 'HashRef', builder => '_vars_builder', lazy => 1 );

extends 'MalariaGEN::AGV::Tool';

around BUILDARGS => sub {
   my $orig = shift;
   my $class = shift;
   $_{steps} = _process_steps($_{steps});
   $_{locals} = _process_locals($_{locals});
   $class->$orig(@_);
};

sub _process_locals {
  my $locals = shift;
  my @result = map {
    return () unless defined $_;
    return MalariaGEN::AGV::Workflow::Var->new(name => $_, type => "local") unless blessed($_);
    $_->isa("MalariaGEN::AGV::Workflow::Var") && $_->type eq "local" or die "wrong local class or type $_";
    return $_; 
  } @{$locals || []};
  return wantarray ? @result : \@result;
}

sub _steps_by_id_builder {
  my $steps = $_->steps;
  my %result = ();
  foreach my $step (@$steps) {
    die "multiple steps with the same id " . $step->id if (exists($result{$step->id}));
    $result{$step->id} =$step;
  }
  return \%result;
}

sub step {
  my ($self,$name) = @_;
  $name or die "you must specify and non empty step name";
  if ($name eq "0" or $name =~ /[1-9]\d*/) {
    my $steps = $self->steps;
    return $steps->[$name] if (scalar(@$steps) > $name);
  }
  return $self->_steps_by_id->{$name};
}

sub _vars_builder {
  my $self = shift;
  my ($locals,$inputs,$outputs,$steps) = ($self->locals,$self->inputs,$self->outputs,$self->steps);
  my @in_vars = map {
    return MalariaGEN::AGV::Workflow::Var->new(name => $_->name, type => 'in');
  } @{$inputs || []}; 
  my @out_vars = map {
    return MalariaGEN::AGV::Workflow::Var->new(name => $_->name, type => 'out');
  } @{$outputs || []};
  my @step_vars = map { _step_vars($_) } @{$steps || []}; 
  my @result = (@$locals,@in_vars,@out_vars,@step_vars);
  my %result = ();
  foreach my $var (@result) {
    my @names = $var->names;
    foreach my $name (@names) {
      die "more than one variable with the same name $name" if (exists ($result{$name})); 
      $result{$name} = $var;
    }
  }
  return wantarray ? %result : \%result;
}

sub _step_vars {
  my $step = shift;
  my @result = ();
  foreach my $in ($step->tool->inputs) {
    my $name = $in->name;
    push @result, MalariaGEN::AGV::Workflow::Var->(name => $name, step => $step, type => 'in');
  }
  foreach my $out ($step->tool->outputs) {
    my $name = $out->name;
    push @result, MalariaGEN::AGV::Workflow::Var->(names => $name, step => $step, type => 'out');
  }
  return wantarray ? @result : \@result;
}

sub _process_steps {
  my $steps = shift;
  return wantarray ? () : [] unless defined $steps;
  my $next_step_num = 0;
  my @result = map { 
    $_ = MalariaGEN::AGV::Workflow::Step->new(%$_) unless blessed($_);
    $_->isa('MalariaGEN::AGV::Workflow::Step') or die "step has wrong class " . ref($_);
    $_->_set_num($next_step_num++);
    return $_;
  } @{$steps || []};
  return wantarray ? @result : \@result;
}

sub _deps_builder {
  my $self = shift;
  my @res = ();
  my $steps = $self->steps;
  my $step_count = scalar(@$steps);
  foreach (my $i = 0; $i < $step_count; $i++) {
    my %step_deps = ();
    my $in_map = $steps->[$i]->inputs;
    foreach my $v (values %$in_map) {
      next unless ($v =~ /^\$:/);
      next unless $1 ne "0";
      $step_deps{$1}++;
    }
    push @res , [ keys(%step_deps) ];
  }
  return wantarray ? @res : \@res;
}

sub _rule_graph {
  my ($self) = @_;
}

sub TTS {
   return Text::Template->new(TYPE => 'STRING', SOURCE => shift);
}

sub script_template {
  my ($self,%args) = @_;
  my $steps = $self->steps;
  my @rules = $self->_step_rule(@$steps);
  my $header = ".PHONY: " . join(" ",@rules) . "\n";
  my %tfn = (); # tfn = temporal file name
  my $rules = "";
  foreach my $step (@$steps) {
     my $rule = $self->_step_rule($step);
     my @rule_deps = map { $self->_step_rule($_) } $self->_step_deps($step);
     my $sub_command = $self->_sub_command_template($step,\%tfn,%args);
     $rules .= "\n$rule: " . join(" ",@rule_deps) . "\n" . $sub_command . "\n";
  }
  return TTS($header . $rules);
}

sub _step_rule {
  my ($self,@steps) = @_;
  my @rules = map { "rule_" . $_->id  } @steps;
  return wantarray ? @rules : \@rules;
}

sub _sub_command_template {
 my ($self,$step,$tfn,%args) = @_;
 $step->resolve();

 my $intr = $step->interpreter(%args);
 my $result = ""; 
 if ($intr eq '$SHELL') {
   $result = $step->script(%args);  
 }
 else {
   my $tmp_provider = $args{job} || $args{engine} or die "cannot find a temporal file provider\n";
   my $sfile = $tmp_provider = $tmp_provider->tempfile;
   my $sfh = IO::File->new($sfile,"w");
   print $sfh $step->script(%args);
   $sfh->close();
   $result = $step->command(%args);
   $result =~ s/\$(\S+)/\$\($1\)/;
 }
 my @lines = split(/\n+/,$result);
 $result = join("\n",map { "\t" . $_  } @lines);
 return $result;
}

1;

__END__


{
  inputs => {},
  outputs => {},
  vars => { id1 => { type => 'file' } },
  steps => [
   { id => "start", tool => "xxx", inputs => { in => "$::id1" }, outputs => { out => "$out1", out2 => "$start::out2" }}
  ],
  
}
