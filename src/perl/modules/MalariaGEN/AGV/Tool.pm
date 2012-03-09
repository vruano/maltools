package MalariaGEN::AGV::Tool;

use strict;
use warnings;
use Scalar::Util qw(blessed);

use MalariaGEN::AGV::Job;
use MalariaGEN::AGV::ToolData;
use Text::Template;

sub new {
  my ($class,%data) = @_;
  my %self = %data;
  my $self = \%self;
  no strict 'refs';
  $self->{inputs} = ${$class . "::INPUTS"} unless exists($self->{inputs});
  $self->{outputs} = ${$class . "::OUTPUTS"} unless exists($self->{outputs});
  $self->{inputs} = $class->_process_data($self->{inputs},mode => 'in');
  $self->{outputs} = $class->_process_data($self->{outputs},mode => 'out');
  return bless $self,$class;
}

sub job {
  my ($self,%args) = @_;
  $args{cpu_ratio} = $self->calculate_cpu_ratio(%args) unless exists $args{cpu_ratio};
  $args{memory} = $self->calculate_memory(%args) unless exists $args{memory};
  $args{cpu_time} = $self->calculate_cpu_time(%args) unless exists $args{cpu_time};
  return MalariaGEN::AGV::Job->new(tool => $self,%args);
}

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 2000;
}

sub calculate_cpu_time {
  return 30 * 60;
}

sub name {
  my $self = shift;
  my $class = ref($self);
  $class =~ s/MalariaGEN::AGV::Tools//;
  my @cname = split(/::/,$class);
  @cname = grep { lc(substr($_,0,1)) . substr($_,1) } @cname;
  return join("-",@cname);
}


sub tool_by_name {
  my ($class,$name) = @_;
  my @cname = split(/\-+/,$name);
  my $sname = join("::",map { uc(substr($_,0,1)) . substr($_,1) } @cname);
  my $cls = "MalariaGEN::AGV::Tools::" . $sname;
  eval "require $cls" and return $cls->new();
  $sname = join("",map { uc(substr($_,0,1)) . substr($_,1) } @cname);
  $cls = "MalariaGEN::AGV::Tools::" . $sname;
  eval "require $cls" and return $cls->new();
  die "not such a tool '$name' and class '$sname'";
}

sub script {
  my ($self,%args) = @_;
  $args{T} = \$self unless exists $args{T};
  my $result = $self->script_template(%args)->fill_in(HASH => {%args}) or die $Text::Template::ERROR;
  return $result;
}

sub script_template {
  die "yet not implemented";
}

sub command {
  my ($self,%args) = @_;
  $args{T} = \$self unless exists $args{T};
  return $self->command_template(%args)->fill_in(HASH => {%args});
}


sub command_template {
  return Text::Template->new(TYPE => 'STRING', SOURCE => '{$T->interpreter} {$S}');
}

sub interpreter {
  die "yet not implemented";
}


sub input {
  my ($self,$name) = @_;
  my $result = $self->inputs()->{$name};
  return $result;
}

sub output {
  my ($self,$name) = @_;
  return $self->outputs()->{$name};
}

sub data {
  my ($self,%args) = @_;
  my $in = $args{as} || 'HASH';
  delete $args{as};
  my %data = %{ $self->_data('inputs',%args) };
  %data = (%data,% {$self->_data('outputs',%args) });
  if ($in eq "HASH") {
    return wantarray ? %data : \%data;
  }
  elsif ($in eq "ARRAY") {
    return wantarray ? values(%data) : [ values(%data) ];
  }
  elsif ($in eq "REF") {
    my @datas = values(%data);
    return $#datas < 0 ? undef : $datas[0];
  }
  else {
    die "cannot handle return type '$in'\n";
  }
  
}

sub inputs {
  my ($self,%args) = @_;
  return $self->_data('inputs',%args);
}

sub add_input {
  my ($self,@args) = @_;
  if (scalar(@args) == 1) {
    $self->_check_input_override($args[0]);
    $self->{inputs}->{$args[0]->name} = $args[0];
  }
  else {
    my %args = @args;
    my $in = MalariaGEN::AGV::ToolData->new(%args);  
    $self->_check_input_override($in);
    $self->{inputs}->{$in->name} = $in;
  }
}


sub _check_input_override {
   my ($self,$new_in) = @_;
   my $old_in = $self->input($new_in->name) or return;
   if ($old_in->type ne $new_in->type) {
     die "trying to override input " . $old_in->name . " with different type " . $new_in->type . "\n";
   }
}

sub outputs {
  my ($self,%args) = @_;
  return $self->_data('outputs',%args);
}


sub _process_data {
  my ($self,$data,%args) = @_;
  my %result = defined($data) ? %$data : ();
  foreach my $k (keys %$data) {
     my $v = $data->{$k};
     next if blessed($v);
     $v->{name} = $k;
     $v->{mode} = $args{mode} if !exists($v->{mode}) && exists($args{mode});
     $data->{$k} = MalariaGEN::AGV::ToolData->new(%$v);
  }
  return wantarray ? %$data : $data;
}

sub _data {
  my ($self,$name,%args) = @_;
  my %result = %{$self->{$name}};
  %result = map { ($_ => $result{$_}) } grep { 
      my $no_pass = 0;
      foreach my $k (keys %args) {
        if ($k eq "type") {
          my $t = $result{$_}->type;
          next if $t->extends_type($args{$k});
          $no_pass = 1;
          last;
        }
        else {
          next if $result{$_}->$k  eq $args{$k};
          $no_pass =1;
          last;
        }
      }
      !$no_pass;
  } keys(%result); 
  return wantarray ? %result : \%result;
}

1;
