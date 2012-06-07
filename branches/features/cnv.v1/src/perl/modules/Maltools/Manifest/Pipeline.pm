package Maltools::Manifest::Pipeline;

use Moose;

require Maltools::Config;
require Maltools::Manifest::Config;
require Maltools::Manifest::Pipeline::Sample;
require Maltools::DataFreeze;

extends 'Maltools::Manifest';

has '_config' => (is => 'ro' , isa => 'Maltools::Manifest::Config', default => sub { Maltools::Config->new() });

sub name {
  return $_[0]->get('name');
}

sub get_path {
  my $self = shift;
  my @names = @_;
  my $paths = $self->get('paths') || {};
  my @result = map { $paths->{$_} } @names;
  if ($#names == 0 && !wantarray) {
    return $result[0];
  }
  else {
    return wantarray ? @result : \@result;
  }
}

override '_set_variable_hash' => sub {
    my $self = shift;
    my $hash = shift;
    if (!exists $hash->{basedir}) {
        $self->{basedir} = ".";
    }
    return super();
};

sub get_parameter {
  my $self = shift;
  return $self->get('parameters',@_) || undef;
}

sub get_sample_ids {
  my $self = shift;
  my $samples = $self->_get_samples;
  
  if (ref($samples) eq "ARRAY") {
    return @$samples;
  }
  elsif (ref($samples) eq "HASH") {
    return keys %$samples;
  }
  else {
    die "invalid sample section neither and array nor hash";
  }
}
   
sub get_sample_names {
  my $self = shift;
  my $samples = $self->_get_samples;
 
  if (ref($samples) eq "ARRAY") {
    return @$samples;
  }
  elsif (ref($samples) eq "HASH") {
    return map { $samples->{$_}->{oxfordCode} || $_ }  keys %$samples;
  }
  else {
    die "invalid sample section neither and array nor hash";
  }
}

sub get_reference {
  my $self = shift;
  my $sequence_file = $self->get_path('reference') || $self->_config()->get_reference_sequences(); 
  return Maltools::Reference->new(sequence_file => $sequence_file);
}

sub get_candidate_snps {
  my $self = shift;
  return $self->get_path('candidateSnps') || $self->_config()->get('reference/candidateSNPs/list');
}

sub get_sample {
  my $self = shift;
  my $samples = $self->_get_samples;
  if (ref($samples) eq "ARRAY") {
    $samples = { map { $_ => { oxfordCode => $_ } }  @$samples };
  }
  if (ref($samples) ne "HASH") {
    die "invalid sample section neither and array nor hash";
  }
  my @result = map { Maltools::Manifest::Pipeline::Sample->new(manifest => $self, %{$samples->{$_}}) } @_;
  if ($#_ == 0) {
    return $#result == -1 ? undef : $result[0];
  }
  else {
    return wantarray ? @result : \@result;
  }
}
  
sub get_all_samples {
  my $self = shift;
  my $result = scalar($self->get_sample(@_));
  return wantarray ? @$result : $result;
}

sub _get_samples {
    my $self = shift;
    return $self->get('samples') || {};
}

sub _set_samples {
    my $self = shift;
    $self->get()->{samples} = $_[0];
}

sub auto_complete {
    my $self = shift;
    my %options = @_;
    my $freeze = $options{freeze} || $self->get_variable('freeze');
    $freeze = Maltools::DataFreeze->new(file => $freeze) if $freeze && ref($freeze) eq "";
    my $samtrak = $options{samtrak} || $self->get_variable('samtrak');
    if ($freeze) {
        $self->_auto_complete_from_freeze($freeze);
    }
    elsif ($samtrak) {
        $self->_auto_complete_from_samtrak($samtrak);
    }
}

sub _auto_complete_from_samtrak {
    my $self = shift;
    my $samtrak = shift;
    my $samples = $self->_get_samples();
    if (ref($samples) eq "ARRAY") {
        $samples = { map { $_ => { oxfordCode => $_ }} @$samples }; 
    }
    elsif (ref($samples) ne "HASH") {
        die "invalid samples specification neither a ARRAY or HASH";
    }
    
    my $sm = $samtrak->manager_class_of('Sample');
    my $fm = $samtrak->manager_class_of('File');
    
    my %missing = ();
    my %found = ();
    foreach my $key (keys %$samples) {
        my $val = $samples->{$key};
        my $ox_code = $val->{oxfordCode};
        my $sample_id = $val->{solarisId};
        my $sample = ($sample_id ? $sm->get_by_id($sample_id) : undef) ||
          ($ox_code ? $sm->get_by_ox_code($ox_code) : undef);
        $ox_code = $sample->code();
        $sample_id = $sample->id;
    
        unless ($sample) {
            warn "sample with Oxford code '$ox_code' cannot be found!!!";
            $missing{$key} = $val;
            next;
        }

        my $files = $fm->get_files(query => [ 'samples.id' => $sample->id ,
               'lanes.use_flag' => 1,
                name => { like => '%_nonhuman%' }],
              with_objects => [qw(samples lanes)], multi_many_ok => 1) || [];
        $files = [ $files ] if blessed($files);
        
        unless ($#$files >= 0) {
          warn "sample with Oxford code '$ox_code' does not have usable sequencing data available";
          $missing{$key} = $val;
        }
        
        my @files = map { _map_freeze_file($_) } @$files;
        $found{$key} = {
          solarisId => $sample_id,
          oxfordCode => $ox_code,
          dataSources => \@files
        };
    }
    $self->_set_samples(\%found);
    $self->get()->{missingSamples} = \%missing;
}

sub _auto_complete_from_freeze {
    my $self = shift;
    my $freeze = shift;
    my $samples = $self->_get_samples();
    if (ref($samples) eq "ARRAY") {
        $samples = { map { $_ => { oxfordCode => $_ }} @$samples }; 
    }
    elsif (ref($samples) ne "HASH") {
        die "invalid samples specification neither a ARRAY or HASH";
    }
    
    my %missing = ();
    my %found = ();
    foreach my $key (keys %$samples) {
        my $val = $samples->{$key};
        my $ox_code = $val->{oxfordCode} or die "some samples do not indicate their name/ox_code!";
        my $f_data = $freeze->find_sample($ox_code);
        unless ($f_data) {
            warn "sample with Oxford code '$ox_code' missing in freeze, won't be proccessed";
            $missing{$key} = $val;
            next;    
        }
        unless ($f_data->{lanes} && scalar(@{$f_data->{lanes}}) > 0) {
            warn "sample with Oxford code '$ox_code' missing lane information and won't be processed";
            $missing{$key} = $val;
            next
        }
        my @usable_lanes = grep { $_->{files} && scalar(@{$_->{files}}) > 0 } @{$f_data->{lanes}};
        unless ($#usable_lanes >= 0) {
            $missing{$key} = $val;
            warn "sample with Oxford code '$ox_code' does not have lanes with files and won't be processed";
            next;        
        }
        my @files = (map { @{$_->{files}} } @{$f_data->{lanes}});
        @files =  grep {$_->{name} =~ /nonhuman/ } @files;
        unless ($#files >= 0) {
            $missing{$key} = $val;
            warn "sample with Oxford code '$ox_code' does not have usable files with non-human data and won't be processed";
            next;        
        }
        @files = map { _map_freeze_file($_) } @files;
        $found{$key} = {
          solarisId => $f_data->{sample_id},
          oxfordCode => $ox_code,
          dataSources => \@files
        };
    }
    $self->_set_samples(\%found);
    $self->get()->{missingSamples} = \%missing;
}
    
sub _map_freeze_file {
  my $freeze_file = shift;
  my $url = "";
  my $filesystem = $freeze_file->{filesystem} || "file";
  if ($filesystem =~ /irods/) {
    $url .= "irods://";
  }
  else {
    $url .= "file://";
  }
  my $full_name = ($freeze_file->{path} || "/") . ($freeze_file->{name} || "");
  $full_name =~ s/\/+/\//g;
  $url .= $full_name;
  my $result = { url => $url };
  $result->{md5} = $freeze_file->{md5} if exists $freeze_file->{md5};
  return $result;
}


1;
