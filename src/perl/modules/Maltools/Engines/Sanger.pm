package Maltools::Engines::Sanger;

use Moose;
extends 'Maltools::Engine';

use IO::Pipe;
use Maltools::Engines::Local;
use Maltools::Engines::Scratch::Array;
use Maltools::Tool;
use Cwd qw(getcwd realpath);
use File::Temp;
use File::Copy::Recursive qw(rcopy);
use Maltools::Config qw(sanger_config);
use File::Spec::Functions qw(catfile file_name_is_absolute);
use File::Basename qw(basename);
use IO::File;
use File::Path qw(remove_tree make_path);
use POSIX;

has '+name' => ( default => 'Sanger');
has 'root' => ( is => 'ro', isa => 'Str' , lazy => 1, default => sub { $_[0]->scratch()->tempdir });
has 'temp_root' => (is => 'ro', isa => 'Str', lazy => 1, default => sub { my $result = catfile($_[0]->root,"tmp"); make_path($result); return $result; });
has 'scratch' => (is => 'ro', isa => 'Maltools::Engines::ScratchArray', 
    default => sub { Maltools::Engines::Sanger::Scratch::Array->new() });
has '_sc' => ( is => 'ro' , lazy => 1, default => sub { sanger_config });
has '_in_lsf' => ( is => 'ro', isa  => 'Bool', init_arg => undef, default => sub { defined $ENV{LSB_JOBID} });

sub scratch_excluded {
  my ($self,$file) = @_;
  return 1 if $file eq "-";
  return 1 if $file =~ /^\s*\/dev/;
  return 1 if $file =~ /^\s*\/lustre/;
  return 0;
}

sub new_script_file {
  my ($self,$cwd) = @_;
  return $self->scratch->tempfile(template=>'script_XXXX',dir => $cwd);
}

sub run_job {
  my ($self,$job,%options) = @_;
  
  # delegates into local engine if this is already a LSF job
  if ($self->_in_lsf) { 
    my $local = Maltools::Engines::Local->new(root => $self->root, temp_root => $self->temp_root);
    return $local->run_job($job,%options);
  }

  # prepare and copy inputs to scratch area, prepare outputs (creates directory structure).
  my $mvs = {};
  $self->_prepare_data($job,$mvs);

  # composes the executing launching or executing command. 
  my $bsub_cmd = $self->bsub_command($job, %options);

  # runs the command: 
  print STDERR "Sending job for execution to LSF ...\n";
  print STDERR "\tCommand: $bsub_cmd\n\tLSFOUT: ";
  my $retries = $self->running_options->{retries};
  my $code = 0;
  my $uses_standard_io = $job->uses_standard_io;
  local $ENV{TMPDIR} = $job->wd;
  do {
    my $bsub_output;
    $bsub_output = IO::Pipe->new();
    $bsub_output->reader(" ($bsub_cmd > /dev/stdout) >& /dev/stderr || echo \"<<<<ERROR>>>>\" \$?");
    while (my $output_line = <$bsub_output>) {
       if ($output_line =~ /^<<<<ERROR>>>> (\d+)/) {
         $code = $1;
         next;
       }
       elsif ($uses_standard_io) {
         print STDOUT $output_line if $uses_standard_io;
       }
       else {
         print STDERR "\tLSFOUT: $output_line";
       }
    }
    print STDERR "Done" . ($code != 0 ? " with errors!!!" : ".") . "\n";
    $bsub_output->close(); 
  } while ($code != 0 && $retries-- > 0);
  $job->return_code($code);
  eval { $self->_export_outputs($job,$mvs); };
  my $return = ($code == 0) && !$@;
  unless ($return) {
    print STDERR  "Check '" . catfile($job->wd, 'error.txt') . "' content for details\n";
    $job->error_message($@);
  }
  else {
    #eval { remove_tree($job->wd); };
  }
  return $return;
}

# prepare the job working directory, create the script file, composes and return the executing command.
sub bsub_command {
  my ($self,$job,%options) = @_;
  $job->_set_engine($self);
  my $tool = $job->tool;
  my $cwd = $self->scratch->tempdir(dir => $self->root);
  $job->_set_wd($cwd);
  my $uses_standard_io = $job->uses_standard_io;
  my $sfn = $self->new_script_file($cwd);
  my $sfn_simple = basename($sfn);
  my $procs = $options{procs} || 1;
  my %job_args = (
     J => \$job, T => \$tool, E => \$self, D => $cwd,
     S => \$sfn_simple,  MAKE => "lsmake" , SH => "tcsh", P => $procs); 
##
  my $cmd = $job->command(%job_args);
  # translate standard progrma place-holders.
  $cmd =~ s/^\$MAKE/lsmake/;
  $cmd =~ s/^\$SHELL/tcsh/;

  my @management_arguments = $self->_bsub_management_arguments($job,$uses_standard_io);
  my @memory_arguments = $self->_bsub_command_memory_arguments($job);
  my @temp_arguments = $self->_bsub_command_temp_arguments($job);
  my @cpu_arguments = $self->_bsub_command_cpu_arguments($job);
  
  my $sfh = IO::File->new($sfn,"w");
  my $script = $job->script(%job_args);
  print $sfh $script;
  $sfh->close();
  
  my $result = "";
  if ($cmd =~ /^lsmake (.*)/) {
    $result = join(" ","lsmake","-C",$job->wd,$1);
  }
  else {
    $result = join(" ","bsub",
      ($uses_standard_io ? "-I" : "-K -o output.txt") , "-e", "error.txt",
      "-cwd",$job->wd,
      @management_arguments,
      @memory_arguments,
      @temp_arguments,
      @cpu_arguments,
      $cmd);
  }
  $result .= " | bsub-remove" if $uses_standard_io;
  return $result;
}

sub _bsub_command_memory_arguments {
  my ($self,$job) = @_;
  my $memory = $self->running_options()->{memory} || $job->memory;
  return ("-M",$memory * 1000,"-R","'select[mem>${memory}] rusage[mem=${memory}]'");
}

sub _bsub_command_temp_arguments {
  my ($self,$job) = @_;
  my $tmp_required = $job->tmp_required();
  return () unless $tmp_required; 
  return ("-R","'select[tmp>$tmp_required] rusage[tmp=$tmp_required]'");
}

sub _bsub_command_cpu_arguments {
  my ($self,$job,) = @_;
  my $cpu_ratio = $self->running_options()->{load} || $job->cpu_ratio;
  my $cpu_count = $self->running_options()->{procs} || floor($cpu_ratio);
  unless ($self->running_options()->{procs}) {
    $cpu_count++ if $cpu_count == 0 || $cpu_ratio - $cpu_count >= 0.5;
  }
  # No parallel jobs:
  return () unless $cpu_count > 1;
  # Multi-threaded jobs that cannot run across hosts:
  $job->cpu_count($cpu_count);
  return ("-n",$cpu_count,"-R","'span[hosts=1] span[ptile=$cpu_count]'") unless $job->is_distributed;
  # Others:
  return ("-n",$cpu_count);
}

# Wall time thresholds in seconds.
our $THR_FACTOR = 0.75; # multipled to the hard limit to get a risk threshold.
our $HHOUR_THR = floor(30 * 60 * $THR_FACTOR);
our $HOUR_THR = floor(60 * 60 * $THR_FACTOR);
our $HALF_DAY_THR = floor(60 * 60 * 12 * $THR_FACTOR);
our $DAY_THR = 2 * $HALF_DAY_THR;
our $TWO_DAYS_THR = 2 * $DAY_THR;

sub _bsub_management_arguments {
  my ($self,$job,$uses_standard_io) = @_;
  if ($self->running_options()->{queue}) {
    return ('-q',$self->running_options()->{queue});
  }
  my $cpu_ratio = $job->cpu_ratio;
  my $cpu_count = floor($cpu_ratio);
  $cpu_count++ if $cpu_count == 0 || $cpu_ratio - $cpu_count >= 0.5;
  my $mem = $job->memory;
  my $wall_time = $job->wall_time;
  # Decrease the number of CPUs for non distributable jobs to a maximum of 12.
  if (! $job->is_distributed && $cpu_count > 12) {
    $wall_time *= ($cpu_count/12);
    $cpu_count = 12;
  }
  if ($wall_time <= $HHOUR_THR && $cpu_count <= 8 && $mem <= 20000) {
    return ('-q','small');
  }
  if ($wall_time <= $HOUR_THR && $cpu_count <= 8 && $mem <= 20000) { 
    return ('-q','normal');
  }
  elsif ($wall_time <= $TWO_DAYS_THR && $cpu_count <= 12 && $mem <= 20000) {
    return ('-q','long');
  }
  elsif (!$uses_standard_io && $mem <= 20000) {
    return ('-q','basement');
  }
  elsif ($wall_time < $TWO_DAYS_THR) {
    return ('-q','hugemem');
  }
  else {
    die "there is no good queue for required memory ${mem}MB, cpu count $cpu_count, wall time estimate " . ($wall_time / 60 * 60 ). " minutes, a " . ($uses_standard_io ? "" : "non-") . "interactive job that is " 
    . ($job->is_distributed() ? "": "not") . " distributable"; 
  }
}

sub tempdir {
  my ($self,@args) = @_;
  $self->scratch->tempdir(dir => $self->temp_root,@args);
}

sub tempfile {
  my ($self,%args) = @_;
  $self->scratch->tempfile(dir => $self->temp_root,%args);
}

sub _real_path {
  my ($self,$job,$in,$value) = @_;
  unless (defined $value) {
    $value = $job->data($in->name) or return undef;
  }
  if (ref($value) eq "ARRAY") {
     return [ map { $self->_real_path($job,$in,$_) } @$value ];
  }
  my @files = $in->type->file_names($value);
  my @real_files = map { realpath($_) || $_ } @files;
  foreach (my $i = 0; $i < scalar(@files); $i++) {
     if ($files[$i] ne $real_files[$i]) {
        $value = $in->type->relocate_file_on_value($value,$files[$i],$real_files[$i]);
        last;
     }
  }
  return $value;
}

sub _scratch_path {
  my ($self,$job,$mvs,$in,$value) = @_;
  unless (defined $value) {
    $value = $job->data($in->name) or return undef;
  }
  if (ref($value) eq "ARRAY") {
    return [ map { $self->_scratch_path($job,$mvs,$in,$_) } @$value ];
  }
  my @files = $in->type->file_names($value);
  return undef unless $#files >= 0;
  my $scratch_value = undef;
  my $scratch = $self->scratch();
  foreach my $file (@files) {
    if ($self->scratch_excluded($file)) {
       next;
    }
    elsif (defined $scratch_value) {
       my $scratch_file = $scratch->to_scratch(file_or_dir => $file, make_path => 1, close_to => $scratch_value);
       $$mvs{$scratch_file} = $file if defined $mvs;
       print STDERR "\t$file <<<>>> $scratch_file\n";
    }
    else {
       my $scratch_file = $scratch->to_scratch(file_or_dir => $file, make_path => 1, close_to => $scratch_value);
       $$mvs{$scratch_file} = $file if defined $mvs;
       print STDERR "\t$file <<<>>> $scratch_file\n";
       $scratch_value = $in->type->relocate_file_on_value($value,$file,$scratch_file);
    }
  }
  return $scratch_value || $value; 
}

sub _prepare_data {
  my ($self,$job,$mvs) = @_;
  my $tool = $job->tool;
  my %file_data = map { ($_->name => $_ ) } grep { $_->type->has_files } values %{$tool->data()};
  
  print STDERR "Mapping data to Scratch area ...\n";
  foreach my $data (values(%file_data)) {
    my $real_value = $self->_real_path($job,$data);
    $job->set_data($data->name => $real_value);
    my $scratch_value = $self->_scratch_path($job,($data->mode =~ /out/ ? $mvs : undef),$data);
    $job->set_data($data->name => $scratch_value);
  }
  print STDERR "... Done.\n";
}


# copy the outputs back out in the warehouse, nfs or whatever.
sub _export_outputs {
  my ($self,$job,$mvs) = @_;
  my $message_shown = 0;
  my $error_count = 0;
  foreach my $sfile (keys %$mvs) {
    my $ofile = $mvs->{$sfile}; 
    print STDERR "Copying outputs back from Scratch area ...\n" unless $message_shown++;
    if (-e $sfile) {
      if ($sfile ne $ofile) {
        print STDERR "\t$ofile <<<< $sfile\n";
        rcopy($sfile,$ofile);
      }
    }
    else {
      $error_count++;
      print STDERR "\t$ofile <-X- $sfile !!!! MISSING !!!!\n";
    }
  }
  print STDERR "... Done" . ($error_count ? "with $error_count error(s)" : "") ."\n" if $message_shown;
}

sub import_resource {
  my ($self,$what) = @_;
  if (-e $what) {
    return $self->scratch->to_scratch(file_or_dir => $what);
  }
  else {
    return $what;
  }
}

1;
