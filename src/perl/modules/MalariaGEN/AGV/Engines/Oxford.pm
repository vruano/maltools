package MalariaGEN::AGV::Engines::Oxford;

use Moose;

extends 'MalariaGEN::AGV::Engine';

use IO::Pipe;
use MalariaGEN::AGV::Engines::Local;
use MalariaGEN::AGV::Tool;
use Cwd qw(getcwd realpath);
use File::Temp;
use File::Copy::Recursive qw(rcopy);
use MalariaGEN::AGV::Config qw(scratch_config);
use File::Spec::Functions qw(catfile file_name_is_absolute);
use File::Basename qw(basename);
use IO::File;
use File::Path qw(remove_tree make_path);
use POSIX;
use MalariaGEN::AGV::Engines::Scratch;

has '+name' => ( default => 'Oxford');
has 'root' => ( is => 'ro', isa => 'Str' , lazy => 1, default => sub { $_[0]->scratch()->tempdir });
has 'temp_root' => (is => 'ro', isa => 'Str', lazy => 1, default => sub { my $result = catfile($_[0]->root,"tmp"); make_path($result); return $result; });
has 'scratch' => (is => 'ro', isa => 'MalariaGEN::AGV::Engines::Scratch', 
    default => sub { scratch_config()->scratch_instance() });
has '_in_ge' => ( is => 'ro', isa  => 'Bool', init_arg => undef, default => sub { defined $ENV{JOB_ID} });

sub scratch_excluded {
  my ($self,$file) = @_;
  return 1 if $file eq "-";
  return 1 if $file =~ /^\s*\/dev/;
  return 1 if $file =~ /^\s*\/data\/malariagen/;
  return 0;
}

sub new_script_file {
  my ($self,$cwd) = @_;
  return $self->scratch->tempfile(template=>'script_XXXX',dir => $cwd);
}

sub new_wrapper_file {
  my ($self,$cwd) = @_;
  return $self->scratch->tempfile(template=>'ge_script_XXXX',dir => $cwd);
}

sub run_job {
  my ($self,$job,%options) = @_;
  
  if ($self->_in_ge) { 
    my $local = MalariaGEN::AGV::Engines::Local->new(root => $self->root, temp_root => $self->temp_root);
    return $local->run_job($job,%options);
  }

  # prepare and copy inputs to scratch area, prepare outputs (creates directory structure).
  my $mvs = {};
  my $uses_standard_io = $self->_prepare_inputs($job);
  $uses_standard_io += $self->_prepare_outputs($job,$mvs);

  # composes the executing launching or executing command. 
  my $qsub_cmd = $self->qsub_command($job, %options, uses_standard_io => $uses_standard_io);

  # runs the command: 
  print STDERR "Sending job for execution to SGE ...\n";
  print STDERR "\tCommand: $qsub_cmd\n";
  my $retries = $self->running_options->{retries};
  my $code = 0;
  $| = 1;
  do {
    my $qsub_output;
    $qsub_output = IO::Pipe->new();
    my $redirect = $uses_standard_io ? "" : "2>&1";
    $qsub_output->reader("$qsub_cmd $redirect || echo \"<<<<ERROR>>>>\" \$?");
    while (my $output_line = <$qsub_output>) {
       if ($output_line =~ /^<<<<ERROR>>>> (\d+)/) {
         $code = $1;
         next;
       }
       elsif ($uses_standard_io) {
         print STDOUT $output_line;
       }
       else {
         print STDERR "\tSGEOUT: $output_line";
       }
    }
    print STDERR "Done" . ($code != 0 ? " with errors!!!" : ".") . "\n";
    $qsub_output->close(); 
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
sub qsub_command {
  my ($self,$job,%options) = @_;
  $job->_set_engine($self);
  my $tool = $job->tool;
  my $cwd = $self->scratch->tempdir(dir => $self->root);
  $job->_set_wd($cwd);
  my $uses_standard_io = $options{uses_standard_io};
  my $sfn = $self->new_script_file($cwd);
  my $wfn = $self->new_wrapper_file($cwd);
  my $sfn_simple = basename($sfn);
  my $procs = $options{procs} || 1;
  my %job_args = (
     J => \$job, T => \$tool, E => \$self, D => $cwd,
     S => \$sfn_simple,  MAKE => "lsmake" , SH => "tcsh", P => $procs); 
##
  my $cmd = $job->command(%job_args);
  # translate standard progrma place-holders.
  $cmd =~ s/^\$MAKE/make/;
  $cmd =~ s/^\$SHELL/tcsh/;

  my @management_arguments = $self->_qsub_management_arguments($job,$uses_standard_io);
  my @memory_arguments = $self->_qsub_command_memory_arguments($job);
  my @temp_arguments = $self->_qsub_command_temp_arguments($job);
  my @cpu_arguments = $self->_qsub_command_cpu_arguments($job);
  
  my $sfh = IO::File->new($sfn,"w");
  my $script = $job->script(%job_args);
  print $sfh $script;
  $sfh->close();
  
  my $result = "";
  if ($cmd =~ /^make (.*)/) {
    $result = join(" ","make","-C",$job->wd,$1);
  }
  else {
    $result = join(" ",$uses_standard_io ? "qsh" : "qsub -o output.txt" , "-e", "error.txt",
      ($uses_standard_io ? () : ("-sync y")),"-N",$job->name,
      "-wd",$job->wd,"-V",
      @management_arguments,
      @memory_arguments,
      @cpu_arguments,
      $wfn);
  }

  my $wfh = IO::File->new($wfn,"w");
  print $wfh $cmd,"\n";
  $wfh->close();

  return $result;
}

sub _qsub_command_memory_arguments {
  my ($self,$job) = @_;
  my $memory = $self->running_options()->{memory} || $job->memory;
  return ("-l mem_free=${memory}");
}

sub _qsub_command_temp_arguments {
  my ($self,$job) = @_;
  return ();
}

sub _qsub_command_cpu_arguments {
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
  return ("-pe smp $cpu_count");
}

# Wall time thresholds in seconds.

our $THR_FACTOR = 0.75; # multipled to the hard limit to get a risk threshold.
our $NORMAL_THR = floor(2* 60 * 60 * $THR_FACTOR); 
our $MEMORY_THR = 10000;

sub _qsub_management_arguments {
  my ($self,$job,$uses_standard_io) = @_;
  if ($self->running_options()->{queue}) {
    return ('-q',$self->running_options()->{queue});
  }
  my $cpu_ratio = $job->cpu_ratio;
  my $cpu_count = floor($cpu_ratio);
  $cpu_count++ if $cpu_count == 0 || $cpu_ratio - $cpu_count >= 0.5;
  my $mem = $job->memory;
  my $wall_time = $job->wall_time;
  if ($mem > $MEMORY_THR) {
    return('-q','hugemem.q');
  }
  elsif ($wall_time <= $NORMAL_THR) { 
    return ('-q','normal.q');
  }
  else {
    return ('-q','long.q');
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

sub _prepare_inputs {
  my ($self,$job) = @_;
  my $uses_standard_io = 0;
  my $tool = $job->tool;
  my $message_shown = 0;
  my %file_inputs = map { ($_->name => $_ ) } grep { $_->type->has_files } values %{$tool->inputs()};
  my $scratch = $self->scratch();
  foreach my $in (values(%file_inputs)) {
    my $value = $job->input($in->name) or next;
    my @files = $in->type->file_names($value);
    my @real_files = map { $_ eq "-" ? "-" : (realpath($_) || $_) } @files;
    foreach (my $i = 0; $i < scalar(@files); $i++) {
       if ($files[$i] ne $real_files[$i]) {
         $value = $in->type->relocate_file_on_value($value,$files[$i],$real_files[$i]);
         @files = $in->type->file_names($value);
         $job->set_input($in->name => $value);
         last;
       }
    }
    my $svalue = undef;
    foreach my $ofile (@files) {
      if ($ofile eq "-") {
        die "std error in-out not allowed in Oxford engine";
        $uses_standard_io = 1;
        next;
      }
      my $file = $ofile;
      my $sfile = $self->scratch_excluded($file)  ? $file : 
         $scratch->to_scratch(file_or_dir => $file, make_path => 1, $svalue ? (close_to => (ref($svalue) eq "ARRAY" ? $$svalue[0] : $svalue)) : ());
      
      print STDERR "Copying inputs to Scratch area ...\n" unless $message_shown++;
      print STDERR "\t$file >>>> $sfile\n" unless $file eq $sfile;
      next if $sfile eq $ofile;
      $svalue = $svalue || $in->type->relocate_file_on_value($value,$ofile,$sfile);
    }
    $job->set_input($in->name => $svalue) if defined $svalue;
  }
  print STDERR "... Done.\n" if $message_shown;
  return $uses_standard_io;
}

sub _prepare_outputs {
  my ($self,$job,$mvs) = @_;
  my $uses_standard_io = 0;
  my $tool = $job->tool;
  my %file_inputs = $tool->outputs(type => "file");
  my $scratch = $self->scratch();
  foreach my $in (values(%file_inputs)) {
    my $value = $job->output($in->name) or next;
    my @files = $in->type->file_names($value);
    my @real_files = map { $_ eq "-" ? "-" : (realpath($_) || $_) } @files;
    foreach (my $i = 0; $i < scalar(@files); $i++) {
       if ($files[$i] ne $real_files[$i]) {
         $value = $in->type->relocate_file_on_value($value,$files[$i],$real_files[$i]);
         @files = $in->type->file_names($value);
         $job->set_output($in->name => $value);
         last;
       }
    }
    my $svalue  = undef;
    foreach my $ofile (@files) {
      if ($ofile eq "-") {
        die "std error in-out not allowed in Oxford engine";
        $uses_standard_io = 1;
        next;
      }
      my $file = $ofile; 
      my $sfile = $self->scratch_excluded($file) ? $file : 
           $scratch->to_scratch(file_or_dir => $file, make_path => 1, $svalue ? (close_to => (ref($svalue) eq "ARRAY" ? $$svalue[0] : $svalue)) : ());
      $mvs->{$sfile} = $file;
      next if $sfile eq $ofile;
      $svalue = $svalue || $in->type->relocate_value($value,$ofile,$sfile);
    }
    $job->set_output($in->name => $svalue) if defined $svalue;
  }
  return $uses_standard_io;
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
      print STDERR "-e $sfile: " . (-e $sfile) . " stats: " . `stat $sfile`, "\n";
      print STDERR "\t$ofile <-X- $sfile !!!! MISSING !!!!\n";
    }
  }
  print STDERR "... Done" . ($error_count ? " with $error_count error(s)" : "") ."\n" if $message_shown;
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
