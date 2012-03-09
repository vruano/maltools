package MalariaGEN::AGV::Commands::Testbed;

use Moose;
extends 'MalariaGEN::AGV::Command';

use Template;
use JSON::XS;

use POSIX;
use File::Spec::Functions qw(catfile);
use MalariaGEN::AGV::Config qw(reference_config);
use MalariaGEN::AGV::Config::Reader;
use MalariaGEN::AGV::Commands::Reference;
use MalariaGEN::AGV::Reference;
use MalariaGEN::AGV::DataFreeze;
use Getopt::Long qw(:config pass_through);
use Samtrak::Site;
use Scalar::Util qw(looks_like_number);
use Cwd qw(realpath);
use IO::File;

has freeze => (is => 'rw', isa => 'Maybe[MalariaGEN::AGV::DataFreeze]', default => undef);


sub help_summary {
   return "manages and execute testbeds to evaluate pipelines";
}

sub help_text {
   return "testbed:\n\t" .
          $_[0]->help_summary ."\n".
          "Syntaxis:\n" .
          $_[0]->cl_name . " testbed [sub-command] [arguments]\n" .
          "Subcommands:\n" .
          "\tcreate <descr.json> <directory> creates the testbed given " . 
          "\t\t\t" . "its description json document and the name of the destination directory\n";
}

sub reference_command {
  my ($self,$config) = @_;
  my $refCommand = MalariaGEN::AGV::Commands::Reference->new();
  if (exists $config->{resources} && exists $config->{resources}->{reference}) {
     $refCommand->reference_file($config->{resources}->{reference});
  }
  if (ref($config->{sequences}) eq "ARRAY") {
     $refCommand->sequences($config->{sequences});
  }
  return $refCommand;
}

sub reference {
  my ($self,$config) = @_;
  my %args =  ();
  unless (exists $config->{resources} && exists $config->{resources}->{reference}) {
    return MalariaGEN::AGV::Config->new()->get_reference();
  }
  if (exists $config->{resources} && exists $config->{resources}->{reference}) {
    $args{sequence_file} = $config->{resources}->{reference}
  }
  if (exists $config->{resources} && exists $config->{resources}->{referenceIndex}) {
    $args{index_file} = $config->{resources}->{referenceIndex};
  }

  return MalariaGEN::AGV::Reference->new(%args);
}

sub execute {
  my ($self,$samtrak_site,$params) = @_;
  local @ARGV = @{$params->{arguments}};
  my $subcommand = shift @ARGV;

  if ($subcommand eq "create") {
    return $self->create(@ARGV);
  }
  elsif ($subcommand eq "run") {
    return $self->run_testbed(@ARGV);
  }
  else {
    return $self->error_return("unknown subcommand '$subcommand'");
  }
}

sub run_testbed {
  my ($self,@args) = @_;
  local @ARGV = (@args);
  my $threads = 200;
  my $testbed = ".";
  my $commands_only = 0;
  my @make_options = ();
  GetOptions("n=i" => \$threads,"testbed|dir=s" => \$testbed,"commands-only|co!" => \$commands_only,"make-options|mo=s@" => \@make_options);
  $testbed = shift(@ARGV) unless defined $testbed;
  $testbed = realpath($testbed);
  -f catfile($testbed,"MANIFEST.json") or return $self->error_return("there is no testbed deployed in directory '$testbed'");
  -f catfile($testbed,"Makefile") or return $self->error_return("testbed at '$testbed' lacks Makefile");
  my @targets = map { $_ eq "all" ? $_ : catfile($_,"DONE.ts") } @ARGV;
  push @targets, 'all' unless $#targets >= 0;
  my $co_flag = $commands_only ? "-n" : "";
  my $make_options = join(" ",@make_options);
  system("cd $testbed; make -k -j $threads $co_flag $make_options " . join(" ",@targets));
  if ($?) {
    return $self->error_return("error ocurred while running the test-bed");
  }
  else {
    return $self->ok_return();
  }
}

sub create {
  my ($self,@args) = @_;
  my $config = undef;
  my $dir = undef;
  my $freeze = undef;
  local @ARGV = @args;
  GetOptions("manifest|m=s" => \$config, "dir=s" => \$dir, "freeze=s" => \$freeze);
  $config ||= shift;
  $dir ||= shift;
  $self->freeze(MalariaGEN::AGV::DataFreeze->new(file => $freeze)) if defined $freeze;
  
  $config or return $self->error_return("you must indicate a manifest file");
  -e $config or return $self->error_return("cannot find configuration file '$config'");
  -f $config or return $self->error_return("the configuration file provided '$config' is not a regular file");
  $dir or return $self->error_return("you must indicate a testbed destination directory");
  return $self->error_return("the destination '$dir' does not seem to be a directory") if -e $dir && ! -d $dir;
  mkdir($dir) unless -e $dir;  
  mkdir(catfile($dir,'intervals')) unless -e catfile($dir,'intervals');
  
  my $template = Template->new({
     DEBUG => 1,
     EVAL_PERL => 1,
     ABSOLUTE => 1,
     INCLUDE_PATH => MalariaGEN::AGV::Config->template_dir("crosses"),
  });

  my $cross = $self->read_cross_config($config) or return $self->error_return("could not read cross configuration");
  ref($cross) eq "HASH" or return $self->error_return($cross);

  my %vars = (
    cross => $cross,
  );

  my $make_content = "";
  $template->process("Makefile.tt",\%vars,\$make_content) or return $self->error_return("Error in template: " . $template->error);
  $self->write_into_file(catfile($dir,'Makefile'),\$make_content);
  foreach my $mapping_name (keys %{$cross->{mappings}}) {
    mkdir(catfile($dir,$mapping_name));
    mkdir(catfile($dir,$mapping_name,'vf-out'));
    $vars{mapping} = $cross->{mappings}->{$mapping_name};
    $vars{mapping_name} = $mapping_name;
    $vars{mapping_dir} = catfile($dir,$mapping_name);
    my $vf_manifest_content = "";
    $template->process("vf-manifest.tt",\%vars,\$vf_manifest_content) or return $self->error_return("Error in template: " . $template->error);
    $self->write_into_file(catfile($dir,$mapping_name,'vf-manifest.json'),\$vf_manifest_content);
  }

  my $cfh = IO::File->new(catfile($dir,"MANIFEST.json"),'w');
  print $cfh JSON::XS->new->pretty->encode($cross);
  $cfh->close();
  return $self->ok_return();
}

sub read_cross_config {
  my ($self,$filename) = @_;
  my $reader = MalariaGEN::AGV::Config::Reader->new();
  $reader->set_vars(%{MalariaGEN::AGV::Config->get('vars')});
  $reader->load($filename);
  
  return $self->complete_config($reader->get());
}

sub complete_config {
  my ($self,$config) = @_;

  $config = $self->complete_config_sequences($config);
  $config = $self->complete_config_intervals($config);
  $config = $self->complete_config_samples($config);
  return $config;
}


sub complete_config_intervals {
   my ($self,$config) = @_;
   if (exists $config->{intervals}) {
     $config->{intervals} = $self->process_config_intervals($config->{intervals},$config);
   }
   else {
     $config->{intervals} = $self->process_config_intervals($self->default_config_intervals($config),$config);
   }
   return $config;
}

sub complete_config_samples {
  my ($self,$config)  = @_;
  if (exists $config->{samples}) {
    $config->{samples} = $self->process_config_samples($config->{samples},$config);
  }
  else {
    die "you must indicate targeted samples";
  }
  return $config;
}

sub process_config_samples {
  my ($self,$sample_config,$config) = @_;

  my $all_list = $sample_config->{all} or die "you must indicate samples within the 'all' sublist";

  if (ref($all_list) eq "ARRAY") {
    $all_list = { map { $_ => 1 } @$all_list };
  }

  unless (ref($all_list) eq "HASH") {
    die "cannot handle a 'all' sample sublist that is not a hash or a id array";
  }

  foreach my $name (keys %$all_list) {
     my $value = $$all_list{$name};
     $value = 0 unless $value;
     my $info = $self->compose_sample_info($name,$value,$config);
     $$all_list{$name} = $info;
  }
  my $missing_list = { map { $_ => $$all_list{$_} } (grep { $$all_list{$_}->{missing} } keys %$all_list)  };
  $all_list = { map { $_ => $$all_list{$_} } (grep  { !exists($$missing_list{$_})} keys %$all_list ) };    

  $sample_config->{missing} = $missing_list;
  $sample_config->{all} = $all_list;

  return $sample_config;
}

sub compose_sample_info {
  my ($self,$name,$value,$config) = @_;
  if (ref($value) eq "") {
    unless ($value) {
      return $self->compose_sample_info($name,{},$config);
    }
    elsif (looks_like_number($value) && $value == 1) {
      return $self->compose_sample_info($name,{},$config);
    }
    else {
      return $self->compose_sample_info($name,{ id => $value },$config);
    }
  }
  elsif (ref($value) eq "ARRAY") {
    return $self->compose_sample_info($name,{ id => $name, lanes => $value });
  }
  elsif (ref($value) eq "HASH") {
    my $result;
    if ((exists($$value{id}) && exists($$value{lanes})) || $$value{missing}) {
      $result = $value;
    }
    else {
      my $id = $$value{id} || $name;
      my $lanes = $self->compose_sample_lanes($id,$$value{lanes},$config);
      $result = { id => $id, lanes => $lanes, scalar(keys %$lanes) == 0? (missing => 1):() };
    }
    $self->_check_sample_info_against_freeze($result) if defined $self->freeze;
    return $result;
  }
}

sub _check_sample_info_against_freeze {
  my ($self,$info) = @_;
  my $freeze = $self->freeze;
  my $id = $info->{id};
  my $sample_freeze = $freeze->find_sample($id);

  unless (defined $sample_freeze) {
    print STDERR "Sample '$id' missing in freeze!\n";
  }

  my %lane_and_files = $freeze->find_sample_lane_and_files($id);

  my $freeze_count = scalar(keys %lane_and_files);
  my $info_count = $info->{missing} ? 0 : scalar(keys(%{$info->{lanes}}));
  if ($info->{missing} && $freeze_count != 0) {
    delete $info->{missing};
    $info->{lanes} = {%lane_and_files};
    print STDERR "Missing in database but present in freeze '$id' adding them\n";
  }
  elsif ($info_count != $freeze_count) {
    print STDERR "Mismatching number of lanes for sample '$id' ($info_count vs $freeze_count)\n";
    $info->{missing} = 1;
    $info->{lanes} = {%lane_and_files};
  }
  else {
    my @lanes_diff = grep { ! exists ($lane_and_files{$_}) } keys (%{$info->{lanes}});
    if (scalar(@lanes_diff) > 0) {
      print STDERR "Different lane names between db and freeze for sample '$id'\n";
      $info->{lanes} = {%lane_and_files};
    }
    else {
      my @files_diff = grep { $lane_and_files{$_} ne $info->{lanes}->{$_} } keys (%lane_and_files);
      if ($#files_diff >= 0) {
        print STDERR "Different file names between db and freeze for sample '$id' and lanes '" . join(",", @files_diff);
        $info->{lanes} = {%lane_and_files};
      }
    }
  }
}

sub compose_sample_lanes {
  my ($self,$id,$value,$config) = @_;
  if (ref($value) eq "") {
    unless ($value) {
       return $self->sample_lane_and_files_from_id($id,$config);
    }
    else {
      return { $self->lane_name_and_file_from_scalar($value,$config) };
    }
  }
  elsif (ref($value) eq "ARRAY") {
    return { map { %{$self->compose_sample_lanes($_,$config)} } @$value };
  }
  elsif (ref($value) eq "HASH") {
    my %result = ();
    foreach my $key (keys %$value) {
      my $v = $$value{$key};
      next unless $v;
      $result{$key} = $v eq "1" ? $self->file_from_lane_name($key) : $v; 
    }
    return \%result;
  }
}

sub lane_name_and_file_from_scalar {
  my ($self,$value,$config) = @_;

  if ($value =~ /^\w+:\/\/.+$/i) {
    return ( $self->lane_name_from_url($value,$config) => $value );
  }
  else {
    my $file_name = $self->file_from_lane_name($value,$config);
    return $file_name ? ( $value => $file_name ) : ();
  }
}


sub sample_lane_and_files_from_id {
  my ($self,$id,$config) = @_;
  my $sm = $self->samtrak->manager_class_of('Sample');
  my $fm = $self->samtrak->manager_class_of('File');
  my $sample = ($id =~ /^\d+$/) ? $sm->get_by_id($id) : $sm->get_by_ox_code($id);
  $sample or die "could not find a sample with identifier '$id'";
  my $files = $fm->get_files(query => [ 'samples.id' => $sample->id ,
          'lanes.use_flag' => 1,
          name => { like => '%_nonhuman%' }], with_objects => [qw(samples lanes)], multi_many_ok => 1) || [];
  $files = [ $files ] if blessed($files);
  unless ($#$files >= 0) {
    warn "sample '"  .  join(",",$id,$sample->id) . "'  does not seem to have usable sequencing data available";
    return {};  
  }
  @$files = grep { $_->is_current } @$files;
  unless  ($#$files >= 0) {
    warn "sample '$id' does not seem to have current usable sequencing data available";
    return {}; 
  }
  my %lanes = (map { $self->lane_name_from_file_object($_) => $_->full_name } @$files); 
  return \%lanes;
}

sub lane_name_from_file_object {
  my ($self,$file) = @_;
  return $self->lane_name_from_file_name($file->name);
}

sub lane_name_from_file_name {
  my ($self,$fname) = @_;
  $fname =~ /(\d+)_(\d+)_nonhuman(#(\d+))?\.bam$/ or die "unexpected lane file name format '$fname'";
  my ($run,$lane,$tag) = ($1,$2,$4);
  if (defined($tag)) {
    return $run . "." . $lane . "#" . $tag;
  }
  else {
    return $run . "." . $lane;
  }
}

sub lane_name_from_url {
  my ($self,$url) = @_;
  if ($url =~ /^\w+:\/\/(.+)$/) {
    return $self->lane_name_from_file_name($1);
  }
  else {
    return $self->lane_name_from_file_name($url);
  }
}

sub file_from_lane_name {
  my ($self,$lane_name) = @_;
  my $fm = $self->samtrak->manager_class_of('File');
  $lane_name =~ /^(\d+)\.(\d+)(#(\d+))$?/;
  my ($run,$lane,$tag) = ($1,$2,$3);
  my $file_name = join("",$run,'.',$lane,'_nonhuman',(defined $tag) ? ('#' . $tag):(),'.bam');
  my $files = $fm->get_files(query => [
           filesystem => { like => 'irods' },
           type => { like => 'bam' },
           is_current => 1,
           name => $file_name ],
         with_objects => [qw(samples lanes lanes.lib samples.tag_values)], multi_many_ok => 1) || [];
  $files = [ $files ] unless ref($files) eq "ARRAY";
  unless ($#$files >= 0) {
     warn "could not find file '$file_name' for a lane '$lane_name'";
     return undef;
  }
  elsif ($#$files > 0) {
     warn "more than one file '$file_name' for lane '$lane_name' returning first available";
     return $$files[0]->full_name;
  }
  else {
      return $$files[0]->full_name;
  }
}



sub complete_config_sequences {
  my ($self,$config) = @_;
   
  if (exists $config->{sequences}) {
    $config->{sequences} = $self->process_config_sequences($config->{sequences},$config);
  }
  else {
    $config->{sequences} = $self->process_config_sequences($self->default_config_sequences($config),$config);
  }
  return $config;
}

sub process_config_intervals {
  my ($self,$intervals,$config) = @_;
  unless (defined $intervals) {
    return $self->default_config_intervals($config);
  }
  elsif (ref($intervals) eq "") {
    my $result;
      my $refCommand = $self->reference_command($config);
    if ($intervals =~ /^(\d+)n?$/) {
      $result = $refCommand->calculate_intervals($1,1000000000);
    }
    elsif ($intervals =~ /^(\d+)s?$/) {
      $result = $refCommand->calculate_intervals(1,$1);
    }
    else {
      die "cannot process interval spec string $intervals"
    }
    return { map { $$_[3] => "$$_[0]:$$_[1]-$$_[2]" } @$result };
  }
  elsif (ref($intervals) eq "ARRAY") {
    my $index = 0;
    return { map { ("interval_" . $index++) => $$_ } @$intervals };
  }
  elsif (ref($intervals) eq "HASH") {
    if (exists($intervals->{number}) || exists($intervals->{size})) {
       my $refCommand = $self->reference_command($config);
       my $number = $intervals->{number} || 1;
       my $size = $intervals->{size} || 100000000000;
       return { map { $$_[3] => "$$_[0]:$$_[1]-$$_[2]" } @{$refCommand->calculate_intervals($number,$size)} };
    }
    else {
       return $intervals;
    }
  }
  else {
    die "cannot handle interval specification with ref-type " . ref($intervals);
  }
}

sub process_config_sequences {
  my ($self,$sequences,$config) = @_;

  my $result;
  my $ref = $self->reference($config);
  unless (defined $sequences) {
    $result = $self->default_config_sequences($config);
  }
  elsif (ref($sequences) eq "") {
    my $ref = $self->reference();
    if ($sequences eq "all") {
       $result = $ref->chromosomes();
    }
    elsif ($sequences eq "nuclear") {
       $result = $ref->chromosomes();
       $result = [ grep { $_ !~ /^MT/ && $_ !~ /PLAST/ && $_ !~ /APICO/ && $_ !~ /MITO/ } @$result ];
    }
  }
  elsif (ref($sequences) eq "ARRAY") {
    $result = $sequences;
  }
  return $result;
}


sub default_config_sequences {
  my ($self,$config) = @_;
  my $ref = $self->reference($config);
  return $ref->chromosomes;
}

sub default_config_intervals {
  my ($self,$config) = @_;
  my $ref = $self->reference_command($config);
  my $intervals = $ref->default_intervals();
  my %result = (map { $$_[3] => "$$_[0]:$$_[1]-$$_[2]" } @$intervals);
  return \%result;
}


sub write_into_file {
  my ($self,$filename,$content_ref) = @_;
  my $fh = IO::File->new($filename,'w')
    or return "could not open '$filename' to write";
  print $fh $$content_ref;
  $fh->close();
}


1;
