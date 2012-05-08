package Maltools::Commands::Query;

use Moose;
extends 'Maltools::Command';

use strict;
use warnings;

use Getopt::Long qw(:config require_order);

use IO::File;
use DateTime;
use DateTime::Format::ISO8601;
use File::Spec::Functions qw(catfile);
use base 'Maltools::Command';
use Digest::MD5;
use Scalar::Util qw(blessed);
use Maltools::Config qw();

has 'data_config' => (is => 'ro', lazy => 1 , default => sub { Maltools::Config::data_config() });

sub help_summary {
   return 'execute some query relevant to the AGV pipeline';
}

sub help_text {
   my $self = shift;
   my $prog = $self->cl_name;
   my $summary = $self->help_summary;
   return <<EOM
Command:
  
  query - $summary

Summary:

  This command can execute some common queries against the sample tracking system and 
  other support files in order to obtain information that might be used for some
  other components of the pipeline. 

  It does have a set of predefined queries that might be extended by developing the code.
  Unfortunatelly there is no mechanism yet to perform arbitrary queries.

Syntax:

  $prog query query-name [-o output-file] [query arguments]

Supported Queries:

  usable-sequence-files [-s sample] [-t taxon] 

    returns the table of available files containing sequencing that has been marked as usable
    i.e. has pass some for of QC. The query can be constrain to a certain sample or taxon. 

  incoming-lanes

    lists usable lanes
  
  group-samples 
  
    lists samples that belong to certain group

  mapped-samples 

    ???

  genotyping-group

    groups of samples to be genotyped together

  good-samples

    list of samples that we think are good based on genotype data.

  samtrak-env [shell-name]

    returns the environment setup required to connect to the samtrak database.

    It will output a file that should be sourced to set up the environment
    appropriate to use Samtrak API to access SolarisDB tracking information. 

    Due to the fact that different shell program set-up variables in a different 
    way, you should provide what kind of shell you plan tu use it with. Supported
    shells are (ba)sh or (t)csh. If none is provided this query will try to guess
    it, but this is not reliable and in fact cannot be bullet-proof so please
    provide one.

EOM
}

sub execute {
  my ($self,$samtrak_site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output_file;
  GetOptions("o=s" => \$output_file);
  $params->{arguments} = \@ARGV;
  $params->{ofh} = (defined $output_file) ? IO::File->new($output_file,"w") : \*STDOUT;
  my $query_name = shift @{$params->{arguments}}
     or return $self->error_return("you need to specify a query name");
  if ($query_name eq "incoming-lanes") {
     return $self->incoming_lanes($samtrak_site,$params);
  }
  elsif ($query_name eq "samtrak-env") {
     return $self->samtrak_env($samtrak_site,$params);
  }
  elsif ($query_name eq "sequenced-samples") {
     return $self->sequenced_samples($samtrak_site,$params);
  }
  elsif ($query_name eq "group-samples") {
     return $self->group_samples($samtrak_site,$params);
  }
  elsif ($query_name eq "incoming-files") {
     return $self->incoming_files($samtrak_site,$params);
  }
  elsif ($query_name eq "mapped-samples") {
     return $self->mapped_samples($samtrak_site,$params);
  }
  elsif ($query_name eq "genotyping-groups") {
     return $self->genotyping_groups($samtrak_site,$params);
  }
  elsif ($query_name eq "good-samples") {
     return $self->good_samples($samtrak_site,$params);
  }
  elsif ($query_name eq "usable-sequence-files") {
     return $self->usable_sequence_files($samtrak_site,$params);
  }
  else {
     return $self->error_return("unknown query name $query_name");
  }
}

sub samtrak_env {
  my ($self,$samtrak_site,$params) = @_;
  my $config = Maltools::Config->new();
  my $db = $config->get('samtrak/db');
  my $shell = shift (@{$params->{arguments}}) || $ENV{SHELL} || "sh";
  my $is_csh = $shell =~ /csh$/;
  if (ref($db) eq "HASH") {
    my %db = %$db;
    if ($is_csh) {
        print "setenv SAMTRAK_DBTYPE default\n";
        print "setenv SAMTRAK_DBNAME $db{schema}\n";
        print "setenv SAMTRAK_DBUSER $db{username}\n";
        print "setenv SAMTRAK_DBPASS $db{password}\n" if $db{password};
        print "setenv SAMTRAK_DBHOST $db{host}\n";
        print "setenv SAMTRAK_DBPORT $db{port}\n";
    }
    else {
        print "export SAMTRAK_DBTYPE=default\n";
        print "export SAMTRAK_DBNAME=$db{schema}\n";
        print "export SAMTRAK_DBUSER=$db{username}\n";
        print "export SAMTRAK_DBPASS=$db{password}\n" if $db{password};
        print "export SAMTRAK_DBHOST=$db{host}\n";
        print "export SAMTRAK_DBPORT=$db{port}\n";
    }
  }
}

sub sequenced_samples {
  my ($self,$samtrak,$params) = @_;
  my $sm = $samtrak->manager_class_of('Sample');
  my $samples = $sm->get_samples(query => [ 'dna_sample.sample2lanes.use_flag' => 1 ],
           with_objects => [qw(tag_values tag_values.tag dna_sample dna_sample.sample2lanes taxon biosample_location)], multi_many_ok => 1 );
  print join("\t",qw(id ox_code taxon biosample_location)),"\n";
  foreach my $s (@$samples) {
     print join("\t",$s->id,$s->code,$s->taxon->scientific_name, !defined($s->biosample_location) ? 'NULL': $s->biosample_location->name),"\n";
  }
  return $self->ok_return;
}

sub usable_sequence_files {
  my ($self,$samtrak,$params) = @_;
  my $args = $params->{arguments};
  local @ARGV = @$args;
  my @sample = (); 
  my @taxon = ();
  GetOptions("sample|s=s" => \@sample, "taxon|t=s" => \@taxon);
  my @constraints = (use_flag => 1);
  if ($#sample >= 0) {
     push @constraints, ($self->sample_constraints($samtrak,@sample));
  }
  if ($#taxon >= 0) {
     push @constraints, ($self->taxon_constraints($samtrak,@taxon));
  }
  my $slm = $samtrak->manager_class_of('Sample2Lane');
  my $lanes = $slm->get_sample2lanes(query => \@constraints, with_objects =>[qw(dna_sample dna_sample.sample dna_sample.sample.taxon)]);
  my $ofh = $params->{ofh} || \*STDOUT; 
  print $ofh join("\t",qw(lane_id tag_id sample_id taxon_id run_num lane_num forward_file reverse_file
                            sample_name taxon_name)),"\n";
  return $self->file_print_from_lanes($ofh,@$lanes); 
}

sub sample_constraints {
  my $self = shift;
  my $samtrak = shift;
  my @sample = @_;
  my @code_likes = ();
  my @codes = ();
  my @ids = ();
  foreach my $s (@sample) {
    my $escaped = 0;
    $s =~ s/\%/\\\%/g and $escaped = 1;
    $s =~ s/\_/\\\_/g and $escaped = 1;
    $s =~ s/\*/\%/g and $escaped = 1;
    $s =~ s/\?/\_/g and $escaped = 1;
    if ($escaped) {
       push @code_likes, $s;
    }
    elsif ($s =~ /^\s*(\d+)\s*$/) {
       push @ids, $s;
    }
    else {
       push @codes, $s;
    }
  }
  if ($#codes >= 0 || $#code_likes >= 0) {
    my $tvm = $samtrak->manager_class_of('TagValue');
    
    my @constraints = ();
    push @constraints,   'value' => \@codes  if $#codes >= 0;
    push @constraints,  map { ('value' => { like => $_ })  } @code_likes;
    
    my $sample_ids = $tvm->get_objects(
       query => ['tag.abbreviation' => 'ox_code', 
            or =>  \@constraints],
           select => 'sample_tag_values.sample_id', 
       group_by => 'sample_tag_values.sample_id', require_objects => ['sample_tag_values'])  || [];
    push @ids, map { $_->id } @$sample_ids;
  }
  my %ids = map { $_ => 1 } @ids;
  @ids = keys (%ids);
  return ('dna_sample.sample.id' => \@ids);
}

sub taxon_constraints {
  my $self = shift;
  my $samtrak = shift;
  my @taxon = @_;
  my @ids = ();
  my @names = ();
  my @name_likes = ();
  foreach my $t (@taxon) {
    my $escaped = 0;
    $t =~ s/\%/\\\%/g and $escaped = 1;
    $t =~ s/\_/\\\_/g and $escaped = 1;
    $t =~ s/\*/\%/g and $escaped = 1;
    $t =~ s/\?/\_/g and $escaped = 1;
    if ($escaped) {
       push @name_likes, $t;
    }
    elsif ($t =~ /^\s*(\d+)\s*$/) {
       push @ids, $t;
    }
    else {
       push @names, $t;
    }
  }
  if ($#names >= 0 || $#name_likes >= 0) {
    my $tm = $samtrak->manager_class_of('Taxon');
    
    my @constraints = ();
    push @constraints, 'ncbi_taxid' => \@names if $#names >= 0;
    push @constraints, map { ('scientific_name' => { like => $_ })  } @name_likes;
    push @constraints, map { ('nickname' => { like => $_ }) } @name_likes;
    
    my $taxon_ids = $tm->get_objects(query => [ or =>  \@constraints], select => 'id') || [];
    push @ids, map { $_->id } @$taxon_ids;
  }
  my %ids = map { $_ => 1 } @ids;
  @ids = keys (%ids);
  return ('dna_sample.sample.taxon_id' => \@ids);
}

sub file_print_from_lanes {
  my $self = shift;
  my $ofh = shift;
  my @lanes = @_;
  foreach my $l (@lanes) {
    my $f_file = $l->f_file;
    $f_file .= "_nonhuman" unless $f_file =~ /nonhuman/ || $l->run < 5000;
    my $r_file = $l->r_file;
    $r_file .= "_nonhuman" unless $r_file =~ /nonhuman/ || $l->run < 5000;
    my $lane_file_name =  "/seq/" .$l->run . "/" . $l->run . "_" . $l->lane . "_nonhuman.bam";
    print $ofh join("\t",$l->lane_id,$l->tag_id || "NA",$l->sample_id,
          $l->dna_sample->sample->taxon->id,$l->run,$l->lane,
          $lane_file_name,$lane_file_name,$l->ox_code,$l->species),"\n";
  } 
  return $self->ok_return;
}


sub incoming_lanes {
  my ($self,$samtrak,$params) = @_;
  my $slm = $samtrak->manager_class_of('Sample2Lane');
  my $lanes = $slm->get_sample2lanes(query => [ 'use_flag' => 1 ], 
           with_object => [qw(dna_sample dna_sample.sample dna_sample.sample.taxon)]);
  my $ofh = $params->{ofh};
  print $ofh join("\t",qw(lane_id tag_id sample_id taxon_id run_num lane_num forward_file reverse_file
                            sample_name taxon_name)),"\n";
  foreach my $l (@$lanes) {
    my $f_file = $l->f_file;
    $f_file .= "_nonhuman" unless $f_file =~ /nonhuman/ || $l->run < 5000;
    my $r_file = $l->r_file;
    $r_file .= "_nonhuman" unless $r_file =~ /nonhuman/ || $l->run < 5000;
    my $lane_file_name =  "/seq/" .$l->run . "/" . $l->run . "_" . $l->lane . "_nonhuman.bam";
    
    print $ofh join("\t",$l->lane_id,$l->tag_id || "NA",$l->sample_id,
          $l->dna_sample->sample->taxon->id,$l->run,$l->lane,
          $lane_file_name,$lane_file_name,$l->ox_code,$l->species),"\n";
  } 
  return $self->ok_return;
}

sub group_samples {
  my ($self,$samtrak,$params) = @_;
  local @ARGV = @{$params->{arguments}};
  my $query_sample = undef;
  my $no_header = 0;
  my $time = undef;
  GetOptions('time=s' => \$time, 'sample=s' => \$query_sample, 'H' => \$no_header);
  my @target_sample_sets = $self->_resolve_sample_sets($samtrak,@ARGV);
  $#target_sample_sets >= 0 or return $self->error_return("there is no sample set that matches query criteria");
  my $fm = $samtrak->manager_class_of('File');
  my @sample_ids = map { map { $_->sample_id } @{$_->items} } @target_sample_sets;
  my %samples_outputs = ();
  my $files = $fm->get_files(query => ['sample.id' => \@sample_ids,
          filesystem => { like => 'irods' },
          type => { like => 'bam' },
          #is_current => 1, 
          creation_date => { le => (defined($time) ? DateTime::Format::ISO8601->parse_datetime($time) : DateTime->now()) },
          'lanes.use_flag' => 1,
          name => { like => '%_nonhuman%' }], with_objects => [qw(samples lanes lanes.lib samples.tag_values)], multi_many_ok => 1);
  print "ox_code\n" unless ($no_header);
  foreach my $file (@$files) {
    my $sample = $file->samples; $sample = $sample->[0] if ref($sample) eq "ARRAY";
    my $sample_code = $sample->attribute('ox_code')->get_externalized();
    print $sample_code,"\n" unless exists $samples_outputs{$sample_code};
    $samples_outputs{$sample_code}++;
  }
  return $self->ok_return();
}

sub _refine_sample_ids {
  my ($self,$samtrak,$sample_ids,$query_sample) = @_;
  return unless $query_sample;

  my $sm = $samtrak->manager_class_of('Sample');
  my $candidate_samples = $sm->get_samples(query => [ 'tag_values.value' => { like => $query_sample } ],
       require_objects => ['tag_values']) || [];
  $candidate_samples = [ $candidate_samples ] if blessed($candidate_samples);
  if ($#$sample_ids < 0) {
    push @$sample_ids, map { $_->id } @$candidate_samples;
  }
}

sub _resolve_sample_sets {
  my ($self,$samtrak,@ds_names) = @_;
  my @target_sample_sets = ();
  my $ssm = $samtrak->manager_class_of('SampleSet');
  if ($#ds_names < 0) {
    @target_sample_sets = ();
    #@target_sample_sets = @{$ssm->get_sample_sets(query => [ name => { like => 'AGV%' } ], with_objects => ['items'])};
  }
  else {
    foreach my $ss_name (@ds_names) {
       my $ss = undef;
       my @query_args = ();
       if ($ss_name =~ /^[A-Z]+$/) {
          push @query_args, query => [name => {like => 'AGV-%-' . $ss_name . '-%'}];
          push @query_args, sort_by => 'name desc';
       }
       elsif ($ss_name =~ /^AGV-\d{4}-[A-Z]+-\d{4}$/) {
          push @query_args, query => [name => $ss_name];
       }
       else {
          $ss_name =~ s/\*/%/g;
          $ss_name =~ s/\?/_/g;
          push @query_args, query => [name => {like => $ss_name}];
          push @query_args, sort_by => 'name desc';
       }
       my $groups = $ssm->get_sample_sets(@query_args,with_objects => [qw(items)]) || [];
       $groups = [ $groups ] unless ref($groups) eq "ARRAY";
       push @target_sample_sets, @$groups;
    }
  }
  return wantarray ? @target_sample_sets : \@target_sample_sets;
}

sub incoming_files {
  my ($self,$samtrak,$params) = @_;
  require DateTime::Format::ISO8601;
  local @ARGV = @{$params->{arguments}};
  my $run = undef;
  my $lane = undef;
  my $query_sample = undef;
  my $no_header = 0;
  my $files_only = 0;
  my $samples_only = 0;
  my $time = undef;
  GetOptions('time|t=s' => \$time, 'run=i' => \$run,'lane=i' => \$lane,
             'sample=s' => \$query_sample, 'H' => \$no_header, 'f' => \$files_only,'s' => \$samples_only);
  
  $time = defined $time ? DateTime::Format::ISO8601->parse_datetime($time) : DateTime->now();
  $self->error_return("you cannot request files and sample only output at the sam time") if ($files_only && $samples_only);
  my @ds_names = @ARGV;
  my @target_sample_sets = $self->_resolve_sample_sets($samtrak,@ARGV);
  #$#target_sample_sets >= 0 or return $self->error_return("there is no sample set that matches query criteria");
  my $fm = $samtrak->manager_class_of('File');
  my @sample_ids = map { map { $_->sample_id } @{$_->items} } @target_sample_sets;
  $self->_refine_sample_ids($samtrak,\@sample_ids,$query_sample);
  my @constraints = ();
  push @constraints, 'samples.id' => \@sample_ids if $#sample_ids >= 0;
  push @constraints, 'lanes.f_run_id' => $run if defined $run;
  push @constraints, 'lanes.lane' => $lane if defined $lane;
  my %samples_outputs = ();
  my $files = $fm->get_files(query => [@constraints, # 'samples.id' => \@sample_ids, 
          filesystem => { like => 'irods' }, 
          type => { like => 'bam' }, 
          #is_current => 1, 
          'lanes.use_flag' => 1,
          creation_date => { le => $time }, 
          name => { like => '%_nonhuman%' }], with_objects => [qw(samples lanes lanes.lib samples.tag_values)], multi_many_ok => 1);
  unless ($no_header) {
    print join("\t",qw(file_id sample_id lane_id lib_id file_name ox_code run_name lane_index tag_num lib_name prep_method creation_date md5)),"\n" unless $files_only || $samples_only;
    print "file_name\n" if $files_only;
    print "ox_code\n" if $samples_only;
  }
  foreach my $file (@$files) {
    my $lane = $file->lanes; $lane = $lane->[0] if ref($lane) eq "ARRAY";
    my $lib = $lane->lib;
    my $tag_num = $file->name =~ /\#(\d+)/ ? $1 : '';
    my $sample = $file->samples; $sample = $sample->[0] if ref($sample) eq "ARRAY";
    my $sample_code = $sample->attribute('ox_code')->get_externalized();
    next if ($query_sample && ($sample_code ne $query_sample && $query_sample ne $sample->id));  
    print join("\t", map { defined $_ ?  $_ : "" } ($file->id,$sample->id,$lane->id,$lib->id,$file->full_name,
           $sample_code,$lane->f_run_id,$lane->lane,$tag_num,$lib->lib_name,$lib->prep_method,$file->creation_date,$file->md5)),"\n" unless $files_only || $samples_only;
    print $file->full_name,"\n" if $files_only; 
    print $sample_code,"\n" if $samples_only && !exists $samples_outputs{$sample_code};
    $samples_outputs{$sample_code}++;
  }
  return $self->ok_return();
}

sub mapped_samples {
  my ($self,$samtrak,$params) = @_;
  my $qcfsm = $samtrak->manager_class_of('QCFlagStat');
  my $qcfs = $qcfsm->get_qcflagstats(query => [ 'sample.taxon.scientific_name' => {like => 'anoph%'} , current => 1],
        require_objects => [qw(sample sample.taxon sample.biosample_location)],
        sort_by => [qw(sample.taxon_id sample.biosample_location_id)]);
  my $ofh = $params->{ofh};
  print $ofh join("\t",qw(id ox_code taxon_id scientific_name biosample_location_id biosample_location_name)),"\n";
  foreach my $qc (@$qcfs) {
    my $t = $qc->sample->taxon;
    my $l = $qc->sample->biosample_location;
    print $ofh join("\t",$qc->sample_id,$qc->ox_code,$t->id,$t->scientific_name,$l->id,$l->name),"\n";
  }
  return $self->ok_return; 
}

sub genotyping_groups {
  my ($self,$samtrak,$params) = @_;
  my $qcfsm = $samtrak->manager_class_of('QCFlagStat');
  my $qcfs = $qcfsm->get_qcflagstats(query => [ 'sample.taxon.scientific_name' => {like => 'anoph%'} , current => 1],
        require_objects => [qw(sample sample.taxon sample.biosample_location sample.biosample_location.country)],
        sort_by => [qw(sample.taxon_id sample.biosample_location_id)]);
  my $ofh = $params->{ofh};
  print $ofh join("\t",qw(code taxon_scientific_name total per-country lmd digest samples)),"\n";
  my $current_spc = _gg_spc_from_qc($qcfs->[0]);
  my $current_loc = _gg_loc_from_qc($qcfs->[0]); 
  my @qcs = ();
  my %locs = ();
  for (my $i = 0; $i < scalar(@$qcfs); $i++) {
    my $qc = $qcfs->[$i];
    my $spc = _gg_spc_from_qc($qc);
    my $loc = _gg_loc_from_qc($qc);
    $locs{$loc}++;
    #if ($i == (scalar(@$qcfs) - 1) || $spc ne $current_spc || $loc ne $current_loc) {
    if ($i == (scalar(@$qcfs) - 1) || $spc ne $current_spc) {
      my ($digest,$lmd) = _gg_digest_and_lmd(@qcs);
      my @sname_bits = split(/\s+/,$qcfs->[$i-1]->sample->taxon->scientific_name);
    #  print $ofh join("\t",$current_spc . "-" . $current_loc, join(" ",@sname_bits[0 .. 1]), 
    #                $qcfs->[$i-1]->sample->biosample_location->name, scalar(@qcs), $lmd, $digest, join ("; ",map {$_->ox_code} @qcs )),"\n";
      print $ofh join("\t",$current_spc, join(" ",@sname_bits[0 .. 1]), scalar(@qcs),
                    join("; ",map { $_ . ":" . $locs{$_} } sort keys %locs), $lmd, $digest, join ("; ",map {$_->ox_code} @qcs )),"\n";
      @qcs = ();
      $current_spc = $spc;
      $current_loc = $loc;
    }
    push @qcs, $qc;
  }
  return $self->ok_return; 
}

sub good_samples {
    my ($self,$samtrak,$params) = @_;
    local @ARGV;
    @ARGV = @{$params->{arguments}};
    my $group = shift @ARGV or return $self->error_return("you need to speficy a genotyping group");
    my $dc = $self->data_config;
    my $god = $dc->genotyping_out_dir(group => $group);
    -d $god or $self->error_return("you need to calculate genotypes fro group '$group' before performing this query. Please use 'agv pipeline genotyping $group' for that");
    my $sample_tab_file = catfile($god,'sample-table.tab');
    unless (-f $sample_tab_file) {
       print STDERR "Calculating good samples list...\n";
       my $engine = $self->resolve_engine;
       my $tool = Maltools::Tool->tool_by_name('vcf-goodSamples') or return $self->error_return("good sample calculation error","could not find a tool with id 'vcf-goodSamples'");
       my $job = $tool->job(inputs => { in => catfile($god,'genotypes.vcf')}, outputs => { out => $sample_tab_file} );
       unless ($engine->run_job($job)) {
          return $self->error_return("error executing vcf-goodSamples tool with return code " . $job->return_code . " and message: " . ($job->error_message || 'no-message'));
       }
       -f $sample_tab_file or return $self->error_return("did not succeed in running the vcf-goodSamples tool");
    }
    my $sample_tab_fh = IO::File->new($sample_tab_file,'r') or return $self->error_return("could not genotyping sample table file '$sample_tab_file'");
    while (my $line = <$sample_tab_fh>) {
       next if $line =~ /^#/;
       my ($name,$class) = split(/\t/,$line);
       next unless uc($class) eq "GOOD";
       print $name,"\n";
    }
    $sample_tab_fh->close();
    return $self->ok_return;
}

sub _gg_digest_and_lmd {
  my $latest = '1970-01-01';
  my $md5 = Digest::MD5->new;
  my @qcs = sort { $a->id <=> $b->id } @_;
  foreach my $qc (@qcs) {
    $md5->add($qc->id);
    $md5->add($qc->ox_code);
    my $td = $qc->added;
    $md5->add($td);
    my $tds = "" . $td;
    $latest = $tds if ($tds gt $latest); 
  }
  return ($md5->b64digest,$latest);
}

sub _gg_spc_from_qc {
  my $qc = shift;
  my $name = $qc->sample->taxon->scientific_name;
  my @name_bits = split(/\s+/,uc($name));
  my $result = substr($name_bits[0],0,1) . substr($name_bits[1],0,1);
  return $result;
}

sub _gg_loc_from_qc {
  my $qc = shift;
  my $iso2 = $qc->sample->biosample_location->country->iso2;
  return $iso2;
}

1;
