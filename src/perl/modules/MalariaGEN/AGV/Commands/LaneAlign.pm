package MalariaGEN::AGV::Commands::LaneAlign;

use base 'MalariaGEN::AGV::Command';

use Cwd qw(realpath);
use strict;
use warnings;

use MalariaGEN::AGV::Config qw(sanger_config reference_config jobs_config);
use MalariaGEN::AGV::Tool;
use Getopt::Long;
use IO::File;

sub help_summary {
   return 'aligns reads in a sam/bam input file with the reference';
}

sub help_text {
   return "lane-align\n\t" .
          $_[0]->help_summary . "\n" .
          "Syntaxis:\n" .
          $_[0]->cl_name . " lane-align single-lane-query\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $lane_id = undef;
  my $lane_idx = undef;
  my $lane_name = undef;
  my $run = undef;
  my $file_id = undef;
  my $output = undef;
  my $tag_num = undef;
  GetOptions("output|out|o=s" => \$output, "tag-num|tn=i" => \$tag_num, "lane-id|li=i" => \$lane_id, 
            "lane-idx|lx=i" => \$lane_idx,"run-num|run=i" => \$run, "lane-name|ln=s" => \$lane_name, "file-id|fi=i" => \$file_id);
  if (defined $lane_name && $lane_name =~ s/\#(\d+)//) {
    $tag_num = $1;
  }
  $output or return $self->error_return("you must indicate and output file name");
  my @query_where = (filesystem => { like => 'irods' }, type => { like => 'bam'},'lanes.use_flag' => 1, name => { like => '%_nonhuman%'});
  push @query_where ,id => $file_id if defined $file_id;
  push @query_where , name => { like => '%#' . $tag_num . '%'} if defined $tag_num;
  push @query_where ,'lanes.f_run_id' => $run if defined $run;
  push @query_where ,'lanes.id' => $lane_id if defined $lane_id;
  push @query_where ,'lanes.lane' => $lane_idx if defined $lane_idx;
  push @query_where , \"concat(t3.f_run_id,'.',t3.lane) = '$lane_name'" if defined $lane_name;
  my $fm = $site->manager_class_of('File');
  my $files = $fm->get_files(query => \@query_where, with_objects => [qw(lanes samples lanes.lib samples.tag_values)], multi_many_ok => 1, sort_by => 'creation_date desc');
  unless (defined $files || scalar(@$files) > 0) {
    return $self->error_return("there is no such a lane/file that complies with the selecting conditions");
  }
  else {
    my $file = $files->[0];
    print join("\t",qw(file_id sample_id lane_id lib_id file_name ox_code run_name lane_index lib_name prep_method creation_date md5)),"\n";
    my $lane = $file->lanes; $lane = $lane->[0] if ref($lane) eq "ARRAY";
    my $lib = $lane->lib;
    my $sample = $file->samples; $sample = $sample->[0] if ref($sample) eq "ARRAY";
    my $sample_code = $sample->attribute('ox_code')->get_externalized();
    print STDERR join("\t", map { defined $_ ?  $_ : "" } ($file->id,$sample->id,$lane->id,$lib->id,$file->full_name,
           $sample_code,$lane->f_run_id,$lane->lane,$lib->lib_name,$lib->prep_method,$file->creation_date,$file->md5)),"\n";
    my $tool = MalariaGEN::AGV::Tool->tool_by_name('pipeline-laneAlign');
    my $job = $tool->job(inputs => { ( $file->full_name =~ /^\/seq/ ? "in_irods" : "in")  => $file->full_name, ref => reference_config()->file_name(), rg_id => join(".",$lane->f_run_id,$lane->lane . (defined $tag_num ? '#' . $tag_num : '' ),$lane->id), 
                                rg_sm => $sample_code, rg_lb => $lib->lib_name, rg_cn => "Wellcome Trust Sanger Institute", 
                            rg_ds => "Plasmodium Genome Variation Project, sample " . $sample_code  . " prep. method " . $lib->prep_method,
                            pg_cl => "agv lane-align " . join(" ",@{$params->{arguments}}),
                            pg_id => "agv_lane_align",
                            pg_vn => $MalariaGEN::AGV::Config::VERSION,
                            pg_pn => "PGV Pipeline Lane aligner",
                            pg_ds => "Obtain sequencing Raw BAM file in iRODS '" . $file->full_name . "' and realign it based in our parameters" }, 
                         outputs => { out => $output  }); 
    my $engine = $self->resolve_engine;
    $engine->run_job($job) or return $self->error_return("error executing the pipeline-laneAlign tool with exit code " . $job->return_code . " and message: " . $job->error_message);
    return $self->ok_return;
  }
}


1;

__DATA__
