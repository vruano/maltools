use strict;
use warnings;

use JSON::XS;
use Samtrak::DB;

use Samtrak::Data::Sample::Manager;
use Samtrak::Data::Lane::Manager;
use Samtrak::Data::Tag::Manager;
use Samtrak::Data::TagValue::Manager;

my $sm = "Samtrak::Data::Sample::Manager";
my $lm = "Samtrak::Data::Lane::Manager";
my $tvm = "Samtrak::Data::TagValue::Manager";
my $tm = "Samtrak::Data::Tag::Manager";

my $json = JSON::XS->new();


my $ox_code_tag = $tm->get_by_abbreviation('ox_code');


while (my $line = <STDIN>) {
   $json->incr_parse($line);
}



my $mani = $json->incr_parse();


my %samples = ();

my $missing = $mani->{samples}->{missing};
my $all = $mani->{samples}->{all};
foreach my $s (keys %$missing) {
   my $ss = $sm->get_by_ox_code($s,with_objects => ['tag_values','tag_values.tag']);
   $samples{$s} = [$s,$ss->id,join(",",@{$ss->attribute('ox_src_code')->get()}),other_similar($ss),0];
}

foreach my $s (keys %$all) {
   my $lanes = $all->{$s}->{lanes};
   my $ss = $sm->get_by_ox_code($s,with_objects => ['tag_values','tag_values.tag']);
   my $ll = { map { lane_id_to_idobj($_) } keys (%$lanes) };
   $samples{$s} = [$s,$ss->id, join(",",@{$ss->attribute('ox_src_code')->get()}),other_similar($ss),int(scalar(keys %$lanes)),map { lane_to_output($ll->{$_}) } sort keys (%$ll)];


}

foreach my $s (sort keys %samples) {
  print STDOUT join("\t",@{$samples{$s}}),"\n";
}

sub lane_id_to_idobj {
  my $id = shift;
  $id =~ /^(\d+)\.(\d+)$/ or die "cannot process lane name '$id'";
  my ($run,$idx) = ($1,$2);
  my $lanes = $lm->get_lanes(query => [ f_run_id => $run, lane => $idx], with_objects => ['lib','lib.seqsamples']);
  die "could not find lane '$id' id Solaris" if ($#$lanes < 0);
  die "too many entires for lane '$id' in Solaris" if ($#$lanes > 0);
  return ($id , $lanes->[0]);
}

sub other_similar {
  my $ss = shift;
  my $id = $ss->attribute('ox_code')->get();
  my $pattern = $id;
  $pattern =~ s/\-(\S+)$//;
  $pattern = "\%$pattern\%";
  my $tvs = $tvm->get_tag_values(query => [ tag_id => $ox_code_tag->id, value => { like => $pattern } ]);
  my $tvs_codes = [ grep { $id ne $_ } map { $_->value } @$tvs];
  grep { $id eq $_ } map { $_->value } @$tvs or die "no self-found";
  return join(",",@$tvs_codes);
}

sub lane_to_output {
  my $lane = shift;
  return $lane->f_run_id . "." . $lane->lane . " (" . $lane->lib->lib_name . ")";
}
