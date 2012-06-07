package Maltools::Commands::AnnotateSampleTree;

use strict;
use warnings;
use Moose;

use Maltools::Tool;
use Getopt::Long qw(:config no_ignore_case);
use File::Copy qw(copy);
use Samtrak::Site;
use Bio::TreeIO;
use URI::Escape;

extends 'Maltools::Command';



sub hidden {
  return 1;
}

has '+engine_name' => ( default => 'local' );

sub help_summary {
   return 'reformat a sample tree adding attributes values to it';
}

sub help_text {
   return "sample-attr-tree\n\t" .
          $_[0]->help_summary . "\n" .
          "Synopsis:\n" .
          $_[0]->cl_name . " sample-attr-tree [-i input-tree | < input-tree] [-o output-tree | > output-tree ] [-I input-format] [-O output-format] [-d (property|suffix) ] ATTR_ID1 ATTR_ID2 ... \n\n" .
          "Options:\n" .
          "\t-I/-O\tvalid formats are 'newick', 'nhx', 'nexus' and 'phyloxml' for either and 'tabtree' for output only\n".
          "\t-d\tindicates how to add the attributes to the node:\n" .
          "\t\t'property' will attempt to use the output format natural way to add tags, attributes or properties to the node\n"  .
          "\t\t'suffix' will add the value into the leave id, if only one attribute is requested it just adds the value, if there are several it will add ATTR_ID1:VALUE_ATTR_I2:VALUE\n";
}

sub execute {
  my ($self,$site,$params) = @_;
  local @ARGV;
  @ARGV = @{$params->{arguments}};
  my $output = "-";
  my $input = "-";
  my $decoration_type = "suffix";
  my $output_format = "newick";
  my $input_format = "newick";
  my @attributes = ();
  GetOptions("i|in|input=s" => \$input, "out|output|o=s" => \$output, "I|input-format|in-format=s" => \$input_format, "O|output-format|out-format=s" => \$output_format,
             "d|decoration-type|d-type=s" => \$decoration_type);
  @attributes = @ARGV;
  $decoration_type eq "suffix" || $decoration_type eq "property" or return $self->error_return("wrong decoration type '$decoration_type' must be either 'suffix' or 'property'");
  grep { $_ eq $input_format } (qw(newick nhx phyloxml nexus)) or return $self->error_return("not supported input format '$input_format'"); 
  grep { $_ eq $output_format } (qw(newick nhx phyloxml nexus tabtree)) or return $self->error_return("not supported output format '$output_format'"); 
  #TODO check for non-sense combinations of output format and decoration types. eg. property and newick.
 
  my $samtrak = Samtrak::Site->new(); # handles connection to the solaris db.

  my $sam = $samtrak->manager_class_of('Sample::Attribute');
  if (scalar(@attributes) == 0) { # get all.
    @attributes = @{$sam->get_all};  
  }
  else {
    my @attributes_tmp = ();
    for (my $i = 0; $i <= $#attributes; $i++) {
      my $candidate = $sam->get_by_abbreviation($attributes[$i]) || $sam->get_by_id_or_nickname($attributes[$i]);
      unless (defined $candidate) {
        print STDERR "Warning: unknown attribute " . $attributes[$i] . " will be ignored\n";
      }
      else {
        push @attributes_tmp , $candidate;
      }
    }
    @attributes = @attributes_tmp;
  }

  my $in = $input eq "-" ? Bio::TreeIO->new(-fh => \*STDIN, -format => $input_format) : Bio::TreeIO->new(-file => $input, format => $input_format) 
     or return $self->error_return("could not open input tree file '$input'");
  my $out = $output eq "-" ? Bio::TreeIO->new(-fh => \*STDOUT, -format => $output_format) : Bio::TreeIO->new(-file => ">$output", format => $output_format)
     or return $self->error_return("could not open output tree file '$output'");

  my %samples_data = (); # will cache each sample attributes as these are needed.
  my $sm = $samtrak->manager_class_of('Sample');
  while (my $tree = $in->next_tree) {
    $tree = $self->annotate_tree($tree,$decoration_type,$sm,\%samples_data,\@attributes);
    $tree->set_tag_value('rooted','false') if ($output_format eq 'phyloxml' && !$tree->get_tag_values('rooted'));
    $out->write_tree($tree);
  }
  $in->close();
  if ($output_format eq "phyloxml") {
    $out->_print('</phyloxml>');
    $out->flush();
  }
  $out->close();
  return $self->ok_return();
}

sub annotate_tree {
  my ($self,$tree,$how,$sm,$samples_data,$attributes) = @_;
  my @leaves = grep { $_->is_Leaf } $tree->get_nodes;
  $self->fetch_samples($sm,$samples_data,@leaves);
  foreach my $leaf (@leaves) {
    $self->annotate_leaf($leaf,$how,$samples_data,$attributes);
  }
  return $tree;
}

sub annotate_leaf {
  my ($self,$node,$how,$sample_data,$attributes) = @_;
  my $sample = $$sample_data{$node->id};
  return unless defined $sample;# no data nothing to do here.
  if ($how eq 'property') {
    foreach my $attr (@$attributes) {
      my $values = $attr->values_of($sample);
      next if (scalar(@$values) == 0);
      $node->set_tag_value($attr->abbreviation,map { $attr->externalize($_) } @$values);
    }
  }
  elsif ($#$attributes == 0) {
    my $attr = $$attributes[0];
    my $values = $attr->values_of($sample);
    if (scalar(@$values) > 0) {
     $node->id(join("=",$node->id(),join("&",map { escape_tree_chars($_) } map { $attr->externalize($_) }  @$values)));
    }
  }
  else {
     my @suffixes = ();
     foreach my $attr (@$attributes) {
       my $values = $attr->values_of($sample);
       if (scalar(@$values) > 0) {
         push @suffixes , $attr->abbreviation . "=" . join("&",map { escape_tree_chars($_) } map { $attr->externalize($_) } @$values);
       }
     }
     $node->id(join("\$",$node->id,@suffixes));
  }
}

sub escape_tree_chars {
  my $v = shift;
  uri_escape($v,"():;[],\"\'_%\$\\=\&");
}

sub fetch_samples {
  my ($self,$sm,$samples_data,@leaves) = @_;
  my %by_id = map { $_ => undef } map { $_->id } @leaves;
  my @ids = grep { !exists($$samples_data{$_}) } keys %by_id; # filter out retrieved ones.
  return unless scalar(@ids) > 0;
  %by_id = map { $_ => undef } @ids;
  my @ids_ref = map { \$_ } map { "'$_'" } @ids;
  my ($subsql,$binds) = $sm->get_objects_sql(query => ['tag_values.tag.abbreviation' => \'\'ox_code\'', 'tag_value.value' => \@ids_ref],
                                             select => 'id', require_objects => ['tag_values','tag_values.tag']);
  my $samples = $sm->get_samples(query => [id => [\$subsql]], 
      with_objects => [qw(contributed_sample 
                          contributed_sample.contribution 
                          contributed_sample.contribution.sender 
                          contributed_sample.contribution.owner
                          tag_values tag_values.tag dna_sample)]);
  foreach my $s (@$samples) { 
    $by_id{$s->attribute('ox_code')->get()} = $s;
  }
  my @missing_ids = grep { !defined $by_id{$_} }  keys %by_id; 
  foreach my $id (@missing_ids) {
    print STDERR "Warning: could not find sample with ox_code '$id', will be left alone";
  }
  foreach my $id (keys %by_id) {
    print STDERR "Info: found data for ox_code $id\n";
    $$samples_data{$id} = $by_id{$id};
  }
}


1;
