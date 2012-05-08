package Maltools::Vcf::Expression;

use strict;
use warnings;
use Moose;
use Carp;

has 'vcf' => (is => 'ro', isa => 'Vcf', required => 1);
has 'expr' => (is => 'ro', isa => 'Str', required => 1);
has 'code' => (is => 'ro', isa => 'Str', lazy => 1, builder => '_code_builder');


sub _code_builder {
  my $self = shift;
  my $expr = $self->expr;
  my $vcf = $self->vcf;
  # Infos:
  my $infos = $vcf->{header}->{INFO};
  foreach my $id (keys %$infos) {
    $expr =~ s/\@\{$id\}/split(\/,\/,\$_info_->{$id})/g;
    $expr =~ s/\b$id\b/\$_info_->{$id}/g;
  }
  my $formats = $vcf->{header}->{FORMAT};
  foreach my $id (keys %$formats) {
    $expr =~ s/\b$id\((\S+)\)/\$_gt_->{'$1'}->{$id}/g;
  }
  my $filters = $vcf->{header}->{FILTER}; 
  foreach my $id (keys %$filters) {
    $expr =~ s/\b$id\b/grep { \$_ eq '$id' } @\$_filter_/g;
  }
  $expr =~ s/\$v\b/\$_var_/g;
  return $expr;
}

sub evaluate {
  my $self = shift;
  my $code = $self->code;
  my $_var_ = shift;
  my $_info_ = $_var_->{INFO};
  my $_filter_ = $_var_->{FILTER};
  my $_gt_ = $_var_->{gtypes};
  my $result = eval $code;
  $@ and die "issues evaluation expression code $code on variant";
  return $result;
}




1;
