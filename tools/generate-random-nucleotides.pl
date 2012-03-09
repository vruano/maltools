#!/usr/bin/env perl -w

use strict;
use warnings;

use Getopt::Long;


my @NUCS = (qw(A C G T));
my @seq_names = ();

my @seq_length = ();

GetOptions("sn=s" => \@seq_names, "sl=s" => \@seq_length);

@seq_names = ('seq1') if (scalar(@seq_names) == 0);
@seq_length = (100) if (scalar(@seq_length) == 0);

my $init_length = scalar(@seq_length);
while (scalar(@seq_length) < scalar(@seq_names)) {
  my $idx = (scalar(@seq_length) + 1) % $init_length; 
  push @seq_length, $seq_length[$idx];
}

my $seq_count = scalar(@seq_names);
for (my $i = 0; $i < $seq_count; $i++) {
  print STDOUT ">" . $seq_names[$i] . "\n";
  my $slen = $seq_length[$i];
  my $seq = "";
  for (my $j = 0; $j < $slen; $j++) {
     $seq .= $NUCS[int(rand(4))];
  }
  my $line_start = 0;
  while ($line_start < $slen) {
     print STDOUT substr($seq,$line_start,(($line_start + 60) <= $slen) ? 60 : ($slen - $line_start)), "\n";
     $line_start += 60;
  }
}



