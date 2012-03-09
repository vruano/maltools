use strict;
use warnings;

while (my $line = <STDIN>) {
  chomp $line;
  my ($chr,$pos,$ref,$cand,$typ,$cas1,$cas2) = split(/\t/,$line);
  my $cas = join(",",grep { $_ } ($cas1,$cas2));
  my @cas = grep { $_ !~ /\*/ } split(/,/,$cas);
  unless (@cas) {
   print "NoShow\t1.000\n";
   next;
  }
  elsif (@cas == 1) {
   print cause(split(/:/,$cas[0])) . "\t1.000\n";
   next;
  }
  #print STDERR "$line\n";
  #print STDERR "XXXXX",scalar(@cas),join("|",@cas),"\n";
  my %cas = (map { split(/:/,$_) } @cas);
  print STDERR $line,"\n" if (grep { $_ =~ /:/ } keys %cas);
  my $total = 0;
  $total += $_ foreach (values(%cas));
  print join ("\n", map { cause($_,$line) . "\t" . sprintf ("%.4f",$cas{$_}/$total) } keys (%cas)  ),"\n";
}

sub cause {
   my $code = shift;
   my $line = shift;
   if ($code eq "C") {
     return "CQ";
   }
   elsif ($code eq "X") {
     return "SQ";
   }
   elsif ($code eq "M") {
     return "MMQ";
   }
   elsif ($code eq "D") {
     return "MaxDepth";
   }
   elsif ($code eq "d") {
     return "MinDepth";
   }
   elsif ($code eq "W") {
     return "SnpDensity";
   }
   elsif ($code eq "Q") {
     return "MapQual";
   }
   elsif ($code eq "G") {
     return "AroundIndel";
   }
   else {
     return $code;
   }
}

