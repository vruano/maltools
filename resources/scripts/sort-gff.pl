=pod

sort-gff.pl - sorts features within sequences (chromosomes) so that these are listed in ascending order based in start, stop and type.


=cut

use strict;
use warnings;

my $line = <STDIN>;
while (defined $line) {
 print STDERR "LINE: $line\n";
 if ($line !~ /^\s*#/) {
   my ($seq,$db,$type,$start,$stop,@others) = split(/\t/,$line);
   my $current_seq = $seq;
   my @features = ();
   while (1) {
     push @features,[$start,$stop,$type,$line];
     $line = <STDIN>;
     last unless defined $line && $line !~ /^#/; 
     ($seq,$db,$type,$start,$stop,@others) = split(/\t/,$line);
     last unless $seq eq $current_seq;
   }
   @features = sort { $a->[0] <=> $b->[0] 
                   || $a->[1] <=> $b->[1] 
                   || $a->[2] cmp $b->[2] } @features;
   foreach my $f (@features) {
     print $f->[3];
   }
 }
 else {
   print $line;
   $line = <STDIN>;
 }
}
