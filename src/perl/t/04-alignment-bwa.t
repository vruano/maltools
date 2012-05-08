
use Maltools::Alignment::BWA qw(bwa);


bwa(in => shift, out => shift, ref => shift, sort => 1, index => 1);
