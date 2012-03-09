package Samtrak::Data::QCFlagStat;

use strict;
use warnings;

use base 'Samtrak::DB::Object';

__PACKAGE__->meta->setup(
  table => 'QC_flagstat',
  description => 'Represents all Physical samples including DNA and cryoisolates.',
  default_cascade_save => 1,
  columns => [
    id => { type => 'serial', primary_key => 1, not_null => 1,
            description => 'Unique numeric indentifer for flag_stat entries' },
    sample_id => { type => 'integer', not_null => 1,
            description => 'numeric identifier to the related sample' },
    ox_code => { type => 'varchar', not_null => 1, length => 15,
            description => 'Oxford Code for the related sample' },
    QC_failure => { type => 'integer', not_null => 1 },
    duplicates => { type => 'integer', not_null => 1 },
    total_reads => { type => 'integer', not_null => 1 },
    mapped_reads => { type => 'integer', not_null => 1 },
    paired_reads => { type => 'integer', not_null => 1 },
    mapped_as_pair => { type => 'integer', not_null => 1 },
    read1 => { type => 'integer', not_null => 1 },
    read2 => { type => 'integer', not_null => 1 },
    singleton_mapped => { type => 'integer', not_null => 1 },
    mate_pairs_mapped => { type => 'integer', not_null => 1 },
    mate_seq_mismap => { type => 'integer', not_null => 1 },
    added => { type => 'timestamp', not_null => 1 },
    current => { type => 'boolean' , not_null => 1 },
  ],

  foreign_keys => [
    sample => {
      deprecated => 'now using an alternative tag-value property with abbrevition <emphasis>storage_location</emphasis>',
      description => 'Not in use any more',
      class => 'Samtrak::Data::Sample',
      key_columns => { sample_id => 'id' },
      relationship_type => 'many to one',
    },
  ],

#  relationships => [],
);

1;
