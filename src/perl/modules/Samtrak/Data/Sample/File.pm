package Samtrak::Data::Sample::File;

use strict;
use warnings;

use base 'Samtrak::DB::Object';

__PACKAGE__->meta->setup(
  table => 'sample2file',
  description => 'Implements the relationship between files and samples',
  default_cascade_save => 1,
  columns => [
    id => { type => 'serial', primary_key => 1, not_null => 1,
            description => 'Unique numeric indentifer for sample-file entries' },
    sample_id => { type => 'integer', not_null => 1,
            description => 'numeric identifier to the related sample' },
    file_id => { type => 'integer', not_null => 1,
            description => 'reference identifier to the related file' },
    lane_id => { type => 'integer', not_null => 1,
            description => 'reference identifier to the related lane' },
  ],

  foreign_keys => [
    sample => {
      description => 'reference to the corresponding sample',
      class => 'Samtrak::Data::Sample',
      key_columns => { sample_id => 'id' },
      relationship_type => 'many to one',
    },

    lane => {
      description => 'reference to the corresponding lane',
      class => 'Samtrak::Data::Lane',
      key_columns => { lane_id => 'id' },
      relationship_type => 'many to one',
    },
    
    file => {
      description => 'reference to the corresponding file',
      class => 'Samtrak::Data::File',
      key_columns => { file_id => 'id' },
      relationship_type => 'many to one',
    },
    
  ],

);

1;

