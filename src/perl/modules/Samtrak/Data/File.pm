package Samtrak::Data::File;

use strict;
use warnings;
use File::Spec::Functions qw(catfile);

use base 'Samtrak::DB::Object';

__PACKAGE__->meta->setup(
  table => 'file',
  description => 'Represents all Physical samples including DNA and cryoisolates.',
  default_cascade_save => 1,
  columns => [
    id => { type => 'serial', primary_key => 1, not_null => 1,
            description => 'Unique numeric indentifer for file entries' },
    name => { type => 'varchar', length => 45, not_null => 1,
            description => 'base file name without the containing directory path' },
    path => { type => 'varchar', length => 255, not_null => 1,
            description => 'enclosing directory name' },
    type => { type => 'varchar', length => 45, not_null => 1,
            description => 'file type bam, refseq, snpomatic ... etc' },
    filesystem => { type => 'varchar', length => 45, not_null => 1,
            description => 'volume, logical space, device or service containg the file such as fuse, nfs, ftp irods'},
    descr => { alias => 'description', type => 'varchar', length => 45, 
            description => 'human readable description of the content' },
    md5 => { type => 'varchar', length => 45, not_null => 1,
             description => 'MD5 hex digest of the contents of the file' },
    creation_date => { type => 'datetime', not_null => 1,
             description => 'when the file was created' },
    algorithm => { type => 'varchar', length => 45,
             description => 'method or piece of software that generated the file' },
    ref_seq => { type => 'integer' ,
             description => 'reference to the file containing the reference sequence to use in analysing the contents of this file' },
    record_added => { type => 'datetime', not_null => 1,
             description => 'when this record was added to the table' },
    is_current =>  { type => 'boolean', not_null => 1,
             description => 'flags the latest updated on the same file' },
  ],

  foreign_keys => [
    reference => {
      class => 'Samtrak::Data::File',
      key_columns => { ref_seq => 'id' },
      relationship_type => 'many to one',
    },
  ],
  relationships => [
    samples => {
      map_class => 'Samtrak::Data::Sample::File',
      map_to => 'sample',
      map_from => 'file',
      type => 'many to many',
    },
    lanes => {
      map_class => 'Samtrak::Data::Sample::File',
      map_to => 'lane',
      map_from => 'file',
      type => 'many to many',
    },
  ],

);


sub full_name {
  my $self = shift;
  catfile($self->path,$self->name);
}

1;

__END__

  `id` int(10) unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(45) NOT NULL COMMENT 'file name without path',
  `path` varchar(255) NOT NULL COMMENT 'full file path without name',
  `type` varchar(45) NOT NULL COMMENT 'file type eg bam, refseq, snpomatic',
  `filesystem` varchar(45) NOT NULL COMMENT 'fuse, nfs, ftp',
  `descr` varchar(45) DEFAULT NULL,
  `size` bigint(20) unsigned NOT NULL,
  `md5` varchar(45) NOT NULL COMMENT 'md5 checksum ',
  `creation_date` datetime NOT NULL,
  `algorithm` varchar(45) DEFAULT NULL COMMENT 'algorithm used to create file',
  `ref_seq` int(10) unsigned DEFAULT NULL COMMENT 'FK file.id for reference DNA sequence files',
  `record_added` datetime NOT NULL COMMENT 'when this record was added',
  `is_current` tinyint(1) NOT NULL DEFAULT '1' COMMENT 'is this record curretn, 1= true',
  PRIMARY KEY (`id`),


+---------------+---------------------+------+-----+---------+----------------+
| Field         | Type                | Null | Key | Default | Extra          |
+---------------+---------------------+------+-----+---------+----------------+
| id            | int(10) unsigned    | NO   | PRI | NULL    | auto_increment |
| name          | varchar(45)         | NO   | MUL | NULL    |                |
| path          | varchar(255)        | NO   |     | NULL    |                |
| type          | varchar(45)         | NO   |     | NULL    |                |
| filesystem    | varchar(45)         | NO   |     | NULL    |                |
| descr         | varchar(45)         | YES  |     | NULL    |                |
| size          | bigint(20) unsigned | NO   |     | NULL    |                |
| md5           | varchar(45)         | NO   |     | NULL    |                |
| creation_date | datetime            | NO   |     | NULL    |                |
| algorithm     | varchar(45)         | YES  |     | NULL    |                |
| ref_seq       | int(10) unsigned    | YES  | MUL | NULL    |                |
| record_added  | datetime            | NO   |     | NULL    |                |
| is_current    | tinyint(1)          | NO   |     | 1       |                |
+---------------+---------------------+------+-----+---------+----------------+
13 rows in set (0.00 sec)

mysql> describe sample2file
    -> ;
+-----------+------------------+------+-----+---------+----------------+
| Field     | Type             | Null | Key | Default | Extra          |
+-----------+------------------+------+-----+---------+----------------+
| id        | int(10) unsigned | NO   | PRI | NULL    | auto_increment |
| sample_id | int(11)          | NO   | MUL | NULL    |                |
| file_id   | int(10) unsigned | NO   | MUL | NULL    |                |
| lane_id   | int(10) unsigned | YES  | MUL | NULL    |                |
+-----------+------------------+------+-----+---------+----------------+
4 rows in set (0.00 sec)

