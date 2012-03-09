package MalariaGEN::AGV::ToolData::fastq;

use strict;
use warnings;
use Moose;

extends 'MalariaGEN::AGV::ToolData::file';

has '+has_files' => ( default => 1 );

1;
