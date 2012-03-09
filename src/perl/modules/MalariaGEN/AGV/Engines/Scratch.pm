package MalariaGEN::AGV::Engines::Scratch;

use Moose::Role;

requires 'to_scratch';

requires 'from_scratch';

requires 'in_scratch';

requires 'quota';

requires 'tempdir';

requires 'tempfile';

1;
