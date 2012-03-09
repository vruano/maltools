package MalariaGEN::AGV::Pipelines::Probgen::Manifest;

use Moose;
extends 'MalariaGEN::AGV::Manifest::Pipeline';

has 'interval_maxsize' => (is => 'ro', isa => 'Int', lazy => 1, builder => '_build_interval_maxsize');

has 'interval_mincount' => (is => 'ro', isa => 'Int', lazy => 1, builder => '_build_interval_mincount'); 

override 'auto_complete' => sub {
    my $self = shift;
    my $result = super();
    my $ref = $self->get_reference()
      or die "could not resolve reference genome location";
    my $intervals = $self->_process_interval_spec($self->get('intervals'));
    if (exists($intervals->{url})) {
        my $url = $intervals->{url};
        $url =~ '^((\w+)://)?(\S+)$';
        my $proto = defined $2 ? $2 : 'file';
        my $file = $3;
        die "yet unsupported url protocol '$proto'" unless $proto eq 'file';
        $ref = $ref->clone(interval_file => $file);
    }
    if (exists($intervals->{minCount}) || exists($intervals->{maxSize})) {
        my @intervals = $ref->create_intervals(min_number => $intervals->{minCount} || 1,
                               max_size => $intervals->{maxSize} || 1000000000);
        @intervals = map { join("",$$_[0],':',$$_[1],'-',$$_[2]) } @intervals;
        $intervals = $self->_process_interval_spec(\@intervals);
    }
    $self->content->{'intervals'} = $intervals;
};

sub get_intervals {
    return $_[0]->get('intervals');
}

sub _process_interval_spec {
    my $self = shift;
    my $intervals = shift || { maxSize => 10000, minCount => 1};
    my $reftype = ref($intervals);
    if (ref($intervals) eq "") {
        if ($intervals =~ /^d+$/) {
            $intervals = { minCount => int($intervals) };
        }
        elsif ($intervals =~ /^(\d+)bp$/) {
            $intervals = { maxSize => int($1) };
        }
        elsif ($intervals =~ '^\w+://') {
            $intervals = { url => $intervals }
        }
        else {
            die "cannot handle intervals spec '$intervals'";
        }
    }
    elsif (ref($intervals) eq "ARRAY") {
        my $key = 'interval0000001';
        my @bad = grep { ! $_ =~ /\S:\d+-\d+/ } @$intervals;
        die "invalid interval specification '$bad[0]'" if @bad;
        $intervals = { map { ($key++ => $_) } @$intervals};
    }
    elsif (ref($intervals) eq "HASH") {
        my @bad = grep { ! $_ =~ /\S:\d+-\d+/ } values %$intervals;
        die "invalid interval specification '$bad[0]'" if @bad;        
    }
    else {
        die "cannot handle interval specification reftype '" . ref($intervals) . "'";
    }
    return $intervals;
}

1;