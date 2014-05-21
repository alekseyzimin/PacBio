package UnionFind;

use Exporter;
use strict;
use warnings;

BEGIN {
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    @ISA = ('Exporter');
    $VERSION = 0.0.1;
    @EXPORT = qw();
    %EXPORT_TAGS = ( );
    @EXPORT_OK = @EXPORT;
}
our @EXPORT_OK;

sub new {
    my ($class) = @_;
    return bless({}, $class);
}

sub __find {
    my ($tree) = @_;
    my @reroot;

    while(my $p = $$tree{PARENT}) {
        push(@reroot, $tree);
        $tree = $p;
    }

    foreach my $t (@reroot) {
        $$t{PARENT} = $tree;
    }

    return $tree;
}

sub union {
    my ($self, $obj1, $obj2) = @_;
    my ($tree1, $tree2) = ($$self{$obj1}, $$self{$obj2});
    my ($root1, $root2);

    if(defined($tree1)) {
        $root1 = __find($tree1);
    } else {
        $$self{$obj1} = $root1 = $tree1 = { SELF => $obj1, 
                                            RANK => 0, PARENT => undef};
    }

    return if($obj1 eq $obj2);

    if(defined($tree2)) {
        $root2 = __find($tree2);
    } else {
        $$self{$obj2} = $root2 = $tree2 = { SELF => $obj2,
                                            RANK => 0, PARENT => undef};
    }

    if($$root1{RANK} > $$root2{RANK}) {
        $$root2{PARENT} = $root1;
    } elsif($$root1{RANK} < $$root2{RANK}) {
        $$root1{PARENT} = $root2;
    } elsif($$root1{SELF} ne $$root2{SELF}) {
        $$root2{PARENT} = $root1;
        $$root1{RANK}++;
    }
}

sub find {
    my ($self, $obj) = @_;
    my $tree = $$self{$obj};

    if(!defined($tree)) {
        $$self{$obj} = { SELF => $obj, RANK => 0, PARENT => undef };
        return $obj;
    }

    my $root = __find($tree);
    return $$root{SELF};
}

sub keys {
    my ($self) = @_;

    return keys %$self;
}

sub components {
    my $self = shift @_;
    my %comps;
    my %avoidh;
    
    if(@_) {
        foreach my $e (@_) {
            my $root = $$self{$e} || next;
            $avoidh{__find($root)->{SELF}} = 1;
        }
    }

    while(my ($obj, $tree) = each %$self) {
        my $root = __find($tree);
        next if $avoidh{$$root{SELF}};
        my $comp = $comps{$$root{SELF}};
        $comp = $comps{$$root{SELF}} = [] unless defined($comp);
        push(@$comp, $obj);
    }

    return values(%comps);
}

1; 

=head1 NAME

UnionFind - Union find data structure

=head1 SYNOPSIS

  use UnionFind;

  my $uf = UnionFind::new();
  $uf->union("a", "b");
  $uf->union("b", "c");
  $uf->union("d", "e");
  $uf->find("c");       # -> "a"
  $uf->components();    # -> [["a", "b", "c"], ["d", "e"]]

=head1 DESCRIPTION

Implements the unionfind data structure with union by rank and path
compression algorithm. It is good at finding connected components in
undirected graphs.

=head1 USAGE

=over

=item new

Create and return a new empty unionfind data structure. Every time an
object is queried (via union or find), the object is added to the data
structure if not already present.

=item union(obj1, obj2)

Union the trees containing obj1 and obj2.

=item find(obj)

Find the root of the tree containing obj.

=item keys

Returns an array of all objects added to the unionfind data structure.

=item components

Returns an array reference of array reference. Each internal array is
a connected component.

=back

=head1 AUTHOR

Guillaume Marcais <gus@math.umd.edu>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Guillaume Marcais

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.6 or,
at your option, any later version of Perl 5 you may have available.

=cut
