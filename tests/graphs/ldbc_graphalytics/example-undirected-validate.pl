#!/usr/bin/env perl

# Generate the list of edges to validate for test_graphlog.cpp Graphlog::EdgeReader

use strict;
use warnings;

open(my $f, "/tmp/log.txt"); # create a log (input.txt) in the generator of the form [src: <id>, dst: <id>, weight: <value>]
my $matches = 0;

while(my $line = <$f>){
    if($line =~ /\[src: (\d+), dst: (\d+), weight: (-?[\d\.]+)\]/){
        print("validate_edge(reader, $1, $2, $3);\n");
        $matches ++;
    }
}

close($f);

print("matches: $matches\n");