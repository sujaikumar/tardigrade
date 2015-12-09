#!/usr/bin/env perl

# sujai.kumar@ed.ac.uk

use strict;
use warnings;

my ($line, $curr_contigid, $curr_taxclass, $curr_start, $curr_end,
    @F,    $prev_contigid, $prev_taxclass, $prev_start, $prev_end);

$prev_contigid = "";

while (<>) {
    $line = $_;
    @F = split /\t/, $line;
    $curr_contigid = $F[0];
    $curr_taxclass = $1 if $line =~ /taxclass=(.*)/;
    $curr_start    = $F[3];
    $curr_end      = $F[4];
    next if $curr_taxclass eq "NotSure";
    #print "INPUT\t$curr_contigid\t$curr_taxclass\t$curr_start\t$curr_end\n";
    if ($curr_contigid ne $prev_contigid) {
        $prev_contigid = $curr_contigid;
        $prev_taxclass = $curr_taxclass;
        $prev_start    = $curr_start;
        $prev_end      = $curr_end;
    }
    elsif ($curr_taxclass ne $prev_taxclass) {
        print "$curr_contigid\t$prev_start\t$curr_end\n";
        $prev_taxclass = $curr_taxclass;
        $prev_start    = $curr_start;
        $prev_end      = $curr_end;
    }
    elsif ($curr_taxclass eq $prev_taxclass) {
        $prev_taxclass = $curr_taxclass;
        $prev_start    = $curr_start;
        $prev_end      = $curr_end;
    }
}
