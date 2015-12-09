#!/usr/bin/env perl

# sujai.kumar@ed.ac.uk

while (<>) {
    @F = split /\t/;
    if (/,Metazoa,/) {
        $proteins{$F[0]}{metazoa} += $F[11];
    } else {
        $proteins{$F[0]}{nonmetazoa} += $F[11];
    }
}

foreach $id (keys %proteins) {
    if (defined ( $proteins{$id}{metazoa}    )) { $metazoa_score    = $proteins{$id}{metazoa}    } else { $metazoa_score = 0 }
    if (defined ( $proteins{$id}{nonmetazoa} )) { $nonmetazoa_score = $proteins{$id}{nonmetazoa} } else { $nonmetazoa_score = 0 }
    if    (   $metazoa_score/($metazoa_score + $nonmetazoa_score) >= 0.9) { print $id . "\tMetazoa\n" }
    elsif ($nonmetazoa_score/($metazoa_score + $nonmetazoa_score) >= 0.9) { print $id . "\tNonMetazoa\n" }
    else  { print $id . "\tNotSure\n" }
}

