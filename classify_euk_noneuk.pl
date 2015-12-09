#!/usr/bin/env perl

# sujai.kumar@ed.ac.uk

while (<>) {
    @F = split /\t/;
    if (/,Eukaryota,/) {
        $proteins{$F[0]}{euk} += $F[11];
    } else {
        $proteins{$F[0]}{noneuk} += $F[11];
    }
}

foreach $id (keys %proteins) {
    if (defined ( $proteins{$id}{euk}    )) { $euk_score    = $proteins{$id}{euk}    } else { $euk_score = 0 }
    if (defined ( $proteins{$id}{noneuk} )) { $noneuk_score = $proteins{$id}{noneuk} } else { $noneuk_score = 0 }
    if    (   $euk_score/($euk_score + $noneuk_score) >= 0.9) { print $id . "\teuk\n" }
    elsif ($noneuk_score/($euk_score + $noneuk_score) >= 0.9) { print $id . "\tNoneuk\n" }
    else  { print $id . "\tNotSure\n" }
}

