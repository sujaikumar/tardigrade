use strict;
use warnings;
my %tax;
my %hit;
my %unc;
my %name;
my %exp;
open DIAMOND, "zcat nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.gz |";
while (<DIAMOND>) {
    if (/^(\S+)\t(\S+).*\t.*,(.*?,.*?,.*?),.*?,.*?,$/) {
        if (not defined $hit {$1}) {
            $hit {$1} = $2;
            $tax {$1} = $3;
        }
    }        
}
close DIAMOND;
open NAMES, "</exports/blast_db/uniref90.fasta";
while (<NAMES>) {
    if (/>(\S+)\s+(.*)/) {
        $name {$1} = $2;
    }
}
close NAMES;
open EXP, "<edi_abundance.tsv";
while (<EXP>) {
    if (/^(\S+).*\t(\S+\t\S+)$/) {
        if (not defined $exp {$1}) {
            $exp {$1} = $2;
        }
    }        
}
close EXP;
open UNC, "<nHd.2.3.1.aug.proteins.fasta.tg.default.maker.proteins.final.fasta.1e-5.tophit.txt";
while (<UNC>) {
    if (/^(\S+)\t(\S+)/) {
        if (not defined $unc {$1}) {
            $unc {$1} = $2;
        }
    }        
}
close UNC;

####

while (<>) {
    chomp;
    print $_;
    print "\t" . ( defined $exp {"$_"} ? $exp {"$_"} : "NO_EXPRESSION" );
    print "\t" . ( defined $unc {"$_"} ? $unc {"$_"} : "NO_UNC_HIT" );
    print "\t" . ( defined $hit {"$_"} ? $hit {"$_"} : "NO_UNIREF90_HIT" );
    if (defined $hit {"$_"}) {
        my $uniref = $hit {$_};
        print "\t" . $name { $uniref } . "\t" . $tax { $_ } . "\n";
    }
    else {
        print "\t.\t.\n";
    }
}
