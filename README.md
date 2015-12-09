# tardigrade
Scripts and relevant processed data files for  Boothby et al 2015 and Koutsovoulos et al 2015 tardigrade genome papers

## PacBio verification of euk-noneuk junctions

Also see my comment at http://www.igs.umaryland.edu/labs/hotopp/2015/12/05/quick-look-at-the-two-manuscripts-on-tardigrade-lgt/

Note: This does not take into account our additional RNAseq data or our genomic (lack of) coverage data as described in http://biorxiv.org/content/early/2015/12/01/033464 . For now, I just wanted to see how many eukaryote-noneukaryote adjacent gene pairs from UNC's own list of genes could be verified using their own PacBio data.

Because I did not have access to the UNC team's assignment of genes as Metazoan/Bacterial/Fungal etc, I had to redo the taxonomic assignment of all genes so it is possible that my final number will be slightly different from theirs. Rather than pick gene pairs at random, I selected *all* instances where a eukaryote-like gene was next to a non-eukaryote like gene on the same scaffold.

Only 294 gene pairs matched these requirements (euk-noneuk adjacent on same scaffold), and of those, only 10 were fully spanned by a PacBio scaffold.

Repeating this for Metazoa-Nonmetazoa pairs, 713 such pairs remain (some of the euk-euk pairs become meta-nonmeta, and so the number of meta-nonmeta pairs went up). Of those, only 26 were verified by PacBio scaffolds.

Details of scripts/commands are below (I hope I've not made any mistakes as it is 430am! If someone could double check my numbers I'd be very grateful).

In summary, although UNC claim that the PacBio data verify their assembly, the reality (assuming my analysis is correct) is that only 10 euk-noneuk (or 26 metazoa-nonmetazoa) gene pairs are actually verified by their own data.

### Overview of steps

1. Classify protein predictions as being eukaryote / Noneukaryote / NotSure

    a. blastp protein fasta against **uniref90** database 1e-5 (this is a very lenient threshold, because we want to be sure we're not missing hits just because uniref doesn't have sequences from closely related species, can change this)

    b. Look up uniref90 hits and assign taxonomy ID. Create taxon hierarchy for each hit

    c. Classify each protein/gene as **eukaryote** / **Noneukaryote** depending on whether the sum of bitscores of all hits to eukaryotes is >90% sum of bitscores to all hits (and vice-versa). Use **NotSure** for all the others

2. Pull out adjacent genes that are eukaryote-Noneukaryote from GFF file

    a. Pull out mRNA coordinates from GFF file

    b. Append taxonomy classification eukaryote Noneukaryote NotSure to each GFF entry

    c. Pull out outer coordinates of adjacent eukaryote-Noneukaryote genes (NotSures allowed in between)
    
3. How many eukaryote-Noneukaryote regions are verified by PacBio scaffolds

    a. Use UNCs PacBio/UNC mapping file to extract UNC assembly scaffold regions that have pacbio scaffolds mapping

    b. Use bedtools to see how many euk-noneuk regions (from step 2 above) are covered completely by pacbio verified UNC regions (3a)
    
### Implementation

1. Classify protein predictions as being eukaryote / Noneukaryote / NotSure

    a. blastp protein fasta against **uniref90** database 1e-5 (this is a very lenient threshold, because we want to be sure we're not missing hits just because uniref doesn't have sequences from closely related species, can change this)

            wget http://weatherby.genetics.utah.edu/seq_transf/tg.default.maker.proteins.final.fasta.gz
            gunzip tg.default.maker.proteins.final.fasta.gz
            
            # we used Diamond blastp for speed
            /exports/software/diamond/diamond-v0.7.9/diamond blastp \
                -d /exports/blast_db/uniref90.0.79 \
                -q tg.default.maker.proteins.final.fasta \
                -a tg.default.maker.proteins.final.fasta.uniref90.blastp \
                -e 1e-5 -t /dev/shm -c 8 -p 32
            # -t -c -p are performance parameters
    
    b. Look up uniref90 hits and assign taxonomy ID (needs uniref100.taxlist). Create taxon hierarchy for each hit
    
            perl daa_to_tagc.pl \
                /exports/blast_db/uniref100.taxlist \
                tg.default.maker.proteins.final.fasta.uniref90.blastp.daa
            
            paste \
                tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc \
                <(cut -f13 tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc | cut -f1 -d ";" \
                | perl taxid_parents_list.pl -) \
            >tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc.taxonhierarchy

                
    c. Classify each protein/gene as **eukaryote** / **Noneukaryote** depending on whether the sum of bitscores of all hits to eukaryotes is >90% sum of bitscores to all hits (and vice-versa). The 90% threshold can be changed. Use **NotSure** for all the others

            perl classify_euk_noneuk.pl \
            < tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc.taxonhierarchy \
            > tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.euknoneuk

2. Pull out adjacent genes that are eukaryote-Noneukaryote from GFF file

    a. Pull out mRNA coordinates from GFF file

            wget http://weatherby.genetics.utah.edu/seq_transf/tg.default.final.gff.gz
            zgrep -P "\tmRNA\t" tg.default.final.gff.gz >tg.default.final.mRNA.gff
    
    b. Append taxonomy classification eukaryote Noneukaryote NotSure to GFF entry
    
            perl -e '
                open CLASS, "<tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.classified";
                while (<CLASS>) { $class{$1} = $2 if /^(\S+)\t(.*euk)/ }
            
                open GFF, "<tg.default.final.mRNA.gff";
                while (<GFF>) {
                    chomp;
                    if (/ID=(.*?);/ and exists $class{$1}) {
                        print $_ . "taxclass=$class{$1}\n";
                    } else {
                        print $_ . "taxclass=NotSure\n";        
                    }
                }
            ' | sort -k1,1V -k4,4n >tg.default.final.mRNA.euknoneuk.gff
    
    c. Pull out outer coordinates of adjacent eukaryote-Noneukaryote genes (NotSures allowed in between). i.e. if there is a pair of euk-noneuk genes next to each other on same UNC scaffold, then make a bed file for each such pair with region from start of gene 1 to end of gene 2.

            perl adjacent_difftaxa_regions.pl \
            < tg.default.final.mRNA.euknoneuk.gff \
            > UNC.gff.euknoneuk.bed
            
            wc -l UNC.gff.euknoneuk.bed
            294
    
3. How many eukaryote-Noneukaryote regions are verified by PacBio scaffolds

    a. Use UNCs PacBio/UNC mapping file to extract UNC assembly scaffold regions that have pacbio scaffolds mapping
    
            wget http://weatherby.genetics.utah.edu/seq_transf/pacbio/h_dujardini_illvspb_assemblies.m5.gz
            zcat h_dujardini_illvspb_assemblies.m5.gz \
            | perl -lne '@F=split/\s+/;print join("\t",@F[0..15])' \
            | cut -f1,3,4 | perl -plne 's/\|/_/;s/\/\S+//;' \
            >UNC_pacbio_verified_regions_ALL.bed

    b. Use bedtools to see how many euk-noneuk regions (from 2 above) are covered completely by pacbio verified UNC regions (3a)

            /exports/software/bedtools/bedtools-2.25.0/bin/intersectBed -wa -f 1 \
                -a UNC.gff.euknoneuk.bed \
                -b UNC_pacbio_verified_regions_ALL.bed | uniq | wc -l
            10
