# tardigrade
Scripts and relevant processed data files for  Boothby et al 2015 and Koutsovoulos et al 2015 tardigrade genome papers
Results:
-----

Using UNC's own gene model GFF file, only **713** adjacent metazoa-nonmetazoa pairs of genes could be found. Of these, only **26** are fully spanned by UNC's own PacBio scaffolds (using their own mapping file).

If this analysis is repeated using eukaryote-noneukaryote as the classification, the number of adjacent euk-noneuk pairs drops to **294**, of which only **10** are fully spanned by the PacBio scaffolds

I'm putting the scripts and intermediate data files at github.com/sujaikumar/tardigrade so that they can be scrutinised to make sure I haven't done anything silly. I'd be grateful for comments on these analyses before we put them in to our manuscript. Thanks!

Steps:
-----

###Overview

1. Classify protein predictions as being Metazoa / NonMetazoa / NotSure

    a. blastp protein fasta against **uniref90** database 1e-5 (this is a very lenient threshold, because we want to be sure we're not missing hits just because uniref doesn't have sequences from closely related species, can change this)

    b. Look up uniref90 hits and assign taxonomy ID. Create taxon hierarchy for each hit

    c. Classify each protein/gene as **Metazoa** / **NonMetazoa** depending on whether the sum of bitscores of all hits to Metazoa is >90% sum of bitscores to all hits (and vice-versa). The 90% threshold can be changed. Use **NotSure** for all the others

2. Pull out adjacent genes that are Metazoa-NonMetazoa from GFF file

    a. Pull out mRNA coordinates from GFF file

    b. Append taxonomy classification Metazoa NonMetazoa NotSure to GFF entry

    c. Pull out outer coordinates of adjacent Metazoa-NonMetazoa genes (NotSures allowed in between)
    
3. How many Metazoa-NonMetazoa regions are verified by PacBio scaffolds

    a. Use UNCs PacBio/UNC mapping file to extract UNC assembly scaffold regions that have pacbio scaffolds mapping

    b. Use bedtools to see how many meta-nonmeta regions (from 2 above) are covered completely by pacbio verified UNC regions (3a)
    
### Implementation

1. Classify protein predictions as being Metazoa / NonMetazoa / NotSure

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

                
    c. Classify each protein/gene as **Metazoa** / **NonMetazoa** depending on whether the sum of bitscores of all hits to Metazoa is >90% sum of bitscores to all hits (and vice-versa). The 90% threshold can be changed. Use **NotSure** for all the others

            perl classify_meta_nonmeta.pl \
            < tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc.taxonhierarchy \
            > tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.classified

2. Pull out adjacent genes that are Metazoa-NonMetazoa from GFF file

    a. Pull out mRNA coordinates from GFF file

            wget http://weatherby.genetics.utah.edu/seq_transf/tg.default.final.gff.gz
            zgrep -P "\tmRNA\t" tg.default.final.gff.gz >tg.default.final.mRNA.gff
    
    b. Append taxonomy classification Metazoa NonMetazoa NotSure to GFF entry
    
            perl -e '
                open CLASS, "<tg.default.maker.proteins.final.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.classified";
                while (<CLASS>) { $class{$1} = $2 if /^(\S+)\t(.*Metazoa)/ }
            
                open GFF, "<tg.default.final.mRNA.gff";
                while (<GFF>) {
                    chomp;
                    if (/ID=(.*?);/ and exists $class{$1}) {
                        print $_ . "taxclass=$class{$1}\n";
                    } else {
                        print $_ . "taxclass=NotSure\n";        
                    }
                }
            ' | sort -k1,1V -k4,4n >tg.default.final.mRNA.metanonmeta.gff
    
    c. Pull out outer coordinates of adjacent Metazoa-NonMetazoa genes (NotSures allowed in between). i.e. if there is a pair of meta-nonmeta genes next to each other on same UNC scaffold, then make a bed file for each such pair with region from start of gene 1 to end of gene 2.


            # if tg.default.final.mRNA.metanonmeta.gff has two adjacent genes like this:
            
            perl adjacent_difftaxa_regions.pl \
            < tg.default.final.mRNA.metanonmeta.gff \
            > UNC.gff.metanonmeta.bed
            
            wc -l UNC.gff.metanonmeta.bed
            713
    
3. How many Metazoa-NonMetazoa regions are verified by PacBio scaffolds

    a. Use UNCs PacBio/UNC mapping file to extract UNC assembly scaffold regions that have pacbio scaffolds mapping
    
            wget http://weatherby.genetics.utah.edu/seq_transf/pacbio/h_dujardini_illvspb_assemblies.m5.gz
            zcat h_dujardini_illvspb_assemblies.m5.gz \
            | perl -lne '@F=split/\s+/;print join("\t",@F[0..15])' \
            | cut -f1,3,4 | perl -plne 's/\|/_/;s/\/\S+//;' \
            >UNC_pacbio_verified_regions_ALL.bed

    b. Use bedtools to see how many meta-nonmeta regions (from 2 above) are covered completely by pacbio verified UNC regions (3a)

            /exports/software/bedtools/bedtools-2.25.0/bin/intersectBed -wa -f 1 \
                -a UNC.gff.metanonmeta.bed \
                -b UNC_pacbio_verified_regions_ALL.bed | uniq | wc -l
            26

