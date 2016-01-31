Hard non-eukaryotic fHGT candidates
===================

Goal: To identify which of the 23,021 genes in the Edinburgh gene annotation set nHd.2.3.1 are hard HGT candidates

1. Best blast hits to non-eukaryotic sequences
2. On same scaffold as eukaryotic sequences
3. Gene expression above 0.1 TPM

### Overview

1. Classify protein predictions as being eukaryote / Noneukaryote / NotSure
    a. blastp protein fasta against **uniref90** database 1e-5 (this is a very lenient threshold, because we want to be sure we're not missing hits just because uniref doesn't have sequences from closely related species, can change this)
    b. Look up uniref90 hits and assign taxonomy ID. Create taxon hierarchy for each hit
    c. Classify each protein/gene as **eukaryote** / **Noneukaryote** depending on whether the sum of bitscores of all hits to eukaryotes is >90% sum of bitscores to all hits (and vice-versa). Use **NotSure** for all the others

2. Pull out noneuk genes that are on same scaffold as euk genes using GFF file
    a. Pull out mRNA coordinates from GFF file
    b. Append taxonomy classification eukaryote Noneukaryote NotSure to each GFF entry
    c. Pull out genes classified as noneuk that are on the same scaffold as a euk gene

3. Use gene expression level to further select genes
    a. Append expression level from Kallisto abundance.tsv results file to gene list
    b. Select genes with gene expression >0.1 TPM

### Implementation

1. Classify protein predictions as being eukaryote / Noneukaryote / NotSure

    a. blastp protein fasta against **uniref90** database 1e-5 (this is a very lenient threshold, because we want to be sure we're not missing hits just because uniref doesn't have sequences from closely related species, can change this)

            # we used Diamond blastp for speed
            /exports/software/diamond/diamond-v0.7.9/diamond blastp \
                -d /exports/blast_db/uniref90.0.79 \
                -q nHd.2.3.1.aug.proteins.fasta \
                -a nHd.2.3.1.aug.proteins.fasta.uniref90.blastp \
                -e 1e-5 -t /dev/shm -c 8 -p 32
            # -t -c -p are performance parameters
    
    b. Look up uniref90 hits and assign taxonomy ID (needs uniref100.taxlist). Create taxon hierarchy for each hit
    
            perl daa_to_tagc.pl \
                /exports/blast_db/uniref100.taxlist \
                nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa
            # this perl script takes the very large Diamond blastp output .daa file
            # and returns a tabular blast format file with the taxid as the last column

            paste \
                nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc \
                <(cut -f13 nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc | cut -f1 -d ";" \
                | perl taxid_parents_list.pl -) \
            >nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy

            # taxid_parents_list.pl takes a taxid and returns the full NCBI taxonomy path to the root
            # paste simply appends that to the blast output as an extra column
                
    c. Classify each protein/gene as **eukaryote** / **Noneukaryote** depending on whether the sum of bitscores of all hits to eukaryotes is >90% sum of bitscores to all hits (and vice-versa). The 90% threshold can be changed. Use **NotSure** for all the others

            perl classify_euk_noneuk.pl \
            < nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy \
            > nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.euknoneuk

            # this file has two columns - col1=protein name, col2=euk/Noneuk/NotSure

2. Pull out noneuk genes that are on same scaffold as euk genes using GFF file

    a. Pull out mRNA coordinates from GFF file

            grep -P "\tmRNA\t" nHd.2.3.1.aug.gff > nHd.2.3.1.aug.mRNA.gff
    
    b. Append taxonomy classification eukaryote Noneukaryote NotSure to GFF entry
    
            perl -e '
                open CLASS, "zcat nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.euknoneuk.gz |";
                while (<CLASS>) { $class{$1} = $2 if /^(\S+)\t(.*euk)/ }
            
                open GFF, "<nHd.2.3.1.aug.mRNA.gff";
                while (<GFF>) {
                    chomp;
                    if (/ID=(.*?);/ and exists $class{$1}) {
                        print $_ . ";taxclass=$class{$1}\n";
                    } else {
                        print $_ . ";taxclass=NotSure\n";        
                    }
                }
            ' | sort -k1,1V -k4,4n > nHd.2.3.1.aug.mRNA.euknoneuk.gff

    c. Pull out genes classified as noneuk that are on the same scaffold as a euk gene

            fgrep.pl \
                -f nHd.cov2cov.ge_10.names.txt \
                nHd.2.3.1.aug.mRNA.euknoneuk.gff \
            | grep -v NotSure | perl -plne 's/^(\S+).*=(\S+)/$1\t$2/' | sort | uniq | cut -f1 | uniq -c | awk '$1>1' \
            | awk '{print $2}' | fgrep.pl -f - nHd.2.3.1.aug.mRNA.euknoneuk.gff | grep "Noneuk" \
            | perl -plne 's/.*ID=//'| cut -f1 -d ";" >216_euknoneuk_same_scaffold.ge_10.geneIDs.txt

            # Explanation:
            # fgrep.pl
            #     is like fgrep, but more efficient for large numbers of search terms)
            #     only selects those genes from nHd.2.3.1.aug.mRNA.euknoneuk.gff which are on scaffolds with
            #     read coverage >=10 (i.e. 23,021 genes - 41 genes = 22,980 genes)
            # grep -v NotSure
            #     removes genes we can't assign as euk/noneuk unequivocally (leaving 12,747 genes).
            # perl -plne 's/^(\S+).*=(\S+)/$1\t$2/'
            #     creates two columns for each gene: scaffold_name\teuk/Noneuk
            # sort | uniq | cut -f1 | uniq -c | awk '$1>1' | awk '{print $2}'
            #     prints names of 162 scaffolds that have BOTH euk and noneuk genes on them
            #     (scaffolds with only one kind of gene will be discarded by awk '$1>1' step)
            # fgrep.pl .... grep "Noneuk"
            #     pulls out those GFF mRNA entries which are Noneuk on scaffolds with both euk+noneuk genes
            # perl -plne 's/.*ID=//'| cut -f1 -d ";"
            #     extracts the IDs of these genes (count 216) and prints them to a file

3. Use gene expression level to further select genes

    a. Append expression level from Kallisto edi_abundance.tsv results file (and best blast hits to UNC proteins/uniref90) to gene list

            perl append_abundance_blast_hits.pl \
              216_euknoneuk_same_scaffold.ge_10.geneIDs.txt \
            > 216_euknoneuk_same_scaffold.ge_10.geneIDs.table.txt
            
    b. Select genes with gene expression >0.1 TPM
    
            awk '$3>0.1' 216_euknoneuk_same_scaffold.ge_10.geneIDs.table.txt | wc -l
            196
            
        i.e. 196 genes were hard HGT candidates on the basis of the criteria above (blast hits 
        unequivocally matched noneuk, present on same scaffold as unequivocally euk genes, and
        with RNAseq abundance >0.1 TPM)


Eukaryotic non-metazoan fHGT candidates
===================

Goal: To identify which of the 23,021 genes in the Edinburgh gene annotation set nHd.2.3.1 are euk non-metazoan HGT candidates

(similar to Hard HGT above, but instead of non-euk genes on same scaffold as a eukaryotic scaffold, we look for eukaryotic-nonmetazoan candidates that are on the same scaffold as a metazoan gene)

1. Best blast hits to eukaryotic non-metazoan sequences
2. On same scaffold as metazoan sequences
3. Gene expression above 0.1 TPM

### Overview

1. Classify protein predictions as being Metazoa / Eukaryotic-Nonmetazoa

    Steps a. and b. and c. are the same as for hard HGT candidates

    a. blastp protein fasta against **uniref90** database 1e-5 (this is a very lenient threshold, because we want to be sure we're not missing hits just because uniref doesn't have sequences from closely related species, can change this)
    b. Look up uniref90 hits and assign taxonomy ID. Create taxon hierarchy for each hit
    c. Classify each protein/gene as **eukaryote** / **Noneukaryote** depending on whether the sum of bitscores of all hits to eukaryotes is >90% sum of bitscores to all hits (and vice-versa). Use **NotSure** for all the others
    d. Extract only the genes that are Eukaryotic, and reclassify these as metazoa-nonmetazoa (similar to steps a, b, c)    

2. Pull out euk-nonmetazoa genes that are on same scaffold as metazoa genes using GFF file
    a. Pull out mRNA coordinates from GFF file (same as above)
    b. Append taxonomy classification Euk-Nonmetazoa / Metazoa / NotSure (or not relevant as they are noneukaryotic) to each GFF entry
    c. Pull out genes classified as euk-nonmetazoa that are on the same scaffold as a metazoa gene

3. Use gene expression level to further select genes
    a. Append expression level from Kallisto abundance.tsv results file and other info (best uniref90 blast hit) to gene list
    b. Select genes with gene expression >0.1 TPM

### Implementation

1. Classify protein predictions as being Metazoa / Eukaryotic-Nonmetazoa

    Steps a. and b. and c. are the same as for hard HGT candidates

    a. blastp protein fasta against **uniref90** database 1e-5 (this is a very lenient threshold, because we want to be sure we're not missing hits just because uniref doesn't have sequences from closely related species, can change this)

    b. Look up uniref90 hits and assign taxonomy ID. Create taxon hierarchy for each hit

    c. Classify each protein/gene as **eukaryote** / **Noneukaryote** depending on whether the sum of bitscores of all hits to eukaryotes is >90% sum of bitscores to all hits (and vice-versa). Use **NotSure** for all the others
    
    d. Extract only the genes that are Eukaryotic for sure, and reclassify these as metazoa-nonmetazoa (similar to steps a, b, c)
    
            zgrep ",Eukaryota," nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.gz \
            | fgrep.pl -v -f <(zgrep Noneuk nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.euknoneuk.gz | cut -f1) - \
            | perl classify_meta_nonmeta.pl >nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.noNoneuk.metanonmeta

        A previous version of this analysis (as reported in the second version http://biorxiv.org/content/early/2015/12/13/033464 of the biorxiv paper http://dx.doi.org/10.1101/033464) selected only the Eukaryota blast hits from uniref90 but did not remove proteins that had been marked as unequivocally non Eukaryotic . Some non-eukaryotic proteins can have low quality hits to Eukaryotic sequences in uniref which resulted in 24 nHd.2.3.1. proteins being double counted in the list of non-eukaryotic proteins as well as the list of eukaryotic non-metazoan proteins.
        
2. Pull out euk-nonmetazoa genes that are on same scaffold as metazoa genes using GFF file

    a. Pull out mRNA coordinates from GFF file (same as above)

    b. Append taxonomy classification Euk-Nonmetazoa / Metazoa / NotSure (or not relevant as they are noneukaryotic) to each GFF entry

            perl -e '
                open CLASS, "nHd.2.3.1.aug.proteins.fasta.uniref90.blastp.daa.tagc.taxonhierarchy.noNoneuk.metanonmeta";
                while (<CLASS>) { $class{$1} = $2 if /^(\S+)\t(.*Metazoa)/ }
            
                open GFF, "<nHd.2.3.1.aug.mRNA.gff";
                while (<GFF>) {
                    chomp;
                    if (/ID=(.*?);/ and exists $class{$1}) {
                        print $_ . ";taxclass=$class{$1}\n";
                    } else {
                        print $_ . ";taxclass=NotSure\n";        
                    }
                }
            ' | sort -k1,1V -k4,4n > nHd.2.3.1.aug.mRNA.noNoneuk.metanonmeta.gff

    c. Pull out genes classified as euk-nonmetazoa that are on the same scaffold as a metazoa gene

            fgrep.pl \
                -f nHd.cov2cov.ge_10.names.txt \
                nHd.2.3.1.aug.mRNA.noNoneuk.metanonmeta.gff \
            | grep -v NotSure | perl -plne 's/^(\S+).*=(\S+)/$1\t$2/' | sort | uniq | cut -f1 | uniq -c | awk '$1>1' \
            | awk '{print $2}' | fgrep.pl -f - nHd.2.3.1.aug.mRNA.noNoneuk.metanonmeta.gff | grep "NonMetazoa" \
            | perl -plne 's/.*ID=//'| cut -f1 -d ";" >385_noNoneuk_metanonmeta_same_scaffold.geneIDs.txt

        (The version of this analysis at http://biorxiv.org/content/early/2015/12/13/033464 reported 24 extra proteins - [Supplemental_File_9_409_EukaryoticNonMetazoa_HGT_candidates_in_the_Edinburgh_assembly](http://biorxiv.org/highwire/filestream/9185/field_highwire_adjunct_files/8/033464-9.txt) which had been double counted - see Note at 1.d above)
        
3. Use gene expression level to further select genes

    a. Append expression level from Kallisto edi_abundance.tsv results file (and best blast hits to UNC proteins/uniref90) to gene list

            perl append_abundance_blast_hits.pl \
              385_noNoneuk_metanonmeta_same_scaffold.ge_10.geneIDs.txt \
            > 385_noNoneuk_metanonmeta_same_scaffold.ge_10.geneIDs.table.txt
            
    b. Select genes with gene expression >0.1 TPM
    
            awk '$3>0.1' 385_noNoneuk_metanonmeta_same_scaffold.geneIDs.table.txt | wc -l
            369
            
        i.e. 369 genes were soft HGT candidates on the basis of the criteria above (blast hits 
        unequivocally matched eukaryotic nonMetazoa, present on same scaffold as unequivocally metazoan genes, and
        with RNAseq abundance >0.1 TPM)

