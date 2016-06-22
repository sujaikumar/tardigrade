## Goal

To determine what percentage of genes are putative HGTs in the Arakawa et al 2016 Tardigrade genome sequencing PNAS letter
http://www.pnas.org/content/113/22/E3057.full

## Overview

Arakawa et al 2016 sequenced a single Hypsibius dujardini individual after extensive steps to remove bacteria on it's surface.
They used genomic read coverage and the presence of RNA seq reads (from 35 separate libraries) to estimate the amount of contamination
in the Boothby et al 2015 genome.

Arakawa et al report that:

> 7,135 contigs (31.7% of all contigs), including the longest 11 contigs, are contaminated under a rather
> conservative estimate (less than ×1 coverage in genomic reads, or without even a single hit of RNA-Seq reads
> from any of the 35 datum). Only 1,771 putative HGT genes remain after removing these contigs, and their percentage
> within the genome, 4.47%, is in line with other eukaryotic genomes.

This 4.47% figure is reported in the PNAS letter but is nowhere to be found in the bioRxiv detailed description.

Using their own sequence coverage table, I found that the number of contaminant contigs with coverage < 1 was 7499 contigs
(not 7135 as reported by them).

These 7499 contigs account for 5648 of the 6663 putative HGT genes in Boothby et al 2015. i.e. the remaining 
putative HGT are (6663 - 5648 =) 1015 genes. 1015 out of Boothby et al's 39532 genes is 2.57%, not 4.47%.

I'm happy to be corrected if I've made any mistakes in pulling out contigs with coverage <1 from the Arakawa table.

## Detailed Methods

1. Arakawa et al provide genomic coverage data in a PDF table at
https://figshare.com/articles/Validation_of_assembled_contigs_of_Hypsibius_dujardini_genome_of_Boothby_et_al_using_ultra_low_input_sequencing/2300251/1

2. Convert this pdf to a text file:

        pdftohtml report.pdf report.html

  pdftohtml creates several files. actual data is in [reports.html](reports.html) (available in this repo)

        awk 'NR>=82 && NR <= 115643' reports.html \
        | perl  -pne 's/<br>\n/\t/; s/<hr>//; s/<.*?>//g; s/Page\s+\d+\s+//; s/\s+/\n/g; s/,//g' \
        | grep '\S' | paste - - - - - - \
        | awk '{print $1$6"\t"$2"\t"$3"\t"$4"\t"$5}' >arakawa.cov.txt

  [arakawa.cov.txt](arakawa.cov.txt) available in this repo.
  
  Apologies for the ugliness of this hack - but it was the quickest way to get the data into a parsable text file from the PDF.
  You can check that it is accurate - all the contigs are there and each row has 5 columns:
  ScaffoldName, Length, Mapped bases, Mean coverage, Standard deviation

3. Get all scaffolds/contigs with mean coverage < 1 (the tardigrade genome has coverage ~100 X)

        awk '$4<1' arakawa.cov.txt | wc -l
        7499

4. Calculate the span of these 7499 scaffolds:

        awk '{if ($4<1){s+=$2}}END{print s}' arakawa.cov.txt
        70228131

  A nice result is that this 70.2 Mbp span almost exactly matches the 68.9 Mbp span that we called contaminant based on
  our sequencing in Table 1 of Koutsovoulos et al 2016 http://www.pnas.org/content/113/18/5053.full . If anything, we
  have been too generous and not called tan extra 1.3 Mbp as contamination, and we note in our paper that there is
  some remaining contamination in our version of the assembly.

5. Get the genes from the 6663 putative HGTs from Boothby et al that are on these 7499 contigs

        awk '$4<1' arakawa.cov.txt | cut -f1 | fgrep.pl -f - 6663_HGT_genes.ids -d "-" | wc -l
        5648

  Explanation: fgrep.pl (available in this repo) can pull out matches from a file by splitting the file to be searched using a 
  delimiter (-d "=" in this case)

6. Calculate putative HGT percentage

    Remove 5648 genes which come from contaminant contigs, from the set of 6663 putative HGTs.
    Calculate percentage out of 39,532 gene predictions:
    
        (6663 - 5648)/39532 = 0.0256754 = 2.57%
