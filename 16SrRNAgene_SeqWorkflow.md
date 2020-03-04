Raw Sequences Analysis of 16S rRNA gene (v4 region) from Apple Root-zone Soil

# Analysis of 16S Miseq Data

Here we used the USEARCH pipeline (v10.0.240_i86linux64) for pre-processing raw sequences data and UPARSE method for OTU clustering (Edgar 2013). Additional analyses were conducted using QIIME1 (1.9.1).

raw sequence data stored on HPCC:
/mnt/research/ShadeLab/Sequence/raw_sequence/AppleReplant/20171113_16S-V4_PE/

Moved/copy apple replant sequences (45 soil samples (F01-F45)) to working space for analysis:
ShadeLab/WorkingSpace/Bintarti/16S_45samples_aprep/20171113_16S-V4_PE/

# Part I: Clustering

usearch v10.0.240_i86linux64, 16.3Gb RAM, 4 cores
(C) Copyright 2013-17 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

## 1) Merge Paired End Reads
```
# decompress the reads
gunzip *.gz

# make directory called "mergedfastq1"
mkdir mergedfastq1

# merge paired end reads
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq1/merged.fq -fastq_merge_maxee 1.0 -tabbedout mergedfastq1/merged.report.txt -alnout mergedfastq1/merged_aln.txt

# -fastq_maxdiffs 10: Allow 10 max differences in overlap region

### Output ###

 2414383  Pairs (2.4M)
   1951063  Merged (2.0M, 80.81%)
    694512  Alignments with zero diffs (28.77%)
    455483  Too many diffs (> 10) (18.87%)
      7837  No alignment found (0.32%)
         0  Alignment too short (< 16) (0.00%)
         0  Exp.errs. too high (max=1.0) (0.00%)
      6572  Staggered pairs (0.27%) merged & trimmed
    246.79  Mean alignment length
    253.05  Mean merged length
      0.49  Mean fwd expected errors
      1.42  Mean rev expected errors
      0.16  Mean merged expected errors
```
## 2) Check Sequence Quality of Merged Seqs
```
mkdir fastq_info
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq1/merged.fq -output fastq_info/eestats.txt

### output ###
1951063 reads, max len 466, avg 253.1

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    1926510( 98.7%)    1949362( 99.9%)    1951052(100.0%)
   100    1884745( 96.6%)    1939990( 99.4%)    1950730(100.0%)
   150    1853669( 95.0%)    1927180( 98.8%)    1949097( 99.9%)
   200    1827941( 93.7%)    1918095( 98.3%)    1947689( 99.8%)
   250    1777104( 91.1%)    1891340( 96.9%)    1938656( 99.4%)
   300        420(  0.0%)        563(  0.0%)        658(  0.0%)
   350         72(  0.0%)        152(  0.0%)        245(  0.0%)
   400         16(  0.0%)         47(  0.0%)        113(  0.0%)
   450          1(  0.0%)          4(  0.0%)          6(  0.0%)
```
## 3) Filter and Truncate the Merged Seqs
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter merged.fq -fastq_maxee 1 -fastq_trunclen 250 -fastqout filtered_merged.fq

### output ###
100.0% Filtering, 96.9% passed
   1951063  Reads (2.0M)                    
      5317  Discarded reads length < 250
     54406  Discarded reads with expected errs > 1.00
   1891340  Filtered reads (1.9M, 96.9%)
```
## 4) Dereplicate Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_merged.fq -fastqout uniques_filtered_merged.fastq -sizeout

### output ###
1891340 seqs, 758945 uniques, 544833 singletons (71.8%)
Min size 1, median 1, max 6947, avg 2.49
```
## 5) Remove Singeltons
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize uniques_filtered_merged.fastq -fastqout nosigs_uniques_filtered_merged.fastq -minsize 2

### output ###
Sorting 214112 sequences
```
## 6) Precluster Sequences
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_fast nosigs_uniques_filtered_merged.fastq -centroids_fastq denoised_nosigs_uniques_filtered_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

### output ###
Seqs  214112 (214.1k)
  Clusters  70493 (70.5k)
  Max size  13158 (13.2k)
  Avg size  19.1
  Min size  2
Singletons  0, 0.0% of seqs, 0.0% of clusters
   Max mem  924Mb
      Time  55.0s
Throughput  3892.9 seqs/sec.
```
## 7) Closed Reference-based OTU Picking Using SILVA_132 Database
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global denoised_nosigs_uniques_filtered_merged.fastq -id 0.97 -db /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -strand plus -uc ref_seqs.uc -dbmatched SILVA_closed_reference.fasta -notmatchedfq failed_closed.fq

### output ###
100.0% Searching, 46.4% matched

Closed reference OTU picking:
Pick OTUs based on the Silva database in the home directory.
Produce some output files - ref_seqs.uc (pre-clustered), SILVA_closed_reference.fasta will be the matched ones, and failed_closed.fq will be used in de novo OTU picking
```
## 8) De novo OTU picking
```
# sort by size
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -sortbysize failed_closed.fq -fastaout sorted_failed_closed.fq

# cluster de novo
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus sorted_failed_closed.fq -minsize 2 -otus denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up

### output ###
13731 OTUs, 12958 chimeras
```
## 9) Combine the Rep Sets Between De novo and SILVA Reference-based OTU Picking
```
cat SILVA_closed_reference.fasta denovo_otus.fasta > FULL_REP_SET.fna
```
## 10) Map 'FULL_REP_SET.fna' Back to Pre-dereplicated Sequences and Make OTU Tables
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -usearch_global merged.fq -db FULL_REP_SET.fna -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

### output ###
1803423 / 1951063 mapped to OTUs (92.4%) 
```
# Part II: Switch to QIIME 1.9.1

## 1) Assign taxonomy to SILVA_132_QIIME_release with uclust 
```
# assignment method: uclust
# --uclust_min_consensus_fraction: [default:0.51]
# --uclust_similarity: [default:0.9]
# --uclust_max_accepts: [default: 3]

assign_taxonomy.py -i FULL_REP_SET.fna -o taxonomy -r /mnt/home/bintarti/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna -t /mnt/home/bintarti/SILVA_132_QIIME_release/taxonomy/16S_only/97/consensus_taxonomy_7_levels.txt
```
## 2) Convert OTU_table.txt. to OTU_table.from_txt_json.biom
```
biom convert -i OTU_table.txt -o OTU_table.biom --table-type="OTU table" --to-json
```
## 3) Add taxonomy to OTU table
```
biom add-metadata -i OTU_table.biom -o OTU_table_tax.biom --observation-metadata-fp=taxonomy/FULL_REP_SET_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```
## 4) Filter non-bacteria/archaea
```
filter_taxa_from_otu_table.py -i OTU_table_tax.biom -o OTU_table_tax_filt.biom -n D_4__Mitochondria,D_3__Chloroplast

# summarize OTU table

biom summarize-table -i OTU_table_tax.biom -o OTU_table_tax_sum.txt

biom summarize-table -i OTU_table_tax_filt.biom -o OTU_table_tax_filt_sum.txt
```
## 5) Rarefaction
```
single_rarefaction.py -d 27716 -o OTU_rarefied.biom -i OTU_table_tax_filt.biom

# summarize OTU table
biom summarize-table -i OTU_rarefied.biom -o OTU_rarefied_sum.txt
```
## 6) Convert and add taxonomy
```
biom convert -i OTU_rarefied.biom -o OTU_rarefied.txt --header-key taxonomy --to-tsv

biom convert -i OTU_table_tax_filt.biom -o OTU_table_tax_filt.txt --header-key taxonomy --to-tsv
```









