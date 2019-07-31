# Analysis of ITS_apple_replant_45samples

```
######## Modified from GLBRC/Bonito's Lab Workflow for ITS Sequences###########
```
```
# all the analyses results are stored in /mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/Workflow2_ITS1_OpenReferenceClust/
# ITS raw sequence data stored on HPCC 
# the analyses were conducted in new centos7 HPCC system (ssh dev-intel18)
# quality filtering was conducted using USEARCH10
# primer/adapter removal was performed using cutadapt 1.17 with Python 2.7.15.
# trimming was conducted with length of 200 bp 
# OTU table was built using both open reference of UNITE V.7.2 and de novo clustering with 97% of similarity
# Classifying was conducted using "consensus_taxonomy.txt" file produced by CONSTAX
# further analyses were conducted using QIIME 1.9.1 pipeline (you have to install it in your home directory because it cannot be loaded either in old centos6 or new centos7 HPCC system)
# unfortunately, taxonomy file (consensus_taxonomy.txt) produced from CONSTAX is biom unfriendly so it cannot be added into OTU table biom file using command: biom add-metadata
# add taxonomy to OTU table can be conducted/modified using R software, Excel, and Phyloseq 

apple replant soil ITS:  45 total files
Copied apple replant ITS sequences to working space for analysis
ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/raw_reads/
```

# QUALITY CHECKING AND PRE-FILTERING

a) count read numbers
```
for fastq in ../rawreads/*.fastq
do wc -l $fastq
done > reads_raw.counts
```
b) produce reads quality graphs using FastQC
```
mkdir stats

cat ../rawreads/*R1_001.fastq > raw_reads_R1.fastq; cat ../rawreads/*R2_001.fastq > raw_reads_R2.fastq

/mnt/home/bintarti/FastQC/fastqc raw_reads_R1.fastq raw_reads_R2.fastq -o stats && rm -rf raw_reads_R1.fastq raw_reads_R2.fastq
```

## MERGE PAIRED READS 
```
mkdir mergedfastq

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_mergepairs ../rawreads/*R1*.fastq -reverse ../rawreads/*R2*.fastq -fastq_maxdiffs 10 -fastq_minmergelen 50 -relabel @ -fastqout mergedfastq/merged.fastq
```
### output
```
6818052  Pairs (6.8M)
   4601240  Merged (4.6M, 67.49%)
   2563437  Alignments with zero diffs (37.60%)
   1854285  Too many diffs (> 10) (27.20%)
    234054  No alignment found (3.43%)
         0  Alignment too short (< 16) (0.00%)
    128473  Merged too short (< 50)
    578828  Staggered pairs (8.49%) merged & trimmed
    202.81  Mean alignment length
    281.03  Mean merged length
      1.04  Mean fwd expected errors
      0.85  Mean rev expected errors
      0.23  Mean merged expected errors
```
## CHECK SEQUENCE QUALITY OF MERGED SEQS USING USEARCH AND VSEARCH
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 mergedfastq/merged.fastq -output stats_eestats2_USEARCH.txt

/mnt/home/bintarti/vsearch-2.8.4-linux-x86_64/bin/vsearch -fastq_stats mergedfastq/merged2.fastq -log stats_results_VSEARCH.txt
```
### output
```
4601240 reads, max len 484, avg 281.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4555392( 99.0%)    4597024( 99.9%)    4601216(100.0%)
   100    4272714( 92.9%)    4421917( 96.1%)    4440687( 96.5%)
   150    4117291( 89.5%)    4371651( 95.0%)    4438364( 96.5%)
   200    3944832( 85.7%)    4245254( 92.3%)    4360838( 94.8%)
   250    3770636( 81.9%)    4091828( 88.9%)    4238271( 92.1%)
   300     677838( 14.7%)     794362( 17.3%)     867429( 18.9%)
   350     173782(  3.8%)     218912(  4.8%)     256250(  5.6%)
   400      15972(  0.3%)      22710(  0.5%)      29544(  0.6%)
   450       1319(  0.0%)       2299(  0.0%)       3547(  0.1%)
```
## REMOVE PRIMER AND ADAPTERS WITH CUTADAPT, CHECK THE STATS
```
#upgrade pip
pip install --upgrade pip

#install latest version of cutadapt 1.17
pip install --user cutadapt
```
```
CS1-ITS1 (fwd): 5’- CTTGGTCATTTAGAGGAAGTAA – 3’ (EMP/Smith and Peay 2014)
CS2-ITS2 (rev): 5’- GCTGCGTTCTTCATCGATGC – 3’ (EMP/Smith and Peay 2014)
Reverse complement of adapter-reverse primer: GCATCGATGAAGAACGCAGC

/mnt/home/bintarti/.local/bin/cutadapt -g CTTGGTCATTTAGAGGAAGTAA -a GCATCGATGAAGAACGCAGC -f fastq -n 2 --discard-untrimmed --match-read-wildcards -o cut_merged.fastq mergedfastq/merged.fastq > cut_adpt_results.txt

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_eestats2 cut_merged.fastq -output cutdapt_eestats2_USEARCH.txt

#output: 4600167 reads, max len 463, avg 239.0

Length         MaxEE 0.50         MaxEE 1.00         MaxEE 2.00
------   ----------------   ----------------   ----------------
    50    4384372( 95.3%)    4435805( 96.4%)    4440797( 96.5%)
   100    4243581( 92.2%)    4411100( 95.9%)    4438871( 96.5%)
   150    4046648( 88.0%)    4289781( 93.3%)    4367229( 94.9%)
   200    3862516( 84.0%)    4138727( 90.0%)    4256775( 92.5%)
   250    1450277( 31.5%)    1606508( 34.9%)    1693344( 36.8%)
   300     203478(  4.4%)     248866(  5.4%)     286097(  6.2%)
   350      39477(  0.9%)      52476(  1.1%)      64854(  1.4%)
   400       3066(  0.1%)       4801(  0.1%)       6918(  0.2%)
   450          0(  0.0%)          1(  0.0%)          1(  0.0%)
```

## FILTERING, TRIMMING, QUALITY CHECK
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastq_filter cut_merged.fastq -fastq_maxee 1 -fastq_trunclen 200 -fastq_maxns 0 -fastaout filtered_cut_merged.fa -fastqout filtered_cut_merged.fastq

###output:
   4600167  Reads (4.6M)                    
    328218  Discarded reads length < 200
        76  Discarded read with > 0 Ns
    133191  Discarded reads with expected errs > 1.00
   4138682  Filtered reads (4.1M, 90.0%)

#quality check
/mnt/home/bintarti/FastQC/fastqc filtered_cut_merged.fastq
```
# CLUSTERING AND DENOISING (OTUs vs. ESV)

## DEREPLICATE
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -fastx_uniques filtered_cut_merged.fastq -fastaout derep_filtered_cut_merged.fasta -sizeout

#output: 4138682 seqs, 547383 uniques, 372398 singletons (68.0%)
Min size 1, median 1, max 401391, avg 7.56
```
###########################################################################
```
## GENERATING Exact Sequence Variants (ESV) ALSO CALLED 0-RADIUS OTUs (DENOISING-MAKING ZOTU)

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -unoise3 derep_filtered_cut_merged.fasta -tabbedout unoise_zotus.txt -zotus zotus.fasta

#output:
00:02 190Mb   100.0% 3017 amplicons, 1097757 bad (size >= 8) 
00:39 202Mb   100.0% 3012 good, 5 chimeras                  
00:39 202Mb   100.0% Writing zotus

/mnt/home/bintarti/python_scripts-master/fasta_number.py zotus.fasta OTU_ > zotus_numbered.fasta

## CONSTRUCT ZOTU TABLE-MAP READS BACK TO zotus_numbered.fasta AS DATABASE FILE

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/merged.fastq -db zotus_numbered.fasta -strand plus -id 0.97 -otutabout otu_table_ITS_UNOISE.txt

#output:00.0% Searching merged.fastq, 92.4% matched
4250658 / 4601240 mapped to OTUs (92.4%)  
```
###########################################################################

## OPEN REFERENCE-BASED OTU PICKING AGAINST UNITE FUNGAL ITS DATABASE (v. 7.2) at 97% IDENTITY
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global derep_filtered_cut_merged.fasta -id 0.97 -db /mnt/research/ShadeLab/UNITE_v7.2/sh_refs_qiime_ver7_97_s_01.12.2017.fasta  -strand plus -uc ref_seqs.uc -dbmatched UNITE_reference.fasta -notmatched UNITE_failed_closed.fq

# output: 100.0% Searching, 57.3% matched
```
### SORT BY SIZE AND CLUSTERING-MAKING OTU TABLE at 97% SEQUENCE SIMILARITY (de novo otu-picking)-UPARSE ON SEQUENCE THAT FAILED TO HIT REFERENCE/UNITE
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64  -sortbysize UNITE_failed_closed.fq -fastaout sorted_UNITE_failed_closed.fq

/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -cluster_otus  sorted_UNITE_failed_closed.fq -minsize 2 -otus DENOVO_otus.fasta -uparseout uparse_otus.txt -relabel OTU_

#output: 100.0% 2763 OTUs, 524 chimeras
```
## Combine the rep sets between de novo and reference-based OTU picking
```
cat UNITE_reference.fasta DENOVO_otus.fasta > REPRESENTATIVE_SET.fna

/mnt/home/bintarti/python_scripts-master/fasta_number.py REPRESENTATIVE_SET.fna OTU_ > NUMBERED_REP_SET.fasta
```
## CONSTRUCT OTU TABLE (MAP READS BACK TO REPRESENTATIVE SEQUENCES AS DATABASE FILE)
```
/mnt/research/rdp/public/thirdParty/usearch10.0.240_i86linux64 -usearch_global mergedfastq/merged.fastq -db NUMBERED_REP_SET.fasta -strand plus -id 0.97 -uc OTU_map.uc -otutabout OPEN_REF_OTU_TABLE_ITS.txt

#output: 4283105 / 4373305 mapped to OTUs (97.9%)   
```
## TAXONOMIC CLASSIFICATION USING CONSTAX
```
please refer to how "Running CONSTAX on the MSU HPCC on lab guru : https://my.labguru.com/knowledge/documents/330

# YOU HAVE TO RUN CONSTAX WITHIN THE OLD CENTOS 6 HPCC SYSTEM!
# add "NUMBERED_REP_SET.fasta" into: mnt/home/bintarti/CONSTAX_hpcc/otus
# change the PATH in: mnt/home/bintarti/CONSTAX_hpcc/config
# output directories: outputs, taxonomy_assignments, training_files
# move the output directories to: /mnt/research/ShadeLab/WorkingSpace/Bintarti/ITS_apple_replant/ITS_45_aprep/Workflow2_ITS1_OpenReferenceClust/CONSTAX_OTU
# Use the taxonomy file : outputs/consensus_taxonomy.txt
```
## CONVERT .TXT FILE INTO .BIOM FILE - CONSTAX output is not BIOM friendly, thus do not use this command: biom add_metadata
```
biom convert -i OPEN_REF_OTU_TABLE_ITS.txt -o OPEN_REF_OTU_TABLE_ITS.biom --table-type="OTU table" --to-json

#output file: OPEN_REF_OTU_TABLE_ITS.biom
```
## ALIGN SEQUENCES USING MUSCLE (in CENTOS 6)
```
module load MUSCLE/3.8.31
muscle -in NUMBERED_REP_SET.fasta -out OTUS_aligned.fasta -maxiters 2 -diags1
```

# START THE ANALYSES USING QIIME 1.9.1 PIPELINE

```
# should install QIIME 1.9.1 using miniconda 2 (Python2.7) environment
# download miniconda 2 and put it in home directory
```
##install miniconda 
```
bash Miniconda2-latest-Linux-x86_64
```
##create qiime1 environment and install qiime
```
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
```
##activate qiime1 environment and test qiime installation
```
source activate qiime1
print_qiime_config.py -t
```
## FILTER EXCESS GAPS FROM ALIGNMENT USING QIIME 1.9.1
```
filter_alignment.py -i OTUS_aligned.fasta -o OTUS_filtered_alignment
```
## MAKE PHYLOGENY WITH FASTREE
```
make_phylogeny.py -i OTUS_filtered_alignment/OTUS_aligned_pfiltered.fasta -o OTUS_rep_set.tre
```
## Rarefy OTU table to lowest sequencing depth
```
# summarize OTU table
biom summarize-table -i OPEN_REF_OTU_TABLE_ITS.biom -o OTU_table_sum.txt

# rarefaction
single_rarefaction.py -d 56240 -o rarefied/OTU_single_rare.biom -i OPEN_REF_OTU_TABLE_ITS.biom

# summarize rarefied OTU table
biom summarize-table -i rarefied/OTU_single_rare.biom -o rarefied/OTU_single_rare_sum.txt

# Convert rarefied OTU table file: .biom into .txt for further analyses on R
biom convert -i rarefied/OTU_single_rare.biom -o rarefied/OTU_single_rare.txt --to-tsv

# Convert OPEN_REF_OTU_TABLE_ITS.biom file into .txt for further analyses on R
biom convert -i OPEN_REF_OTU_TABLE_ITS.biom -o rarefied/OPEN_REF_OTU_TABLE_ITS.txt --to-tsv
```