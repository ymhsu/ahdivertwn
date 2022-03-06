#!/bin/bash
# Program:
# Make the complete vcf by merging rows of headers (starting with ##) and the main part of allele information of 31 accessions or their subsets

file="../data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/vcf_name_list_gd"



for i in $(cat "$file")
do

echo ${i}

cat <(gunzip -c ../data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/31_varieties_cultivated_aln_passed_geno.vcf.gz | grep "##") ../data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/var_${i}_merged_filtered_geno_marker_20p_list > ../data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/var_${i}_merged_filtered_geno_marker_20p_m.vcf

done
