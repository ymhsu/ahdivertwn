#!/bin/bash
# Program:
# Make the complete vcf by merging rows of headers (starting with ##) and the genotyping result of 284 accessions using 14 KASP markers







echo "mergeing header"

cat <(gunzip -c ../data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/31_varieties_cultivated_aln_passed_geno.vcf.gz | grep "##") ../data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/genotype_core_pop_temp_vcf > ../data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/284_var_KASP_result_14markers.vcf
