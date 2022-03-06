library(pacman)
Packages <- c("tidyverse", "ggrepel", "ggsci", "poppr", "ape", "vcfR", "FactoMineR", "factoextra", "ggpubr")
p_load(Packages, character.only = TRUE)


#calculating the transition to transversion ratio of SNP
SNP_mutation_raw <- tibble(REF = rep(c("A", "C", "G", "T"), each = 3), ALT = c("C", "T", "G", "A", "G", "T", "C", "T", "A", "A", "G", "C"), type = rep(c("transversion", "transversion", "transition"), 4)) 
variety_list <- read_delim("/data/projects/peanut_work/data/snp-calling_samtools/snp-calling/variety_list", delim = "\t", col_names = "var")
variety_list$var
SNP_mutation_var <- vector("list", length = length(variety_list$var))
SNP_mutation_var_table <- vector("list", length = length(variety_list$var))

for (i in seq_along(SNP_mutation_var)) {
  x <- read_delim(str_c("/data/projects/peanut_work/data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/", variety_list$var[[i]], "_SNP_mutation_raw_data"), delim = "\t", col_names = c("REF", "ALT")) %>%
    group_by(REF, ALT) %>%
    summarise(count = n()) %>%
    ungroup()
  x_2 <- SNP_mutation_raw %>%
    left_join(x) %>%
    group_by(type) %>%
    summarise(count = sum(count)) %>%
    spread(key = "type", value = "count") %>%
    mutate(ts_tv_ratio = transition/transversion)
  SNP_mutation_var[[i]] <- x_2$ts_tv_ratio
  #SNP_mutation_var_table[[i]] <- x_2 
}



ts_tv_ratio_31_var <- variety_list %>%
  mutate(ts_tv_ratio = round(unlist(SNP_mutation_var), 2))

ts_tv_ratio_31_var_m <- ts_tv_ratio_31_var %>%
  arrange(-ts_tv_ratio) %>%
  mutate(var = if_else(var == "TN9", "TNS9", if_else(var == "Penghu1", "PH1", var)))

write_delim(ts_tv_ratio_31_var_m, "/data/projects/peanut_work/analysis/ts_tv_ratio_31_var_aln", delim = "\t", col_names = FALSE)

#create SNPs that can distinguish 31 accessions
#import SNP data of 31 accessions 
SNP_31var <- read_delim("./data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/31_varieties_merged_filtered_geno_gotten_modified", delim = "\t", col_names = TRUE)

#extract the genotype of two alleles of each identified markers of 31 accessions (keep biallelic variants)
SNP_31_var_homo_backbone <- SNP_31var %>%
  mutate(length_alt = str_length(ALT)) %>%
  filter(length_alt == 1) %>%
  mutate(Chr = str_replace(Chr, "arahy.Tifrunner.gnm2.Arahy.", "Chr")) %>%
  select(-ID, -QUAL, -FILTER, -INFO, -FORMAT, -length_alt) %>%
  gather(key = "var", value = "var_genotype", c(colnames(SNP_31var)[10:40])) %>%
  distinct() %>%
  mutate(var_genotype1 = str_sub(var_genotype, 1, 1), var_genotype2 = str_sub(var_genotype, 3, 3)) 

#keep SNPs that can distinguish 31 accessions
SNP_31_var_homo_hetero_noNA <- SNP_31_var_homo_backbone %>%
  filter(var_genotype1 != "." & var_genotype2 != ".") %>%
  #filter(var_genotype1 == var_genotype2) %>%
  mutate(var_genotype1 = as.double(var_genotype1), var_genotype2 = as.double(var_genotype2)) %>%
  mutate(var_genotype = var_genotype1 + var_genotype2) %>%
  group_by(Chr, POS) %>%
  mutate(sum_var_genotype = sum(var_genotype), count_group = n()) %>%
  filter(sum_var_genotype != 0) %>%
  filter(sum_var_genotype/count_group != 2) %>%
  ungroup() %>%
  select(Chr, POS, REF, ALT, var, var_genotype)

#make the table of SNP info with all SNPs that can distinguish 31 accessions (no info of 31 accessions)
SNP_31_var_allele_uniq <- SNP_31_var_homo_hetero_noNA %>%
  select(Chr, POS, REF, ALT) %>%
  distinct()

#calculate the sum of NA of each marker
SNP_31_var_marker_score <- SNP_31_var_homo_backbone %>%
  select(Chr, POS, var) %>%
  left_join(SNP_31_var_homo_hetero_noNA) %>% 
  select(-REF, -ALT) %>%
  spread(key = "var", value = "var_genotype") %>%
  select(-Chr, -POS) %>%
  mutate(geno_score = rowSums(., na.rm = TRUE)) %>%
  select(geno_score) 

#make the table of allele information of markers that can differentiate 31 accessions (with the allele info of 31 accessions)
SNP_31_var_homo_hetero_list <- SNP_31_var_homo_backbone %>%
  select(Chr, POS, var) %>%
  left_join(SNP_31_var_homo_hetero_noNA) %>% 
  select(-REF, -ALT) %>%
  spread(key = "var", value = "var_genotype") %>%
  #drop_na() %>%
  #select(-Chr, -POS) %>%
  mutate(x = rowSums(is.na(.))) %>%
  mutate(geno_score = SNP_31_var_marker_score$geno_score) %>%
  filter(62-2*x != geno_score) %>%
  left_join(SNP_31_var_allele_uniq) %>%
  gather(key = "var", value = "var_genotype", c(colnames(SNP_31var)[10:40])) %>%
  select(-x, -geno_score) %>%
  #drop_na() %>%
  #mutate(var_genotype = if_else(var_genotype == 0, str_c(REF, REF), if_else(var_genotype == 2, str_c(ALT, ALT), str_c(REF, ALT)))) %>%
  mutate(marker_name = str_c(Chr, "_", POS)) %>%
  mutate(var_genotype = if_else(var_genotype == 1, 0.5, if_else(var_genotype == 2, 1, var_genotype))) %>%
  spread(key = "var", value = "var_genotype") %>%
  select(-Chr, -POS, -REF, -ALT) 

#create the function of establishing the table based on the tolerance of missing values (accessions lacking data)
produce_table_df_per_NA <- function(a, data){
  x <- data %>%
    mutate(x = rowSums(is.na(.))) %>%
    #drop_na() %>%
    filter(x <= a) %>%
    select(-x)
  
  x2 <- x %>%
    select(-marker_name) %>%
    t() %>%
    as_tibble() 
  
  colnames(x2) <- x$marker_name
  varname_table <- tibble(var = colnames(x)[-1])
  
  x3 <- bind_cols(varname_table, x2) 
  x3
}

#create SNP info with the tolerance of missing (3 accessions > 10%, 6 accessions > 20%)
SNP_31_var_marker_20p_NA_f <- produce_table_df_per_NA(6, SNP_31_var_homo_hetero_list)
SNP_31_var_marker_10p_NA_f <- produce_table_df_per_NA(3, SNP_31_var_homo_hetero_list)


write_csv(SNP_31_var_marker_20p_NA_f, "./analysis/SNP_31_var_marker_20p_NA_f_aln.csv", col_names = TRUE)
write_csv(SNP_31_var_marker_10p_NA_f, "./analysis/SNP_31_var_marker_10p_NA_f_aln.csv", col_names = TRUE)

#genetic diversity parameters (He, PIC, MAF)
#create the table for the following grouping used in the calculation of parameters of genetic diversity

#all
peanut_pop_gene_diversity_all <- tibble(var = c("E01001", "E01004", "HL1", "HL2", "India_song", "NS011001", "PH1", "PI109839", "PI118480", "PI118989", "PI145681", "PI153169",
                                            "PI155112", "PI203396", "PI259717", "PI314817", "PI338337", "PI565455", "PI599592", "Red", "TC1", "TN14", "TN15", 
                                            "TN16", "TN17", "TN18", "TNS9", "TNG10", "TNG7", "Vietnam", "Xiamen"),
                                    pop = rep("all", 31))


#botanical varieties
peanut_pop_gene_diversity_botanical <- tibble(var = c("E01001", "E01004", "HL1", "HL2", "India_song", "NS011001", "PH1", "PI109839", "PI118480", "PI118989", "PI145681", "PI153169",
                                       "PI155112", "PI203396", "PI259717", "PI314817", "PI338337", "PI565455", "PI599592", "Red", "TC1", "TN14", "TN15", 
                                       "TN16", "TN17", "TN18", "TNS9", "TNG10", "TNG7", "Vietnam", "Xiamen"),
                                    pop = c("VA", "VA", "SP", "VA", "SP", "VR_F", "VR_F", "VR_F", "VA", "VA", "VR_F", "SP",
                                    "VA", "VR_F", "SP", "VA", "VA", "SP", "VR_F", "VA", "VR_F", "SP", "SP", "VA", "VA", "SP", "SP", "SP", "SP", "SP", "SP"))

#global and local subsets
peanut_pop_gene_diversity_col <- tibble(var = c("E01001", "E01004", "HL1", "HL2", "India_song", "NS011001", "PH1", "PI109839", "PI118480", "PI118989", "PI145681", "PI153169",
                                            "PI155112", "PI203396", "PI259717", "PI314817", "PI338337", "PI565455", "PI599592", "Red", "TC1", "TN14", "TN15", 
                                            "TN16", "TN17", "TN18", "TNS9", "TNG10", "TNG7", "Vietnam", "Xiamen"),
                                    pop = c("global", "global", "local", "local", "global", "local", "local", "global", "global", "global", "global", "global",
                                            "global", "global", "global", "global", "global", "global", "global", rep("local", 10), "global", "global"))

#centers of origin and local cultivars
peanut_pop_gene_diversity_col_2 <- tibble(var = c("PI109839", "PI118480", "PI118989", "PI153169",
                                                "PI155112", "PI203396", "PI314817", "PI338337", "TC1", "TN14", "TN15", 
                                                "TN16", "TN17", "TN18", "TNS9", "TNG10", "TNG7", "HL1", "HL2", "PH1"),
                                        pop = c(rep("center", 8), rep("cultivar", 12)))


#create the function for creating the table with parameters of genetic diversity

the_genetic_diversity_pop_calculation <- function(data, data2){
x <- data %>%
  gather(key = "marker", value = "var_genotype", c(2:3475)) %>%
  mutate(var = if_else(var == "Taichung1", "TC1", if_else(var == "red", "Red", if_else(var == "TN9", "TNS9", if_else(var == "Penghu1", "PH1", var))))) %>%
  left_join(data2) %>%
  mutate(REF = if_else(var_genotype == 0, 2, if_else(var_genotype == 0.5, 1, if_else(var_genotype == 1, 0, var_genotype)))) %>%
  mutate(ALT = if_else(var_genotype == 1, 2, if_else(var_genotype == 0.5, 1, if_else(var_genotype == 0, 0, var_genotype)))) %>%
  mutate(REF_homo = if_else(REF == 2, 1, 0), ALT_homo = if_else(ALT == 2, 1, 0)) %>%
  group_by(pop, marker) %>%
  summarise(REF_sum = sum(REF, na.rm = TRUE), ALT_sum = sum(ALT, na.rm = TRUE), REF_homo_sum = sum(REF_homo, na.rm = TRUE), ALT_homo_sum = sum(ALT_homo, na.rm = TRUE)) %>%
  mutate(allele_sum = REF_sum + ALT_sum, REF_freq = REF_sum/allele_sum, ALT_freq = ALT_sum/allele_sum, REF_homo_freq = REF_homo_sum/(allele_sum/2), ALT_homo_freq = ALT_homo_sum/(allele_sum/2)) %>%
  mutate(major_allele = if_else(REF_freq > ALT_freq, REF_freq, ALT_freq), minor_allele = if_else(REF_freq > ALT_freq, ALT_freq, REF_freq)) %>%
  mutate(Ho = 1 - REF_homo_freq - ALT_homo_freq, He = 1 - (major_allele^2 + minor_allele^2), PIC = 1 - (major_allele^2 + minor_allele^2) - 2*major_allele^2*minor_allele^2) %>%
  mutate(Ho = if_else(REF_homo_sum + ALT_homo_sum == allele_sum/2, 0, Ho)) %>%
  mutate(Fis = 1-(Ho/He)) %>%
  summarise(count = n(), mean_Fis = mean(Fis, na.rm = TRUE), min_Ho = min(Ho), max_Ho = max(Ho), mean_Ho = mean(Ho), median_Ho = median(Ho), min_He = min(He), max_He = max(He), mean_He = mean(He), median_He = median(He), min_PIC = min(PIC), max_PIC = max(PIC), mean_PIC = mean(PIC), median_PIC = median(PIC), min_major_allele = min(major_allele), max_major_allele = max(major_allele), mean_major_allele = mean(major_allele), median_major_allele = median(major_allele)) 
x
}

#test the table with the info of genetic diversity using SNP which can distinguish 31 accessions (6 accessions at most for the missing values)
bind_rows(the_genetic_diversity_pop_calculation(SNP_31_var_marker_20p_NA_f, peanut_pop_gene_diversity_all),
the_genetic_diversity_pop_calculation(SNP_31_var_marker_20p_NA_f, peanut_pop_gene_diversity_botanical),
the_genetic_diversity_pop_calculation(SNP_31_var_marker_20p_NA_f, peanut_pop_gene_diversity_col))



#create the list in vcf format for later analysis (genetic distance and PCA) 
SNP_31_var_marker_list <- SNP_31_var_homo_hetero_list %>%
  select(marker_name) %>%
  mutate(status = "marker")

SNP_31_var_marker_list_20p <- tibble(
  marker_name = colnames(SNP_31_var_marker_20p_NA_f)[-1],
  status = "marker")

#create vcf file for the following analysis
create_marker_vcf <- function(data, data2){
x <- data %>%
  gather(key = "var", value = "genotype", c(10:40)) %>%
  mutate(genotype = str_sub(genotype, 1, 3)) %>%
  spread(key = var, value = genotype) %>%
  mutate(Chr = str_remove(Chr, "arahy.Tifrunner.gnm2.Arahy.")) %>%
  mutate(marker_name = str_c("Chr", Chr, "_", POS)) %>%
  left_join(data2) %>%
  drop_na() %>%
  select(-marker_name, -status) %>%
  mutate(Chr = str_c("arahy.Tifrunner.gnm2.Arahy.", Chr)) %>%
  distinct()

colnames(x)[[1]] <- "#CHROM"

x
}


#with SNPs only having 20% of missing values at most (6 accessions) based on all accessions
var_31_merged_filtered_geno_marker_20p_list <- create_marker_vcf(SNP_31var, SNP_31_var_marker_list_20p) %>%
  mutate(pop = "31") %>%
  split(.$pop) %>%
  map(. %>% select(-pop))



#We then used the raw vcf file to create 5 subsets of raw vcf files
#establish the subset of three botanical varieties and two origins for the calculation of genetic distance
produce_var_subset_fc <- function(data, data2){
  x <- data %>%
    gather(key = "var", value = "genotype", c(10:40)) %>%
    left_join(data2) %>%
    split(.$pop) %>%
    map(. %>% select(-pop)) %>%
    map(. %>% spread(key = var, value = genotype))
  x
}

#create the subset of the raw vcf files 
var_bot_merged_filtered_geno_marker_20p_list <- produce_var_subset_fc(var_31_merged_filtered_geno_marker_20p_list[[1]], peanut_pop_gene_diversity_botanical)

var_ori_merged_filtered_geno_marker_20p_list <- produce_var_subset_fc(var_31_merged_filtered_geno_marker_20p_list[[1]], peanut_pop_gene_diversity_col)

#merge 5 subsets of of vcf files into 1 list
vcf_raw_subset_append_list <- append(append(var_31_merged_filtered_geno_marker_20p_list, var_bot_merged_filtered_geno_marker_20p_list), var_ori_merged_filtered_geno_marker_20p_list)
vcf_subset_name <- c(31, names(var_bot_merged_filtered_geno_marker_20p_list), names(var_ori_merged_filtered_geno_marker_20p_list))

paths_vcf_raw <-
  str_c(
    "./data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/",
    "var_",
    vcf_subset_name,
    "_merged_filtered_geno_marker_20p_list"
  )

pwalk(list(
  vcf_raw_subset_append_list,
  paths_vcf_raw,
  delim = "\t",
  col_names = TRUE
),
write_delim)

#run bash (adding headers on the top of raw vcf, including 31 accessions and their subsets)
system(paste("cd script/", "&& bash make_vcf_final_genetic_diversity.sh", sep = " "))

#import vcf files (the total one and 5 subsets)
list_subset_vcf <- vector("list", length = length(vcf_subset_name))

for (i in seq_along(vcf_subset_name)) {
  list_subset_vcf[[i]] <- read.vcfR(str_c("./data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/var_", vcf_subset_name[[i]], "_merged_filtered_geno_marker_20p_m.vcf"), convertNA = TRUE)
}

#merge vcf files with the vcf file including 31 accessions and calculate their genetic distance
#list_vcf_genetic_distance <- append(peanut.VCF, list_subset_vcf)
gd_var <- vector("list", length(list_subset_vcf))

for (i in seq_along(list_subset_vcf)) {
  x <- vcfR2genlight(list_subset_vcf[[i]])
  gd_matrix <- as.matrix(bitwise.dist(x))
  gd_sum <- sum(rowSums(gd_matrix))
  n_row <- nrow(gd_matrix)
  gd_var[[i]] <- gd_sum/((n_row*n_row-n_row)/2)/2
}

#merge the information as the Table 
parameters_of_g_diversity_31_var <- bind_rows(the_genetic_diversity_pop_calculation(SNP_31_var_marker_20p_NA_f, peanut_pop_gene_diversity_all),
                                              the_genetic_diversity_pop_calculation(SNP_31_var_marker_20p_NA_f, peanut_pop_gene_diversity_botanical),
                                              the_genetic_diversity_pop_calculation(SNP_31_var_marker_20p_NA_f, peanut_pop_gene_diversity_col)) %>%
  mutate(genetic_distance = unlist(gd_var)) 

write_delim(parameters_of_g_diversity_31_var, "./analysis/parameters_of_g_diversity_31_var_Table2", delim = "\t", col_names = TRUE)

#create the table of genetic distance of 31 accessions (Table S2)
dist_table <- bitwise.dist(gl.peanut)
dist_publisehd_test <- round(dist_table, 3)

write_delim(as_tibble(as.matrix(dist_publisehd_test)), "./analysis/dist_publisehd_test_aln", delim = "\t", col_names = TRUE)

##plot phylogenetic tree (Fig. 1)
peanut.VCF <- list_subset_vcf[[1]]

peanut_pop.data <- tibble(AccessID = c("E01001", "E01004", "HL1", "HL2", "India_song", "NS011001", "PH1", "PI109839", "PI118480", "PI118989", "PI145681", "PI153169",
                                       "PI155112", "PI203396", "PI259717", "PI314817", "PI338337", "PI565455", "PI599592", "Red", "TC1", "TN14", "TN15", 
                                       "TN16", "TN17", "TN18", "TNS9", "TNG10", "TNG7", "Vietnam", "Xiamen"),
                          Country = c("later_defined"),
                          State = c("VA", "VA", "SP", "VA", "SP", "VR", "RN", "VR", "VA", "VA", "RN", "SP",
                                    "VA", "RN", "SP", "VA", "VA", "SP", "RN", "VA", "VR", "SP", "SP", "VA", "VA", "SP", "SP", "SP", "SP", "SP", "SP"))
gl.peanut <- vcfR2genlight(peanut.VCF)
ploidy(gl.peanut) <- 2
pop(gl.peanut) <- factor(peanut_pop.data$State, levels = c("SP", "VA", "VR", "RN"))

tree_peanut <- aboot(gl.peanut , tree = "upgma", distance = bitwise.dist, sample = 1000, showtree = F, cutoff = 50, quiet = T)


#https://link.springer.com/content/pdf/10.1007%2F978-0-387-35100-1_4.pdf
#make high-resolution figures
#https://www.r-bloggers.com/2013/03/high-resolution-figures-in-r/
#par(mar = c(bottom, left, top, right)
tiff("/data/projects/peanut_work/analysis/The_distance_of_peanut_20p_NA_04032022_v2_low.tiff", width = 6, height = 4, units = 'in', res = 300)
jpeg("/data/projects/peanut_work/analysis/The_distance_of_peanut_20p_NA_04032022_v2_low.jpeg", width = 6, height = 4, units = 'in', res = 300)
par(mar = c(2, 0.5, 0, 0))
cols_peanut <- c(pal_npg("nrc")(10))[c(2,3,4,1)]
plot.phylo(tree_peanut, cex = 0.8, font = 2, tip.color =  cols_peanut[pop(gl.peanut)])

nodelabels(tree_peanut$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.6, font = 3, xpd = TRUE)

#tricks of adjusting ticks 
#https://r-charts.com/base-r/axes/#rotate-ticks
axis(side = 1, pos = 29, at = c(0, 0.01, 0.02), labels = c("", "", ""), tck = 0)
text(0.01, 28, "0.04", cex = 0.7)
axis(side = 1, at = 0, pos = 29, labels = c(""), tck = 0.01)
axis(side = 1, at = 0, pos = 29, labels = c(""), tck = -0.01)
axis(side = 1, at = 0.02, pos = 29, labels = c(""), tck = 0.01)
axis(side = 1, at = 0.02, pos = 29, labels = c(""), tck = -0.01)


legend(0, 32.5, legend = c("SP", "VA", "VR", "RN"), fill = cols_peanut, border = TRUE, bty = "n", cex = 0.8, horiz = TRUE)

#change the position of title
#https://stackoverflow.com/questions/20355410/adjust-plot-title-main-position/20355606
title(xlab = list("Genetic distance (Fraction of different loci)", cex = 0.8), line = 0)

dev.off()





#perform PCA using 3,474 SNPs
peanut.pca <- glPca(gl.peanut, nf = 3)
barplot(100*peanut.pca$eig/sum(peanut.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

peanut_31_var_PCA_result_geno <- as_tibble(peanut.pca$scores) %>%
  mutate(var = rownames(peanut.pca$scores), pop = gl.peanut$pop) %>%
  mutate(type = "Genotypic evaluation (3,474 SNPs)")


peanut_31_var_PCA_result_geno$pop <- pop(gl.peanut)

peanut.pca$eig
sum(peanut.pca$eig[2])/sum(peanut.pca$eig)

#PCA of phenotypic data of 31 var
peanut_31_var_pop_pheno_PCA <- peanut_pop.data %>%
  select(var = AccessID, pop = State)
peanut_31_var_phenotype_2016 <- read_delim("./data/peanut_31_var_phenotype_2016", delim = "\t", col_names = TRUE)

peanut_31_var_phenotype_2016_m <- peanut_31_var_phenotype_2016 %>%
  left_join(peanut_31_var_pop_pheno_PCA) %>%
  drop_na() %>%
  mutate(market_type = factor(pop, levels = c("SP", "VA", "VR", "RN")))

res.pca_31_var <- PCA(peanut_31_var_phenotype_2016_m[,-c(1,12:13)], graph = FALSE)


peanut_31_var_PCA_result_pheno <- as_tibble(res.pca_31_var$ind$coord) %>%
  select(PC1 = Dim.1, PC2 = Dim.2, PC3 = Dim.3) %>%
  mutate(var = peanut_31_var_phenotype_2016_m$var, pop = peanut_31_var_phenotype_2016_m$market_type) %>%
  mutate(pop = factor(pop, levels = c("SP", "VA", "VR", "RN"))) %>%
  mutate(type = "Phenotypic evaluation (8 agronomic traits)") 

#create vectors of explained variance for the latter figures
explained_var_geno_3PC <- c(str_c("PC1 (", "24.2", "%)"), str_c("PC2 (", "20.8", "%)"), str_c("PC3 (", 8.2, "%)"))
explained_var_pheno_3PC <- c(str_c("PC1 (", "33.0", "%)"), str_c("PC2 (", "23.9", "%)"), str_c("PC3 (", 16.3, "%)"))

PCA_result_geno_pheno <- list(peanut_31_var_PCA_result_geno, peanut_31_var_PCA_result_pheno)
PCA_plot_raw_PC1_2 <- vector("list", length = length(PCA_result_geno_pheno))
PCA_plot_raw_PC1_3 <- vector("list", length = length(PCA_result_geno_pheno))

for (i in seq_along(PCA_result_geno_pheno)) {
  PCA_plot_raw_PC1_2[[i]] <- PCA_result_geno_pheno[[i]] %>%
    ggplot(aes(x=PC1, y=PC2))
  
  PCA_plot_raw_PC1_3[[i]] <- PCA_result_geno_pheno[[i]] %>%
    ggplot(aes(x=PC1, y=PC3))
}



PCA_plot_31_var <- function(data, a){
  data +
  scale_color_manual(values = cols_peanut) +
  geom_point(aes(colour=pop), size=3) +
  #stat_ellipse(aes(colour = pop), level = 0.9, size = 1, show.legend = FALSE) +
  geom_label_repel(aes(label = var, color = pop),
                   box.padding   = 0.3, 
                   point.padding = 0.4,
                   segment.color = 'grey50', show.legend = FALSE) +
  #geom_text(aes(label=var), hjust=0, vjust=0, size = 3) + 
  facet_wrap(~type, nrow = 1, scales = "free") +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  labs(x = a[1], y = a[2]) +
  theme(
    #axis.ticks = element_blank(),
    strip.text.x = element_text(
      colour = "black",
      face = "bold",
      size = 16
    ),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.position="bottom"
  ) + 
  theme(axis.title.y = element_text(color = "black"), axis.title.x = element_text(color = "black"),  plot.margin = margin(0,0.5,0,0.5, "cm")) 
}


#create PC1/PC2 and PC1/PC3 for 3,474 SNPs and 8 phenotypic traits of 31 accessions
PCA_plot_insilico_SNP_31_var_1 <- PCA_plot_31_var(PCA_plot_raw_PC1_2[[1]], explained_var_geno_3PC[c(1:2)])
PCA_plot_phenotype_31_var_1 <- PCA_plot_31_var(PCA_plot_raw_PC1_2[[2]], explained_var_pheno_3PC[c(1:2)])

PCA_plot_insilico_SNP_31_var_2 <- PCA_plot_31_var(PCA_plot_raw_PC1_3[[1]], explained_var_geno_3PC[c(1,3)])
PCA_plot_phenotype_31_var_2 <- PCA_plot_31_var(PCA_plot_raw_PC1_3[[2]], explained_var_pheno_3PC[c(1,3)])

PCA_merge_31_var_1 <- ggarrange(PCA_plot_insilico_SNP_31_var_1, PCA_plot_insilico_SNP_31_var_2, nrow = 1, ncol = 2) +
  theme(plot.margin = margin(0.1,0,0,0, "cm"), plot.background = element_rect(fill = "white", color = "white"))

PCA_merge_31_var_2 <- ggarrange(PCA_plot_phenotype_31_var_1, PCA_plot_phenotype_31_var_2, nrow = 1) +
  theme(plot.margin = margin(0.1,0,0,0, "cm"), plot.background = element_rect(fill = "white", color = "white"))

ggsave("./analysis/PCA_merge_31_var.tiff", PCA_merge_31_var_1, width = 300, height = 180, units = c("mm"), dpi = 600)
ggsave("./analysis/PCA_merge_31_var_2.tiff", PCA_merge_31_var_2, width = 300, height = 180, units = c("mm"), dpi = 600)
ggsave("./analysis/PCA_merge_31_var.jpeg", PCA_merge_31_var_1, width = 300, height = 180, units = c("mm"), dpi = 200)
ggsave("./analysis/PCA_merge_31_var_2.jpeg", PCA_merge_31_var_2, width = 300, height = 180, units = c("mm"), dpi = 200)


#KASP marker development based on SNP calling using the reference genome of two diploid ancestors
SNP_31var_old <- read_delim("./data/snp-calling_samtools/snp-calling/ancestor_genome_snp_calling/31_varieties_merged_filtered_geno_gotten_modified", delim = "\t", col_names = TRUE)
KASP_design_coordinate <- read_delim("./data/snp-calling_samtools/snp-calling/KASP_design_coordinate", delim = "\t", col_names = c("marker_coordinate", "type"))

SNP_31_var_homo_backbone_old <- SNP_31var_old %>%
  mutate(length_alt = str_length(ALT)) %>%
  filter(length_alt == 1) %>%
  mutate(Chr = str_replace(Chr, "arahy.Tifrunner.gnm2.Arahy.", "Chr")) %>%
  select(-ID, -QUAL, -FILTER, -INFO, -FORMAT, -length_alt) %>%
  gather(key = "var", value = "var_genotype", c(colnames(SNP_31var)[10:40])) %>%
  distinct() %>%
  mutate(var_genotype1 = str_sub(var_genotype, 1, 1), var_genotype2 = str_sub(var_genotype, 3, 3))

#only take homozygous alleles
SNP_31_var_homo_hetero_noNA_old <- SNP_31_var_homo_backbone_old %>%
  filter(var_genotype1 != "." & var_genotype2 != ".") %>%
  filter(var_genotype1 == var_genotype2) %>%
  mutate(var_genotype1 = as.double(var_genotype1), var_genotype2 = as.double(var_genotype2)) %>%
  mutate(var_genotype = var_genotype1 + var_genotype2) %>%
  group_by(Chr, POS) %>%
  mutate(sum_var_genotype = sum(var_genotype), count_group = n()) %>%
  filter(sum_var_genotype != 0) %>%
  filter(sum_var_genotype/count_group != 2) %>%
  ungroup() %>%
  select(Chr, POS, REF, ALT, var, var_genotype)

SNP_31_var_allele_uniq_old <- SNP_31_var_homo_hetero_noNA_old %>%
  select(Chr, POS, REF, ALT) %>%
  distinct()

SNP_31_var_marker_score_old <- SNP_31_var_homo_backbone_old %>%
  select(Chr, POS, var) %>%
  left_join(SNP_31_var_homo_hetero_noNA_old) %>% 
  select(-REF, -ALT) %>%
  spread(key = "var", value = "var_genotype") %>%
  select(-Chr, -POS) %>%
  mutate(geno_score = rowSums(., na.rm = TRUE)) %>%
  select(geno_score) 

SNP_31_var_homo_hetero_list_old <- SNP_31_var_homo_backbone_old %>%
  select(Chr, POS, var) %>%
  left_join(SNP_31_var_homo_hetero_noNA_old) %>% 
  select(-REF, -ALT) %>%
  spread(key = "var", value = "var_genotype") %>%
  #drop_na() %>%
  #select(-Chr, -POS) %>%
  mutate(x = rowSums(is.na(.))) %>%
  mutate(geno_score = SNP_31_var_marker_score_old$geno_score) %>%
  filter(62-2*x != geno_score) %>%
  left_join(SNP_31_var_allele_uniq_old) %>%
  gather(key = "var", value = "var_genotype", c(colnames(SNP_31var)[10:40])) %>%
  select(-x, -geno_score) %>%
  #drop_na() %>%
  #mutate(var_genotype = if_else(var_genotype == 0, str_c(REF, REF), if_else(var_genotype == 2, str_c(ALT, ALT), str_c(REF, ALT)))) %>%
  mutate(marker_name = str_c(Chr, "_", POS)) %>%
  mutate(var_genotype = if_else(var_genotype == 1, 0.5, if_else(var_genotype == 2, 1, var_genotype))) %>%
  spread(key = "var", value = "var_genotype") %>%
  select(-Chr, -POS, -REF, -ALT) 

#create the list of markers from SNPs without any missing value
nona_31var_PIC_old <- produce_table_df_per_NA(0, SNP_31_var_homo_hetero_list_old) %>%
  gather(key = "marker", value = "var_genotype", c(2:1231)) %>%
  mutate(var = if_else(var == "Taichung1", "TC1", if_else(var == "red", "Red", if_else(var == "TN9", "TNS9", if_else(var == "Penghu1", "PH1", var))))) %>%
  mutate(REF = if_else(var_genotype == 0, 2, if_else(var_genotype == 0.5, 1, if_else(var_genotype == 1, 0, var_genotype)))) %>%
  mutate(ALT = if_else(var_genotype == 1, 2, if_else(var_genotype == 0.5, 1, if_else(var_genotype == 0, 0, var_genotype)))) %>%
  mutate(REF_homo = if_else(REF == 2, 1, 0), ALT_homo = if_else(ALT == 2, 1, 0)) %>%
  group_by(marker) %>%
  summarise(REF_sum = sum(REF, na.rm = TRUE), ALT_sum = sum(ALT, na.rm = TRUE), REF_homo_sum = sum(REF_homo, na.rm = TRUE), ALT_homo_sum = sum(ALT_homo, na.rm = TRUE)) %>%
  mutate(allele_sum = REF_sum + ALT_sum, REF_freq = REF_sum/allele_sum, ALT_freq = ALT_sum/allele_sum) %>%
  mutate(major_allele = if_else(REF_freq > ALT_freq, REF_freq, ALT_freq), minor_allele = if_else(REF_freq > ALT_freq, ALT_freq, REF_freq)) %>%
  mutate(PIC = 1 - (major_allele^2 + minor_allele^2) - 2*major_allele^2*minor_allele^2) %>%
  select(marker, PIC)

#create the raw vcf file including SNPs without missing values 
nona_31var_SNP_old_vcf_raw <- SNP_31var_old %>%
  mutate(marker = str_c(Chr, "_", POS)) %>%
  left_join(nona_31var_PIC_old) %>%
  drop_na() %>%
  select(-marker, -PIC)

#create the raw vcf file including 29 SNPs for KASP marker development
nona_31var_KASP_old_vcf_raw <- SNP_31var_old %>%
  mutate(marker_coordinate = POS) %>%
  #mutate(marker_coordinate = as.double(marker_coordinate)) %>%
  left_join(KASP_design_coordinate) %>%
  drop_na() %>%
  select(-marker_coordinate, -type)

colnames(nona_31var_SNP_old_vcf_raw)[[1]] <- "#CHROM"
colnames(nona_31var_KASP_old_vcf_raw)[[1]] <- "#CHROM"

write_delim(nona_31var_SNP_old_vcf_raw, "./data/snp-calling_samtools/snp-calling/ancestor_genome_snp_calling/nona_31var_SNP_old_vcf_raw.tsv", delim = "\t", col_names = TRUE)
write_delim(nona_31var_KASP_old_vcf_raw, "./data/snp-calling_samtools/snp-calling/ancestor_genome_snp_calling/nona_31var_KASP_old_vcf_raw", delim = "\t", col_names = TRUE)

#create vcf file of homozygous SNPs without missing value (Supplementary File S1) and 29 SNP markers for KASP development
system(paste("cd data/snp-calling_samtools/snp-calling/ancestor_genome_snp_calling/", "&& cat header_ancient_genome_vcf nona_31var_SNP_old_vcf_raw.tsv > nona_31var_SNP_ancient.vcf", sep = " "))
system(paste("cd data/snp-calling_samtools/snp-calling/ancestor_genome_snp_calling/", "&& cat header_ancient_genome_vcf nona_31var_KASP_old_vcf_raw > 31_var_KASP_29markers_F.vcf", sep = " "))

peanut_31_var_KASP.VCF <- read.vcfR("./data/snp-calling_samtools/snp-calling/ancestor_genome_snp_calling/31_var_KASP_29markers_F.vcf", convertNA = TRUE)

gl_31var_KASP.peanut <- vcfR2genlight(peanut_31_var_KASP.VCF)
ploidy(gl_31var_KASP.peanut) <- 2
pop(gl_31var_KASP.peanut) <- factor(peanut_pop.data$State, levels = c("SP", "VA", "VR", "RN"))

peanut_31var_KASP.pca <- glPca(gl_31var_KASP.peanut, nf = 3)
barplot(100*peanut_31var_KASP.pca$eig/sum(peanut_31var_KASP.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

peanut_31var_KASP.pca.scores <- as_tibble(peanut_31var_KASP.pca$scores) %>%
  mutate(var = rownames(peanut_31var_KASP.pca$scores), pop = gl_31var_KASP.peanut$pop) %>%
  mutate(type = "Genotypic evaluation (29 SNPs)")


explained_var_geno_KASP_2PC <- c(str_c("PC1 (", round(sum(peanut_31var_KASP.pca$eig[1])/sum(peanut_31var_KASP.pca$eig)*100, 1), "%)"), str_c("PC2 (", round(sum(peanut_31var_KASP.pca$eig[2])/sum(peanut_31var_KASP.pca$eig)*100, 1), "%)"))

PCA_plot_raw_PC1_2_KASP <- peanut_31var_KASP.pca.scores %>%
  ggplot(aes(x=PC1, y=PC2))

PCA_plot_31_var_29_insilico_SNP <- PCA_plot_31_var(PCA_plot_raw_PC1_2_KASP, explained_var_geno_KASP_2PC[c(1:2)]) 

ggsave("./analysis/PCA_plot_31_var_29_insilico_SNP.tiff", PCA_plot_31_var_29_insilico_SNP, width = 240, height = 192, units = c("mm"), dpi = 600)

ggsave("./analysis/PCA_plot_31_var_29_insilico_SNP.jpeg", PCA_plot_31_var_29_insilico_SNP, width = 240, height = 192, units = c("mm"), dpi = 200)
