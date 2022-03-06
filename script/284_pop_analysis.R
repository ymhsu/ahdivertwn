#load necessary libraries
library(pacman)
Packages <- c("tidyverse", "ggrepel", "ggsci", "poppr", "ape", "vcfR", "FactoMineR", "factoextra", "ggpubr")
p_load(Packages, character.only = TRUE)

#generate the vcf files of 284 germplasms (the information of REF and ALT is not correct since they are not important for the following analysis)
genotype_core_pop <- read_csv("./data/genotype_core_pop_284.csv", col_names = TRUE) 

genotype_core_pop_backbone <- genotype_core_pop %>%
  select(marker) %>%
  mutate(Chr_raw = str_sub(marker, 1, 1), Chr_n = str_sub(marker, 2, 3)) %>%
  mutate(Chr_n = as.double(Chr_n)) %>%
  mutate(Chr_n = if_else(Chr_raw == "B", Chr_n + 10, Chr_n)) %>%
  arrange(Chr_n) %>%
  mutate(Chr_n = as.character(Chr_n)) %>%
  mutate(Chr_n = if_else(nchar(Chr_n)==1, str_c("0", Chr_n), Chr_n)) %>%
  mutate(CHROM = str_c("arahy.Tifrunner.gnm2.Arahy.", Chr_n)) %>%
  mutate(POS = str_sub(marker, 5, 20)) %>%
  mutate(ID = ".", REF = "T", ALT = "C", QUAL = 250, FILTER = "PASS", INFO = "MQ0F=0;MQ=60;VDB=0;SGB=-0.693147;DP=2611;DP4=0,2499,0,111;AN=52;AC=2", FORMAT = "GT:PL")

genotype_core_pop_temp <- genotype_core_pop_backbone %>%
  left_join(genotype_core_pop) %>%
  select(-marker, -Chr_raw, -Chr_n) %>%
  mutate(`284` = if_else(`284` == "H", "1/0", `284`)) 


peanut_pop_info <- tibble(
  pop_name = c(str_c(c(1:79), "_VA"), str_c(c(80:169), "_RN"), str_c(c(170:218), "_VR"), str_c(c(219:284), "_SP")))


colnames(genotype_core_pop_temp) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO","FORMAT", peanut_pop_info$pop_name)


write_delim(genotype_core_pop_temp, "./data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/genotype_core_pop_temp_vcf", delim = "\t", col_names = TRUE)


#produce vcf files for PCA
system(paste("cd script/", "&& bash make_vcf_14_KASP_284_pop.sh", sep = " "))

peanut_pop.VCF <- read.vcfR("./data/snp-calling_samtools/snp-calling/cultivated_genome_aln_snp_calling/284_var_KASP_result_14markers.vcf", convertNA = TRUE)

peanut_core_pop.data <- peanut_pop_info %>%
  select(AccessID = pop_name) %>%
  mutate(Country = c("later_defined"), market_type = c(rep("VA", 79), rep("RN", 90), rep("VR", 49), rep("SP", 66)))

gl_284.peanut <- vcfR2genlight(peanut_pop.VCF)
ploidy(gl_284.peanut) <- 2
pop(gl_284.peanut) <- peanut_core_pop.data$market_type

#PCA for 14 KASP genotyping
peanut_284.pca <- glPca(gl_284.peanut, nf = 3)
peanut_284.pca.scores <- as.data.frame(peanut_284.pca$scores)

#calculate the contribution of 2 PC for the total variance of SNP
peanut_284.pca$eig[[1]]/sum(peanut_284.pca$eig)
peanut_284.pca$eig[[2]]/sum(peanut_284.pca$eig)

peanut_284_PCA_result_1 <- as_tibble(peanut_284.pca$scores) %>%
  mutate(var = rownames(peanut_284.pca$scores), pop = pop(gl_284.peanut)) %>%
  mutate(pop = factor(pop, levels = c("SP", "VA", "VR", "RN"))) %>%
  mutate(type = "Genotypic evaluation (14 KASP markers)")

#perform PCA on phenotypic data from the experiment in the fall of 2016
peanut_pop_phenotype_2016 <- read_delim("./data/peanut_pop_phenotype_2016", delim = "\t", col_names = TRUE)

peanut_pop_phenotype_2016_m <- peanut_pop_phenotype_2016 %>%
  select(-pod_length, -pod_width) %>%
  drop_na() %>%
  mutate(market_type = str_sub(var, 1, 2)) %>%
  mutate(market_type = factor(market_type, levels = c("SP", "VA", "VR", "RN")))

res.pca <- PCA(peanut_pop_phenotype_2016_m[,-c(1,10,11)], graph = FALSE)

peanut_284_PCA_result_2 <- as_tibble(res.pca$ind$coord) %>%
  select(PC1 = Dim.1, PC2 = Dim.2, PC3 = Dim.3) %>%
  mutate(var = peanut_pop_phenotype_2016_m$name_v2, pop = peanut_pop_phenotype_2016_m$market_type) %>%
  mutate(pop = factor(pop, levels = c("SP", "VA", "VR", "RN"))) %>%
  mutate(type = "Phenotypic evaluation (8 agronomic traits)")

res.pca$eig

cols_peanut_284 <- c(pal_npg("nrc")(10))[c(2,3,4,1)]


PCA_figure_function <- function(data, a) {
  x <- data %>%
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(colour=pop), size=2) +
    stat_ellipse(aes(colour = pop), level = 0.9, size = 1, show.legend = FALSE) +
    scale_color_manual(values = cols_peanut_284) +
    facet_wrap(~type, nrow = 1, scales = "free") +
    #geom_text(aes(label=var), hjust=0, vjust=0, size = 3) + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    theme_bw() +
    labs(x = str_c("PC1 ", "(", a[1], "%)"), y = str_c("PC2 ", "(", a[2], "%)")) +
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
    )  + 
    theme(axis.title.y = element_text(color = "black"), axis.title.x = element_text(color = "black"),  plot.margin = margin(0,0.5,0,0.5, "cm")) 
  x
}

PCR_plot_phenotype <- PCA_figure_function(peanut_284_PCA_result_2, c(28.4, 24.1))
PCR_plot_KASP <- PCA_figure_function(peanut_284_PCA_result_1, c(36.9, 18.3))

PCA_merge_KASP_phenotype <- ggarrange(PCR_plot_KASP, PCR_plot_phenotype, nrow = 1) +
  theme(plot.margin = margin(0.1,0,0,0, "cm"), plot.background = element_rect(fill = "white", color = "white"))

ggsave("./analysis/PCA_merge_KASP_phenotype.jpeg", PCA_merge_KASP_phenotype, width = 300, height = 180, units = c("mm"), dpi = 200)
ggsave("./analysis/PCA_merge_KASP_phenotype.tiff", PCA_merge_KASP_phenotype, width = 300, height = 180, units = c("mm"), dpi = 600)

#cluster analysis based on the result of 14 KASP markers
gl_284_raw.peanut <- vcfR2genlight(peanut_pop.VCF)
ploidy(gl_284_raw.peanut) <- 2
pop(gl_284.peanut) <- peanut_core_pop.data$market_type

test_cluster <- find.clusters(gl_284_raw.peanut, n.pca = 100, choose.n.clust = FALSE, stat = "BIC")
BIC_plot <- tibble(BIC = test_cluster$Kstat) %>%
  mutate(groups = c(1:n())) %>%
  filter(groups <= 6) %>%
  ggplot() + 
  geom_point(mapping = aes(groups, BIC)) + 
  geom_line(mapping = aes(groups, BIC)) + 
  labs(y = "BIC", x = "Number of groups (K)") +
  theme_bw() +
  theme(
    #panel.background = element_rect(colour = "black", linetype = 1),
    #strip.background = element_rect(colour = "black", linetype = 1),
    strip.text = element_text(
      colour = "black",
      face = "bold",
      size = 16
    ),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    legend.position="none"
  ) + scale_x_continuous(breaks = c(1:6)) +
  theme(axis.title.y = element_text(color = "black"), axis.title.x = element_text(color = "black"),  plot.margin = margin(0,0.5,0,0.5, "cm")) 

ggsave("./analysis/BIC_plot.tiff", BIC_plot, width = 300, height = 180, units = c("mm"), dpi = 600)
ggsave("./analysis/BIC_plot.jpeg", BIC_plot, width = 300, height = 180, units = c("mm"), dpi = 100)

grp_284_peanut_K2 <- find.clusters(gl_284_raw.peanut, n.pca = 100, choose.n.clust = TRUE, stat = "BIC")
grp_284_peanut_K3 <- find.clusters(gl_284_raw.peanut, n.pca = 100, choose.n.clust = TRUE, stat = "BIC")
grp_284_peanut_K4 <- find.clusters(gl_284_raw.peanut, n.pca = 100, choose.n.clust = TRUE, stat = "BIC")
grp_284_peanut_K5 <- find.clusters(gl_284_raw.peanut, n.pca = 100, choose.n.clust = TRUE, stat = "BIC")
grp_284_peanut_K6 <- find.clusters(gl_284_raw.peanut, n.pca = 100, choose.n.clust = TRUE, stat = "BIC")

plot(grp_284_peanut_K4$Kstat)
grp_284_peanut_K4$grp

produce_post_prob_df_K <- function(data, a){
  
  pop(gl_284_raw.peanut) <- data$grp
  dapc_result <- dapc(gl_284_raw.peanut, n.pca = 3, n.da = 2)
  
  results <- as_tibble(dapc_result$posterior) %>%
    mutate(pop = pop(gl_284_raw.peanut), indNames = rownames(dapc_result$posterior), market_type = peanut_core_pop.data$market_type) %>%
    gather(key = "Assigned_Pop", value = "Posterior_membership_probability", c(1:a[1])) %>%
    arrange(indNames) %>%
    select(Original_Pop = pop, Sample = indNames, market_type, Assigned_Pop,Posterior_membership_probability) %>%
    mutate(clusters = str_c("K = ", a[1]), market_type = factor(market_type, levels = c("SP", "VA", "VR", "RN")))
  
  results
}
post_prob_df_K_2 <- produce_post_prob_df_K(grp_284_peanut_K2, 2) %>%
  mutate(Assigned_Pop_n = if_else(Assigned_Pop == 2, 1, 2))

post_prob_df_K_3 <- produce_post_prob_df_K(grp_284_peanut_K3, 3) %>%
  mutate(Assigned_Pop_n = if_else(Assigned_Pop == 3, 1, if_else(Assigned_Pop == 1, 3, 2)))

post_prob_df_K_4 <- produce_post_prob_df_K(grp_284_peanut_K4, 4) %>%
  mutate(Assigned_Pop_n = if_else(Assigned_Pop == 3, 2, if_else(Assigned_Pop == 2, 3, if_else(Assigned_Pop == 4, 4, 1))))

post_prob_df_K_f <- bind_rows(post_prob_df_K_2, post_prob_df_K_3, post_prob_df_K_4) %>%
  mutate(Assigned_Pop_n = as.character(Assigned_Pop_n))

write_delim(post_prob_df_K_f, "analysis/post_prob_df_K_f", delim = "\t", col_names = TRUE)

dapc_284_plot_scluster <- ggplot(post_prob_df_K_f, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop_n)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual(values = cols_peanut_284) +
  facet_grid(clusters~market_type, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 6)) +
  labs(y = "Posterior membership probability", x = "284 peanut accessions") +
  theme_bw() +
  theme(
    #panel.background = element_rect(colour = "black", linetype = 1),
    #strip.background = element_rect(colour = "black", linetype = 1),
    strip.text = element_text(
      colour = "black",
      face = "bold",
      size = 16
    ),
    legend.text = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 4, angle = 90),
    legend.position="none"
  ) + 
  theme(axis.title.y = element_text(color = "black"), axis.title.x = element_text(color = "black"),  plot.margin = margin(0,0.5,0,0.5, "cm")) 

ggsave("./analysis/dapc_284_plot_scluster.tiff", dapc_284_plot_scluster, width = 300, height = 180, units = c("mm"), dpi = 600)
ggsave("./analysis/dapc_284_plot_scluster.jpeg", dapc_284_plot_scluster, width = 300, height = 180, units = c("mm"), dpi = 100)