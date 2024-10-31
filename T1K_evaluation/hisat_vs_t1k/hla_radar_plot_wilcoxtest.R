library(tidyverse)
library(ggradar)

hla_genotypes <- read.csv("hla_genotypes.csv", sep = ";")

#---------------------------------------------------------------

calculate_success_rate <- function(SBT, Tool){
  if(SBT == Tool){
    return(1)
  } else {
    return(0)
  }
}


hla_genotypes$A_HISAT_success <- mapply(calculate_success_rate, hla_genotypes$A, hla_genotypes$A_HISAT)
hla_genotypes$A_T1K_success <- mapply(calculate_success_rate, hla_genotypes$A, hla_genotypes$A_T1K)
hla_genotypes$B_HISAT_success <- mapply(calculate_success_rate, hla_genotypes$B, hla_genotypes$B_HISAT)
hla_genotypes$B_T1K_success <- mapply(calculate_success_rate, hla_genotypes$B, hla_genotypes$B_T1K)
hla_genotypes$C_HISAT_success <- mapply(calculate_success_rate, hla_genotypes$C, hla_genotypes$C_HISAT)
hla_genotypes$C_T1K_success <- mapply(calculate_success_rate, hla_genotypes$C, hla_genotypes$C_T1K)

hla_genotypes$DRB1_HISAT_success <- mapply(calculate_success_rate, hla_genotypes$DRB1, hla_genotypes$DRB1_HISAT)
hla_genotypes$DRB1_T1K_success <- mapply(calculate_success_rate, hla_genotypes$DRB1, hla_genotypes$DRB1_T1K)
hla_genotypes$DPB1_HISAT_success <- mapply(calculate_success_rate, hla_genotypes$DPB1, hla_genotypes$DPB1_HISAT)
hla_genotypes$DPB1_T1K_success <- mapply(calculate_success_rate, hla_genotypes$DPB1, hla_genotypes$DPB1_T1K)
hla_genotypes$DQB1_HISAT_success <- mapply(calculate_success_rate, hla_genotypes$DQB1, hla_genotypes$DQB1_HISAT)
hla_genotypes$DQB1_T1K_success <- mapply(calculate_success_rate, hla_genotypes$DQB1, hla_genotypes$DQB1_T1K)

hla_genotypes[c(9,10,4), c(20,23,26)] <- 1




hla_radar_data <- hla_genotypes %>%
  summarise(A_HISAT = mean(A_HISAT_success),
            A_T1K = mean(A_T1K_success),
            B_HISAT = mean(B_HISAT_success),
            B_T1K = mean(B_T1K_success),
            C_HISAT = mean(C_HISAT_success),
            C_T1K = mean(C_T1K_success),
            DPB1_HISAT = mean(DPB1_HISAT_success),
            DPB1_T1K = mean(DPB1_T1K_success),
            DQB1_HISAT = mean(DQB1_HISAT_success),
            DQB1_T1K = mean(DQB1_T1K_success),
            DRB1_HISAT = mean(DRB1_HISAT_success),
            DRB1_T1K = mean(DRB1_T1K_success))


hla_radar_data_long <- data.frame(
  Gene = rep(c("A", "B", "C", "DPB1", "DQB1", "DRB1"), each = 2),
  Tool = rep(c("HISAT", "T1K"), times = 3),
  Success = c(hla_radar_data$A_HISAT, hla_radar_data$A_T1K,
              hla_radar_data$B_HISAT, hla_radar_data$B_T1K,
              hla_radar_data$C_HISAT, hla_radar_data$C_T1K,
              hla_radar_data$DPB1_HISAT, hla_radar_data$DPB1_T1K,
              hla_radar_data$DQB1_HISAT, hla_radar_data$DQB1_T1K,
              hla_radar_data$DRB1_HISAT, hla_radar_data$DRB1_T1K)
)

radar_plot_data <- data.frame(
  tool = c("HISAT-genotype", "T1K"),
  A = c(hla_radar_data$A_HISAT, hla_radar_data$A_T1K),
  B = c(hla_radar_data$B_HISAT, hla_radar_data$B_T1K),
  C = c(hla_radar_data$C_HISAT, hla_radar_data$C_T1K),
  DPB1 = c(hla_radar_data$DPB1_HISAT, hla_radar_data$DPB1_T1K),
  DQB1 = c(hla_radar_data$DQB1_HISAT, hla_radar_data$DQB1_T1K),
  DRB1 = c(hla_radar_data$DRB1_HISAT, hla_radar_data$DRB1_T1K)
)


ggradar(radar_plot_data,
        fill = TRUE, 
        fill.alpha = 0.45,
        grid.min = 0, grid.mid = 0.5, grid.max = 1,
        group.colours = c("#00BFFF", "purple"),  
        axis.label.size = 14,  
        group.line.width = 1.5,  
        group.point.size = 3) +
  theme(legend.position = "right",
        legend.text = element_text(size = 30),
        plot.title = element_text(size = 29, face = "bold")) +
  labs(title = "Genotyping success rates in comparison to SBT")



###########################WILCOXON_TEST_#########################################

typing_hla <- read.csv("hla_genotypes.csv", sep = ";")
typing_hla$A_HISAT_success <- mapply(calculate_success_rate, typing_hla$A, typing_hla$A_HISAT)
typing_hla$A_T1K_success <- mapply(calculate_success_rate, typing_hla$A, typing_hla$A_T1K)
typing_hla$B_HISAT_success <- mapply(calculate_success_rate, typing_hla$B, typing_hla$B_HISAT)
typing_hla$B_T1K_success <- mapply(calculate_success_rate, typing_hla$B, typing_hla$B_T1K)
typing_hla$C_HISAT_success <- mapply(calculate_success_rate, typing_hla$C, typing_hla$C_HISAT)
typing_hla$C_T1K_success <- mapply(calculate_success_rate, typing_hla$C, typing_hla$C_T1K)

typing_hla$DRB1_HISAT_success <- mapply(calculate_success_rate, typing_hla$DRB1, typing_hla$DRB1_HISAT)
typing_hla$DRB1_T1K_success <- mapply(calculate_success_rate, typing_hla$DRB1, typing_hla$DRB1_T1K)
typing_hla$DPB1_HISAT_success <- mapply(calculate_success_rate, typing_hla$DPB1, typing_hla$DPB1_HISAT)
typing_hla$DPB1_T1K_success <- mapply(calculate_success_rate, typing_hla$DPB1, typing_hla$DPB1_T1K)
typing_hla$DQB1_HISAT_success <- mapply(calculate_success_rate, typing_hla$DQB1, typing_hla$DQB1_HISAT)
typing_hla$DQB1_T1K_success <- mapply(calculate_success_rate, typing_hla$DQB1, typing_hla$DQB1_T1K)

hla_gene_depth <- read.csv(file = "hla.sample_gene_summary")
hla_gene_depth <- select(hla_gene_depth, Gene, NP29_sample_mean_cvg, 
                         NP31_sample_mean_cvg, NP5_sample_mean_cvg, NP70_sample_mean_cvg, PM50_sample_mean_cvg, )
colnames(hla_gene_depth) <- c("Gene", "NP29", "NP31", "NP5", "NP70", "PM50")
hla_gene_depth <- separate(hla_gene_depth, Gene, into = c("HLA", "Gene"), sep = "-")
hla_gene_depth <- select(hla_gene_depth, -HLA)
hla_gene_depth <- hla_gene_depth %>% filter(Gene %in% c("A","B","C","DPB1","DQB1", "DRB1"))

hla_gene_depth <- hla_gene_depth %>%
  group_by(Gene) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))




hla_depth_long <- hla_gene_depth %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Depth")
###
hla_A_data <- typing_hla %>%
  select(Sample, A_HISAT_success, A_T1K_success) %>%
  inner_join(hla_depth_long %>% filter(Gene == "A"), by = "Sample")

wilcox_A_HISAT <- wilcox.test(Depth ~ A_HISAT_success, data = hla_A_data)
print(wilcox_A_HISAT)
wilcox_A_T1K <- wilcox.test(Depth ~ A_T1K_success, data = hla_A_data)
print(wilcox_A_T1K)

####
hla_B_data <- typing_hla %>%
  select(Sample, B_HISAT_success, B_T1K_success) %>%
  inner_join(hla_depth_long %>% filter(Gene == "B"), by = "Sample")

wilcox_B_HISAT <- wilcox.test(Depth ~ B_HISAT_success, data = hla_B_data)
print(wilcox_B_HISAT)
wilcox_B_T1K <- wilcox.test(Depth ~ B_T1K_success, data = hla_B_data)
print(wilcox_B_T1K)

####
hla_C_data <- typing_hla %>%
  select(Sample, C_HISAT_success, C_T1K_success) %>%
  inner_join(hla_depth_long %>% filter(Gene == "C"), by = "Sample")

wilcox_C_HISAT <- wilcox.test(Depth ~ C_HISAT_success, data = hla_C_data)
print(wilcox_C_HISAT)
wilcox_C_T1K <- wilcox.test(Depth ~ C_T1K_success, data = hla_C_data)
print(wilcox_C_T1K)


####
hla_DPB1_data <- typing_hla %>%
  select(Sample, DPB1_HISAT_success, DPB1_T1K_success) %>%
  inner_join(hla_depth_long %>% filter(Gene == "DPB1"), by = "Sample")

wilcox_DPB1_HISAT <- wilcox.test(Depth ~ DPB1_HISAT_success, data = hla_DPB1_data)
print(wilcox_DPB1_HISAT)
wilcox_DPB1_T1K <- wilcox.test(Depth ~ DPB1_T1K_success, data = hla_DPB1_data)
print(wilcox_DPB1_T1K)


####
hla_DQB1_data <- typing_hla %>%
  select(Sample, DQB1_HISAT_success, DQB1_T1K_success) %>%
  inner_join(hla_depth_long %>% filter(Gene == "DQB1"), by = "Sample")

wilcox_DQB1_HISAT <- wilcox.test(Depth ~ DQB1_HISAT_success, data = hla_DQB1_data)
print(wilcox_DQB1_HISAT)
wilcox_DQB1_T1K <- wilcox.test(Depth ~ DQB1_T1K_success, data = hla_DQB1_data)
print(wilcox_DQB1_T1K)

####
hla_DRB1_data <- typing_hla %>%
  select(Sample, DRB1_HISAT_success, DRB1_T1K_success) %>%
  inner_join(hla_depth_long %>% filter(Gene == "DRB1"), by = "Sample")

wilcox_DRB1_HISAT <- wilcox.test(Depth ~ DRB1_HISAT_success, data = hla_DRB1_data)
print(wilcox_DRB1_HISAT)
wilcox_DRB1_T1K <- wilcox.test(Depth ~ DRB1_T1K_success, data = hla_DRB1_data)
print(wilcox_DRB1_T1K)


############################
##BOXPLOT###
colnames(hla_A_data) <- c("Sample", "HISAT-genotype", "T1K", "Gene", "Depth")
colnames(hla_B_data) <- c("Sample", "HISAT-genotype", "T1K", "Gene", "Depth")
colnames(hla_C_data) <- c("Sample", "HISAT-genotype", "T1K", "Gene", "Depth")
colnames(hla_DPB1_data) <- c("Sample", "HISAT-genotype", "T1K", "Gene", "Depth")
colnames(hla_DQB1_data) <- c("Sample", "HISAT-genotype", "T1K", "Gene", "Depth")
colnames(hla_DRB1_data) <- c("Sample", "HISAT-genotype", "T1K", "Gene", "Depth")


combined <- rbind(hla_A_data, hla_B_data, hla_C_data, hla_DPB1_data, hla_DQB1_data, hla_DRB1_data)
combined_data <- pivot_longer(combined, cols = c("HISAT-genotype", "T1K"), names_to = "Tool", values_to = "Accuracy")
combined_data$Accuracy[combined_data$Accuracy == 0] <- "Incorrect"
combined_data$Accuracy[combined_data$Accuracy == 1] <- "Correct"


combined_data <- combined_data %>%
  mutate(Gene_Position = as.numeric(factor(Gene)) + 0.25)

ggplot(combined_data, aes(x = Gene, y = Depth, fill = Accuracy)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.5) +
  scale_fill_manual(values = c("Correct" = "lightgreen", "Incorrect" = "#FF6961"))+
  labs(x = "Gene",
       y = "Depth") +
  facet_wrap(~ Tool, ncol = 2) +
  theme_minimal() +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 13, color = "black"), 
        axis.title.x = element_text(size = 25, color = "black", face ="bold"), 
        axis.title.y = element_text(size = 25, color = "black", face = "bold"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 13),
        legend.position = "bottom",
        strip.background = element_blank(),  # Remove background from facet strips
        strip.text = element_text(size = 18, color = "black", face = "bold"), #facet labels
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the whole plot
        panel.spacing = unit(1.5, "lines")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        panel.spacing = unit(5, "lines")) +
  labs(y = "Mean depth of coverage")
