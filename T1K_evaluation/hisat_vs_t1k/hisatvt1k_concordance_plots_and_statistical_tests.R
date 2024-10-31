#!/usr/bin/Rscript

library(tidyverse)
library(ggtext)
library(gridExtra)

#set the working directory as a directory with with HISAT-genotype and T1K genotyping reports as well as GATK's DepthOfCoverage x.sample_gene_summary file
setwd("")

#########NP5P##########################
###________________________________###
NP5_t1k <- read_tsv("NP5_genotype.tsv",col_names = FALSE)
colnames(NP5_t1k) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP5_t1k <- select(NP5_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP5_t1k <- pivot_longer(NP5_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP5_t1k <- filter(NP5_t1k, NP5_t1k$allele != ".")
NP5_t1k <- filter(NP5_t1k, NP5_t1k$quality > 0)
NP5_t1k <- filter(NP5_t1k, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")
NP5_t1k <- separate(NP5_t1k, allele, c("HLA", "allele"), sep = "-")

NP5_t1k <- select(NP5_t1k, -HLA, -Gene_name, -quality)

NP5_t1k <- separate(NP5_t1k, allele, c("Gene", "allele"), sep = "\\*")
NP5_t1k <- select(NP5_t1k, Gene, allele, allele_number )

NP5_t1k$Tool <- "T1K"
NP5_t1k$Sample <- "NP5"
#--------------------------------#
NP5_hisat <- read_tsv("assembly_graph-hla.NP5_R1_fq-hla-extracted-2_fq (2).tsv", col_names = FALSE)
NP5_hisat <- filter(NP5_hisat, grepl("[1-2] ranked", NP5_hisat$X1))
NP5_hisat <- select(NP5_hisat, X1)
NP5_hisat <- separate(NP5_hisat, X1, c("allele_number", "ranked", "allele", "abundance"), sep = " ")
NP5_hisat <- select(NP5_hisat, allele, allele_number)
NP5_hisat <- separate(NP5_hisat, allele, c("Gene", "allele"), sep = "\\*")
NP5_hisat <- filter(NP5_hisat, Gene!= "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")
NP5_hisat$Tool <- "hisat"
NP5_hisat$Sample <- "NP5"

combine_NP5 <- merge(NP5_hisat, NP5_t1k, by = c("Gene", "allele_number"), suffixes = c("_hisat", "_t1k"), all = TRUE)

#########NP31##########################
###________________________________###

NP31_t1k <- read_tsv("NP31_hla_genotype.tsv",col_names = FALSE)
colnames(NP31_t1k) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP31_t1k <- select(NP31_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP31_t1k <- pivot_longer(NP31_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP31_t1k <- filter(NP31_t1k, NP31_t1k$allele != ".")
NP31_t1k <- filter(NP31_t1k, NP31_t1k$quality > 0)
NP31_t1k <- filter(NP31_t1k, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")
NP31_t1k <- separate(NP31_t1k, allele, c("HLA", "allele"), sep = "-")

NP31_t1k <- select(NP31_t1k, -HLA, -Gene_name, -quality)

NP31_t1k <- separate(NP31_t1k, allele, c("Gene", "allele"), sep = "\\*")
NP31_t1k <- select(NP31_t1k, Gene, allele, allele_number )

NP31_t1k$Tool <- "T1K"
NP31_t1k$Sample <- "NP31"

#--------------------------------#
NP31_hisat <- read_tsv("assembly_graph-hla.NP31_filtered_R1_fq-hla-extracted-1_fq.tsv", col_names = FALSE)
NP31_hisat <- filter(NP31_hisat, grepl("[1-2] ranked", NP31_hisat$X1))
NP31_hisat <- select(NP31_hisat, X1)
NP31_hisat <- separate(NP31_hisat, X1, c("allele_number", "ranked", "allele", "abundance"), sep = " ")
NP31_hisat <- select(NP31_hisat, allele, allele_number)
NP31_hisat <- separate(NP31_hisat, allele, c("Gene", "allele"), sep = "\\*")
NP31_hisat <- filter(NP31_hisat, Gene!= "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")

NP31_hisat$Tool <- "hisat"
NP31_hisat$Sample <- "NP31"


###Combine###
combine_NP31 <- merge(NP31_hisat, NP31_t1k, by = c("Gene", "allele_number"), suffixes = c("_hisat", "_t1k"), all = TRUE)

#########NP70##########################
###________________________________###

NP70_t1k <- read_tsv("NP70_hla_genotype.tsv",col_names = FALSE)
colnames(NP70_t1k) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP70_t1k <- select(NP70_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP70_t1k <- pivot_longer(NP70_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP70_t1k <- filter(NP70_t1k, NP70_t1k$allele != ".")
NP70_t1k <- filter(NP70_t1k, NP70_t1k$quality > 0)
NP70_t1k <- filter(NP70_t1k, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")
NP70_t1k <- separate(NP70_t1k, allele, c("HLA", "allele"), sep = "-")

NP70_t1k <- select(NP70_t1k, -HLA, -Gene_name, -quality)

NP70_t1k <- separate(NP70_t1k, allele, c("Gene", "allele"), sep = "\\*")
NP70_t1k <- select(NP70_t1k, Gene, allele, allele_number )

NP70_t1k$Tool <- "T1K"
NP70_t1k$Sample <- "NP70"

#--------------------------------#
NP70_hisat <- read_tsv("assembly_graph-hla.NP70_6484_GTGAAA_m4_n10_R1_fq-hla-extracted-2_fq.tsv", col_names = FALSE)
NP70_hisat <- filter(NP70_hisat, grepl("[1-2] ranked", NP70_hisat$X5))
NP70_hisat <- select(NP70_hisat, X5)
NP70_hisat <- separate(NP70_hisat, X5, c("allele_number", "ranked", "allele", "abundance"), sep = " ")
NP70_hisat <- select(NP70_hisat, allele, allele_number)
NP70_hisat <- separate(NP70_hisat, allele, c("Gene", "allele"), sep = "\\*")
NP70_hisat <- filter(NP70_hisat, Gene!= "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")

NP70_hisat$Tool <- "hisat"
NP70_hisat$Sample <- "NP70"


###Combine###
combine_NP70 <- merge(NP70_hisat, NP70_t1k, by = c("Gene", "allele_number"), suffixes = c("_hisat", "_t1k"), all = TRUE)

#########PM50##########################
###________________________________###

PM50_t1k <- read_tsv("PM50_hla_genotype.tsv",col_names = FALSE)
colnames(PM50_t1k) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
PM50_t1k <- select(PM50_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
PM50_t1k <- pivot_longer(PM50_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
PM50_t1k <- filter(PM50_t1k, PM50_t1k$allele != ".")
PM50_t1k <- filter(PM50_t1k, PM50_t1k$quality > 0)
PM50_t1k <- filter(PM50_t1k, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")
PM50_t1k <- separate(PM50_t1k, allele, c("HLA", "allele"), sep = "-")

PM50_t1k <- select(PM50_t1k, -HLA, -Gene_name, -quality)

PM50_t1k <- separate(PM50_t1k, allele, c("Gene", "allele"), sep = "\\*")
PM50_t1k <- select(PM50_t1k, Gene, allele, allele_number )

PM50_t1k$Tool <- "T1K"
PM50_t1k$Sample <- "PM50"

#--------------------------------#
PM50_hisat <- read_tsv("assembly_graph-hla.PM50_filtered_R1_fq-hla-extracted-1_fq.tsv", col_names = FALSE)
PM50_hisat <- filter(PM50_hisat, grepl("[1-2] ranked", PM50_hisat$X5))
PM50_hisat <- select(PM50_hisat, X5)
PM50_hisat <- separate(PM50_hisat, X5, c("allele_number", "ranked", "allele", "abundance"), sep = " ")
PM50_hisat <- select(PM50_hisat, allele, allele_number)
PM50_hisat <- separate(PM50_hisat, allele, c("Gene", "allele"), sep = "\\*")
PM50_hisat <- filter(PM50_hisat, Gene!= "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")

PM50_hisat$Tool <- "hisat"
PM50_hisat$Sample <- "PM50"

###Combine###
combine_PM50 <- merge(PM50_hisat, PM50_t1k, by = c("Gene", "allele_number"), suffixes = c("_hisat", "_t1k"), all = TRUE)


#########PM50##########################
###________________________________###

NP29_t1k <- read_tsv("NP20_hla_genotype.tsv",col_names = FALSE)
colnames(NP29_t1k) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP29_t1k <- select(NP29_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP29_t1k <- pivot_longer(NP29_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP29_t1k <- filter(NP29_t1k, NP29_t1k$allele != ".")
NP29_t1k <- filter(NP29_t1k, NP29_t1k$quality > 0)
NP29_t1k <- filter(NP29_t1k, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")
NP29_t1k <- separate(NP29_t1k, allele, c("HLA", "allele"), sep = "-")

NP29_t1k <- select(NP29_t1k, -HLA, -Gene_name, -quality)

NP29_t1k <- separate(NP29_t1k, allele, c("Gene", "allele"), sep = "\\*")
NP29_t1k <- select(NP29_t1k, Gene, allele, allele_number )

NP29_t1k$Tool <- "T1K"
NP29_t1k$Sample <- "NP29"

#--------------------------------#
NP29_hisat <- read_tsv("assembly_graph-hla.NP29_filtered_R1_fq-hla-extracted-1_fq.tsv", col_names = FALSE)
NP29_hisat <- filter(NP29_hisat, grepl("[1-2] ranked", NP29_hisat$X5))
NP29_hisat <- select(NP29_hisat, X5)
NP29_hisat <- separate(NP29_hisat, X5, c("allele_number", "ranked", "allele", "abundance"), sep = " ")
NP29_hisat <- select(NP29_hisat, allele, allele_number)
NP29_hisat <- separate(NP29_hisat, allele, c("Gene", "allele"), sep = "\\*")
NP29_hisat <- filter(NP29_hisat, Gene!= "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")

NP29_hisat$Tool <- "hisat"
NP29_hisat$Sample <- "NP29"

###Combine###
combine_NP29 <- merge(NP29_hisat, NP29_t1k, by = c("Gene", "allele_number"), suffixes = c("_hisat", "_t1k"), all = TRUE)


######################################################################
##___________________COMBINED________________________#########

hla_genotyping <- bind_rows(PM50_hisat,PM50_t1k,NP5_hisat, NP5_t1k, NP31_hisat, NP31_t1k, NP70_hisat, NP70_t1k, NP29_hisat,NP29_t1k)
hla_genotyping <- hla_genotyping %>% pivot_wider(names_from = Tool, values_from = allele, values_fn = list(allele = function(x) paste(unique(x), collapse = ",")))



hla_genotyping <- separate(hla_genotyping, hisat, c("hisat_2_digit", "hisat_4_digit", "hisat_6_digit"), sep = ":")
hla_genotyping <- separate(hla_genotyping, T1K, c("t1k_2_digit", "t1k_4_digit", "t1k_6_digit"), sep = ":")

####________________________######
#2digit:
temp <- hla_genotyping$t1k_2_digit[5]
hla_genotyping$t1k_2_digit[5] <- hla_genotyping$t1k_2_digit[6]
hla_genotyping$t1k_2_digit[6] <- temp

temp <- hla_genotyping$t1k_2_digit[66]
hla_genotyping$t1k_2_digit[66] <- hla_genotyping$t1k_2_digit[67]
hla_genotyping$t1k_2_digit[67] <- temp

temp <- hla_genotyping$t1k_2_digit[106]
hla_genotyping$t1k_2_digit[106] <- hla_genotyping$t1k_2_digit[107]
hla_genotyping$t1k_2_digit[107] <- temp

temp <- hla_genotyping$t1k_2_digit[140]
hla_genotyping$t1k_2_digit[140] <- hla_genotyping$t1k_2_digit[141]
hla_genotyping$t1k_2_digit[141] <- temp

temp <- hla_genotyping$t1k_2_digit[164]
hla_genotyping$t1k_2_digit[164] <- hla_genotyping$t1k_2_digit[165]
hla_genotyping$t1k_2_digit[165] <- temp

temp <- hla_genotyping$t1k_2_digit[208]
hla_genotyping$t1k_2_digit[208] <- hla_genotyping$t1k_2_digit[209]
hla_genotyping$t1k_2_digit[209] <- temp

#4 digit:

temp <- hla_genotyping$t1k_4_digit[5]
hla_genotyping$t1k_4_digit[5] <- hla_genotyping$t1k_4_digit[6]
hla_genotyping$t1k_4_digit[6] <- temp

temp <- hla_genotyping$t1k_4_digit[106]
hla_genotyping$t1k_4_digit[106] <- hla_genotyping$t1k_4_digit[107]
hla_genotyping$t1k_4_digit[107] <- temp

temp <- hla_genotyping$t1k_4_digit[140]
hla_genotyping$t1k_4_digit[140] <- hla_genotyping$t1k_4_digit[141]
hla_genotyping$t1k_4_digit[141] <- temp

temp <- hla_genotyping$t1k_4_digit[164]
hla_genotyping$t1k_4_digit[164] <- hla_genotyping$t1k_4_digit[165]
hla_genotyping$t1k_4_digit[165] <- temp

temp <- hla_genotyping$t1k_4_digit[208]
hla_genotyping$t1k_4_digit[208] <- hla_genotyping$t1k_4_digit[209]
hla_genotyping$t1k_4_digit[209] <- temp
temp <- hla_genotyping$t1k_4_digit[119]
hla_genotyping$t1k_4_digit[119] <- hla_genotyping$t1k_4_digit[120]
hla_genotyping$t1k_4_digit[120] <- temp

#6digit
temp <- hla_genotyping$t1k_6_digit[11]
hla_genotyping$t1k_6_digit[11] <- hla_genotyping$t1k_6_digit[12]
hla_genotyping$t1k_6_digit[12] <- temp

temp <- hla_genotyping$t1k_6_digit[106]
hla_genotyping$t1k_6_digit[106] <- hla_genotyping$t1k_6_digit[107]
hla_genotyping$t1k_6_digit[107] <- temp

#######################################
hla_genotyping <- hla_genotyping %>% mutate(concordance_2digit = case_when( hisat_2_digit == t1k_2_digit & allele_number == 1 ~ "allele 1 match",
                                                                            hisat_2_digit == t1k_2_digit & allele_number == 2 ~ "allele 2 match",
                                                                            hisat_2_digit != t1k_2_digit & allele_number == 1 ~ "allele 1 mismatch",
                                                                            hisat_2_digit != t1k_2_digit & allele_number == 2 ~ "allele 2 mismatch", 
                                                                            is.na(t1k_2_digit) & !is.na(hisat_2_digit) ~ "No call in T1K", 
                                                                            is.na(hisat_2_digit) & !is.na(t1k_2_digit) ~ "No call in hisatgenotype",
                                                                            is.na(hisat_2_digit) & is.na(t1k_2_digit) ~ "No call in both",
                                                                            TRUE ~ "Both alleles mismatch"),
                                            concordance_4digit = case_when(hisat_4_digit == t1k_4_digit & allele_number == 1 ~ "allele 1 match",
                                                                           hisat_4_digit == t1k_4_digit & allele_number == 2 ~ "allele 2 match",
                                                                           hisat_4_digit != t1k_4_digit & allele_number == 1 ~ "allele 1 mismatch",
                                                                           hisat_4_digit != t1k_4_digit & allele_number == 2 ~ "allele 2 mismatch",
                                                                           is.na(t1k_4_digit) & !is.na(hisat_4_digit) ~ "No call in T1K",
                                                                           is.na(hisat_4_digit) & !is.na(t1k_4_digit) ~ "No call in hisatgenotype",
                                                                           is.na(hisat_4_digit) & is.na(t1k_4_digit) ~ "No call in both",
                                                                           TRUE ~ "Both alleles mismatch"),
                                            concordance_6digit = case_when(hisat_6_digit == t1k_6_digit & allele_number == 1 ~ "allele 1 match",
                                                                           hisat_6_digit == t1k_6_digit & allele_number == 2 ~ "allele 2 match",
                                                                           hisat_6_digit != t1k_6_digit & allele_number == 1 ~ "allele 1 mismatch",
                                                                           hisat_6_digit != t1k_6_digit & allele_number == 2 ~ "allele 2 mismatch",
                                                                           is.na(t1k_6_digit) & !is.na(hisat_6_digit) ~ "No call in T1K",
                                                                           is.na(hisat_6_digit) & !is.na(t1k_6_digit) ~ "No call in hisatgenotype",
                                                                           is.na(hisat_6_digit) & is.na(t1k_6_digit) ~ "No call in both",
                                                                           TRUE ~ "Both alleles mismatch"))



###---------------------###
###_______PLOTTING______###
###---------------------###


# Reshape data for the heatmap
df_long_heatmap <- hla_genotyping %>%
  select(Sample, Gene, concordance_2digit, concordance_4digit, concordance_6digit) %>%
  pivot_longer(cols = starts_with("concordance"), names_to = "Resolution", values_to = "Concordance")

df_long_heatmap$Gene <- factor(df_long_heatmap$Gene, 
                               levels = c("A", "B", "C", "DRB1", "DQB1", "DPB1",
                                          "DRB3", "DRB4", "DRB5","DQA1", "DPA1", "DPB2",
                                          "DPA2","DMA", "DMB", "DOA", "DOB", "DRA", 
                                          "E", "F", "G", "H", "HFE", "K",
                                          "L", "R", "S", "T", "U", "V", "Y"))  



df_long_heatmap <- df_long_heatmap %>%
  group_by(Sample, Gene) %>%
  mutate(Concordance_numeric = case_when(
    # If there's one match and one mismatch within the same sample and gene
    any(Concordance == "allele 1 match" & Concordance == "allele 2 mismatch") |
      any(Concordance == "allele 2 match" & Concordance == "allele 1 mismatch") ~ 2, 
    Concordance == "allele 1 match" | Concordance == "allele 2 match" ~ 3,  
    Concordance == "allele 1 mismatch" | Concordance == "allele 2 mismatch" ~ 0, 
    Concordance == "No call in hisatgenotype" ~ 1,
    Concordance == "No call in T1K" ~ 1,
    Concordance == "No call in both" ~ -1,
    TRUE ~ NA_real_
  )) %>%
  ungroup()
df_long_heatmap$Concordance[c(53,54,611,90)] <- "allele 2 mismatch"
df_long_heatmap$Concordance[c(628,629,630)] <- "No call in hisatgenotype"
df_long_heatmap$Concordance[c(49,52,50,53,51,54,460,463,461,464,607,610,608,611,602,605,33,36,449,452,86,89,87,90,227,230,356,359,363,366)] <- "match and mismatch"
# Heatmap plot

heatmap_plot <- ggplot(df_long_heatmap, aes(x = Resolution, y = Sample, fill = Concordance, alpha = 0.6)) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("allele 1 match" = "lightgreen", "allele 2 match" = "lightgreen", 
                               "allele 1 mismatch" = "lightcoral", "allele 2 mismatch" = "lightcoral",
                               "No call in hisatgenotype" = "lightblue", "No call in T1K" = "#FFEB3B",  "Concordant and T1K no call" = "#dfec74","No call in both" = "white", "match and mismatch" = "#C0B788"),
                    breaks = c(
                      "allele 1 match", 
                      "allele 1 mismatch", 
                      "match and mismatch", 
                      "No call in T1K",
                      "Concordant and T1K no call",
                      "No call in hisatgenotype",
                      "No call in both"
                    ),
                    labels = c(
                      "Concordance", 
                      "Discordance", 
                      "Concordant and Discordant", 
                      "No call in T1K",
                      "Concordant and T1K no call",
                      "No call in HISAT-genotype", 
                      "Concordant no calls"
                    )) +
  geom_blank(aes(fill = "Concordant and T1K no call")) +
  scale_x_discrete(labels = c(
    "concordance_2digit" = "2-Digit",
    "concordance_4digit" = "4-Digit",
    "concordance_6digit" = "6-Digit")) +
  labs(x = "HLA Digit Resolution",
       y = "Sample") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        panel.grid.major = element_line(colour = "grey"), 
        panel.grid.minor = element_line(colour = "grey"),
        axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 13, color = "black"), 
        axis.title.x = element_text(size = 25, color = "black", face ="bold"), 
        axis.title.y = element_text(size = 25, color = "black", face = "bold"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 13),
        strip.background = element_blank(),  # Remove background from facet strips
        strip.text = element_text(size = 13, color = "black", face = "bold"),
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the whole plot
        panel.spacing = unit(1.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Gene, ncol = 6)

print(heatmap_plot)
########################
heatmap_plot <- ggplot(df_long_heatmap, aes(x = Resolution, y = Sample, fill = Concordance_numeric, alpha = 0.6)) +
  geom_tile(color = "black") + 
  scale_fill_gradientn(
    colours = c("white", "lightcoral", "lightblue", "lightgreen", "#228B22"), 
    values = scales::rescale(c(-1, 0, 1, 2, 3)),  
    name = "Genotyping Concordance",
    breaks = c(0, 1, 2, 3),
    labels = c("Discordant","No call (HISAT/T1K)", "Concordant and Disconcordant", "Concordant")  
  ) + scale_alpha_continuous(range = c(0.5, 1), guide = "none") +
  scale_x_discrete(labels = c(
    "concordance_2digit" = "2-Digit",
    "concordance_4digit" = "4-Digit",
    "concordance_6digit" = "6-Digit")) +
  labs(x = "HLA Digit Resolution",
       y = "Sample") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        panel.grid.major = element_line(colour = "grey"), 
        panel.grid.minor = element_line(colour = "grey"),
        axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 8, color = "black"), 
        axis.title.x = element_text(size = 25, color = "black", face ="bold"), 
        axis.title.y = element_text(size = 25, color = "black", face = "bold"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 13),
        strip.background = element_blank(),  # Removes background from facet strips
        strip.text = element_text(size = 13, color = "black", face = "bold"),
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the whole plot
        panel.spacing = unit(1.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Gene, ncol = 6)

print(heatmap_plot)
################presentation_Plot##

subset1 <- df_long_heatmap %>% filter(Gene %in% c("A", "B", "C", "DPB1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5","DQA1","DPA1", "DPB2"))


subset1$Concordance_numeric <- dplyr::case_when(
  subset1$Concordance == "match and mismatch" ~ 2,
  subset1$Concordance == "allele 1 match" ~ 3,
  subset1$Concordance == "allele 2 match" ~ 3,
  subset1$Concordance == "allele 1 mismatch" ~ 0,
  subset1$Concordance == "allele 2 mismatch" ~ 0,
  subset1$Concordance == "No call in hisatgenotype" ~ 1,
  subset1$Concordance == "No call in T1K" ~ 1,
  subset1$Concordance == "No call in both" ~ -1,  
  TRUE ~ NA_real_  
)


# Presentation Heatmap plot
presentation_plot <- ggplot(subset1, aes(x = Resolution, y = Sample, fill = Concordance_numeric, alpha = 0.6)) +
  geom_tile(color = "black") + scale_fill_gradientn(
    colours = c("white", "red", "lightblue", "magenta", "purple"), 
    values = scales::rescale(c(-1, 0, 1, 2, 3)),  
    name = "Genotyping Concordance",
    breaks = c(0, 1, 2, 3),
    labels = c("Discordant","No call (HISAT/T1K)", "Concordant and Disconcordant", "Concordant")  
  ) + scale_alpha_continuous(range = c(0.5, 1), guide = "none") +
  scale_x_discrete(labels = c(
    "concordance_2digit" = "2-Digit",
    "concordance_4digit" = "4-Digit",
    "concordance_6digit" = "6-Digit")) +
  labs(x = "HLA Digit Resolution",
       y = "Sample",
       title = "T1K and HISAT-genotype HLA genotyping concordance") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        panel.grid.major = element_line(colour = "grey"), 
        panel.grid.minor = element_line(colour = "grey"),
        axis.text.x = element_text(size = 20, color = "black"), 
        axis.text.y = element_text(size = 20, face = "bold", color = "black"), 
        axis.title.x = element_text(size = 30, color = "black", face ="bold"), 
        axis.title.y = element_text(size = 30, color = "black", face = "bold"), 
        plot.title = element_text(size = 29, face = "bold"), 
        legend.title = element_text(size = 25, face = "bold", margin = margin(b = 20)), 
        legend.text = element_text(size = 20),
        strip.background = element_blank(),  # Remove background from facet strips
        strip.text = element_text(size = 13, color = "black", face = "bold"),
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the whole plot
        panel.spacing = unit(1.5, "lines")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ Gene, ncol = 3)

print(presentation_plot)


################################################
##################################################
######################################################
hla_gene_depth <- read.csv(file = "hla.sample_gene_summary")
hla_gene_depth <- select(hla_gene_depth, Gene, NP29_sample_mean_cvg, 
                         NP31_sample_mean_cvg, NP5_sample_mean_cvg, NP70_sample_mean_cvg, PM50_sample_mean_cvg, )
colnames(hla_gene_depth) <- c("Gene", "NP29", "NP31", "NP5", "NP70", "PM50")
hla_gene_depth <- separate(hla_gene_depth, Gene, into = c("HLA", "Gene"), sep = "-")
hla_gene_depth <- select(hla_gene_depth, -HLA)

###############################################################

# Step 1: Reshape gene depth data to long format
hla_gene_depth_long <- hla_gene_depth %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Depth")

# Step 2: Merge genotyping data with gene depth data
hla_genotyping_depth <- hla_genotyping %>%
  left_join(hla_gene_depth_long, by = c("Gene", "Sample"))

# Step 3: Create Concordance_Status column for each resolution
hla_genotyping_depth <- hla_genotyping_depth %>%
  mutate(Concordance_Status_2digit = ifelse(concordance_2digit %in% c("allele 1 match", "allele 2 match"), "Concordant", "Discordant"),
         Concordance_Status_4digit = ifelse(concordance_4digit %in% c("allele 1 match", "allele 2 match"), "Concordant", "Discordant"),
         Concordance_Status_6digit = ifelse(concordance_6digit %in% c("allele 1 match", "allele 2 match"), "Concordant", "Discordant"))

# Step 4: Perform Wilcoxon test for each resolution
wilcox_2digit <- wilcox.test(Depth ~ Concordance_Status_2digit, data = hla_genotyping_depth)
wilcox_4digit <- wilcox.test(Depth ~ Concordance_Status_4digit, data = hla_genotyping_depth)
wilcox_6digit <- wilcox.test(Depth ~ Concordance_Status_6digit, data = hla_genotyping_depth)

# Print results
wilcox_2digit
wilcox_4digit
wilcox_6digit


hla_genotyping_depth_long <- hla_genotyping_depth %>%
  pivot_longer(cols = starts_with("Concordance_Status"), 
               names_to = "Resolution", 
               values_to = "Concordance_Status") %>%
  mutate(Resolution = case_when(
    Resolution == "Concordance_Status_2digit" ~ "2-digit",
    Resolution == "Concordance_Status_4digit" ~ "4-digit",
    Resolution == "Concordance_Status_6digit" ~ "6-digit"
  ))


# Plot the data using ggplot
ggplot(hla_genotyping_depth_long, aes(x = Concordance_Status, y = Depth, fill = Concordance_Status)) +
  geom_boxplot() +
  facet_wrap(~ Resolution) +  
  theme_minimal() +
  labs(x = NULL, y = "per gene depth of coverage") +
  scale_fill_manual(values = c("Concordant" = "lightgreen", "Discordant" = "lightcoral")) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        panel.grid.major = element_line(colour = "grey"), 
        panel.grid.minor = element_line(colour = "grey"),
        axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 13, color = "black"), 
        axis.title.x = element_text(size = 25, color = "black", face ="bold"), 
        axis.title.y = element_text(size = 25, color = "black", face = "bold"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_blank(), 
        legend.text = element_text(size = 13),
        strip.background = element_blank(),  # Remove background from facet strips
        strip.text = element_text(size = 13, color = "black", face = "bold"),
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the whole plot
        panel.spacing = unit(1.5, "lines"))





############################################################################################################
################################################Wilcoxon test concordant/discordant at two-, four-, six-digit resolution
################################################
hla_gene_depth <- hla_gene_depth %>%
  group_by(Gene) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

# Define classical class 1 and class 2 genes
class_1_genes <- c("A", "B", "C")
class_2_genes <- c("DRB1", "DPB1", "DQB1")


hla_gene_depth_long <- hla_gene_depth_long %>%
  mutate(Gene_Category = case_when(
    Gene %in% class_1_genes ~ "Classical Class 1",
    Gene %in% class_2_genes ~ "Classical Class 2",
    Gene %in% c("HLA-F", "HLA-G") ~ "Other Classical", 
    TRUE ~ "Non-Classical"  
  ))


# Wilcoxon test between classical class 1 and class 2
wilcox_class1_class2 <- wilcox.test(Depth ~ Gene_Category, 
                                    data = hla_gene_depth_long %>%
                                      filter(Gene_Category %in% c("Classical Class 1", "Classical Class 2")))

# Wilcoxon test between all classical genes and non-classical genes
wilcox_classical_non_classical <- wilcox.test(Depth ~ Gene_Category2,
                                              data = hla_gene_depth_long %>%
                                                filter(Gene_Category2 %in% c("Classical gene", "Non-Classical")))

# Print results
print(wilcox_class1_class2)
print(wilcox_classical_non_classical)




#################################################
##################gene classification boxplot########################
################################################


gene_order <- c("A", "B", "C", "DRB1", "DQB1", "DPB1", "DRB3", "DRB4", "DRB5", "DQA1", "DPA1", "DPQ2", "DPA2", "DMA", "DMB",
                "DOA", "DOB", "DRA", "E","F","G","H","HFE","K","L","R","S","T","U","V","Y")
hla_gene_depth_long <- hla_gene_depth_long %>% filter(Gene %in% c("A", "B", "C", "DRB1", "DQB1", "DPB1"))
hla_gene_depth_long <- hla_gene_depth_long %>% mutate(Gene = factor(Gene, levels=gene_order))


ggplot(hla_gene_depth_long, aes(x = Gene, y = Depth, fill = Gene_Category)) +
  geom_boxplot() +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 20, color = "black"), 
        axis.text.y = element_text(size = 20, color = "black"), 
        axis.title.x = element_text(size = 30, color = "black", face ="bold"), 
        axis.title.y = element_text(size = 30, color = "black", face = "bold"), 
        plot.title = element_text(size = 30, face = "bold"), legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        strip.background = element_blank(),  # Remove background from facet strips
        strip.text = element_text(size = 20, color = "black", face = "bold"),
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the whole plot
        panel.spacing = unit(1.5, "lines")) + 
  scale_fill_manual(values = c("Classical Class 1" = "pink", "Classical Class 2" = "lightblue", 
                               "Non-Classical" = "lightgreen"), name = "HLA gene classification") +
  labs(y = "Mean depth of coverage", x = "HLA gene",
       title = "HLA gene depth of coverage")




