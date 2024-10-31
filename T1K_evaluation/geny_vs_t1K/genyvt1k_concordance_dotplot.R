#!/usr/bin/Rscript

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggtext)
######################################################################
##___________________NP5_________________________#########

NP5_t1k_kir <- read_tsv(file = "NP5_kir_t1k.tsv", col_names = FALSE)
colnames(NP5_t1k_kir) <- c("Gene", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP5_t1k_kir <- select(NP5_t1k_kir, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)

NP5_t1k_kir <- pivot_longer(NP5_t1k_kir, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP5_t1k_kir <- filter(NP5_t1k_kir, NP5_t1k_kir$allele != ".")
NP5_t1k_kir <- filter(NP5_t1k_kir, NP5_t1k_kir$quality > 0 | NP5_t1k_kir$quality > 0)
NP5_t1k_kir <- filter(NP5_t1k_kir, Gene != "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")
NP5_t1k_kir <- NP5_t1k_kir[, c(1:3)]
NP5_t1k_kir <- select(NP5_t1k_kir, -allele_number)
NP5_t1k_kir <- NP5_t1k_kir %>%
  mutate(allele = str_extract_all(allele, "(?<=\\*)(\\d+)(?=\\D|$)")) %>%
  rowwise() %>%
  mutate(allele = paste(allele, collapse = ",")) %>%
  ungroup()
NP5_t1k_kir$allele <- sub(",.*", "", NP5_t1k_kir$allele)
NP5_t1k_kir$Sample <- "NP5"
NP5_t1k_kir$Tool <- "T1K"
#-------------------------------------------------------------------------------------------------------------------------------
NP5_geny_kir <- read_tsv(file = "NP5_geny_kir.tsv", col_names = FALSE)
NP5_geny_kir <- select( NP5_geny_kir, -X2, -X3)
colnames(NP5_geny_kir) <- "allele"
NP5_geny_kir <- NP5_geny_kir[-c(1:8),]

NP5_geny_kir <- NP5_geny_kir %>% 
  separate(allele, into = c("Prefix", "Gene", "Allele_Abundance"), sep = " ", extra = "merge") %>% 
  separate(Allele_Abundance, into = c("allele", "Abundance"), sep = " => ") %>% select(-Prefix)
NP5_geny_kir <- head(NP5_geny_kir, -1)
NP5_geny_kir <- select(NP5_geny_kir, -Abundance)
NP5_geny_kir$Sample <- "NP5"
NP5_geny_kir$Tool <- "Geny"
NP5_geny_kir$allele <- substr(NP5_geny_kir$allele, 1, 3)

######################################################################
##___________________PM50_________________________#########
PM50_t1k <- read_tsv(file = "PM50_kir_t1k.tsv", col_names = FALSE)
colnames(PM50_t1k) <- c("Gene", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
PM50_t1k <- select(PM50_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)

PM50_t1k <- pivot_longer(PM50_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
PM50_t1k <- filter(PM50_t1k, PM50_t1k$allele != ".")
PM50_t1k <- filter(PM50_t1k, PM50_t1k$quality > 0 | PM50_t1k$quality > 0)
PM50_t1k <- filter(PM50_t1k, Gene != "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")
PM50_t1k <- PM50_t1k[, c(1:3)]
PM50_t1k <- select(PM50_t1k, -allele_number)
PM50_t1k <- PM50_t1k %>%
  mutate(allele = str_extract_all(allele, "(?<=\\*)(\\d+)(?=\\D|$)")) %>%
  rowwise() %>%
  mutate(allele = paste(allele, collapse = ",")) %>%
  ungroup()
PM50_t1k$allele <- sub(",.*", "", PM50_t1k$allele)
PM50_t1k$Sample <- "PM50"
PM50_t1k$Tool <- "T1K"
#------------------------------------------------
PM50_geny <- read_tsv(file = "PM50_geny.tsv", col_names = FALSE)
PM50_geny <- select( PM50_geny, -X2, -X3)
colnames(PM50_geny) <- "allele"
PM50_geny <- PM50_geny[-c(1:8),]

PM50_geny <- PM50_geny %>% 
  separate(allele, into = c("Prefix", "Gene", "Allele_Abundance"), sep = " ", extra = "merge") %>% 
  separate(Allele_Abundance, into = c("allele", "Abundance"), sep = " => ") %>% select(-Prefix)
PM50_geny$allele <- gsub(".*\\*(\\d{3}).*", "\\1", PM50_geny$allele)

PM50_geny <- head(PM50_geny, -1)
PM50_geny <- select(PM50_geny, -Abundance)

PM50_geny$Sample <- "PM50"
PM50_geny$Tool <- "Geny"
PM50_geny$allele <- substr(PM50_geny$allele, 1, 3)


######################################################################
##___________________NP31_________________________#########

NP31_t1k <- read_tsv(file = "NP31_kir_t1k.tsv", col_names = FALSE)
colnames(NP31_t1k) <- c("Gene", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP31_t1k <- select(NP31_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)

NP31_t1k <- pivot_longer(NP31_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP31_t1k <- filter(NP31_t1k, NP31_t1k$allele != ".")
NP31_t1k <- filter(NP31_t1k, NP31_t1k$quality > 0 | NP31_t1k$quality > 0)
NP31_t1k <- filter(NP31_t1k, Gene != "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")
NP31_t1k <- NP31_t1k[, c(1:3)]
NP31_t1k <- select(NP31_t1k, -allele_number)
NP31_t1k <- NP31_t1k %>%
  mutate(allele = str_extract_all(allele, "(?<=\\*)(\\d+)(?=\\D|$)")) %>%
  rowwise() %>%
  mutate(allele = paste(allele, collapse = ",")) %>%
  ungroup()
NP31_t1k$allele <- sub(",.*", "", NP31_t1k$allele)
NP31_t1k$Sample <- "NP31"
NP31_t1k$Tool <- "T1K"
#------------------------------------------------
NP31_geny <- read_tsv(file = "NP31_geny.tsv", col_names = FALSE)
NP31_geny <- select( NP31_geny, -X2, -X3)
colnames(NP31_geny) <- "allele"
NP31_geny <- NP31_geny[-c(1:8),]

NP31_geny <- NP31_geny %>% 
  separate(allele, into = c("Prefix", "Gene", "Allele_Abundance"), sep = " ", extra = "merge") %>% 
  separate(Allele_Abundance, into = c("allele", "Abundance"), sep = " => ") %>% select(-Prefix)
NP31_geny$allele <- gsub(".*\\*(\\d{3}).*", "\\1", NP31_geny$allele)

NP31_geny <- head(NP31_geny, -1)
NP31_geny <- select(NP31_geny, -Abundance)
NP31_geny$Sample <- "NP31"
NP31_geny$Tool <- "Geny"
NP31_geny$allele <- substr(NP31_geny$allele, 1, 3)


######################################################################
##___________________NP70_________________________#########

NP70_t1k <- read_tsv(file = "NP70_kir_genotype.tsv", col_names = FALSE)
colnames(NP70_t1k) <- c("Gene", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP70_t1k <- select(NP70_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)

NP70_t1k <- pivot_longer(NP70_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP70_t1k <- filter(NP70_t1k, NP70_t1k$allele != ".")
NP70_t1k <- filter(NP70_t1k, NP70_t1k$quality > 0 | NP70_t1k$quality > 0)
NP70_t1k <- filter(NP70_t1k, Gene != "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")
NP70_t1k <- NP70_t1k[, c(1:3)]
NP70_t1k <- select(NP70_t1k, -allele_number)
NP70_t1k <- NP70_t1k %>%
  mutate(allele = str_extract_all(allele, "(?<=\\*)(\\d+)(?=\\D|$)")) %>%
  rowwise() %>%
  mutate(allele = paste(allele, collapse = ",")) %>%
  ungroup()
NP70_t1k$allele <- sub(",.*", "", NP70_t1k$allele)
NP70_t1k$Sample <- "NP70"
NP70_t1k$Tool <- "T1K"
#------------------------------------------------
NP70_geny <- read_tsv(file = "NP70_geny_kir.tsv", col_names = FALSE)
NP70_geny <- select( NP70_geny, -X2, -X3)
colnames(NP70_geny) <- "allele"
NP70_geny <- NP70_geny[-c(1:8),]

NP70_geny <- NP70_geny %>% 
  separate(allele, into = c("Prefix", "Gene", "Allele_Abundance"), sep = " ", extra = "merge") %>% 
  separate(Allele_Abundance, into = c("allele", "Abundance"), sep = " => ") %>% select(-Prefix)
NP70_geny$allele <- gsub(".*\\*(\\d{3}).*", "\\1", NP70_geny$allele)

NP70_geny <- head(NP70_geny, -1)
NP70_geny <- select(NP70_geny, -Abundance)
NP70_geny$Sample <- "NP70"
NP70_geny$Tool <- "Geny"
NP70_geny$allele <- substr(NP70_geny$allele, 1, 3)


######################################################################
##___________________NP29_________________________#########

NP29_t1k <- read_tsv(file = "NP29_kir_t1k.tsv", col_names = FALSE)
colnames(NP29_t1k) <- c("Gene", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP29_t1k <- select(NP29_t1k, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)

NP29_t1k <- pivot_longer(NP29_t1k, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP29_t1k <- filter(NP29_t1k, NP29_t1k$allele != ".")
NP29_t1k <- filter(NP29_t1k, NP29_t1k$quality > 0 | NP29_t1k$quality > 0)
NP29_t1k <- filter(NP29_t1k, Gene != "MICA" & Gene != "MICB" & Gene != "TAP1" & Gene != "TAP2")
NP29_t1k <- NP29_t1k[, c(1:3)]
NP29_t1k <- select(NP29_t1k, -allele_number)
NP29_t1k <- NP29_t1k %>%
  mutate(allele = str_extract_all(allele, "(?<=\\*)(\\d+)(?=\\D|$)")) %>%
  rowwise() %>%
  mutate(allele = paste(allele, collapse = ",")) %>%
  ungroup()
NP29_t1k$allele <- sub(",.*", "", NP29_t1k$allele)
NP29_t1k$Sample <- "NP29"
NP29_t1k$Tool <- "T1K"
#------------------------------------------------
NP29_geny <- read_tsv(file = "NP29_geny.tsv", col_names = FALSE)
NP29_geny <- select( NP29_geny, -X2, -X3)
colnames(NP29_geny) <- "allele"
NP29_geny <- NP29_geny[-c(1:8),]

NP29_geny <- NP29_geny %>% 
  separate(allele, into = c("Prefix", "Gene", "Allele_Abundance"), sep = " ", extra = "merge") %>% 
  separate(Allele_Abundance, into = c("allele", "Abundance"), sep = " => ") %>% select(-Prefix)
NP29_geny$allele <- gsub(".*\\*(\\d{3}).*", "\\1", NP29_geny$allele)

NP29_geny <- head(NP29_geny, -1)
NP29_geny <- select(NP29_geny, -Abundance)
NP29_geny$Sample <- "NP29"
NP29_geny$Tool <- "Geny"
NP29_geny$allele <- substr(NP29_geny$allele, 1, 3)


######################################################################
##___________________COMBINED________________________#########

kir_genotyping <- bind_rows(PM50_geny,PM50_t1k,NP5_geny_kir, NP5_t1k_kir,NP29_geny, NP29_t1k, NP31_geny, NP31_t1k, NP70_geny, NP70_t1k)
kir_genotyping <- kir_genotyping %>% pivot_wider(names_from = Tool, values_from = allele, values_fn = list(allele = function(x) paste(unique(x), collapse = ",")))

kir_genotyping <- kir_genotyping %>%
  mutate(Geny_alleles = str_split(Geny, ","),
    T1K_alleles = str_split(T1K, ","))
geny_split <- kir_genotyping %>% separate_rows(Geny, sep = ",") %>%
  mutate(Geny = ifelse(Geny == "NA", NA, Geny))

t1k_split <- kir_genotyping %>% separate_rows(T1K, sep = ",") %>% 
  mutate(T1K = ifelse(T1K == "NA", NA, T1K))

kir_genotyping <- full_join(geny_split, t1k_split, by = c("Gene", "Sample"))

final_kir_genotyping <- select(kir_genotyping, Gene, Sample, Geny.x, T1K.y, Geny.y, T1K.x)

deleted_rows <- c(4,7,11,12,14,17,20,22,35,36,57,60,62,66,69,71,72,74,76,90,93)
final_kir_genotyping <- final_kir_genotyping[-deleted_rows, ]
#final_kir_genotyping <- select(final_kir_genotyping, -Geny.y, -T1K.x)

##########___________________________##########

compare_alleles <- function(geny, t1k) {
  # Split alleles and trim whitespace
  geny_alleles <- str_split(geny, ",\\s*")[[1]] %>% sort() %>% unique()
  t1k_alleles <- str_split(t1k, ",\\s*")[[1]] %>% sort() %>% unique()
  
  # Check for missing values
  if (all(is.na(geny_alleles) | geny_alleles == "")) {
    return("No call in Geny")
  } else if (all(is.na(t1k_alleles) | t1k_alleles == "")) {
    return("No call in T1K")
  }
  
  # Compare alleles
  geny_count <- length(geny_alleles)
  t1k_count <- length(t1k_alleles)
  
  # Check for agreement

  if (geny_count == t1k_count && all(geny_alleles == t1k_alleles)) {
    return("Concordance")
  }  
  if (length(intersect(geny_alleles, t1k_alleles)) > 0) {
    return("Partial Concordance")
  }
  return("Disconcordance")
}

# Apply the comparison function to your data frame
final_kir_genotyping <- final_kir_genotyping %>%
  mutate(Comparison = mapply(compare_alleles, Geny.y, T1K.x, SIMPLIFY = TRUE))

final_kir_genotyping <- select(final_kir_genotyping, -Geny.y, -T1K.x)
# Manually correct specific rows based on row indices
final_kir_genotyping <- final_kir_genotyping %>%
  mutate(Comparison = ifelse(row_number() %in% c(4,9,10,14,26,41,47,52,55,71), "Concordance", Comparison))
#NUCLEOTIDE MISMATCHES
final_kir_genotyping$Mismatch <- c(3, 3, 0, 0, 1, 0, 0, 0, 0, 0, 
                                   3, 2, 1, 0, 2, 0, 0, 0, 0, 0, 
                                   0, 1, 0, 6, 0, 0, 1, 33, 33,
                                   5, 5, 5, 5, 0, 3, 0, 0, 2, 2, 
                                   2, 0, 3, 0, 0, 0, 0, 0, 5, 5, 
                                   0, 3, 0, 0, 0, 0, 6, 4, 4, 0,
                                   0, 0, 0, 0, 0, 0, 1, 1, 1, 3, 
                                   3, 0, 3, 0, 0, 0, 0, 0)

#-----------------------------------------
final_kir_genotyping$Comparison <- factor(final_kir_genotyping$Comparison, 
                                          levels = c("Concordance", "Partial Concordance", "Disconcordance", "No call in Geny", "No call in T1K", "No call in both"))

ggplot(final_kir_genotyping, aes(x = Gene, y = Sample, color = Comparison, size = Mismatch, alpha = 0.75)) +
  geom_point() +  # Add point size as a function of mismatch
  scale_size_continuous(range = c(3, 9), name = "Mismatch Count") +  
  scale_color_manual(values = c("Concordance" = "green", "Partial Concordance" = "orange", "Disconcordance" = "red", 
                                "No call in Geny" = "cyan", "No call in T1K" = "purple")) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        panel.grid.major = element_line(colour = "grey"), 
        panel.grid.minor = element_line(colour = "grey"),
        axis.text.x = element_text(size = 20, color = "black"), 
        axis.text.y = element_text(size = 20, color = "black"), 
        axis.title.x = element_text(size = 25, color = "black", face ="bold"), 
        axis.title.y = element_text(size = 25, color = "black", face = "bold"), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        title = element_text(size = 30, face = "bold"),
        strip.background = element_blank(),  # Remove background from facet strips
        strip.text = element_text(size = 13, color = "black", face = "bold"),
        plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),  # Border around the whole plot
        panel.spacing = unit(1.5, "lines")) +
  labs(x = "Gene", y = "Sample", title = "T1K and Geny KIR genotyping concordance") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.spacing = unit(1.5, "lines"))


ggsave("compare_t1k_vs_geny_1.png", x, width = 12, height =8, bg = "white")
