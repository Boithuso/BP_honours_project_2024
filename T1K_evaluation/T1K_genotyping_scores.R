library(tidyverse)
#There is no difference between leaving the -s minimum aligment similarity at the default of 0.8 and changing it to 0.6
NP5_genotype <- read_tsv("NP5_genotype.tsv",col_names = FALSE)
colnames(NP5_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP5_genotype <- select(NP5_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP5_genotype$sample <- "NP5"
NP5_genotype <- pivot_longer(NP5_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP5_genotype <- filter(NP5_genotype, NP5_genotype$allele != ".")
NP5_genotype <- filter(NP5_genotype, NP5_genotype$quality > 0 | NP5_genotype$quality > 0)
NP5_genotype <- filter(NP5_genotype, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")


#--------------------------------------------------------------------

NP70_genotype <- read_tsv("NP70_genotype.tsv",col_names = FALSE)
colnames(NP70_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP70_genotype <- select(NP70_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP70_genotype$sample <- "NP70"
NP70_genotype <- pivot_longer(NP70_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP70_genotype <- filter(NP70_genotype, NP70_genotype$allele != ".")
NP70_genotype <- filter(NP70_genotype, NP70_genotype$quality > 0 | NP70_genotype$quality > 0)
NP70_genotype <- filter(NP70_genotype, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")

#------------------------------------------------------------
PM50_hla_genotype <- read_tsv("PM50_hla_genotype.tsv", col_names = FALSE)
colnames(PM50_hla_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
PM50_hla_genotype <- select(PM50_hla_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
PM50_hla_genotype$sample <- "PM50"
PM50_hla_genotype <- pivot_longer(PM50_hla_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
PM50_hla_genotype <- filter(PM50_hla_genotype, PM50_hla_genotype$allele != ".")
PM50_hla_genotype <- filter(PM50_hla_genotype, PM50_hla_genotype$quality > 0 | PM50_hla_genotype$quality > 0)
PM50_hla_genotype <- filter(PM50_hla_genotype, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")

#---------------------------------------------------------------

NP31_hla_genotype <- read_tsv("NP31_hla_genotype.tsv", col_names = FALSE)
colnames(NP31_hla_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP31_hla_genotype <- select(NP31_hla_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP31_hla_genotype$sample <- "NP31"
NP31_hla_genotype <- pivot_longer(NP31_hla_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP31_hla_genotype <- filter(NP31_hla_genotype, NP31_hla_genotype$allele != ".")
NP31_hla_genotype <- filter(NP31_hla_genotype, NP31_hla_genotype$quality > 0 | NP31_hla_genotype$quality > 0)
NP31_hla_genotype <- filter(NP31_hla_genotype, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")


#-------------------------------------------------------------
NP29_hla_genotype <- read_tsv("NP20_hla_genotype.tsv", col_names = FALSE)
colnames(NP29_hla_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP29_hla_genotype <- select(NP29_hla_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP29_hla_genotype$sample <- "NP29"
NP29_hla_genotype <- pivot_longer(NP29_hla_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP29_hla_genotype <- filter(NP29_hla_genotype, NP29_hla_genotype$allele != ".")
NP29_hla_genotype <- filter(NP29_hla_genotype, NP29_hla_genotype$quality > 0 | NP29_hla_genotype$quality > 0)
NP29_hla_genotype <- filter(NP29_hla_genotype, Gene_name != "MICA" & Gene_name != "MICB" & Gene_name != "TAP1" & Gene_name != "TAP2")



HLA_typing <- rbind(NP5_genotype, NP70_genotype, PM50_hla_genotype, NP31_hla_genotype, NP29_hla_genotype)
HLA_typing <- HLA_typing %>% group_by(Gene_name)


ggplot(HLA_typing, mapping = aes(x =sample, y = Gene_name, fill = quality)) + geom_tile() +
  scale_fill_gradient(low = "#2c7bb6", high = "maroon") + labs(x= "Sample", y= "Gene name") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), 
        axis.title.x = element_text(size = 25, color = "black"),  axis.title.y = element_text(size = 25, color = "black"))



#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

NP5_kir_genotype <- read_tsv("NP5_kir_genotype.tsv",col_names = FALSE)
colnames(NP5_kir_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP5_kir_genotype <- select(NP5_kir_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP5_kir_genotype$sample <- "NP5"
NP5_kir_genotype <- pivot_longer(NP5_kir_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP5_kir_genotype <- filter(NP5_kir_genotype, NP5_kir_genotype$allele != ".")
NP5_kir_genotype <- filter(NP5_kir_genotype, NP5_kir_genotype$quality > 0 | NP5_kir_genotype$quality > 0)




#-----------------------------------------------------------------

NP70_kir_genotype <- read_tsv("NP70_kir_genotype.tsv",col_names = FALSE)
colnames(NP70_kir_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP70_kir_genotype <- select(NP70_kir_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP70_kir_genotype$sample <- "NP70"
NP70_kir_genotype <- pivot_longer(NP70_kir_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP70_kir_genotype <- filter(NP70_kir_genotype, NP70_kir_genotype$allele != ".")
NP70_kir_genotype <- filter(NP70_kir_genotype, NP70_kir_genotype$quality > 0 | NP70_kir_genotype$quality > 0)

#--------------------------------------------------------------
PM50_kir_genotype <- read_tsv("PM50_kir_genotype.tsv",col_names = FALSE)
colnames(PM50_kir_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
PM50_kir_genotype <- select(PM50_kir_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
PM50_kir_genotype$sample <- "PM50"
PM50_kir_genotype <- pivot_longer(PM50_kir_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
PM50_kir_genotype <- filter(PM50_kir_genotype, PM50_kir_genotype$allele != ".")
PM50_kir_genotype <- filter(PM50_kir_genotype, PM50_kir_genotype$quality > 0 | PM50_kir_genotype$quality > 0)

write.csv(PM50_kir_genotype, "output.csv", row.names=FALSE, quote=FALSE) 
#------------------------------------------------------------------
NP31_kir_genotype <- read_tsv("NP31_kir_genotype.tsv",col_names = FALSE)
colnames(NP31_kir_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP31_kir_genotype <- select(NP31_kir_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP31_kir_genotype$sample <- "NP31"
NP31_kir_genotype <- pivot_longer(NP31_kir_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP31_kir_genotype <- filter(NP31_kir_genotype, NP31_kir_genotype$allele != ".")
NP31_kir_genotype <- filter(NP31_kir_genotype, NP31_kir_genotype$quality > 0 | NP31_kir_genotype$quality > 0)


NP29_kir_genotype <- read_tsv("NP29_kir_genotype.tsv",col_names = FALSE)
colnames(NP29_kir_genotype) <- c("Gene_name", "num_of_allele", "allele_1", "abundance_1", "quality_1", "allele_2", "abundance_2", "quality_2", "secondary_allele")
NP29_kir_genotype <- select(NP29_kir_genotype, -abundance_1, - abundance_2, -secondary_allele, -num_of_allele)
NP29_kir_genotype$sample <- "NP29"
NP29_kir_genotype <- pivot_longer(NP29_kir_genotype, cols = c(allele_1, quality_1, allele_2, quality_2), names_to = c(".value", "allele_number"), names_pattern = "(allele|quality)_(1|2)")
NP29_kir_genotype <- filter(NP29_kir_genotype, NP29_kir_genotype$allele != ".")
NP29_kir_genotype <- filter(NP29_kir_genotype, NP29_kir_genotype$quality > 0 | NP29_kir_genotype$quality > 0)


KIR_typing <- rbind(NP5_kir_genotype, NP70_kir_genotype, PM50_kir_genotype, NP31_kir_genotype, NP29_kir_genotype)


ggplot(KIR_typing, mapping = aes(x =sample, y = Gene_name, fill = quality)) + geom_tile() +
  scale_fill_gradient(low = "#2c7bb6", high = "maroon") + labs(x= "Sample", y= "Gene name") +
   theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 20, color = "black"), axis.text.y = element_text(size = 20, color = "black"), 
        axis.title.x = element_text(size = 25, color = "black"),  axis.title.y = element_text(size = 25, color = "black"))

