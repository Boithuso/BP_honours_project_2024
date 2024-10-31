#install.packages(tidyverse)
library(tidyverse)

#below set the working directory to the location of multiQC TSV files

setwd("")
seq_counts <- read_tsv(file = "fastqc_sequence_counts_plot.tsv")

seq_counts$Samples <- c("PM50_Bx21_R2", "PM50_Bx21_R1", "NP99_R2", "NP99_R1", "NP70_Bx5_R2","NP70_Bx5_R1", "NP5_Bx5_R2", "NP5_Bx5_R1", "NP31_Bx21_R2", "NP31_Bx21_R1", "NP29_Bx5_R2", "NP29_Bx5_R1", "CD52_R2", "CD52_R1", "CD4_R2", "CD4_R1")
seq_counts$KIR_genotype <- c("Bx21","Bx21", "","","Bx5","Bx5","Bx5","Bx5","Bx21","Bx21","Bx5","Bx5", "","","","")
seq_counts <- filter(seq_counts, KIR_genotype == "Bx5" | KIR_genotype == "Bx21")

seq_counts <- mutate(seq_counts, `Total Reads` = `Unique Reads` + `Duplicate Reads`)
seq_counts <- mutate(seq_counts, `Read number` = case_when(grepl("R1", seq_counts$Samples)~ "Read 1", grepl("R2", seq_counts$Samples)~ "Read 2"))

#TOTAL READ COUNTS BOXPLOT-----------------------:
Total_reads_boxplot <- ggplot(seq_counts, mapping = aes(x = KIR_genotype, y = `Total Reads`, fill = `KIR_genotype`)) + 
  geom_boxplot() + facet_wrap(~`Read number`) + 
  scale_fill_manual(values = c("red", "darkblue")) + 
  labs(y = expression("Total Number of Reads" ~(10^7)), title = "Total Read Counts", x = "KIR genotype") +
  scale_y_continuous(labels = scales::label_number(scale = 1e-7)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 16, color = "black"), axis.title.y = element_text(size = 16, color = "black"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14))
  
ggsave(filename="Total_Reads_Boxplot.pdf", plot = Total_reads_boxplot, width = 10, height = 6.5, dpi = 1200)
ggsave(filename="Total_Reads_Boxplot.png", plot = Total_reads_boxplot, width = 13, height = 8, dpi = 1200)

seq_counts_long <- pivot_longer(seq_counts, cols = c("Unique Reads", "Duplicate Reads"), names_to = "Reads", values_to = "Counts")

Sequence_counts <- ggplot(seq_counts_long, mapping = aes(x = Samples, y = Counts, fill = Reads)) + geom_bar(stat = "identity") + 
  coord_flip() + scale_y_continuous(expand = c(0, 0), labels = scales::label_number(scale = 1e-7)) + 
  labs(y = expression("Counts" ~(10^6)), x = "Samples", title = "Sequence Counts") +
  scale_fill_manual(values = c("skyblue", "violet"), labels = c("Duplicate Reads", "Unique Reads")) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 16, color = "black"), axis.title.y = element_text(size = 16, color = "black"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14))
ggsave(filename="Unique_and_Duplicate_seq_counts.pdf", plot = Sequence_counts, width = 10, height = 6.5, dpi = 1200)


#PER BASE SEQUENCE QUALITY--------------------------------------------------------------------------------------------
per_base_seq_qual <- read_tsv(file = "fastqc_per_base_sequence_quality_plot.tsv")
per_base_seq_qual$Sample <- c("CD4_R1", "CD4_R2", "CD52_R1", "CD52_R2", "NP29_Bx5_R1", "NP29_Bx5_R2", "NP31_Bx21_R1", "NP31_Bx21_R2", "NP5_Bx5_R1", "NP5_Bx5_R2", "NP70_Bx5_R1","NP70_Bx5_R2", "NP99_R1", "NP99_R2", "PM50_Bx21_R1", "PM50_Bx21_R2")
per_base_seq_qual <- per_base_seq_qual[-(c(1,2,3,4,13,14)),]
long_per_base_qual <- pivot_longer(per_base_seq_qual, cols = -Sample, names_to = "Position", values_to = "Quality score")
long_per_base_qual <- mutate(long_per_base_qual, `Read_number` = case_when(grepl("R1", long_per_base_qual$Sample)~ "Read 1", grepl("R2", long_per_base_qual$Sample)~ "Read 2"))
long_per_base_qual <- mutate(long_per_base_qual, `KIR genotype` = case_when(grepl("Bx5", long_per_base_qual$Sample)~ "Bx5", grepl("Bx21", long_per_base_qual$Sample)~ "Bx21"))
long_per_base_qual$Phenotype <- sub("_R[12]$", "", long_per_base_qual$Sample)
long_per_base_qual$Position <- as.numeric(as.character(long_per_base_qual$Position))

labele_values <- data.frame()


Plot_per_base_seq_qual <- ggplot(long_per_base_qual, mapping = aes(x = Position, y = `Quality score`, color= `KIR genotype`, group = Phenotype)) + 
  geom_line() + facet_grid(~`Read_number`) + ylim(0, 40) +
  geom_text(data = max_values, aes(x = max_x, y = max_y, label = round(max_y, 3)), vjust = -1.5, color = "maroon", size = 3.5) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 14, color = "black"), 
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 16, color = "black"), 
        axis.title.y = element_text(size = 16, color = "black"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) + 
  labs(x = "Position (bp)", y = "Phred Quality Score", title = "Per base quality") + 
  scale_color_manual(values = c("red", "darkblue"))
  
Plot_per_base_seq_qual

ggsave(filename="Per_Base_Sequence_Quality.pdf", plot = Plot_per_base_seq_qual, width = 10, height = 6.5, dpi = 1200)


#PER SEQUENCE QUALITY---------------------------------
per_seq_qual <- read_tsv(file = "fastqc_per_sequence_quality_scores_plot.tsv")
per_seq_qual$Sample <- c("CD4_R1", "CD4_R2", "CD52_R1", "CD52_R2", "NP29_Bx5_R1", "NP29_Bx5_R2", "NP31_Bx21_R1", "NP31_Bx21_R2", "NP5_Bx5_R1", "NP5_Bx5_R2", "NP70_Bx5_R1","NP70_Bx5_R2", "NP99_R1", "NP99_R2", "PM50_Bx21_R1", "PM50_Bx21_R2")
per_seq_qual <- per_seq_qual[-(c(1,2,3,4,13,14)),]

per_seq_qual_long <- pivot_longer(per_seq_qual, cols = -Sample, names_to = "Mean Sequence Quality (Phred) Score", values_to = "Number of sequences")
per_seq_qual_long <- mutate(per_seq_qual_long, `Read number` = case_when(grepl("R1", per_seq_qual_long$Sample)~ "Read 1", grepl("R2", per_seq_qual_long$Sample)~ "Read 2"))
per_seq_qual_long <- mutate(per_seq_qual_long, `KIR genotype` = case_when(grepl("Bx5", per_seq_qual_long$Sample)~ "Bx5", grepl("Bx21", per_seq_qual_long$Sample)~ "Bx21"))
per_seq_qual_long$Phenotype <- sub("_R[12]$", "", per_seq_qual_long$Sample)
per_seq_qual_long$`Mean Sequence Quality (Phred) Score`<- as.numeric(as.character(per_seq_qual_long$`Mean Sequence Quality (Phred) Score`))

max_values <- per_seq_qual_long %>% group_by(`Read number`) %>% summarise(max_y = max(`Number of sequences`), max_x = `Mean Sequence Quality (Phred) Score`[which.max(`Number of sequences`)])

Plot_per_seq_qual <- ggplot(per_seq_qual_long, mapping = aes(x = `Mean Sequence Quality (Phred) Score`, y = `Number of sequences`, colour = `KIR genotype`, group = Phenotype)) +
  geom_line() + facet_grid(~`Read number`) + 
  scale_y_continuous(expand = c(0, 0), labels = scales::label_number(scale = 1e-6), limits = c(0, 26000000)) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 16, color = "black"), axis.title.y = element_text(size = 16, color = "black"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) + 
  scale_color_manual(values = c("red", "darkblue")) + labs(y = expression("Number of sequences" ~(10^6)), title = "Per sequence quality") 
  
Plot_per_seq_qual
ggsave(filename="Per_Sequence_Quality_Plot.pdf", plot = Plot_per_seq_qual, width = 10, height = 6.5, dpi = 1200)


#ADAPTER CONTENT--------------------------------------------------------------------------

adapter_content <- read_tsv(file = "fastqc_adapter_content_plot.tsv")
adapter_content$Sample <- c("NP70_Bx5_R1-illumina_universal_adapter", "NP70_Bx5_R2-illumina_universal_adapter", "NP99_R2-illumina_universal_adapter")
adapter_content <- adapter_content[-3,]
long_adapter_content <- pivot_longer(adapter_content, cols = -Sample, names_to = "Position (bp)", values_to = "Percentage of sequences")
long_adapter_content <- mutate(long_adapter_content, `Read number` = case_when(grepl("R1", long_adapter_content$Sample)~ "Read 1", grepl("R2", long_adapter_content$Sample)~ "Read 2"))
long_adapter_content$`Position (bp)` <- as.numeric(as.character(long_adapter_content$`Position (bp)`))

long_adapter_content <- mutate(long_adapter_content, Legend = "NP70")

max_values <- long_adapter_content %>% group_by(`Read number`) %>% summarise(max_y = max(`Percentage of sequences`), max_x = `Position (bp)`[which.max(`Percentage of sequences`)])

Adapter_content_plot <- ggplot(long_adapter_content, mapping = aes(x= `Position (bp)`, y = `Percentage of sequences`, color = Legend)) + 
  geom_line(size = 0.7) + 
  facet_grid(~`Read number`) + 
  geom_text(data = max_values, aes(x = max_x, y = max_y, label = round(max_y, 3)), vjust = -1.5, color = "maroon", size = 3.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(long_adapter_content$`Percentage of sequences`) + 0.03)) + 
  scale_x_continuous(limits = c(0, 91)) +
  scale_color_manual(values = c("NP70" = "darkblue")) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 16, color = "black"), axis.title.y = element_text(size = 16, color = "black"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)) +
  labs(y = "Percentage of sequences (%)", title = "Adapter content")

Adapter_content_plot

ggsave(filename="Adapter_content_plot.pdf", plot = Adapter_content_plot, width = 10, height = 6.5, dpi = 1200)
