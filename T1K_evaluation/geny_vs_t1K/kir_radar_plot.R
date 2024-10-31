#install.packages("tidyverse")
#install.packages(ggradar)
library(tidyverse)
library(ggradar)


#------------------------------------------------------


radar_data <- c(40, 0.8, 1, 0, 0, 0.6, 0.6, 0.8, 1, 1, 0.8, 0.4, 0.8, 0.6, 1, 1, 1, 40, 1, 1, 0.4, 0.2, 0.8, 0.4, 0.4, 1, 1, 0.4, 0, 0.2, 0.2, 0.8, 0.4, 0.8)
radar_plot <- data.frame(matrix(data = radar_data, ncol = 17, nrow = 2, byrow = TRUE))
colnames(radar_plot) <- c("Tools", "KIR2DL4", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR2DP1", "KIR2DL1", "KIR2DL3", "KIR2DS4", "KIR3DL1", "KIR2DL2", "KIR2DL5", "KIR2DS2", "KIR2DS3", "KIR2DS1", "KIR2DS5", "KIR3DS1")
radar_plot[, 1] <- c("T1K", "Geny")

#colour of lines:
colours <- c("cyan", "purple")





#------------------------------------------------

# colour by gene function
kir_genes <- c("KIR2DL4", "KIR3DL2", "KIR3DL3", "KIR3DP1", "KIR2DP1", "KIR2DL1", "KIR2DL3", "KIR2DS4", "KIR3DL1", "KIR2DL2", "KIR2DL5", "KIR2DS2", "KIR2DS3", "KIR2DS1", "KIR2DS5", "KIR3DS1")
gene_functions <- c("framework", "framework", "framework", "pseudogene", "pseudogene", "inhibitory", "inhibitory", "activating", "inhibitory", "inhibitory", "inhibitory", "activating", "activating", "activating", "activating", "activating")
label_colours <- c("inhibitory" = "red", "activating" = "green", "framework" = "blue", "pseudogene" = "black")





#________________________________#


# Calculate positions for the custom labels
n_vars <- length(kir_genes)
angles <- seq(0, 2 * pi, length.out = n_vars + 1)[-1]

# Convert angles to x and y coordinates
label_data <- data.frame(
  x = cos(angles) * 1.3,  # Position slightly outside the radar plot (1.2)
  y = sin(angles) * 1.3,
  label = kir_genes,
  category = gene_functions
)

# Add custom colored axis labels





#_________________________________________#

kir_genes <- c("KIR3DP1", "KIR3DL3", "KIR3DL2", "KIR2DL4", "KIR3DS1", "KIR2DS5", "KIR2DS1", "KIR2DS3", "KIR2DS2", "KIR2DL5", "KIR2DL2", "KIR3DL1", "KIR2DS4", "KIR2DL3", "KIR2DL1", "KIR2DP1")
gene_functions <- c("pseudogene", "framework", "framework", "framework", "activating", "activating", "activating", "activating", "activating", "inhibitory", "inhibitory", "inhibitory", "activating", "inhibitory", "inhibitory", "pseudogene")

# Set colors for each category
label_colours <- c("inhibitory" = "red", "activating" = "green", "framework" = "blue", "pseudogene" = "black")



# Add custom colored axis labels using manual positions
# First, you need to calculate the angle positions for each label
n_vars <- length(kir_genes)
angles <- seq(0, 2 * pi, length.out = n_vars + 1)[-1]


#________________________________#


# Calculate positions for the custom labels
n_vars <- length(kir_genes)
angles <- seq(0, 2 * pi, length.out = n_vars + 1)[-1]

# Convert angles to x and y coordinates
label_data <- data.frame(
  x = cos(angles) * 1.3,  # Position slightly outside the radar plot (1.2)
  y = sin(angles) * 1.3,
  label = kir_genes,
  category = gene_functions
)




##################################
label_names <- c("inhibitory" = "Inhibitory KIR", 
                 "activating" = "Activating KIR", 
                 "framework" = "Framework KIR", 
                 "pseudogene" = " KIR Pseudogene")




################################
ggradar(radar_plot, 
        group.colours = colours,
        fill = TRUE, 
        fill.alpha = 0.6, 
        group.point.size = 4,
        axis.label.size = 0) + 
  geom_text(data = label_data, aes(x = x, y = y, label = label, color = category), size = 6) +
  scale_color_manual(values = label_colours, labels = label_names) +  # Set custom labels
  theme(legend.position = "right",
        legend.text = element_text(size = 18),
        title = element_text(size = 35, face ="bold")) +
  guides(color = guide_legend(title = "Gene Function"), 
         fill = guide_legend(title = "Tool")) +
  labs(title = "Genotyping success rates in comparison to SSP-qPCR")
