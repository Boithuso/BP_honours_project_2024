#install.packages("tidyverse")
library(tidyverse)
#set the working directory as the directory with the pidstat memory consumption files
setwd("")


memory_consumption <- function(file_path, sample_name, time_increment) {
  #cleaning up the pidstat command output
  data <- read_tsv(file = file_path, colnames = FALSE)
  data <- filter(data, !startsWith(X1, "Linux"))
  data <- data[seq(4, nrow(data), by = 4), ]
  data <- separate(data, col = "X1", into c("Time", "UID", "PID", "minorfaults",
                                            "majorfaults", "VSZ", "RSS", "%MEM", "Command"), sep = "\\s+")
  data <- select(data, Time, RSS)
  
  #setting the time and fixing units:
  data$Time <- seq(from = time_increment, by = time_increment, length.out = nrow(data))
  data$RSS <- as.numeric(data$RSS) / 1000
  data$Time <- data$Time / 3600
  
  #make an origin at 0,0
  new_row <- data.frame(Time = 0, RSS = 0)
  data <- rbind(new_row, data)
  
  #calculate the maximum RAM
  max_rss <- data %>% summarise(max_y = max(RSS), max_x = Time[which.max(RSS)])
  data$sample <- sample_name
  
  return(list(data = data, max_rss = max_rss))
  
  
}


combined <- rbind(NP5, NP29, NP31, PM50, NP70)


max_rss_values <- combined %>% group_by(sample) %>% 
  summarise(max_y = max(RSS), max_x = Time[which.max(RSS)]) %>% 
  mutate(max_x = max_x + (1.1 * row_number()))

mean_time <- mean(max_time)
mean_rss <- combined %>%
  group_by(Time) %>%
  summarise(mean_RSS = mean(RSS, na.rm = TRUE)) 
mean_rss <- filter(mean_rss, Time <=  mean_time)

overall_mean_time <- mean(combined$Time)
mean_rss_line <- data.frame(Time = overall_mean_time, mean_RSS = mean(mean_rss$mean_RSS, na.rm = TRUE))

y_limit <- max(max_rss_values$max_y) + 0.5
max_time <- max(max_time)
#########################################################################
ggplot(combined, mapping = aes(x= Time, y = RSS, colour = sample)) + 
  geom_line(linewidth = 0.8) + 
  geom_line(data = mean_rss, aes(x = Time, y = mean_RSS, colour = "Mean"), color = "black", linetype = "dashed", linewidth = 0.75) +
  scale_x_continuous(breaks = seq(0, max_time, by = 1), limits = c(0, max_time)) +
  scale_y_continuous(limits = c(0, y_limit)) +
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "grey20"), 
        axis.text.x = element_text(size = 25, color = "black"), 
        axis.text.y = element_text(size = 25, color = "black"), 
        axis.title.x = element_text(size = 30, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 30, color = "black", face = "bold"), 
        plot.title = element_text(size = 16.5, face = "bold"), legend.title = element_text(size = 30), 
        legend.text = element_text(size = 25),
        panel.grid.major = element_line(colour = "grey"), 
        panel.grid.minor = element_line(colour = "grey")) + 
  labs(x = "Elapsed Time (hours)", y = "Memory usage (MB)") +
  geom_label(data = max_rss_values, 
             aes(x = max_x, y = max_y, label = round(max_y, 2), color = sample), 
             vjust = -0.5, size = 5, fontface = "bold", 
             show.legend = FALSE) +
  scale_color_manual(values = c("NP5" = "lightpink", "NP29" = "red", "NP31" = "purple", "PM50" = "lightblue", "NP70" = "green3","Mean" = "black")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5)))



