library(RBioFormats)
library(peakRAM)
library(bench)
library(magick)
library(dplyr)
library(ggplot2)
library(scales)
library(xlsx)
source("../functions.R")

####
# import ####
####

# get results
all_results <- read.xlsx("registration_speed_benchmark_other_methods.xlsx", sheetIndex = 1)
all_results <- reshape2::melt(all_results, measure.vars = c("Time_sec", "PeakmemoryMiB"))

# all_results$Dataset <- paste0(all_results$Dataset, "\n (", all_results$Step, ")")

####
# visualize time ####
####

# visualize time 
all_results_time <- all_results[all_results$variable %in% "Time_sec",]
breaks <- c(5,10,200,400,800)
all_results_time$value <- transform_breaks_time(all_results_time$value, breaks)
tmp <- all_results_time$Method
unitmp <- unique(all_results_time$Method)
unitmp <- c(unitmp[grepl("PASTE", unitmp)], "SANTO", "STalign", "wsireg", "VoltRon")
tmp <- factor(tmp, levels = unitmp)
all_results_time$Method <- tmp
tmp <- all_results_time$Step
tmp <- factor(tmp, levels = c("Align. + Cells/Spots", "Cells/Spots", "Molecules", "Image (Align.)", "Image (Trans.)"))
all_results_time$Step <- tmp
tmp <- all_results_time$Dataset
unitmp <- unique(all_results_time$Dataset)
unitmp <- c(unitmp[grepl("DLPFC", unitmp)], unitmp[!grepl("DLPFC", unitmp)])
tmp <- factor(tmp, levels = unitmp)
all_results_time$Dataset <- tmp
g1 <- ggplot(mapping = aes(x = Dataset, y = value, fill = Method), data = all_results_time) + 
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.05), width = 1.0) + 
  facet_grid(~Step, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    limits = c(0, length(breaks)),
    breaks = 0:length(breaks),   # <-- choose your custom break points
    labels = as.character(c(0,breaks))
  ) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(legend.position = "top") + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE)) + 
  ylab("Time (seconds)") + xlab("")

####
# visualize memory ####
####

# visualize time 
all_results_memory <- all_results[all_results$variable %in% "PeakmemoryMiB",]
all_results_memory$value <- all_results_memory$value * (1024^2) # output is MB, convert to byte first
breaks <- c(1,1024, 1024^2, 1024^3, 4*1024^3)
tmp <- all_results_time$Method
unitmp <- unique(all_results_time$Method)
unitmp <- c(unitmp[grepl("PASTE", unitmp)], "SANTO", "STalign", "wsireg", "VoltRon")
tmp <- factor(tmp, levels = unitmp)
all_results_time$Method <- tmp
tmp <- all_results_memory$Step
tmp <- factor(tmp, levels = c("Align. + Cells/Spots", "Cells/Spots", "Molecules", "Image (Align.)", "Image (Trans.)"))
all_results_memory$Step <- tmp
tmp <- all_results_memory$Dataset
unitmp <- unique(all_results_memory$Dataset)
unitmp <- c(unitmp[grepl("DLPFC", unitmp)], unitmp[!grepl("DLPFC", unitmp)])
tmp <- factor(tmp, levels = unitmp)
all_results_memory$Dataset <- tmp
g2 <- ggplot(mapping = aes(x = Dataset, y = value, fill = Method), data = all_results_memory) + 
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single", padding = 0.05), width = 1.0) + 
  facet_grid(~Step, scales = "free_x", space = "free_x") +
  scale_y_continuous(
    trans = "log2",
    breaks = c(1, 1024, 1024^2, 1024^3, 4*1024^3),
    labels = c("1B", "1KB", "1MB", "1GB", "4GB"),
    limits = c(1, 4*1024^3)
  ) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(legend.position = "none") +
  ylab("Memory Consumption (Bytes)") + xlab("") +
  scale_x_discrete(expand = expansion(add = 0.4)) + 
  theme(panel.spacing.x = unit(0.5, "lines"))

g1 / g2
# ggsave("time_memory.pdf", plot = last_plot(), device = "pdf",
#        height = 8, width = 15)