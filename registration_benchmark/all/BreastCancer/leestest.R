library(ggplot2)
library(reshape2)
library(dplyr)

#### 
# selected genes ####
####

# colors
# "#440154"
# "#21908C"
# "#FDE725"

datax <- read.table("leestest.txt", row.names = 1, header = TRUE)
datax <- data.frame(datax, method = rownames(datax))
datax <- reshape2::melt(datax)
tmp <- datax$method 
tmp <- factor(tmp, levels = c("janesick", "stalign", "wsireg", "voltron"))
# levels(tmp) <- c("janesick", "stalign", "wsireg", "voltron")
datax$method <- tmp
ggplot(datax, aes(x = method, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_classic() +
  xlab("") + ylab("") + 
  theme(legend.position = "top", legend.title = element_blank()) +
  scale_fill_manual(values = c("#440154", "#21908C")) +
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# ggsave("figures/breastcancer_leestest.pdf", plot = last_plot(), device = "pdf", width = 4, height = 5)

#### 
# all genes ####
####

# datag <- data.frame()
# for(file in list.files("data/")[1:2]){
#   cur_data <- read.table("data/wsireg.txt", row.names = 2)
#   cur_data <- cur_data[,-1]
#   cur_data <- data.frame(cur_data, group = gsub(".txt", "", file))
#   datag <- rbind(datag, cur_data)
# }
# cur_data <- datag

cur_data <- read.table("data/wsireg.txt", row.names = 2)
cur_data <- cur_data[,-1]
cur_data <- data.frame(cur_data, group = "wsireg")
cur_data2 <- read.table("data/stalign.txt", row.names = 2)
cur_data2 <- cur_data2[,-1]
cur_data2 <- data.frame(cur_data2, group = "stalign")
cur_data3 <- read.table("data/janesick.txt", row.names = 2)
cur_data3 <- cur_data3[,-1]
cur_data3 <- data.frame(cur_data3, group = "janesick")
cur_data4 <- read.table("data/voltron.txt", row.names = 2)
cur_data4 <- cur_data4[,-1]
cur_data4 <- data.frame(cur_data4, group = "VoltRon")
datag <- rbind(cur_data, cur_data2, cur_data3, cur_data4)

ggplot(datag, aes(x = `Lee.s.L.statistic`, fill = group, col = group)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.7,
    bins = 30
  ) +
  geom_density(
    linewidth = 0.4, alpha = 0.7
  ) +
  theme_classic()
