# library(VoltRon)
library(ggpubr)
library(magick)

# subset
Vis <- importVisium("../../../../../data/10X_Xenium_Visium/Visium/", sample_name = "VisiumR1", resolution_level = "hires")
Vis_subset <- subset(Vis, interactive = TRUE)
saveRDS(VRBlock_subset$subset_info_list[[1]], file = "Visium_vis_subset.rds")

# image
g1 <- magick::image_ggplot(vrImages(Vis))
ggsave("figures/breastcancer_vis.pdf", plot = g1, device = "pdf", width = 6, height = 7)