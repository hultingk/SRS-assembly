########
## SCRIPT NAME: 02_pcoa.R
## AUTHOR: Katherine Hulting
## PURPOSE: Visualize PCoA of community trajectories
## PRODUCTS: figureS1.pdf : supplementary figure of PCoA visualization
#########


# loading libraries
librarian::shelf(tidyverse, vegan, ape, BiodiversityR, glmmTMB)

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))


# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.factor(srs_data_wider$time)
srs_data_wider$unique_id <- as.factor(srs_data_wider$unique_id)
srs_data_wider$year <- as.factor(srs_data_wider$year)

# patch data
patch_info <- srs_data_wider %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info <- srs_data_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))


# Jaccard distance matrix
jaccard_dist <- vegdist(sp_info, method = "jaccard")
jaccard_dist_all <- as.matrix(jaccard_dist)
jaccard_dist_all_df <- as.data.frame(jaccard_dist_all)

jaccard_dist_all_df <- cbind(patch_info, jaccard_dist_all_df) # merge with patch info

# pcoa
pcoa_all <- pcoa(jaccard_dist)

# add axes to patch and time
pcoa_axes <- pcoa_all$vectors[,c(1,2)]
pcoa_axes <- cbind(patch_info, pcoa_axes)

#### plotting ####
pcoa_axes_plot <- pcoa_axes %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-", remove = F)
pcoa_axes_plot$time <- as.numeric(pcoa_axes_plot$time)

# separating into group of blocks with duplicate rectanglar patches and group of blocks with duplicate winged patches
duplicate_wing <- pcoa_axes_plot %>%
  filter(block %in% c("08", "52", "53N", "53S", "75E"))
duplicate_rectangle <- pcoa_axes_plot %>%
  filter(block %in% c("10", "54N", "54S", "57", "75W"))

# PCOA plot for blocks with duplicate winged patch
duplicate_wing_names <- c( # assinging patch rep letters to patch types
  'B'="Connected",
  'C'="Winged",
  'D'="Rectangular",
  'E'="Winged"
)
duplicate_wing_pcoa <- duplicate_wing %>%
  mutate(patch_rep = factor(patch_rep, levels = c("B", "D", "C", "E"))) %>%
  ggplot(aes(Axis.1, Axis.2, color = time)) +
  geom_point(size = 1.5) +
  geom_path(aes(Axis.1, Axis.2, group = patch_type, color = time), linewidth = 1) +
  facet_grid(block~patch_rep, labeller = labeller(patch_rep = duplicate_wing_names)) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        strip.text.x= element_text(size = 22),
        strip.text.y= element_text(size = 22),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = "none") +
  labs(color='Time',
       x = "PCoA 1",
       y = "PCoA 2") 
duplicate_wing_pcoa


# PCOA plot for blocks with duplicate rectangular patches
duplicate_rectangle_names <- c( # assinging patch rep letters to patch types
  'B'="Connected",
  'C'="Rectangular",
  'D'="Winged",
  'E'="Rectangular"
)
duplicate_rectangle_pcoa <- duplicate_rectangle %>%
  mutate(patch_rep = factor(patch_rep, levels = c("B", "D", "C", "E"))) %>%
  ggplot(aes(Axis.1, Axis.2, color = time)) +
  geom_point(size = 1.5) +
  geom_path(aes(Axis.1, Axis.2, group = patch_type, color = time), linewidth = 1) +
  facet_grid(block~patch_rep, labeller = labeller(patch_rep = duplicate_rectangle_names)) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        strip.text.x= element_text(size = 22),
        strip.text.y= element_text(size = 22),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(color='Time',
       x = "PCoA 1",
       y = "PCoA 2") 
duplicate_rectangle_pcoa

# putting plots together
figureS1 <- cowplot::plot_grid(duplicate_wing_pcoa, duplicate_rectangle_pcoa, rel_widths = c(1, 1.1))
figureS1

# exporting
# pdf(file = file.path("plots", "figureS1.pdf"), width = 22, height =12.5)
# figureS1
# dev.off()




