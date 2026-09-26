### Repeating analyses using calender year as predictor
#loading libraries 
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, performance, cowplot, kableExtra, car, scales) # Install missing packages and load needed libraries

source(here::here(file.path("scripts", "00_functions.R")))

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

#### Temporal variability: CTA segment length ####
# pivot to wider format
srs_data_year_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
srs_data_year_wider$time <- as.numeric(srs_data_year_wider$time)
srs_data_year_wider$unique_id <- as.factor(srs_data_year_wider$unique_id)
srs_data_year_wider$year <- as.factor(srs_data_year_wider$year)

# patch data
patch_info_year <- srs_data_year_wider %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_year <- srs_data_year_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_year <- vegdist(sp_info_year, method = "jaccard")

# defining trajectories
srs_trajectory_year <- defineTrajectories(jaccard_dist_year, sites = patch_info_year$unique_id, surveys = patch_info_year$time)

# segment lengths of trajectories between consectutive years
segment_lengths_year <- trajectoryLengths(srs_trajectory_year)
segment_lengths_year <- segment_lengths_year %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)

# creating time info to join with segment lengths - some surveys were not consecutive years
time_surveys_year <- patch_info_year %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) # removing first survey for sites created in 2007

# joining with segment lengths
segment_lengths_year <- cbind(segment_lengths_year, time_surveys_year)
segment_lengths_year$dispersal_mode <- "All Species"
segment_lengths_year <- segment_lengths_year %>% # making calendar year a number for analysis
  mutate(year_numeric = as.numeric(sub("^20", "", year)))
segment_lengths_year$patch_type <- as.factor(segment_lengths_year$patch_type)

# modeling
# quadratic
m_length_quad_year <- glmmTMB(distance ~ patch_type * year_numeric + patch_type * I(year_numeric^2) + (1|block/patch),
                                data = segment_lengths_year)
# model fit
summary(m_length_quad_year)
plot(simulateResiduals(m_length_quad_year))

# posthoc
m_length_posthoc_year <- emmeans(m_length_quad_year, ~ patch_type*year_numeric + patch_type * I(year_numeric^2), at = list(year_numeric = c(12)))
m_length_pairs_year <- pairs(m_length_posthoc_year, simple = "patch_type")
m_length_pairs_year

# model predictions for plotting
m_length_predict_year <- ggpredict(m_length_quad_year, terms = c("year_numeric [all]", "patch_type"))
m_length_predict_year <- as.data.frame(m_length_predict_year)
m_length_predict_year$dispersal_mode <- "All Species"
m_length_predict_year$linetype <- "solid" # adding line type

m_length_predict_year$year <- m_length_predict_year$x + 2000
segment_lengths_year$year_plotting <- segment_lengths_year$year_numeric + 2000

# plotting
all_segment_plot_year <- m_length_predict_year %>%
  ggplot() +
  geom_point(aes(year_plotting, distance, color = patch_type), size = 6, alpha = 0.15, data = segment_lengths_year) +
  geom_ribbon(aes(x = year, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(year, predicted, color = group), linewidth = 3.5) +
  theme_minimal(base_size = 32) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Calendar year") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  #guides(fill=guide_legend(ncol=1)) +
  #guides(color=guide_legend(ncol=1)) +
  #scale_y_continuous(limits = c(0.19, 0.38), labels = label_number(accuracy = 0.01)) +
  theme(axis.text = element_text(size = 20), 
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(legend.position = "right")
all_segment_plot_year


#### Directionality ####
## directionality broken into time periods 
###### directionality in first 12 years ###
srs_data_year_wider <- srs_data_year_wider %>% # making calendar year a number for analysis
  mutate(year_numeric = as.numeric(sub("^20", "", year)))

sp_info_1_12_year <- srs_data_year_wider %>%
  filter(year_numeric <= 12)

patch_info_1_12_year <- sp_info_1_12_year %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_1_12_year <- sp_info_1_12_year %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_1_12_year <- vegdist(sp_info_1_12_year, method = "jaccard")

# defining trajectories
srs_trajectory_1_12_year <- defineTrajectories(jaccard_dist_1_12_year, sites = patch_info_1_12_year$unique_id, surveys = patch_info_1_12_year$time)

# directionality
segment_direction_1_12_year <- trajectoryDirectionality(srs_trajectory_1_12_year)
segment_direction_1_12_year <- data.frame(segment_direction_1_12_year)
segment_direction_1_12_year <- segment_direction_1_12_year %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = segment_direction_1_12_year)

###### directionality in second 12 years ###
sp_info_13_24_year <- srs_data_year_wider %>%
  filter(year_numeric >= 13) %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  filter(block != "54N") %>% # 54N only has 2 years of surveys after year 12 -- not enough to calculate directionality
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) 

patch_info_13_24_year <- sp_info_13_24_year %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_13_24_year <- sp_info_13_24_year %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year", "block", "patch", "patch_type"))

# Jaccard distance matrix
jaccard_dist_13_24_year <- vegdist(sp_info_13_24_year, method = "jaccard")

# defining trajectories
srs_trajectory_13_24_year <- defineTrajectories(jaccard_dist_13_24_year, sites = patch_info_13_24_year$unique_id, surveys = patch_info_13_24_year$time)

# directionality
segment_direction_13_24_year <- trajectoryDirectionality(srs_trajectory_13_24_year)
segment_direction_13_24_year <- data.frame(segment_direction_13_24_year)
segment_direction_13_24_year <- segment_direction_13_24_year %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = segment_direction_13_24_year)

#### putting all together
segment_direction_all_year <- rbind(
  segment_direction_1_12_year,
  segment_direction_13_24_year
)
segment_direction_all_year$dispersal_mode <- "All Species"

# modeling
m.direction_year <- glmmTMB(directionality ~ patch_type * time + (1|block),
                              data = segment_direction_all_year)
summary(m.direction_year)
plot(simulateResiduals(m.direction_year))
performance::r2(m.direction_year)

# posthoc tests
m.direction.posthoc_year <- emmeans(m.direction_year, ~ patch_type*time)
m.direction_pairs_year <- pairs(m.direction.posthoc_year, simple = "patch_type")
m.direction_pairs_year
m.direction_pairs2_year <- pairs(m.direction.posthoc_year, simple = "time")
m.direction_pairs2_year

# predictions for plotting
m.direction.predict_year <- ggpredict(m.direction_year, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.predict_year$dispersal_mode <- "All Species"

# plotting
direction_plot_year <- m.direction.predict_year %>%
  ggplot() +
  geom_jitter(aes(x = time, y = directionality, color = patch_type),
             data = segment_direction_all_year, alpha = 0.2, size = 5.5,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = m.direction.predict_year, width = 0, linewidth = 3,  position = position_dodge(width = 0.7)) +
  #scale_y_continuous(limits = c(0.32, 0.38), labels = label_number(accuracy = 0.01)) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  geom_point(aes(x = x, y = predicted, fill = group), size = 7.5, 
             data = m.direction.predict_year,  position = position_dodge(width = 0.7),
             colour="black", pch=21, stroke = 2)+ 
  labs(title = NULL,
       x = "Time Period",
       y = "Trajectory directionality") +
  scale_x_discrete(
    labels = c(
      "Year 1-12" = "Year 2000-2012",
      "Year 13-24" = "Year 2013-2024"
    )) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                    labels = c("Connected", "Rectangular", "Winged"), 
                    name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                     labels = c("Connected", "Rectangular", "Winged"), 
                     name = "Patch Type") +
  theme(axis.text = element_text(size = 18)) +
  theme(legend.position = "right") 
direction_plot_year


#### Convergence/divergence: Spatial beta diversity over time ####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
convergence_jaccard_year <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

year_info <- srs_data %>%
  count(block, time, year) %>%
  dplyr::select(-n) %>%
  mutate(year_numeric = as.numeric(sub("^20", "", year)))
  

# removing same patch type comparisons and time 0 (only for 52 and 57)
convergence_jaccard_year <- convergence_jaccard_year %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "All Species") %>%
  left_join(year_info, by = c("block", "time"))
convergence_jaccard_year$patch_pair <- as.factor(convergence_jaccard_year$patch_pair)

#  model
m.converge_quad_year <- glmmTMB(jaccard ~ patch_pair * year_numeric + patch_pair * I(year_numeric^2) + (1|block),
                                  data = convergence_jaccard_year)
## model checking
summary(m.converge_quad_year)

# model checking
plot(simulateResiduals(m.converge_quad_year))

## posthoc comparisons
m.converge_posthoc_year <- emmeans(m.converge_quad_year, ~ patch_pair*year_numeric+ patch_pair * I(year_numeric^2))
m.converge_pairs_year <- pairs(m.converge_posthoc_year, simple = "patch_pair")
m.converge_pairs_year

# model predictions
m.converge.predict_year <- ggpredict(m.converge_quad_year, terms=c("year_numeric [all]", "patch_pair [all]"), back_transform = T)
m.converge.predict_year <- as.data.frame(m.converge.predict_year)
m.converge.predict_year$dispersal_mode <- "All Species"


m.converge.predict_year$year <- m.converge.predict_year$x + 2000
convergence_jaccard_year$year_plotting <- convergence_jaccard_year$year_numeric + 2000

# plotting
convergence_plot_year <- m.converge.predict_year %>%
  ggplot() +
  geom_point(aes(year_plotting, jaccard, color = patch_pair), size = 6, alpha = 0.07, data = convergence_jaccard_year) +
  geom_ribbon(aes(x = year, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(year, predicted, color = group), linewidth = 3.5) +
  theme_minimal(base_size = 28) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), labels = c(expression("Connected"%<->%"Rectangular"), 
                                                                            expression("Connected"%<->%"Winged"),
                                                                            expression("Rectangular"%<->%"Winged")), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), labels = c(expression("Connected"%<->%"Rectangular"), 
                                                                             expression("Connected"%<->%"Winged"),
                                                                             expression("Rectangular"%<->%"Winged")), name = "Patch Comparison") +
  xlab("Calendar year") +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(legend.position = "right") 
convergence_plot_year


#### exporting plots ####
pdf(file = file.path("plots", "calendar_year_length.pdf"), width = 14, height = 8)
all_segment_plot_year
dev.off()

pdf(file = file.path("plots", "calendar_year_direction.pdf"), width = 12, height = 7)
direction_plot_year
dev.off()

pdf(file = file.path("plots", "calendar_year_convergence.pdf"), width = 15.5, height = 8)
convergence_plot_year
dev.off()



