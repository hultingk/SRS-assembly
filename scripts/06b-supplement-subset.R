### Repeating analyses with only blocks that were active for the full time series 
#loading libraries 
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, performance, cowplot, kableExtra, car, scales) # Install missing packages and load needed libraries

source(here::here(file.path("scripts", "00_functions.R")))

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))

# subsetting only blocks that were active the full time series
srs_data_subset <- srs_data %>%
  filter(block %in% c('08', '10', "53N", '53S', "54S"))


#### Temporal variability: CTA segment length ####
# pivot to wider format
srs_data_subset_wider <- srs_data_subset %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
srs_data_subset_wider$time <- as.numeric(srs_data_subset_wider$time)
srs_data_subset_wider$unique_id <- as.factor(srs_data_subset_wider$unique_id)
srs_data_subset_wider$year <- as.factor(srs_data_subset_wider$year)

# patch data
patch_info_subset <- srs_data_subset_wider %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_subset <- srs_data_subset_wider %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_subset <- vegdist(sp_info_subset, method = "jaccard")

# defining trajectories
srs_trajectory_subset <- defineTrajectories(jaccard_dist_subset, sites = patch_info_subset$unique_id, surveys = patch_info_subset$time)

# segment lengths of trajectories between consectutive years
segment_lengths_subset <- trajectoryLengths(srs_trajectory_subset)
segment_lengths_subset <- segment_lengths_subset %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)

# creating time info to join with segment lengths - some surveys were not consecutive years
time_surveys_subset <- patch_info_subset %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) # removing first survey for sites created in 2007

# joining with segment lengths
segment_lengths_subset <- cbind(segment_lengths_subset, time_surveys_subset)
segment_lengths_subset$dispersal_mode <- "All Species"
segment_lengths_subset$s.time <- as.numeric(scale(segment_lengths_subset$time)) # scaling time
segment_lengths_subset$patch_type <- as.factor(segment_lengths_subset$patch_type)

# modeling
# quadratic
m_length_quad_subset <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                         data = segment_lengths_subset)

# model fit
summary(m_length_quad_subset)
plot(simulateResiduals(m_length_quad_subset))
#check_model(m_length_quad)
performance::r2(m_length_quad_subset)

# posthoc
m_length_posthoc_subset <- emmeans(m_length_quad_subset, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0)))
m_length_pairs_subset <- pairs(m_length_posthoc_subset, simple = "patch_type")
m_length_pairs_subset

# model predictions for plotting
m_length_predict_subset <- ggpredict(m_length_quad_subset, terms = c("s.time [all]", "patch_type"))
m_length_predict_subset <- as.data.frame(m_length_predict_subset)
m_length_predict_subset$dispersal_mode <- "All Species"
m_length_predict_subset$linetype <- "solid" # adding line type
scaled_time_key_subset <- segment_lengths_subset %>% # creating key of scaled times to join to predictions for easy visualization
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m_length_predict_subset <- m_length_predict_subset %>%
  left_join(scaled_time_key_subset, by = c("x" = "s.time"))

# plotting
all_segment_plot_subset <- m_length_predict_subset %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 6, alpha = 0.15, data = segment_lengths_subset) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(time, predicted, color = group), linewidth = 3.5) +
  theme_minimal(base_size = 32) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Years since site creation") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  #guides(fill=guide_legend(ncol=1)) +
  #guides(color=guide_legend(ncol=1)) +
  #scale_y_continuous(limits = c(0.19, 0.38), labels = label_number(accuracy = 0.01)) +
  theme(axis.text = element_text(size = 20), 
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(legend.position = "top")
all_segment_plot_subset



#### Directionality ####
## directionality broken into time periods 
###### directionality in first 12 years ###
sp_info_1_12_subset <- srs_data_subset_wider %>%
  filter(time <= 12)

patch_info_1_12_subset <- sp_info_1_12_subset %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_1_12_subset <- sp_info_1_12_subset %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_1_12_subset <- vegdist(sp_info_1_12_subset, method = "jaccard")

# defining trajectories
srs_trajectory_1_12_subset <- defineTrajectories(jaccard_dist_1_12_subset, sites = patch_info_1_12_subset$unique_id, surveys = patch_info_1_12_subset$time)

# directionality
segment_direction_1_12_subset <- trajectoryDirectionality(srs_trajectory_1_12_subset)
segment_direction_1_12_subset <- data.frame(segment_direction_1_12_subset)
segment_direction_1_12_subset <- segment_direction_1_12_subset %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = segment_direction_1_12_subset)

###### directionality in second 12 years ###
sp_info_13_24_subset <- srs_data_subset_wider %>%
  filter(time >= 13) %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  filter(block != "54N") %>% # 54N only has 2 years of surveys after year 12 -- not enough to calculate directionality
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) 

patch_info_13_24_subset <- sp_info_13_24_subset %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_13_24_subset <- sp_info_13_24_subset %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year", "block", "patch", "patch_type"))

# Jaccard distance matrix
jaccard_dist_13_24_subset <- vegdist(sp_info_13_24_subset, method = "jaccard")

# defining trajectories
srs_trajectory_13_24_subset <- defineTrajectories(jaccard_dist_13_24_subset, sites = patch_info_13_24_subset$unique_id, surveys = patch_info_13_24_subset$time)

# directionality
segment_direction_13_24_subset <- trajectoryDirectionality(srs_trajectory_13_24_subset)
segment_direction_13_24_subset <- data.frame(segment_direction_13_24_subset)
segment_direction_13_24_subset <- segment_direction_13_24_subset %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = segment_direction_13_24_subset)

#### putting all together
segment_direction_all_subset <- rbind(
  segment_direction_1_12_subset,
  segment_direction_13_24_subset
)
segment_direction_all_subset$dispersal_mode <- "All Species"

# modeling
m.direction_subset <- glmmTMB(directionality ~ patch_type * time + (1|block),
                       data = segment_direction_all_subset)
summary(m.direction_subset)
plot(simulateResiduals(m.direction_subset))
performance::r2(m.direction_subset)

# posthoc tests
m.direction.posthoc_subset <- emmeans(m.direction_subset, ~ patch_type*time)
m.direction_pairs_subset <- pairs(m.direction.posthoc_subset, simple = "patch_type")
m.direction_pairs_subset
m.direction_pairs2_subset <- pairs(m.direction.posthoc_subset, simple = "time")
m.direction_pairs2_subset

# predictions for plotting
m.direction.predict_subset <- ggpredict(m.direction_subset, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.predict_subset$dispersal_mode <- "All Species"

# plotting
direction_plot_subset <- m.direction.predict_subset %>%
  ggplot() +
  # geom_jitter(aes(x = time, y = directionality, color = patch_type),
  #            data = dispersal_mode_direction_1, alpha = 0.2, size = 5.5,
  #            position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = m.direction.predict_subset, width = 0, linewidth = 3,  position = position_dodge(width = 0.7)) +
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
             data = m.direction.predict_subset,  position = position_dodge(width = 0.7),
             colour="black", pch=21, stroke = 2)+ 
  labs(title = NULL,
       x = NULL,
       y = "Trajectory directionality") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                    labels = c("Connected", "Rectangular", "Winged"), 
                    name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                     labels = c("Connected", "Rectangular", "Winged"), 
                     name = "Patch Type") +
  theme(axis.text = element_text(size = 18)) +
  theme(legend.position = "none") 
direction_plot_subset



#### Convergence/divergence: Spatial beta diversity over time ####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
convergence_jaccard_subset <- srs_data_subset %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
convergence_jaccard_subset <- convergence_jaccard_subset %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "All Species")
convergence_jaccard_subset$s.time <- as.numeric(scale(convergence_jaccard_subset$time)) # scaling time
convergence_jaccard_subset$patch_pair <- as.factor(convergence_jaccard_subset$patch_pair)

#  model
m.converge_quad_subset <- glmmTMB(jaccard ~ patch_pair * s.time + patch_pair * I(s.time^2) + (1|block),
                           data = convergence_jaccard_subset)
## model checking
summary(m.converge_quad_subset)

# model checking
plot(simulateResiduals(m.converge_quad_subset))

## posthoc comparisons
m.converge_posthoc_subset <- emmeans(m.converge_quad_subset, ~ patch_pair*s.time+ patch_pair * I(s.time^2))
m.converge_pairs_subset <- pairs(m.converge_posthoc_subset, simple = "patch_pair")
m.converge_pairs_subset

# model predictions
m.converge.predict_subset <- ggpredict(m.converge_quad_subset, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge.predict_subset <- as.data.frame(m.converge.predict_subset)
m.converge.predict_subset$dispersal_mode <- "All Species"

# creating time key for easy visualization
scaled_time_key_subset <- convergence_jaccard_subset %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

# plotting
convergence_plot_subset <- m.converge.predict_subset %>%
  left_join(scaled_time_key_subset, by = c("x" = "s.time")) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.07, data = convergence_jaccard_subset) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(time, predicted, color = group), linewidth = 3.5) +
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
  xlab("Years since site creation") +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 16),
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(legend.position = "right") 
convergence_plot_subset



