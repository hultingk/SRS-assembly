########
## SCRIPT NAME: 04b_cta_directionality.R
## AUTHOR: Katherine Hulting
## PURPOSE: Calculate trajectory directionality, repeat within dispersal mode groups
## PRODUCTS: tableS7.html, tableS8.html, table S9.html, figure4.pdf
#########

### loading libraries
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, 
                 ggeffects, AICcmodavg, multcomp, multcompView, here, kableExtra) # Install missing packages and load needed libraries

source(here::here(file.path("scripts", "04_cta_segments.R")))

########################
#### ALL SPECIES ####
########################
## directionality broken into time periods 
###### directionality in first 12 years ###
sp_info_1_12 <- srs_data_wider %>%
  filter(time <= 12)

patch_info_1_12 <- sp_info_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_1_12 <- sp_info_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
jaccard_dist_1_12 <- vegdist(sp_info_1_12, method = "jaccard")

# defining trajectories
srs_trajectory_1_12 <- defineTrajectories(jaccard_dist_1_12, sites = patch_info_1_12$unique_id, surveys = patch_info_1_12$time)

# directionality
segment_direction_1_12 <- trajectoryDirectionality(srs_trajectory_1_12)
segment_direction_1_12 <- data.frame(segment_direction_1_12)
segment_direction_1_12 <- segment_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = segment_direction_1_12)


###### directionality in second 12 years ###
sp_info_13_24 <- srs_data_wider %>%
  filter(time >= 13) %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  filter(block != "54N") %>% # 54N only has 2 years of surveys after year 12 -- not enough to calculate directionality
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) 
  

patch_info_13_24 <- sp_info_13_24 %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
sp_info_13_24 <- sp_info_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year", "block", "patch", "patch_type"))

# Jaccard distance matrix
jaccard_dist_13_24 <- vegdist(sp_info_13_24, method = "jaccard")

# defining trajectories
srs_trajectory_13_24 <- defineTrajectories(jaccard_dist_13_24, sites = patch_info_13_24$unique_id, surveys = patch_info_13_24$time)

# directionality
segment_direction_13_24 <- trajectoryDirectionality(srs_trajectory_13_24)
segment_direction_13_24 <- data.frame(segment_direction_13_24)
segment_direction_13_24 <- segment_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = segment_direction_13_24)


#### putting all together
segment_direction_all <- rbind(
  segment_direction_1_12,
  segment_direction_13_24
)
segment_direction_all$dispersal_mode <- "All Species"

###### modeling #####
m.direction <- glmmTMB(directionality ~ patch_type * time + (1|block),
                       data = segment_direction_all)
summary(m.direction)
plot(simulateResiduals(m.direction))
performance::r2(m.direction)
anova.direction <- Anova(m.direction, type = "III")

# posthoc tests
m.direction.posthoc <- emmeans(m.direction, ~ patch_type*time)
m.direction_pairs <- pairs(m.direction.posthoc, simple = "patch_type")
m.direction_pairs
m.direction_pairs2 <- pairs(m.direction.posthoc, simple = "time")
m.direction_pairs2
# percent change from decade 1 to decade 2
(-0.035256)/0.368895 * 100 # -9.557191% decrease in directionality
confint(m.direction)
(-0.047920779)/(0.359860618) *100 # -13.31648%
(-0.022590546)/(0.377929254) *100 # -5.977454%






###############################
###### ANIMAL #####
###############################
animal_data$time <- as.numeric(animal_data$time)
# first 12 years
animal_1_12 <- animal_data %>%
  filter(time <= 12)

# patch info
animal_patch_info_1_12 <- animal_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
animal_1_12 <- animal_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
animal_jaccard_dist_1_12 <- vegdist(animal_1_12, method = "jaccard")
# defining trajectories
animal_trajectory_1_12 <- defineTrajectories(animal_jaccard_dist_1_12, sites = animal_patch_info_1_12$unique_id, surveys = animal_patch_info_1_12$time)

# directionality
animal_direction_1_12 <- trajectoryDirectionality(animal_trajectory_1_12)
animal_direction_1_12 <- data.frame(animal_direction_1_12)
animal_direction_1_12 <- animal_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = animal_direction_1_12)


# second half of succession
animal_13_24 <- animal_data %>%
  filter(time >= 13) %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  filter(block != "54N") %>% # 54N only has 2 years of surveys after year 12 -- not enough to calculate directionality
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) 

# patch info
animal_patch_info_13_24 <- animal_13_24 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
animal_13_24 <- animal_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year", "block", "patch", "patch_type"))

# Jaccard distance matrix
animal_jaccard_dist_13_24 <- vegdist(animal_13_24, method = "jaccard")
# defining trajectories
animal_trajectory_13_24 <- defineTrajectories(animal_jaccard_dist_13_24, sites = animal_patch_info_13_24$unique_id, surveys = animal_patch_info_13_24$time)

# directionality
animal_direction_13_24 <- trajectoryDirectionality(animal_trajectory_13_24)
animal_direction_13_24 <- data.frame(animal_direction_13_24)
animal_direction_13_24 <- animal_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = animal_direction_13_24)

#### putting all together
animal_direction_all <- rbind(
  animal_direction_1_12,
  animal_direction_13_24
)
animal_direction_all$dispersal_mode <- "Animal"

###### modeling #####
m.animal_direction <- glmmTMB(directionality ~ patch_type * time + (1|block),
                       data = animal_direction_all)
summary(m.animal_direction)
plot(simulateResiduals(m.animal_direction))
performance::r2(m.animal_direction)
anova.animal.direction <- Anova(m.animal_direction, type = "III")

# posthoc tests
m.animal_direction.posthoc <- emmeans(m.animal_direction, ~ patch_type*time)
m.animal_direction_pairs <- pairs(m.animal_direction.posthoc, simple = "patch_type")
m.animal_direction_pairs
m.animal_direction_pairs2 <- pairs(m.animal_direction.posthoc, simple = "time")
m.animal_direction_pairs2



##################################
###### GRAVITY #####
###################################

gravity_data$time <- as.numeric(gravity_data$time)
# first 12 years
gravity_1_12 <- gravity_data %>%
  filter(time <= 12)

# patch info
gravity_patch_info_1_12 <- gravity_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
gravity_1_12 <- gravity_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
gravity_jaccard_dist_1_12 <- vegdist(gravity_1_12, method = "jaccard")
# defining trajectories
gravity_trajectory_1_12 <- defineTrajectories(gravity_jaccard_dist_1_12, sites = gravity_patch_info_1_12$unique_id, surveys = gravity_patch_info_1_12$time)

# directionality
gravity_direction_1_12 <- trajectoryDirectionality(gravity_trajectory_1_12)
gravity_direction_1_12 <- data.frame(gravity_direction_1_12)
gravity_direction_1_12 <- gravity_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = gravity_direction_1_12)


# second half of succession
gravity_13_24 <- gravity_data %>%
  filter(time >= 13) %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  filter(block != "54N") %>% # 54N only has 2 years of surveys after year 12 -- not enough to calculate directionality
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) 

# patch info
gravity_patch_info_13_24 <- gravity_13_24 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
gravity_13_24 <- gravity_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year", "block", "patch", "patch_type"))

# Jaccard distance matrix
gravity_jaccard_dist_13_24 <- vegdist(gravity_13_24, method = "jaccard")
# defining trajectories
gravity_trajectory_13_24 <- defineTrajectories(gravity_jaccard_dist_13_24, sites = gravity_patch_info_13_24$unique_id, surveys = gravity_patch_info_13_24$time)

# directionality
gravity_direction_13_24 <- trajectoryDirectionality(gravity_trajectory_13_24)
gravity_direction_13_24 <- data.frame(gravity_direction_13_24)
gravity_direction_13_24 <- gravity_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = gravity_direction_13_24)

#### putting all together
gravity_direction_all <- rbind(
  gravity_direction_1_12,
  gravity_direction_13_24
)
gravity_direction_all$dispersal_mode <- "Gravity"

###### modeling #####
m.gravity_direction <- glmmTMB(directionality ~ patch_type * time + (1|block),
                              data = gravity_direction_all)
summary(m.gravity_direction)
plot(simulateResiduals(m.gravity_direction))
performance::r2(m.gravity_direction)
anova.gravity.direction <- Anova(m.gravity_direction, type = "III")

# posthoc tests
m.gravity_direction.posthoc <- emmeans(m.gravity_direction, ~ patch_type*time)
m.gravity_direction_pairs <- pairs(m.gravity_direction.posthoc, simple = "patch_type")
m.gravity_direction_pairs
m.gravity_direction_pairs2 <- pairs(m.gravity_direction.posthoc, simple = "time")
m.gravity_direction_pairs2


#################################
###### WIND #####
#################################

wind_data$time <- as.numeric(wind_data$time)
# first 12 years
wind_1_12 <- wind_data %>%
  filter(time <= 12)

# patch info
wind_patch_info_1_12 <- wind_1_12 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
wind_1_12 <- wind_1_12 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
wind_jaccard_dist_1_12 <- vegdist(wind_1_12, method = "jaccard")
# defining trajectories
wind_trajectory_1_12 <- defineTrajectories(wind_jaccard_dist_1_12, sites = wind_patch_info_1_12$unique_id, surveys = wind_patch_info_1_12$time)

# directionality
wind_direction_1_12 <- trajectoryDirectionality(wind_trajectory_1_12)
wind_direction_1_12 <- data.frame(wind_direction_1_12)
wind_direction_1_12 <- wind_direction_1_12 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 1-12") %>%
  rename(directionality = wind_direction_1_12)


# second half of succession
wind_13_24 <- wind_data %>%
  filter(time >= 13) %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), sep = "-") %>%
  filter(block != "54N") %>% # 54N only has 2 years of surveys after year 12 -- not enough to calculate directionality
  mutate(unique_id = paste(block, patch, patch_type, sep = "-")) 

# patch info
wind_patch_info_13_24 <- wind_13_24 %>%
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
wind_13_24 <- wind_13_24 %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year", "block", "patch", "patch_type"))

# Jaccard distance matrix
wind_jaccard_dist_13_24 <- vegdist(wind_13_24, method = "jaccard")
# defining trajectories
wind_trajectory_13_24 <- defineTrajectories(wind_jaccard_dist_13_24, sites = wind_patch_info_13_24$unique_id, surveys = wind_patch_info_13_24$time)

# directionality
wind_direction_13_24 <- trajectoryDirectionality(wind_trajectory_13_24)
wind_direction_13_24 <- data.frame(wind_direction_13_24)
wind_direction_13_24 <- wind_direction_13_24 %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  mutate(time = "Year 13-24") %>%
  rename(directionality = wind_direction_13_24)

#### putting all together
wind_direction_all <- rbind(
  wind_direction_1_12,
  wind_direction_13_24
)
wind_direction_all$dispersal_mode <- "Wind"

###### modeling #####
m.wind_direction <- glmmTMB(directionality ~ patch_type * time + (1|block),
                              data = wind_direction_all)
summary(m.wind_direction)
plot(simulateResiduals(m.wind_direction))
performance::r2(m.wind_direction)
anova.wind.direction <- Anova(m.wind_direction, type = "III")

# posthoc tests
m.wind_direction.posthoc <- emmeans(m.wind_direction, ~ patch_type*time)
m.wind_direction_pairs <- pairs(m.wind_direction.posthoc, simple = "patch_type")
m.wind_direction_pairs
m.wind_direction_pairs2 <- pairs(m.wind_direction.posthoc, simple = "time")
m.wind_direction_pairs2




###############################
###### TABLES #####
###############################
# Anova table
anova.direction_df <- as.data.frame(anova.direction)
anova.direction_df <- anova.direction_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "All Species")

anova.animal.direction_df <- as.data.frame(anova.animal.direction)
anova.animal.direction_df <- anova.animal.direction_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Animal-Dispersed")

anova.gravity.direction_df <- as.data.frame(anova.gravity.direction)
anova.gravity.direction_df <- anova.gravity.direction_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed")

anova.wind.direction_df <- as.data.frame(anova.wind.direction)
anova.wind.direction_df <- anova.wind.direction_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Wind-Dispersed")

m.direction_anova_all <- rbind(
  anova.direction_df, anova.animal.direction_df, anova.gravity.direction_df, anova.wind.direction_df
)

rename_variable_anova <- tibble(model_term = c("patch_type", "time", "patch_type:time"),
                                Variable = c("Patch Type", "Time Period", "Patch Type:Time Period"))

m.direction_anova_all <- m.direction_anova_all %>%
  filter(model_term != "(Intercept)") %>%
  left_join(rename_variable_anova, by = "model_term") %>%
  dplyr::select(`Dispersal mode`, Variable, Chisq, Df, `Pr(>Chisq)`) %>%
  rename(p.value = `Pr(>Chisq)`, df = Df)

tableS7 <- m.direction_anova_all %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.direction_anova_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.direction_anova_all), extra_css = "padding-bottom: 5px;")
tableS7

# exporting
#save_kable(tableS7, file = file.path("tables", "tableS7.html"))

# emmeans posthoc tables
# creating dataframes of results
# all species
m.direction_pairs_df <- as.data.frame(m.direction_pairs)
m.direction_pairs_df <- m.direction_pairs_df %>%
  mutate(`Dispersal mode` = "All Species") %>%
  mutate(Time = time) 
# animal dispersed
m.animal_direction_pairs_df <- as.data.frame(m.animal_direction_pairs)
m.animal_direction_pairs_df <- m.animal_direction_pairs_df %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") %>%
  mutate(Time = time) 
# gravity dispersed
m.gravity_direction_pairs_df <- as.data.frame(m.gravity_direction_pairs)
m.gravity_direction_pairs_df <- m.gravity_direction_pairs_df %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed") %>%
  mutate(Time = time) 
# wind dispersed
m.wind_direction_pairs_df <- as.data.frame(m.wind_direction_pairs)
m.wind_direction_pairs_df <- m.wind_direction_pairs_df %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") %>%
  mutate(Time = time) 

m_direction_table_all <- rbind(
  m.direction_pairs_df, m.animal_direction_pairs_df, m.gravity_direction_pairs_df, m.wind_direction_pairs_df
)

tableS8 <- m_direction_table_all %>% 
  dplyr::select(`Dispersal mode`, Time, contrast, estimate, SE, df, z.ratio, p.value) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2)) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_direction_table_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_direction_table_all), extra_css = "padding-bottom: 5px;")
tableS8

# exporting
save_kable(tableS8, file = file.path("tables", "tableS8.html"))


# creating dataframes of results
# all species
m.direction_pairs_df2 <- as.data.frame(m.direction_pairs2)
m.direction_pairs_df2 <- m.direction_pairs_df2 %>%
  mutate(`Dispersal mode` = "All Species") %>%
  mutate(`Patch Type` = patch_type) 
# animal dispersed
m.animal_direction_pairs_df2 <- as.data.frame(m.animal_direction_pairs2)
m.animal_direction_pairs_df2 <- m.animal_direction_pairs_df2 %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") %>%
  mutate(`Patch Type` = patch_type) 
# gravity dispersed
m.gravity_direction_pairs_df2 <- as.data.frame(m.gravity_direction_pairs2)
m.gravity_direction_pairs_df2 <- m.gravity_direction_pairs_df2 %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed") %>%
  mutate(`Patch Type` = patch_type) 
# wind dispersed
m.wind_direction_pairs_df2 <- as.data.frame(m.wind_direction_pairs2)
m.wind_direction_pairs_df2 <- m.wind_direction_pairs_df2 %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") %>%
  mutate(`Patch Type` = patch_type) 

m_direction_table_all2 <- rbind(
  m.direction_pairs_df2, m.animal_direction_pairs_df2, m.gravity_direction_pairs_df2, m.wind_direction_pairs_df2
)

tableS9 <- m_direction_table_all2 %>% 
  dplyr::select(`Dispersal mode`, `Patch Type`, contrast, estimate, SE, df, z.ratio, p.value) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_direction_table_all2), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_direction_table_all2), extra_css = "padding-bottom: 5px;")
tableS9

# exporting
save_kable(tableS9, file = file.path("tables", "tableS9.html"))







###############################
###### PLOTTING #####
###############################

# predictions 
m.direction.predict <- ggpredict(m.direction, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.predict$dispersal_mode <- "All Species"
m.direction.animal.predict <- ggpredict(m.animal_direction, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.animal.predict$dispersal_mode <- "Animal"
m.direction.gravity.predict <- ggpredict(m.gravity_direction, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.gravity.predict$dispersal_mode <- "Gravity"
m.direction.wind.predict <- ggpredict(m.wind_direction, terms=c("time [all]", "patch_type [all]"), back_transform = T)
m.direction.wind.predict$dispersal_mode <- "Wind"

# # FACET BY ROWS - total and animal together and gravity and wind together
predict_direction_1 <- rbind(
  m.direction.predict, m.direction.animal.predict
)
predict_direction_2 <- rbind(
  m.direction.wind.predict, m.direction.gravity.predict
)
# making sure factors are in the right order
predict_direction_1$dispersal_mode <- factor(predict_direction_1$dispersal_mode, levels = c("All Species", "Animal"))
predict_direction_2$dispersal_mode <- factor(predict_direction_2$dispersal_mode, levels = c("Wind", "Gravity"))

# joining together data points
dispersal_mode_direction_1 <- rbind(
  segment_direction_all, animal_direction_all
)
dispersal_mode_direction_2 <- rbind(
  wind_direction_all, gravity_direction_all
)

# making sure factors are in the right order
dispersal_mode_direction_1$dispersal_mode <- factor(dispersal_mode_direction_1$dispersal_mode, levels = c("All Species", "Animal"))
dispersal_mode_direction_2$dispersal_mode <- factor(dispersal_mode_direction_2$dispersal_mode, levels = c("Wind", "Gravity"))

# two faceted plots
# first set of plots
direction_predict_plot_1 <- predict_direction_1 %>%
  ggplot() +
  geom_jitter(aes(x = time, y = directionality, color = patch_type),
             data = dispersal_mode_direction_1, alpha = 0.2, size = 5.5,
             position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = predict_direction_1, width = 0, linewidth = 3,  position = position_dodge(width = 0.7)) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("All Species" = "(A) All species", "Animal" = "(B) Animal-dispersed"))) +
   scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  geom_point(aes(x = x, y = predicted, fill = group), size = 7.5, 
             data = predict_direction_1,  position = position_dodge(width = 0.7),
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
direction_predict_plot_1

# second set of plots
direction_predict_plot_2 <- predict_direction_2 %>%
  ggplot() +
  geom_jitter(aes(x = time, y = directionality, color = patch_type),
             data = dispersal_mode_direction_2, alpha = 0.2, size = 5.5,
             position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = predict_direction_2, width = 0, linewidth = 3,  position = position_dodge(width = 0.7)) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("Gravity" = "(C) Gravity-dispersed", "Wind" = "(D) Wind-dispersed"))) +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  geom_point(aes(x = x, y = predicted, fill = group), size = 7.5, 
             data = predict_direction_2,  position = position_dodge(width = 0.7),
             colour="black", pch=21, stroke = 2)+ 
  labs(title = NULL,
       x = "Time period",
       y = "Trajectory directionality") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                    labels = c("Connected", "Rectangular", "Winged"), 
                    name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), 
                     labels = c("Connected", "Rectangular", "Winged"), 
                     name = "Patch Type") +
  theme(axis.text = element_text(size = 18)) +
  theme(legend.position = "none") 
direction_predict_plot_2

# get legend
pL <- predict_direction_2 %>%
  ggplot() +
  geom_jitter(aes(x = time, y = directionality, color = patch_type), 
              data = dispersal_mode_direction_2, alpha = 0.2, size = 5.5, 
              position = position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.7)) +
  geom_errorbar(aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), color = "black",
                data = predict_direction_2, width = 0, linewidth = 2.5,  position = position_dodge(width = 0.7)) +
  geom_point(aes(x = x, y = predicted, fill = group), size = 7.5, 
             data = predict_direction_2,  position = position_dodge(width = 0.7),
             colour="black", pch=21, stroke = 2)+ 
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom") +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), labels = c("Connected", "Rectangular", "Winged"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), labels = c("Connected", "Rectangular", "Winged"), name = "Patch Type")
l <- get_legend(pL)



# put together
figure3 <- cowplot::plot_grid(direction_predict_plot_1, direction_predict_plot_2, l,
                                        ncol = 1, nrow = 3, rel_heights = c(1, 1.1, 0.15),
                                        label_size = 20, label_x = 0.2, label_y = 0.95)
figure3

# # exporting
# pdf(file = file.path("plots", "figure3.pdf"), width = 12.5, height = 13)
# figure3
# dev.off()

