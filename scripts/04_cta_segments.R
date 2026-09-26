########
## SCRIPT NAME: 04a_cta_segments.R
## AUTHOR: Katherine Hulting
## PURPOSE: Calculate interannual trajectory distances, repeat within dispersal mode groups
## PRODUCTS: tableS4.html, tableS5.html, tableS6.html, figure3.pdf
#########

### community trajectory analysis - segment lengths
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, performance, cowplot, kableExtra, car, scales) # Install missing packages and load needed libraries

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))


########################
#### ALL SPECIES ####
########################
# pivot to wider format
srs_data_wider <- srs_data %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format


# make factor
srs_data_wider$time <- as.numeric(srs_data_wider$time)
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

# defining trajectories
srs_trajectory <- defineTrajectories(jaccard_dist, sites = patch_info$unique_id, surveys = patch_info$time)


# segment lengths of trajectories between consectutive years
segment_lengths <- trajectoryLengths(srs_trajectory)
segment_lengths <- segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)

# creating time info to join with segment lengths - some surveys were not consecutive years
time_surveys <- patch_info %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) # removing first survey for sites created in 2007

# joining with segment lengths
segment_lengths <- cbind(segment_lengths, time_surveys)
segment_lengths$dispersal_mode <- "All Species"
segment_lengths$s.time <- as.numeric(scale(segment_lengths$time)) # scaling time
segment_lengths$patch_type <- as.factor(segment_lengths$patch_type)


###### MODELS ######
# segment lengths model 
# linear
m_length <- glmmTMB(distance ~ patch_type * s.time + (1|block/patch),
                    data = segment_lengths)
# quadratic
m_length_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                    data = segment_lengths)
# null
m_length_null <- glmmTMB(distance ~ 1 + (1|block/patch), # null model
                         data = segment_lengths)

# AIC comparison
length.aic <- list(m_length, m_length_quad, m_length_null)
length.aic.table <- aictab(length.aic) # quadratic much better fit
length.aic.table

# model fit
summary(m_length_quad)
plot(simulateResiduals(m_length_quad))
#check_model(m_length_quad)
performance::r2(m_length_quad)
anova.length <- Anova(m_length_quad, type = "III")


# posthoc
m_length_posthoc <- emmeans(m_length_quad, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0)))
m_length_pairs <- pairs(m_length_posthoc, simple = "patch_type")
m_length_pairs

# percent change in segment length from year 2-22 (20 years)
# time 2 = -1.5077466
# time 22 = 1.5664900
# intercept + time_estimate (time) + time^2_estimate (time)
# at time 2
0.241442 + -0.017123*(-1.5077466) + 0.014279*((-1.5077466)^2) # 0.2997196
# at time 22
0.241442 + -0.017123*(1.5664900) + 0.014279*((1.5664900)^2) # 0.2496581

# percent change = time 22 - time 2 / time 2 * 100
(0.2496581 - 0.2997196)/0.2997196 * 100 #-16.70278 % decrease

# 95% CI, percent change
confint(m_length_quad)
# lower 95% CI
# at time 1
0.214760382 + -0.025648096*(-1.5077466) + 0.006113128*((-1.5077466)^2) # 0.2673282
# at time 21
0.214760382 + -0.025648096*(1.5664900) + 0.006113128*((1.5664900)^2) # 0.1895838
# percent change = time 21 - time 1 / time 1 * 100
(0.1895838 - 0.2673282)/0.2673282 * 100 #-29.082 % 

# upper 95% CI
# at time 1
0.268123543 + -0.008598510*(-1.5077466) + 0.022445636*((-1.5077466)^2) # 0.3321136
# at time 21
0.268123543 + -0.008598510*(1.5664900) + 0.022445636*((1.5664900)^2) # 0.3097332
# percent change = time 21 - time 1 / time 1 * 100
(0.3097332 - 0.3321136)/0.3321136 * 100 #-6.738779 % 




#####################
#### ANIMAL ####
#####################

##### animal CTA segments #####
## animal dispersed 
animal_data <- srs_data %>%
  filter(dispersal_mode == "Animal") %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

  
# make factor
animal_data$unique_id <- as.factor(animal_data$unique_id)

# patch data
animal_patch_info <- animal_data %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
animal_sp_info <- animal_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
animal_jaccard_dist <- vegdist(animal_sp_info, method = "jaccard")

# defining trajectories
animal_srs_trajectory <- defineTrajectories(animal_jaccard_dist, sites = animal_patch_info$unique_id, surveys = animal_patch_info$time)

# segment lengths of trajectories between consectutive years
animal_segment_lengths <- trajectoryLengths(animal_srs_trajectory)
animal_segment_lengths <- animal_segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)


# creating time info to join with segment lengths - some surveys were not consecutive years
animal_time_surveys <- animal_patch_info %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), remove = F) %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) %>% # removing first survey for sites created in 2007
  dplyr::select(!c("block", "patch", "patch_type"))


# joining with segment lengths
animal_segment_lengths <- cbind(animal_segment_lengths, animal_time_surveys)
animal_segment_lengths$dispersal_mode <- "Animal"
animal_segment_lengths$s.time <- as.numeric(scale(as.numeric(animal_segment_lengths$time))) # scaling time
animal_segment_lengths$patch_type <- as.factor(animal_segment_lengths$patch_type)


###### MODELS ######
# segment lengths model 
# linear animal dispersed
m_length_animal <- glmmTMB(distance ~ patch_type * s.time + (1|block/patch),
                                data = animal_segment_lengths)
# quadratic animal dispersed
m_length_animal_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                         data = animal_segment_lengths)
# null animal dispersed
m_length_animal_null <- glmmTMB(distance ~ 1 + (1|block/patch),
                           data = animal_segment_lengths)
# AIC comparison
length.aic.animal <- list(m_length_animal, m_length_animal_quad, m_length_animal_null)
length.aic.animal.table <- aictab(length.aic.animal) # null is best fit
length.aic.animal.table

# model fit
summary(m_length_animal_null)
plot(simulateResiduals(m_length_animal_null))
#check_model(m_length_animal_null)
performance::r2(m_length_animal_null)
anova.animal.length <- Anova(m_length_animal_null, type = "III")










#####################
#### GRAVITY ####
#####################

##### gravity CTA segments #####
gravity_data <- srs_data %>%
  filter(dispersal_mode == "Gravity") %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
gravity_data$unique_id <- as.factor(gravity_data$unique_id)

# patch data
gravity_patch_info <- gravity_data %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
gravity_sp_info <- gravity_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))
# Jaccard distance matrix
gravity_jaccard_dist <- vegdist(gravity_sp_info, method = "jaccard")

# defining trajectories
gravity_srs_trajectory <- defineTrajectories(gravity_jaccard_dist, sites = gravity_patch_info$unique_id, surveys = gravity_patch_info$time)

# segment lengths of trajectories between consectutive years
gravity_segment_lengths <- trajectoryLengths(gravity_srs_trajectory)
gravity_segment_lengths <- gravity_segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)


# creating time info to join with segment lengths - some surveys were not consecutive years
gravity_time_surveys <- gravity_patch_info %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), remove = F) %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) %>% # removing first survey for sites created in 2007
  dplyr::select(!c("block", "patch", "patch_type"))


# joining with segment lengths
gravity_segment_lengths <- cbind(gravity_segment_lengths, gravity_time_surveys)
gravity_segment_lengths$dispersal_mode <- "Gravity"
gravity_segment_lengths$s.time <- as.numeric(scale(as.numeric(gravity_segment_lengths$time))) # scaling time
gravity_segment_lengths$patch_type <- as.factor(gravity_segment_lengths$patch_type)


###### MODELS ######
# segment lengths model
# linear gravity dispersed
m_length_gravity <- glmmTMB(distance ~ patch_type * s.time + (1|block/patch),
                           data = gravity_segment_lengths)
# quadratic gravity dispersed
m_length_gravity_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                                data = gravity_segment_lengths)
# null gravity dispersed
m_length_gravity_null <- glmmTMB(distance ~ 1 + (1|block/patch),
                                data = gravity_segment_lengths)

# AIC comparison
length.aic.gravity <- list(m_length_gravity, m_length_gravity_quad, m_length_gravity_null)
length.aic.gravity.table <- aictab(length.aic.gravity) # quadratic is best fit
length.aic.gravity.table

# model fit
summary(m_length_gravity_quad)
plot(simulateResiduals(m_length_gravity_quad))
#check_model(m_length_gravity_quad)
performance::r2(m_length_gravity_quad)
anova.gravity.length <- Anova(m_length_gravity_quad, type = "III")


# posthoc
m_length_gravity_posthoc <- emmeans(m_length_gravity_quad, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0)))
m_length_gravity_pairs <- pairs(m_length_gravity_posthoc, simple = "patch_type")
m_length_gravity_pairs





#####################
#### WIND ####
#####################

##### wind CTA segments #####
## wind dispersed 
wind_data <- srs_data %>%
  filter(dispersal_mode == "Wind") %>%
  dplyr::count(unique_id, time, year, sppcode) %>%
  pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) # wide format

# make factor
wind_data$unique_id <- as.factor(wind_data$unique_id)

# patch data
wind_patch_info <- wind_data %>% 
  arrange(unique_id, time) %>%
  dplyr::select(unique_id, time, year)

# species matrix
wind_sp_info <- wind_data %>%
  arrange(unique_id, time) %>%
  mutate(unique_id_year = paste(unique_id, time, year, sep = "-")) %>%
  column_to_rownames("unique_id_year") %>%
  dplyr::select(!c("unique_id", "time", "year"))

# Jaccard distance matrix
wind_jaccard_dist <- vegdist(wind_sp_info, method = "jaccard")

# defining trajectories
wind_srs_trajectory <- defineTrajectories(wind_jaccard_dist, sites = wind_patch_info$unique_id, surveys = wind_patch_info$time)

# segment lengths of trajectories between consectutive years
wind_segment_lengths <- trajectoryLengths(wind_srs_trajectory)
wind_segment_lengths <- wind_segment_lengths %>%
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch", "patch_type")) %>%
  pivot_longer(cols = S1:S22, names_to = "time", values_to = "distance") %>%
  mutate(time = as.numeric(sub("S", "", time))) %>%
  filter(!is.na(distance)) %>%
  dplyr::select(!time)


# creating time info to join with segment lengths - some surveys were not consecutive years
wind_time_surveys <- wind_patch_info %>%
  separate(unique_id, into = c("block", "patch", "patch_type"), remove = F) %>%
  filter(year != 2001) %>% # removing first year for sites created in 2000
  filter(time!= 0) %>% # removing first survey for sites created in 2007
  dplyr::select(!c("block", "patch", "patch_type"))


# joining with segment lengths
wind_segment_lengths <- cbind(wind_segment_lengths, wind_time_surveys)
wind_segment_lengths$dispersal_mode <- "Wind"
wind_segment_lengths$s.time <- as.numeric(scale(as.numeric(wind_segment_lengths$time))) # scaling time
wind_segment_lengths$patch_type <- as.factor(wind_segment_lengths$patch_type)


###### MODELS ######
# segment lengths model 
# linear wind dispersed
m_length_wind <- glmmTMB(distance ~ patch_type * s.time + (1|block/patch),
                            data = wind_segment_lengths)
# quadratic wind dispersed
m_length_wind_quad <- glmmTMB(distance ~ patch_type * s.time + patch_type * I(s.time^2) + (1|block/patch),
                                 data = wind_segment_lengths)
# null wind dispersed
m_length_wind_null <- glmmTMB(distance ~ 1 + (1|block/patch),
                                 data = wind_segment_lengths)

# AIC comparison
length.aic.wind <- list(m_length_wind, m_length_wind_quad, m_length_wind_null)
length.aic.wind.table <- aictab(length.aic.wind) # quadratic is best fit
length.aic.wind.table

# model fit
summary(m_length_wind_quad)
plot(simulateResiduals(m_length_wind_quad))
#check_model(m_length_wind_quad)
performance::r2(m_length_wind_quad)
anova.wind.length <- Anova(m_length_wind_quad, type = "III")


# posthoc
m_length_wind_posthoc <- emmeans(m_length_wind_quad, ~ patch_type*s.time + patch_type * I(s.time^2), at = list(s.time = c(0)))
m_length_wind_pairs <- pairs(m_length_wind_posthoc, simple = "patch_type")
m_length_wind_pairs




########################
#### TABLES ####
########################
# AICc table
model.names <- tibble(Modnames = c("Mod1", "Mod2", "Mod3"),
                      Model = c("Linear", "Quadratic", "Null"),
                      `Model Formula` = c(" ", " ", " ")) # filling in manually after

length.aic.table.df <- as.data.frame(length.aic.table)
length.aic.table.df$`Dispersal Mode` <- "All Species"

length.aic.animal.table.df <- as.data.frame(length.aic.animal.table)
length.aic.animal.table.df$`Dispersal Mode` <- "Animal-Dispersed"

length.aic.gravity.table.df <- as.data.frame(length.aic.gravity.table)
length.aic.gravity.table.df$`Dispersal Mode` <- "Gravity-Dispersed"

length.aic.wind.table.df <- as.data.frame(length.aic.wind.table)
length.aic.wind.table.df$`Dispersal Mode` <- "Wind-Dispersed"


# all together
length.aic.table.all <- rbind(
  length.aic.table.df, length.aic.animal.table.df, length.aic.gravity.table.df, length.aic.wind.table.df
)

length.aic.table.all <- length.aic.table.all %>%
  left_join(model.names, by = "Modnames") %>%
  dplyr::select(`Dispersal Mode`, Model, `Model Formula`, K, LL, AICc, Delta_AICc, Cum.Wt) %>%
  rename(`Cumulative Weight` = Cum.Wt, `Delta AICc` = Delta_AICc)

tableS1 <- length.aic.table.all %>% 
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(length.aic.table.all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(length.aic.table.all), extra_css = "padding-bottom: 5px;")
tableS1

# exporting
# save_kable(tableS1, file = file.path("tables", "tableS1.html"))







# anova table
anova.length_df <- as.data.frame(anova.length)
anova.length_df <- anova.length_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "All Species") %>%
  mutate(`Top Model` = "Quadratic")

anova.animal.length_df <- as.data.frame(anova.animal.length)
anova.animal.length_df <- anova.animal.length_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") %>%
  mutate(`Top Model` = "Null")


anova.gravity.length_df <- as.data.frame(anova.gravity.length)
anova.gravity.length_df <- anova.gravity.length_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed")%>%
  mutate(`Top Model` = "Quadratic")


anova.wind.length_df <- as.data.frame(anova.wind.length)
anova.wind.length_df <- anova.wind.length_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") %>%
  mutate(`Top Model` = " Quadratic ")

# putting all together
m.length_anova_all <- rbind(
  anova.length_df, anova.animal.length_df, anova.gravity.length_df, anova.wind.length_df
)

rename_variable_anova <- tibble(model_term = c("(Intercept)", "patch_type", "s.time", "I(s.time^2)", "patch_type:s.time", "patch_type:I(s.time^2)"),
                                Variable = c("Intercept", "Patch Type", "Time", "Time^2", "Patch Type:Time", "Patch Type:Time^2"))

m.length_anova_all <- m.length_anova_all %>%
  #filter(model_term != "(Intercept)") %>%
  left_join(rename_variable_anova, by = "model_term") %>%
  dplyr::select(`Dispersal mode`, `Top Model`, Variable, Chisq, Df, `Pr(>Chisq)`) %>%
  rename(p.value = `Pr(>Chisq)`, df = Df)

tableS2 <- m.length_anova_all %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.length_anova_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.length_anova_all), extra_css = "padding-bottom: 5px;")
tableS2

# exporting
#save_kable(tableS2, file = file.path("tables", "tableS2.html"))


# emmeans posthoc tables
# creating dataframes of results
# all species
m_length_pairs_df <- as.data.frame(m_length_pairs)
m_length_pairs_df <- m_length_pairs_df %>%
  mutate(`Dispersal mode` = "All Species") 
 
# gravity dispersed
m_length_gravity_pairs_df <- as.data.frame(m_length_gravity_pairs)
m_length_gravity_pairs_df <- m_length_gravity_pairs_df %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed") 

# wind dispersed
m_length_wind_pairs_df <- as.data.frame(m_length_wind_pairs)
m_length_wind_pairs_df <- m_length_wind_pairs_df %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") 

# putting all together
m_length_table_all <- rbind(
  m_length_pairs_df, m_length_gravity_pairs_df, m_length_wind_pairs_df
)

tableS3 <- m_length_table_all %>% 
  dplyr::select(`Dispersal mode`, contrast, estimate, SE, df, z.ratio, p.value) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2)) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_length_table_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_length_table_all), extra_css = "padding-bottom: 5px;")
tableS3

# exporting
# save_kable(tableS3, file = file.path("tables", "tableS3.html"))




####################
#### PLOTS ####
###################

# model predictions
# All species predictions
m_length_predict <- ggpredict(m_length_quad, terms = c("s.time [all]", "patch_type"))
m_length_predict <- as.data.frame(m_length_predict)
m_length_predict$dispersal_mode <- "All Species"
m_length_predict$linetype <- "solid" # adding line type
scaled_time_key <- segment_lengths %>% # creating key of scaled times to join to predictions for easy visualization
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m_length_predict <- m_length_predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))

# Animal dispersed predictions - generating predictions from quadratic model for plotting purposes only - showing trends as dashed line to emphasize that null model was top model
m.animal_segments.predict <- ggpredict(m_length_animal_quad, terms=c("s.time [all]", "patch_type [all]"), back_transform = T)
m.animal_segments.predict <- as.data.frame(m.animal_segments.predict)
m.animal_segments.predict$dispersal_mode <- "Animal"
m.animal_segments.predict$linetype <- "dashed" # adding line type
animal_time_key <- animal_segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m.animal_segments.predict <- m.animal_segments.predict %>%
  left_join(animal_time_key, by = c("x" = "s.time"))

# gravity dispersed predictions
m.gravity_segments.predict <- ggpredict(m_length_gravity_quad, terms=c("s.time [all]", "patch_type [all]"), back_transform = T)
m.gravity_segments.predict <- as.data.frame(m.gravity_segments.predict)
m.gravity_segments.predict$dispersal_mode <- "Gravity"
m.gravity_segments.predict$linetype <- "solid" # adding line type
gravity_time_key <- gravity_segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m.gravity_segments.predict <- m.gravity_segments.predict %>%
  left_join(gravity_time_key, by = c("x" = "s.time"))

# wind dispersed predictions
m.wind_segments.predict <- ggpredict(m_length_wind_quad, terms=c("s.time [all]", "patch_type [all]"), back_transform = T)
m.wind_segments.predict <- as.data.frame(m.wind_segments.predict)
m.wind_segments.predict$dispersal_mode <- "Wind"
m.wind_segments.predict$linetype <- "solid" # adding line type
wind_time_key <- wind_segment_lengths %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
m.wind_segments.predict <- m.wind_segments.predict %>%
  left_join(wind_time_key, by = c("x" = "s.time"))


# FACET BY ROWS - total and animal together and gravity and wind together
# joining together predictions
predict_segments_1 <- rbind(
  m_length_predict, m.animal_segments.predict
)

predict_segments_2 <- rbind(
  m.wind_segments.predict, m.gravity_segments.predict
)

predict_segments_1$dispersal_mode <- factor(predict_segments_1$dispersal_mode, levels = c("All Species", "Animal"))
predict_segments_2$dispersal_mode <- factor(predict_segments_2$dispersal_mode, levels = c("Wind", "Gravity"))

predict_segments_1$time <- as.numeric(as.character(predict_segments_1$time))
predict_segments_2$time <- as.numeric(as.character(predict_segments_2$time))

# joining together data points
# segment_lengths <- segment_lengths %>%
#   dplyr::select(-soil_moisture, -year_since_fire)
# putting data together for plotting
dispersal_mode_segments_1 <- rbind(
  segment_lengths, animal_segment_lengths
)
dispersal_mode_segments_1$dispersal_mode <- factor(dispersal_mode_segments_1$dispersal_mode, levels = c("All Species", "Animal"))
dispersal_mode_segments_1$time <- as.numeric(as.character(dispersal_mode_segments_1$time))

dispersal_mode_segments_2 <- rbind(
  wind_segment_lengths, gravity_segment_lengths
)
dispersal_mode_segments_2$dispersal_mode <- factor(dispersal_mode_segments_2$dispersal_mode, levels = c("Wind", "Gravity"))
dispersal_mode_segments_2$time <- as.numeric(as.character(dispersal_mode_segments_2$time))

# two faceted plots
# first set of plots
segments_plot_1 <- predict_segments_1 %>%
  ggplot() +
  #geom_point(aes(time, distance, color = patch_type), size = 3, alpha = 0.11, data = dispersal_mode_segments_1) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(time, predicted, color = group, linetype = linetype), linewidth = 3) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("All Species" = "(A) All species", "Animal" = "(B) Animal-dispersed"))) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_linetype_manual(values = c('solid','longdash'), guide = "none") +
  xlab(NULL) +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  scale_y_continuous(limits = c(0.16, 0.35), breaks = c(0.20, 0.30), labels = label_number(accuracy = 0.01)) +
  theme(axis.text = element_text(size = 18)) +
  theme(legend.position = "none") 
segments_plot_1

# second set of plots
segments_plot_2 <- predict_segments_2 %>%
  ggplot() +
  #geom_point(aes(time, distance, color = patch_type), size = 3, alpha = 0.11, data = dispersal_mode_segments_2) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("Gravity" = "(C) Gravity-dispersed", "Wind" = "(D) Wind-dispersed"))) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  xlab("Years since site creation") +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  scale_y_continuous(limits = c(0.15, 0.45), labels = label_number(accuracy = 0.01)) +
  theme(axis.text = element_text(size = 18)) +
  theme(legend.position = "none") 
segments_plot_2

# get legend
pL <- predict_segments_2 %>%
  ggplot() +
  #geom_point(aes(time, distance, color = patch_type), size = 4, alpha = 0.3, data = dispersal_mode_segments_2) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.4) +
  geom_line(aes(time, predicted, color = group), linewidth = 2.5) +
  theme_minimal(base_size = 26) +
  theme(legend.position = "bottom") +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
l <- get_legend(pL)

# put together
figure2 <- cowplot::plot_grid(segments_plot_1, segments_plot_2, l,
                                      ncol = 1, nrow = 3, rel_heights = c(1, 1.1, 0.15),
                                      label_size = 20, label_x = 0.2, label_y = 0.95)
figure2


# exporting
# pdf(file = file.path("plots", "figure2.pdf"), width = 12.5, height = 13)
# figure2
# dev.off()



# individual plot
all_segment_plot <- m_length_predict %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 6, alpha = 0.15, data = segment_lengths) +
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
all_segment_plot


pdf(file = file.path("plots", "all_segment_plot.pdf"), width = 12, height = 12)
all_segment_plot
dev.off()



all_segment_plot <- segment_lengths %>%
  ggplot(aes(time, distance, color = patch_type, fill = patch_type)) +
  geom_point(aes(time, distance, color = patch_type), size = 6, alpha = 0.15, data = segment_lengths) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2, linewidth = 4) +
  #geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
 # geom_line(aes(time, predicted, color = group), linewidth = 3.5) +
  theme_minimal(base_size = 36) +
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
  scale_y_continuous(limits = c(0.15, 0.45), labels = label_number(accuracy = 0.01)) +
  theme(axis.text = element_text(size = 20), 
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(legend.position = "top")
all_segment_plot

pdf(file = file.path("plots", "all_segment_plot.pdf"), width = 12, height = 12)
all_segment_plot
dev.off()



# dispersal mode plots
# FACET BY ROWS - total and animal together and gravity and wind together
# joining together predictions
predict_segments_mode <- rbind(
  m.animal_segments.predict, m.wind_segments.predict, m.gravity_segments.predict
)


predict_segments_mode$dispersal_mode <- factor(predict_segments_mode$dispersal_mode, levels = c("Animal", "Wind", "Gravity"))
predict_segments_mode$linetype <- factor(predict_segments_mode$linetype, levels = c("solid", "dashed"))

predict_segments_mode$time <- as.numeric(as.character(predict_segments_mode$time))

# joining together data points
# segment_lengths <- segment_lengths %>%
#   dplyr::select(-soil_moisture, -year_since_fire)
# putting data together for plotting
dispersal_mode_segments_points <- rbind(
  animal_segment_lengths, wind_segment_lengths, gravity_segment_lengths
)
dispersal_mode_segments_points$dispersal_mode <- factor(dispersal_mode_segments_points$dispersal_mode, levels = c("Animal", "Wind", "Gravity"))
dispersal_mode_segments_points$time <- as.numeric(as.character(dispersal_mode_segments_points$time))

# two faceted plots
# first set of plots
segments_plot_mode <- predict_segments_mode %>%
  ggplot() +
  geom_point(aes(time, distance, color = patch_type), size = 3, alpha = 0.1, data = dispersal_mode_segments_points) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(time, predicted, color = group, linetype = linetype), linewidth = 2.5) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("Animal" = "(A) Animal-dispersed", "Wind" = "(B) Wind-dispersed", "Gravity" = "(C) Gravity-dispersed"))) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_linetype_manual(values = c('solid','longdash'), guide = "none") +
  xlab(NULL) +
  ylab(expression(atop("Trajectory distance", paste("between consecutive surveys")))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
 # scale_y_continuous(limits = c(0.16, 0.35), breaks = c(0.20, 0.30), labels = label_number(accuracy = 0.01)) +
  theme(axis.text = element_text(size = 18)) +
  theme(legend.position = "none") 
segments_plot_mode


pdf(file = file.path("plots", "segments_plot_mode_plot.pdf"), width = 14, height = 6)
segments_plot_mode
dev.off()

