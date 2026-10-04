########
## SCRIPT NAME: 03_convergence.R
## AUTHOR: Katherine Hulting
## PURPOSE: Calculate convergence/divergence between patch type communities across time, repeat within dispersal mode groups
## PRODUCTS: figure2.pdf, tableS1.html, tableS2.html, tableS3.html
#########

librarian::shelf(tidyverse, vegan, ape, BiodiversityR, glmmTMB, AICcmodavg, 
                 DHARMa, emmeans, car, ggeffects, performance, cowplot, here, kableExtra, scales)

source(here::here(file.path("scripts", "02_pcoa.R")))
source(here::here(file.path("scripts", "00_functions.R")))

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))


#########################
#### ALL SPECIES ####
#########################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
convergence_jaccard <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
convergence_jaccard <- convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "All Species")
convergence_jaccard$s.time <- as.numeric(scale(convergence_jaccard$time)) # scaling time
convergence_jaccard$patch_pair <- as.factor(convergence_jaccard$patch_pair)

##### Models #####
# linear model
m.converge <- glmmTMB(jaccard ~ patch_pair * s.time + (1|block),
              data = convergence_jaccard)

# quadratic model
m.converge_quad <- glmmTMB(jaccard ~ patch_pair * s.time + patch_pair * I(s.time^2) + (1|block),
                         data = convergence_jaccard)
# null model
m.converge_null <- glmmTMB(jaccard ~ 1 + (1|block), # null model
                         data = convergence_jaccard)
# AIC comparison
converge.aic <- list(m.converge, m.converge_quad, m.converge_null)
converge.aic.table <- aictab(converge.aic) # quadratic much better fit
converge.aic.table

## model checking
summary(m.converge_quad)




# percent change in dissimilarity from year 1-21 (20 years)
# time 1 = -1.519113
# time 21 = 1.436766
# intercept + time_estimate (time) + time^2_estimate (time)
# at time 1
0.384906 + 0.013563*(-1.519113) + 0.005123*((-1.519113)^2) # 0.3761246
# at time 21
0.384906 + 0.013563*(1.436766) + 0.005123*((1.436766)^2) # 0.4149682

# percent change = time 21 - time 1 / time 1 * 100
(0.4149682 - 0.3761246)/0.3761246 * 100 #10.32732 % increase

# 95% CI, percent change
confint(m.converge_quad)
# lower 95% CI
# at time 1
0.3607879039 + 0.0081246938*(-1.519113) + -0.0003483755*((-1.519113)^2) # 0.3476416
# at time 21
0.3607879039 + 0.0081246938*(1.436766) + -0.0003483755*((1.436766)^2) # 0.371742
# percent change = time 21 - time 1 / time 1 * 100
(0.371742 - 0.3476416)/0.3476416 * 100 #6.932542 % increase

# upper 95% CI
# at time 1
0.4090249487 + 0.0190010704*(-1.519113) + 0.0105941226*((-1.519113)^2) # 0.4046083
# at time 21
0.4090249487 + 0.0190010704*(1.436766) + 0.0105941226*((1.436766)^2) # 0.4581945
# percent change = time 21 - time 1 / time 1 * 100
(0.4046083 - 0.4581945)/0.4581945 * 100 #11.69508 % increase


# emmeans(m.converge_quad, ~s.time+I(s.time^2), at = list(s.time = c(-1.5320869, 1.4326678)),
#         type = "response")
# (0.382 - 0.398) / 0.398 * 100 # 8.28877 % increase over 20 years

# Anova
anova.converge <- Anova(m.converge_quad, type = "III")
# model checking
plot(simulateResiduals(m.converge_quad))
#check_model(m.converge_quad)
performance::r2(m.converge_quad)


## posthoc comparisons
m.converge_posthoc <- emmeans(m.converge_quad, ~ patch_pair*s.time+ patch_pair * I(s.time^2))
m.converge_pairs <- pairs(m.converge_posthoc, simple = "patch_pair")
m.converge_pairs

m.converge_posthoc
# % increase in dissimilarity from (connected-winged) to (connected-rectangular)
(0.385-0.359)/0.359 * 100 # connected patches are %7.24234 more similar to winged patches than rectangular patches across time
#95% CI
0.02548 +1.96 *0.00541 # upper = 0.0360836
0.0360836/0.359 * 100 # upper = 10.05114%
0.02548 -1.96 *0.00541 # lower = 0.0148764
0.0148764/0.359 *100 # lower = 4.143844%

# winged patches are closer in composition to connected patches than rectangular patches
(0.359-0.375)/0.375 *100 # winged patches are 4.266667 % less similar to rectangular patches and closer to connected patches
-0.01590 +1.96 *0.00499 # upper = -0.0061196
-0.0061196/0.375 *100 #-1.631893
-0.01590 -1.96 *0.00499 # upper = -0.0256804
-0.0256804/0.375 *100 #-6.848107






###################################
#### ANIMAL DISPERSED SPECIES ####
###################################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
animal_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Animal") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
animal_convergence_jaccard <- animal_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Animal")

# scaling time
animal_convergence_jaccard$s.time <- as.numeric(scale(animal_convergence_jaccard$time))
animal_convergence_jaccard$patch_pair <- as.factor(animal_convergence_jaccard$patch_pair)

##### Models #####
# animal dispersed linear
m.converge_animal_linear <- glmmTMB(jaccard ~ patch_pair*s.time + (1|block),
                                  data = animal_convergence_jaccard)

# animal dispersed quadratic
m.converge_animal_quad <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                             data = animal_convergence_jaccard)

# animal dispersed null
m.converge_animal_null <- glmmTMB(jaccard ~ 1 + (1|block),
                                    data = animal_convergence_jaccard)

# AIC comparison
converge.aic.animal <- list(m.converge_animal_linear, m.converge_animal_quad, m.converge_animal_null)
converge.aic.animal.table <- aictab(converge.aic.animal) # linear within 2 delta AICc and more parsimonious model
converge.aic.animal.table

# model checking
summary(m.converge_animal_linear)
plot(simulateResiduals(m.converge_animal_linear))
#check_model(m.converge_animal_linear)
anova.animal.converge <- Anova(m.converge_animal_linear, type = "III")

# posthoc
m.converge_animal_posthoc <- emmeans(m.converge_animal_linear, ~ patch_pair*s.time)
m.converge_animal_pairs <- pairs(m.converge_animal_posthoc, simple = "patch_pair")
m.converge_animal_pairs





###################################
#### GRAVITY DISPERSED SPECIES ####
###################################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
gravity_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Gravity") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
gravity_convergence_jaccard <- gravity_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Gravity")


# scaling time
gravity_convergence_jaccard$s.time <- as.numeric(scale(gravity_convergence_jaccard$time))
gravity_convergence_jaccard$patch_pair <- as.factor(gravity_convergence_jaccard$patch_pair)

##### Models #####
# gravity dispersed linear
m.converge_gravity_linear <- glmmTMB(jaccard ~ patch_pair*s.time + (1|block),
                                    data = gravity_convergence_jaccard)

# gravity dispersed quadratic
m.converge_gravity_quad <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                                  data = gravity_convergence_jaccard)

# gravity dispersed null
m.converge_gravity_null <- glmmTMB(jaccard ~ 1 + (1|block),
                                  data = gravity_convergence_jaccard)

# AIC comparison
converge.aic.gravity <- list(m.converge_gravity_linear, m.converge_gravity_quad, m.converge_gravity_null)
converge.aic.gravity.table <- aictab(converge.aic.gravity) # quadratic better fit
converge.aic.gravity.table

# model checking
summary(m.converge_gravity_quad)
plot(simulateResiduals(m.converge_gravity_quad))
#check_model(m.converge_gravity_quad)
anova.gravity.converge <- Anova(m.converge_gravity_quad, type = "III")

# posthoc
m.converge_gravity_posthoc <- emmeans(m.converge_gravity_quad, ~ patch_pair*s.time + patch_pair * I(s.time^2))
m.converge_gravity_pairs <- pairs(m.converge_gravity_posthoc, simple = "patch_pair")
m.converge_gravity_pairs






###################################
#### WIND DISPERSED SPECIES ####
###################################

##### Convergence/divergence between patch types based on distance #####
# iterate over blocks, for each patch pair within a block, compute jaccard dissimilarity for each year
# splitting into blocks, applying function, putting back together
wind_convergence_jaccard <- srs_data %>%
  filter(dispersal_mode == "Wind") %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block) %>%
  group_split() %>%
  lapply(compute_convergence_jaccard) %>%
  bind_rows() # putting together into a dataframe

# removing same patch type comparisons and time 0 (only for 52 and 57)
wind_convergence_jaccard <- wind_convergence_jaccard %>%
  filter(!patch_pair %in% c("Rectangular-Rectangular", "Winged-Winged")) %>%
  filter(time != 0) %>%
  mutate(dispersal_mode = "Wind")

# scaling time
wind_convergence_jaccard$s.time <- as.numeric(scale(wind_convergence_jaccard$time))
wind_convergence_jaccard$patch_pair <- as.factor(wind_convergence_jaccard$patch_pair)

##### Models #####
# wind dispersed linear
m.converge_wind_linear <- glmmTMB(jaccard ~ patch_pair*s.time + (1|block),
                                     data = wind_convergence_jaccard)

# wind dispersed quadratic
m.converge_wind_quad <- glmmTMB(jaccard ~ patch_pair*s.time + patch_pair*I(s.time^2) + (1|block),
                                   data = wind_convergence_jaccard)

# wind dispersed null
m.converge_wind_null <- glmmTMB(jaccard ~ 1 + (1|block),
                                   data = wind_convergence_jaccard)

# AIC comparison
converge.aic.wind <- list(m.converge_wind_linear, m.converge_wind_quad, m.converge_wind_null)
converge.aic.wind.table <- aictab(converge.aic.wind) # quadratic better fit
converge.aic.wind.table

# model checking
summary(m.converge_wind_quad)
plot(simulateResiduals(m.converge_wind_quad))
#check_model(m.converge_wind_quad)
anova.wind.converge <- Anova(m.converge_wind_quad, type = "III")

# posthoc
m.converge_wind_posthoc <- emmeans(m.converge_wind_quad, ~ patch_pair*s.time + patch_pair * I(s.time^2))
m.converge_wind_pairs <- pairs(m.converge_wind_posthoc, simple = "patch_pair")
m.converge_wind_pairs




#########################
#### TABLES ####
#########################
# AICc table
model.names <- tibble(Modnames = c("Mod1", "Mod2", "Mod3"),
                      Model = c("Linear", "Quadratic", "Null"),
                      `Model Formula` = c(" ", " ", " ")) # filling in manually after

converge.aic.table.df <- as.data.frame(converge.aic.table)
converge.aic.table.df$`Dispersal Mode` <- "All Species"

converge.aic.animal.table.df <- as.data.frame(converge.aic.animal.table)
converge.aic.animal.table.df$`Dispersal Mode` <- "Animal-Dispersed"

converge.aic.gravity.table.df <- as.data.frame(converge.aic.gravity.table)
converge.aic.gravity.table.df$`Dispersal Mode` <- "Gravity-Dispersed"

converge.aic.wind.table.df <- as.data.frame(converge.aic.wind.table)
converge.aic.wind.table.df$`Dispersal Mode` <- "Wind-Dispersed"

# all together
converge.aic.table.all <- rbind(
  converge.aic.table.df, converge.aic.animal.table.df, converge.aic.gravity.table.df, converge.aic.wind.table.df
)

converge.aic.table.all <- converge.aic.table.all %>%
  left_join(model.names, by = "Modnames") %>%
  dplyr::select(`Dispersal Mode`, Model, `Model Formula`, K, LL, AICc, Delta_AICc, Cum.Wt) %>%
  rename(`Cumulative Weight` = Cum.Wt, `Delta AICc` = Delta_AICc)

tableS10 <- converge.aic.table.all %>% 
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(converge.aic.table.all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(converge.aic.table.all), extra_css = "padding-bottom: 5px;")
tableS10

# exporting
# save_kable(tableS10, file = file.path("tables", "tableS10.html"))





# Anova table 
anova.converge_df <- as.data.frame(anova.converge)
anova.converge_df <- anova.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "All Species") %>%
  mutate(`Top Model` = "Quadratic")

anova.animal.converge_df <- as.data.frame(anova.animal.converge)
anova.animal.converge_df <- anova.animal.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") %>%
  mutate(`Top Model` = "Linear")

anova.gravity.converge_df <- as.data.frame(anova.gravity.converge)
anova.gravity.converge_df <- anova.gravity.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed") %>%
  mutate(`Top Model` = "Quadratic")

anova.wind.converge_df <- as.data.frame(anova.wind.converge)
anova.wind.converge_df <- anova.wind.converge_df %>%
  rownames_to_column("model_term") %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") %>%
  mutate(`Top Model` = "Quadratic")

m.converge_anova_all <- rbind(
  anova.converge_df, anova.animal.converge_df, anova.gravity.converge_df, anova.wind.converge_df
)

rename_variable_anova <- tibble(model_term = c("(Intercept)", "patch_pair", "s.time", "I(s.time^2)", "patch_pair:s.time", "patch_pair:I(s.time^2)"),
                                Variable = c("Intercept", "Patch Pair", "Time", "Time^2", "Patch Pair:Time", "Patch Pair:Time^2"))


m.converge_anova_all <- m.converge_anova_all %>%
  #filter(model_term != "(Intercept)") %>%
  left_join(rename_variable_anova, by = "model_term") %>%
  dplyr::select(`Dispersal mode`, `Top Model`, Variable, Chisq, Df, `Pr(>Chisq)`) %>%
  rename(p.value = `Pr(>Chisq)`, df = Df)

tableS11 <- m.converge_anova_all %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.converge_anova_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.converge_anova_all), extra_css = "padding-bottom: 5px;")
tableS11

# exporting
# save_kable(tableS11, file = file.path("tables", "tableS11.html"))


# emmeans posthoc tables
# creating dataframes of results
# all species
m.converge_pairs_df <- as.data.frame(m.converge_pairs)
m.converge_pairs_df <- m.converge_pairs_df %>%
  mutate(`Dispersal mode` = "All Species") 

# animal dispersed
m.converge_animal_pairs_df <- as.data.frame(m.converge_animal_pairs)
m.converge_animal_pairs_df <- m.converge_animal_pairs_df %>%
  mutate(`Dispersal mode` = "Animal-Dispersed") 

# gravity dispersed
m.converge_gravity_pairs_df <- as.data.frame(m.converge_gravity_pairs)
m.converge_gravity_pairs_df <- m.converge_gravity_pairs_df %>%
  mutate(`Dispersal mode` = "Gravity-Dispersed")

# wind dispersed
m.converge_wind_pairs_df <- as.data.frame(m.converge_wind_pairs)
m.converge_wind_pairs_df <- m.converge_wind_pairs_df %>%
  mutate(`Dispersal mode` = "Wind-Dispersed") 

# putting all together
m.converge_emmeans_all <- rbind(
  m.converge_pairs_df, m.converge_animal_pairs_df, m.converge_gravity_pairs_df, m.converge_wind_pairs_df
)

# table S3
tableS12 <- m.converge_emmeans_all %>% 
  dplyr::select(`Dispersal mode`, contrast, estimate, SE, df, z.ratio, p.value) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m.converge_emmeans_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m.converge_emmeans_all), extra_css = "padding-bottom: 5px;")
tableS12

# exporting
# save_kable(tableS12, file = file.path("tables", "tableS12.html"))




#########################
#### PLOTS ####
#########################
## All species dispersed plot
# model predictions
m.converge.predict <- ggpredict(m.converge_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge.predict <- as.data.frame(m.converge.predict)
m.converge.predict$dispersal_mode <- "All Species"

## animal dispersed plot
# model predictions
m.converge_animal.predict <- ggpredict(m.converge_animal_linear, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_animal.predict <- as.data.frame(m.converge_animal.predict)
m.converge_animal.predict$dispersal_mode <- "Animal"

## gravity dispersed plot
# model predictions
m.converge_gravity.predict <- ggpredict(m.converge_gravity_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_gravity.predict <- as.data.frame(m.converge_gravity.predict)
m.converge_gravity.predict$dispersal_mode <- "Gravity"

## wind dispersed plot
# model predictions
m.converge_wind.predict <- ggpredict(m.converge_wind_quad, terms=c("s.time [all]", "patch_pair [all]"), back_transform = T)
m.converge_wind.predict <- as.data.frame(m.converge_wind.predict)
m.converge_wind.predict$dispersal_mode <- "Wind"


# FACET BY ROWS - total and animal together and gravity and wind together
# joining together predictions
predict_converge_1 <- rbind(
  m.converge.predict, m.converge_animal.predict
)
predict_converge_2 <- rbind(
  m.converge_wind.predict, m.converge_gravity.predict
)


# creating key of scaled times to join to predictions for easy visualization
scaled_time_key <- convergence_jaccard %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))

# joining with time
predict_converge_1 <- predict_converge_1 %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))
predict_converge_2 <- predict_converge_2 %>%
  left_join(scaled_time_key, by = c("x" = "s.time"))

predict_converge_1$dispersal_mode <- factor(predict_converge_1$dispersal_mode, levels = c("All Species", "Animal"))
predict_converge_2$dispersal_mode <- factor(predict_converge_2$dispersal_mode, levels = c("Wind", "Gravity"))

# joining together data points
dispersal_mode_convergence_1 <- rbind(
  convergence_jaccard, animal_convergence_jaccard
)
dispersal_mode_convergence_1$dispersal_mode <- factor(dispersal_mode_convergence_1$dispersal_mode, levels = c("All Species", "Animal"))

dispersal_mode_convergence_2 <- rbind(
  wind_convergence_jaccard, gravity_convergence_jaccard
)
dispersal_mode_convergence_2$dispersal_mode <- factor(dispersal_mode_convergence_2$dispersal_mode, levels = c("Wind", "Gravity"))

# first set of plots
converge_plot_1 <- predict_converge_1 %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3, alpha = 0.05, data = dispersal_mode_convergence_1) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("All Species" = "(A) All species", "Animal" = "(B) Animal-dispersed"))) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(color = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab(NULL) +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  guides(fill=guide_legend(ncol=1)) +
  guides(color=guide_legend(ncol=1)) +
  #ylim(0.25, 0.55) +
  # scale_y_continuous(limits = c(0.25, 0.47), labels = label_number(accuracy = 0.01)) +
  theme(axis.text = element_text(size = 16)) +
  theme(legend.position = "none") 
converge_plot_1

# second set of plots
converge_plot_2 <- predict_converge_2 %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 3, alpha = 0.05, data = dispersal_mode_convergence_2) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  geom_line(aes(time, predicted, color = group), linewidth = 3) +
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("Gravity" = "(C) Gravity-dispersed", "Wind" = "(D) Wind-dispersed"))) +
  theme_minimal(base_size = 26) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Comparison") +
  xlab("Years since site creation") +
  ylab(expression(paste("Spatial ", beta, " diversity (Jaccard)"))) +
  guides(fill=guide_legend(ncol=1)) +
  #ylim(0.25, 0.55) +
  # scale_y_continuous(limits = c(0.32, 0.55), labels = label_number(accuracy = 0.01)) +
  guides(color=guide_legend(ncol=1)) +
  theme(axis.text = element_text(size = 16)) +
  theme(legend.position = "none") 
converge_plot_2

# get legend
pL <- m.converge_wind.predict %>%
  left_join(scaled_time_key, by = c("x" = "s.time")) %>%
  ggplot() +
  geom_point(aes(time, jaccard, color = patch_pair), size = 4, alpha = 0.05, data = wind_convergence_jaccard) +
  geom_ribbon(aes(x = time, ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.5) +
  geom_line(aes(time, predicted, color = group), linewidth = 1.5) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "bottom") +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), labels = c(expression("Connected"%<->%"Rectangular"), 
                                                                            expression("Connected"%<->%"Winged"),
                                                                            expression("Rectangular"%<->%"Winged")), name = "Patch Comparison") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), labels = c(expression("Connected"%<->%"Rectangular"), 
                                                                             expression("Connected"%<->%"Winged"),
                                                                             expression("Rectangular"%<->%"Winged")), name = "Patch Comparison")
l <- get_legend(pL)

# put together
figure4 <- cowplot::plot_grid(converge_plot_1, converge_plot_2, 
                                      ncol = 1, nrow = 2, rel_heights = c(1, 1.1),
                                      label_size = 20, label_x = 0.2, label_y = 0.95)
figure4

# # exporting
# pdf(file = file.path("plots", "figure4.pdf"), width = 12.5, height = 12.5)
# figure4
# dev.off()
# #
# # #exporting legend seperately
# pdf(file = file.path("plots", "figure4_legend.pdf"), width = 11.5, height = 1.5)
# plot(l)
# dev.off()

