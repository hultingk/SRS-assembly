########
## SCRIPT NAME: 05_supplement.R
## AUTHOR: Katherine Hulting
## PURPOSE: Supplemental exploratory plots: species richness over time, proportion of dispersal mode, turnover
## PRODUCTS: figureS2.pdf, figureS3.pdf, figureS6.pdf
#########

# loading libraries
librarian::shelf(tidyverse, vegan, ecotraj, glmmTMB, DHARMa, emmeans, ggeffects, 
                 AICcmodavg, performance, cowplot, mgcViz) # Install missing packages and load needed libraries

source(here::here(file.path("scripts", "00_functions.R"))) # loading functions

# loading data
srs_data <- read_csv(file = file.path("data", "L1_wrangled", "srs_plant_all.csv"))



##### compositional changes #####
# calculating raw changes in composition: # of gains, # of losses, # of species staying the same from one year to the next
species_changes <- srs_data %>%
  count(block, patch, patch_type, unique_id, year, time, sppcode) %>%
  group_by(block, patch) %>%
  group_split() %>%
  lapply(compute_composition_change) %>%
  bind_rows() %>% # putting together into a dataframe
  rownames_to_column("unique_id") %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type", "time2"), sep = "-") %>%
  dplyr::select(-time2) # removing duplicate time column

# removing the first year 
species_changes <- species_changes %>%
  filter(!time %in% c("0")) %>% # removing time 0 - only recorded for 52 and 57, first year of sampling so no turnover from previous year
  filter(!(time == 1 & block %in% c("08", "10", "53N", "53S", "54N", "54S", "75E", "75W"))) %>% # first year of sampling for these blocks so no turnover from previous year
  pivot_longer(5:8, names_to = "type", values_to = "change")


#### modeling # of species that persisted between years ####
stayed_present <- species_changes %>%
  filter(type %in% c("stayed_present"))

stayed_present$s.time <- as.numeric(scale(stayed_present$time))
stayed_present$patch_type <- as.factor(stayed_present$patch_type)
stayed_present$block <- as.factor(stayed_present$block)
stayed_present$patch_rep <- as.factor(stayed_present$patch_rep)

m_present <-gam(change ~ patch_type +
                  s(s.time, by = patch_type, k = 5)  + 
                  s(block, bs = "re") + s(block, patch_rep, bs="re"),
                family = nb(),
                data = stayed_present)

m_present_summary <- summary(m_present)
m_present_summary
gam.check(m_present)
plot(simulateResiduals(m_present))

# posthoc
m_present_posthoc <- emmeans(m_present, ~patch_type)
m_present_posthoc
m_present_pairs <- pairs(m_present_posthoc)
m_present_pairs

#### modeling # of species that changed between years ####
changed_total <- species_changes %>%
  filter(type %in% c("changed_total"))

changed_total$s.time <- as.numeric(scale(changed_total$time))
changed_total$patch_type <- as.factor(changed_total$patch_type)
changed_total$block <- as.factor(changed_total$block)
changed_total$patch_rep <- as.factor(changed_total$patch_rep)

m_change <-gam(change ~ patch_type +
                 s(s.time, by = patch_type, k = 5)  + 
                 s(block, bs = "re") + s(block, patch_rep, bs="re"),
               family = nb(),
               data = changed_total)
m_change_summary <- summary(m_change)
gam.check(m_change)
plot(simulateResiduals(m_change))


# posthoc
m_change_posthoc <- emmeans(m_change, ~patch_type)
m_change_pairs <- pairs(m_change_posthoc)
m_change_pairs


#### model summary table ####
# table for model summary of amount of change 
# dataframe of parametric predictors
m_change_table_patch <- as.data.frame(m_change_summary$p.table) %>%
  tibble::rownames_to_column("Variable") %>%
  rename(
    estimate = Estimate,
    SE = `Std. Error`,
    p.value = `Pr(>|z|)`
  )

# dataframe of smoothing predictors
rename_change_table <- tibble(Term = c("(Intercept)", "s(s.time):patch_typeConnected", "s(s.time):patch_typeRectangular", "s(s.time):patch_typeWinged"),
                              Variable = c("Intercept", "Time:Connected Patch", "Time:Rectangular Patch", "Time:Winged Patch"))

m_change_table_smooth <- as.data.frame(m_change_summary$s.table) %>%
  tibble::rownames_to_column("Term") %>%
  filter(!Term %in% c("s(block)", "s(block,patch_rep)")) %>%
  rename(
    edf = edf,
    p.value = `p-value`
  ) %>%
  left_join(rename_change_table, by = "Term") %>%
  dplyr::select(Variable, edf, Ref.df, Chi.sq, p.value)

# table of parametric predictors
table_change_patch <- m_change_table_patch %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_change_table_patch), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_change_table_patch), extra_css = "padding-bottom: 5px;") %>%
  row_spec(3, extra_css = "border-bottom: 5px double;")
table_change_patch

# exporting
#save_kable(table_change_patch, file = file.path("tables", "table_change_patch.html"))

# table of smoothing predictors
table_change_smooth <- m_change_table_smooth %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_change_table_smooth), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_change_table_smooth), extra_css = "padding-bottom: 5px;")
table_change_smooth

# exporting
#save_kable(table_change_smooth, file = file.path("tables", "table_change_smooth.html"))


# table for model summary of # of persisting species
# dataframe of parametric predictors
m_present_table_patch <- as.data.frame(m_present_summary$p.table) %>%
  tibble::rownames_to_column("Variable") %>%
  rename(
    estimate = Estimate,
    SE = `Std. Error`,
    p.value = `Pr(>|z|)`
  )

# dataframe of smoothing predictors
rename_present_table <- tibble(Term = c("(Intercept)", "s(s.time):patch_typeConnected", "s(s.time):patch_typeRectangular", "s(s.time):patch_typeWinged"),
                              Variable = c("Intercept", "Time:Connected Patch", "Time:Rectangular Patch", "Time:Winged Patch"))

m_present_table_smooth <- as.data.frame(m_present_summary$s.table) %>%
  tibble::rownames_to_column("Term") %>%
  filter(!Term %in% c("s(block)", "s(block,patch_rep)")) %>%
  rename(
    edf = edf,
    p.value = `p-value`
  ) %>%
  left_join(rename_present_table, by = "Term") %>%
  dplyr::select(Variable, edf, Ref.df, Chi.sq, p.value)

# table of parametric predictors
table_present_patch <- m_present_table_patch %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_present_table_patch), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_present_table_patch), extra_css = "padding-bottom: 5px;") %>%
  row_spec(3, extra_css = "border-bottom: 5px double;")
table_present_patch

# exporting
#save_kable(table_present_patch, file = file.path("tables", "table_present_patch.html"))

# table of smoothing predictors
table_present_smooth <- m_present_table_smooth %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2), target = 1) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_present_table_smooth), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_present_table_smooth), extra_css = "padding-bottom: 5px;")
table_present_smooth

# exporting
#save_kable(table_present_smooth, file = file.path("tables", "table_present_smooth.html"))


#### emmeans summary table ####
m_change_pairs_df <- as.data.frame(m_change_pairs)
m_change_pairs_df <- m_change_pairs_df %>%
  mutate(`Response` = "Number of gains and losses") 

m_present_pairs_df <- as.data.frame(m_present_pairs)
m_present_pairs_df <- m_present_pairs_df %>%
  mutate(`Response` = "Number of persisting species") 

# putting all together
m_change_table_all <- rbind(
  m_change_pairs_df, m_present_pairs_df
)

change_table_all <- m_change_table_all %>% 
  dplyr::select(`Response`, contrast, estimate, SE, df, t.ratio, p.value) %>%
  kbl(digits = 3) %>%
  kable_classic(full_width = T) %>%
  kable_styling(html_font = "Times New Roman",
                font_size = 16) %>%
  collapse_rows(columns = c(1, 2)) %>%
  row_spec(0, extra_css = "border-bottom: 5px double;") %>%
  row_spec(1:nrow(m_change_table_all), extra_css = "border-bottom: 1px solid;") %>%
  row_spec(0:nrow(m_change_table_all), extra_css = "padding-bottom: 5px;")
change_table_all

# exporting
#save_kable(change_table_all, file = file.path("tables", "emmeans_change_table_all.html"))



#### plotting Figure S2 ####
# changed total plot
predict_changed <- ggpredict(m_change, terms = c("s.time [all]", "patch_type"))
# adding time as non-scaled
changed_time_key <- changed_total %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
predict_changed <- predict_changed %>%
  left_join(changed_time_key, by = c("x" = "s.time"))

changed_plot <- predict_changed %>%
  ggplot(aes(x = time, y = predicted, color = group)) +
  geom_point(aes(time, change, color = patch_type), data = changed_total, size = 4, alpha = 0.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_line(size = 3) +
  xlab("Years since site creation") +
  ylab(expression(atop("Number of gains and losses", paste("between consecutive surveys")))) +
  theme_minimal(base_size = 22) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text = element_text(size = 16)) +
  ylim(0, 150) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  theme(legend.position = "none")
changed_plot


# stayed present plot
predict_present <- ggpredict(m_present, terms = c("s.time [all]", "patch_type"))
# adding time as non-scaled
present_time_key <- stayed_present %>%
  count(time, s.time) %>%
  dplyr::select(-n) %>%
  mutate(s.time = round(s.time, 2))
predict_present <- predict_present %>%
  left_join(present_time_key, by = c("x" = "s.time"))

present_plot <- predict_present %>%
  ggplot(aes(x = time, y = predicted, color = group)) +
  geom_point(aes(time, change, color = patch_type), data = stayed_present, size = 4, alpha = 0.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_line(size = 3) +
  xlab("Years since site creation") +
  ylab(expression(atop("Number of species consistent", paste("between consecutive surveys")))) +
  theme_minimal(base_size = 22) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(axis.text = element_text(size = 16)) +
  ylim(0, 150) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
present_plot

# all together
figureS2 <- cowplot::plot_grid(changed_plot, present_plot, rel_widths = c(1, 1.44),
                               labels = c("(A)", "(B)"), label_size = 16)
figureS2


# exporting
pdf(file = file.path("plots", "figureS2.pdf"), width = 14, height = 5.5)
figureS2
dev.off()


# individual present plot
present_plot <- predict_present %>%
  ggplot(aes(x = time, y = predicted, color = group)) +
  geom_point(aes(time, change, color = patch_type), data = stayed_present, size = 6, alpha = 0.15) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  geom_line(size = 3.5) +
  xlab("Years since site creation") +
  ylab(expression(atop("Number of species persisting", paste("between consecutive surveys")))) +
  theme_minimal(base_size = 32) +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        panel.grid.major = element_line(linetype = 2, linewidth = 0.7, color = "grey85"), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black", linewidth = 0.7),
        strip.text.x = element_text(hjust = -0.05)) +
  ylim(30, 150) +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  theme(axis.text = element_text(size = 20), 
        legend.text = element_text(size = 26),
        legend.title = element_text(size = 26),
        panel.background = element_rect(fill = "transparent", color = NA), # Inside axes
        plot.background = element_rect(fill = "transparent", color = NA)) +
  theme(legend.position = "none")
present_plot

pdf(file = file.path("plots", "present_plot.pdf"), width = 12.3, height = 11.5)
present_plot
dev.off()




#### species richness over time ####
srs_richness <- srs_data %>%
  dplyr::count(unique_id, time) 

# species richness by dispersal mode
srs_richness_dispersal <- srs_data %>%
  dplyr::count(unique_id, time, dispersal_mode) %>%
  pivot_wider(names_from = dispersal_mode, values_from = n) %>%
  left_join(srs_richness, by = c("unique_id", "time")) %>% # joining with total richness
  rename(`All Species` = n) %>% # renaming the count of total richness as "all species"
  pivot_longer(cols = c("Animal", "Gravity", "Wind", "All Species"), names_to = "dispersal_mode", # pivoting longer to making column of dispersal mode
               values_to = "richness") %>%
  separate(unique_id, into = c("block", "patch_rep", "patch_type"), sep = "-") # seperating unique ID

# ordering factors
srs_richness_dispersal$dispersal_mode <- factor(srs_richness_dispersal$dispersal_mode, levels = c("All Species", "Animal", "Wind", "Gravity"))

## species richness plot
figureS3 <- srs_richness_dispersal %>%
  ggplot(aes(time, richness, color = patch_type, fill = patch_type)) +
  geom_point(alpha = 0.1, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2, linewidth = 2) + # allowing quadratic line
  facet_wrap(~dispersal_mode, scales = "free", labeller = as_labeller(c("All Species" = "(A) All species", "Animal" = "(B) Animal-dispersed", "Gravity" = "(C) Gravity-dispersed", "Wind" = "(D) Wind-dispersed"))) +
  theme_minimal(base_size = 22) +
  theme(axis.text = element_text(size = 18)) +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "darkgrey", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05)) +
  xlab("Years since site creation") +
  ylab("Species richness") +
  scale_fill_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type") +
  scale_color_manual(values = c("#5389A4", "#CC6677", "#DCB254"), name = "Patch Type")
figureS3

# exporting
# pdf(file = file.path("plots", "figureS3.pdf"), width = 10.5, height = 9)
# figureS3
# dev.off()


#### proportion of dispersal mode over time ####
srs_dispersal_prop <- srs_richness_dispersal %>%
  pivot_wider(names_from = dispersal_mode, values_from = richness) %>%
  mutate(prop_animal = Animal/`All Species`, # calculating proportion of each dispersal mode by dividing the sp richnes of that dispersal mode by the total species richness at that time
         prop_gravity = Gravity/`All Species`,
         prop_wind = Wind/`All Species`) %>%
  dplyr::select(-`All Species`, -Animal, -Gravity, -Wind) %>%
  pivot_longer(cols = c("prop_animal", "prop_gravity", "prop_wind"), names_to = "dispersal_mode", values_to = "proportion")

# proportion plot over time
figureS4 <- srs_dispersal_prop %>%
  ggplot(aes(time, proportion, color = dispersal_mode, fill = dispersal_mode)) +
  geom_point(alpha = 0.1, size = 3) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), alpha = 0.2, linewidth = 2) +
  facet_wrap(~patch_type, scales = "free", labeller = as_labeller(c("Connected" = "(A) Connected", "Rectangular" = "(B) Rectangular", "Winged" = "(C) Winged"))) +
  theme_minimal(base_size = 18) +
  theme(axis.text = element_text(size = 16)) +
  theme(panel.border = element_rect(colour = "darkgrey", fill=NA, linewidth=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "darkgrey", linewidth = 0.5),
        strip.text.x = element_text(hjust = -0.05)) +
  theme(legend.position = "bottom") +
  ylim(0.03, 0.59) +
  xlab("Years since site creation") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#E1BE6A", "#8c510a", "#40B0A6"), name = "Dispersal mode", labels = c("Animal", "Gravity", "Wind")) +
  scale_color_manual(values = c("#E1BE6A", "#8c510a", "#40B0A6"), name = "Dispersal mode", labels = c("Animal", "Gravity", "Wind"))
figureS4

# exporting
# pdf(file = file.path("plots", "figureS4.pdf"), width = 10, height = 5)
# figureS4
# dev.off()



