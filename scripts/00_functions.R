########
## SCRIPT NAME: 00_functions.R
## AUTHOR: Katherine Hulting
## PURPOSE: Script containing custom functions used throughout analysis.
## PRODUCTS: NA
#########

# loading libraries
librarian::shelf(tidyverse, vegan) # Install missing packages and load needed libraries

### compute jaccard dissimilarity between patch pairs for each year
compute_convergence_jaccard <- function(df) {
  # convert data to correct format
  df_wide <- df %>% 
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    arrange(unique_id, time) %>%  # Ensure years are sorted properly
    mutate(unique_id = paste(unique_id, time, sep = "-")) %>% # create unique ID
    dplyr::select(-block, -patch_type, -patch, -year, -time) %>% # drop unnecessary columns
    column_to_rownames("unique_id") # make unique ID the rownames
  
  jac <- vegdist(df_wide, method = "jaccard") # jaccard dissimilarity
  jac_df <- as.data.frame(as.table(as.matrix(jac))) # make dataframe
  names(jac_df) <- c("unique_id1", "unique_id2", "jaccard") # set names
  
  jac_df <- jac_df %>%
    separate(unique_id1, into = c("block1", "patch_rep1", "patch_type1", "time1"), remove = T) %>% # seperate unique ID of patch 1 into its parts
    separate(unique_id2, into = c("block2", "patch_rep2", "patch_type2", "time2"), remove = T) %>% # seperate unique ID of patch 2 into its parts
    mutate(unique_id1 = paste(block1, patch_rep1, patch_type1, sep = "-")) %>% # recreate unique ID without time
    mutate(unique_id2 = paste(block2, patch_rep2, patch_type2, sep = "-")) %>% # recreate unique ID without time
    filter(!(unique_id1 == unique_id2)) %>% # remove observations comparing patch to itself
    filter(time1 == time2) %>% # only keeping observations that are comparing patches in the same year
    mutate(patch_pair = paste(patch_type1, patch_type2, sep = "-")) %>% # creating variable definining which patch types are being compared
    dplyr::select(-block2, -time2, -patch_rep1, -patch_rep2, -patch_type1, -patch_type2) %>% # dropping irrelevant variables
    rename(block = block1, time = time1) %>% # renaming 
    mutate(time = as.numeric(time)) %>% # making numeric
    mutate(patch_pair = dplyr::case_when( # making sure patch type comparison character strings are consistent
      patch_pair %in% c("Connected-Rectangular", "Rectangular-Connected") ~ "Connected-Rectangular",
      patch_pair %in% c("Connected-Winged", "Winged-Connected") ~ "Connected-Winged",
      patch_pair %in% c("Rectangular-Winged", "Winged-Rectangular") ~ "Rectangular-Winged",
      patch_pair %in% c("Center-Connected", "Connected-Center") ~ "Center-Connected",
      patch_pair %in% c("Center-Rectangular", "Rectangular-Center") ~ "Center-Rectangular",
      patch_pair %in% c("Center-Winged", "Winged-Center") ~ "Center-Winged",
      .default = patch_pair
    ))
  
  return(jac_df) # returning dataframe
}



#### calculate raw composition change between consecutive years ####
compute_composition_change <- function(df) {
  
  df_wide <- df %>% # convert data to correct format
    pivot_wider(names_from = sppcode, values_from = n, values_fill = 0) %>% # wide format
    arrange(unique_id, time) %>%  # ensure years are sorted properly
    mutate(unique_id = paste(unique_id, time, sep = "-")) %>%
    dplyr::select(-block, -patch_type, -patch, -year) %>% # remove unneeded columns
    column_to_rownames("unique_id")
  
  # calculate losses and gains between consecutive years
  turnover <- df_wide %>%
    arrange(time) %>% # ensure years are sorted properly
    mutate(
      gains = rowSums(across(-time, ~ (.x == 1 & lag(.x) == 0)), na.rm = TRUE), # calculate sp gains, where 0 -> 1
      losses = rowSums(across(-time, ~ (.x == 0 & lag(.x) == 1)), na.rm = TRUE), # calculate sp losses, 1 -> 0
      stayed_present = rowSums(across(-time, ~ (.x == 1 & lag(.x) == 1)), na.rm = TRUE), # calculate sp staying consistent, 1 -> 1
      changed_total = gains + losses # total turnover
    ) %>%
    dplyr::select(time, gains, losses, stayed_present, changed_total)
  
  
  return(turnover)
}
