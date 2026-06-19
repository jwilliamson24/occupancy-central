## =================================================
##
## Title: multinom-data-manipulation
## Author: Jasmine Williamson
## Date Created: 6/16/2026
##
## Description: Looked at raw counts for age classes per plot, and looks like there arent
## enough to justify doing a separate multinom for each spp/age class (30 indivs).
## Also created count dataframes for multinomial model, and fixed missing obs values from 2023.
##
## =================================================


## settings -----------------------------------------------------------------------------------------------

  # rm(list=ls())
  setwd("/Users/jasminewilliamson/Library/CloudStorage/OneDrive-Personal/Documents/Academic/OSU/Git/occupancy-central/data")
  
  # library(unmarked)
  # library(ggplot2)
  # library(stats)
  # library(MASS)
  library(tidyverse)
  library(dplyr)


## load data ----------------------------------------------------------------------------------------

  subplot <- read.csv("data/covariate matrices/subplot.complete.csv")
  sals <- read.csv("data/occupancy/sals.complete.csv", 
                   colClasses = c(landowner="factor", stand="character", trt="factor",
                                  obs="factor", subplot="factor", recap="factor",
                                  pass="factor", spp="factor", cover_obj="factor", 
                                  substrate="factor", age_class="factor"))


## df with all plots/occasions ----------------------------------------------------------------------

  #creating occupancy df with non-detections
  df.new <- subplot[,c(1,3,5,6,7,9)] # siteID, stand, subplots, year, date, treatment
  df.new$subplot <- as.factor(df.new$subplot)
  df.new <- df.new[order(df.new$site_id, df.new$subplot),] # reorder
  
  # add 3 passes for each subplot
  df.new <- df.new %>%
    slice(rep(1:n(), each = 3)) %>%  # duplicate each row 3 times
    group_by(site_id, subplot) %>%
    mutate(pass = 1:3) %>%
    ungroup()
  df.new$pass <- as.factor(df.new$pass)
  df.new <- as.data.frame(df.new)


# merge sals with full subplot/pass df -------------------------------------------------------------
  
  sals <- sals %>%
    dplyr::filter(recap != 1)
  sals.new <- sals[,c(1,11:18)] 
  sals.new <- filter(sals.new, spp %in% c("OSS", "ENES"))
  sals.new$detect = 1  

  
  # merge
  
  #add in sites with no detections by merging with df that has all site/subplot combos listed
  df.merge <- full_join(df.new,sals.new,by=c("site_id","subplot","pass"))
  df.merge$detect <- ifelse(is.na(df.merge$detect), 0, df.merge$detect) #make NA's = 0

  
#### raw counts of sals ---------------------------------------------------------------------------
  
  sal_counts <- df.merge %>%
    filter(detect == 1) %>%
    group_by(spp, age_class) %>%
    summarise(count = n(), .groups = "drop")

  sal_counts_wide <- sal_counts %>%
    unite("spp_age", spp, age_class, sep = "_") %>%
    pivot_wider(names_from = spp_age, values_from = count, values_fill = 0)
  sal_counts_wide <- sal_counts_wide[-6] #removing enes unknown age class
  sal_counts_wide <- sal_counts_wide[-7] #removing oss unknown age class
  
  sum(df.merge$detect[df.merge$spp == "OSS"], na.rm = TRUE)
  sum(df.merge$detect[df.merge$spp == "ENES"], na.rm = TRUE)
  
  
#### fix observer data ------------------------------------------------------------------------

## observer data for sites without sals was not inputted into the master data sheet
## it was only on the "salamanders" tab for each sally found
## below i organized the obs data that we do have, downloaded that csv,
## then i went back to paper datasheets and inputted remaining obs info into that csv
## the person who filled out datasheet was V1, then other two were randomly assigned
## i only did this for part of 2023, and then decided it actually didnt matter, and 
## wrote code below to fill in randomly


# get usable df
    
    # subset
    df.merge.new <- as.data.frame(df.merge[, c("site_id", "date", "stand", "trt", "year", "subplot", "pass", "obs")])
    
    # condense indiv sal replicates per pass, and widen pass reps
    df.new.wide <- df.merge.new %>% 
      group_by(site_id, stand, year, subplot, pass, date, trt) %>%
      summarise(obs = first(obs), .groups = "drop") %>%
      pivot_wider(names_from = pass, values_from = obs, names_prefix = "V") 
    
    # convert columns to character =, then change NA to blank cell for csv readability
    df.new.wide <- df.new.wide %>%
      mutate(across(starts_with("V"), as.character))
    
    df.new.wide[is.na(df.new.wide)] <- ""    
    
    write.csv(df.new.wide, "data/abundance/df-merge-for-obs.csv", row.names = FALSE)
    
# here is where i inputted some 2023 data in the above csv, gave up, and uploaded partially fixed df below
    
    part.fixed <- read.csv("data/abundance/df-merge-for-obs.csv")
    
# code from claude for inputting random obs assignments for rest of df
    
    # create function for filling techs
    fill_techs <- function(row, techs) {
      present <- row[!is.na(row) & row != ""]      # which initials are already filled in this row?
      missing <- setdiff(techs, present)           # which of the 3 techs are NOT yet in this row?
      row[is.na(row) | row == ""] <- missing[seq_len(sum(is.na(row) | row == ""))]  # fill blanks with missing names
      row                                          # return the completed row
    }
    
    techs_2023 <- c("JW", "JB", "BZ")
    techs_2024 <- c("JW", "RMM", "SLG")
    
    # make sure V1, V2, V3 are characters, not factors
    part.fixed <- part.fixed %>% 
      rowwise() %>%
      mutate(
        across(c(V1, V2, V3), as.character)
      ) %>%
      ungroup()
    
    # loop through df
    for (i in seq_len(nrow(part.fixed))) {
      yr <- part.fixed$year[i]                              # which year 
      techs <- if (yr == 2023) techs_2023 else techs_2024   # select correct techs list
      row_vals <- c(part.fixed$V1[i], part.fixed$V2[i], part.fixed$V3[i]) # view current row values
      filled <- fill_techs(row_vals, techs)                 # run fill techs function
      part.fixed$V1[i] <- filled[1]   # write new rows back into df
      part.fixed$V2[i] <- filled[2]
      part.fixed$V3[i] <- filled[3]
    }
    
    
# take back to long format
    
    obs.fixed.long <- part.fixed %>%
      pivot_longer(
        cols = c(V1, V2, V3),
        names_to = "pass",
        values_to = "obs"
      ) %>%
      mutate(pass = as.factor(readr::parse_number(pass)))

# merge with df.merge from way above to have fixed obs values
    
    df.merge.fixed <- df.merge %>%
      left_join(
        obs.fixed.long %>% 
          select(site_id, subplot, year, pass, obs_filled = obs) %>%
          mutate(subplot = as.factor(subplot)),
        by = c("site_id", "subplot", "year", "pass")
      ) %>%
      mutate(obs = coalesce(obs, obs_filled)) %>%
      select(-obs_filled)
    
    
#### make counts dataframe with updated df ------------------------------------------------------------------
    
    # oss
    
    # subsest oss and merge with df pass-level df to fill in long format
    sals.oss <- subset(sals.new, spp=="OSS")
    df.merge.o <- full_join(df.new,sals.oss,by=c("site_id","subplot","pass"))
    df.merge.o$detect <- ifelse(is.na(df.merge.o$detect), 0, df.merge.o$detect) #make NA's = 0
    
    # merge with fixed df to fill in missing obs
    df.merge.o <- df.merge.o %>%
      left_join(
        obs.fixed.long %>% 
          select(site_id, subplot, year, pass, obs_filled = obs) %>%
          mutate(subplot = as.factor(subplot)),
        by = c("site_id", "subplot", "year", "pass")
      ) %>%
      mutate(obs = coalesce(obs, obs_filled)) %>%
      select(-obs_filled)
    
    # summarize count and detections per site/subplot/pass so the rows are unique (no repliate rows)
    df.counts.o <- df.merge.o %>%
      group_by(site_id, stand, trt, year, subplot, pass) %>%
      summarise(
        count = sum(detect, na.rm = TRUE),
        detect = as.integer(sum(detect, na.rm = TRUE) > 0),
        obs = first(obs),
        .groups = "drop"
      ) 
    
    
    write.csv(df.counts.o, "data/abundance/pass-level-counts-o.csv", row.names = FALSE)
    
    
    
    # enes
    
    # subsest enes and merge with df pass-level df to fill in long format
    sals.enes <- subset(sals.new, spp=="ENES")
    df.merge.e <- full_join(df.new,sals.enes,by=c("site_id","subplot","pass"))
    df.merge.e$detect <- ifelse(is.na(df.merge.e$detect), 0, df.merge.e$detect) #make NA's = 0
    
    # merge with fixed df to fill in missing obs
    df.merge.e <- df.merge.e %>%
      left_join(
        obs.fixed.long %>% 
          select(site_id, subplot, year, pass, obs_filled = obs) %>%
          mutate(subplot = as.factor(subplot)),
        by = c("site_id", "subplot", "year", "pass")
      ) %>%
      mutate(obs = coalesce(obs, obs_filled)) %>%
      select(-obs_filled)
    
    # summarize count and detections per site/subplot/pass so the rows are unique (no repliate rows)
    df.counts.e <- df.merge.e %>%
      group_by(site_id, stand, trt, year, subplot, pass) %>%
      summarise(
        count = sum(detect, na.rm = TRUE),
        detect = as.integer(sum(detect, na.rm = TRUE) > 0),
        obs = first(obs),
        .groups = "drop"
      ) 
    
    
    write.csv(df.counts.e, "data/abundance/pass-level-counts-e.csv", row.names = FALSE)
    
    
    
    
    
    
    ######## these last two dfs are the ones that load into the multinomial model script ########
    
    
    
    
    
    
    
    
    
    
    
    