# function to get zones from study data

recall_binding <- function(study_data, recall_data){
  complete_zone_recall <-  left_join(recall_data, study_data,
                                     by = c("subject_id", "block", "item_num")) %>%
    select(subject_id, block, item_num, event, trial_num, euc_distance)
  return(complete_zone_recall)
}

# Functions to compute overall lag-CRPs

## Function to calculate transitions between each item for overall lag-CRP
calculate_transitions <- function(study_data, recall_data){
  
  # define block_length and all subject_ids
  complete_recall_trans <- data.frame()
  # max_block <- max(study_attn_data$block_recode)
  # block_length <- study_attn_data %>%
  #   group_by(block) %>%
  #   summarise(block_length = max(study_attn_data$trial_num)) %>%
  #   ungroup() # correct variable names
  usable_subjects <- unique(study_data$subject_id)
  
  # loop subjects  
  for (this_subject in usable_subjects){
    
    # Instead of looping from 1 to max_block, get the blocks available for this subject
    max_subject_block <- unique(recall_data$block_recode[recall_data$subject_id == this_subject])
    
    # loop block number
    for (this_block_num in max_subject_block){
      subject_block <- recall_data[
        recall_data$subject_id == this_subject &
          recall_data$block_recode == this_block_num, 
      ]
      # Skip this iteration if the block is empty
      if (nrow(subject_block) == 0) next
      
      # check if the block is "H" or "P"
      if (!is.na(unique(subject_block$block)[1]) && unique(subject_block$block)[1] == "H"){
        # if it's "H", the block length is 18
        block_length <- 18
      } else if(!is.na(unique(subject_block$block)[1]) && unique(subject_block$block)[1] == "P"){
        # if it's "P", the block length is 13
        block_length <- 13
      } 
      
      # get all possible items for that block (non-cumulative)
      all_items <- seq(1, block_length, 1)
      
      # filter the current block for the current subject
      this_block <- recall_data %>%
        filter(subject_id == this_subject & block_recode == this_block_num) %>%
        mutate(actual_lags = lead(as.numeric(item_num)) - item_num) %>%
        mutate(possible_lags =  rep(NA, block_length)) %>%
        
        # calculate distance between recalled items
        mutate(abs_lag_distance = abs(lead(as.numeric(euc_distance)) - euc_distance))
      
      # iterate through each item in the block
      for (this_item in 1:block_length){
        
        # ensure both current and proceeding item are not NAs
        continue =  !is.na(this_block$item_num[this_item]) & !is.na(this_block$item_num[this_item+1])
        
        if (continue){# if Item_num is not NA
          
          curr_item_num <- this_block$item_num[this_item] # get the curr item's item number
          all_items <- all_items[which(all_items != curr_item_num)] # remove this item from the list of all items
          all_pos_trans <- all_items - curr_item_num # get all possible transitions, i.e. difference between curr item num and all the item numbers still in the list
          
          if (length(all_pos_trans) != 0) { # if there are possible transitions list them within the dataframe cell
            this_block$possible_lags[this_item] <- list(all_pos_trans)
          } 
          
          else {
            
            this_block$possible_lags[this_item] <- list(0) # if there are no possible transitions
          } 
        } else {
          this_block$possible_lags[this_item] <- list(0) # if there are no possible transitions
        }
        
      }
      
      complete_recall_trans <- bind_rows(complete_recall_trans, this_block)
      
    }
    
  }  
  return(complete_recall_trans)
}


# loops over subjects and blocks - so CRP is calculated per block per subject
CRP_calculation <- function(complete_recall_trans){
  complete_CRP <- data.frame() # dataframe to store calculated CRPs 
  # dataframe that will get used to calculate CRP at various possible lags
  
  # block_length <- max(study_attn_data$trial_num)
  # max_block <- max(complete_recall_trans$block_recode)  
  usable_subjects <- unique(complete_recall_trans$subject_id)
  
  # loop subjects  
  for (this_subject in usable_subjects){
    
    # Instead of looping from 1 to max_block, get the blocks available for this subject
    max_subject_block <- unique(complete_recall_trans$block_recode[complete_recall_trans$subject_id == this_subject])
    
    # loop block number
    for (this_block_num in max_subject_block){
      # create subject_block for block_length
      subject_block <- complete_recall_trans[
        complete_recall_trans$subject_id == this_subject &
          complete_recall_trans$block_recode == this_block_num, 
      ]
      # check if the block is "H" or "P"
      if (unique(subject_block$block)[1] == "H"){
        # if it's "H", the block length is 18
        block_length <- 18
      } else if (unique(subject_block$block)[1] == "P") {
        # if it's "P", the block length is 13
        block_length <- 13
      }
      
      # create lag_dataframe based on different block_lengths
      lag_dataframe <- data.frame(lags =c(seq(-(block_length - 1), -1, 1), seq(1,(block_length - 1),1))) %>%
        mutate(lags = as.character(lags)) # generaing all possible lags for each item
      
      # get the recall data for each subject's each block
      this_block <- complete_recall_trans %>% 
        filter(subject_id == this_subject & block_recode == this_block_num)
      
      if (dim(this_block)[1] == 0){
        next
      }
      
      # unlisting all the possible lags for each item, and getting counts for each possible lag
      possible_lag_counts <- as.vector(unlist(this_block$possible_lags)) 
      possible_lag_counts <- as.data.frame(table(possible_lag_counts)) %>%
        dplyr::rename(lags = possible_lag_counts, possible_counts = Freq)%>%
        filter(lags != 0) %>%
        mutate(lags = as.character(lags)) 
      possible_lag_counts <- full_join(lag_dataframe, possible_lag_counts, by = "lags") %>%
        mutate(subject_id = this_subject, block = this_block_num) %>%
        select(subject_id, block, lags, possible_counts)
      
      # counting number of actual lags for a particular value
      actual_lag_counts <- this_block %>%
        dplyr::count(actual_lags) %>%
        filter(!is.na(actual_lags)) %>% # removing actual lag of NA
        dplyr::rename(lags = actual_lags, actual_counts = n) %>%
        as.data.frame()%>%
        mutate(lags = as.character(lags))
      actual_lag_counts <- full_join(lag_dataframe, actual_lag_counts, by = "lags") %>%
        mutate(subject_id = this_subject, block = this_block_num) %>%
        select(subject_id, block, lags, actual_counts)  
      
      # calculating CRP as actual/possible transitions. Any NA values are replaced with 0 (either actual or both NA)
      this_block_CRP <- possible_lag_counts %>%
        full_join(actual_lag_counts, by = c("subject_id", "block", "lags")) %>%
        mutate(CRP = actual_counts / possible_counts) %>%
        mutate(CRP = replace_na(CRP, 0)) %>%
        mutate(lags = as.integer(lags))
      
      
      complete_CRP <- bind_rows(complete_CRP, this_block_CRP)
    }
  }
  return(complete_CRP)
}


## Function to calculate percentile rank

temp_percentile_rank <- function(actual, possible) {
  # If fewer than 2 possible transitions, return NA
  if (length(possible) < 2) {
    return(NA)
  }
  
  # Sort possible transitions from largest to smallest
  possible <- sort(possible, decreasing = TRUE)
  
  # Get indices of transitions that match the actual transition
  matches <- which(possible == actual)
  
  if (length(matches) > 0) {
    # Convert R's 1-based indices to zero-based (like Python)
    rank <- mean(matches - 1)  
    # Convert rank to the proportion of possible transitions that were more distant
    ptile_rank <- rank / (length(possible) - 1)
  } else {
    ptile_rank <- NA
  }
  
  return(ptile_rank)
}


# Functions to compute event type-based lag-CRPs using both "common denominator" approach

# Creating event transition column
zone_based_category <- function(complete_recall_trans){
  
  # define an empty data frame
  zone_based_data <- data.frame()
  
  # sort all unique subjects
  usable_subjects <- unique(complete_recall_trans$subject_id)
  
  # Loop through each subject, each block
  for (this_subject in usable_subjects){
    
    # Instead of looping from 1 to max_block, get the blocks available for this subject
    max_subject_block <- unique(complete_recall_trans$block_recode[complete_recall_trans$subject_id == this_subject])
    
    for (this_block_num in max_subject_block){
      
      # create subject_block for block_length
      this_subject_block <- complete_recall_trans[
        complete_recall_trans$subject_id == this_subject &
          complete_recall_trans$block_recode == this_block_num, 
      ]
      
      # filter the current block for the current subject
      this_block <- this_subject_block %>%
        mutate(event_transition = 
                 case_when(
                   # If event is NA or the lead is NA, we make this NA (or "unclassified", etc.)
                   is.na(event) | is.na(lead(event)) ~ NA_character_,
                   # if event matches the lead
                   lead(event) == event ~ "within",
                   # otherwise
                   TRUE ~ "between"
                 )
        )
      # Append these rows to zone_based_data
      zone_based_data <- bind_rows(zone_based_data, this_block)
    }  
  }  
  return(zone_based_data)
}

# Creating event transition column
zone_based_encoding <- function(encoding_data){
  
  # define an empty data frame
  zone_based_encoding <- data.frame()
  
  # sort all unique subjects
  usable_subjects <- unique(encoding_data$subject_id)
  
  # Loop through each subject, each block
  for (this_subject in usable_subjects){
    
    # Instead of looping from 1 to max_block, get the blocks available for this subject
    max_subject_block <- unique(encoding_data$block[encoding_data$subject_id == this_subject])
    
    for (this_block_letter in max_subject_block){
      
      # create subject_block for block_length
      this_subject_block <- encoding_data[
        encoding_data$subject_id == this_subject &
          encoding_data$block == this_block_letter, 
      ]
      
      # filter the current block for the current subject
      this_block <- this_subject_block %>%
        mutate(event_transition = 
                 case_when(
                   # If event is NA or the lead is NA, we make this NA (or "unclassified", etc.)
                   is.na(event) | is.na(lead(event)) ~ NA_character_,
                   # if event matches the lead
                   lead(event) == event ~ "within",
                   # otherwise
                   TRUE ~ "between"
                 )
        ) %>%
        
        # calculates the distance for each lag
        mutate(abs_lag_distance = abs(lead(as.numeric(euc_distance)) - euc_distance))
      
      # Append these rows to zone_based_encoding data frame
      zone_based_encoding <- bind_rows(zone_based_encoding, this_block)
    }  
  }  
  return(zone_based_encoding)
}


# Function for computing actual transitions 
actual_lag_finder <- function(this_block, lag_dataframe) {
  actual_lag_counts <- this_block %>%
    count(actual_lags) %>%
    filter(!is.na(actual_lags)) %>%
    dplyr::rename(lags = actual_lags, actual_counts = n) %>%
    as.data.frame()%>%
    mutate(lags = as.character(lags))
  
  actual_lag_counts <- full_join(lag_dataframe, actual_lag_counts, by = "lags")
  
  
  return(actual_lag_counts)
}


# Function for computing possible transitions 
possible_lag_finder <- function (possible_lag_column, lag_dataframe){
  possible_lag_counts <- as.vector(unlist(possible_lag_column)) 
  possible_lag_counts <- as.data.frame(table(possible_lag_counts)) %>%
    dplyr::rename(lags = possible_lag_counts, possible_counts = Freq)%>%
    filter(lags !=0) %>%
    mutate(lags = as.character(lags)) 
  possible_lag_counts <- full_join(lag_dataframe, possible_lag_counts, by = "lags") 
  
  return(possible_lag_counts)
  
}

# Function for computing CRP
CRP_calculator <- function(possible_lag_counts, actual_lag_counts){
  
  this_block_CRP <- possible_lag_counts %>%
    full_join(actual_lag_counts, by = c("Subject_ID", "BlockNumber", "lags", "zone")) %>%
    mutate(CRP = actual_counts / possible_counts) %>%
    mutate(CRP = replace_na(CRP, 0)) %>%
    mutate(lags = as.integer(lags))
  
  return(this_block_CRP)
  
}


# Function for common-denominator approach - CRP
common_deno_CRP <- function(zone_based_data){
  
  complete_CRP_within <- data.frame()
  complete_CRP_between <- data.frame()
  
  usable_subjects <- unique(zone_based_data$subject_id)
  
  # loop subjects  
  for (this_subject in usable_subjects){
    
    # Instead of looping from 1 to max_block, get the blocks available for this subject
    max_subject_block <- unique(zone_based_data$block_recode[zone_based_data$subject_id == this_subject])
    
    for (this_block_num in max_subject_block){
      
      # create subject_block for block_length
      this_subject_block <- zone_based_data[
        zone_based_data$subject_id == this_subject &
          zone_based_data$block_recode == this_block_num, 
      ]
      # check if the block is "H" or "P"
      if (unique(this_subject_block$block)[1] == "H"){
        # if it's "H", the block length is 18
        block_length <- 18
      } else if (unique(this_subject_block$block)[1] == "P") {
        # if it's "P", the block length is 13
        block_length <- 13
      }
      
      # create lag_dataframe based on different block_lengths
      lag_dataframe <- data.frame(lags =c(seq(-(block_length - 1), -1, 1), seq(1,(block_length - 1),1))) %>%
        mutate(lags = as.character(lags)) # generaing all possible lags for each item
      
      # get the recall data for each subject's each block
      this_block <- zone_based_data %>% 
        filter(subject_id == this_subject & block_recode == this_block_num)
      
      if (dim(this_block)[1] == 0){
        next
      }
      
      # unlisting all the possible lags for each item, and getting counts for each possible lag
      possible_lag_counts <- as.vector(unlist(this_block$possible_lags)) 
      possible_lag_counts <- as.data.frame(table(possible_lag_counts)) %>%
        dplyr::rename(lags = possible_lag_counts, possible_counts = Freq)%>%
        filter(lags != 0) %>%
        mutate(lags = as.character(lags)) 
      possible_lag_counts <- full_join(lag_dataframe, possible_lag_counts, by = "lags") %>%
        mutate(subject_id = this_subject, block = this_block_num) %>%
        select(subject_id, block, lags, possible_counts)
      
      # Compute actual lag counts and CRP separately for "within" and "between"
      for (transition_type in c("within", "between")) {
        this_subblock <- this_block %>% filter(event_transition == transition_type)
        
        if (nrow(this_subblock) == 0) next
        
        actual_counts <- this_subblock %>%
          count(actual_lags) %>%
          rename(lags = actual_lags, actual_counts = n) %>%
          mutate(lags = as.character(lags)) %>%
          full_join(lag_dataframe, by = "lags") %>%
          mutate(subject_id = this_subject, block = this_block_num) %>%
          select(subject_id, block, lags, actual_counts)
        
        this_CRP <- possible_lag_counts %>%
          full_join(actual_counts, by = c("subject_id", "block", "lags")) %>%
          mutate(CRP = actual_counts / possible_counts) %>%
          mutate(CRP = replace_na(CRP, 0), lags = as.integer(lags), event_transition = transition_type)
        
        if (transition_type == "within") {
          complete_CRP_within <- bind_rows(complete_CRP_within, this_CRP)
        } else if (transition_type == "between") {
          complete_CRP_between <- bind_rows(complete_CRP_between, this_CRP)
        }
      }
    }
  }
  
  return(list(within = complete_CRP_within, between = complete_CRP_between))
  
}

run_paired_ttest <- function(subdata) {
  
  if (length(unique(subdata$event_transition)) < 2) {
    return(NULL)
  }
  
  wide <- subdata %>% 
    select(subject_id, lags, mean_CRP, event_transition) %>% 
    pivot_wider(
      names_from = event_transition,
      values_from = mean_CRP
    ) %>% 
    filter(!is.na(within) & !is.na(between))
  
  model <- t.test(wide$within, wide$between, paired = TRUE)
  
  # Group means and SDs
  mean_within <- mean(wide$within, na.rm = TRUE)
  sd_within <- sd(wide$within, na.rm = TRUE)
  mean_between <- mean(wide$between, na.rm = TRUE)
  sd_between <- sd(wide$between, na.rm = TRUE)
  
  # Cohen'd 
  cohen_value <- cohens_d(wide$within, wide$between, paired = TRUE)
  
  # CI
  ci_lower <- model$conf.int[1]
  ci_upper <- model$conf.int[2]
  
  # summarise
  data.frame(
    lag = unique(subdata$lags),
    cond1 = "within",
    cond2 = "between",
    mean_within = mean_within,
    sd_within = sd_within, 
    mean_between = mean_between,
    sd_between = sd_between,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    df = unname(model$parameter),
    t = unname(model$statistic),
    p_raw = model$p.value,
    cohens_d = cohen_value$Cohens_d[1])
}


# Function for common-denominator approach - temporal cluster
common_temp_percentile_rank <- function(actual_within, actual_between, possible) {
  # If fewer than 2 possible transitions, return NA
  if (length(possible) < 2) {
    return(NA)
  }
  
  # Sort possible transitions from largest to smallest
  possible <- sort(possible, decreasing = TRUE)
  
  # Get indices of transitions that match the actual transition
  matches_within <- which(possible == actual_within)
  matches_between <- which(possible == actual_between)
  
  if (length(matches_within) > 0) {
    # Convert R's 1-based indices to zero-based (like Python)
    rank_within <- mean(matches_within - 1)  
    # Convert rank to the proportion of possible transitions that were more distant
    ptile_rank_within <- rank_within / (length(possible) - 1)
  } else {
    ptile_rank_within <- NA
  }
  
  if (length(matches_between) > 0) {
    # Convert R's 1-based indices to zero-based (like Python)
    rank_between <- mean(matches_between - 1)  
    # Convert rank to the proportion of possible transitions that were more distant
    ptile_rank_between <- rank_between / (length(possible) - 1)
  } else {
    ptile_rank_between <- NA
  }
  
  return(list(
    ptile_rank_within = ptile_rank_within,
    ptile_rank_between = ptile_rank_between
  ))
}

