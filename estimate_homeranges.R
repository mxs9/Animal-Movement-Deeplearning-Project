#You need to install of these libraries in order to run the functions
library(tidyverse); 
library(lubridate); 
library(ctmm); 
library(parallel); 
library(foreach); 
library(gridExtra)
library(dplyr)
library(stringr)
doParallel::registerDoParallel(4)

#Examples of how to run all the functions are at the bottom of the script
compare_overlap <- function(before_model, after_model){
  #Takes in two lists. 
  #Before_model takes a list comtaining two fitted models of the individuals before encounter.
  #After_model is the list of fitted models after encounter
  #Returns a data.frame that lets you compare values easier
  before_overlap <- overlap(before_model)
  after_overlap <- overlap(after_model)
  low_comp <- c(before_overlap[2], after_overlap[2])
  est_comp <- c(before_overlap[6], after_overlap[6])
  high_comp <- c(before_overlap[10], after_overlap[10])
  return(data.frame( rows = c("before","after"), low = low_comp, est = est_comp, high = high_comp, 
                     row.names = "rows"))
}
----------------------------------------------
  #The following is a script that runs all the AKDE analysis and saves all the data into a folder it creates
  #This splits the data EXACTLY before and after the encounter time, it does not factor in "transition times"
  #To factor in transitions for on encounter, use run_akde_with_data instead
  #You need to input the following
  #name1: the name of individual one (make sure it's spelled exactly)
  #name2: the name of individual two (make sure it's spelled exaclty)
  #study: the study the pair is from
  #date_encounter: the date of encounter as a lubridate datetime (use ymd_hms("2016-04-23 14:30:00") to create them)
run_akde_analysis <- function(name1, name2, study, date_encounter){
    folder_name <- str_c(name1, "_", name2, "_Homeranges", sep = "")
    #Change the folder_dir string to the folder where you want all your results to be stored
    #Do not add any slashes to the end of the string
    #You ONLY need to change this variable
    folder_dir <- "C:/workhard/Movement Seminar/MovementSeminar/Code_Mingxi"
    folder_path <- str_c(folder_dir, folder_name, sep = "\\")
    dir.create(folder_path)
    
    id1 <- movement_data_id_all4 %>% filter(study.id == study, individual.local.identifier ==name1)
    id2 <- movement_data_id_all4 %>% filter(study.id == study,individual.local.identifier == name2)
    
    before <- list(filter(id1, timestamp < date_encounter),
                   filter(id2, timestamp < date_encounter)) %>% map(., as.telemetry)
    
    after <- list(filter(id1, timestamp > date_encounter),
                  filter(id2, timestamp > date_encounter)) %>% map(., as.telemetry)
    
    save(before, file = str_c(folder_path, "\\", folder_name, "_before_data.rda"))
    save(after, file = str_c(folder_path, "\\", folder_name, "_after_data.rda"))
    
    #Sets the projections to be the same as before. Change both to after is you want them standardized to after
    projection(before) <- median(before)
    projection(after) <- median(before)
    
    #Print the number of entries in each split dataset, gives you a sense of scale
    print(map_dbl(before, nrow))
    print(map_dbl(after, nrow))
    
    GUESS <- lapply(Hedley_Diabla_before[1:2], function(b) ctmm.guess(b,interactive=FALSE) )
    FITS_Hedley_Diabla_before <- lapply(1:2, function(i) ctmm.select(Hedley_Diabla_before[[i]],GUESS[[i]]) )
    
    akde_all <- function(data){
      guess <- map(data, ctmm.guess, interactive = FALSE)
      fit <- map2(data, guess, ctmm.select)
      akde_final <- map2(data, fit, akde)
      return(list(fit = fit, akde = akde_final))
    }
    
    seq_tele <- list(before = before, after = after)
    model_list <- foreach ( i = seq_along(seq_tele), .packages = c("ctmm","lubridate","tidyverse")) %dopar% 
      {akde_all(seq_tele[[i]])}
    
    before_model <- model_list[[1]]
    after_model <- model_list[[2]]
    
    save(before_model, file = str_c(folder_path, "\\", folder_name, "_before_model.rda", sep = ""))
    save(after_model, file = str_c(folder_path, "\\", folder_name, "_after_model.rda", sep = ""))
    
    #Overlap analysis table
    overlap_analysis <- compare_overlap(before_model$fit, after_model$fit)
    png(filename = str_c(folder_path, "\\", folder_name,"_overlap_table.png", sep = ""), 
        height = 50*nrow(overlap_analysis), width = 200*ncol(overlap_analysis))
    grid.table(overlap_analysis)
    dev.off()
    #Save Overlap Table as dataframe
    save(overlap_analysis, file = str_c(folder_path, "\\", folder_name,"_overlap_table.rda", sep = ""))
    
    pair_names <- names(before_model$akde)
    #Make the effective home range table for reporting
    eff_home_range <- data.frame(x = c(before_model$akde[[1]]$DOF.area[1], 
                                       after_model$akde[[1]]$DOF.area[1]),
                                 y = c(before_model$akde[[2]]$DOF.area[1], 
                                       after_model$akde[[2]]$DOF.area[1]))
    colnames(eff_home_range) <- pair_names
    rownames(eff_home_range) <- c("Before","After")
    
    png(filename = str_c(folder_path, "\\", folder_name,"_effective_homerange_table.png", sep = ""), 
        height = 50*nrow(eff_home_range), width = 200*ncol(eff_home_range))
    grid.table(eff_home_range)
    dev.off()
    
    save(eff_home_range, file = str_c(folder_path, "\\", folder_name,"_effective_homerange_table.rda", sep = ""))
    
    extents <- rbind(ctmm::extent(before_model$fit[[1]]), ctmm::extent(before_model$fit[[2]]),
                     ctmm::extent(before_model$akde[[1]]), ctmm::extent(before_model$akde[[2]]),
                     ctmm::extent(after_model$fit[[1]]), ctmm::extent(after_model$fit[[2]]),
                     ctmm::extent(after_model$akde[[1]]), ctmm::extent(after_model$akde[[2]]))
    
    
    x_range <- c((min(extents$x) - 500 %#% "meters"),(max(extents$x) + 500 %#% "meters"))
    y_range <- c((min(extents$y) - 500 %#% "meters"),(max(extents$y) + 500 %#% "meters"))
    
    sorted_names <- sort(c(name1, name2), decreasing = F)
    #Before encounter home range estimates
    png(filename = str_c(folder_path, "\\", folder_name,"_before.png", sep = ""),
        width = 13, height = 6, units = "in", res = 75)
    plot(before, col = color(before, by = "individual"), UD = before_model$akde,
         ylim = y_range, xlim = x_range, level.UD = 0.95,
         main = str_c(name1, " and ", name2, " before encounter: ",date_encounter))
    legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex =2.5) 
    
    dev.off()
    #After encounter home range estimates
    png(filename = str_c(folder_path, "\\", folder_name,"_after.png", sep = ""), 
        width = 13, height = 6, units = "in", res = 75)
    plot(after, col = color(after, by = "individual"), UD = after_model$akde, 
         ylim = y_range, xlim = x_range, 
         main = str_c(name1, " and ", name2, " after encounter: ",date_encounter))
    legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex = 2.5) 
    dev.off()
    
  }
---------------------------------------------
  #This runs that same AKDE analysis as in the function above, expect you give it the data list to run the analysis on
  #This allows you to split and discard some data how you like
  #It does not work with pairs with multiple encounter times 
  #You must input
  #before_list: a list containing the data frames of both individual's before data. It must be names corresponding to the names of the pairs
  #after_list: same as before list except it contains the data after the encounter
  #date_encounter: the date of the encounter, formatted as a lubridate datetime (use ymd_hms() function)
  run_akde_with_data <- function(before_list, after_list, date_encounter){
    pair_names <- names(before_list)
    name1 <- pair_names[1]
    name2 <- pair_names[2]
    folder_name <- str_c(name1, "_", name2, "_Homeranges", sep = "")
    #Change the folder_dir string to the folder where you want all your results to be stored
    #Do not add any slashes to the end of the string
    #You ONLY need to change this variable
    folder_dir <- "C:/workhard/Movement Seminar/MovementSeminar/Code_Mingxi"
    folder_path <- str_c(folder_dir, folder_name, sep = "/")
    dir.create(folder_path)
    
    before <- before_list %>% map(., as.telemetry)
    
    after <-  after_list %>% map(., as.telemetry)
    
    #Sets the projections to be the same as before. Change both to after is you want them standardized to after
    projection(before[1]) <- median(before)
    projection(before[2]) <- median(before)
    projection(after[1]) <- median(before)
    projection(after[2]) <- median(before)
    
    #Print the number of entries in each split dataset, gives you a sense of scale
    print(map_dbl(before, nrow))
    print(map_dbl(after, nrow))
    
    
    akde_all <- function(data){
      guess <- map(data, ctmm.guess, interactive = FALSE)
      fit <- map2(data, guess, ctmm.select)
      akde_final <- map2(data, fit, akde)
      return(list(fit = fit, akde = akde_final))
    }
    
    seq_tele <- list(before = before, after = after)
    model_list <- foreach ( i = seq_along(seq_tele), .packages = c("ctmm","lubridate","tidyverse")) %dopar% 
      {akde_all(seq_tele[[i]])}
    
    before_model <- model_list[[1]]
    after_model <- model_list[[2]]
    
    save(before_model, file = str_c(folder_path, "/", folder_name, "_before_model.rda", sep = ""))
    save(after_model, file = str_c(folder_path, "/", folder_name, "_after_model.rda", sep = ""))
    
    #Overlap analysis table
    overlap_analysis <- compare_overlap(before_model$fit, after_model$fit)
    png(filename = str_c(folder_path, "/", folder_name,"_overlap_table.png", sep = ""), 
        height = 50*nrow(overlap_analysis), width = 200*ncol(overlap_analysis))
    grid.table(overlap_analysis)
    dev.off()
    #Save Overlap Table as dataframe
    save(overlap_analysis, file = str_c(folder_path, "/", folder_name,"_overlap_table.rda", sep = ""))
    
    pair_names <- names(before_model$akde)
    #Make the effective home range table for reporting
    eff_home_range <- data.frame(x = c(before_model$akde[[1]]$DOF.area[1], 
                                       after_model$akde[[1]]$DOF.area[1]),
                                 y = c(before_model$akde[[2]]$DOF.area[1], 
                                       after_model$akde[[2]]$DOF.area[1]))
    colnames(eff_home_range) <- pair_names
    rownames(eff_home_range) <- c("Before","After")
    
    png(filename = str_c(folder_path, "/", folder_name,"_effective_homerange_table.png", sep = ""), 
        height = 50*nrow(eff_home_range), width = 200*ncol(eff_home_range))
    grid.table(eff_home_range)
    dev.off()
    
    save(eff_home_range, file = str_c(folder_path, "/", folder_name,"_effective_homerange_table.rda", sep = ""))
    
    extents <- rbind(ctmm::extent(before_model$fit[[1]]), ctmm::extent(before_model$fit[[2]]),
                     ctmm::extent(before_model$akde[[1]]), ctmm::extent(before_model$akde[[2]]),
                     ctmm::extent(after_model$fit[[1]]), ctmm::extent(after_model$fit[[2]]),
                     ctmm::extent(after_model$akde[[1]]), ctmm::extent(after_model$akde[[2]]))
    
    
    x_range <- c((min(extents$x) - 500 %#% "meters"),(max(extents$x) + 500 %#% "meters"))
    y_range <- c((min(extents$y) - 500 %#% "meters"),(max(extents$y) + 500 %#% "meters"))
    
    sorted_names <- sort(pair_names, decreasing = F)
    #Before encounter home range estimates
    png(filename = str_c(folder_path, "/", folder_name,"_before.png", sep = ""),
        width = 13, height = 6, units = "in", res = 75)
    plot(before, col = color(before, by = "individual"), UD = before_model$akde,
         ylim = y_range, xlim = x_range, level.UD = 0.95,
         main = str_c(name1, " and ", name2, " before encounter: ",date_encounter))
    legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex =2.5) 
    
    dev.off()
    #After encounter home range estimates
    png(filename = str_c(folder_path, "/", folder_name,"_after.png", sep = ""), 
        width = 13, height = 6, units = "in", res = 75)
    plot(after, col = color(after, by = "individual"), UD = after_model$akde, 
         ylim = y_range, xlim = x_range, 
         main = str_c(name1, " and ", name2, " after encounter: ",date_encounter))
    legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex = 2.5) 
    dev.off()
    
  }

#Examples-------------------------------

#make sure you have the movement_data_id_all4 dataset loaded in before trying these examples

#run_akde_analysis example

#make datetime to use as date_encounter
#Apethe_Merimela_encounter_date <- ymd_hms("2016-08-12 03:00:00")
#run_akde_analysis("Apethe", "Merimela", "Walton", Apethe_Merimela_encounter_date)

---------------------------------------
  #run_akde_with_data example
  #Apethe_Merimela_encounter_date <- ymd_hms("2012-05-26 03:10:37")
  BORR_ORRR_encounter_date <-ymd_hms("2019/7/11  23:02:06")
#Get all data on both pairs
BORR_data <- movement_data_id_all4 %>% filter(individual.local.identifier == "BORR", study.id == "Bertreaux")
ORRR_data <- movement_data_id_all4 %>% filter(individual.local.identifier == "ORRR", study.id == "Bertreaux")

# #Split before data sets
# Apethe_before <- Apethe_data %>% filter(timestamp < Apethe_Merimela_encounter_date)
# Merimela_before <- Merimela_data %>% filter(timestamp < Apethe_Merimela_encounter_date)
# #Split after data sets
# Apethe_after <- Apethe_data %>% filter(timestamp > Apethe_Merimela_encounter_date)
# Merimela_after <- Merimela_data %>% filter(timestamp > Apethe_Merimela_encounter_date)

#Split before data sets
BORR_before <- BORR_data %>% filter(timestamp < ymd_hms("2019/7/11  23:02:06"))
ORRR_before <- ORRR_data %>% filter(timestamp < ymd_hms("2019/7/11  23:02:06"))
#Split after data sets
BORR_after <- BORR_data %>% filter(timestamp > ymd_hms("2019/7/11  23:02:06"))
ORRR_after <- ORRR_data %>% filter(timestamp > ymd_hms("2019/7/11  23:02:06"))

#Make lists to input into function (MAKE SURE IT"S NAMES)
before_list <- list(BORR = BORR_before, ORRR = ORRR_before)
after_list <- list(BORR = BORR_after, ORRR = ORRR_after)
#Run Function!
run_akde_with_data(before_list, after_list, BORR_ORRR_encounter_date)


