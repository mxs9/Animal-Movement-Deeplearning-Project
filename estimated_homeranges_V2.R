#code for estimating the homeranges
#overlap of UDS
#overlap of model fitting
#making plots of home ranges


#Author
#first version by Maria Nikolaitchik 05/2022
#second version by Qianru Liao --08/2022

#You need to install of these libraries in order to run the functions
library(tidyverse); library(lubridate); library(ctmm); library(parallel); library(foreach); library(gridExtra)
doParallel::registerDoParallel(4)

#Examples of how to run all the functions are at the bottom of the script
compare_overlap <- function(before_model, after_model){
  #Takes in two lists. 
  #Before_model takes a list comtaining two fitted models of the individuals before encounter.
  #After_model is the list of fitted models after encounter
  #Returns a data.frame that lets you compare values easier
  before_overlap <- overlap(before_model)
  after_overlap <- overlap(after_model)
  low_comp <- c(before_overlap$CI[2], after_overlap$CI[2]) #change it to CI --Qianru
  est_comp <- c(before_overlap$CI[6], after_overlap$CI[6])#change it to CI --Qianru
  high_comp <- c(before_overlap$CI[10], after_overlap$CI[10])#change it to CI --Qianru
  return(data.frame( rows = c("before","after"), low = low_comp, est = est_comp, high = high_comp, 
                     row.names = "rows"))
}
#----------------------------------------------
#The following is a script that runs all the AKDE analysis and saves all the data into a folder it creates
#This splits the data EXACTLY before and after the encounter time, it does not factor in "transition times"
#To factor in transitions for on encounter, use run_akde_with_data instead
#You need to input the following
#name1: the name of individual one (make sure it's spelled exactly)
#name2: the name of individual two (make sure it's spelled exaclty)
#study: the study the pair is from
#date_encounter: the date of encounter as a lubridate datetime (use ymd_hms("2016-04-23 14:30:00") to create them)
run_akde_analysis = function(name1, name2, study, date_encounter){
  #folder_name <- str_c(name1, "_", name2, "_Homeranges",date_encounter, sep = "")
  folder_name <- str_c(name1, "_", name2, "_Homeranges", sep = "")
  #Change the folder_dir string to the folder where you want all your results to be stored
  #Do not add any slashes to the end of the string
  #You ONLY need to change this variable
  folder_dir <- "C:/Users/MINGXI SONG/Documents/Code_Mingxi"
  folder_path <- str_c(folder_dir, folder_name, sep = "/")
  dir.create(folder_path)
  
  id1 <- movement_data_id_all4 %>% filter(study.id == study, individual.local.identifier ==name1)
  id2 <- movement_data_id_all4 %>% filter(study.id == study,individual.local.identifier == name2)
  
  before <- list(filter(id1, timestamp < date_encounter),
                 filter(id2, timestamp < date_encounter)) %>% map(., as.telemetry)
  
  after <- list(filter(id1, timestamp > date_encounter),
                filter(id2, timestamp > date_encounter)) %>% map(., as.telemetry)
  
  save(before, file = str_c(folder_path, "/", folder_name, "_before_data.rda"))
  save(after, file = str_c(folder_path, "/", folder_name, "_after_data.rda"))
  
  #Sets the projections to be the same as before. Change both to after is you want them standardized to after
  ctmm:::projection(before[1]) <- median(before)   #correct the projection --Qianru
  ctmm:::projection(before[2]) <- median(before)   #correct the projection --Qianru
  ctmm:::projection(after[1]) <- median(before)   #correct the projection --Qianru
  ctmm:::projection(after[2]) <- median(before)   #correct the projection --Qianru
  
  #Print the number of entries in each split dataset, gives you a sense of scale
  print(map_dbl(before, nrow))
  print(map_dbl(after, nrow))
  
  akde_all <- function(data){
    guess <- map(data, ctmm.guess, interactive = FALSE)
    fit <- map2(data, guess, ctmm.select)
    #akde_final <- map2(data, fit, akde)
    return(fit)
    #, akde = akde_final))
  }
  
  seq_tele <- list(before = before, after = after)
  model_list <- foreach ( i = seq_along(seq_tele), .packages = c("ctmm","lubridate","tidyverse")) %dopar% 
    {akde_all(seq_tele[[i]])}
  
  before_model_fit <- model_list[[1]]
  after_model_fit <- model_list[[2]]
  
  save(before_model_fit, file = str_c(folder_path, "/", folder_name, "_before_model.rda", sep = ""))
  save(after_model_fit, file = str_c(folder_path, "/", folder_name, "_after_model.rda", sep = ""))
  
  #fit overlap analysis table
  overlap_analysis_fit <- compare_overlap(before_model_fit, after_model_fit)
  png(filename = str_c(folder_path, "/", folder_name,"_overlap_table_fit.png", sep = ""),
      height = 50*nrow(overlap_analysis_fit), width = 200*ncol(overlap_analysis_fit))
  grid.table(overlap_analysis_fit)
  dev.off()
  #Save Overlap Table as dataframe
  save(overlap_analysis_fit, file = str_c(folder_path, "/", folder_name,"_overlap_table_fit.rda", sep = ""))
  
  #UDS overlap analysis table ---------------added by Qianru
  
  UDS_before <- akde(seq_tele[[1]],before_model_fit)
  UDS_after <- akde(seq_tele[[2]],after_model_fit)
  overlap_analysis_akde <- compare_overlap(UDS_before, UDS_after)
  png(filename = str_c(folder_path, "/", folder_name,"_overlap_table_UDS.png", sep = ""), 
      height = 50*nrow(overlap_analysis_akde), width = 200*ncol(overlap_analysis_akde))
  grid.table(overlap_analysis_akde)
  dev.off()
  
  #Save UDS Overlap Table as dataframe ---------------added by Qianru
  save(overlap_analysis_akde, file = str_c(folder_path, "/", folder_name,"_overlap_table_UDS.rda", sep = ""))
  
  pair_names <- c(name1, name2)
  #Make the effective home range table for reporting  ---------------edited by Qianru
  eff_home_range <- data.frame(x = c(UDS_before[[1]]$DOF.area[1], 
                                     UDS_after[[1]]$DOF.area[1]),
                               y = c(UDS_before[[2]]$DOF.area[1], 
                                     UDS_after[[2]]$DOF.area[1]))
  colnames(eff_home_range) <- pair_names
  rownames(eff_home_range) <- c("Before","After")
  
  png(filename = str_c(folder_path, "/", folder_name,"_effective_homerange_table.png", sep = ""), 
      height = 50*nrow(eff_home_range), width = 200*ncol(eff_home_range))
  grid.table(eff_home_range)
  dev.off()
  
  save(eff_home_range, file = str_c(folder_path, "/", folder_name,"_effective_homerange_table.rda", sep = ""))
  
  extents <- rbind(ctmm::extent(before_model_fit[[1]]), ctmm::extent(before_model_fit[[2]]),
                   ctmm::extent(UDS_before[[1]]), ctmm::extent(UDS_before[[2]]),
                   ctmm::extent(after_model_fit[[1]]), ctmm::extent(after_model_fit[[2]]),
                   ctmm::extent(UDS_after[[1]]), ctmm::extent(UDS_after[[2]]))
  
  
  x_range <- c((min(extents$x) - 500 %#% "meters"),(max(extents$x) + 500 %#% "meters"))
  y_range <- c((min(extents$y) - 500 %#% "meters"),(max(extents$y) + 500 %#% "meters"))
  
  sorted_names <- sort(c(name1, name2), decreasing = F)
  #Before encounter home range estimates
  png(filename = str_c(folder_path, "/", folder_name,"_before.png", sep = ""),
      width = 13, height = 6, units = "in", res = 75)
  plot(before, col = color(before, by = "individual"), UD = UDS_before,
       ylim = y_range, xlim = x_range, level.UD = 0.95,
       main = str_c(name1, " and ", name2, " before encounter: ",date_encounter))
  legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex =2.5) 
  
  dev.off()
  #After encounter home range estimates
  png(filename = str_c(folder_path, "/", folder_name,"_after.png", sep = ""), 
      width = 13, height = 6, units = "in", res = 75)
  plot(after, col = color(after, by = "individual"), UD = UDS_after, 
       ylim = y_range, xlim = x_range, 
       main = str_c(name1, " and ", name2, " after encounter: ",date_encounter))
  legend("bottomleft", legend = sorted_names, col = c("red","blue"), pch = 1, cex = 2.5) 
  dev.off()
  
}


#Examples-------------------------------

#make sure you have the movement_data_id_all4 dataset loaded in before trying these examples

#run_akde_analysis example

#make datetime to use as date_encounter
#Tex_Mirinha_encounter_date1 <- ymd_hms("2016-08-12 03:00:00")
#run_akde_analysis("Tex", "Mirinha", "Gustavo", Tex_Mirinha_encounter_date1)
BVOB_ORRR_encounter_date1 <- ymd_hms("2019-07-04 7:37:00")
run_akde_analysis("BVOB", "ORRR", "Bertreaux", BVOB_ORRR_encounter_date1)