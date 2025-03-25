library(iterators)
library(parallel)
library(foreach)
library(doParallel)
registerDoParallel(4)
library(ctmm)
library(lubridate)
library(geosphere)
library(dplyr)

#overlap_range function
overlap_range<-function(start1,end1,start2,end2){
  f1<-max(start1,start2)
  f2<-min(end1,end2)
  do_overlap= f1<f2
  if(do_overlap){
    return(c(f1,f2))
  }else{
    return(NA)
  }
}
#distance function
distance_time<-function(pair1,pair2){
  distance<-c()
  #pair1<-Dog_ctmm$pair1
  #pair2<-Dog_ctmm$pair2
  start1<-pair1$timestamp%>%min()
  end1<-pair1$timestamp%>%max()
  start2<-pair2$timestamp%>%min()
  end2<-pair2$timestamp%>%max()
  range<-overlap_range(start1,end1,start2,end2)
  range_start<-range[1]%>%as.numeric()
  range_end<-range[2]%>%as.numeric()
  #if (is.na(range)==FALSE){
  t_range<-seq(range_start,range_end,600)
  Guess_pair1 <- ctmm.guess(pair1, interactive = FALSE)
  model_pair1 <- ctmm.select(pair1, Guess_pair1,trace=2)
  Guess_pair2 <- ctmm.guess(pair2, interactive = FALSE)
  model_pair2 <- ctmm.select(pair2, Guess_pair2,trace=2)
  s1 <- predict(model_pair1, t=t_range,data = pair1,complete=TRUE,res=1)
  s2 <- predict(model_pair2, t=t_range,data = pair2,complete=TRUE,res=1)
  for (i in 1:length(t_range)){
    distance[i]<-geosphere::distHaversine(c(s1$longitude[i],s1$latitude[i]), c(s2$longitude[i],s2$latitude[i]))
  }
  result_distance<-data.frame(s1$timestamp,distance)
  return (result_distance)
}

#load good pair data
#remember to edit the file path
good_pair<-read.csv("./Code_Qianru/Data/good_pair.csv")

#input the columns you plan to run this time
#change the name 'good_pair_?'
#change good_pair[c(?:?),]
good_pair_400=good_pair[c(300:400),]

#run to calculate the distance!
#change the name 'distance_good_pair_?'
#change the name 'good_pair_?'
library(dplyr)
distance_good_pair_400<-list()
for (i in 1:nrow(good_pair_400)){
  pair1<-subset(movement_data_id_all4,individual.local.identifier==good_pair_400$pair1[i])%>%as.telemetry()
  pair2<-subset(movement_data_id_all4,individual.local.identifier==good_pair_400$pair2[i])%>%as.telemetry()
  distance_good_pair_400[[i]] <- distance_time(pair1,pair2)
  save(distance_good_pair_400, file = "distance_good_pair_400.rda")
}


#next step: filter out the pairs with distance <100 meters
#change the name: 'result_?'
#change the name 'distance_good_pair_?'
#change the name 'good_pair_?'

result_400 = foreach(i=1:length(distance_good_pair_400),.combine=rbind,.packages = c("base","parallel","doParallel","foreach","iterators","dplyr"))%dopar% {
  result=distance_good_pair_400[[i]][distance_good_pair_400[[i]]$distance<=100 & sum(distance_good_pair_400[[i]]$distance<=100)<4,]%>%
    data.frame%>%
    slice_min(distance)%>%
    mutate(pair1=good_pair_400$pair1[i],pair2=good_pair_400$pair2[i],study=good_pair_400$study[i],species=good_pair_400$species[i])
  return(result)
}

#save the data as an csv file
#change the file path
write.csv(result_266,file="./Code_Qianru/Data/result_266.csv",row.names = FALSE,quote = FALSE)

