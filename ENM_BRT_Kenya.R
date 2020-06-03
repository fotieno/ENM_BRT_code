library(maps)
library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(maptools)
library(mapdata)
library(gbm)
library(dismo)
library(ResourceSelection)
library(dplyr)
library(tidyverse)
library(SDMPlay)
library(spatstat)
library(xlsx)
library(MASS)
library(spatialEco)
library(pROC)
library(randomForest)
library(boot)
library(Hmisc)
library(verification)
library(pdp)
library(caret)
library(ROCR)
memory.limit(size=40000)
inputs<- ("C:/RF_data/Inputs") #data path
setwd(inputs)

#Read spatialdata
################
studarea_bnd <- readOGR("historicalCentWestSouthSubs_bnd.shp")
kenya1 <- readOGR("kenya1.shp")
projcrs <- crs(kenya1)
maskedbnd<- readOGR("historicalWestSouthSubsLess5kmbuff.shp")
plot(maskedbnd)

#read occurrences-prescence CSV
#-------------------------------
anthrax_occurrences <- read.csv("occurrence69Points.csv",header=TRUE, sep=",")
anthrax_occurrences <- anthrax_occurrences[,3:4]
head(anthrax_occurrences)

#accessing and stacking the output files 
#---------------------------------------

setwd("C:\\RF_data\\Inputs\\VIF\\BRT") #Tiffs path
tifFiles <- Sys.glob ('*.tif')
predictors_anthrax <-stack(tifFiles,quick
                           =F)

#loop run steps
##################
set.seed(100)
mylist<-list()
ivlist <- list()
testlist <- list()
trainlist <- list()
qs.list <- list()
nrep= 100
for (r in 1:nrep){
  setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2")
  nam <- print(paste("run", r))
  
  #generate random points
  #-------------------------
  anthrax_pseudos<-spsample(maskedbnd,n=69,"random")
  anthrax_pseudos_df<-data.frame(anthrax_pseudos)
  anthrax_pseudos_df<-anthrax_pseudos_df[,1:2]
  points(anthrax_pseudos)
  
  #training(75%) and Testing (25%) 
  #---------------------------------
  #Presence data
  anthrax_occurrences$species <- c("anthrax")
  trainIndex <- createDataPartition(anthrax_occurrences$species, p=0.75, list = FALSE)
  presence_Train <- anthrax_occurrences[trainIndex,]
  presence_Test <- anthrax_occurrences[-trainIndex,]
  presence_Train$species <- NULL
  presence_Test$species <- NULL
  
  
  #pseudo-absence data
  anthrax_pseudos_df$species <- c("anthrax")
  trainIndex_a <- createDataPartition(anthrax_pseudos_df$species, p=0.75, list = FALSE)
  abs_Train <- anthrax_pseudos_df[trainIndex_a,]
  abs_Test <- anthrax_pseudos_df[-trainIndex_a,]
  abs_Train$species <- NULL
  abs_Test$species <- NULL
  
  #extracting train presence and absence data
  #----------------------------------------------
  file1<-raster::extract(predictors_anthrax,presence_Train,method='simple')
  outcome<-rep(1,dim(file1)[[1]])
  file1<-cbind(outcome,file1)
  
  file2<-raster::extract(predictors_anthrax,abs_Train,method='simple')
  outcome<-rep(0,dim(file2)[[1]])
  file2<-cbind(outcome,file2)
  
  train_data<-as.data.frame(rbind(file1,file2))
  train_data<-na.omit(train_data)
  train_data2<-train_data
  colnames(train_data2) <- c('outcome','Organic_carbon','Cattle_density','Calcic_Vertisols','Vegetation_index','Haplic_Calcisols','Haplic_Vertisols','Relative_humidity','Soil_moisture','Rain_wettest_month','Temp_seasonality','Soil_clay','Longest_dry_season','Soil_pH','Drought_severity','Evapotranspiration','Silt','Slope','Soil_texture')
  getwd()
  write.csv(train_data2,paste0("train","_", r, ".csv"),row.names=TRUE,col.names=FALSE)
  trainlist[[r]] <- train_data2
  
  
  #extracting test presence and absence data
  #----------------------------------------------
  pre_testdata<-raster::extract(predictors_anthrax,presence_Test,method='simple')
  outcome<-rep(1,dim(pre_testdata)[[1]])
  pre_testdata<-cbind(outcome,pre_testdata)
  
  abs_testdata<-raster::extract(predictors_anthrax,abs_Test,method='simple')
  outcome<-rep(0,dim(abs_testdata)[[1]])
  abs_testdata<-cbind(outcome,abs_testdata)
  
  test_data<-as.data.frame(rbind(pre_testdata,abs_testdata))
  test_data<-na.omit(test_data)
  test_data2 <- test_data
  colnames(test_data2) <- c('outcome','Organic_carbon','Cattle_density','Calcic_Vertisols','Vegetation_index','Haplic_Calcisols','Haplic_Vertisols','Relative_humidity','Soil_moisture','Rain_wettest_month','Temp_seasonality','Soil_clay','Longest_dry_season','Soil_pH','Drought_severity','Evapotranspiration','Silt','Slope','Soil_texture')
  write.csv(test_data2,paste0("test","_", r, ".csv"),row.names=TRUE,col.names=FALSE)
  testlist[[r]] <- test_data2
  
  #fitting model and selecting variables
  #######################################
  col<-ncol(train_data)
  brt_step2 <-gbm.step(data = train_data2, gbm.x = c(2:col), gbm.y = 1,family = "bernoulli", tree.complexity = 5,learning.rate = 0.001, bag.fraction = 0.5,max.trees = 2500,n.folds = 10)
 
  #obtaining and plotting variable influence
  #--------------------------------------
  ivlist[[r]]<-data.frame(summary(brt_step2))

  
  #storing brt objects
  #-----------------------------
  brt_step.bs <- try(gbm.step(data=train_data2[sample(NROW(train_data), NROW(train_data), replace=T),], gbm.x =c(2:col),
                              gbm.y = 1, family = "bernoulli", tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.5,max.trees = 2500,n.folds = 10,
                              verbose=TRUE, silent=FALSE, plot.main=TRUE))
  if (class(brt_step.bs) == "try-error") next
  qs.list[[r]] <- brt_step.bs
  cat("This is replicate number ", r, "\n")
  save(qs.list, file="qs.list.Rdata")
  rm(brt_step.bs)
  
  #storing run AUCs
  #----------------
  AUC<-brt_step2$cv.statistics$discrimination.mean
  mylist<-append(mylist,AUC)
}
#create CV AUC excel file
#--------------------------
write.xlsx(unlist(mylist),"C:\\RF_data\\Inputs\\VIF\\runs_brt2\\AUC.xlsx")


#Variable Influence--------------
#Individual run variable influences
ibv <- list()
for (r in 1:length(ivlist)) {
  ivl <- ggplot(ivlist[[r]], aes(x = reorder(ivlist[[r]][["var"]], ivlist[[r]][["rel.inf"]]), y = ivlist[[r]][["rel.inf"]])) + 
    geom_bar(stat = "identity",color="white", fill="dodgerblue3") +
    coord_flip() +
    ylab("Relative Influence") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=10, colour="black"),
          axis.text.y = element_text(size=10, colour="black"),
          axis.title.x = element_text(size=11),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5, vjust = -1))
  ibv[[r]] <- ivl
  plot(ibv[[r]])
}


#Averaged run variable influences
iv <- NULL
for (r in 1:length(ivlist)) {
  temp_iv <- data.frame(ivlist[[r]][["var"]],ivlist[[r]][["rel.inf"]], col = rep(r:r, each = 18)) 
  iv <- rbind(iv, temp_iv)
  
}
names(iv) <- c("var", "rel.inf", "Run")
a_iv <- iv %>% group_by(var) %>% summarise(mean_rel.inf = mean(rel.inf), sd_rel.inf = sd(rel.inf)) 
names(a_iv) <- c("var", "mean_rel_inf", "sd_rel_inf")
a_iv$var <- c('Calcic Vertisols','Cattle Density','Drought Severity', 'Evapotranspiration', 'Haplic Calcisols', 'Haplic Vertisols', 'Longest Dry Season','Soil Organic Carbon', 'Rain Wettest Month', 'Relative Humidity','Silt', 'Slope', 'Soil Clay', 'Soil Moisture', 'Soil pH', 'Soil Texture','Temperature Seasonality', 'Enhanced Vegetation Index')
ivy <- ggplot(data=a_iv, aes(x=reorder(var, mean_rel_inf), y=mean_rel_inf)) +
  geom_bar(stat = "identity",color="white", fill="dodgerblue3") +
  geom_errorbar(aes(ymin=mean_rel_inf-sd_rel_inf, ymax=mean_rel_inf+sd_rel_inf), width=.2,
                position=position_dodge(.9),size=0.7, color="black") +
  coord_flip() +
  ylab("Relative Influence") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=11, colour="black"),
        axis.text.y = element_text(size=11, colour="black"),
        axis.title.x = element_text(size=11),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = -1))
plot(ivy)
getwd()
tiff('influence_plot.tif',units="in", width=6, height=6, res=300)
plot(ivy)
dev.off()

#BRT partial plots with CI--------------------------
#Code Layout for PDPs
#1. PDP for each run for a variable
#2. Combined PDPS for all runs for a variable
#3. Averaged PDPs for all runs for a variable

#The bottom of this code section is grid arrange for all averaged PDPs for each variable


#Cattle_Density
yty <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt <- partial(qs.list[[r]], pred.var = "Cattle_density", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
            plot = FALSE, train = df)
  rb <- ggplot(rt, aes(x = Cattle_density, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=6, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty[[r]] <- rb
  plot(yty[[r]])
}

ry <- NULL
temp_ry <- NULL
for(r in 1:length(yty)) {
    temp_ry <- data.frame(x = yty[[r]][["data"]][["Cattle_density"]], y=yty[[r]][["data"]][["yhat"]], 
                          col = rep(r:r, length.out = NA))
    ry <- rbind(ry, temp_ry) 
}
names(ry) <- c("Cattle_Density", "P_Probability", "Run")
  by <- ggplot(ry, aes(x = Cattle_Density,y = P_Probability, group=Run, color = factor(Run))) + 
        geom_smooth() +
        guides(color=guide_legend(title = "Run")) +  ggtitle("Cattle Density")+
        xlab("animals/90km^2") +
        ylab("P(Probability)") +
        theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  plot(by)

ry_averaged <- ry
ry_averaged$Run <- 1
by_av <- ggplot(ry_averaged, aes(x=Cattle_Density,y=P_Probability)) + 
  geom_smooth() + ggtitle("Cattle Density")+
  xlab("animals/90km^2") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av)

#Organic_carbon
yty_1 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_1 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_1 <- partial(qs.list[[r]], pred.var = "Organic_carbon", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                plot = FALSE, train = df_1)
  rb_1 <- ggplot(rt_1, aes(x = Organic_carbon, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_1[[r]] <- rb_1
  plot(yty_1[[r]])
}

ry_1 <- NULL
for(r in 1:length(yty_1)) {
  temp_ry_1 <- data.frame(x = yty_1[[r]][["data"]][["Organic_carbon"]], y=yty_1[[r]][["data"]][["yhat"]], 
                          col = rep(r:r, length.out = NA))
  ry_1 <- rbind(ry_1, temp_ry_1) 
}
names(ry_1) <- c("Organic_carbon", "P_Probability", "Run")
by_1 <- ggplot(ry_1, aes(x = Organic_carbon,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) +  ggtitle("Soil Organic Carbon")+
  xlab("kg/m^3") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_1)

ry_averaged_1 <- ry_1
ry_averaged_1$Run <- 1
by_av_1 <- ggplot(ry_averaged_1, aes(x=Organic_carbon,y=P_Probability)) + 
  geom_smooth() + ggtitle("Soil Organic Carbon")+
  xlab("kg/m^3") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_1)


#Calcic_Vertisols
yty_2 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_2 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_2 <- partial(qs.list[[r]], pred.var = "Calcic_Vertisols", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_2)
  rb_2 <- ggplot(rt_2, aes(x = Calcic_Vertisols, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_2[[r]] <- rb_2
  plot(yty_2[[r]])
}

ry_2 <- NULL
for(r in 1:length(yty_2)) {
  temp_ry_2 <- data.frame(x = yty_2[[r]][["data"]][["Calcic_Vertisols"]], y=yty_2[[r]][["data"]][["yhat"]], 
                          col = rep(r:r, length.out = NA))
  ry_2 <- rbind(ry_2, temp_ry_2) 
}
names(ry_2) <- c("Calcic_Vertisols", "P_Probability", "Run")
by_2 <- ggplot(ry_2, aes(x = Calcic_Vertisols,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Calcic Vertisols")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_2)

ry_averaged_2 <- ry_2
ry_averaged_2$Run <- 1
by_av_2 <- ggplot(ry_averaged_2, aes(x=Calcic_Vertisols,y=P_Probability)) + 
  geom_smooth() + ggtitle("Calcic Vertisols")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_2)


#Vegetation_index
yty_3 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_3 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_3 <- partial(qs.list[[r]], pred.var = "Vegetation_index", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_3)
  rb_3 <- ggplot(rt_3, aes(x = Vegetation_index, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_3[[r]] <- rb_3
  plot(yty_3[[r]])
}

ry_3 <- NULL
for(r in 1:length(yty_3)) {
  temp_ry_3 <- data.frame(x = yty_3[[r]][["data"]][["Vegetation_index"]], y=yty_3[[r]][["data"]][["yhat"]], 
                          col = rep(r:r, length.out = NA))
  ry_3 <- rbind(ry_3, temp_ry_3) 
}
names(ry_3) <- c("Vegetation_index", "P_Probability", "Run")
by_3 <- ggplot(ry_3, aes(x = Vegetation_index,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Enhanced Vegetation Index")+
  xlab("index") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_3)

ry_averaged_3 <- ry_3
ry_averaged_3$Run <- 1
by_av_3 <- ggplot(ry_averaged_3, aes(x=Vegetation_index,y=P_Probability)) + 
  geom_smooth() + ggtitle("Enhanced Vegetation Index")+
  xlab("index") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_3)


#Haplic_Calcisols
yty_4 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_4 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_4 <- partial(qs.list[[r]], pred.var = "Haplic_Calcisols", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_4)
  rb_4 <- ggplot(rt_4, aes(x = Haplic_Calcisols, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_4[[r]] <- rb_4
  plot(yty_4[[r]])
}

ry_4 <- NULL
for(r in 1:length(yty_4)) {
  temp_ry_4 <- data.frame(x = yty_4[[r]][["data"]][["Haplic_Calcisols"]], y=yty_4[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_4 <- rbind(ry_4, temp_ry_4) 
}
names(ry_4) <- c("Haplic_Calcisols", "P_Probability", "Run")
by_4 <- ggplot(ry_4, aes(x = Haplic_Calcisols,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Haplic Calcisols")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_4)

ry_averaged_4 <- ry_4
ry_averaged_4$Run <- 1
by_av_4 <- ggplot(ry_averaged_4, aes(x=Haplic_Calcisols,y=P_Probability)) + 
  geom_smooth() + ggtitle("Haplic Calcisols")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_4)


#Haplic_Vertisols
yty_5 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_5 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_5 <- partial(qs.list[[r]], pred.var = "Haplic_Vertisols", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_5)
  rb_5 <- ggplot(rt_5, aes(x = Haplic_Vertisols, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_5[[r]] <- rb_5
  plot(yty_5[[r]])
}

ry_5 <- NULL
for(r in 1:length(yty_5)) {
  temp_ry_5 <- data.frame(x = yty_5[[r]][["data"]][["Haplic_Vertisols"]], y=yty_5[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_5 <- rbind(ry_5, temp_ry_5) 
}
names(ry_5) <- c("Haplic_Vertisols", "P_Probability", "Run")
by_5 <- ggplot(ry_5, aes(x = Haplic_Vertisols,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Haplic Vertisols")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_5)

ry_averaged_5 <- ry_5
ry_averaged_5$Run <- 1
by_av_5 <- ggplot(ry_averaged_5, aes(x=Haplic_Vertisols,y=P_Probability)) + 
  geom_smooth() + ggtitle("Haplic Vertisols")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_5)

#Relative_humidity
yty_6 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_6 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_6 <- partial(qs.list[[r]], pred.var = "Relative_humidity", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_6)
  rb_6 <- ggplot(rt_6, aes(x = Relative_humidity, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_6[[r]] <- rb_6
  plot(yty_6[[r]])
}

ry_6 <- NULL
for(r in 1:length(yty_6)) {
  temp_ry_6 <- data.frame(x = yty_6[[r]][["data"]][["Relative_humidity"]], y=yty_6[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_6 <- rbind(ry_6, temp_ry_6) 
}
names(ry_6) <- c("Relative_humidity", "P_Probability", "Run")
by_6 <- ggplot(ry_6, aes(x = Relative_humidity,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Relative Humidity")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_6)

ry_averaged_6 <- ry_6
ry_averaged_6$Run <- 1
by_av_6 <- ggplot(ry_averaged_6, aes(x=Relative_humidity,y=P_Probability)) + 
  geom_smooth() + ggtitle("Relative Humidity")+
  xlab("%") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_6)


#Soil_moisture
yty_7 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_7 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_7 <- partial(qs.list[[r]], pred.var = "Soil_moisture", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_7)
  rb_7 <- ggplot(rt_7, aes(x = Soil_moisture, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_7[[r]] <- rb_7
  plot(yty_7[[r]])
}

ry_7 <- NULL
for(r in 1:length(yty_7)) {
  temp_ry_7 <- data.frame(x = yty_7[[r]][["data"]][["Soil_moisture"]], y=yty_7[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_7 <- rbind(ry_7, temp_ry_7) 
}
names(ry_7) <- c("Soil_moisture", "P_Probability", "Run")
by_7 <- ggplot(ry_7, aes(x = Soil_moisture,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Soil Moisture")+
  xlab("m^3/m^3") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_7)

ry_averaged_7 <- ry_7
ry_averaged_7$Run <- 1
by_av_7 <- ggplot(ry_averaged_7, aes(x=Soil_moisture,y=P_Probability)) + 
  geom_smooth() + ggtitle("Soil Moisture")+
  xlab("m^3/m^3") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_7)


#Rain_wettest_month
yty_8 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_8 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_8 <- partial(qs.list[[r]], pred.var = "Rain_wettest_month", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_8)
  rb_8 <- ggplot(rt_8, aes(x = Rain_wettest_month, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_8[[r]] <- rb_8
  plot(yty_8[[r]])
}

ry_8 <- NULL
for(r in 1:length(yty_8)) {
  temp_ry_8 <- data.frame(x = yty_8[[r]][["data"]][["Rain_wettest_month"]], y=yty_8[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_8 <- rbind(ry_8, temp_ry_8) 
}
names(ry_8) <- c("Rain_wettest_month", "P_Probability", "Run")
by_8 <- ggplot(ry_8, aes(x = Rain_wettest_month,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() + 
  guides(color=guide_legend(title = "Run")) + ggtitle("Rainfall Wettest Month")+
  xlab("mm") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_8)

ry_averaged_8 <- ry_8
ry_averaged_8$Run <- 1
by_av_8 <- ggplot(ry_averaged_8, aes(x=Rain_wettest_month,y=P_Probability)) + 
  geom_smooth() + ggtitle("Rainfall Wettest Month")+
  xlab("mm") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_8)


#Temp_seasonality
yty_9 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_9 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_9 <- partial(qs.list[[r]], pred.var = "Temp_seasonality", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_9)
  rb_9 <- ggplot(rt_9, aes(x = Temp_seasonality, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_9[[r]] <- rb_9
  plot(yty_9[[r]])
}

ry_9 <- NULL
for(r in 1:length(yty_9)) {
  temp_ry_9 <- data.frame(x = yty_9[[r]][["data"]][["Temp_seasonality"]], y=yty_9[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_9 <- rbind(ry_9, temp_ry_9) 
}
names(ry_9) <- c("Temp_seasonality", "P_Probability", "Run")
by_9 <- ggplot(ry_9, aes(x = Temp_seasonality,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Temp Seasonality")+
  xlab("°C") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_9)

ry_averaged_9 <- ry_9
ry_averaged_9$Run <- 1
by_av_9 <- ggplot(ry_averaged_9, aes(x=Temp_seasonality,y=P_Probability)) + 
  geom_smooth() + ggtitle("Temp Seasonality")+
  xlab("°C") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_9)



#Soil_clay
yty_10 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_10 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_10 <- partial(qs.list[[r]], pred.var = "Soil_clay", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                  plot = FALSE, train = df_10)
  rb_10 <- ggplot(rt_10, aes(x = Soil_clay, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_10[[r]] <- rb_10
  plot(yty_10[[r]])
}

ry_10 <- NULL
for(r in 1:length(yty_10)) {
  temp_ry_10 <- data.frame(x = yty_10[[r]][["data"]][["Soil_clay"]], y=yty_10[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_10 <- rbind(ry_10, temp_ry_10) 
}
names(ry_10) <- c("Soil_clay", "P_Probability", "Run")
by_10 <- ggplot(ry_10, aes(x = Soil_clay,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Soil Clay Content")+
  xlab("mass fraction(%)") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_10)

ry_averaged_10 <- ry_10
ry_averaged_10$Run <- 1
by_av_10 <- ggplot(ry_averaged_10, aes(x=Soil_clay,y=P_Probability)) + 
  geom_smooth() + ggtitle("Soil Clay Content")+
  xlab("mass fraction(%)") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_10)



#Longest_dry_season
yty_11 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_11 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_11 <- partial(qs.list[[r]], pred.var = "Longest_dry_season", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                   plot = FALSE, train = df_11)
  rb_11 <- ggplot(rt_11, aes(x = Longest_dry_season, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_11[[r]] <- rb_11
  plot(yty_11[[r]])
}

ry_11 <- NULL
for(r in 1:length(yty_11)) {
  temp_ry_11 <- data.frame(x = yty_11[[r]][["data"]][["Longest_dry_season"]], y=yty_11[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_11 <- rbind(ry_11, temp_ry_11) 
}
names(ry_11) <- c("Longest_dry_season", "P_Probability", "Run")
by_11 <- ggplot(ry_11, aes(x = Longest_dry_season,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Longest Dry Season")+
  xlab("months") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_11)

ry_averaged_11 <- ry_11
ry_averaged_11$Run <- 1
by_av_11 <- ggplot(ry_averaged_11, aes(x=Longest_dry_season,y=P_Probability)) + 
  geom_smooth() + ggtitle("Longest Dry Season")+
  xlab("months") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_11)


#Soil_pH
yty_12 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_12 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_12 <- partial(qs.list[[r]], pred.var = "Soil_pH", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                   plot = FALSE, train = df_12)
  rb_12 <- ggplot(rt_12, aes(x = Soil_pH, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_12[[r]] <- rb_12
  plot(yty_12[[r]])
}

ry_12 <- NULL
for(r in 1:length(yty_12)) {
  temp_ry_12 <- data.frame(x = yty_12[[r]][["data"]][["Soil_pH"]], y=yty_12[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_12 <- rbind(ry_12, temp_ry_12) 
}
names(ry_12) <- c("Soil_pH", "P_Probability", "Run")
by_12 <- ggplot(ry_12, aes(x = Soil_pH,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() + ggtitle("Soil pH") +
  guides(color=guide_legend(title = "Run")) + 
  xlab("pH") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_12)

ry_averaged_12 <- ry_12
ry_averaged_12$Run <- 1
by_av_12 <- ggplot(ry_averaged_12, aes(x=Soil_pH,y=P_Probability)) + 
  geom_smooth() + ggtitle("Soil pH") +
  xlab("pH") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_12)


#Drought_severity
yty_13 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_13 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_13 <- partial(qs.list[[r]], pred.var = "Drought_severity", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                   plot = FALSE, train = df_13)
  rb_13 <- ggplot(rt_13, aes(x = Drought_severity, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_13[[r]] <- rb_13
  plot(yty_13[[r]])
}

ry_13 <- NULL
for(r in 1:length(yty_13)) {
  temp_ry_13 <- data.frame(x = yty_13[[r]][["data"]][["Drought_severity"]], y=yty_13[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_13 <- rbind(ry_13, temp_ry_13) 
}
names(ry_13) <- c("Drought_severity", "P_Probability", "Run")
by_13 <- ggplot(ry_13, aes(x = Drought_severity,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Drought Severity")+
  xlab("index") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_13)

ry_averaged_13 <- ry_13
ry_averaged_13$Run <- 1
by_av_13 <- ggplot(ry_averaged_13, aes(x=Drought_severity,y=P_Probability)) + 
  geom_smooth() +
  xlab("index") + ggtitle("Drought Severity")+
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_13)


#Evapotranspiration
yty_14 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_14 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_14 <- partial(qs.list[[r]], pred.var = "Evapotranspiration", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                   plot = FALSE, train = df_14)
  rb_14 <- ggplot(rt_14, aes(x = Evapotranspiration, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_14[[r]] <- rb_14
  plot(yty_14[[r]])
}

ry_14 <- NULL
for(r in 1:length(yty_14)) {
  temp_ry_14 <- data.frame(x = yty_14[[r]][["data"]][["Evapotranspiration"]], y=yty_14[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_14 <- rbind(ry_14, temp_ry_14) 
}
names(ry_14) <- c("Evapotranspiration", "P_Probability", "Run")
by_14 <- ggplot(ry_14, aes(x = Evapotranspiration,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Evapotranspiration")+
  xlab("mm") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_14)

ry_averaged_14 <- ry_14
ry_averaged_14$Run <- 1
by_av_14 <- ggplot(ry_averaged_14, aes(x=Evapotranspiration,y=P_Probability)) + 
  geom_smooth() +
  xlab("mm") + ggtitle("Evapotranspiration")+
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_14)


#Silt
yty_15 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_15 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_15 <- partial(qs.list[[r]], pred.var = "Silt", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                   plot = FALSE, train = df_15)
  rb_15 <- ggplot(rt_15, aes(x = Silt, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_15[[r]] <- rb_15
  plot(yty_15[[r]])
}

ry_15 <- NULL
for(r in 1:length(yty_15)) {
  temp_ry_15 <- data.frame(x = yty_15[[r]][["data"]][["Silt"]], y=yty_15[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_15 <- rbind(ry_15, temp_ry_15) 
}
names(ry_15) <- c("Silt", "P_Probability", "Run")
by_15 <- ggplot(ry_15, aes(x = Silt,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Silt")+
  xlab("mass fraction(%)") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_15)

ry_averaged_15 <- ry_15
ry_averaged_15$Run <- 1
by_av_15 <- ggplot(ry_averaged_15, aes(x=Silt,y=P_Probability)) + 
  geom_smooth() +
  xlab("mass fraction(%)") + ggtitle("Silt")+
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_15)


#Slope
yty_16 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_16 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_16 <- partial(qs.list[[r]], pred.var = "Slope", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                   plot = FALSE, train = df_16)
  rb_16 <- ggplot(rt_16, aes(x = Slope, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_16[[r]] <- rb_16
  plot(yty_16[[r]])
}

ry_16 <- NULL
for(r in 1:length(yty_16)) {
  temp_ry_16 <- data.frame(x = yty_16[[r]][["data"]][["Slope"]], y=yty_16[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_16 <- rbind(ry_16, temp_ry_16) 
}
names(ry_16) <- c("Slope", "P_Probability", "Run")
by_16 <- ggplot(ry_16, aes(x = Slope,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Slope")+
  xlab("degrees") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_16)

ry_averaged_16 <- ry_16
ry_averaged_16$Run <- 1
by_av_16 <- ggplot(ry_averaged_16, aes(x=Slope,y=P_Probability)) + 
  geom_smooth() + ggtitle("Slope")+
  xlab("degrees") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_av_16)


#Soil_texture
yty_17 <- list()
for (r in 1:length(qs.list)) {
  print(paste("pdp", r))
  df_17 <- data.frame(qs.list[[r]][["gbm.call"]][["dataframe"]])
  rt_17 <- partial(qs.list[[r]], pred.var = "Soil_texture", n.trees = qs.list[[r]][['n.trees']], recursive = FALSE, ice = TRUE, center = TRUE, na.rm = TRUE,
                   plot = FALSE, train = df_17)
  rb_17 <- ggplot(rt_17, aes(x = Soil_texture, y = yhat)) + 
    geom_smooth() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.position="none",
          axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour = "black"),
          axis.text.x = element_text(size=7, colour="black"),
          axis.text.y = element_text(size=7, colour="black"),
          axis.title.x = element_text(size=7),
          axis.title.y = element_text(size=7),
          plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
  yty_17[[r]] <- rb_17
  plot(yty_17[[r]])
}

ry_17 <- NULL
for(r in 1:length(yty_17)) {
  temp_ry_17 <- data.frame(x = yty_17[[r]][["data"]][["Soil_texture"]], y=yty_17[[r]][["data"]][["yhat"]], col = rep(r:r, length.out = NA))
  ry_17 <- rbind(ry_17, temp_ry_17) 
}
names(ry_17) <- c("Soil_texture", "P_Probability", "Run")
by_17 <- ggplot(ry_17, aes(x = Soil_texture,y = P_Probability, group=Run, color = factor(Run))) + 
  geom_smooth() +
  guides(color=guide_legend(title = "Run")) + ggtitle("Soil Texture")+
  xlab("factor") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
plot(by_17)

ry_averaged_17 <- ry_17
ry_averaged_17$Run <- 1
by_av_17 <- ggplot(ry_averaged_17, aes(x=Soil_texture,y=P_Probability)) + 
  geom_smooth() + ggtitle("Soil Texture")+
  xlab("factor") +
  ylab("P(Probability)") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position="none",
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        axis.text.x = element_text(size=7, colour="black"),
        axis.text.y = element_text(size=7, colour="black"),
        axis.title.x = element_text(size=7),
        axis.title.y = element_text(size=7),
        plot.title = element_text(hjust = 0.5, vjust = -1,size=7))
        
plot(by_av_17)


#creating mean of the predictions and plot

#Landscape_predictions-----------------------------------------
predictors_anthrax2 <- predictors_anthrax
names(predictors_anthrax2) <- c('Organic_carbon','Cattle_density','Calcic_Vertisols','Vegetation_index','Haplic_Calcisols','Haplic_Vertisols','Relative_humidity','Soil_moisture','Rain_wettest_month','Temp_seasonality','Soil_clay','Longest_dry_season','Soil_pH','Drought_severity','Evapotranspiration','Silt','Slope','Soil_texture')
predictors_anthrax2

#Training Runs to predict onto the landscape
for (r in 1:length(qs.list)) {
  setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2")
  print(paste("pred", r))

}

#Mean of Predictions
setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2")
tifFiles <- Sys.glob('*.tif') 
predicted<-stack(tifFiles, quick=F)
mean_narm = function(x,...){mean(x,na.rm=TRUE)}
Predict_mean<- do.call(overlay, c(predicted, fun = mean_narm))
setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\mean")
writeRaster(Predict_mean,filename='mean_predictBRT2', format="GTiff", overwrite=TRUE)
meanPred<-raster("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\mean\\mean_predictBRT2.tif")
plot(meanPred)

#Upper 97.5% of predictions
setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\mean")
tifFiles <- Sys.glob('*.tif') 
predicted<-stack(tifFiles, quick=F)
uci_narm = function(x,...){quantile(x, probs = c(0.975),na.rm=TRUE)}
Predict_uci<- do.call(overlay, c(predicted, fun = uci_narm))
setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\CI")
writeRaster(Predict_uci,filename='uci_predictBRT2', format="GTiff", overwrite=TRUE)
uciPred<-raster("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\CI\\uci_predictBRT2.tif")
plot(uciPred)

#Lower 2.5% of predictions
setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\mean")
tifFiles <- Sys.glob('*.tif') 
predicted<-stack(tifFiles, quick=F)
lci_narm = function(x,...){quantile(x, probs = c(0.025),na.rm=TRUE)}
Predict_lci<- do.call(overlay, c(predicted, fun = lci_narm))
setwd("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\CI")
writeRaster(Predict_lci,filename='lci_predictBRT2', format="GTiff", overwrite=TRUE)
lciPred<-raster("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\CI\\lci_predictBRT2.tif")
plot(lciPred)



#Test_Data_Evaluation-------------------
#Model evaluation using testdata 


evallist<- list()
auclist <- list()
brt_p_2list <- list()
brt_predlist <- list()
for (r in 1:length(qs.list)) {
  print(paste("testpred", r))
  
  brt_prediction <- predict(qs.list[[r]], testlist[[r]], n.trees = qs.list[[r]]["n.trees"], type = "response")
  calc.deviance(obs=testlist[[r]]$outcome, pred=brt_prediction, calc.mean=TRUE)
  eval_data <- cbind(testlist[[r]]$outcome, brt_prediction)
  presences <- eval_data[eval_data[,1]==1, 2]
  abscenses <- eval_data[eval_data[,1]==0, 2]
  evaluation <- evaluate(p=presences, a=abscenses)
  evaluation
  t_AUC <- evaluation@auc
  auclist <- append(auclist, t_AUC)
  evallist[[r]] <- evaluation
  brt_predlist[[r]] <- brt_prediction
  brt_p_2 <- data.frame(brt_prediction)
  brt_p_2list[[r]] <- brt_p_2
  
}

write.xlsx(unlist(auclist), ("C:\\RF_data\\Inputs\\VIF\\runs_brt2\\AUC\\AUC.xlsx"))
