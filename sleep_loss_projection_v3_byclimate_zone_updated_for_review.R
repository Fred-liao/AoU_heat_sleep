#Code for Figure 5 and Supplement Table 12 and 13
#Projection of the Sleep Loss in the Future by Climate Zones
#Jiawen Liao
#USC, 2025-01
#For External Review Purpose

library(data.table)
library(readxl)
library(ggplot2)
library(lubridate)
library(MonteCarlo)
library(plyr)
library(ggpubr)
library(ggh4x)
#package for plot
library(gridExtra)
library(grid)

#packages for netCDF processing
library(RNetCDF)
library(ncdf4) # Network Common Data Form v4, for loading climate data
library(fields) # For plotting map-like images
library(maps) # For getting coastlines and such
library(sf) #for spatial join
library(tidyverse) 
library(stringi)
rm(list = ls())

# setwd("/Volumes/Macintosh HD - Data/AoU/SSP_CMIP6_predictions/grid_with_climate_zones/")

### Step 1: Preprocessing the CMIP6 Downscaled netCDF predictions and map them to each climate zones ####

#read in geo-processed data grid with climate zone
grid_with_CZ = fread("grid_0.25_cliamte_zone.csv")
table(grid_with_CZ$BA_Climate)
grid_with_CZ[BA_Climate %in% c("Hot-Humid","Hot-Dry"),BA_Climate2 := "Hot"]
grid_with_CZ[BA_Climate %in% c("Cold","Very Cold"),BA_Climate2 := "Cold"]
grid_with_CZ[BA_Climate %in% c("Mixed-Humid","Mixed-Dry"),BA_Climate2 := "Mixed"]
grid_with_CZ[BA_Climate %in% c("Marine"),BA_Climate2 := "Marine"]
grid_with_CZ2 = grid_with_CZ[!is.na(BA_Climate2),.(lon,lat,BA_Climate2)]

### Step 2: Read in monthly TNN anomaly netCDF files and process them ####
#data from https://climateknowledgeportal.worldbank.org/country/united-states/climate-data-projections

#a function to read in .nc file, merge with grid with CZ classificaiton and calculate CZ average tnn anomaly
inFile_list = list.files("../raw/",pattern = ".nc",full.names = T)
#function
CZ_tnn_calcualtion = function(inFile,month){
   
   #message("Processing ", inFile)
   if(stringi::stri_sub(stringi::stri_split_fixed(inFile,"/")[[1]][8], 1,7)  != "anomaly") {message("not correct file"); return(NULL)}
   
   
   #get SSP senario, prediction period 
   stat_type = stri_split_fixed(stri_split_fixed(inFile,"/")[[1]][8], "_")[[1]][5]
   ssp_senario = stri_split_fixed(stri_split_fixed(inFile,"/")[[1]][8], "_")[[1]][3]
   period_period = stri_sub(stri_split_fixed(stri_split_fixed(inFile,"/")[[1]][8], "_")[[1]][6], 1,9) 
   
   file.nc <- nc_open(inFile)
   #Air Temperature Data Shape: Lon/Lat/Julian Date, units are Kelvin 
   tnn_anomaly <- ncvar_get(file.nc, 'anomaly-tnn-monthly-mean') 
   #Lat/Lon Data Shape: Lon or Lat length
   lat <- ncvar_get(file.nc, 'lat') #721
   lon <- ncvar_get(file.nc, 'lon') #1440
   nc_close(file.nc)
   #Account for reversed latitude
   
   var_name = paste("tnn_anomaly",stat_type,ssp_senario,period_period,sep = "_")
   
   grid = expand.grid(lon,lat)
   dim(grid)
   df <- reshape2::melt(tnn_anomaly[,,month], na.rm = TRUE)
   df = df$value
   
   grid2 = cbind.data.frame(grid,df)
   grid2 = as.data.table(grid2)
   setnames(grid2, old = c('Var1','Var2',"df"), new = c("lon",'lat',"value"))
   grid2[,lat := as.numeric(lat)]
   grid2[,lon := as.numeric(lon)]
   #merge with CZ grid
   grid3 = merge.data.table(grid_with_CZ2, grid2, by = c("lon","lat"), all.x = T, all.y = F)
   #calculate climate zone mean
   CZ_cal = grid3[,.(value_CZ = mean(value)), by = "BA_Climate2"]
   CZ_cal[,stat := stat_type]
   CZ_cal[,ssp := ssp_senario]
   CZ_cal[,period := period_period]
   CZ_cal[,month := month]
   return(CZ_cal)
   
}


in_df = expand.grid(inFile_list,c(1:12))
in_df = as.data.table(in_df)
colnames(in_df) = c('inFile','month')
in_df[,inFile := as.character(inFile)]
#run the function to read in data
tnn_prediction_CZ = as.data.table(mdply(.data = in_df, .fun = CZ_tnn_calcualtion, .progress = "text" )) 

#filelocation = "/Users/Jiawen/Google Drive/Research/Postdoc-USC/AllofUS/Data/Projection/"

write.csv(tnn_prediction_CZ, file = paste(filelocation,"projected_tnn_monthly_anomlay_mean_by_Climate_Zones.csv", sep = "/"),row.names = F)

### Step 3: Summary of the TNN anomaly by SSP senarios and by climate zones ####

# data_al_by_CZ = fread("paste(filelocation,"projected_tnn_monthly_anomlay_mean_by_Climate_Zones.csv", sep = "/"))

data_al_by_CZ =tnn_prediction_CZ

data_al_by_CZ = data_al_by_CZ[,.(BA_Climate2,value_CZ,stat,ssp,period,month)]
data_al_by_CZ_wide = dcast.data.table(data_al_by_CZ, BA_Climate2 + ssp+period + month ~ stat, value.var = "value_CZ")
data_al_by_CZ_wide[,ssp := tstrsplit(ssp,"-")[3]]

data_al_by_CZ = data_al_by_CZ[,.(BA_Climate2,value_CZ,stat,ssp,period,month)]
data_al_by_CZ_wide = dcast.data.table(data_al_by_CZ, BA_Climate2 + ssp+period + month ~ stat, value.var = "value_CZ")
data_al_by_CZ_wide[,ssp := tstrsplit(ssp,"-")[3]]

#write.csv(data_al_by_CZ_wide, file = "/Users/Jiawen/Google Drive/Research/Postdoc-USC/AllofUS/Data/Projection/projected_tnn_annomaly_by_climate_zone_long.csv")


#Read in NTA and sleep duration reduction Exposure-response information, stratified by climate zones and month
NTA_stratified = read_excel("All_results_v4_cz_month_stra.xlsx"
                            ,sheet = "NTA_zone_month_stratify")
NTA_stratified = as.data.table(NTA_stratified)
NTA_stratified[,`...1` := NULL]
NTA_stratified_sub = NTA_stratified[stratify_var  == 'month']
NTA_stratified_sub[, month_levels := tstrsplit(exp,"month)")[2]]
#NTA_stratified_sub[,month_levels2 := NULL]
NTA_stratified_sub[,month_levels2 := lubridate::month(as.numeric(month_levels), label = TRUE, abbr = TRUE)]
NTA_stratified_sub = NTA_stratified_sub[c_variable_select  == "tmmn_anomaly"]
NTA_stratified_sub = NTA_stratified_sub[,.(climate_zone_selected, month_levels2,month_levels,exp_var,outcome,Estimate, `Std. Error`,  LCI,  UCI,`Pr(>|t|)`)]

#merge NTA in future with NTA - sleep associations with selected outcomes
NTA_stratified_selected_outcome = NTA_stratified_sub[outcome == "minute_asleep"]
#by zone, change to wide data
unique(NTA_stratified_selected_outcome$climate_zone_selected)

#change to wide data
NTA_stratified_sub2 = copy(NTA_stratified_sub)

NTA_stratified_sub2[,LCI := Estimate - 1.96*`Std. Error`]
NTA_stratified_sub2[,UCI := Estimate + 1.96*`Std. Error`]

#only total sleep time
NTA_cz_month_wide = dcast.data.table(NTA_stratified_sub2[exp_var == "tmmn_anomaly" & outcome == "minute_asleep"], month_levels + month_levels2 + climate_zone_selected ~ outcome,
                                     value.var = c("Estimate","LCI","UCI","Pr(>|t|)") ) 
NTA_cz_month_wide[,month_levels := as.numeric(month_levels)]
NTA_cz_month_wide = NTA_cz_month_wide[order(climate_zone_selected,month_levels)]

### Step 4: Merge exposure-response (ER) and predicted NTA (exposure) and Monte Carlo Simulation by SSP and Zones ####

data_all2 = merge.data.table(data_al_by_CZ_wide, NTA_cz_month_wide , by.x = c('month',"BA_Climate2"),
                               by.y = c("month_levels","climate_zone_selected")  )

data_all2 = data_all2[order(period,month)]

##Monte Carlo Simulation Function by climate zones and month
monte_carlo_sleep_loss = function(month,BA_Climate2,ssp,period,median,temp_se  , Estimate_minute_asleep , effect_SD ){
      NTA_simu = rnorm(n = 10000, mean = median, sd = temp_se)
      Effect_simu = rnorm(n = 10000, mean = Estimate_minute_asleep, sd = effect_SD)
      NTA_effect_simu = NTA_simu*Effect_simu
      effect_median = median(NTA_effect_simu)
      effect_p25 = as.numeric(quantile(NTA_effect_simu,0.025))
      effect_p975 = as.numeric(quantile(NTA_effect_simu,0.975))
      temp = cbind.data.frame(effect_median,effect_p25,effect_p975)
      return(temp)
      
      
}


#Runing the MC Function
SSP_simu_all2 = as.data.table(mdply(.data = data_all2_sub, .fun =monte_carlo_sleep_loss, .progress = "text" ))

#calcualte average sleep loss, total time in hour
SSP_simu_all2_yearly = SSP_simu_all2[,.(effect_median_yearly = mean(effect_median)*365.25/60 ,
                                        effect_p975_yearly = mean(effect_p975)*365.25/60 ,
                                        effect_p25_yearly = mean(effect_p25)*365.25/60 ), by = c("ssp","period","BA_Climate2")]

SSP_simu_all2_yearly = SSP_simu_all2_yearly[order(ssp,period)]


#reverse the number, change negative number to postive since its sleep duration loss
SSP_simu_all2_yearly[,effect_median_yearly_rev := -effect_median_yearly]
SSP_simu_all2_yearly[,effect_p25_yearly_rev := -effect_p25_yearly]
SSP_simu_all2_yearly[,effect_p975_yearly_rev := -effect_p975_yearly]

#plot yearly results
year_plot = ggplot(SSP_simu_all2_yearly, aes(x =period , y = effect_median_yearly_rev) )+
      geom_line(aes(color = ssp, group = ssp))+
      geom_point(aes(color = ssp))+
      scale_color_manual(values = c("#4393c3","#f4a582","#d6604d","#b2182b") )+
      scale_fill_manual(values = c("#4393c3","#f4a582","#d6604d","#b2182b"))+
      geom_hline(yintercept = 0, linetype = "dashed")+
      geom_ribbon(aes(ymin = effect_p975_yearly_rev, ymax =effect_p25_yearly_rev , group = ssp, fill = ssp), alpha = 0.1)+
      labs(x = "Year Period", y = "Sleep Duration\n Reduction by year (h)",fill = "SSP Senario",color = "SSP Senario")+
      theme_classic()+
      facet_grid(. ~ BA_Climate2 ) +
      theme(axis.title.x = element_text(size=12, angle=0, face="bold", vjust=1.1),
            axis.title.y = element_text(size=12, angle=90, face="bold", vjust=-0),
            axis.text.y = element_text(size=12, angle=0, face="bold"),
            axis.text.x = element_text(size=10, angle=90),
            strip.background = element_part_rect(side = "b"),
            strip.text = element_text(size = 10),
            strip.placement = "inside",
            axis.ticks.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            legend.text = element_text(size=10),
            legend.title = element_text(size=10,face = "bold"),
            #legend.title=element_blank(),
            text = element_text(family = "sans"),
            legend.position="none")

#plot yearly results
#plot monthly results, using minutes loss per day?
SSP_simu_all2[,effect_median_month_rev := -effect_median*30.4375 / 60]
SSP_simu_all2[,effect_p975_month_rev := -effect_p975*30.4375 / 60]
SSP_simu_all2[,effect_p25_month_rev := -effect_p25*30.4375 / 60]
SSP_simu_all3 = copy(SSP_simu_all2)


#monthly plot
month_plot = ggplot(SSP_simu_all3, aes(x =month , y = effect_median_month_rev) )+
      geom_line(aes(color = ssp, group = ssp))+
      geom_point(aes(color = ssp))+
      scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))+
      scale_color_manual(values = c("#4393c3","#f4a582","#d6604d","#b2182b") )+
      scale_fill_manual(values = c("#4393c3","#f4a582","#d6604d","#b2182b"))+
      geom_ribbon(aes(ymin = effect_p975_month_rev, ymax =effect_p25_month_rev , group = ssp, fill = ssp), alpha = 0.1)+
      geom_hline(yintercept = 0, linetype = "dashed")+
      labs(x = "Month", y = "Sleep Duration\n Reduction (h/month)",fill = "SSP Senario",color = "SSP Senario")+
      theme_classic()+
      facet_grid(period ~ BA_Climate2 ) +
      theme(axis.title.x = element_text(size=12, angle=0, face="bold", vjust=1.1),
            axis.title.y = element_text(size=12, angle=90, face="bold", vjust=-0),
            axis.text.y = element_text(size=12, angle=0, face="bold"),
            axis.text.x = element_text(size=10, angle=0),
            strip.background = element_part_rect(side = "b"),
            strip.text = element_text(size = 10),
            axis.ticks.y = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line.x = element_line(),
            axis.line.y = element_line(),
            legend.text = element_text(size=10),
            legend.title = element_text(size=10,face = "bold"),
            #legend.title=element_blank(),
            text = element_text(family = "sans"),
            legend.position="bottom")


y_m_plot = arrangeGrob(arrangeGrob(year_plot,
                        top = grid::textGrob("A)", x = 0, hjust = 0,gp = gpar(fontsize= 14,col  = "#951515",fontface  = "bold"))),
            arrangeGrob(month_plot,
                        top = grid::textGrob("B)", x = 0, hjust = 0,gp = gpar(fontsize= 14,col  = "#951515",fontface  = "bold")))
            ,nrow = 2 ,heights = c(1,2))


plot(y_m_plot)

#export the figure
#Figure 5
ggsave(y_m_plot, file = "fig5_prediction_v3.jpeg",
       device = "jpeg",
       width = 10,height = 9)

#version 3 is by climate zone
#supplement tables 12 and table 13
write.csv(SSP_simu_all3,file = "TST_loss_NTA_prediction_by_month_v3_1.csv")
write.csv(SSP_simu_all2_yearly, file = "TST_loss_NTA_prediction_by_year_v3_1.csv")



