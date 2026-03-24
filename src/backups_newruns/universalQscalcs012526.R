################Draft Universal I metrics. Christine Albano, 1/20/26

###################decisions to discuss:
################summary by water year -- do all agree?
###############threshold for degree days -- using 10C;  for dry days -- using 2 mm -- do all agree? may need to consider adjusting based on climatologies/how data are summarized spatially
#################smoothing -- using 5 day for GDD calcs -- do all agree?
#################climatological baseline -- 1996-2025 -- do all agree?
################## do cool season drought metrics make sense for SC???
################## calculate snow metrics for ID/SC?
##################how to calculate UTCI and wet bulb globe temperature
#################for large areas -- use spatial mean, min, median, or max? -- currently using mean but have other stats for NV areas
###################duration - climdex (# days with duration > 6 days) vs. annual max duration (simpler, easier)

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(zoo)
library(ggplot2)
library(trend)
'%!in%' <- function(x, y) !(x %in% y)

########## a few climate index resources
###https://www.ecad.eu//indicesextremes/indicesdictionary.php#2
###https://cran.r-project.org/web/packages/ClimInd/ClimInd.pdf
##https://rdrr.io/cran/ClimInd/f/README.md
###https://etccdi.pacificclimate.org/list_27_indices.shtml


setwd("C:\\Users\\albano\\OneDrive - Desert Research Institute\\W2L\\UniversalSurvey\\")

areas<-c("NV\\UpperHumboldt_","NV\\MiddleHumboldt_","NV\\LowerHumboldt_","SC\\Allendale_Box_","SC\\EPSCOR_AOI_Allendale_","ID\\Deary_", "ID\\Bovill_", "ID\\ElkRiver_","ID\\Kendrick_", "ID\\Juliaetta_", "ID\\Troy_")##"NV\\HumboldtHUC6_",
areas2<-c("UpperHumboldt_","MiddleHumboldt_","LowerHumboldt_","Allendale_Box_","EPSCOR_AOI_Allendale_","Deary_", "Bovill_", "ElkRiver_","Kendrick_", "Juliaetta_", "Troy_")##"HumboldtHUC6_",

dtasources<-c("gridmet","nclim",'era5land')

for (i in 1: length(areas)){
  for (j in 1: length(dtasources)){
    if(dtasources[j]=="nclim"){
      vars<-c("system.index","precip", "tmax", "tmin")}
    if(dtasources[j]=="gridmet"){
      vars<-c("system.index","bi" ,"erc","eto","etr","fm100", "fm1000",
                     "pr","rmax" ,"rmin", "sph","srad","th","tmmn", "tmmx","vpd","vs")}
    if(dtasources[j]=="era5land"){
      vars<-c('system.index','temperature_2m','dewpoint_temperature_2m','snow_depth_water_equivalent', 
              'snow_depth', 'snowfall_sum','snowmelt_sum','potential_evaporation_sum',
              'runoff_sum','u_component_of_wind_10m','v_component_of_wind_10m',
              'total_precipitation_sum','temperature_2m_min','temperature_2m_max',
              'surface_latent_heat_flux_sum','surface_net_thermal_radiation_sum',
              'surface_net_solar_radiation_sum','surface_sensible_heat_flux_sum',
              'surface_solar_radiation_downwards_sum','surface_thermal_radiation_downwards_sum')}
    
    #######select water years to analyze based on data source 
    if(dtasources[j]!="gridmet"){wateryearselect<-c(1952:2025)}
    if(dtasources[j]=="gridmet"){wateryearselect<-c(1980:2025)}

########################################################################################
#############################data prep##################################################
########################################################################################

#####################import daily data

  df<-read.csv(paste0("DailyData\\",areas[i],dtasources[j],".csv"))%>% select(all_of(vars))
  df$date<-ymd(substr(df$system.index, 1,8))
  df$year<-year(df$date)
  df$month<-month(df$date)
  df$wateryear<-ifelse(df$month>9, df$year+1, df$year)
  df$yday<-yday(df$date)
  df$wyday<-ifelse(df$year==df$wateryear & leap_year(df$year), df$yday+92, ifelse(df$year==df$wateryear, df$yday+91, ifelse(leap_year(df$wateryear-1), df$yday-274,df$yday-273)))
  df$mday<-mday(df$date)
  df$season2<-ifelse(df$month %in% 4:9, "growing", "cool")
  df$season4<-ifelse(df$month %in% 3:5, "MAM", ifelse(df$month %in% 6:8,"JJA", ifelse(df$month %in% 9:11, "SON", "DJF")))
  if(dtasources[j]=="nclim"){
  names(df)[2:4]<-c("pr","tmmx", "tmmn")
  }
  if(dtasources[j]=="era5land"){
    names(df)[c(8,12,13,14)]<-c("eto","pr","tmmn", "tmmx")
    df[,c(4:8,12)]<-df[,c(4:8,12)]*1000
  }
  if(dtasources[j]!="nclim"){
  df$tmmn<-df$tmmn-273.15
  df$tmmx<-df$tmmx-273.15}
  
  df$drydays<-ifelse(df$pr<=2,1,0)##dry days defined as p < 2mm
 


##################calculate quantiles
  q<-c(0,0.01, 0.05, 0.1,0.25,0.5,0.75,0.9,0.95, 0.99,1)
  quibble <- function(x, q = c(0,0.01, 0.05,0.1,0.25, 0.5, 0.75,0.9,0.95, 0.99,1), dropNA = TRUE) {
    tibble(x = quantile(x, q, na.rm = dropNA), q = q)
  }


############ calculate daily climos
  if(dtasources[j]=="nclim"){
    daily9625climoquants<- df %>% filter(wateryear %in% 1996:2025) %>% group_by(month, mday) %>%
      reframe(across(2:4, ~ quibble(.x, dropNA = TRUE))) %>%
      pivot_longer(3:5)}
  
  if(dtasources[j]=="gridmet"){
  daily9625climoquants<- df %>% filter(wateryear %in% 1996:2025) %>% group_by(month, mday) %>%
    reframe(across(2:17, ~ quibble(.x, dropNA = TRUE))) %>%
    pivot_longer(3:18)}
  
  if(dtasources[j]=="era5land"){
    daily9625climoquants<- df %>% filter(wateryear %in% 1996:2025) %>% group_by(month, mday) %>%
      reframe(across(2:20, ~ quibble(.x, dropNA = TRUE))) %>%
      pivot_longer(3:21)}
  
  daily9625climoquants$quantile<-daily9625climoquants$value$q
  
  daily9625climoquants$value<-daily9625climoquants$value$x
  
  daily9625climoquants$name<-paste0("dailyclimo9625",daily9625climoquants$name)
  
  dfdaily90climos<-left_join(df, daily9625climoquants %>% 
                               filter(quantile==0.9, name %in% paste0("dailyclimo9625",c("tmmx","tmmn","pr"))) %>% 
                               pivot_wider(names_from=name, values_from=value))


  
  ############ calculate monthly climos
  monthly9625climoquants<-df %>% filter(wateryear %in% 1996:2025) %>% 
    group_by(wateryear,month) %>% 
    summarize(pr=sum(pr), tmmn=mean(tmmn), tmmx=mean(tmmx)) %>% 
    group_by(month) %>% 
    reframe(across(2:4, ~ quibble(.x, dropNA = TRUE))) %>%
    pivot_longer(2:4)
  
  monthly9625climoquants$quantile<-monthly9625climoquants$value$q
  monthly9625climoquants$value<-unname(monthly9625climoquants$value$x)
  monthly9625climoquants$name<-paste0("monthlyclimo9625",monthly9625climoquants$name)
  
  # ggplot(monthly9625climoquants %>% 
  #        filter(name=="monthlyclimo9625pr", quantile==0.5), aes(x=month, y=value))+
  #        geom_point()+geom_smooth(se=FALSE)
  # 
  
  
  ############ calculate annual climos

  annual9625climoquants<-df %>% filter(wateryear %in% 1996:2025) %>% 
    group_by(wateryear) %>% 
    summarize(annpr=sum(pr), anntmmn=mean(tmmn), anntmmx=mean(tmmx)) %>% 
    reframe(across(2:4, ~ quibble(.x, dropNA = TRUE))) %>%
    pivot_longer(1:3)
  
  annual9625climoquants$quantile<-annual9625climoquants$value$q
  annual9625climoquants$value<-unname(annual9625climoquants$value$x)
  
  if(dtasources[j]=="era5land"){
  annual9625snow<-df %>% filter(wateryear %in% 1996:2025) %>% 
    group_by(wateryear) %>% 
    summarize(annsnowfall=sum(snowfall_sum),
              annswe=mean(snow_depth_water_equivalent),
              annsnowdepth=mean(snow_depth),
              annpeakswe=max(snow_depth_water_equivalent)) %>% 
    reframe(across(2:5, ~ quibble(.x, dropNA = TRUE)))%>%
    pivot_longer(1:4)
  
  annual9625snow$quantile<-as.numeric(annual9625snow$value$q)
  annual9625snow$value<-as.numeric(annual9625snow$value$x)
  annual9625climoquants<-rbind(annual9625climoquants,annual9625snow)
  }

##################daily smoothed data

### 5 day, centered
  df5day<-df
  if(dtasources[j]=="nclim"){
    df5day[,2:4]<-rollmean(df[,2:4],5, align="center", na.pad=TRUE)}
  if(dtasources[j]=="gridmet"){
    df5day[,2:17]<-rollmean(df[,2:17],5, align="center", na.pad=TRUE)}
  if(dtasources[j]=="era5land"){
    df5day[,2:20]<-rollmean(df[,2:20],5, align="center", na.pad=TRUE)}

### 10 day, centered
  df10day<-df
  if(dtasources[j]=="nclim"){
    df10day[,2:4]<-rollmean(df[,2:4],10, align="center", na.pad=TRUE)}
  if(dtasources[j]=="gridmet"){
    df10day[,2:17]<-rollmean(df[,2:17],10, align="center", na.pad=TRUE)}
  if(dtasources[j]=="era5land"){
    df10day[,2:20]<-rollmean(df[,2:20],10, align="center", na.pad=TRUE)}
  


##############################################################################################################
#TEMPERATURE#############################################################################################################
#############################################################################################################

###magnitude#########################################################################
####Based on your experience in this community since you have lived or worked here, has it gotten hotter or cooler overall?
####Since you’ve lived or worked here, has it gotten hotter or cooler during the summer months?
  


    ### mean annual water year temperature
    annualtemp<-df %>% filter(wateryear %in% wateryearselect) %>% group_by(wateryear) %>% summarize(meananntemp=mean((tmmn+tmmx)/2,na.rm=TRUE))
    annualtemp<-left_join(annualtemp, df %>% filter(wateryear %in% wateryearselect) %>% group_by(wateryear) %>% summarize(meanminanntemp=mean(tmmn,na.rm=TRUE)))
    annualtemp<-left_join(annualtemp,df %>% filter(wateryear %in% wateryearselect) %>% group_by(wateryear) %>% summarize(meanmaxanntemp=mean(tmmx,na.rm=TRUE)))
    
    
    ### mean summer temperature, June-Aug
    annualtemp<-left_join(annualtemp,df %>% filter(wateryear %in% wateryearselect, season4=="JJA") %>% group_by(wateryear) %>% summarize(maxminjjatemp=max(tmmn,na.rm=TRUE)))
    annualtemp<-left_join(annualtemp,df %>% filter(wateryear %in% wateryearselect, season4=="JJA") %>% group_by(wateryear) %>% summarize(maxmaxjjatemp=max(tmmx,na.rm=TRUE)))
    annualtemp<-left_join(annualtemp,df %>% filter(wateryear %in% wateryearselect, season4=="JJA") %>% group_by(wateryear) %>% summarize(meanjjatemp=mean((tmmn+tmmx)/2,na.rm=TRUE)))
    
    
    ### Wet bulb globe temperature (use cloud cover/solar rad to calculate); 
    #???
    
    
    ### UTCI; number of days > daily 90/95th pcntl; 
    
    # df$utci<-utci(
    #   (df$tmmn+df$tmmx)/2,
    #   (df$rmax+d$rmin)/2,
    #   df$ws,
    #   tmrt #sunshine duration??radiation temperature, Celsius
    # )

###frequency#############################################################################
####If you think it has gotten hotter or cooler during the summer, how frequently is it hotter?
    
    ### number of days with temp > answer to #8
    ## need P data to calculate this
    
    ####number of frost days (tmin<0)
    annualtemp<-left_join(annualtemp, df %>%
                            filter(wateryear %in% wateryearselect) %>%
                            group_by(wateryear) %>%
                            filter(tmmn< 0) %>%
                            summarize(frostdays=n()))
    
    ####number of summer days (tmax>25)
    annualtemp<-left_join(annualtemp, df %>% 
                            filter(wateryear %in% wateryearselect) %>% 
                            group_by(wateryear) %>% 
                            filter(tmmx>25) %>% 
                            summarize(summerdays=n())) 
    
    ####number of days with daily tmmx greater than 90ptile
    annualtemp<-left_join(annualtemp, dfdaily90climos %>% 
                            filter(wateryear %in% wateryearselect) %>% 
                            group_by(wateryear) %>% 
                            filter(tmmx>dailyclimo9625tmmx) %>% 
                            summarize(tmaxgt90thdays=n())) 
    
    ####number of days with daily tmmn greater than 90ptile
    annualtemp<-left_join(annualtemp, dfdaily90climos %>% 
                            filter(wateryear %in% wateryearselect) %>% 
                            group_by(wateryear) %>% 
                            filter(tmmn>dailyclimo9625tmmn) %>% 
                            summarize(tmingt90thdays=n())) 

###timing#################################################################################
###When does it start getting hotter in the spring?
###When does it start getting cooler in the fall?

    gdd<-df5day %>% 
      mutate(tmean=(tmmn+tmmx)/2) %>% 
      mutate(tmeangt10=ifelse(tmean>=10, tmean,0)) %>% 
      group_by(year) %>%
      mutate(gdd10c=cumsum(tmeangt10))
    
    
    ####start of growing season  -- 5 day mean > 10 C
    annualtemp<-left_join(annualtemp,gdd %>% filter(year %in% wateryearselect) %>% group_by(year) %>% filter(tmeangt10 > 10) %>%
      slice_head(n = 1) %>%
      ungroup() %>% select (wateryear,gsstartdoy=yday))
    
    
    
    ####end of growing season  -- 5 day mean > 10 C -- this might not quite get it!!--maybe need a longer smooth? other suggestions?
    annualtemp<-left_join(annualtemp,gdd %>% group_by(wateryear) %>% filter(tmeangt10 > 10) %>%
      slice_tail(n = 1) %>%
      ungroup() %>% select(wateryear,gsenddoy=yday))


###duration##############################################################################
####Do heat conditions last for shorter or longer than you remember?
        
    ###growing season length
    annualtemp<-annualtemp %>% mutate(gslength=gsenddoy-gsstartdoy)
    
    ### consecutive number of days with temp > answer to #8; 
    ##need p data to calculate
    
    ####number of days with daily tmmx greater than 90ptile for 6 days or more
    dfdaily90climos$tmxclimodiff<-ifelse(dfdaily90climos$tmmx-dfdaily90climos$dailyclimo9625tmmx<0,0,1)
    dfdaily90climos$tmxdur<-NA
    for(n in 2: nrow(dfdaily90climos)){
      dfdaily90climos$tmxdur[n]<-ifelse(dfdaily90climos$tmxclimodiff[n]==0,0,dfdaily90climos$tmxdur[n-1]+1)
    }
    annualtemp<-left_join(annualtemp, dfdaily90climos%>% 
                            filter(tmxdur>5) %>% 
                            group_by(wateryear, tmxdur) %>% 
                            summarize(n()) %>% 
                            mutate(ndays=ifelse(tmxdur==6,6,1)) %>% 
                            group_by(wateryear) %>% 
                            summarize(durtmaxgt90th=sum(ndays)))
    
    annualtemp$durtmaxgt90th<-ifelse(is.na(annualtemp$durtmaxgt90th),0,annualtemp$durtmaxgt90th)

####################################################################################################
####################DROUGHT##########################################################################
#############################################################################################################


###magnitude############################################################################
###Based on your experience in this community since you have lived or worked here, has it gotten wetter or drier overall?
###Since you have lived or worked here, has it become wetter or drier during the winter months?
    
    ###Annual water year precipitation, water deficit, eto
    if(dtasources[j]=="nclim"){
      annualdrought<-df %>% filter(wateryear %in% wateryearselect) %>% 
        group_by(wateryear) %>% 
        summarize(annualprecip=sum(pr,na.rm=TRUE))
    }
    if(dtasources[j]=="gridmet"){
    annualdrought<-df %>% filter(wateryear %in% wateryearselect) %>% 
      group_by(wateryear) %>% 
      summarize(annualprecip=sum(pr,na.rm=TRUE),
                annualeto=sum(eto, na.rm=TRUE),
                annualpwd=annualprecip-annualeto)}
    
    if(dtasources[j]=="era5land"){
      annualdrought<-df %>% filter(wateryear %in% wateryearselect) %>% 
        group_by(wateryear) %>% 
        summarize(annualprecip=sum(pr,na.rm=TRUE),
                  annualeto=sum(eto, na.rm=TRUE),
                  annualpwd=annualprecip-annualeto,
                  annualsnowf=sum(snowfall_sum, na.rm=TRUE),
                  peakswe=max(snow_depth_water_equivalent, na.rm=TRUE))}

    ###Annual winter precipitation, Oct-March; 
    annualdrought<-left_join(annualdrought, df %>% 
                               filter(wateryear %in% wateryearselect, season2=="cool") %>% 
                               group_by(wateryear) %>% 
                               summarize(coolseasonprecip=sum(pr,na.rm=TRUE)))
    
###frequency#############################################################################
####If you think it has gotten wetter or drier during the winter, how frequently is it drier?

    # ### Annual number of dry days (precip < 2 mm) 
    # annualdrought<-left_join(annualdrought, df %>% 
    #                            filter(wateryear %in% wateryearselect,pr<2) %>% 
    #                            group_by(wateryear) %>% 
    #                            summarize(drydays2mm=n()))
    
    ### Winter (Oct-March) number of dry days (precip < 2 mm)
    annualdrought<-left_join(annualdrought, df %>% 
                               filter(wateryear %in% wateryearselect, season2=="cool",pr<2) %>% 
                               group_by(wateryear) %>% 
                               summarize(coolseasondrydays2mm=n()))
    
    if(dtasources[j]=="era5land" & areas[i] %!in% c("SC\\Allendale_Box_","SC\\EPSCOR_AOI_Allendale_")){
    ### Days with snowfall
    annualdrought<-left_join(annualdrought, df %>% 
                               filter(wateryear %in% wateryearselect,snowfall_sum>=2) %>% 
                               group_by(wateryear) %>% 
                               summarize(snowdays2mm=n()))
    annualdrought$snowdays2mm<-ifelse(is.na(annualdrought$snowdays2mm),0,annualdrought$snowdays2mm)}
    

###timing#################################################################################
###Has there been a change in when snowfall/rainfall starts during the year?

    ## https://climatology.tamu.edu/research/Rainy-and-Dry-Season-RADS.html-- another method....
    # https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL090350 -- using this but changed 20 mm criterion to 10 mm b/c not met in many years in NV. Also using 2 mm for dry days instead of 1 mm
    
    ####first Julian wet day (>1 mm) starting from October 1  
    ####of 2 consecutive days receiving at least 20 mm not followed
    ####by a long dry spell (at least 10 days) receiving less than 5 mm in the following 20 days
    annualdrought<-left_join(annualdrought, df %>% 
                               mutate(pr20day=rollsum(pr,20, align="left", na.pad=TRUE),
                                      dryday10=rollsum(drydays,10, align="left", na.pad=TRUE),
                                      pr2day=rollsum(pr,2, align="right", na.pad=TRUE)) %>% 
                               filter(pr2day>=10, pr20day>=5, dryday10<10, drydays!=1) %>% 
                               group_by(wateryear) %>% 
                               slice_head(n = 1)%>%
                               ungroup() %>% 
                               select(wateryear,preciponsetdoy=yday))
    
    ### precip center of mass
    annualdrought<-left_join(annualdrought,df %>% 
                               filter(wateryear %in% wateryearselect) %>% 
                               mutate(wpr=pr*yday) %>% 
                               group_by(wateryear) %>% 
                               summarize(precipCOMdoy=sum(wpr)/sum(pr)))
    
    #######snow vars for era5land######
    if(dtasources[j]=="era5land"& areas[i] %!in% c("SC\\Allendale_Box_","SC\\EPSCOR_AOI_Allendale_")){
    ### snowonset
    
    annualdrought<-left_join(annualdrought, df %>% 
      mutate(sf20day=rollsum(snowfall_sum,20, align="left", na.pad=TRUE),
             dryday10=rollsum(drydays,10, align="left", na.pad=TRUE),
             sf2day=rollsum(snowfall_sum,2, align="right", na.pad=TRUE)) %>% 
      filter(sf2day>=10, sf20day>=5, dryday10<10, drydays!=1) %>% 
      group_by(wateryear) %>% 
      slice_head(n = 1)%>%
      ungroup() %>% 
      select(wateryear,snowonsetdoy=yday))

    ### swe center of mass
    annualdrought<-left_join(annualdrought,df %>% 
      filter(wateryear %in% wateryearselect) %>% 
      mutate(wswe=snow_depth_water_equivalent*yday) %>% 
      group_by(wateryear) %>% 
      summarize(sweCOMdoy=sum(wswe)/sum(snow_depth_water_equivalent)))

    
    #####Has there been a change in how fast the snowpack melts during the year?
    ### Day of Snow Disappearance
    ###snowmeltcomplete - last day of rolling mean of 10-day snowmelt> 0.0001 mm -- note that using snowmelt vs. Swe is about the same result
    ###sensitive to choice of threshold....not sure how best to justify....can't use 0 because values always > 0
    annualdrought<-left_join(annualdrought,df10day %>%
                               filter(wateryear %in% wateryearselect,snowmelt_sum>0.0001) %>%
                               group_by(wateryear) %>%
                               slice_tail(n = 1)%>%
                               ungroup() %>%
                               select(wateryear,snowmeltcompletedoy=yday))

    #####Has there been a change in when snowmelt starts during the year?
    ### Timing of Peak SWE
    annualdrought<-left_join(annualdrought,df %>% 
                               filter(wateryear %in% wateryearselect) %>% 
                               group_by(wateryear) %>% 
                               slice_max(snow_depth_water_equivalent,n = 1, with_ties = FALSE)%>%
                               ungroup() %>%
                               select(wateryear,peakswedoy=yday))
    
    ### snowmeltstart
    ###sensitive to choice of threshold....not sure how best to justify....can't use 0 because values always > 0
    annualdrought<-left_join(annualdrought,df10day %>%
                               filter(wateryear %in% wateryearselect,snowmelt_sum>0.0001) %>%
                               group_by(wateryear) %>%
                               slice_head(n = 1)%>%
                               ungroup() %>%
                               select(wateryear,snowmeltstartdoy=yday))
    
    }
    
#######duration##############################################################################
####Do drier conditions last shorter or longer than you remember?
    df$drydaydur<-NA
    df$drydaydur[1]<-0
    for(n in 2: nrow(df)){
      df$drydaydur[n]<-ifelse(df$drydays[n]==0,0,df$drydaydur[n-1]+1)
    }

    ### Water year maximum consecutive number of dry days
    annualdrought<-left_join(annualdrought, df %>% 
                            group_by(wateryear) %>% 
                            summarize(maxdrydur= max(drydaydur, na.rm=TRUE)))
    
    ### Winter maximum consecutive number of dry days
    annualdrought<-left_join(annualdrought, df %>% 
                               filter(season2=="cool") %>% 
                               group_by(wateryear) %>% 
                               summarize(maxcoolseasondrydur= max(drydaydur, na.rm=TRUE)))
    
    ###days with SWE> 0
    ###sensitive to choice of threshold....not sure how best to justify....can't use 0 because values always > 0
    ###SWE duration 
    if(dtasources[j]=="era5land" & areas[i] %!in% c("SC\\Allendale_Box_","SC\\EPSCOR_AOI_Allendale_")){
    annualdrought<-left_join(annualdrought, df %>% 
                               filter(snow_depth_water_equivalent>0.0001) %>% 
                               group_by(wateryear) %>% 
                               summarize(sweduration= n()))}


###############################################################################################################
###calculate trends###########################################################################################
#########################################################################################################
    
    anndata<-left_join(annualtemp, annualdrought)
    write.csv(anndata, paste0(dtasources[j], "_",areas2[i],"AnnualMetrics_zonalmean.csv"),row.names=FALSE)
    
    anndatalong<-anndata %>% pivot_longer(2:ncol(anndata))
    slopes<-as.data.frame(matrix(nrow=0, ncol=4, data=NA))
    names(slopes)<-c("sen_slope","p_value","TimePeriod","datasource")
    
    if(dtasources[j]!="gridmet"){
    startyears<-c(1955, 1960, 1965, 1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)
    }else{
    startyears<-c(1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)}
    
    for(n in 1: length(startyears)){
      startyear<-startyears[n]
      slopesn<-anndatalong %>% 
        group_by(name) %>% 
        filter(wateryear>=startyear, !is.na(value)) %>% 
        summarize(sen_slope = sens.slope(value)$estimates,
                  p_value = mk.test(value)$p.value)
      slopesn$TimePeriod<-paste0(startyear, "-2025")
      slopesn$datasource<-dtasources[j]
      slopes<-rbind(slopes, slopesn)
    }
    if(j==1){slopes2<-slopes}else{slopes2<-rbind(slopes2, slopes)}
  }
  slopes2wide<-slopes2 %>% pivot_wider(names_from=datasource, values_from=c(sen_slope, p_value))
  write.csv(slopes2wide, paste0(areas2[i],"Trends_zonalmean.csv"), row.names=FALSE)
 }
