### Download updated/current sea level (sea surface height, ssh) information and plot #### 
# format for SSH info for dfa or other analysis
# calculate mean spring sea level for dfa analysis

rm(list=ls())
library(tidyverse)
library(ggplot2)

# file locations
home_file = getwd()
# where to put ssh downloads
raw_data = paste0(home_file,'/Raw_SSH_Data/')
# where to put processed data
processed_data =  paste0(home_file,'/Data/')

# choose data type meantrend or intannvar
# intannvar used in Tolimier & Haltuch 2023
data_type = 'intannvar' #'meantrend' - which tide gauge data to use.

### Scrape sea level data data #############################################################
# note: these stations are sub-setted for the DFA analysis. 
# only 16 are used

STATIONS_0 = c(
     "9440910_TokePoint", ########### WA
     "9443090_NeahBay",
     "9444090_PortAngeles",
     "9444900_PortTownsend",
     "9447130_Seattle",
     "9449424_CherryPoint",
     "9449880_FridayHarbor",
     "9431647_PortOrford", ########## OR
     "9432780_Charleston",
     "9435380_SouthBeachOR",
     "9437540_Garibaldi",
     "9439011_Hammond",
     "9439040_Astoria",
     '9410170_SanDiego', ############ CA
     '9410230_LaJolla', 
     '9410580_NewportBeach',
     '9410660_LosAngeles',
     '9410840_SantaMonica',
     '9411270_RinconIsland',
     '9411340_SantaBarbara',
     '9412110_PortSanLuis',
     '9413450_Monterey',
     '9414290_SanFrancisco',
     '9414523_RedwoodCity',
     '9414750_Alameda',
     '9415020_PointReyes',
     '9415144_PortChicago',
     '9416841_ArenaCove',
     '9418767_NorthSpit',
     '9419750_CrescentCity'
)

STATIONS = paste0(STATIONS_0,"_",data_type,".csv")

#### Download updated sea level timeseries  ##############################################
new.data=FALSE # just a trigger to turn on and off
if(new.data==TRUE){
     
     for(i in 1:length(STATIONS)){
          print(i)
          station = substring(STATIONS[i], 1,7)
          URL = paste0("https://tidesandcurrents.noaa.gov/sltrends/data/",station, "_",data_type,".csv")
          DESTFILE = paste0(raw_data, unlist(strsplit(STATIONS[i],"_"))[2],"_",station,"_",data_type,".csv")
          download.file(url=URL, destfile=DESTFILE)
     } # end i
} # end if

############### RELOAD AND PROCESS SSH DATA into Wide Data frame for DFA analysis #################

# get list all files
ssh_files0  = dir(raw_data)
# get list of just correct data_type (intannvar or meantrend)
# inannvar used in Tolimieri & Haltuch 2023
ssh_files = ssh_files0[grep(data_type,ssh_files0)]

# combine all sea level files into one file
for(k in 1:length(ssh_files)){
  site = ssh_files[k]
  print(paste("k = ", k, "site = ",site))
  fname = paste0(raw_data,ssh_files[k])
  # load csv
  d1 = data.frame(read.csv(fname,header=FALSE, skip=1))
  # name columns
  if(dim(d1)[2]==8){
     colnames(d1)[c(1,2,3)]=c('year','month','ssh')}else{
     colnames(d1)[c(1,2,3)]=c('year','month','ssh')}
  d2 = d1[,c('year', 'month','ssh')]
  d2$site = unlist(strsplit(site,"_"))[1]
  if(k==1){SSH = d2}else{SSH = rbind(SSH,d2)}
} # end k

SSH$id = paste(SSH$year,SSH$month,sep="_")
# save out raw data
write.csv(SSH, paste0(processed_data, "Sea_level_data_combined_monthly.csv"), row.names = FALSE)
saveRDS(SSH, file = paste0(processed_data, "Sea_level_data_combined_monthly.rds"))

#### get spring average ##########################

ssh_spring_long = SSH %>% 
  filter(month %in% c(4:6)) %>%
  group_by(year,site) %>%
  summarise(st_dev = sd(ssh, na.rm = TRUE),
            ssh = mean(ssh, na.rm=TRUE)) %>%
  select(year, site, ssh, st_dev)

# save out spring average in long format (useful for ggplot figures)
write.csv(ssh_spring_long, paste0(processed_data, "Sea_level_data_spring_long.csv"), row.names = FALSE)
saveRDS(ssh_spring_long, file = paste0(processed_data, "Sea_level_data_spring_long.rds"))

# reshape to wide format for MARSS/DFA ###########
ssh_spring_wide = ssh_spring_long %>%
  select(year, site, ssh) %>%
   pivot_wider( ., values_from = ssh, names_from = year)

write.csv(ssh_spring_wide, paste0(processed_data, "Sea_level_data_spring_wide.csv"), row.names = FALSE)
saveRDS(ssh_spring_wide, file = paste0(processed_data, "Sea_level_data_spring_wide.rds"))


#### FIGURES #################################################
#### set up site selection and plot order ####################

location = levels(factor(ssh_spring_wide$site))

# set plot order for ssh graph
site_order = data.frame(location)
site_order$order = NA
site_order$order[site_order$location=="NeahBay"] <-16
site_order$order[site_order$location=="TokePoint"] <-15
site_order$order[site_order$location=="Astoria"] <- 14
site_order$order[site_order$location=="SouthBeachOR"] <- 13
site_order$order[site_order$location=="Charleston"] <- 12
site_order$order[site_order$location=="PortOrford"] <- 11
site_order$order[site_order$location=="CrescentCity"] <- 10
site_order$order[site_order$location=="NorthSpit"] <- 9
site_order$order[site_order$location=="PointReyes"] <- 8
site_order$order[site_order$location=="SanFrancisco"] <- 7
site_order$order[site_order$location=="Monterey"] <- 6  
site_order$order[site_order$location=="PortSanLuis"] <- 5
site_order$order[site_order$location=="SantaMonica"] <- 4
site_order$order[site_order$location=="LosAngeles"] <- 3
site_order$order[site_order$location=="LaJolla"] <- 2
site_order$order[site_order$location=="SanDiego"] <- 1
site_order <- site_order[order(site_order$order, decreasing =TRUE),]
# drop sites not on list above
site_order <- na.omit(site_order)
# add spaces between capitals
site_order$site = gsub("([a-z])([A-Z])","\\1 \\2",site_order$location)
site_order$site = stringr::str_remove(site_order$site, " OR")
site_order

write.csv(site_order, 'Site_Order.csv', row.names = FALSE)

######## plot sea level ########################################################
# get spaces in names
ssh_spring_long$location = site_order$site[ match(ssh_spring_long$site, site_order$location) ]
# set plot order
ssh_spring_long$location = factor(ssh_spring_long$location, levels = site_order$site)
# 
# get rid of any NA sites
df = ssh_spring_long %>% filter(location %in% site_order$site)
min_yr = 1925
max_yr = max(df$year, na.rm=TRUE)

graphics.off()
png(paste0("Figure_Spring_sea_level_",data_type,"_",min_yr,"-",max_yr,".png"), units='in',width=6, height=8, res=300)

ssh_plot <- ggplot(df, aes(x=year, y=ssh)) +
  geom_line() + 
  xlab("Year")+xlim(min_yr,max_yr) +
  ylab("Interannual variation in sea level (m)") +
  facet_wrap( facets = df$location,  nrow = 8, ncol=2) + 
  theme_bw()

ssh_plot

dev.off()


































