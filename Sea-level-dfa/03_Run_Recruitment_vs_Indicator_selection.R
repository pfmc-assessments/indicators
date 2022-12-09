

# load packages ####
library(bayesdfa)
library(tidyverse)
library(ggplot2)
library(MuMIn)

# set file locations
home_file = getwd()
data_loc = paste0(home_file, '/Data/')
results_loc = paste0(home_file, '/Results/')

# Bring in fish data and combine with DFA results #################################
# fish
fish = data.frame(read.csv( paste0(data_loc, 'Assessment_Data_w_StDev_2022_12.csv')))
fish$year = substring(fish$Label, nchar(fish$Label)-3, nchar(fish$Label))
x = grep("RecrDev", fish$Label)
RecDev = fish[x,] %>%
  rename(RecDev = Value) %>% 
  select(year,RecDev)
# DFA
dfa = readRDS( "Data_Sea_Level_Indices_1925_2022.rds")
# combine
dfa$year = as.numeric(dfa$year)
RecDev$year = as.numeric(RecDev$year)
dfile = full_join(RecDev,dfa)

# set analysis years #############################################################
fish_years = 1975:2020
dfile = dfile %>% 
  filter(year %in% fish_years)
write.csv( dfile, paste0(data_loc,'Data_combined_4dredge_',min(fish_years),'-',max(fish_years),'.csv'), row.names = FALSE)

# Set up model equation ##########################################################
# get df names
cn1 = colnames(dfile %>% select(contains("DF")))
# cn1 = c('DF1','DF2','DF3','DF4','DF5') # 
cn2 = paste0("I(",cn1,"^2)")
cn = c(cn1,cn2)
# make equation
for(i in 1:length(cn)){
  x1 = cn[i]
  if(i==length(cn)){x2 = x1}else{x2 = paste(x1,"+")}
  if(i==1){eq = x2}else{eq = paste(eq,x2)}
}
# quick look
eq

## Run DREDGE #####################################################################
# make a formula 
form1 = as.formula(paste0('RecDev ~ ', eq))

# fit full model; step 1 for 'dredge'
fit = lm(form1, data=dfile, na.action = na.fail)

#### run dredge ####

mtable = dredge(fit, m.lim = c(NA, 5),
                subset= # require both linear and quadratic
                  dc(DF1, I(DF1^2)) &&
                  dc(DF2, I(DF2^2)) &&
                  dc(DF3, I(DF3^2)) &&
                  dc(DF4, I(DF4^2)) &&
                  dc(DF5, I(DF5^2)),
                extra = list(R2 = function(x)
                  summary(x)$r.squared[[1]],
                  F = function(y)
                    summary(y)$fstatistic[[1]]),
                trace=2 )

mtable2 = subset(mtable, delta<2)
mtable2
mtable4 = subset(mtable, delta<4)
mtable4
mtable =data.frame(mtable)
write.csv(mtable, paste0(results_loc, "AIC-table-dredge-", min(fish_years),'-',max(fish_years), '.csv') )

# fit just DF1 - northern sea level index

fit_df1 = lm(RecDev ~ DF1, data=dfile)
summary(fit_df1)
anova(fit_df1)












