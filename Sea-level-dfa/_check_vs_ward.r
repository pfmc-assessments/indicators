#### check dfs vs Eric Ward model ################

# load packages for 'bayesdfa'
library(bayesdfa)
library(rstan)
library(rstudioapi)
library(tidyverse)
library(ggplot2)

# options
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

# set file locations
home_file = getwd()
data_loc = paste0(home_file, '/Data/')
results_loc = paste0(home_file, '/Results/')


# eric ward values # shorter sl level (through 2020)

EW_loc = "C:/Users/nick.tolimieri/Documents/GitHub/Sablefish-SSH-2022/Bayes-SHH-Output/"
df_ew <- data.frame(read.csv(paste0(EW_loc,"Data_for_dredge.csv"), header = TRUE))

# current values

df_new = data.frame(read.csv( 
  paste0(data_loc,'Data_combined_4dredge_',min(fish_years),'-',max(fish_years),'.csv'),header=TRUE))

df_ew = df_ew %>% filter(year %in% df_new$year)  

plot(df_ew$DF1, df_new$DF1, 
     xlab = "Eric Ward bayesdfa",
     ylab = "Nick T bayesdfa")
# small changes due to longer time-frame for the NT dfa (additional years)


################
dfile = df_ew
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
mtable  = data.frame(mtable)

