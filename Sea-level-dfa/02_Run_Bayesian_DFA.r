# Run bayesian dynamic factor analysis to obtain northern sea level index for sablefish assessment ####

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

##### load data files ###############################################################
# for plotting etc
site_order = data.frame(read.csv('Site_Order.csv',header=TRUE))
# the ssh data
ssh_x = readRDS( paste0(data_loc, 'Sea_level_data_spring_wide.rds'))

##### prep data for DF #########################################################
# set years to include
years1 = colnames(ssh_x) 
years1 = as.numeric(years1[years1 != 'site'])
# reset min year here for shorter dfa time series ####
dfa_years = as.character(1925:max(years1))
# subset years
ssh = as.matrix(ssh_x[,dfa_years])
rownames(ssh) <- ssh_x$site
# subset to 16 sites
ssh = ssh[site_order$location,]
# get range for labeling files
yrmin = min(dfa_years)
yrmax = max(dfa_years)

# archive analyzed data
write.csv(ssh, paste0(data_loc,"Sea_level_for_DFA_",yrmin,'-',yrmax,".csv"), row.names = TRUE)
saveRDS(ssh, file = paste0(data_loc,"Sea_level_for_DFA_",yrmin,'-',yrmax,".rds"))

###################################################################################################
#### Run Bayes DFA ################################################################################
###################################################################################################

## Model comparison ####
# fit multiple models and compare to get best fit model
# note model selection in Tolimieri & Haltuch was based of
# of original MARSS (non Bayesian) DFA analysis in the 
# 2019 sablefish stock assessment

# very, very long run time
run_model_selection = FALSE

if(run_model_selection == TRUE){ # begin if
  
  R_un = 1:length(rownames(ssh))
  
  model_comparison = fit_dfa_trends(
    y = ssh,
    kmin = 1, # min trends 
    kmax = 5, # max trends
    compare_normal = TRUE,
    variance = c('equal','unequal')
  )
  
  model_comparison$summary
  model_comparison$best_model

saveRDS( model_comparison, file = paste0( results_loc, "BayesDFA_Model_Comparison_", yrmin,'-',yrmax,'.rds'))  
    
} # end if


#### Run just one model ###########################################################################
# run the Tolimieri & Haltuch 2023 Bayesian DFA

run_dfa = FALSE # trigger to run or not run model 

if(run_dfa == TRUE){ # begin if
  # R_eq = rep(1,length(rownames(ssh)))
  R_un = 1:length(rownames(ssh))
  bfit_5 <- fit_dfa( y = ssh, 
                 num_trends = 5,
                 varIndx = R_un,
                 chains = 3,
                 iter = 3000,
                 scale = "zscore",
                 data_shape = 'wide'
                 # warmup = warmups,
                 # thin=1,
                 # estimate_nu = TRUE,
                 # estimate_trend_ar = TRUE,
                 #control = list(adapt_delta = 0.99,
                 #max_treedepth = 30
  )
  
  saveRDS(bfit_5 , paste0(results_loc,"Bayesian-DFA-SSH-DF5-",yrmin,'-',yrmax,"-fit.rds"))
  # rotate loadings and trends for saving out
  rotated_bfit_5 <- rotate_trends(bfit_5)
  saveRDS(rotated_bfit_5, paste0(results_loc,"Bayesian-DFA-SSH-DF5-",yrmin,'-',yrmax,"-rotated.rds"))  
  
}# end if

####### Bayesian figures ###############################
# reload model instead of re-running bayes dfa; select model
model_name = paste0("Bayesian-DFA-SSH-DF5-",yrmin,'-',yrmax,"-fit.rds")
bmod = readRDS(paste0(results_loc,model_name))

# rotate loadings and trends
df_rotated <- rotate_trends(bmod)

# some prelim stuff for prettier plots #################
# reorder to plot north-top south-bottm

#get and add regions to various sections
regions = rownames(bmod$orig_data)
colnames(df_rotated$Z_rot) <- regions
rownames(df_rotated$Z_rot_mean) <- regions
rownames(df_rotated$Z_rot_median) <- regions

# reorder for nicer plots
region_order = rev(rownames(bmod$orig_data))
df_rotated$Z_rot = df_rotated$Z_rot[,region_order,]
df_rotated$Z_rot_mean = df_rotated$Z_rot_mean[region_order,]
df_rotated$Z_rot_median = df_rotated$Z_rot_median[region_order,]
bmod$orig_data = bmod$orig_data[region_order,]

# for file naming
yrmin = min(colnames(bmod$orig_data))
yrmax = max(colnames(bmod$orig_data))

## plot fitted model ####
graphics.off()
png( paste0("Figure_Fitted_Model_",yrmin,'-',yrmax,'.png') )
plot_fitted(bmod, names = rownames(bmod$orig_data))+ theme_bw() 
dev.off()

## plot rotated loadings ####
graphics.off()
png( paste0("Figure_Rotated_Loadings_",yrmin,'-',yrmax,'.png') )
plot_loadings(df_rotated)
dev.off()

## plot rotated trends ####
graphics.off()
png( paste0("Figure_Rotated_Trends_",yrmin,'-',yrmax,'.png') )
plot_trends(df_rotated)
dev.off()

## Output rotated trends for stock assessment ####

trends = t(df_rotated$trends_mean)
colnames(trends) <- c(paste0('DF',1:5))
uCI = t(df_rotated$trends_upper)
colnames(uCI) <- c(paste0('uCI',1:5))
lCI = t(df_rotated$trends_lower)
colnames(lCI) <- c(paste0('lCI',1:5))
df.out <- data.frame(cbind(trends,uCI,lCI))
df.out$year = yrmin:yrmax
df_out <-  df.out %>% select(year,
                             DF1, DF2, DF3, DF4, DF5, 
                             lCI1, lCI2, lCI3, lCI4, lCI5,
                             uCI1, uCI2, uCI3, uCI4, uCI5
                             )
write.csv(df_out, paste0('Data_Sea_Level_Indices_',yrmin,"_",yrmax,'.csv'), row.names = FALSE)
saveRDS(df_out, file = paste0('Data_Sea_Level_Indices_',yrmin,"_",yrmax,'.rds'))

#########################################################################################
######## Run a DFA with one dynamic factor ##############################################
# run a one (1) DF version for complarison

run_one = FALSE # trigger; TRUE or FALSE

if(run_one == TRUE){ # begin if

  R_un = 1:length(rownames(ssh))
  bfit_1 <- fit_dfa( y = ssh, 
                   num_trends = 1,
                   varIndx = R_un,
                   chains = 4,
                   iter = 5000,
                   scale = "zscore",
                   data_shape = 'wide'
                   # warmup = warmups,
                   # thin=1,
                   # estimate_nu = TRUE,
                   # estimate_trend_ar = TRUE,
                   #control = list(adapt_delta = 0.99,
                   #max_treedepth = 30)
  )

  saveRDS(bfit_1 , paste0(results_loc,"Bayesian-DFA-SSH-DF1-",yrmin,'-',yrmax,"-fit.rds"))
  # rotate loadings and trends for saving out
  rotated_bfit_1 <- rotate_trends(bfit_1)
  saveRDS(rotated_bfit_1, paste0(results_loc,"Bayesian-DFA-SSH-DF1-",yrmin,'-',yrmax,"-rotated.rds"))  
}

m1 = readRDS(paste0(results_loc,"Bayesian-DFA-SSH-DF1-",yrmin,'-',yrmax,"-fit.rds"))
r1 = rotate_trends(m1)
plot_fitted(m1)
plot_loadings(r1)
plot_trends(r1)

# Northern sea level index from 1-trend and 5-trend model
# just for comparison

plot(r1$trends_mean, df_rotated$trends_mean[1,],
     xlab="1-trend DFA", ylab='5-trend DFA')
x = as.vector(r1$trends_mean)
y = as.vector(df_rotated$trends_mean[1,])
cor1 = cor.test(x,y)
cor1


####################################################################################################
###  MARSS DFA for comparison ######################################################################
# just for fun 
# settings below are mostly manual


run_marss_dfa = FALSE # trigger to run or not

if(run_marss_dfa == TRUE){ # begin if

library(MARSS)

iterations = 2000
cntl.list = list(maxit=2000)
model.list = list( m = 5, R = "diagonal and unequal" )

marss_fit = MARSS(ssh, model=model.list, method='BFGS',
             demean=TRUE, z.score=TRUE, form="dfa", control=cntl.list)
# save out model
saveRDS(marss_fit, paste0(results_loc,"MARSS_dfa.rds"))

# varimax rotation

Z.est = coef(marss_fit, type="matrix")$Z
H.inv <- 1
if (ncol(Z.est) > 1){H.inv <- varimax(coef(marss_fit, type = "matrix")$Z)$rotmat}
# rotate factor loadings
Z.rot <- Z.est %*% H.inv
# rotate trends
trends.rot <- solve(H.inv) %*% marss_fit$states

# Add CIs to marssMLE object
marss_fit <- MARSSparamCIs(marss_fit)

# Use coef() to get the upper and lower CIs -- factor loadings
Z.low <- coef(marss_fit, type = "Z", what = "par.lowCI")
Z.up <- coef(marss_fit, type = "Z", what = "par.upCI")
Z.rot.up <- Z.up %*% H.inv
Z.rot.low <- Z.low %*% H.inv
df <- data.frame(
  est = as.vector(Z.rot),
  conf.up = as.vector(Z.rot.up),
  conf.low = as.vector(Z.rot.low)
)

# plot rotated marss loadings
df$factor = c(rep('DF1',16),rep('DF2',16),rep('DF3',16),rep('DF4',16),rep('DF5',16))
locs = rownames(Z.rot)
df$site = rep(locs,5)
# set order for ploting

df$site = factor(df$site, levels = rev(site_order$location))

ggplot(df, aes(x=est, y = site)) +
  geom_point() + 
  geom_errorbarh( aes(xmin=conf.low, xmax=conf.up)) +
  annotate("segment", 
           x = 0, xend = 0, 
           y = 1, yend = 17,
           colour = "blue") + 
  facet_wrap(facets = 'factor', ncol = 5)+
  xlab('Loading') + ylab("Site") + 
  theme_bw() 


# plot marss trends ####
marss_trends = data.frame(t(trends.rot)) %>% 
  mutate(year = yrmin:yrmax) %>% 
  rename(DF1 = X1, DF2 = X2, DF3 = X3, DF4 = X4, DF5 = X5) %>%
  pivot_longer(cols=DF1:DF5, values_to = 'trend',names_to = 'factor')


head(marss_trends)

ggplot(marss_trends, aes(x = year, y = trend)) +
  geom_line() + 
  facet_grid(rows='factor')+
  xlab("Year") + ylab("Trend") +
  theme_bw()

# compare MARSS & bayes dfa ####

plot(df_rotated$trends_mean[1,], trends.rot[1,],
     xlab = "Bayes - DF1", ylab='MARSS - DF1')

x = as.vector(df_rotated$trends_mean[1,])
y = as.vector(trends.rot[1,])
cor(x,y)
} # end if



