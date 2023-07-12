### Produces qqplots of model and obs data ##
library(ncdf4)
library(stringr)
library(ggplot2)
library(plotly)
library(stringr)


############################# SET UP ################################
MODELS <- c(
  'MRI-CGCM3',
  'MPI-ESM-LR',
  'ACCESS1-3',
  'IPSL-CM5A-MR'
  # 'PPE_multimodel_mean',
  # 'CMIP_multimodel_mean',
  # 'full_multimodel_mean'
)

DOMAIN <- c(
  #'global'
  'regional',
  'local'
)

DSI_NUMBER <- c(
  3,
  6,
  12
)

COUNTRIES <- list(
  'uk',
  'Scotland',
  'England',
  'Wales',
  'Northern_Ireland'
)

get_unmasked_nonzero_data <- function(cube) {
  data <- ncvar_get(cube)
  data <- data[data < 1000]
  data <- data[data != 0]
  return(data)
}

load_model_cube <- function(country, domain, model, index_number) {
  firstmonth <- str_pad(index_number - 1, width = 2, pad = "0")
  model_path <- paste0('/data/users/mgrant/ukcp/droughts/data/', country, '/',
                       domain,'/', model, '/TS1/monthly/dsi-', index_number, 
                       '/', model, '_1981', firstmonth, '-200011_dsi-', 
                       index_number, '.nc')
  
  model_cube <- nc_open(model_path)
  
  return(model_cube)
}

load_obs_cube <- function(country, index_number, resolution) {
  firstmonth <- str_pad(index_number - 1, width = 2, pad = "0")
  obs_path <- paste0('/data/users/mgrant/ukcp/droughts/data/', country, '/obs/', 
                     resolution, '/dsi-', index_number, 
                     '/1980-2000/obs_', resolution, '_1981', firstmonth, 
                     '-200011_dsi-', index_number, '.nc')
  
  obs_cube <- nc_open(obs_path)
  return(obs_cube)
}

############################## QQ PLOTS ###################################
# Function which creates a qqplot
create_qqplot <- function(domain, model, dsi_number, scale=TRUE) {
  # Capitalise the first letter of the domain for titles
  Domain <- str_to_title(domain)
  
  # Load in the observation data
  obs_cube <- load_obs_cube(dsi_number, '1km')
  obs_data <- get_unmasked_nonzero_data(obs_cube)
  
  # Load in the model data
  model_cube <- load_model_cube(domain, model, dsi_number)
  model_data <- get_unmasked_nonzero_data(model_cube)
  
  if (scale){
    model_data <- scale(model_data)
    obs_data <- scale(obs_data)
  }
  
  # Find the max value between the two so that the axes are even
  max_value <- max(obs_data, model_data)
  qqplot(obs_data, model_data, xlim=c(0, max_value), ylim=c(0, max_value), 
         xlab='Observations', ylab=paste0(Domain, ' Model'),
         main=paste0(model, ' ', Domain, ' vs Observations DSI-', dsi_number))
  abline(0, 1, col=2)
}


# # This tests the function for any changes that need made
# filepath <- c("/scratch/mgrant/UKCP/droughts/plots/distribution_comparisons/qqplots/model_vs_obs/global/MRI-CGCM3/dsi-3/")
# 
# if (!dir.exists(filepath)) {
#   dir.create(filepath, recursive=TRUE)
# }
# 
# png(paste0(filepath, "qqplot_global_MPI-CGCM3-vs-obs_dsi-3.png"))
# 
# create_qqplot('global', 'MRI-CGCM3', 3)
# dev.off()

# Create all the qqplots
# for (number in DSI_NUMBER){
#   obs_cube <- load_obs_cube(number, '1km')
#   obs_data <- get_unmasked_nonzero_data(obs_cube)
#   for (dom in DOMAIN){
#     for (model in MODELS){
#       model_cube <- load_model_cube(dom, model, number)
#       model_data <- get_unmasked_nonzero_data(model_cube)
# 
#       filepath <- paste0("/scratch/mgrant/UKCP/droughts/plots/",
#                          "distribution_comparisons/qqplots/model_vs_obs/",
#                          dom, "/", model, "/dsi-", number, "/scaled_anomalies/")
# 
#       if (!dir.exists(filepath)) {
#         dir.create(filepath, recursive=TRUE)
#       }
# 
#       filename <- paste0('qqplot_', dom, '_', model,
#                          '-vs-obs_dsi-', number, '.png')
# 
#       png(paste0(filepath, filename), width = 800, height = 600)
# 
#       create_qqplot(dom, model, number, scale=TRUE)
#       dev.off()
#       
#       model_cube <- load_model_cube(dom, model, number)
#       model_data <- get_unmasked_nonzero_data(model_cube)
# 
#       filepath <- paste0("/scratch/mgrant/UKCP/droughts/plots/",
#                          "distribution_comparisons/qqplots/model_vs_obs/",
#                          dom, "/", model, "/dsi-", number,
#                          "/absolute_distribution/")
# 
#       if (!dir.exists(filepath)) {
#         dir.create(filepath, recursive=TRUE)
#       }
# 
#       filename <- paste0('qqplot_', dom, '_', model,
#                          '-vs-obs_dsi-', number, '.png')
# 
#       png(paste0(filepath, filename), width = 800, height = 600)
# 
#       create_qqplot(dom, model, number, scale=FALSE)
#       dev.off()
#     }
#   }
# }

###################### QQ PLOTS AND HISTOGRAMS #############################

# Function which creates a figure of local, global, and obs histograms; as well
# as the qqplots comparing local vs obs, global vs obs, and local vs global
create_hist_and_qqplots <- function(model, dsi_number, scale = TRUE) {
  # Load in the observation data
  obs_cube <- load_obs_cube(dsi_number, '1km')
  obs_data <- get_unmasked_nonzero_data(obs_cube)
  
  # Load in the global model data
  global_cube <- load_model_cube('global', model, dsi_number)
  global_data <- get_unmasked_nonzero_data(global_cube)
  
  # Load in the regional model data
  regional_cube <- load_model_cube('regional', model, dsi_number)
  regional_data <- get_unmasked_nonzero_data(regional_cube)
  
  # Load in the local model data
  local_cube <- load_model_cube('local', model, dsi_number)
  local_data <- get_unmasked_nonzero_data(local_cube)
  
  if (scale){
    obs_data <- scale(obs_data)
    local_data <- scale(local_data)
    regional_data <- scale(regional_data)
    global_data <- scale(global_data)
  }
  
  # Find the max value between the two so that the axes are even
  max_total <- max(obs_data, global_data, regional_data, local_data)
  min_total <- min(obs_data, global_data, regional_data, local_data)
  
  max_model <- max(local_data, regional_data, global_data)
  min_model <- min(local_data, regional_data, global_data)
  
  max_obs_local <- max(obs_data, local_data)
  min_obs_local <- min(obs_data, local_data)
  
  max_obs_regional <- max(obs_data, regional_data)
  min_obs_regional <- min(obs_data, regional_data)

  max_obs_global <- max(obs_data, global_data)
  min_obs_global <- min(obs_data, global_data)
  # Create the plots
  par(mfrow=c(2, 4))
  
  # Increase the outer margin size to make room for the main title
  par(oma = c(0.5, 0.5, 5, 0))  # Adjust the last value (2) as needed
  
  hist(obs_data, col='coral1', main='Observed', 
       xlim=c(min_total, max_total), 
       xlab=paste0('Observed DSI-', dsi_number),
       cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
  
  hist(global_data, col='coral1', main='Global', 
       xlim=c(min_total, max_total),
       xlab=paste0('Global DSI-', dsi_number),
       cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
  
  hist(regional_data, col='coral1', main='Regional', 
       xlim=c(min_total, max_total),
       xlab=paste0('Regional DSI-', dsi_number),
       cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
  
  hist(local_data, col='coral1', main='Local', xlim=c(min_total, max_total),
       xlab=paste0('Local DSI-', dsi_number),
       cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
  
  qqplot(obs_data, global_data, xlim=c(min_obs_global, max_obs_global), 
         ylim=c(min_obs_global, max_obs_global), xlab='Observations', 
         ylab='Global Model',
         main='Global vs Observations',
         cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
  abline(0, 1, col=2)
  
  qqplot(obs_data, regional_data, xlim=c(min_obs_regional, max_obs_regional), 
         ylim=c(min_obs_regional, max_obs_regional), xlab='Observations', 
         ylab='Regional Model', 
         main='Regional vs Observations',
         cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
  abline(0, 1, col=2)
  
  qqplot(obs_data, local_data, xlim=c(min_obs_local, max_obs_local), 
         ylim=c(min_obs_local, max_obs_local), xlab='Observations', 
         ylab='Local Model',
         main='Local vs Observations',
         cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
  abline(0, 1, col=2)
  
  # Add a main title
  if (scale){
    main_title <- paste0("1980-2000 UK DSI-", dsi_number, " Scaled Anomalies ",
                         "Distributions for ",
                         "Observations and ", model, " on All Domains")
  } else {
    main_title <- paste0("1980-2000 UK DSI-", dsi_number, " Distributions for ",
                         "Observations and ", model, " on All Domains")
  }
  par(cex.main = 3)
  title(main = main_title, outer = TRUE)
}


## Tests the above function
# filepath <- '/scratch/mgrant/'
# 
# if (!dir.exists(filepath)) {
#  dir.create(filepath, recursive=TRUE)
# }
# 
# png(paste0(filepath, 'qqplot_global_MRI-CGCM3-vs-obs_dsi-3.png'),
#           width = 1200, height = 900)
# 
# create_hist_and_qqplots('MPI-ESM-LR', 6, scale=TRUE)
# 
# dev.off()


######################### CREATE HIST AND QQ-PLOTS ############################
# # Creates the histograms and qqplots for all the models and dsi numbers
# for (number in DSI_NUMBER) {
#   for (model in MODELS) {
#     # print('---------------------------------------')
#     # print('Starting plot for:')
#     # cat('Model: ', model, "\n")
#     # cat('Index: DSI-', number, "\n")
#     # filepath <- paste0('/scratch/mgrant/UKCP/droughts/plots/',
#     #                    'distribution_comparisons/qq_and_hist/model_vs_obs/',
#     #                    'domain_comparison/', model, '/dsi-', number,
#     #                    '/scaled_anomalies/')
#     # 
#     # if (!dir.exists(filepath)) {
#     #   dir.create(filepath, recursive=TRUE)
#     # }
#     # 
#     # filename <- paste0('distribution_comparison_',
#     #                    model, '-vs-obs_dsi-', number, '.png')
#     # 
#     # png(paste0(filepath, filename), width = 1200, height = 900)
#     # 
#     # create_hist_and_qqplots(model, number, scale=TRUE)
#     # dev.off()
#     # print('Plot completed')
#     # # For the distributions
#     filepath <- paste0('/scratch/mgrant/UKCP/droughts/plots/',
#                        'distribution_comparisons/qq_and_hist/model_vs_obs/',
#                        'domain_comparison/', model, '/dsi-', number,
#                        '/absolue_distribution/')
# 
#     if (!dir.exists(filepath)) {
#       dir.create(filepath, recursive=TRUE)
#     }
# 
#     filename <- paste0('distribution_comparison_',
#                        model, '-vs-obs_dsi-', number, '.png')
# 
#     png(paste0(filepath, filename), width = 1200, height = 900)
# 
#     create_hist_and_qqplots(model, number, scale=FALSE)
#     dev.off()
#     
#     # For the scaled anomaly distributions
# 
#     #MRI-CGCM3 dsi-6 problem
#   }
# }

# 
# global_cube <- load_model_cube('global', 'MRI-CGCM3', 3)
# global_data <- get_unmasked_nonzero_data(global_cube)
# hist(scale(global_data), col='coral1', main='Global',
#      xlab=paste0('Global DSI-', 3),
#      cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
# 
# obs_cube <- load_obs_cube(3, '60km')
# obs_data <- get_unmasked_nonzero_data(obs_cube)
# qqplot(scale(obs_data), scale(global_data), xlab='Observations', 
#        ylab='Global Model', main='Global vs Observations',
#        cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
# abline(0, 1, col=2)
# 
# qqplot(obs_data, global_data, xlab='Observations', ylim = c(0, 70), ylab='Global Model',
#        main='Global vs Observations',
#        cex.main=1.8, cex.lab=1.5, cex.axis=1.5)
# abline(lm(sort(global_data) ~ sort(obs_data)), col=2)
