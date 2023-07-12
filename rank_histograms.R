library(stringr)
library(SpecsVerification)
library(parallel)
library(iterators)
library(foreach)
# install.packages("doParallel")
library(doParallel)

source("/net/home/h02/mgrant/ukcp/ukcp_2.0/droughts/qqplots.R")

# CMIP_MODELS <- c(
#     'MRI-CGCM3',
#     'MPI-ESM-LR',
#     'ACCESS1-3',
#     'IPSL-CM5A-MR'
# )

PPE_MEMBERS <- c(
  1,
  4,
  5,
  6,
  7,
  8,
  9,
  10,
  11,
  12,
  13,
  15
)

DSI_NUMBER <- c(
  3,
  6,
  12
)


domain_to_res <- list(
  global = '60km',
  regional = '12km',
  local = '2.2km'
)

ENSEMBLE = c(
  'CMIP',
  'PPE',
  'full'
)

ensemble_length <- list(
  CMIP = 4,
  PPE = 12,
  full = 16
)

COUNTRIES <- list(
  'uk',
  'Scotland',
  'England',
  'Wales',
  'Northern_Ireland'
)

get_ensemble_model_data <- function(country, domain, index_number, ensemble){
  firstmonth <- str_pad(index_number - 1, width = 2, pad = "0")
  file_path <- paste0('/data/users/mgrant/ukcp/droughts/data/', country, '/',
                      domain, '/ensemble_cubes/', ensemble, '/TS1/monthly/dsi-', 
                      index_number, '/', ensemble, '_ensemble_1981', firstmonth, 
                      '-200011_dsi-', index_number, '.nc')
  cube <- nc_open(file_path)
  data <- ncvar_get(cube)
  
  # Remove masked values
  data[data > 10000] <- NA
  
  return(data)
}

get_ensemble_array_at_ij <- function(ensemble_data, i, j){
  data_at_ij = ensemble_data[i, j,,]
  
  
  if (dim(data_at_ij[2])>ensemble_length[[ensemble]]){
    print('Error: Likely wrong dimension indexed')
  }
  
  return(data_at_ij)
}


get_obs_data <- function(country, index_number, resolution){
  obs_cube <- load_obs_cube(country, index_number, resolution)
  data <- ncvar_get(obs_cube)
  data[data > 10000] <- NA
  
  return(data)
}


# Set up to parallelise the code
# Define the number of cores/workers to use
# num_cores <- 6
# 
# # Register parallel backend using doParallel
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)


# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }


find_rank_hist <- function(country, domain, index_number, ensemble){
  # Load in the model ensemble and obs data
  ensemble_data <- get_ensemble_model_data(
    country, domain, index_number, ensemble)
  
  res = domain_to_res[[domain]]
  obs_data <- get_obs_data(country, index_number, res)

  rank_hist = numeric(dim(ensemble_data)[4] + 1)
  for (i in 1:dim(ensemble_data)[1]) {
    for (j in 1:dim(ensemble_data)[2]) {
      obs_at_ij <- obs_data[i, j,]
      model_at_ij = ensemble_data[i, j,,]
      
      if (any(!is.na(model_at_ij))==FALSE){
        next
      }
      
      if (any(!is.na(obs_at_ij))==FALSE){
        next
      }
      
      if (dim(model_at_ij)[2] != ensemble_length[[ensemble]]){
        stop('Error: Likely wrong dimension indexed')
      }
      
      rank_hist_at_ij = Rankhist(model_at_ij, obs_at_ij)
      rank_hist <- rank_hist + rank_hist_at_ij
    }
  }
  return(rank_hist)
}

# 
# find_rank_hist_vectorised <- function(domain, index_number, ensemble){
#   # Load in the model ensemble and obs data
#   ensemble_data <- get_ensemble_model_data(domain, index_number, ensemble)
#   
#   res <- domain_to_res[[domain]]
#   obs_data <- get_obs_data(index_number, res)
#   
#   rank_hist <- matrix(0, nrow = dim(ensemble_data)[1], ncol = dim(ensemble_data)[2])
#   
#   valid_indices <- apply(ensemble_data, c(1, 2), function(x) dim(x)[2] == ensemble_length[[ensemble]])
#   
#   rank_hist <- apply(ensemble_data, c(1, 2), function(x) sum(Rankhist(x, obs_data)))
#   
#   rank_hist <- colSums(rank_hist) + 1
#   
#   return(rank_hist)
# }

# rank_hist <- find_rank_hist('global', 3, 'CMIP')
# rank_hist_vectorised <- find_rank_hist('global', 3, 'CMIP')
# 
# print(rank_hist - rank_hist_vectorised)

# unregister_dopar()


plot_rank_hist <- function(country, domain, index_number, ensemble){
  rank_hist <- find_rank_hist(country, domain, index_number, ensemble)
  
  filepath <- paste0('/scratch/mgrant/UKCP/droughts/plots/',
                     'rank_histograms/', country, '/', domain, '/dsi-', 
                     index_number, '/')
  
  if (!dir.exists(filepath)) {
    dir.create(filepath, recursive=TRUE)
  }
  
  filename <- paste0('rank_histogram_dsi-', index_number, '_', ensemble, 
                     '_ensemble_', domain, '_domain.png')
  
  png(paste0(filepath, filename), width = 600, height = 600)
  
  par(mar=c(5, 5, 3, 1), cex.main=4, cex.lab=2.5, cex.axis=2)
  bp <- barplot(rank_hist, xlab=NA, ylab=NA, col='coral1', axes=FALSE)
  axis(1, at=bp, line=.2, labels=paste(1:length(rank_hist)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  if (ensemble == 'full'){
    ensemble_title = 'Full'
  } else {
    ensemble_title = ensemble
  }
  
  main_title = paste0('Rank Histogram of DSI-', index_number, ' in ', country, 
                      ' for the ', ensemble_title, ' Ensemble \n on the ', 
                      domain, ' Domain')
  par(oma = c(1, 1, 2.2, 0))
  title(main = main_title, outer = TRUE)
  dev.off()
}


plot_multiple_rank_hists <- function(country, index_number){
  
  rank_hist_cmip_global = find_rank_hist(country, 'global', 3, 'CMIP')
  rank_hist_cmip_regional = find_rank_hist(country, 'regional', 6, 'CMIP')
  rank_hist_cmip_local = find_rank_hist(country, 'local', 12, 'CMIP')
  rank_hist_ppe_global = find_rank_hist(country, 'global', 3, 'PPE')
  rank_hist_ppe_regional = find_rank_hist(country, 'regional', 6, 'PPE')
  rank_hist_ppe_local = find_rank_hist(country, 'local', 12, 'PPE')
  rank_hist_full_global = find_rank_hist(country, 'global', 3, 'full')
  rank_hist_full_regional = find_rank_hist(country, 'regional', 6, 'full')
  rank_hist_full_local = find_rank_hist(country, 'local', 12, 'full')
  
  par(mfrow(3, 3))
  bp1 <- barplot(rank_hist_cmip_global, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp1, line=.2, labels=paste(1:length(rank_hist_cmip_global)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_cmip_global), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp2 <- barplot(rank_hist_cmip_regional, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp2, line=.2, labels=paste(1:length(rank_hist_cmip_regional)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_cmip_regional), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp3 <- barplot(rank_hist_cmip_local, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp3, line=.2, labels=paste(1:length(rank_hist_cmip_local)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_cmip_local), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp4 <- barplot(rank_hist_ppe_global, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp4, line=.2, labels=paste(1:length(rank_hist_ppe_global)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_ppe_global), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp5 <- barplot(rank_hist_ppe_regional, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp5, line=.2, labels=paste(1:length(rank_hist_ppe_regional)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_ppe_regional), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp6 <- barplot(rank_hist_ppe_local, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp6, line=.2, labels=paste(1:length(rank_hist_ppe_local)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_ppe_local), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp7 <- barplot(rank_hist_full_global, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp7, line=.2, labels=paste(1:length(rank_hist_full_global)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_full_global), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp8 <- barplot(rank_hist_full_regional, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp8, line=.2, labels=paste(1:length(rank_hist_full_regional)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_full_regional), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  bp9 <- barplot(rank_hist_full_local, xlab=NA, 
                 ylab=NA, col='coral1', axes=FALSE,
                 cex.main=4, cex.lab=2.5, cex.axis=2)
  axis(1, at=bp9, line=.2, labels=paste(1:length(rank_hist_full_local)))
  axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist_full_local), 
       labels=c(0, 20, 40, 60, 80, 100))
  mtext(side=1, text="Rank", line=2)
  mtext(side=2, text="Percentage of Max Bin", line=2)
  
  
  main_title = paste0('Rank Histogram of DSI-', index_number, ' in ', country, 
                      ' for All Ensembles')
  par(oma = c(1, 1, 2.2, 0))
  title(main = main_title, outer = TRUE)
  
}


plot_multiple_rank_hists('uk', 3)

# save_rank_hists <- function(country, domain)
# plot_rank_hist('global', 3, 'CMIP')

# stopCluster(cl)
# plot_all_rank_hists <- function(){
#   foreach(domain=DOMAIN) %dopar% {
#     foreach(number=DSI_NUMBER) %dopar% {
#       foreach(ensemble=ENSEMBLE) %dopar% {
#         plot_rank_hist(domain, number, ensemble)
#       }
#     }
#   }
# }

plot_all_rank_hists <- function(){
  for (country in COUNTRIES){
    for (domain in DOMAIN) {
      for (number in DSI_NUMBER) {
        for (ensemble in ENSEMBLE) {
          plot_rank_hist(country, domain, number, ensemble)
        }
      }
    }
  }
}

# plot_all_rank_hists()


# rank_hist <- find_rank_hist('regional', 3, 'CMIP')
# 
# png('/home/h02/mgrant/ukcp/ukcp_2.0/droughts/test_RH.png', width = 400, height = 400)
# # PlotRankhist(rank_hist)
# 
# par(mar=c(5, 5, 3, 1), cex.lab=1, cex.axis=1)
# bp <- barplot(rank_hist, xlab=NA, ylab=NA, col='coral1', axes=FALSE)
# axis(1, at=bp, line=.2, labels=paste(1:length(rank_hist)))
# axis(2, las=1, at=c(0, .2, .4, .6, .8, 1)*max(rank_hist), 
#      labels=c(0, 20, 40, 60, 80, 100))
# mtext(side=1, text="Rank", line=2)
# mtext(side=2, text="Percentage of Max Bin", line=2)
# 
# # par(mgp = c(1, 5, 4))
# main_title = 'Rank Histogram of DSI-3 for the CMIP5 Ensemble \n on the Regional Domain'
# # par(cex.main = 3)
# # par(mar = c(5, 5, 5, 5))
# par(oma = c(1, 1, 2.2, 0))
# title(main = main_title, outer = TRUE)
# dev.off()
# 
# model_cube <- load_model_cube('global', 'MRI-CGCM3', 3)
# model_data <- ncvar_get(model_cube)
# model_data[model_data > 1000] <- NA
# print(model_data)
# obs_cube <- load_obs_cube(3, '60km')
# 
# model_data <- get_unmasked_nonzero_data(model_cube)
# obs_data <- get_unmasked_nonzero_data(obs_cube)


