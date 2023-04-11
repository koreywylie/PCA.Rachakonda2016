# Principal Components Analyses application,
#   fns. that compress & de-compress subj.- & group-level data
#   as a prerequisite for group Independent Components Analyses or Hierarchical Principal Components Analysis.
#
# Requires:
#   PCA_fns.R : source code for optimized PCA methods,
#                 including ability to mask 3D fMRI data,
#                 simple PCA algorithms,
#                 & a random PCA algorithm parallelized for large datasets (MPOWIT).
#               https://github.com/koreywylie/PCA.Rachakonda2016
#
# Referrences:
#   Erhardt et al. (2011). Comparison of multi-subject ICA methods
#     for analysis of fMRI data. Hum Brain Mapp 32(12), p. 2075-2095.
#


##################################################################
detrend <- function(ts, detrend_rows = T){
  ##################################################################
  # Remove linear trend from time series'
  #    (code from stats::spec.pgram)
  #
  # Input:
  #   ts : time series matrix
  #   detrend_rows : detrend rows of matrix if T, otherwise detrend cols. 
  # Output:
  #   ts : as above, after subtracting global mean at each time point
  #
  
  if (detrend_rows){ ts = t(ts) }
  
  N <- nrow(ts)
  t <- 1L:N - (N + 1)/2
  sumt2 <- N * (N^2 - 1)/12
  for (i in 1L:ncol(ts)) ts[, i] <- ts[, i] - mean(ts[, i], na.rm=T) - 
    sum(ts[, i] * t, na.rm=T) * t/sumt2
  
  if (detrend_rows){ ts = t(ts) }
  return(ts)
} ##################################################################


##################################################################
subj_preprocData <- function(data, 
                             preproc.type = 'rm.voxel.mean', 
                             discard.n = NA,
                             var = NA, var.subset = NA, mask = NA, 
                             verbose = F){
  ##################################################################
  # Pre-processing data prior to subj.-level PCA,
  #   called as part of subj_PCA() & groupICA_BackReconstruction().
  #     Recommended for PCA prior to & part of ICA
  #
  # Input:
  #   data : (options)
  #           1. path to subj. file, w/ data 'var' as [voxels x time] matrix
  #           2. matrix w/ dims. [voxels x time] matrices
  #   preproc.type : pre-processing type (options, or skip if NA)
  #       rm.voxel.mean  : center all voxel time series' to mean = 0
  #       rm.time.mean   : center each timepoint (~vol.) to mean = 0
  #       variance.norm  : detrend all voxels & scale to standard deviation
  #   var  : if 'data' w/ option 1, name of saved variable w/ subj. data
  #   var.subset : if 'var' is a list, name of list element to load
  #   mask : if 'data' is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   discard.n : if 'data' is a Nifti file(s), 
  #                 number of initial vols. to discard to allow for magnet stabilization
  # Output:
  #   data  : [voxels x time] matrix pre-processed as above
  # Requires:
  #   load.data() : data munging fn. from PCA_fns.R
  #   detrend()   : detrending for during variance normalization
  #
  
  if (verbose){  cat('\n......pre-processing subj. data:  ')}
  
  preproc.type = as.character(preproc.type)
  discard.n = as.integer(discard.n)
  data = load.data(data, var, var.subset, mask)
  
  if (!all(is.na(discard.n)) && discard.n > 0){
    if (verbose){  cat(paste0('discarding initial ', discard.n, ' volumes...'))}
    n1 = discard.n+1
    p = ncol(data)
    data = data[,n1:p]
  }
  
  if (all(is.na(preproc.type))){
    return(data)
    
  }else if (!(preproc.type %in% c('rm.voxel.mean', 'rm.time.mean', 'variance.norm'))){
    if (verbose){
      cat('\n  WARNING: could not understand preprocessing type,')
      cat(' defaulting to removing voxel means')
    }
    return(sweep(data, 1, rowMeans(data)))
    
  }else{
    if (preproc.type == 'rm.voxel.mean'){
      if (verbose){  cat('removing voxel means')  }
      return(sweep(data, 1, rowMeans(data)))
      
    }else if (preproc.type == 'rm.time.mean'){
      if (verbose){  cat('removing timepoint means')  }
      return(sweep(data, 2, colMeans(data)))
      
    }else if (preproc.type == 'variance.norm'){
      if (verbose){  cat('detrending voxels & scaling to s.d.')  }
      
      data = detrend(data, detrend_rows=T)
      data = t(apply(data, 1, function(yy) return(yy / sd(yy, na.rm=T))))
      
      if (any(is.na(data))){
        if (verbose){
          NA.any = sum(is.na(data))
          NA.byrow = sum(apply(data,1,function(yy) any(is.na(yy))))
          cat(paste0('\n\nWARNING: Found  ',NA.any,'  NA values in  ',NA.byrow,'/',nrow(data),'  rows in data,'))
          cat(paste0('...setting all NA values to 0...\n\n'))
        }
        data[is.na(data)] = 0
      }
      if (any(is.infinite(data))){
        if (verbose){
          Inf.any = sum(if.infinite(data))
          Inf.byrow = sum(apply(data,1,function(yy) any(is.infinite(yy))))
          cat(paste0('\n\nWARNING: Found  ',Inf.any,'  Infinite values in  ',Inf.byrow,'/',nrow(data),'  rows in data,'))
          cat(paste0('...setting all +/-Inf values to 0...\n\n'))
        }
        data[is.infinite(data)] = 0
      }
      
      return(data)
    }
  }
} ##################################################################


##################################################################
subj_PCA <- function(data, K1 = NA, 
                     preproc.type = 'rm.voxel.mean',
                     discard.n = NA, whitenPCs = F,
                     data_wc = NA, 
                     var = NA, var.subset = NA, mask = NA,
                     prefix = NA, save_path = NA,
                     return.data = F, parallel = F, verbose = F, ...){
  ##################################################################
  # Performances subject-level dimension reduction w/ PCA
  #
  # Input:
  #   data : (options)
  #           1. vector of subjs. saved as .RData files, 
  #               w/ data 'var' as [voxels x time] matrix
  #           2. path to directory containing above
  #           3. list of [voxels x time] matrices
  #           4. array w/ dims. [subjs. x voxels x time]
  #   K1   : PCA subspace dimension,
  #           if NA defaults to min. dim. in data
  #   preproc.type : pre-processing type  (options, or skip if NA)
  #       rm.voxel.mean  : center all voxel time series' to mean = 0
  #       rm.time.mean   : center each timepoint (~vol.) to mean = 0
  #       variance.norm  : detrend all voxels & scale to standard deviation
  #   discard.n : if 'data' is a Nifti file(s), 
  #                 number of initial vols. to discard to allow for magnet stabilization
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   data_wc   : if 'data' w/ option 2, wildcard to filter subj. files
  #   var       : if 'data' w/ option 1, name of saved variable w/ subj. data
  #   var.subset: if 'var' is a list, name of list element to load
  #   mask      : optional binary mask if 'data' is a Nifti file(s),
  #                 or path to such for .RData file(s) to use for back-reconstruction
  #   prefix    : identifying prefix for experiment, 
  #                 attached to saved output filename
  #   save_path   : path to dir. for saved files
  #   return.data : return list of subj. PCA spaces if T,
  #                   or path to dir. w/ saved files if F
  #   ...       : additional inputs passed to specific PCA fns.
  #                 (e.g. S_g=10 for group size in PCA.SubsampledTime())
  # Output: 
  #   saved files w/ *_pca1.RData added as suffix & vars.:
  #     pca$X  : [voxels x K1] matrix of subj. PCA space
  #     pca$U  : [time   x K1] reducing matrix of eigenvectors
  #     pca$Lambda : eigenvalues
  # Requires:
  #   subj_preprocData()   : pre-processes data
  #   load.data()          : data munging from PCA_fns.R
  #   detrend()            : detrending for variance normalization
  #   PCA.Rachakonda2016() : main fn. in PCA_fns.R
  #
  
  if (verbose){  cat('\nCalculating subj.-level PCA spaces...')}
  if (parallel){  library(doParallel)  }
  
  ### Defaults ###
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources
  
  pca_algorithm = NA  # Initialize computation settings for PCA...
  pca_unstacked = NA  # ...as determined by 1st subj. & continued thereafter
  
  ### Inputs ###
  if (is.character(data) && (length(data)==1) && dir.exists(data)){
    stopifnot(length(list.files(data)) > 0)
    data = list.files(data, full.names=T)
    data = data[!dir.exists(data)] # exclude subdirs.
    
    if (!all(is.na(data_wc))){     # filter files
      data = data[grepl(data_wc, data)]
    }
  }
  if (is.vector(data) | is.list(data)){
    S = length(data)
  }else if (is.array(data)){
    S = dim(data)[1]
  }else{
    print('WARNING: Could not understand input data!')
    S = length(data)
  }
  
  if (!all(is.na(var)) && is.character(data) && !any(grepl('.RData$', data))){
    var = NA   # only applicable for data in saved R files
    var.subset = NA
  }
  if (!all(is.na(mask)) && is.character(mask)){
    mask.file = mask
  }else{    # check if mask.file is saved in .RData file
    mask.file = NA
    if (!all(is.na(var)) && is.character(data) && file.exists(data[1]) && 
        grepl('.RData$', data[1])){
      Space.subj_data = new.env()
      load(data[1], envir=Space.subj_data)
      
      if (any(grepl('mask', names(Space.subj_data)))){
        subj.mask.vars = names(Space.subj_data)[grepl('mask', names(Space.subj_data), ignore.case=T)]
        for (mask.var.name in subj.mask.vars){
          mask.var = get(mask.var.name, envir=Space.subj_data)
          
          if (!all(is.na(mask.var)) && is.character(mask.var) &&  # check mask.var for valid nifti vol.
              grepl('\\.(nii|img|hdr)', mask.var) && file.exists(mask.var)){
            mask.file = mask.var          # will be saved for use during post-PCA back-reconstruction, following subj. PCA below
          }
        }
      }
      rm(Space.subj_data)
    }
  }
  
  K1 = as.integer(K1)
  if (verbose){  cat(paste0('\n...reduced dim. K1 = ',K1))  }
  whitenedPCs = whitenPCs = as.logical(whitenPCs)
  if (verbose && whitenPCs){cat(paste0('\n...all Principal Comps. will be whitened by scaling to sqrt. of eigenvalues...'))}
  prefix = as.character(prefix)
  save_path = as.character(save_path)
  if (!all(is.na(save_path)) && !dir.exists(save_path)){
    dir.create(save_path, recursive=T)
  }
  if (return.data){  PCA_s = list() }
  
  
  ### Main fn. ###
  if (parallel){
    ### Get best algorithm & run all subjs. w/ same settings ###
    recommendations = find_PCAtype_best(data[1], k=K1, 
                                        var=var, var.subset=var.subset, mask=mask,
                                        verbose=verbose)
    pca_algorithm = recommendations$PCAtype
    pca_unstacked = recommendations$unstacked
    pca_fns = PCA.Rachakonda2016(NA, return.fns.list=T)
    
    if (verbose && is.character(data)){
      for (s in 1:S){  cat(paste0('\n...PCA compressing (in parallel):  ', basename(data[s])))  }
    }
    PCA_s = foreach(s=icount(S),
                    .export=c('load.data', 'calc_Yi_means', 'zeroMean_Yi', 'whitenedPCs',
                              'detrend', 'subj_preprocData',
                              'PCA.Rachakonda2016', pca_fns,
                              'flatten_img', 'concatTemporal_Yi')) %dopar% {
                                if (is.character(data) && all(file.exists(data))){
                                  ### Create save file name ###
                                  data_files = data[s]  # plural for consistency w/ other fns.
                                  save.fname = basename(data_files)
                                  if (!is.na(prefix) && !grepl(paste0('^',prefix), save.fname)){
                                    save.fname = paste0(prefix, save.fname)
                                  }
                                  save.fname = sub('.nii$', '', save.fname)
                                  save.fname = sub('.RData$', '', save.fname)
                                  save.fname = paste0(save.fname, subj_pca_suffix)
                                  if (!is.na(save_path)){
                                    save.fname = file.path(save_path, save.fname)
                                  }else{
                                    save.fname = file.path(dirname(data_files), save.fname)
                                  }
                                  if (file.exists(save.fname) && !return.data){
                                    cat(paste0('\n.....found complete subj. PCA: ', save.fname,' continuing to next subj...'))
                                    return(invisible(NULL))
                                  }
                                  
                                  ### Load & pre-process data ###
                                  data.preproc = subj_preprocData(data[s], 
                                                                  preproc.type, discard.n,
                                                                  var=var, var.subset=var.subset, mask=mask,
                                                                  verbose=verbose)
                                  
                                }else if (is.list(data)){
                                  ### Create save file name ###
                                  data_files = NA
                                  save.fname = paste0('s',s,subj_pca_suffix)
                                  if (!is.na(prefix) && !grepl(paste0('^',prefix), save.fname)){
                                    save.fname = paste0(prefix, save.fname)
                                  }
                                  if (!is.na(save_path)){
                                    save.fname = file.path(save_path, save.fname)
                                  }else{
                                    save.fname = file.path(dirname(data_files), save.fname)
                                  }
                                  if (file.exists(save.fname) && !return.data){
                                    cat(paste0('\n.....found complete subj. PCA: ', save.fname,' continuing to next subj...'))
                                    return(invisible(NULL))
                                  }
                                  
                                  ### Load & pre-process data ###
                                  data.preproc = subj_preprocData(data[[s]], 
                                                                  preproc.type, discard.n,
                                                                  var=var, var.subset=var.subset, mask=mask,
                                                                  verbose=verbose)
                                }
                                
                                
                                ### Main fn.:  Subj-level PCA ###
                                pca = PCA.Rachakonda2016(data.preproc, K1, 
                                                         var=var, var.subset=var.subset, mask=mask,
                                                         whitenPCs=whitenPCs,
                                                         algorithm=pca_algorithm,
                                                         unstacked=pca_unstacked,
                                                         verbose=verbose, ...)
                                
                                
                                ### Saving ###
                                mask = mask.file
                                if (verbose){  cat(paste0('\n......saving as: ', save.fname))}
                                save(list=c('pca', 'K1', 'preproc.type', 'discard.n', 'whitenedPCs',
                                            'data_files', 'var', 'var.subset', 'mask', 'prefix', 's'), 
                                     file=save.fname)
                                stopifnot(file.exists(save.fname))
                                if (return.data){ return(pca) }
                              }
    
  }else{
    for (s in 1:S){
      if (verbose){  cat(paste0('\n...PCA compressing:  ', basename(data[s])))  }
      
      if (is.character(data) && all(file.exists(data))){
        ### Create save file name ###
        data_files = data[s]  # plural for consistency w/ other fns.
        save.fname = basename(data_files)
        if (!is.na(prefix) && !grepl(paste0('^',prefix), save.fname)){
          save.fname = paste0(prefix, save.fname)
        }
        save.fname = sub('.nii$', '', save.fname)
        save.fname = sub('.RData$', '', save.fname)
        save.fname = paste0(save.fname, subj_pca_suffix)
        if (!is.na(save_path)){
          save.fname = file.path(save_path, save.fname)
        }else{
          save.fname = file.path(dirname(data_files), save.fname)
        }
        if (file.exists(save.fname) && !return.data){
          cat(paste0('\n.....found complete subj. PCA: ', save.fname,' continuing to next subj...'))
          next
        }
        
        ### Load & pre-process data ###
        data.preproc = subj_preprocData(data[s], 
                                        preproc.type, discard.n,
                                        var=var, var.subset=var.subset, mask=mask,
                                        verbose=verbose)
        
      }else if (is.list(data)){
        ### Create save file name ###
        data_files = NA
        save.fname = paste0('s',s,subj_pca_suffix)
        if (!is.na(prefix) && !grepl(paste0('^',prefix), save.fname)){
          save.fname = paste0(prefix, save.fname)
        }
        if (!is.na(save_path)){
          save.fname = file.path(save_path, save.fname)
        }else{
          save.fname = file.path(dirname(data_files), save.fname)
        }
        if (file.exists(save.fname) && !return.data){
          cat(paste0('\n.....found complete subj. PCA: ', save.fname,' continuing to next subj...'))
          next
        }
        
        ### Load & pre-process data ###
        data.preproc = subj_preprocData(data[[s]], 
                                        preproc.type, discard.n,
                                        var=var, var.subset=var.subset, mask=mask,
                                        verbose=verbose)
      }
      
      
      ### Main fn.:  Subj-level PCA ###
      #     Automatically selects PCA algorithm,
      #       will run unstacked in parallel, if nodes available
      
      pca = PCA.Rachakonda2016(data.preproc, K1, 
                               var=var, var.subset=var.subset, mask=mask,
                               whitenPCs=whitenPCs,
                               algorithm=pca_algorithm,
                               unstacked=pca_unstacked,
                               verbose=verbose, ...)
      
      ### Run all subsequent subjs. w/ same settings ###
      pca_algorithm = pca$algorithm  # initialized to NA above...
      pca_unstacked = pca$unstacked  # ...then determined by 1st subj. & continued thereafter
      
      
      ### Saving ###
      mask = mask.file
      if (verbose){  cat(paste0('\n......saving as: ', save.fname))}
      save(list=c('pca', 'K1', 'preproc.type', 'discard.n', 'whitenedPCs',
                  'data_files', 'var', 'var.subset', 'mask', 'prefix', 's'), 
           file=save.fname)
      stopifnot(file.exists(save.fname))
      if (return.data){ PCA_s[[s]] = pca  }
    }
  }
  if (verbose){  cat('\n')  }
  
  ##################################################################
  if (return.data){  return(PCA_s)  }   # list of subj.-specific PCAs
} ##################################################################


##################################################################
get_subjPCA_dat <- function(subj_pca_file, prefix = NA, 
                            return.PCs = T, return.orig.data = F,
                            verbose = F){
  ##################################################################
  # Loads & formats subj-level PCA data,
  #   & applies PCA pre-processing to original data as output of subj_PCA()
  #
  # Requires:
  #   subj_preprocData()   : pre-processes data
  #   load.data()          : data munging fn. from PCA_fns.R
  #   detrend()            : detrending for pre-processing
  #
  
  if (verbose){cat(paste0('\n...loading subj.-level PCA from:\n      ',subj_pca_file))}
  
  Space.subj_pca = new.env()
  load(subj_pca_file, envir=Space.subj_pca)
  stopifnot(all(c('pca', 'K1', 
                  'data_files', 'var', 'var.subset', 
                  'prefix', 'mask') %in% names(Space.subj_pca)))
  
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.subj_pca)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.subj_pca))
  }
  subj_data_file = get('data_files', envir=Space.subj_pca)
  stopifnot(is.character(subj_data_file) 
            && file.exists(subj_data_file))
  var = get('var', envir=Space.subj_pca)
  var.subset = get('var.subset', envir=Space.subj_pca)
  mask = get('mask', envir=Space.subj_pca)
  preproc.type = get('preproc.type', envir=Space.subj_pca)
  discard.n = get('discard.n', envir=Space.subj_pca)
  subjPCA = get('pca', envir=Space.subj_pca)
  
  subjWhitenedPCs = get('whitenedPCs', envir=Space.subj_pca)
  
  ### Get time series info for original data ###
  good.timepoints = NA     # indices of non-scrubbed, non-excluded time points in original time series
  total.timepoints = NA    # total length of original time series
  if (grepl('\\.RData$', subj_data_file)){
    Space.subj_orig = new.env()
    load(subj_data_file, envir=Space.subj_orig)
    stopifnot(var %in% names(Space.subj_orig))
    
    if ('good.timepoints' %in% names(Space.subj_orig)){
      good.timepoints = get('good.timepoints', envir=Space.subj_orig)
    }
    if ('total.timepoints' %in% names(Space.subj_orig)){
      total.timepoints = get('total.timepoints', envir=Space.subj_orig)
    }
  }
  
  ### Get original data used for PCA ###
  data.preproc = NA        # time series used for PCA
  voxel.data.preproc = NA  # voxel-level time series, if different from above
  if (return.orig.data){
    ### Apply PCA pre-processing to data ###
    data.preproc = subj_preprocData(subj_data_file,
                                    preproc.type, discard.n,
                                    var=var, var.subset=var.subset, mask=mask,
                                    verbose=verbose)
    data.preproc = zeroMean_Yi(data.preproc)
    
    ### Check for voxel time series (~denoised, pre-ROI)  ###
    if (grepl('\\.RData$', subj_data_file)){
      Space.subj_orig = new.env()
      load(subj_data_file, envir=Space.subj_orig)
      stopifnot(var %in% names(Space.subj_orig))
      
      if ('mask' %in% names(Space.subj_orig) &&
          !all(is.na(get('mask', envir=Space.subj_orig)))){
        mask = get('mask', envir=Space.subj_orig)
      }
      
      if (all(c('Vx.ts', 'mask') %in% names(Space.subj_orig))){
        voxel.data.preproc = subj_preprocData(get('Vx.ts', envir=Space.subj_orig),
                                              preproc.type, discard.n,
                                              var=NA, var.subset=NA, mask=mask,
                                              verbose=verbose)
        voxel.data.preproc = zeroMean_Yi(voxel.data.preproc)
      }
    }
    
    ### Default, implied time series info for original data ###
    if (all(is.na(good.timepoints))){
      good.timepoints = c(1:ncol(data.preproc))
      if (!all(is.na(discard.n))){
        good.timepoints = good.timepoints + discard.n
      }
    }
    if (all(is.na(total.timepoints))){
      total.timepoints = ncol(data.preproc)
      if (!all(is.na(discard.n))){
        total.timepoints = total.timepoints + discard.n
      }
    }
  }
  
  return.PCs = as.logical(return.PCs)
  if (!return.PCs){  subjPCA$X = NA  }
  
  ##################################################################
  return(list('Y'   = data.preproc,     # pre-processed subj. data used for PCA
              'Y_vx'= voxel.data.preproc, # pre-processed, pre-ROI averaged voxel' t.s.
              'X'   = subjPCA$X,        # subj. Principal Components
              'U_i' = subjPCA$U,        # subj.-level PCA eigenvectors
              'L_i' = subjPCA$Lambda,   # subj.-level PCA eigenvalues
              'whitenedPCs' = subjWhitenedPCs, # if Princ.Comps. scaled to eigenvals.
              'discard.n'       = discard.n,   # number of initial time points dropped
              'good.timepoints' = good.timepoints, # indices of non-scrubbed, non-excluded time points
              'total.timepoints'= total.timepoints, # total number of time points in original t.s.
              'subj_pca_file' = subj_pca_file, # path to subject PCA file
              'prefix'        = prefix,        # file prefix for experiment
              'mask'          = mask))         # path to mask vol.
} ##################################################################



##################################################################
group_PCA <- function(data, K2 = NA, 
                      whitenPCs = F,
                      prefix = NA, save_path = NA,
                      big.data = F, return.data = F, verbose = F, ...){
  ##################################################################
  # Performances group-level dimension reduction w/ PCA
  #
  # Input:
  #   data : (options)
  #           1. vector of subj. PCA files, output from subj_PCA()
  #           2. path to directory containing above
  #           3. list of [voxels x time subspace] matrices
  #           4. array w/ dims. [subjs. x voxels x time subspace]
  #   K2   : reduced group subspace dimension,
  #           if NA defaults to min. dim. in data
  #   whitenPCs   : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   prefix      : identifying prefix for experiment, 
  #                  attached to saved output filename
  #   save_path   : path to dir. for saved files
  #   ...  : additional inputs passed to specific PCA fns.
  #           (e.g. S_g=20 for group size in PCA.SubsampledTime())
  # Output: 
  #   saved files w/ *_groupPCA.RData added as suffix & vars.:
  #     pca$X  : [voxels x K2] matrix of group PCA space
  #     pca$U  : [concatenated ~time x K2] reducing matrix of eigenvectors
  #     pca$Lambda : eigenvalues
  # Requires:
  #   get_subjPCA_dat()    : loads subj.-level PCs
  #   find_PCAtype_best()  : determines best PCA type based on data size & mem. requirements, from PCA_fns.R
  #   get_RAMavailable()   : determines available memory, from PCA_fns.R
  #   PCA.Rachakonda2016() : main PCA fn. from PCA_fns.R
  #
  
  if (verbose){  cat('\nCalculating group-level PCA space...')}
  
  ### Inputs & defaults ###
  subj_pca_suffix  = '_PCA1.RData'       # filesuffix for subj. temporal PCAs from subj_PCA()
  group_pca_suffix = 'groupPCA.RData'    # filesuffix for group-level temporal PCA from group_PCA()
  group_ica_suffix = 'spatialICA.RData'  # filesuffix for group spatial ICA from group_spatialICA()
  br_suffix        = '_ICAsources.RData' # filesuffix for back-reconstructed subj. sources  
  
  prefix = as.character(prefix)
  K2 = as.integer(K2)
  if (verbose){  cat(paste0('\n...reduced dim. K2 = ',K2))  }
  whitenedPCs = whitenPCs = as.logical(whitenPCs)
  if (verbose && whitenPCs){cat(paste0('\n...all Principal Comps. will be whitened by scaling to sqrt. of eigenvalues...'))}
  
  ### Prelims ###
  subj_pca_wc = paste0('*', subj_pca_suffix)
  var = 'pca'      # format of saved output from subj_PCA()
  var.subset = 'X' # name of list element w/ PCs within above 
  algorithm = NA   # PCA algorithm selected based on data size & available memory below
  unstacked = !big.data   # defaults to loading all datasets into memory if True
  mask = NA        # path to mask vol., shared among subjects
  
  if (is.character(data) && all(file.exists(data))){
    if ((length(data)==1) && dir.exists(data)){
      if (is.na(save_path)){
        save_path = data
      }
      data_files = list.files(data, subj_pca_wc, full.names=T)
    }else{
      data_files = data
      if (is.na(save_path)){
        save_path = file.path(unique(dirname(data_files))[1])
      }
    }
    stopifnot(length(data_files) > 0)
    
    ### Load Mask info ###
    Space.subj_pca = new.env()
    load(data_files[1], envir=Space.subj_pca)
    if ('mask' %in% names(Space.subj_pca)){
      mask = get('mask', envir=Space.subj_pca)
    }
    rm(Space.subj_pca)
    
    
    if (big.data){  # avoid loading full dataset into memory
      data.dims = c(length(data_files),  # number of subjs.
                    dim(get_subjPCA_dat(data_files[1], prefix=prefix, 
                                        return.PCs=T, verbose=F)$X))
      recommendations = find_PCAtype_best(data.dims, K2,
                                          verbose=verbose)
      data = data_files
      algorithm = recommendations$PCAtype
      unstacked = recommendations$unstacked
      
    }else{
      data = vector('list', length(data_files))
      if (verbose){
        cat('\n...loading subj.-level PCA data from:')
        cat(paste0('\n      ', file.path(unique(dirname(data_files))[1], subj_pca_wc)))
      }
      for (f in 1:length(data_files)){
        data[[f]] = get_subjPCA_dat(data_files[f], prefix=prefix, 
                                    return.PCs=T, verbose=F)$X
      }
      var = NA
      var.subset = NA
    }

  }else{
    data_files = NA
  }
  
  ### Main fn. ###
  if (verbose && !all(is.na(data_files)) && is.character(data_files)){
    cat(paste0('\n...PCA compressing (as group):'))
    for (s in 1:length(data_files)){  cat(paste0('\n  ', basename(data_files[s])))  }
  }
  pca = PCA.Rachakonda2016(data, K2,       # Automatically selects PCA algorithm, will run in parallel if available
                           var = var, var.subset = var.subset, 
                           whitenPCs = whitenPCs, 
                           algorithm = algorithm, unstacked = unstacked,
                           verbose = verbose, ...)
  
  ### Saving ###
  save.fname = group_pca_suffix
  if (!is.na(prefix) && !grepl(paste0('^',prefix), save.fname)){
    save.fname = paste0(prefix, '_', save.fname)
  }
  if (!is.na(save_path)){
    save.fname = file.path(save_path, save.fname)
  }
  if (verbose){  cat(paste0('\n......saving as: ', save.fname))}
  save(list=c('pca', 'data_files', 'K2', 
              'whitenedPCs', 'prefix', 'mask'), file=save.fname)
  stopifnot(file.exists(save.fname))
  
  if (verbose){  cat('\n\n')  }
  ##################################################################
  if (return.data){  return(pca)  } # group PCA output
} ##################################################################


##################################################################
get_groupPCA_dat <- function(group_pca_file, prefix = NA, 
                             model_orders = NA, return.PCs = T,
                             verbose = F){
  ##################################################################
  # Loads & formats group-level PCA,
  #   output by group_PCA(),
  #     prior to ICA back-reconstruction
  #
  
  if (verbose){cat(paste0('\n...loading group-level PCA data from:\n      ',
                          group_pca_file))}
  
  return.PCs = as.logical(return.PCs)
  if (all(is.na(return.PCs))){  return.PCs = F  }
  
  Space.pca = new.env()
  stopifnot(is.character(group_pca_file) 
            && file.exists(group_pca_file))
  load(group_pca_file, envir=Space.pca)
  stopifnot(all(c('pca', 'data_files', 
                  'K2', 'prefix') %in% names(Space.pca)))
  if (all(is.na(prefix))){
    prefix = get('prefix', envir=Space.pca)
  }else{
    stopifnot(prefix == get('prefix', envir=Space.pca))
  }
  subj_pca_files = get('data_files', envir=Space.pca)
  mask = get('mask', envir=Space.pca)
  groupPCA = get('pca', envir=Space.pca)
  groupWhitenedPCs = get('whitenedPCs', envir=Space.pca)
  K2 = get('K2', envir=Space.pca)
  
  ### Concat subspaces for multi-model order ICA ###
  model_orders = as.integer(model_orders)
  if (!all(is.na(model_orders))){ 
    G = NULL
    L = NULL
    if (return.PCs){  X = NULL  }else{  X = NA  }
    for (k in model_orders){
      G = cbind(G, groupPCA$U[,1:k], deparse.level=0)
      L = c(L, groupPCA$Lambda[1:k])
      if (return.PCs){
        X = cbind(X, groupPCA$X[,1:k], deparse.level=0)
      }
    }
  }else{
    G = groupPCA$U
    L = groupPCA$Lambda
    if (return.PCs){  X = groupPCA$X  }else{  X = NA  }
  }
  
  ##################################################################
  return(list('X' = X,          # group-level Principal Components, in cols.
              'G' = G,          # matrix of eigenvectors for group PCA subspace(s)
              'L' = L,          # group PCA eigenvalues for subspace(s)
              'whitenedPCs' = groupWhitenedPCs, # if Princ.Comps. scaled to eigenvals.
              'num_comps' = K2, # overall dimension of group PCA subspace
              'model_orders' = model_orders,     # dims. of concatenated subspaces
              'subj_pca_files' = subj_pca_files, # subject PCA files used in group PCA
              'prefix'         = prefix,         # file prefix for experiment
              'mask'           = mask))          # path to mask vol.
} ##################################################################
