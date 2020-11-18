# Principal Components Analysis functions
#
# Illustration of methods mentioned in:
#
#   Rachakonda S., Silva R.F., Liu J., and Calhoun V.D. (2016).
#     Memory Efficient PCA Methods for Large Group ICA.
#     Fontiers in Neuroscience 10(17): 1-15.
#
#
# Shared Notation:
#   M : number of subjs.
#   v : number of voxels
#   t : number of time points
#   p : number of subj.-level top components, after subj.-level PCA
#   k : number of group-level top components, after group-level PCA
#   Zi : data for subj. i,                [v by t] matrix
#   Yi : whitened data for subj. i,       [v by p] matrix
#   Y : temp.-concatenated group data,    [v by Mp] matrix
#   X : group-level PCA space,            [v by k] matrix
#

library(doParallel)


#----------------------------- Main fn. ----------------------------------------------------
#   Chooses among featured PCA algorithms, based on available memory vs. size of dataset,
#     runs PCA algorithm,
#       returning PC space (group or otherwise), eigenvectors (if indicated), eigenvalues
#

#############################################
PCA.Rachakonda2016 <- function(Y, k=NA, 
                               return.eigenvectors=T, 
                               algorithm=NA, unstacked=T,
                               verbose=F){
  #############################################
  #
  #     Principal Components Analysis function
  #
  # Illustration of methods mentioned in:
  #
  #   Rachakonda S., Silva R.F., Liu J., and Calhoun V.D. (2016).
  #     Memory Efficient PCA Methods for Large Group ICA.
  #     Fontiers in Neuroscience 10(17): 1-15.
  #
  # Inputs:
  #   Y : temp.-concatenated group data,    [v by Mp] matrix
  #   k : number of group-level top components, after group-level PCA
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #   algorithm           : algorithm PCA algorithm to call, fn. name as a string
  #                           see below for named options
  #   unstacked           : with input algorithm above, if True run as unstacked (low mem., parallelizable),
  #                                                   else if False run as stacked (high mem., fast)
  #
  # Outputs:
  #   X : group-level PCA space,            [v by k] matrix
  #   F : group-level eigenvector basis     [v by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # Dependencies:
  #    (pretty much every function included below)
  #

  
  
  ### PCA algorithms available in "PCA_fns.R" ###
  algorithm.options = c("PCA.EVDtime", "PCA.EVDvoxel", 
                        "PCA.SubsampledTime", "PCA.SubsampledVoxel", 
                        "PCA.MultiPowerIteration")
  
  ### Input/default handling ###
  if (is.character(algorithm) 
      && (algorithm %in% algorithm.options)
      && is.logical(unstacked)){
    if (is.na(unstacked)){  unstacked = T  }
    unstacked = unstacked
    if (verbose){ cat(paste0('\nSelecting input algorithm: ', algorithm))}
    
  }else{
    if (!all(is.na(algorithm))){
      cat(paste0('\nCould not understand requested algorithm: ',algorithm))
    }
    recommendations = find_PCAtype_best(Y, k=k, verbose=verbose)
    algorithm = recommendations$PCAtype
    unstacked = recommendations$unstacked
  }
  
  ##############################################################
  if (algorithm     == "PCA.EVDtime"){
    if (verbose){cat("\n...running PCA.EVDtime...")}
    
    PCA = PCA.EVDtime(Y, k=k, 
                      verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = NA
    
    ##############################################################
  }else if (algorithm     == "PCA.EVDvoxel"){
    if (verbose){cat("\n...running PCA.EVDvoxel...")}
    
    PCA = PCA.EVDvoxel(Y, k=k, 
                       unstacked=unstacked, 
                       return.eigenvectors=return.eigenvectors,
                       verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
    ##############################################################
  }else if (algorithm     == "PCA.SubsampledTime"){
    if (verbose){cat("\n...running PCA.SubsampledTime...")}
    
    PCA = PCA.SubsampledTime(Y, k=k, 
                             unstacked=unstacked, 
                             return.eigenvectors=return.eigenvectors,
                             verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
    ##############################################################  
  }else if (algorithm     == "PCA.SubsampledVoxel"){
    if (verbose){cat("\n...running PCA.SubsampledVoxel...")}
    
    Y.dims = get_data_dims(Y)$Y.dims
    PCA = PCA.SubsampledVoxel(Y, k=k, Y.dims=Y.dims, 
                              unstacked=unstacked, 
                              return.eigenvectors=return.eigenvectors,
                              verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
    ##############################################################
  }else if (algorithm     == "PCA.MultiPowerIteration"){
    if (verbose){cat("\n...running PCA.MultiPowerIteration...")}
    X0 = 'STP'  # intialize w/ Subsampled Time PCA, to speed convergence
    
    PCA = PCA.MultiPowerIteration(Y, k=k, 
                                  X0=X0,
                                  l=5, delta=0.00001, 
                                  unstacked=unstacked,
                                  return.eigenvectors=return.eigenvectors,
                                  verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
  }##############################################################
  
  
  #############################################
  return(PCA)
} #############################################







#----------------------------- Data Munging fns. --------------------------------------------------------------------------
#############################################
load.dataFiles <- function(datafiles, var='Y',
                           y1.inds=NA, y2.inds=NA, y3.inds=NA,
                           verbose=F){
  #############################################
  # Loads saved var. from datafiles, subj. by subj.,
  #   returning slices from y*.inds
  #
  
  library(abind) #asub
  
  if (is.list(datafiles) || is.array(datafiles)){
    return(datafiles)
  }else{
    datafiles = as.character(datafiles)
    
  }
  if (all(dir.exists(datafiles))){
    datafiles = list.files(datafiles, full.names=T)
  }
  stopifnot(all(file.exists(datafiles)))
  if (!is.list(y1.inds)){y1.inds = list(y1.inds)}
  if (!is.list(y2.inds)){y2.inds = list(y2.inds)}
  if (!is.list(y3.inds)){y3.inds = list(y3.inds)}
  
  L = max(length(y1.inds), length(y2.inds), length(y3.inds))        # counter for total number of sets of y*.inds
  if (!all(is.na(y1.inds))){   stopifnot(length(y1.inds) == L)   }  # Ensure same number of sets in each y*.ind
  if (!all(is.na(y2.inds))){   stopifnot(length(y2.inds) == L)   }
  if (!all(is.na(y3.inds))){   stopifnot(length(y3.inds) == L)   }
  
  Y = list()
  for (f in 1:length(datafiles)){
    if (verbose){cat(paste0('\nLoading: ', datafiles[f]))}
    sub = new.env()
    load(datafiles[f], envir=sub)
    stopifnot(exists(var,envir=sub))
    Y_sub = get(var,envir=sub)
    if (verbose){cat('\n')}
    if (verbose){cat(paste0(str(Y_sub)))}
    
    Yi = rep(list(Y_sub), L)
    
    if (!all(is.na(y1.inds))){
      for (l in 1:L){  Yi[[l]] = asub(Yi[[l]], y1.inds[[l]], 1)  }
    }
    if (!all(is.na(y2.inds))){
      for (l in 1:L){  Yi[[l]] = asub(Yi[[l]], y2.inds[[l]], 2)  }
    }
    if (!all(is.na(y3.inds))){
      for (l in 1:L){  Yi[[l]] = asub(Yi[[l]], y3.inds[[l]], 3)  }
    }
    
    Y = c(Y, Yi)
  }
  
  if ((length(datafiles) > 1) && (L > 1)){
    Y_sorted = list()
    for (f in 1:length(datafiles)){
      Y_sorted[[f]] = list()
      for (l in 1:l){
        ind = (f - 1) * L + l
        Y_sorted[[f]][[l]] = Y[[ind]]
      }
    }
  }else{
    Y_sorted = Y # can potentially be either list containing data from a single subj., 
    #                a list of data subsets from a single subj.,
    #                or a list of subjs. w/ 1 dataset for each
  }              # if multiple subjs. & data subsets are requested,
  #                  will be a list of lists,
  #                  with the outer 1st list indexing subjs.
  
  #############################################
  return(Y_sorted)
} #############################################

#############################################
check_stacked <- function(Y){
  #############################################
  # Checks whether or not data is already stacked
  #
  
  if (is.character(Y)){
    stacked = F
  }else if (is.list(Y) && (length(Y)>0) && is.numeric(Y[[1]])){
    stacked = T
  }else if (is.array(Y)){
    stacked = T
  }
  return(stacked)
} #############################################


# 8/15/2020 --kw-- 3d-application of this fn. needs to be debugged & tested
#
#############################################
get_data_dims <- function(Y, unstacked=T, return.data=F, verbose=F){
  #############################################
  # Find dimensions of data, loading 1st subj. if needed
  #
  # OUTPUT:
  #   v : length of 1st dim. ~ number of voxels
  #   M : number of subjs. (2nd dim. is M*p)
  #   p : number of time points/ICs (~2nd dim.)
  #   datafiles : paths to saved files
  #   Y.dims : spatial of data
  #   Y : (optional) if indicated by return.data=T,
  #         loaded data for all subjs.
  #
  #
  
  ### Required fns. ###
  get.spacetime <- function(point){
    v = p = Y.dims = NA
    if (is.matrix(point)){
      Y.dims = NA
      v = nrow(point)
      p = ncol(point)
    }else if (is.array(point)){
      Y.dims.all = dims(point)
      Y.dims = Y.dims.all[-length(Y.dims.all)] # Assume leading dimensions are all spatial...
      v = prod(Y.dims)                         # ...number of voxels is their product...
      p = Y.dims.all[length(Y.dims.all)]       # ...& assume time/component is the last dimension
    }else if (is.vector(point)){
      Y.dims = NA
      v = length(point)                # assume leading dimensions are spatial
      p = 1                             # ...assume time/component is the last dimension
    }
    return(list('v'=v, 'p'=p, 'Y.dims'=Y.dims))
  }
  
  ### Inputs/defaults ###
  stacked = !as.logical(unstacked)
  Y.dims = NA
  
  if (is.character(Y) && dir.exists(Y)){
    datafiles = list.files(path=Y, full.names=T)
    if (stacked){  # load all subjs. into memory at once
      Y = load.dataFiles(datafiles, verbose=verbose)
    }else{
      Y = load.dataFiles(datafiles[1], verbose=verbose)
    }
    M = length(datafiles)
    spacetime = get.spacetime(Y[[1]])
    Y.dims = spacetime$Y.dims
    v = spacetime$v
    p = spacetime$p
  }else if (is.character(Y) && all(file.exists(Y))){
    if (stacked){  # load all subjs. into memory at once
      datafiles = paste(basename(Y), collapse=', ')
      if (verbose){cat(paste0('Loading & stacking: ', datafiles), '...\n')}
      Y = load.dataFiles(Y, verbose=verbose)
      M = length(datafiles)
      spacetime = get.spacetime(Y[[1]])
      Y.dims = spacetime$Y.dims
      v = spacetime$v
      p = spacetime$p
    }else{            # load subjs. one-by-one
      datafiles = Y
      M = length(datafiles)
      Y.dims = dim(load.dataFiles(datafiles[1])[[1]])
      v = Y.dims[1]
      p = Y.dims[2]
    }
  }else{
    if (is.character(Y) && !all(file.exists(Y))){
      cat(paste0('WARNING: could not find file as input:'))
      print(Y)
    }else{
      if (is.list(Y)){
        M = length(Y)
        spacetime = get.spacetime(Y[[1]])
        Y.dims = spacetime$Y.dims
        v = spacetime$v
        p = spacetime$p
      }else if (is.array(Y)){
        M = dim(Y)[1]
        v = dim(Y)[2]
        p = dim(Y)[3]
      }else{  stopifnot(is.list(Y) || is.array(Y))  }
    }
    datafiles = NA
  }
  if (!return.data){  Y = NA  }
  
  #############################################
  return(list('v'=v, 'M'=M, 'p'=p, 
              'Y.dims'=Y.dims,
              'datafiles'=datafiles,
              'Y'=Y))
} #############################################

#############################################
concatTemporal_Yi <- function(Y){
  #############################################
  # Temporally Concatenates subj.-level data Y
  #   where Y is either list of Yi   [v by p] matrices
  #   or Y is array with dim.   [M by v by p]
  #
  
  if (is.matrix(Y)){
    return(Y)  # PASS, nothing to do
  }else if (is.list(Y)){
    # iterates over elements of Y & binds results by column
    return(do.call(cbind, Y))
  }else if (is.array(Y)){
    # re-arranges Y, putting subj. dimension last, iterates over ~voxel dim., binding results by column
    return(t(apply(aperm(Y,c(2,3,1)),1,cbind)))
  }
}

#############################################
zeroMean_Yi <- function(Y){
  #############################################
  # Removes column mean of subj.-level data Y
  #     ~centering every vol. at every timepoint
  #   where Y is either list of Yi   [v by p] matrices
  #   or Y is array with dim.   [M by v by p]
  #
  
  if (is.matrix(Y)){
    Y = sweep(Y, 2, colMeans(Y))
  }else if (is.list(Y)){
    for (i in 1:length(Y)){
      Y[[i]] = sweep(Y[[i]], 2, colMeans(Y[[i]]))
    }
  }else if (is.array(Y)){
    for (i in 1:dim(Y)[1]){
      Y[i,,] = sweep(Y[i,,], 2, colMeans(Y[i,,]))
    }
  }
  return(Y)
}

#############################################
spatialSubsample <- function(dims, depth=2, depth.l=1, verbose=F){
  #############################################
  # Creates indices to subsample data,
  #   taking into spatial structure
  #   using R's column major ordering of indexing
  #     (i.e., rows are filled in 1st,
  #             then columns, then 3rd dim, etc.)
  #
  # INPUT:
  #   dims    : dimensions of data array
  #   depth   : depth of sampling, currently only able to handle depth=2
  #   depth.l : level of the depth to subsample, {1:depth}
  # OUTPUT:
  #   -a list of length depth,
  #       each element a linear T/F indices specifying selected subsamples
  #
  
  stopifnot(is.vector(dims))
  stopifnot(depth.l <= depth)
  
  V = prod(dims)
  D = length(dims)
  C = depth
  C.i = c(1:C)
  sample.ind.mat = matrix(NA,V,D)
  
  #initialize 1st dim. by subdividing its indices
  inds = rep(rep(C.i, times=ceiling(dims[1] / length(C.i)), length.out=dims[1])) #subdivide dim. indices
  for (d in 1:D){
    
    sample.ind.mat[,d] = rep(inds, times=prod(dims) / length(inds)) #...fill col. w/ repeating inds

    if (verbose){
      print(inds)
      print(sample.ind.mat)
    }
    if (d+1<=D){
      inds = rep(rep(C.i, times=ceiling(dims[d+1] / length(C.i)), length.out=dims[d+1])) #subdivide next dim. indices
      inds = rep(inds, each=prod(dims[1:d])) # set up inds for next dimension...
    }
  }
  
  subsample.l = apply(sample.ind.mat==depth.l,1,prod)
  if (verbose){   print(subsample.l)  }
  
  return(as.logical(subsample.l))
}

#############################################
get_eigenValsVecs_fromPCs <- function(Y, X,
                                      unstacked=F, verbose=F){
  #############################################
  # Find eigenvalues & eigenvectors from data Y,
  #   given input Principal Components X
  #
  # Theory: 
  #   If X is the [voxels x k] group PCs, solve def. of PCs in 1st eq. 
  #   for eigenvectors, the [Mp x k] F_, w/o needing generalized inverse:
  #
  #            X =          Y %*% F_
  #   t(Y) %*% X = (t(Y) %*% Y) %*% F_
  #   t(Y) %*% X = F_ %*% Lambda %*% t(F_) %*% F_
  #   t(Y) %*% X = F_ %*% Lambda = F_unscaled
  #
  #     Lambda^2 = Lambda %*% t(F_) %*% F_ %*% Lambda
  #     Lambda^2 = t(F_unscaled) %*% F_unscaled
  #
  #           F_ = F_unscaled %*% Lambda^-1
  #
  
  stacked = !unstacked
  
  data.dims = get_data_dims(Y, T, 
                            return.data=stacked)
  if (verbose){print(str(data.dims))}
  datafiles = data.dims$datafiles
  if (stacked){ Y = data.dims$Y } # loaded data
  M = data.dims$M # number of subjs.
  v = data.dims$v # number of voxels
  p = data.dims$p # dimension of voxels (i.e., time, reduced-PCA T1, etc.)
  
  if (unstacked){
    
    F_unscaled = foreach(sub=iter(datafiles),
                         .export=c('load.dataFiles','zeroMean_Yi'),
                         .combine = rbind) %dopar% {
                           
                           Yi = load.dataFiles(sub, verbose=verbose)[[1]]
                           Yi = zeroMean_Yi(Yi)
                           return( t(Yi) %*% X )
                           
                         }
  }else{
    Y = load.dataFiles(Y)
    Y = zeroMean_Yi(Y)
    Y = concatTemporal_Yi(Y)
    F_unscaled = t(Y) %*% X
    
  }
  
  Lambda = sqrt(diag(t(F_unscaled) %*% F_unscaled))
  F_ = F_unscaled %*% diag(1 / Lambda)
  F_sd = apply(F_,2, sd)   # if needed for normalization
  
  eigs = list('F_' = F_,
              'Lambda' = Lambda,
              'F_sd' = F_sd)
  return(eigs)
}

# #############################################
# get_eigenvalues_fromPCs <- function(Y, X, Lambda=NA,
#                                     unstacked=F, verbose=F){
#   #############################################
#   # Find eigenvectors from data Y,
#   #   given input Principal Components X
#   #
#   # Theory: 
#   #   If X is the [voxels x k] group PCs, solve def. of PCs in 1st eq. 
#   #   for eigenvectors, the [Mp x k] F_, w/o needing generalized inverse:
#   #
#   #            X =          Y %*% F_
#   #   t(Y) %*% X = (t(Y) %*% Y) %*% F_
#   #   t(Y) %*% X = F_ %*% Lambda %*% t(F_) %*% F_
#   #   t(Y) %*% X = F_ %*% Lambda
#   #           F_ = t(Y) %*% X %*% Lambda^-1
#   #
#   
#   stacked = !unstacked
#   
#   data.dims = get_data_dims(Y, T, 
#                             return.data=stacked, verbose=verbose)
#   if (verbose){print(str(data.dims))}
#   datafiles = data.dims$datafiles
#   if (stacked){ Y = data.dims$Y } # loaded data
#   Y = zeroMean_Yi(Y)
#   Y = concatTemporal_Yi(Y)
#   M = data.dims$M # number of subjs.
#   v = data.dims$v # number of voxels
#   p = data.dims$p # dimension of voxels (i.e., time, reduced-PCA T1, etc.)
#   
#   if (unstacked){
#     F_ = foreach(sub=iter(datafiles),
#                  .export=c('load.dataFiles','zeroMean_Yi'),
#                  .combine = rbind) %dopar% {
#                    Yi = load.dataFiles(sub, verbose=verbose)[[1]]
#                    Yi = zeroMean_Yi(Yi)
#                    return( t(Yi) %*% X )
#                  }
#   }else{
#     F_ = t(Y) %*% X
#   }
#   
#   if (!all(is.na(Lambda))){
#     if (isSymmetric.matrix(Lambda)){
#       F_ = F_ %*% solve(Lambda)
#     }else if (is.vector(Lambda)){
#       F_ = F_ %*% diag( 1 / Lambda )
#     }
#   }
#   
#   F_ = apply(F_ ,2, function(yy) return( yy / sd(yy))) # normalize to sd
#   return(F_)
# } 


#------------------------------------ Memory-calculations & algorithm recommendations ----------------------------

#############################################
find_PCAtype_best <- function(data, k=NA, verbose=F){
  #############################################
  # Finds best PCA algorithm,
  #   based on recommendations in Figure 8, Rachakonda et al. (2016)
  #
  
  ### Get data dims ###
  data.dims = get_data_dims(data)
  M = data.dims$M
  v = data.dims$v
  p = data.dims$p
  k = min(v, M*p, k, na.rm=T)
  
  ### Get available RAM ###
  memFree = get_RAMavailable()
  
  ### Decision table in Figure 8 ###
  if (min(v, M*p) <= 10000){
    if (v < M*p){
      PCAtype = 'PCA.EVDvoxel'
      unstacked = T
    }else{
      PCAtype = 'PCA.EVDtime'
      unstacked = F
    }
  }else{
    if (v*M*p > memFree){
      PCAtype = 'PCA.MultiPowerIteration'
      unstacked = T
    }else{
      PCAtype = 'PCA.MultiPowerIteration'
      unstacked = F
    }
  }
  if (verbose){
    cat(paste0('\nRecommended PCA algorithm:  ', PCAtype))
    if (unstacked){cat(' (unstacked)')}else{cat(' (stacked)')}
  }
  
  ### Sanity Check ###
  memInfo = calc_PCA_memRequirements(data, k, pca_type=PCAtype)
  if (length(memInfo) > 1){
    memNeeded = min(memInfo[[1]]$mem, memInfo[[2]]$mem, na.rm=T)
  }else{
    memNeeded = memInfo[[1]]$mem
  }
  
  if (verbose){
    cat(paste0('\n...requiring  ',signif(memNeeded/1024),' Mb  of Memory'))
  }
  if (memFree < memNeeded){
    print(paste0("WARNING:  Memory available (",signif(memFree/1024)," Mb) less than needed for PCA (",signif(memNeeded/1024)," Mb)!"))
  }
  
  #############################################
  return(list('PCAtype' = PCAtype,
              'unstacked' = unstacked))
} #############################################

#############################################
get_RAMavailable <- function(verbose=F){
  #############################################
  # Finds currently available RAM in computer
  #
  
  if (Sys.info()['sysname'] == 'Linux'){
    gc()
    memfree = as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=T)) # returns RAM free, in kB
  }else if (Sys.info()['sysname'] == 'Windows'){
    print("Guessing mem. size for windows...")
    memfree = memory.size() * 1024  # get limits on mem. available for R, in Kb
  }else{
    print("ERROR:  get_RAMavailable() has not been setup for non-linux operating systems currently!")
    memfree = NA
  }
  if (verbose){  cat(paste0('\nCurrently available Memory:  ', memfree))}
  return(memfree)
}

#############################################
calc_PCA_memRequirements <- function(data, k, pca_type='all', verbose=F){
  #############################################
  # Calculate PCA memory requirements,
  #   for various PCA algorithms in 'PCA_fns.R'
  # Based off appendix in Rachakonda et al. (2016). Memory Efficient PCA Methods for Large Group PCA.
  #                           Frontiers in Neuroscience (10)17: p. 1-17.
  #
  # INPUT:
  #   data     : input data
  #   k        : final number of PCs
  #   pca_type : fn. in 'PCA_fns.R' or 'all'
  # OUTPUT: 
  #     list w/ elements indexed by fn. name, containing required RAM in Kb
  # REQUIRES:
  #   get_data_dims() from 'PCA_fns.R'
  #
  
  if (verbose){
    cat(paste0('\nCalculating memory requirements for PCA algorithms...'))
  }
  
  ### Input handling & defaults ###
  if (pca_type == 'all'){
    pca_types = c('PCA.EVDtime', 'PCA.EVDvoxel',
                  'PCA.SubsampledTime', 'PCA.SubsampledVoxel',
                  'PCA.MultiPowerIteration')
  }else{  pca_types = pca_type  }
  bytes = 8 # double precision, numeric data
  if (!exists('l')){  l = 5   }  # default multiplier for MPOWIT
  if (!exists('g')){  g = 20  }  # default number of subjs. / group in STP
  if (!exists('k_prime')){
    k_prime = max(k*5, 500) # dim. of PCA during group update, defaults suggested number in paper
  }
  
  ### Get data dims ###
  data.dims = get_data_dims(data)
  M = data.dims$M
  v = data.dims$v
  p = data.dims$p
  if (!exists('v_prime_a')){  v_prime_a = v / 2  } # all odd voxels
  if (!exists('v_prime_b')){  v_prime_b = v / 2  } # all even voxels
  
  ### Main fn. ###
  memInfo = list()
  m = 0
  for (pcaType in pca_types){
    if (pcaType == 'PCA.EVDtime'){
      mem = v * M * p       # total data size stacked...
      + (M * p)^2           # ...size of [time x time] cov. matrix
      + v * k               # ...size of group PCA space X
      unstacked = F
    }else if (pcaType == 'PCA.EVDvoxel'){
      mem = v * p           # total data size for single subjs. unstacked...
      + v^2                 # ...size of [time x time] cov. matrix
      + v * k               # ...size of group PCA space X
      unstacked = T
    }else if (pcaType == 'PCA.SubsampledTime'){
      g1 = min(g, M)
      mem = v * g1 * p      # total data size for group...
      + (g1 * p)^2          # ...concat. group cov. matrix
      + 2 * v * k_prime     # ...size of group PCA spaces X_g & X_g+1 in eq. (13)
      + g1 * p * k_prime    # ...size of final group PCA space in eq. (15)
      + 2 * k_prime^2       # ...size of Lambda in eq. (14), as block diag. of Lambda_g & Lambda_g+1, for final group PCA est.
      unstacked = T
    }else if (pcaType == 'PCA.SubsampledVoxel'){
      mem = v_prime_a^2 + v_prime_b^2   # Voxel cov. matrices...
      + 2 *v * M * p * k_prime          # ...intermediary projected group PCA spaces
      + (2 * k_prime)^2                 # ...concat. subsampled voxel cov. matrices
      + (2 * k_prime)^2                 # ...size of Lambda in eq. (10)
      unstacked = T
    }else if (pcaType == 'PCA.MultiPowerIteration'){
      Ma = 2 * v * l * k    # X_j & X_j-1, for j=1...
      + (l * k)^2           # ...[Chi_j^T x Chi_j] cov. matrix
      Mb = v * l * k        # ...X_j
      + 3 * (l * k)^2       # ...[Chi_j^T x Chi_j] & [Chi_j-1^T x Chi_j-1] during ortho() & Lambda_j-1 
      mem_stacked = v * M * p       # total data size stacked...
      + max(Ma, Mb)                 # ...largest size of iteration
      mem_unstacked = v * p         # size of single subjs.' data...
      + max(Ma, Mb)                 # ...largest size of iteration
      mem = c(mem_stacked,
              mem_unstacked)
      unstacked = c(F, T)
    }
    mem = mem * bytes / 1024        # return size in Kb
    if (length(mem) == 1){
      m = m + 1
      memInfo[[m]] = list('type' = pcaType,
                          'mem' = mem,
                          'units' = 'Kb',
                          'unstacked' = unstacked)
    }else{
      m = m + 1
      memInfo[[m]] = list('type' = pcaType,
                          'mem' = mem[1],
                          'units' = 'Kb',
                          'unstacked' = unstacked[1])
      m = m + 1
      memInfo[[m]] = list('type' = pcaType,
                          'mem' = mem[2],
                          'units' = 'Kb',
                          'unstacked' = unstacked[2])
    }
    if (verbose){
      cat(paste0('\n...PCA using  ', pcaType))
      if (unstacked[1]){cat(' (unstacked)')}else{cat(' (stacked)')}
      cat(paste0('\n      Memory  =  ',mem[1]/1024,' Mb'))
      if ((length(mem) > 1) && (length(unstacked) > 1)){
        cat(paste0('\n...PCA using  ', pcaType))
        if (unstacked[2]){cat(' (unstacked)')}else{cat(' (stacked)')}
        cat(paste0('\n      Memory  =  ',mem[2]/1024,' Mb'))
      }
    }
  }
  
  #############################################
  return (memInfo)
} #############################################






#----------------------------------- PCA  algorithms  -----------------------------------------------------------
#############################################
PCA.EVDtime <- function(Y, k=NA, verbose=F){
  #############################################
  # Principal Components Analysis via Eigenvalue Decomposition (well, actually ~SVD)
  #   to find group-level PCA space
  #
  # INPUT:
  #   Y : subj.-level data,
  #         as temp.-concatenated data,    [v by Mp] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  # OUTPUT:
  #   X : group-level PCA space,            [v by k] matrix
  #   F : group-level eigenvector basis     [v by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  # REQUIRES:
  #   concatTemporal_Yi() : data munging for arrays & lists
  # EQUATIONS & NOTATION:
  #   C : covariance matrix, where          [Mp by Mp] matrix
  #                   C = (t(Y) %*% Y) / (v-1)                eq. (1)
  #                   C = F %*% Lambda %*% t(F)               eq. (2)       
  #   F_      : group-level rotation matrix, resulting from group-level PCA, [Mp by k]
  #   \Lambda : group-level diagonal matrix of eigenvalues
  #                   X = Y %*% F %*% solve(sqrt(Lambda))     eq. (3)
  #
  
  if (verbose){cat("\nCalculating PCA space using EVD on time*subjs. dimension...")}

  if (is.character(Y) && all(file.exists(Y))){
    Y = load.dataFiles(Y) #loads all subjs. into mem. at once
  }
  Y = zeroMean_Yi(Y)
  Y = concatTemporal_Yi(Y)
  v = nrow(Y)
  Mp = ncol(Y)
  k = min(k, v, Mp, na.rm=T)

  C = (t(Y) %*% Y) / (v-1)  # eq. (1)
  UDV = svd(C, k)
  
  X = Y %*% UDV$u %*% diag( 1 / sqrt(UDV$d[1:k])) # eq. (3)
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'F_'     = UDV$u,       # Eigenvectors
              'Lambda' = UDV$d[1:k])) # Eigenvalues
} #############################################

#############################################
PCA.EVDvoxel <- function(Y, k=NA, 
                         unstacked=F, return.eigenvectors=F,
                         verbose=F){
  #############################################
  # Principal Components Analysis via Eigenvalue Decomposition (well, actually ~SVD)
  #   on voxel by voxel covariance matrix
  #     to find group-level PCA space
  #
  # INPUT:
  #   Y : subj.-level data,
  #         as vector of paths to saved subj. data,
  #         as temp.-concatenated data,    [v by Mp] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #                       (NOTE: requires parallelization cluster set up)
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # OUTPUT:
  #   X   : group-level PCA space,          [v by k] matrix
  #   Chi : group-level eigenvector basis   [v by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # EQUATIONS:
  #   Cv_i : subj.-level voxel-wise cov. matrix, [v by v] matrix
  #   Cv : group-level voxel-wise cov. matrix,   [v by v] matrix
  #                   Cv = (Y %*% t(Y)) / (v-1)
  #                   Cv = sum(Cv_i) 
  #                   Cv = \Chi %*% \Lambda %*% t(\Chi)
  #                   X = \Chi * sqrt(v-1)
  #
  # REQUIRES:
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.dataFiles()    :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #
  # SPEED:   time complexity  O(Mpv^2)
  #
  
  if (verbose){cat("\nCalculating PCA space using EVD on voxel dimension...")}
  
  ### Defaults & Inputs ###
  stacked = check_stacked(Y) || !as.logical(unstacked)
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  
  data.dims = get_data_dims(Y)
  M = data.dims$M
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  stopifnot(exists('M')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = min(v, M*p, k, na.rm=T)
  
  Cv = matrix(0,v,v)
  
  if (unstacked){
    Cv = foreach(sub=iter(datafiles), 
                 .export=c('load.dataFiles',
                           'zeroMean_Yi'), 
                 .combine='+') %dopar% {
                   
            Yi = load.dataFiles(sub, verbose=verbose)[[1]]
            Yi = zeroMean_Yi(Yi)
            return((Yi %*% t(Yi)) / (v-1))
            
          }
  }else{
    Y = load.dataFiles(Y)
    Y = zeroMean_Yi(Y)
    Y = concatTemporal_Yi(Y)
    Cv = (Y %*% t(Y)) / (v-1)
  }
  
  UDV = svd(Cv, k)
  X = UDV$u * sqrt(v-1)
  
  ### Eigenvectors & Eigenvalues ###
  Lambda = UDV$d[1:k]
  if (return.eigenvectors){
    if (all(is.na(datafiles))){  datafiles = Y  }
    F_ = get_eigenValsVecs_fromPCs(datafiles, X, unstacked=unstacked)$F_
  }else{  'F_' = NA  }
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'F_'     = F_,          # Eigenvectors
              'Lambda' = Lambda))     # Eigenvalues
} #############################################

# #############################################
# PCA.EVDvoxel <- function(Y, k=NA, 
#                          unstacked=F, return.eigenvectors=F,
#                          verbose=F){
#   #############################################
#   # Principal Components Analysis via Eigenvalue Decomposition (well, actually ~SVD)
#   #   on voxel by voxel covariance matrix
#   #     to find group-level PCA space
#   #
#   # INPUT:
#   #   Y : subj.-level data,
#   #         as vector of paths to saved subj. data,
#   #         as temp.-concatenated data,    [v by Mp] matrix
#   #         as list of [v by p] matrices
#   #         or as array w/ subjs. indices in 1st dim.
#   #   k : number of group-level top components, after group-level PCA
#   #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
#   #                 in order to minimize memory use, despite very large number of subjs.
#   #                       (NOTE: requires parallelization cluster set up)
#   #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
#   #
#   # OUTPUT:
#   #   X   : group-level PCA space,          [v by k] matrix
#   #   Chi : group-level eigenvector basis   [v by k] matrix
#   #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
#   #
#   # EQUATIONS:
#   #   Cv_i : subj.-level voxel-wise cov. matrix, [v by v] matrix
#   #   Cv : group-level voxel-wise cov. matrix,   [v by v] matrix
#   #                   Cv = (Y %*% t(Y)) / (v-1)
#   #                   Cv = sum(Cv_i) 
#   #                   Cv = \Chi %*% \Lambda %*% t(\Chi)
#   #                   X = \Chi * sqrt(v-1)
#   #
#   # REQUIRES:
#   #   check_stacked()     : data munging
#   #   get_data_dims()     :  "     "
#   #   load.dataFiles()    :  "     "
#   #   zeroMean_Yi()       :  "     "
#   #   concatTemporal_Yi() :  "     "
#   #
#   # SPEED:   time complexity  O(Mpv^2)
#   #
#   
#   if (verbose){cat("\nCalculating PCA space using EVD on voxel dimension...")}
#   
#   ### Defaults & Inputs ###
#   stacked = check_stacked(Y) || !as.logical(unstacked)
#   unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
#   if (unstacked){  library(foreach); library(iterators)  }
#   
#   data.dims = get_data_dims(Y, unstacked, return.data=stacked)
#   M = data.dims$M
#   v = data.dims$v
#   p = data.dims$p
#   datafiles = data.dims$datafiles
#   if (stacked){  Yi = data.dims$Y  }
#   stopifnot(exists('M')); stopifnot(exists('v')); stopifnot(exists('p'))
#   k = min(v, M*p, k, na.rm=T)
#   
#   Cv = matrix(0,v,v)
#   if (unstacked){
#     Cv = foreach(sub=iter(datafiles), .export=c('load.dataFiles',
#                                                 'zeroMean_Yi'), .combine='+') %dopar% {
#       Yi = load.dataFiles(sub, verbose=verbose)[[1]]
#       Yi = zeroMean_Yi(Yi)
#       return((Yi %*% t(Yi)) / (v-1))
#     }
#   }else{
#     stopifnot(exists('Yi'))
#     Yi = zeroMean_Yi(Yi)
#     Yi = concatTemporal_Yi(Yi)
#     Cv = (Yi %*% t(Yi)) / (v-1)
#   }
#   
#   UDV = svd(Cv, k)
#   X = UDV$u * sqrt(v-1)
#   
#   ### Eigenvectors & Eigenvalues ###
#   if (return.eigenvectors){
#     
#     F_ = get_eigenvalues_fromPCs(Y, X, unstacked=unstacked)
#     
#   }else{  F_ = NA  }
#   
#   #############################################
#   return(list('X'      = X,           # Principal Comps.
#               'F_'     = F_,          # Eigenvectors
#               'Lambda' = UDV$d[1:k])) # Eigenvalues
# } #############################################


# 8/15/2020 --kw-- added ability to get eigenvalues by backprojecting PCs onto data
# 8/1/2020 --kw-- Added parallelization, based on example fn. in GIFT toolbox
#                   Good agreement for stacked vs. (parallel) unstacked vs. PCA.EVDtime()
#############################################
PCA.SubsampledTime <- function(Y, k=NA, g=20, parFactor=4,
                               unstacked=F, return.eigenvectors=F,
                               verbose=F){
  #############################################
  # Principal Components Analysis on sub-sample of time dim.,
  #   i.e., w/ temp. concat. PCA, randomly sample subjs. on time dim.
  #     to find group-level PCA space
  #
  # INPUT:
  #   Y : subj.-level data,
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   g : number of subjs. per group, group size
  #   parFactor : parallelization factor, number of worker nodes (only relevant w/ unstacked data)
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # OUTPUT:
  #   X : group-level PCA space,            [v by k] matrix
  #   F_ : group-level eigenvector basis    [Mp by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # REQUIRES:
  #   PCA.EVDtime()
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.dataFiles()    :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #  
  # SPEED:   time complexity  O(v(Mp)^2)
  #
  
  if (verbose){cat("\nCalculating PCA space by subsampling time*subjs. dimension...")}
  
  ### Required fns. ###
  partition_group <- function(subj.count, group.size){ 
    # randomly partitions subjs. into groups with set group size
    s.inds = c(1:subj.count)
    if (group.size == 1){ # return list w/ single subj. index per element
      partition = as.list(sample(subj.count,subj.count))
      s.grouped = unlist(partition)
      s.ungrouped = s.inds[!(s.inds %in% s.grouped)]
    }else if (group.size >= subj.count){ # return list of randomized subj. indices
      partition = list(sample(subj.count,subj.count))
      s.grouped = unlist(partition)
      s.ungrouped = s.inds[!(s.inds %in% s.grouped)]
    }else{
      partition = list()
      parts = subj.count %/% group.size
      s.grouped = c()
      s.ungrouped = s.inds
      for (p in 1:parts){
        partition[[p]] = sample(s.ungrouped, group.size)
        s.grouped = c(s.grouped, partition[[p]])
        s.ungrouped = s.ungrouped[!(s.ungrouped %in% s.grouped)]
      }
      if (length(s.ungrouped) > 0){
        if (length(s.ungrouped) == 1){ #see sample() documentation
          partition[[p+1]] = s.ungrouped 
        }else{
          partition[[p+1]] = sample(s.ungrouped, length(s.ungrouped))
        }
        s.grouped = c(s.grouped, partition[[p+1]])
        s.ungrouped = s.ungrouped[!(s.ungrouped %in% s.grouped)]
      }
    }
    stopifnot(all(s.inds %in% s.grouped))
    stopifnot(!any(s.inds %in% s.ungrouped))
    return(partition)
  } ###################
  calc_GroupSubSampleTimeSVD <- function(data, partition_inds, k_prime,  Xg_0=NA){
    ### Calc. for C_g ###
    if (is.character(data)){
      Yg = load.dataFiles(data[unlist(partition_inds)])
    }else if (is.list(data)){
      s.g = 0
      Yg = list()
      for (i in partition_inds){
        s.g = s.g + 1
        Yg[[s.g]] = data[[i]]
      }
    }else if (is.array(data)){
      Yg = data[partition_inds,,]
    }
    
    ### PCA ###
    Yg = zeroMean_Yi(Yg)
    Yg = concatTemporal_Yi(Yg)
    UVD.Yg = PCA.EVDtime(Yg, k_prime)                # eq. (11) & (12)
    Xg = Yg %*% UVD.Yg$F_                            # eq. (13)
    
    ### Combine Groups ###
    if (all(is.na(Xg_0))){
      return(Xg)    # return PCs, if 1st iteration
    }else{
      stopifnot(exists('Xg_0') && all(is.numeric(Xg_0)))
      
      ### Calc. Cov. of g-1 & g ###
      Cg0g1_g0g0 = (t(Xg_0) %*% Xg_0) / (v-1)        # eq. (14) left
      Cg0g1_g0g1 = (t(Xg_0) %*% Xg)   / (v-1)
      Cg0g1_g1g0 = (t(Xg) %*% Xg_0)   / (v-1)
      Cg0g1_g1g1 = (t(Xg) %*% Xg)     / (v-1)
      Cg0g1 = rbind(cbind(Cg0g1_g0g0, Cg0g1_g0g1),
                    cbind(Cg0g1_g1g0, Cg0g1_g1g1))
      
      UDV.g0g1 = svd(Cg0g1, k_prime)                 # eq. (14) right
      Xg = cbind(Xg_0, Xg) %*% UDV.g0g1$u          # eq. (15) 
    }
    return(Xg)
  } ###################
  
  
  ### Defaults & Inputs ###
  stacked = check_stacked(Y) || !as.logical(unstacked) # Check requested arg. inputs vs. reality of data input Y
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  
  data.dims = get_data_dims(Y, unstacked, return.data=stacked)
  M = data.dims$M
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  if (stacked){  Y = data.dims$Y  }
  stopifnot(exists('M')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = min(v, M*p, k, na.rm=T)
  k_prime = min(max(k*5, 500), # dim. of PCA during group update, defaults suggested number in paper
                v, M*p)        #  limited by rank of matrix
  
  
  ### Subdivide by groups & EVD on each ###
  if (unstacked){
    g_node = ceiling(M / min(M, parFactor))      # number of subjs. on each parallel node
    node.part.list = partition_group(M, g_node)
    numNodes = length(node.part.list)
    Xg = foreach(n.pt.i=icount(numNodes),
                 .export=c('load.dataFiles', 'zeroMean_Yi', 'concatTemporal_Yi',
                           'PCA.EVDtime'), 
                 .combine='cbind') %dopar% {
                   dataSubset = datafiles[unlist(node.part.list[[n.pt.i]])]
                   grp.part.list = partition_group(length(dataSubset), g)
                   Xg_0 = Xg = NA
                   for (pt.i in 1:(length(grp.part.list))){
                     Xg_0 = Xg #re-name prev. iteration
                     Xg = calc_GroupSubSampleTimeSVD(dataSubset, grp.part.list[[pt.i]], 
                                                     k_prime, Xg_0)
                   }
                   return(Xg)
                 }
    UDV.g_final = svd((t(Xg) %*% Xg) / (v-1)) # combine individual subgroup results
    Xg = Xg %*% UDV.g_final$u
  }else{
    grp.part.list = partition_group(M, g)
    Xg_0 = Xg = NA
    for (pt.i in 1:(length(grp.part.list))){
      Xg_0 = Xg  #re-name prev. iteration
      Xg = calc_GroupSubSampleTimeSVD(Y, grp.part.list[[pt.i]],
                                      k_prime, Xg_0)
    }
    UDV.g_final = svd(t(Xg) %*% Xg / (v-1))
  }
  
  ### Group PCA space ###
  X = Xg[,1:k] %*% diag( 1 / sqrt(UDV.g_final$d[1:k])) # ~final paragraph for STP section
  
  ### Eigenvectors & Eigenvalues ###
  Lambda = UDV.g_final$d[1:k]
  if (return.eigenvectors){
    if (all(is.na(datafiles))){  datafiles = Y  }
    F_ = get_eigenValsVecs_fromPCs(datafiles, X, unstacked=unstacked)$F_
  }else{  'F_' = NA  }
  
  # ### Some Interesting PCA theory ###
  # X = Xg_0[,1:k] %*% diag( 1 / apply(Xg_0[,1:k],2,sd))    # An alternat way of calc. X.  Note that the st.dev. of cols. ~ sqrt. of lambdas...
  # Lambda = diag(t(Xg_0[,1:k]) %*% Xg_0[,1:k]) / (v-1)     # ...and the scaled diagonal of the cov. are the eigenvalues!
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'F_'     = F_,          # Eigenvectors
              'Lambda' = Lambda))     # Eigenvalues
} #############################################
#   ### Eigenvectors & Eigenvalues ###
#   Lambda =  UDV.g_final$d[1:k]
#   if (return.eigenvectors){
#     
#     F_ = get_eigenvalues_fromPCs(Y, X,unstacked=unstacked)
#     
#   }else{  F_ = NA  }
#   
#   # ### Some Interesting PCA theory ###
#   # X = Xg_0[,1:k] %*% diag( 1 / apply(Xg_0[,1:k],2,sd))    # An alternat way of calc. X.  Note that the st.dev. of cols. ~ sqrt. of lambdas...
#   # Lambda = diag(t(Xg_0[,1:k]) %*% Xg_0[,1:k]) / (v-1)     # ...and the scaled diagonal of the cov. are the eigenvalues!
#   
#   
#   #############################################
#   return(list('X' = X,
#               'F_' = F_,
#               'Lambda' = UDV.g_final$d[1:k]))
# } #############################################


# 8/1/2020 --kw-- Method for subsampled voxel PCA not included in GIFT,
#                   likely due to relatively poor performance vs. subsampled time PCA
#                     (i.e., 2x data loads, less accurate, potentially larger memory footprint)
# 5/31/2020 --kw-- PCA calc. With spatially equidistant subsampling grid,
#                     with re-scaling parameters partially moved earlier, before final SVD to calc. eigenvalues
#                       no re-scaling eigenvalues by taking sqrt. of final SVD, but still re-scales group-level PCA space X .
#                     Good aggreement in both ouput Lambda's & group-level PCA matrix X compared to PCA.EVDtime().
#                     Appears to be a fundamental trade-off between accurate initial scaling of eigenvalues
#                       & accurate estimation of group PCA space X,
#                       in below, fn. follows paper for accurate X,
#                         lastly re-scales eigenvalues for comparison w/ other PCA fns.
# 6/30/2020 --kw--  Small disagreement between eigenvectors of stacked vs. unstacked versions of fn.,
#                       but all eigenvalues in agreement
#
#############################################
PCA.SubsampledVoxel <- function(Y, k=NA, Y.dims=NA, 
                                unstacked=F, return.eigenvectors=F,
                                verbose=F){
  #############################################
  # Principal Components Analysis on sub-sample of voxels
  #   to find group-level PCA space
  # Currently only supports sampling depth of 2 in single dim.,
  #   resulting in odd/even row sampling of data matrix
  #
  # INPUT:
  #   Y : subj.-level data,
  #         as temp.-concatenated data,    [v by Mp] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   Y.dims    : dims. of Y used for odd/even spatial subsampling,
  #                 defaults to every other row of Y
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #                       (NOTE: requires parallelization cluster set up)
  #  
  # OUTPUT:
  #   X : group-level PCA space,            [v by k] matrix
  #   F_ : group-level eigenvector basis    [Mp by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # REQUIRES:
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.dataFiles()    :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #   spatialSubsample()  : sampling taking into account spatial structure of data
  # EQUATIONS:
  #   complicated, see ref. (6)-(10)
  #
  
  if (verbose){cat("\nCalculating PCA space by subsampling in voxel (or spatial) dimension(s)...")}
  
  ### Defaults & Inputs ###
  stacked = check_stacked(Y) || !as.logical(unstacked)
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  
  data.dims = get_data_dims(Y, unstacked, return.data=stacked)
  M = data.dims$M
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  if (stacked){  Y = data.dims$Y  }
  stopifnot(exists('M')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = min(v, M*p, k, na.rm=T)
  
  
  ### Sub-sample dimensions ###
  if (all(is.na(Y.dims))){ Y.dims=v }else{ Y.dims = as.integer(Y.dims)  }
  odd.i = spatialSubsample(Y.dims,2,1)
  even.i = spatialSubsample(Y.dims,2,2)
  v_prime_a = sum(odd.i)
  v_prime_b = sum(even.i)
  if (verbose){
    cat('\n...Spatial sub-sampling w/ depth 2,    all odd & all even,')
    cat('\n......resulting vector indices of lengths: ')
    cat(paste0(v_prime_a,' & ',v_prime_b,'  /  ',v))
  }
  k_prime_a = min(max(k*5, 500),  # dim. of PCA during group update, defaults suggested number in paper,
                  v_prime_a, M*p)            #    but limited by rank of cov. matrix
  k_prime_b = min(max(k*5, 500),  # dim. of PCA during group update, defaults suggested number in paper,
                  v_prime_b, M*p)            #    but limited by rank of cov. matrix
  k = min(k, v, M*p, k_prime_a+k_prime_b, na.rm=T)
  
  
  ### Main fn. ###
  if (unstacked){  # Calc. cov. matrix...
    Cv_aabb = foreach(sub=iter(datafiles), 
                      .export=c('load.dataFiles', 'zeroMean_Yi'), 
                      .combine='+') %dopar% {
                        Yi = load.dataFiles(sub, y1.inds=list(odd.i, even.i), verbose=verbose)
                        Yi = zeroMean_Yi(Yi)
                        Cv_aa = Yi[[1]] %*% t(Yi[[1]]) / (v_prime_a - 1) # eq (5) unstacked
                        Cv_bb = Yi[[2]] %*% t(Yi[[2]]) / (v_prime_b - 1) # eq (5) unstacked
                        diff.i = sum(odd.i) - sum(even.i)
                        if (diff.i > 0){                                 # pad matrix in order to concatenate & format for return
                          Cv_bb = rbind(Cv_bb, matrix(NA, diff.i, ncol(Cv_bb)))
                          Cv_bb = cbind(Cv_bb, matrix(NA, nrow(Cv_bb), diff.i))
                        }else if (diff.i < 0){                           # pad matrix in order to concatenate & format for return
                          diff.i = abs(diff.i)
                          Cv_aa = rbind(Cv_aa, matrix(NA, diff.i, ncol(Cv_aa)))
                          Cv_aa = cbind(Cv_aa, matrix(NA, nrow(Cv_aa), diff.i))
                        }
                        stopifnot(dim(Cv_aa) == dim(Cv_bb))
                        return(rbind(Cv_aa, Cv_bb))
                      }
    Cv_aa = Cv_aabb[1:v_prime_a, 1:v_prime_a]
    Cv_bb = Cv_aabb[v_prime_a + 1:v_prime_b,
                    1:v_prime_b]
  }else{
    Y = zeroMean_Yi(Y)
    Y = concatTemporal_Yi(Y)
    Ya = Y[odd.i,]
    Yb = Y[even.i,]
    Cv_aa = (Ya %*% t(Ya)) / (v_prime_a - 1)  # eq (5) stacked
    Cv_bb = (Yb %*% t(Yb)) / (v_prime_b - 1)  # eq (5) stacked
  }
  
  ### PCA ###
  UDV.aa = svd(Cv_aa, k_prime_a)  # eq. (6) in Rachakonda et al. 2016
  UDV.bb = svd(Cv_bb, k_prime_b)  # eq. (6) in Rachakonda et al. 2016
  
  if (unstacked){  # Calc. PCs for each set of spatial indices...
    X_ab = foreach(sub=iter(datafiles), 
                   .export=c('load.dataFiles', 'zeroMean_Yi'), 
                   .combine='+') %dopar% {
                     Yi = load.dataFiles(sub)[[1]]
                     Yi = zeroMean_Yi(Yi)
                     Xa = Yi %*% t(Yi[odd.i,])  %*% UDV.aa$u # eq (7) & (8) unstacked, as summation across subjs.
                     Xb = Yi %*% t(Yi[even.i,]) %*% UDV.bb$u # eq (7) & (8) unstacked
                     diff.i = sum(odd.i) - sum(even.i)
                     if (diff.i > 0){                        # pad matrix in order to concatenate & format for return
                       Xb = cbind(Xb, matrix(NA, nrow(Xb), diff.i))                     # 6/29/20 --kw-- check concatenation here
                     }else if (diff.i < 0){
                       diff.i = abs(diff.i)
                       Xa = cbind(Xa, matrix(NA, nrow(Xa), diff.i))
                     }
                     stopifnot(dim(Xa) == dim(Xb))
                     return(rbind(Xa, Xb))
                   }
    Xa = X_ab[1:v, 1:v_prime_a]
    Xb = X_ab[v + 1:v, 1:v_prime_b]
  }else{
    Xa = Y %*% t(Ya) %*% UDV.aa$u # eq (7) & (8) stacked
    Xb = Y %*% t(Yb) %*% UDV.bb$u # eq (7) & (8) stacked
  }
  
  Cv2_aa = (t(Xa) %*% Xa) / (v-1)    # eq. (9) left in Rachakonda et al. 2016
  Cv2_ab = (t(Xa) %*% Xb) / (v-1)    # eq. (9) left in Rachakonda et al. 2016
  Cv2_ba = (t(Xb) %*% Xa) / (v-1)    # eq. (9) left in Rachakonda et al. 2016
  Cv2_bb = (t(Xb) %*% Xb) / (v-1)    # eq. (9) left in Rachakonda et al. 2016
  Cv2 = rbind(cbind(Cv2_aa, Cv2_ab), # eq. (9) left in Rachakonda et al. 2016
              cbind(Cv2_ba, Cv2_bb)) # eq. (9) left in Rachakonda et al. 2016
  
  UDV.v2 = svd(Cv2, k)               # eq. (9) right in Rachakonda et al. 2016
  
  ### Group PCA space ###
  X = cbind(Xa, Xb) %*% UDV.v2$u %*% diag( 1 / sqrt(UDV.v2$d[1:k])) # eq. (10) in Rachakonda et al. 2016
  
  ### Eigenvectors & Eigenvalues ###
  #    NOTE: Normalizes Lambda for agreement w/ other fns. here, simplest to correct for here, allows lambda to be compared w/ other PCA fns.
  #          There apears to be a trade-off between taking the sqrt. previously in eq. (8) vs. re-scaling eigenvalues later below.
  #             For instance, if eqs. (7) & (8) are corrected for duplicating of Y's eigenvalues, X will be inaccurate,
  #             but the eigenvalues below will not need to be corrected here.
  Lambda = sqrt(UDV.v2$d[1:k] / (v-1))    # eq. (8) above essentially squares eigenvalues by doubling during concanation
  if (return.eigenvectors){
    if (all(is.na(datafiles))){  datafiles = Y  }
    F_ = get_eigenValsVecs_fromPCs(datafiles, X, unstacked=unstacked)$F_
  }else{  'F_' = NA  }
  
  # if (return.eigenvectors){
  #   
  #   F_ = get_eigenvalues_fromPCs(Y, X,unstacked=unstacked)
  #   
  # }else{  F_ = NA  }
  
  
  # ### Some interesting PCA theory ###
  # #    NOTE: Normalizes Lambda for agreement w/ other fns. here, simplest to correct for here, allows lambda to be compared w/ other PCA fns.
  # #          There apears to be a trade-off between taking the sqrt. previously in eq. (8) vs. re-scaling eigenvalues later below.
  # #             For instance, if eqs. (7) & (8) are corrected for duplicating of Y's eigenvalues, X will be inaccurate,
  # #             but the eigenvalues below will not need to be corrected here.
  # Lambda = sqrt(UDV.v2$d[1:k] / (v-1))    # eq. (8) above essentially squares eigenvalues by doubling during concanation
  
  
  #############################################
  return(list('X'          = X,
              'F_'         = F_,
              'Lambda'     = Lambda,         #Lambda scaled for agreement with other PCA methods
              'Lambda_eq9' = UDV.v2$d[1:k]   #Lambda from EVD of equation (9) in paper, used to scale X appropriately
  ))
} #############################################


#############################################
PCA.MultiPowerIteration <- function(Y, k=NA, X0=NA, l=5, delta=0.001, 
                                    Y.dims=NA, unstacked=F, return.eigenvectors=F,
                                    verbose=F){
  #############################################
  # PCA based on Multi-Power Iteration (MPOWIT)
  # 
  #   Rachakonda S, Silva RF, Liu J, Calhoun VD. (2016). Memory Efficient PCA Methods for
  #     Large Group ICA. Frontiers in Neuroscience (10): article 17.
  #
  # INPUT:
  #   Y : subj.-level data,
  #         as temp.-concatenated data,    [v by Mp] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   X0 : matrix used to initialize group-level PCA space,
  #           either gaussian (default),
  #           'STP' for subsampled Time PCA,
  #           'SVP' for subsampled Voxel PCA
  #   l : "levels" of k, 
  #         essentially, multiplier used to exceed large number of comps. estimated,
  #         for accurate convergence
  #   delta : tolerance for change in eigenvalues between iterations,
  #             threshold used to define convergenceS
  #   Y.dims : (optional) dims. of Y used for odd/even spatial subsampling by voxel,
  #             defaults to every other row of Y. Unused w/o subsampled voxel initiation of X0
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #               If unstacked=T, will run in parallel using R doParallel library
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # OUTPUT:
  #   X : group-level PCA space,            [v by k] matrix
  #   F : group-level eigenvector basis     [v by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #   delta  : tolerence level
  #
  # REQUIRES:
  #   zeroMean_Yi()       : data munging
  #   concatTemporal_Yi() : data munging for arrays & lists
  #   get_data_dims()     : find number of subjs., voxels, time points, etc.
  #   PCA.SubsampledVoxel()
  #   spatialSubsample()  : sampling, taking into account spatial structure of data
  #   PCA.SubsampledTime()
  #
  # EQUATIONS:
  #   X_0 = G_v,lk,  Chi_0 = (YY^T)X_0,  Lambda_0 = 0          eq (28)
  #   X_j = ortho(Chi_j_1) = Chi_j_1 %*% F %*% L^-1 >= 1       eq (29)
  #   Chi_j = (YY^T)X_j = Y ( X_j^T %*% Y)^T = sum(Yi(X_j^T %*% Yi)^T)  eq (30)
  #   Chi_j^T %*% Chi_j / (v-1) = W_j %*% Lambda_j %*% W_j^T   eq (31)
  #   || Lambda_j  -  Lambda_j_1 || < delta                    eq (32)
  #   X = X_j %*% W_j %*% Lambda_j^-1/2                        eq (33)
  #   
  # SPEED:   scales w/   O(Mpvlkj)   operations
  #
  
  if (verbose){cat("\nCalculating PCA space using Multiple Power Iteration (MPOWIT)...")}
  
  ### Required Internal fns. ###
  L2.norm <- function(Lambda){  sqrt(sum(Lambda^2))  }
  Sup.norm <- function(Lambda){  max(abs(Lambda))    } # Suprenum norm for finite-dim. vectors
  ortho <- function(X){
    F_ = svd(t(X) %*% X)$u
    L = apply(X %*% F_, 2, L2.norm) # get L2 norms of cols.
    return(X %*% F_ %*% diag( 1 / L))
  }
  
  ### Defaults & Input ###
  stacked = check_stacked(Y) || !as.logical(unstacked)
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  
  l = as.integer(l)
  stopifnot(l > 0)
  
  data.dims = get_data_dims(Y, unstacked, return.data=stacked)
  M = data.dims$M
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  if (stacked){  Y = data.dims$Y  }
  stopifnot(exists('M')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = min(v, M*p, k, na.rm=T)
  lk = min(l*k, v, M*p, na.rm=T)
  
  ### Initialize vars. ###
  if (all(is.na(X0))){
    X_j = matrix(rnorm(v*lk), v, lk)            # eq (28) 
    Lambda_j = matrix(0, k, k)                  # eq (28) 
  }else if (is.character(X0)){
    stopifnot(X0 %in% c('SVP', 'STP'))
    if (X0 == 'STP'){
      if (verbose){cat('\n...initializing MPOWIT w/ Subsampled Time PCA...')}
      pca.y = PCA.SubsampledTime(Y, g=4, unstacked=unstacked)
    }else if (X0 == 'SVP'){
      if (verbose){cat('\n...initializing MPOWIT w/ Subsampled Voxel PCA...')}
      if (verbose && !all(is.na(Y.dims))){cat('& spatial dims.: (',paste(Y.dims,collapse=','),')\n')}
      pca.y = PCA.SubsampledVoxel(Y, Y.dims=Y.dims, unstacked=unstacked)
    }
    if (length(pca.y$Lambda) <= lk){
      l0 = lk - length(pca.y$Lambda)
      Lambda_j = c(pca.y$Lambda, rep(0, l0))
      X_j = cbind(pca.y$X, matrix(rnorm(v*l0),v,l0))
    }else{
      Lambda_j = pca.y$Lambda[1:lk]
      X_j = pca.y$X[,1:lk]
    }
  }else if (is.matrix(X0)){
    stopifnot(nrow(X0) == v)
    stopifnot(ncol(X0) == lk)
  }else{  stopifnot(is.matrix(X0))  }
  
  if (verbose){cat('\n...running MPOWIT:')}
  if (unstacked){
    Chi_j = foreach(sub=iter(datafiles), 
                    .export=c('load.dataFiles',
                              'zeroMean_Yi',
                              'concatTemporal_Yi'),
                    .combine='+') %dopar% {
                      Yi = load.dataFiles(sub, verbose=verbose)[[1]]
                      Yi = zeroMean_Yi(Yi)
                      Yi = concatTemporal_Yi(Yi)
                      return(Yi %*% t(t(X_j) %*% Yi)) # eq (30) right hand side, j=0
                    }
  }else{
    Y = zeroMean_Yi(Y)
    Y = concatTemporal_Yi(Y)
    Chi_j = Y %*% t(t(X_j) %*% Y)               # eq (28), j=0
  }
  
  # lambda.norm = L2.norm(Lambda_j)                 # eq (32), L2 norm, using Lambda_0 = 0
  lambda.norm = Sup.norm(Lambda_j)              # eq (32), Suprenum norm
  if (lambda.norm == 0){  lambda.norm = Inf  }  # ensure at least 1 iteration
  
  
  ### Main fn. ###
  j = 0
  while ((lambda.norm > delta) && (j < 1000)){
    j = j + 1
    Lambda_j_1 = Lambda_j
    if (verbose){cat(paste0('\n...  iteration  j=',j,'  '))}
    
    X_j = ortho(Chi_j)                          # eq (29)
    if (unstacked){
      Chi_j = foreach(sub=iter(datafiles), 
                      .export=c('load.dataFiles',
                                'zeroMean_Yi'), 
                      .combine='+') %dopar% {
                        Yi = load.dataFiles(sub)[[1]]
                        Yi = zeroMean_Yi(Yi)
                        return(Yi %*% t(t(X_j) %*% Yi))   # eq (30) right hand side
                      }
    }else{
      Y = zeroMean_Yi(Y)
      Y = concatTemporal_Yi(Y)
      Chi_j = Y %*% t(Y) %*% X_j                # eq (30)
    }
    
    
    UDV.XChi = svd((t(X_j) %*% Chi_j) / (v-1))  # eq (31)
    Lambda_j = UDV.XChi$d[1:k]
    
    # lambda.norm = L2.norm(Lambda_j - Lambda_j_1)  # eq (32)
    lambda.norm = Sup.norm(Lambda_j - Lambda_j_1)  # eq (32)
    if (verbose){cat(paste0('...  || Lambda_j - Lambda_j_1 || =  ', lambda.norm))}
  }
  if (verbose){
    if (lambda.norm >= delta){
      cat(paste0('\n...... stopped at j=',j,', returning last iteration\n'))
    }else{
      cat(paste0('\n......converged at j=',j,'  by criteria   || Lambda_j - Lambda_j_1 || <  ', delta, ' = delta\n'))
    }
  }
  
  ### Group PCA space ###
  W_j = UDV.XChi$u[,1:k]
  X = X_j %*% W_j %*% diag( 1 / sqrt(Lambda_j))   # eq (33)
  X = X %*% diag( 1 / apply(X,2,sd))              # whitten principal components
  
  # # NOTE: interestingly, no other PCA algorithm requires whittening group PCs
  
  ### Eigenvectors & Eigenvalues ###
  if (return.eigenvectors){
    if (all(is.na(datafiles))){  datafiles = Y  }
    F_ = get_eigenValsVecs_fromPCs(datafiles, X, unstacked=unstacked)$F_
  }else{  'F_' = NA  }
  # if (return.eigenvectors){
  #   
  #   F_ = get_eigenvalues_fromPCs(Y, X,unstacked=unstacked)
  #   
  # }else{  F_ = NA  }
  # 
  # # # NOTE: for eigenvectors, project group PCA X onto transpose of data, normalize cols. to L2 norm.
  
  #############################################
  return(list('X'          = X,
              'F_'         = F_,
              'Lambda'     = Lambda_j,
              'delta'      = delta,
              'j'          = j))
} #############################################



