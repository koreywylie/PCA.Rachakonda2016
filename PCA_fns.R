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
#   S : number of subjs.
#   v : number of voxels
#   t : number of time points
#   p : number of subj.-level top components, after subj.-level PCA
#   k : number of group-level top components, after group-level PCA
#   Zi : data for subj. i,                [v by t] matrix
#   Yi : whitened data for subj. i,       [v by p] matrix
#   Y : temp.-concatenated group data,    [v by S*p] matrix
#   X : group-level PCA space,            [v by k] matrix
#
#
# Theory:
#   Principal Components Analysis (PCA) projects the data onto 
#     the subspace defined by the eigenvectors of the cov. matrix.
#     For instance, if Y is the [v x t] data matrix w/ zero mean, then
#
#         t(Y) %*% Y / (v-1) = U %*% Lambda %*% t(U),
#
#     where U is the matrix of eigenvectors in cols., 
#     and Lambda is a diagonal matrix of eigenvalues.
#     The PCA space X is then the pojection of the data onto the eigenvectors
#
#         X = Y %*% U,
#
#     with optional whitening of X by scaling eigenvectors to eigenvalues.
#
#   Eigenvalue Decomposition (EVD) and Singular Value Decomposition (SVD)
#     are used to calculate PCA. EVD is SVD applied to matrices.
#     SVD decomposes a matrix Y using a dual basis with associated singular values
#
#         Y = V %*% Lambda^1/2 %*% U * sqrt(v-1),           (1)
#
#     where V is the [v x t] matrix of voxel eigenvectors,
#     U is the [t x t] matrix of temporal eigenvectors,
#     Lambda^1/2 is a diagonal matrix of the square roots of cov. matrix eigenvalues.
#     Using eq. (1), the cov. matrices with EVD are then
#
#       C = t(Y) %*% Y / (v-1) = U %*% Lambda %*% U         (2)
#       Cv = Y %*% t(Y) / (v-1) = V %*% Lambda %*% V.       (3)
#
#     Additionally, from (1) the whitened PCA matrix X can be calculated
#     using either of the equivalence relations
#
#       V * sqrt(v-1)  =  X  =  Y %*% U %*% Lambda^-1/2.    (4)
#
#     Finally, the eigenvalues can be recovered from the PCA space and original data,
#     by left & right multiplying (3) by X & t(X) respectively, scaling to (v-1)^2,
#     then substituting (1) and either side of (4) into the resulting equation, as
#
#       Lambda = t(X) %*% Y %*% t(Y) %*% X / (v-1)^2.       (5)
#
#     Subsequently, the eigenvectors U can be recovered from the PCA space & data,
#     by re-arranging and scaling (4) right as
#
#       U = t(Y) %*% X %*% Lambda^-1/2 / (v-1).             (6)
#
#
#   Subsampled Time PCA and Multipower Iteration (MPOWIT)
#     are stochastic PCA algorithms, able to run in parallel,
#     used to efficiently estimate PCA spaces with high dimensions.
#     Equations (5) & (6) are used below, ignoring scaling factors that cancel,
#     to estimate the eigenvectors from the data for stochastic PCA algorithms.
#   Subsampled Time PCA samples temporally-concatenated subjects,
#     then combines PCs using PCA. With a single subj. per group, this algorithm
#     MELODIC's incremental group PCA (MIGP, Smith et al. 2014).
#   MPOWIT uses higher powers of the cov. matrix to speed computation,
#     initialized using subsampled Time PCA to speed convergence.
#     


library(doParallel)


#----------------------------- Main fn. ----------------------------------------------------
#   Chooses among featured PCA algorithms, based on available memory vs. size of dataset,
#     runs PCA algorithm,
#       returning PC space (group or otherwise), eigenvectors (if indicated), eigenvalues
#

#############################################
PCA.Rachakonda2016 <- function(Y, k = NA, 
                               var = NA, var.subset = NA, mask = NA,
                               whitenPCs = T,
                               return.eigenvectors = T, 
                               return.fns.list = F,
                               algorithm = NA, unstacked = T,
                               verbose = F, ...){
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
  #   Y : temp.-concatenated group data,    [v by S*p] matrix
  #   k : number of group-level top components, after group-level PCA
  #   var        : if Y is an .RData file(s), name of saved var. w/ data
  #   var.subset : if 'var' is a list, name of list element to load
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs           : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #   algorithm           : algorithm PCA algorithm to call, fn. name as a string
  #                           see below for named options
  #   unstacked           : with input algorithm above, if True run as unstacked (low mem., parallelizable),
  #                                                   else if False run as stacked (high mem., fast)
  #   ...   : additional inputs passed to specific PCA fns.
  #             (e.g. S_g=10 for group size in PCA.SubsampledTime())
  #
  # Outputs:
  #   X : group-level PCA space,            [v by k] matrix
  #   U : group-level eigenvector basis     [p by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # Dependencies:
  #    (pretty much every function included below)
  #
  
  if (verbose){  cat('\nRunning PCA with algorithms from Rachaconda et al. 2016')}
  
  ### PCA algorithms available in "PCA_fns.R" ###
  algorithm.options = c("PCA.EVDtime", 
                        "PCA.EVDvoxel", 
                        "PCA.SubsampledTime",
                        "PCA.MultiPowerIteration")
  if (return.fns.list){  return(algorithm.options)  }
  
  
  ### Input/default handling ###
  if (is.character(algorithm) 
      && (algorithm %in% algorithm.options)
      && is.logical(unstacked)){
    if (is.na(unstacked)){  unstacked = T  }
    unstacked = unstacked
    if (verbose){ cat(paste0('\n...selecting input algorithm:  ', algorithm))}
    
  }else{
    if (!all(is.na(algorithm))){
      cat(paste0('\nCould not understand requested algorithm: ',algorithm))
    }
    recommendations = find_PCAtype_best(Y, k=k, 
                                        var=var, var.subset=var.subset, mask=mask,
                                        verbose=verbose, ...)
    algorithm = recommendations$PCAtype
    unstacked = recommendations$unstacked
  }
  
  #-------------------------------------------------------------
  if (algorithm     == "PCA.EVDtime"){
    if (verbose){cat("  ...running PCA.EVDtime...")}
    
    PCA = PCA.EVDtime(Y, k=k, 
                      var=var, var.subset=var.subset, mask=mask,
                      whitenPCs=whitenPCs,
                      verbose=verbose, ...)
    PCA$algorithm = algorithm
    PCA$unstacked = NA
    
    #-------------------------------------------------------------
  }else if (algorithm     == "PCA.EVDvoxel"){
    if (verbose){cat("  ...running PCA.EVDvoxel...")}
    
    PCA = PCA.EVDvoxel(Y, k=k, 
                       var=var, var.subset=var.subset, mask=mask,
                       unstacked=unstacked, 
                       whitenPCs=whitenPCs,
                       return.eigenvectors=return.eigenvectors,
                       verbose=verbose, ...)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
    #-------------------------------------------------------------
  }else if (algorithm     == "PCA.SubsampledTime"){
    if (verbose){cat("  ...running PCA.SubsampledTime...")}
    
    PCA = PCA.SubsampledTime(Y, k=k, 
                             var=var, var.subset=var.subset, mask=mask,
                             unstacked=unstacked,
                             whitenPCs=whitenPCs,
                             return.eigenvectors=return.eigenvectors,
                             verbose=verbose, ...)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
    #-------------------------------------------------------------
  }else if (algorithm     == "PCA.MultiPowerIteration"){
    if (verbose){cat("  ...running PCA.MultiPowerIteration...")}
    
    X0 = 'STP'  # initialize w/ Subsampled Time PCA, to speed convergence
    PCA = PCA.MultiPowerIteration(Y, k=k, 
                                  X0=X0, l=5, delta=0.000001, 
                                  var=var, var.subset=var.subset, mask=mask,
                                  unstacked=unstacked,
                                  whitenPCs=whitenPCs,
                                  return.eigenvectors=return.eigenvectors,
                                  verbose=verbose, ...)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
  }  #-------------------------------------------------------------
  
  #############################################
  return(PCA)
} #############################################



#----------------------------- Data Munging fns. --------------------------------------------------------------------------
#############################################
load.data <- function(datafiles, 
                      var = NA, var.subset = NA, mask = NA, 
                      verbose = F){
  #############################################
  # Loads saved var. from datafiles, subj. by subj.,
  #   including 3d & 4d nifti volumes applying (optional) mask,
  #   or specified elements if var. is a list.
  # 
  # Input:
  #   datafiles : paths to saved files
  #   var       : named of saved variable to load
  #   var.subset: if 'var' is a list, name of list element to load
  #   mask      : binary mask for nifti vols.
  # Requires:
  #   flatten_img() for nifti volumes
  #   resample_img() if mask is not in space of nifti vol.
  #
  
  if (is.list(datafiles) || is.array(datafiles)){
    return(datafiles)
  }else{
    datafiles = as.character(datafiles)
  }
  if (all(dir.exists(datafiles))){
    datafiles = list.files(datafiles, full.names=T)
    datafiles = datafiles[!dir.exists(datafiles)]
  }
  stopifnot(all(file.exists(datafiles)))
  
  if (!all(is.na(mask)) && is.character(mask)){
    stopifnot(file.exists(mask))
  }
  
  Y = list()
  for (f in 1:length(datafiles)){
    if (verbose){cat(paste0('\nLoading: ', datafiles[f]))}
    if (grepl('*.nii$', datafiles[f])){
      Y_i = flatten_img(datafiles[f], mask)
    }else{
      sub = new.env()
      load(datafiles[f], envir=sub)
      
      if (all(is.na(var)) || !exists(var, envir=sub)){
        if (length(ls(sub)) > 1){
          cat(paste0('WARNING:   data var. not specified or not found. '))
          cat(paste0('Could not find data in ', datafiles[f]))
        }
        stopifnot(length(ls(sub))==1)
        var = ls(sub)
      }
      stopifnot(exists(var, envir=sub))
      Y_i = get(var, envir=sub)
      
      if (is.list(Y_i) && (var.subset %in% names(Y_i))){
        Y_i = Y_i[[var.subset]]
      }
      stopifnot(is.numeric(Y_i))
      
    }
    if (verbose){cat('\n')}
    if (verbose){cat(paste0(str(Y_i)))}
    Y = c(Y, list(Y_i))
    
  }
  if (length(Y)==1){  Y = Y[[1]]  }
  
  #############################################
  return(Y)
} #############################################


#############################################
resample_img <- function(source_img, target_img){
  #############################################
  # Resample/reslice source image (fMRI nifti vol.)
  #   to spatial dimensions of target image
  #
  
  if (is.character(source_img) && file.exists(source_img)){
    source_img = RNifti::readNifti(source_img)
  }
  if (is.character(target_img) && file.exists(target_img)){
    target_img = RNifti::readNifti(target_img, volumes=1)
  }else if (length(dim(target_img)) == 4){
    target_img = RNifti::asNifti(target_img[,,,1], niftiHeader(target_img))
  }
  if (!all(class(source_img)==c("niftiImage","array"))){
    print("WARNING:  source image must be nifti image loaded by RNifti")
    if (is.array(source_img)){
      print("WARNING:  input source_img is an array, will distort data into block of voxels in upper right anterior of vol.")
    }
  }
  if (!all(class(target_img)==c("niftiImage","array"))){
    print("WARNING:  target image must be nifti image loaded by RNifti")
  }
  
  D = min(3, length(dim(target_img)))
  if (all(dim(source_img)[1:D] == dim(target_img)[1:D]) &&
      all(RNifti::pixdim(source_img)[1:D] == RNifti::pixdim(target_img)[1:D]) &&
      all(RNifti::origin(source_img) == RNifti::origin(target_img))){
    return(source_img)  # all dims verified, nothing to do
  }
  
  #############################################
  resampled = RNiftyReg::niftyreg.linear(source_img, 
                                         target_img,
                                         init=diag(4), nLevels=0)
  #############################################
  return(resampled$image)
} #############################################


#############################################
flatten_img <- function(img, mask = NA, img3d = NA){
  #############################################
  # Flattens 4d nifti vol., loaded w/ readNifti()
  #   applying optional mask,
  #     returns matrix as [voxels x time].
  # If mask is not the same dimensions as 4D img,
  #   either 3D ref. image or filename in 'img3d' required to reshape mask. 
  #
  # Reversed by:
  #   unflatten_img()
  # Requires: 
  #   resample_img() if dims. of input differ
  #
  
  if (is.character(img) && file.exists(img)){
    img3d = img
    img = RNifti::readNifti(img)
  }
  
  if (!all(is.na(mask))){
    if (is.character(mask) && file.exists(mask)){
      mask_filename = mask
      mask = RNifti::readNifti(mask)
    }
    D = min(3, length(dim(mask)))
    if (!all(dim(img)[1:D] == dim(mask))){
      if ((length(dim(img)) > 3) && !is.na(img3d)){
        if (is.character(img3d) && file.exists(img3d)){
          img1 = RNifti::readNifti(img3d, volumes=1)
          
        }else if (all(class(img3d)==c("niftiImage","array")) &&
                  (length(dim(img3d)) == 3)){
          img1 = img3d
        }
        mask = resample_img(mask, img1)  # target img must be 3D vol.
        
      }else if (length(dim(img)) > 3){
        print('WARNING:  filename for img not provided, cannot re-shape mask volume!')
        mask = resample_img(mask, img)
      }
    }
    stopifnot(all(dim(img)[1:D] == dim(mask)))
    if (all(c("niftiImage","array") %in% class(mask))){
      stopifnot(all(RNifti::pixdim(img)[1:D] == RNifti::pixdim(mask)))
      stopifnot(all(RNifti::origin(img) == RNifti::origin(mask)))
    }
    if (any(is.na(mask))){  mask[is.na(mask)] = 0   }
    mask = mask != 0
  }else{
    mask = rep(T, prod(dim(img)[1:D], na.rm=T))
  }
  
  if (length(dim(img)) == 4){ 
    # flatten along 4th dim & concat in rows
    flat_img = apply(img, 4, function(x) x[mask])
  }else{
    flat_img = img[mask]
  }
  #############################################
  return(flat_img)
} #############################################


#############################################
unflatten_img <- function(data, mask){
  #############################################
  # Fills in 4d-array with data,
  #   using 3d binary mask,
  #     w/ data formatted as [voxels x ~time] matrix.
  # If mask is nifti, returns 4d-nifti vol.
  #
  # Reversed by:
  #   flatten_img() & appropriate mask
  #

  library(RNifti)
  library(abind)
  library(foreach)
  library(iterators)

  data = as.matrix(data)

  if (is.character(mask) && file.exists(mask)){
    mask_filename = mask
    mask = RNifti::readNifti(mask)
    mask.hdr = RNifti::niftiHeader(mask)
  }else if (all(c("niftiImage","array") %in% class(mask))){
    mask.hdr = RNifti::niftiHeader(mask)
  }else{
    stopifnot(length(dim(mask)) == 3)
    mask.hdr = NA
  }
  if (any(is.na(mask))){
    mask[is.na(mask)] = 0
  }
  if (length(dim(mask)) < 3){
    dim(mask) = c(dim(mask), 1)
    vx.size = round(mean(RNifti::pixdim(mask)))
    RNifti::pixdim(mask) = c(RNifti::pixdim(mask), vx.size)
  }
  mask = mask != 0

  ### Reshape into array ###
  assembly <- function(...){  do.call(abind, list(..., 'along'=4))  }
  img = foreach(y=iapply(data,2),
                .combine=assembly,
                .multicombine=T) %do% {
                  vol = mask
                  vol[mask] = y
                  return(vol)
                }

  ### Format as NifTI-1 4d vol. ###
  if (!all(is.na(mask.hdr))){
    mask.hdr$dim[1] = length(dim(img))
    if (length(dim(img)) >= 4){
      mask.hdr$dim[5] = dim(img)[4]
      mask.hdr$pixdim[5] = 1
    }
    attr(mask.hdr, 'imagedim') = c(mask.hdr$dim[2:5])
    attr(mask.hdr, 'pixdim') = c(mask.hdr$pixdim[2:5])

    img = asNifti(img, reference=mask.hdr)
  }

  #############################################
  return(img)
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


#############################################
get_data_dims <- function(Y, 
                          var = NA, var.subset= NA, mask = NA,
                          unstacked = T, return.data = F, verbose = F){
  #############################################
  # Find dimensions of data, 
  #   loading 1st subj. if needed
  #   finding saved 'var' or applying 'mask' or 'var.subset' if indicated
  #
  # Output:
  #   v : length of 1st dim. ~ number of voxels
  #   S : number of subjs. (2nd dim. is S*p)
  #   p : number of time points/ICs (~2nd dim.)
  #   datafiles : paths to saved files
  #   Y : (optional) if indicated by return.data=T,
  #         loaded data for all subjs.
  #
  
  ### Required fns. ###
  get.spacetime <- function(point){
    v = p = NA
    if (is.list(point)){
      Y.dims = get.spacetime(point[[1]])
      v = Y.dims[[1]]
      p = Y.dims[[2]]
    }else if (is.matrix(point)){
      v = nrow(point)
      p = ncol(point)
    }else if (is.array(point)){
      Y.dims.all = dims(point)
      Y.dims = Y.dims.all[-length(Y.dims.all)] # Assume leading dimensions are all spatial...
      v = prod(Y.dims)                         # ...number of voxels is their product...
      p = Y.dims.all[length(Y.dims.all)]       # ...& assume time/component is the last dimension
    }else if (is.vector(point)){
      v = length(point)                # assume leading dimensions are spatial
      p = 1                             # ...assume time/component is the last dimension
    }
    return(list('v'=v, 'p'=p))
  }
  
  
  ### Inputs/defaults ###
  stacked = !as.logical(unstacked)
  if (is.character(Y) && dir.exists(Y)){
    datafiles = list.files(path=Y, full.names=T)
    datafiles = datafiles[!dir.exists(datafiles)]
    S = length(datafiles)
    
    if (stacked){  # load all subjs. into memory at once
      Y = load.data(datafiles, var, var.subset, mask, verbose=verbose)
      spacetime = get.spacetime(Y[[1]])
    }else{
      Y = load.data(datafiles[1], var, var.subset, mask, verbose=verbose)
      spacetime = get.spacetime(Y)
    }
    
  }else if (is.character(Y) && all(file.exists(Y))){
    if (stacked){  # load all subjs. into memory at once
      datafiles = Y
      if (verbose){cat(paste0('Loading & stacking: ', paste(basename(Y), collapse=', ')), '...\n')}
      Y = load.data(Y, var, var.subset, mask, verbose=verbose)
      S = length(datafiles)
      spacetime = get.spacetime(Y[[1]])
      
    }else{            # load subjs. one-by-one
      datafiles = Y
      S = length(datafiles)
      Y = load.data(datafiles[1], var, var.subset, mask)
      spacetime = get.spacetime(Y)
      
    }
  }else{
    if (is.character(Y) && !all(file.exists(Y))){
      cat(paste0('WARNING: could not find file as input:'))
      print(Y)
      stopifnot(all(file.exists(Y)))
    }else{
      if (is.list(Y)){
        S = length(Y)
        spacetime = get.spacetime(Y[[1]])
        
      }else if (is.array(Y)){
        Y.dims = dim(Y)
        if (length(Y.dims) > 2){
          S = Y.dims[1]
          spacetime = list(prod(Y.dims[2:(length(Y.dims)-1)]),
                           Y.dims[length(Y.dims)])
        }else{
          S = 1
          spacetime = as.list(Y.dims)
          
        }
      }else{  stopifnot(is.list(Y) || is.array(Y))  }
    }
    datafiles = NA
  }
  v = spacetime[[1]]
  p = spacetime[[2]]
  
  if (!return.data){  Y = NA  }
  
  #############################################
  return(list('S'=S, 'v'=v, 'p'=p, 
              'datafiles'=datafiles,
              'Y'=Y))
} #############################################

#############################################
calc_Yi_means <- function(Y){
  #############################################
  # Calculates column means of subj.-level data Y,
  #   where Y is either list of Yi   [v by p] matrices
  #     or Y is array with dim.   [S by v by p]
  #
  
  if (!all(is.finite(unlist(Y)))){
    print('WARNING:  Infinite values found in matrix!  Mean(s) will be +/-Inf')
  }
  if (is.matrix(Y)){
    mY = apply(Y, 2, mean)
  }else if (is.list(Y)){
    mY = lapply(Y, function(yy) apply(yy, 2, mean))
  }else if (is.array(Y)){    # returns list of subj. column vectors
    mY = t(apply(Y, 1, function(yy) apply(yy, 2, mean)))
    mY = lapply(apply(mY,1,list), function(yy) t(matrix(unlist(yy))))
  }
  return(mY)
} #############################################

#############################################
zeroMean_Yi <- function(Y){
  #############################################
  # Removes column mean of subj.-level data Y
  #     ~centering every vol. at every timepoint
  #   where Y is either list of Yi   [v by p] matrices
  #     or Y is array with dim.   [S by v by p]
  #
  
  if (!all(is.finite(unlist(Y)))){
    print('WARNING:  Infinite values found in matrix!  Mean(s) will be +/-Inf')
  }
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
} #############################################

#############################################
concatTemporal_Yi <- function(Y){
  #############################################
  # Temporally Concatenates subj.-level data Y
  #   where Y is either list of Yi   [v by p] matrices
  #   or Y is array with dim.   [S by v by p]
  #
  
  if (is.matrix(Y)){
    return(Y)  # PASS, nothing to do
  }else if (is.list(Y)){
    # iterates over elements of Y & binds results by column
    return(do.call(cbind, Y))
  }else if (is.array(Y)){
    # re-arranges Y, putting subj. dimension last,
    #   iterates over ~voxel dim., binding results by column
    return(t(apply(aperm(Y,c(2,3,1)),1,cbind)))
  }
} #############################################

#############################################
partition_group <- function(subj.count, group.size){ 
  #############################################
  # Randomly partitions subjs. into groups,
  #   with constant group size,
  #     & return list of indice vectors w/ subjs. in each group.
  #
  
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
  
  #############################################
  return(partition)
} #############################################


#############################################
get_eigenValsVecs_fromPCs <- function(Y, X,
                                      var = NA, var.subset = NA, mask = NA,
                                      unstacked = F, verbose = F){
  #############################################
  # Find eigenvalues & eigenvectors from data Y & (optional) variable var,
  #   given input Principal Components X.
  # Used to get eigenvalues
  #   following PCA methods that only calc. PCs
  #
  # Theory: 
  #   If Y is the [v x S*p] data matrix,
  #     & X is the [v x k] group PCs,
  #     solve for the [S*p x k] matrix U of eigenvectors,
  #     using def. of the Principal Components of Y,
  #     & squared eigenvals. from scaling of PCs:
  #
  #   (t(Y) %*% Y) = U %*% Lambda %*% t(U)
  #
  #            X =          Y %*% U
  #   t(Y) %*% X = (t(Y) %*% Y) %*% U
  #   t(Y) %*% X = U %*% Lambda %*% t(U) %*% U
  #   t(Y) %*% X = U %*% Lambda = U_scaled
  #
  #     Lambda^2 = t(Lambda) %*% Lambda
  #     Lambda^2 = t(Lambda) %*% (t(U) %*% U) %*% Lambda
  #     Lambda^2 = t(U_scaled) %*% U_scaled
  #     Lambda = sqrt(diag(t(U_scaled) %*% U_scaled))
  #
  #           U = U_scaled %*% Lambda^-1
  #
  
  stacked = !unstacked
  data.dims = get_data_dims(Y, var, var.subset, mask,
                            unstacked=unstacked, 
                            return.data=stacked)
  if (verbose){print(str(data.dims))}
  datafiles = data.dims$datafiles
  if (stacked){ Y = data.dims$Y } # loaded data
  S = data.dims$S # number of subjs.
  v = data.dims$v # number of voxels
  p = data.dims$p # dimension of voxels (i.e., time, reduced-PCA T1, etc.)
  
  if (unstacked){
    U_scaled = foreach(sub=iter(datafiles),
                       .export=c('load.data','zeroMean_Yi'),
                       .combine = rbind) %dopar% {
                         
                         Yi = load.data(sub, var, var.subset, mask,
                                        verbose=verbose)
                         Yi = zeroMean_Yi(Yi)
                         return( t(Yi) %*% X )
                         
                       }
  }else{
    Y = load.data(Y, var, var.subset, mask)
    Y = zeroMean_Yi(Y)
    Y = concatTemporal_Yi(Y)
    U_scaled = t(Y) %*% X
    
  }
  
  Lambda = sqrt(diag(t(U_scaled) %*% U_scaled))
  U = U_scaled %*% diag(1 / Lambda)
  # U_sd = apply(U,2, sd)   # if needed for normalization
  
  #############################################
  return(list('U'      = U,       # eigenvectors
              'Lambda' = Lambda)) # eigenvalues
} #############################################




#------------------------------------ Memory-calculations & algorithm recommendations ----------------------------
#############################################
find_PCAtype_best <- function(data, k = NA, 
                              var = NA, var.subset = NA, mask = NA,
                              verbose = F, ...){
  #############################################
  # Finds best PCA algorithm,
  #   based on recommendations in Figure 8, Rachakonda et al. (2016)
  #
  
  ### Get data dims ###
  if (is.numeric(data) && (length(data) == 3)){  # avoid loading potentially large datasets
    S = data[1]  # number of subjs.
    v = data[2]  # number of voxels
    p = data[3]  # number of time points / reduced dim.
  }else{
    data.dims = get_data_dims(data, var, var.subset, mask)
    S = data.dims$S
    v = data.dims$v
    p = data.dims$p
  }
  k = min(v, S*p, k, na.rm=T)
  
  ### Get available RAM ###
  memFree = get_RAMavailable()
  
  ### Decision table in Figure 8 ###
  if (min(v, S*p) <= 10000){
    if (v < S*p){
      PCAtype = 'PCA.EVDvoxel'
      unstacked = T
    }else{
      PCAtype = 'PCA.EVDtime'
      unstacked = F
    }
  }else{
    PCAtype = 'PCA.MultiPowerIteration'
    if ((log2(v) + log2(S) + log2(p)) > 9){
      unstacked = T  # large dataset, cannot calculate stacked memory requirements
    }else if (v*S*p > memFree){
      unstacked = T
    }else{
      unstacked = F
    }
  }
  if (verbose){
    cat(paste0('\n      Recommended PCA algorithm:  ', PCAtype))
    if (unstacked){cat(' (unstacked)')}else{cat(' (stacked)')}
  }
  
  ### Reality/Sanity Check ###
  memInfo = calc_PCA_memRequirements(data, k, 
                                     pca_type=PCAtype, 
                                     var=var, var.subset=var.subset, mask=mask, ...)
  if (length(memInfo) > 1){
    memNeeded = min(memInfo[[1]]$mem, memInfo[[2]]$mem, na.rm=T)
  }else{
    memNeeded = memInfo[[1]]$mem
  }
  
  if (verbose){
    cat(paste0('\n        requiring  ',signif(memNeeded/1024),' Mb  of Memory'))
  }
  if (memFree < memNeeded){
    print(paste0("WARNING:  Memory available (",signif(memFree/1024),
                 " Mb) less than needed for PCA (",signif(memNeeded/1024)," Mb)!"))
  }
  
  #############################################
  return(list('PCAtype' = PCAtype,
              'unstacked' = unstacked))
} #############################################

#############################################
get_RAMavailable <- function(verbose = F){
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
} #############################################

#############################################
calc_PCA_memRequirements <- function(data, k, pca_type = 'all', 
                                     var = NA, var.subset = var.subset, mask = NA,
                                     verbose = F, ...){
  #############################################
  # Calculate PCA memory requirements,
  #   for various PCA algorithms in 'PCA_fns.R'
  # Based off appendix in Rachakonda et al. (2016). Memory Efficient PCA Methods for Large Group PCA.
  #                           Frontiers in Neuroscience (10)17: p. 1-17.
  #
  # Input:
  #   data     : input data Y, filepath(s) containing Y, 
  #               or dimensions of data vector formatted as [subjects, voxels, ~time]
  #   k        : final number of PCs
  #   pca_type : fn. in 'PCA_fns.R' or 'all'
  #   var        : if Y is file(s), name of saved var. w/ data
  #   var.subset : if 'var' is a list, name of list elements to load
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  # Output: 
  #     list w/ elements indexed by fn. name, containing required RAM in Kb
  # Requires:
  #   get_data_dims() from 'PCA_fns.R'
  #
  
  if (verbose){
    cat(paste0('\nCalculating memory requirements for PCA algorithms...'))
  }
  
  ### Input handling & defaults ###
  if (pca_type == 'all'){
    pca_types = c('PCA.EVDtime', 
                  'PCA.EVDvoxel',
                  'PCA.SubsampledTime',
                  'PCA.MultiPowerIteration')
  }else{  pca_types = pca_type  }
  bytes = 8 # double precision, numeric data
  if (!exists('l')){  l = 5   }  # default multiplier for MPOWIT
  if (!exists('S_g')){  S_g = 20  }  # default number of subjs. per group in STP
  if (!exists('k_prime') && !exists('k_intermediary')){
    k_prime = max(k*5, 500) # dim. of PCA during group update, defaults suggested number in paper
  }else if (!exists('k_prime') && exists('k_intermediary')){
    k_prime = k_intermediary
  }
  
  ### Get data dims ###
  if (is.numeric(data) && (length(data) == 3)){
    # assume data dims entered & avoid loading potentially large datasets
    S = data[1]
    v = data[2]
    p = data[3]
  }else{
    data.dims = get_data_dims(data, var, var.subset, mask)
    S = data.dims$S
    v = data.dims$v
    p = data.dims$p
  }
  
  ### Main fn. ###
  memInfo = list()
  m = 0
  for (pcaType in pca_types){
    if (pcaType == 'PCA.EVDtime'){
      mem = v * S * p       # total data size stacked...
      + (S * p)^2           # ...size of [time x time] cov. matrix
      + v * k               # ...size of group PCA space X
      unstacked = F
    }else if (pcaType == 'PCA.EVDvoxel'){
      mem = v * p           # total data size for single subjs. unstacked...
      + v^2                 # ...size of [time x time] cov. matrix
      + v * k               # ...size of group PCA space X
      unstacked = T
    }else if (pcaType == 'PCA.SubsampledTime'){
      g1 = min(S_g, S)
      mem = v * g1 * p      # total data size for group...
      + (g1 * p)^2          # ...concat. group cov. matrix
      + 2 * v * k_prime     # ...size of group PCA spaces X_g & X_g+1 in eq. (13)
      + g1 * p * k_prime    # ...size of final group PCA space in eq. (15)
      + 2 * k_prime^2       # ...size of Lambda in eq. (14), as block diag. of Lambda_g & Lambda_g+1, for final group PCA est.
      unstacked = T
    }else if (pcaType == 'PCA.MultiPowerIteration'){
      Ma = 2 * v * l * k    # X_j & X_j-1, for j=1...
      + (l * k)^2           # ...[Chi_j^T x Chi_j] cov. matrix
      Mb = v * l * k        # ...X_j
      + 3 * (l * k)^2       # ...[Chi_j^T x Chi_j] & [Chi_j-1^T x Chi_j-1] during ortho() & Lambda_j-1 
      if ((log2(v) + log2(S) + log2(p)) > 9){
        mem_stacked = NA            # product of dims. exceed R's integer capacity
      }else{
        mem_stacked = v * S * p     # total data size stacked...
        + max(Ma, Mb)               # ...plus size of largest iteration
      }
      mem_unstacked = v * p         # size of single subjs.' data...
      + max(Ma, Mb)                 # ...plus size of largest iteration
      
      ### Include MPOWIT initialization requirements w/ Subsampled Time PCA ###
      g1 = min(S_g, S)
      mem_init = v * g1 * p + (g1 * p)^2 + 2 * v * k_prime 
      + g1 * p * k_prime + 2 * k_prime^2
      
      mem = c(max(mem_init, mem_stacked),
              max(mem_init, mem_unstacked))
      unstacked = c(F, T)
    }
    mem = mem * bytes / 1024        # return size in Kb
    if (length(mem) == 1){
      m = m + 1
      memInfo[[m]] = list('type' = pcaType,
                          'mem' = mem,
                          'units' = 'Kb',
                          'unstacked' = unstacked)
    }else{ # stacked & unstacked MPOWIT
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
PCA.EVDtime <- function(Y, k = NA, 
                        var = NA, var.subset = NA, mask = NA,
                        whitenPCs = T,
                        verbose = F, ...){
  #############################################
  # Principal Components Analysis via Eigenvalue Decomposition (well, actually ~SVD)
  #   to find subj. or group-level PCA space
  #
  # Input:
  #   Y : subjs. or group-level data,
  #         as filename(s)
  #         as temp.-concatenated data,    [v by S*p] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   var  : if Y is file(s), name of saved var. w/ data
  #   var.subset : if 'var' is a list, name of list elements to load
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  # Output:
  #   X : group-level PCA space,            [v by k] matrix
  #   U : group-level eigenvector basis     [S*p by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #   center : center coordinates for data  [S*p by 1] matrix
  #
  # Requires:
  #   load.data()         : data munging
  #   calc_Yi_means()     :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #
  
  if (verbose){cat("\nCalculating PCA space using EVD on time*subjs. dimension...")}
  
  if (is.character(Y) && all(file.exists(Y))){
    Y = load.data(Y, var, var.subset, mask) # loads all subjs. into mem. at once
  }else if (is.list(Y)){
    stopifnot(is.numeric(Y[[1]]))
  }else{  stopifnot(is.numeric(Y))  }
  mY = calc_Yi_means(Y)
  mY = concatTemporal_Yi(mY)
  Y = zeroMean_Yi(Y)
  Y = concatTemporal_Yi(Y)
  v = nrow(Y)
  Sp = ncol(Y)
  k = min(k, v, Sp, na.rm=T)
  
  ### Quick SVD from data ###
  EVD = svd(t(Y), k, 0)    # eigenvectors of data equivalent to eigenvectors of corr. matrix by orthogonality
  EVD$d = EVD$d^2 / (v-1)  # scale to eigenvalues of corr. matrix
  ### Slower EVD from corr. mat. ###
  # C = (t(Y) %*% Y) / (v-1)  # eq. (2), can be costly to calculate for large datasets
  # EVD = svd(C, k)
  
  if (whitenPCs){
    if (verbose){cat('\n...projecting data onto eigenvectors & whitening by scaling to eigenvalues...')}
    X = Y %*% EVD$u %*% diag( 1 / sqrt(EVD$d[1:k]) ) # eq. (4)
  }else{
    if (verbose){cat('\n...projecting data onto eigenvectors...')}
    X = Y %*% EVD$u
  }
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'U'      = EVD$u,       # Eigenvectors
              'Lambda' = EVD$d[1:k],  # Eigenvalues
              'whitenedPCs' = whitenPCs,
              'center' = mY))
} #############################################


#############################################
PCA.EVDvoxel <- function(Y, k = NA, 
                         var = NA, var.subset = NA, mask = NA,
                         whitenPCs = T,
                         unstacked = F, 
                         return.eigenvectors = F,
                         verbose = F, ...){
  #############################################
  # Principal Components Analysis via Eigenvalue Decomposition (well, actually ~SVD)
  #   on voxel by voxel covariance matrix
  #     to find group-level PCA space
  #
  # Input:
  #   Y : subj.-level data,
  #         as vector of paths to saved subj. data,
  #         as temp.-concatenated data,    [v by S*p] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   var  : if Y is file(s), name of saved var. w/ data
  #   var.subset : if 'var' is a list, name of list elements to load
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #                       (NOTE: requires parallelization cluster set up)
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # Output:
  #   X   : group-level PCA space,          [v by k] matrix
  #   U  : group-level eigenvector basis    [S*p by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # Requires:
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.data()         :  "     "
  #   calc_Yi_means()     :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #
  # Speed:   time complexity  O(S*p*v^2)
  #
  
  if (verbose){cat("\nCalculating PCA space using EVD on voxel dimension...")}
  
  ### Defaults & Inputs ###
  stacked = check_stacked(Y) || !as.logical(unstacked)
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  
  data.dims = get_data_dims(Y, var, var.subset, mask)
  S = data.dims$S
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  stopifnot(exists('S')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = min(v, S*p, k, na.rm=T)
  
  Cv = matrix(0,v,v)
  if (unstacked){
    Cv = foreach(sub=iter(datafiles), 
                 .export=c('load.data',
                           'zeroMean_Yi'), 
                 .combine='+') %dopar% {
                   
                   Yi = load.data(sub, var, var.subset, mask)
                   Yi = zeroMean_Yi(Yi)
                   return((Yi %*% t(Yi)) / (v-1))  # eq. (3)
                   
                 }
  }else{
    Y = load.data(Y, var, var.subset, mask)
    Y = zeroMean_Yi(Y)
    Y = concatTemporal_Yi(Y)
    Cv = (Y %*% t(Y)) / (v-1)  # eq. (3)
  }
  
  ### Principal Components ###
  EVD = svd(Cv, k, 0)
  X = EVD$u * sqrt(v-1)        # eq. (4) left
  unwhitenPCs = !whitenPCs
  if (unwhitenPCs){
    if (verbose){cat('un-whitening PCs by scaling to inverse of eigenvalues...')}
    X = X %*% diag(sqrt(EVD$d[1:k]))
  }
  
  ### Eigenvectors & Eigenvalues ###
  Lambda = EVD$d[1:k]
  if (return.eigenvectors){
    if (verbose){cat('finding eigenvectors by projecting PCs onto data...')}
    if (all(is.na(datafiles))){  datafiles = Y  }
    
    U = get_eigenValsVecs_fromPCs(datafiles, X, 
                                  var=var, var.subset=var.subset, mask=mask,
                                  unstacked=unstacked)$U
    
  }else{  'U' = NA  }
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'U'      = U,           # Eigenvectors
              'Lambda' = Lambda,      # Eigenvalues
              'whitenedPCs' = whitenPCs))
} #############################################


#############################################
PCA.SubsampledTime <- function(Y, k = NA, k_intermediary = 500, 
                               S_g = 20, parFactor = NA,
                               var = NA, var.subset = NA, mask = NA,
                               whitenPCs = T,
                               return.eigenvectors = F,
                               verbose = F, ...){
  #############################################
  # Principal Components Analysis on sub-sample of time dim.,
  #   i.e., w/ temp. concat. PCA, randomly sample subjs. on time dim.
  #     to find group-level PCA space
  #
  # Input:
  #   Y : subj.-level data,
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after final group-level PCA
  #   k_intermediary : number of components estimated while combining groups
  #   S_g : number of subjs. per group, i.e. group size
  #   parFactor : Parallelization factor, number of worker nodes,
  #                 used to more quickly combine data from large number of groups.
  #   var  : if Y is file(s), name of saved var. w/ data
  #   var.subset : if 'var' is a list, name of list elements to load
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # Output:
  #   X : group-level PCA space,            [v by k] matrix
  #   U : group-level eigenvector basis     [S*p by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # Requires:
  #   PCA.EVDtime()       : main PCA fn.
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.data()         :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #   partition_group()   :  "     "
  #  
  # Speed:   time complexity  O(v(Sp)^2)
  #
  library(foreach)
  library(iterators)
  
  if (verbose){cat("\nCalculating PCA space by subsampling time*subjs. dimension...")}
  
  ############ Required fns. ##################
  calc_GroupSubSampleTimeSVD <- function(data, partition_inds, k_prime,  
                                         var=NA, var.subset=NA, mask=NA, Xg_0=NA){
    ### Calc. for C_g ###
    if (is.character(data)){
      Yg = load.data(data[unlist(partition_inds)], var, var.subset, mask)
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
    EVD.Yg = PCA.EVDtime(Yg, k_prime)           # eq. (11) & (12) in Rachakonda et al. 2016
    
    if (all(is.na(Xg_0))){
      return(Yg %*% EVD.Yg$U)                   # return PCs, if 1st iteration
    }else{
      ### Combine Groups ###
      stopifnot(exists('Xg_0') && all(is.numeric(Xg_0)))
      Xg = Yg %*% EVD.Yg$U                      # eq. (13) in Rachakonda et al. 2016
      
      ### Quick eigenvectors from SVD ###
      EVD.g0g1 = svd(t(cbind(Xg_0, Xg)), k_prime, 0)
      ### Slower Calc. Cov. of g-1 & g ###
      # Cg0g1_g0g0 = (t(Xg_0) %*% Xg_0) / (v-1)        # eq. (14) left
      # Cg0g1_g0g1 = (t(Xg_0) %*% Xg)   / (v-1)
      # Cg0g1_g1g0 = (t(Xg) %*% Xg_0)   / (v-1)
      # Cg0g1_g1g1 = (t(Xg) %*% Xg)     / (v-1)
      # Cg0g1 = rbind(cbind(Cg0g1_g0g0, Cg0g1_g0g1),
      #               cbind(Cg0g1_g1g0, Cg0g1_g1g1))
      # EVD.g0g1 = svd(Cg0g1, k_prime)                 # eq. (14) right in Rachakonda et al. 2016
      
      return(cbind(Xg_0, Xg) %*% EVD.g0g1$u)    # eq. (15)  in Rachakonda et al. 2016
    }
  } ###################
  
  
  ### Defaults & Inputs ###
  data.dims = get_data_dims(Y, var, var.subset, mask)
  S = data.dims$S
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  stopifnot(exists('S')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = as.integer(k)
  k = min(k, v, S*p, na.rm=T)
  k_intermediary = as.integer(k_intermediary)
  k_intermediary = min(k_intermediary, v, S*p, na.rm=T)  # dim. of PCA during group update
  k_intermediary = max(k, k_intermediary)
  
  S_g = min(as.integer(S_g), S, na.rm=T)  # number of subjs. in each group
  parFactor = max(as.integer(parFactor), 1, na.rm=T)  # number of groups on each parallel node
  if ((parFactor > S) || (ceiling(S / parFactor) < S_g)){
    parFactor = 1  # avoid parallelization if not needed (e.g., small S or S_g)
  }
  
  ### Subdivide by groups & EVD on each ###
  if (parFactor > 1){
    S_g_par = ceiling(S / parFactor)
    node.part.list = partition_group(S, S_g_par)  # partition subjs. to parallel nodes
    numNodes = length(node.part.list)
    
    Xg = foreach(n.pt.i=icount(numNodes),
                 .export=c('load.data', 'calc_Yi_means',
                           'zeroMean_Yi', 'concatTemporal_Yi', 
                           'partition_group', 'PCA.EVDtime'),
                 .combine='cbind') %dopar% {
                   
                   dataSubset = datafiles[unlist(node.part.list[[n.pt.i]])]
                   grp.part.list = partition_group(length(dataSubset), S_g)
                   
                   Xg = NA
                   for (pt.i in 1:(length(grp.part.list))){
                     Xg = calc_GroupSubSampleTimeSVD(dataSubset, grp.part.list[[pt.i]], 
                                                     k_intermediary, 
                                                     var, var.subset, mask, Xg)
                   }
                   return(Xg)
                 }
  }else{
    grp.part.list = partition_group(S, S_g)
    Xg = NA
    for (pt.i in 1:(length(grp.part.list))){
      Xg = calc_GroupSubSampleTimeSVD(Y, grp.part.list[[pt.i]],
                                      k_intermediary, 
                                      var, var.subset, mask, Xg)
    }
  }
  EVD.g_final = svd(t(Xg), k, 0)              # combine individual subgroup results
  EVD.g_final$d = EVD.g_final$d^2 / (v-1)  # scale to eigenvalues of corr. matrix
  Xg = Xg %*% EVD.g_final$u
  
  
  ### Group PCA space ###
  if (whitenPCs){
    if (verbose){cat('projecting data onto eigenvectors & whitening by scaling to eigenvalues...')}
    X = Xg %*% diag( 1 / sqrt(EVD.g_final$d[1:k]) ) # Rachakonda et al. (2016) final paragraph for STP section
  }else{
    X = Xg
  }
  
  ### Eigenvectors & Eigenvalues ###
  Lambda = EVD.g_final$d[1:k]
  if (return.eigenvectors){
    if (verbose){cat('finding eigenvectors by projecting PCs onto data...')}
    if (all(is.na(datafiles))){  datafiles = Y  }
    
    U = get_eigenValsVecs_fromPCs(datafiles, X, 
                                  var=var, var.subset, mask=mask,
                                  unstacked=T)$U
    
  }else{  'U' = NA  }
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'U'      = U,           # Eigenvectors
              'Lambda' = Lambda,      # Eigenvalues
              'whitenedPCs' = whitenPCs))
} #############################################


#############################################
PCA.MultiPowerIteration <- function(Y, k = NA, 
                                    X0 = NA, l = 5, delta = 0.001,
                                    S_g = 20, parFactor = NA,
                                    var = NA, var.subset = NA, mask = NA,
                                    whitenPCs = T,
                                    unstacked = F,
                                    return.eigenvectors = F,
                                    verbose = F, ...){
  #############################################
  # PCA based on Multi-Power Iteration (MPOWIT)
  # 
  #   Rachakonda S, Silva RF, Liu J, Calhoun VD. (2016). Memory Efficient PCA Methods for
  #     Large Group ICA. Frontiers in Neuroscience (10): article 17.
  #
  # Input:
  #   Y : subj.-level data,
  #         as temp.-concatenated data,    [v by S*p] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   X0 : matrix used to initialize group-level PCA space,
  #           either gaussian (default), or
  #           'STP' for Subsampled Time PCA
  #   l : "levels" of k, 
  #         essentially, multiplier used to exceed large number of comps. estimated,
  #         for accurate convergence
  #   delta : tolerance for change in eigenvalues between iterations,
  #             threshold used to define convergence
  #   S_g : number of subjs. per group, if initialized w/ Subsampled Time PCA
  #   parFactor : Parallelization factor, number of worker nodes,
  #                 used to more quickly combine data from large number of groups
  #   var  : if Y is file(s), name of saved var. w/ data
  #   var.subset : if 'var' is a list, name of list elements to load
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #               If unstacked=T, will run in parallel using R doParallel library
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # Output:
  #   X : group-level PCA space,            [v by k] matrix
  #   U : group-level eigenvector basis     [S*p by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #   delta  : tolerence level
  #
  # Requires:
  #   PCA.EVDtime()       : main PCA fn.
  #   PCA.SubsampledTime() : used to initialize PCs to speed convergence
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.data()         :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #   partition_group()   :  "     "
  #
  # Equations:    (Rachakonda et al. 2016)
  #   X_0 = G_v,lk,  Chi_0 = (YY^T)X_0,  Lambda_0 = 0          eq (28)
  #   X_j = ortho(Chi_j_1) = Chi_j_1 %*% U %*% L^-1 >= 1       eq (29)
  #   Chi_j = (YY^T)X_j = Y ( X_j^T %*% Y)^T = sum(Yi(X_j^T %*% Yi)^T)  eq (30)
  #   Chi_j^T %*% Chi_j / (v-1) = W_j %*% Lambda_j %*% W_j^T   eq (31)
  #   || Lambda_j  -  Lambda_j_1 || < delta                    eq (32)
  #   X = X_j %*% W_j %*% Lambda_j^-1/2                        eq (33)
  #   
  # Speed:   scales w/   O(Spvlkj)   operations
  #
  
  if (verbose){cat("\nCalculating PCA space using Multiple Power Iteration (MPOWIT)...")}
  
  ### Required Internal fns. ###
  L2.norm <- function(Lambda){  sqrt(sum(Lambda^2))  }
  Sup.norm <- function(Lambda){  max(abs(Lambda))    } # Suprenum norm for finite-dim. vectors
  ortho <- function(X){
    U = svd(t(X), ncol(X), 0)$u
    L = apply(X %*% U, 2, L2.norm) # get L2 norms of cols.
    return(X %*% U %*% diag( 1 / L))
  }
  
  ### Defaults & Input ###
  stacked = check_stacked(Y) || !as.logical(unstacked)
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  l = as.integer(l)
  stopifnot(l > 0)
  data.dims = get_data_dims(Y, var, var.subset, mask, 
                            unstacked, return.data=stacked)
  S = data.dims$S
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  if (stacked){  Y = data.dims$Y  }
  stopifnot(exists('S')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = min(v, S*p, k, na.rm=T)
  lk = min(l*k, v, S*p, na.rm=T)
  
  ### Initialize vars. ###
  if (all(is.na(X0))){
    X_j = matrix(rnorm(v*lk), v, lk)        # eq (28) 
    Lambda_j = matrix(0, k, k)
  }else if (is.character(X0)){
    stopifnot(X0 %in% c('STP'))
    if (X0 == 'STP'){
      if (verbose){cat('\n...initializing MPOWIT w/ Subsampled Time PCA...')}
      pca.y = PCA.SubsampledTime(Y, k=k, k_intermediary=k, 
                                 S_g=S_g, parFactor=parFactor, 
                                 var=var, var.subset=var.subset, mask=mask, 
                                 whitenPCs=whitenPCs)
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
  
  
  #----------------
  ### Main fn. ###
  #----------------
  if (verbose){cat('\n...running MPOWIT:')}
  
  # lambda.norm = L2.norm(Lambda_j)               # eq (32), L2 norm, using Lambda_0 = 0
  lambda.norm = Sup.norm(Lambda_j)              # eq (32), Suprenum norm
  if (lambda.norm == 0){  lambda.norm = Inf  }  # ensure at least 1 iteration
  
  j = -1
  while ((lambda.norm > delta) && (j < 1000)){
    j = j + 1
    Lambda_j_1 = Lambda_j
    if (verbose){cat(paste0('\n...  iteration  j=',j,'  '))}
    
    if (j > 0){
      X_j = ortho(Chi_j)
    }

    if (unstacked){
      parFactor = max(as.integer(parFactor), 1, na.rm=T)  # number of groups on each parallel node
      S_g_par = ceiling(S / parFactor)
      node.part.list = partition_group(S, S_g_par)  # partition subjs. to parallel nodes
      numNodes = length(node.part.list)
      
      Chi.dims = dim(X_j)
      Chi_j = foreach(n.pt=icount(numNodes),
              .export=c('load.data', 'zeroMean_Yi'),
              .combine='+') %dopar% {
                
                dataSubset = datafiles[unlist(node.part.list[[n.pt]])]
                
                Chi_j = matrix(0, Chi.dims[1], Chi.dims[2])
                for (sub in dataSubset){
                  Yi = load.data(sub, var, var.subset, mask)
                  Yi = zeroMean_Yi(Yi)
                  Chi_j = Chi_j + Yi %*% t(t(X_j) %*% Yi)   # eq (30) right hand side
                }
                return(Chi_j)
              }
    }else{
      Y = zeroMean_Yi(Y)
      Y = concatTemporal_Yi(Y)
      Chi_j = Y %*% t(Y) %*% X_j                # eq (30)
    }
    
    EVD.XChi = svd((t(X_j) %*% Chi_j) / (v-1), lk, 0)  # eq (31)
    Lambda_j = EVD.XChi$d[1:k]
    
    # lambda.norm = L2.norm(Lambda_j - Lambda_j_1)  # eq (32)
    lambda.norm = Sup.norm(Lambda_j - Lambda_j_1)  # eq (32)
    if (verbose && (j > 0)){cat(paste0('...  || Lambda_j - Lambda_j_1 || =  ', lambda.norm))}
  }
  if (verbose){
    if (lambda.norm >= delta){
      cat(paste0('\n...... stopped at j=',j,', returning last iteration\n'))
    }else{
      cat(paste0('\n......converged at j=',j,'  by criteria   || Lambda_j - Lambda_j_1 || < delta = ', delta, '\n'))
    }
  }
  
  ### Group PCA space ###
  W_j = EVD.XChi$u[,1:k]
  if (whitenPCs){
    if (verbose){cat('projecting orthogonalized L2-normed data onto eigenvectors & whitening PCs by scaling to eigenvalues and s.d....')}
    X = X_j %*% W_j %*% diag( 1 / sqrt(Lambda_j) )   # eq (33)
    X = X %*% diag( 1 / apply(X,2,sd) )              # scale principal components, for comparison w/ other PCA methods
  }else{
    if (verbose){cat('projecting orthogonalized L2-normed data onto eigenvectors & de-whitening PCs accordingly...')}
    X = X_j %*% W_j
    X = X %*% diag(sqrt(Lambda_j * (v-1)))           # scale principal components, for comparison w/ other PCA methods
  }
  
  
  ### Eigenvectors & Eigenvalues ###
  if (return.eigenvectors){
    if (verbose){cat('finding eigenvectors by projecting PCs onto data...')}
    if (all(is.na(datafiles))){  datafiles = Y  }
    
    U = get_eigenValsVecs_fromPCs(datafiles, X, 
                                  var=var, var.subset=var.subset, mask=mask,
                                  unstacked=unstacked)$U
  }else{  'U' = NA  }
  
  
  #############################################
  return(list('X'           = X,          # group-level PCA space
              'U'           = U,          # group-level eigenvector basis
              'Lambda'      = Lambda_j,   # group-level eigenvalues
              'whitenedPCs' = whitenPCs,  # PCs scaled to eigenvalues
              'delta'       = delta,      # tolerence level
              'j'           = j))         # stopping iteration algorithm 
} #############################################



