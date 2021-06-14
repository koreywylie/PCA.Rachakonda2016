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
#   Y : temp.-concatenated group data,    [v by Mp] matrix
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
PCA.Rachakonda2016 <- function(Y, k=NA, 
                               var=NA, mask=NA,
                               whitenPCs=T,
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
  #   var  : if Y is an .RData file(s), name of saved var. w/ data
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs           : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #   algorithm           : algorithm PCA algorithm to call, fn. name as a string
  #                           see below for named options
  #   unstacked           : with input algorithm above, if True run as unstacked (low mem., parallelizable),
  #                                                   else if False run as stacked (high mem., fast)
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
    recommendations = find_PCAtype_best(Y, k=k, var=var, mask=mask,
                                        verbose=verbose)
    algorithm = recommendations$PCAtype
    unstacked = recommendations$unstacked
  }
  
  ##############################################################
  if (algorithm     == "PCA.EVDtime"){
    if (verbose){cat("  ...running PCA.EVDtime...")}
    
    PCA = PCA.EVDtime(Y, k=k, var=var, 
                      mask=mask,
                      whitenPCs=whitenPCs,
                      verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = NA
    
    ##############################################################
  }else if (algorithm     == "PCA.EVDvoxel"){
    if (verbose){cat("  ...running PCA.EVDvoxel...")}
    
    PCA = PCA.EVDvoxel(Y, k=k, 
                       var=var, mask=mask,
                       unstacked=unstacked, 
                       whitenPCs=whitenPCs,
                       return.eigenvectors=return.eigenvectors,
                       verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
    ##############################################################
  }else if (algorithm     == "PCA.SubsampledTime"){
    if (verbose){cat("  ...running PCA.SubsampledTime...")}
    
    PCA = PCA.SubsampledTime(Y, k=k, 
                             var=var, mask=mask,
                             unstacked=unstacked,
                             whitenPCs=whitenPCs,
                             return.eigenvectors=return.eigenvectors,
                             verbose=verbose)
    PCA$algorithm = algorithm
    PCA$unstacked = unstacked
    
    ##############################################################  
  }else if (algorithm     == "PCA.MultiPowerIteration"){
    if (verbose){cat("  ...running PCA.MultiPowerIteration...")}
    
    X0 = 'STP'  # initialize w/ Subsampled Time PCA, to speed convergence
    PCA = PCA.MultiPowerIteration(Y, k=k, 
                                  X0=X0, l=5, delta=0.00001, 
                                  var=var, mask=mask,
                                  unstacked=unstacked,
                                  whitenPCs=whitenPCs,
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
load.data <- function(datafiles, 
                      var='pca', mask=NA,
                      verbose=F){
  #############################################
  # Loads saved var. from datafiles, subj. by subj.,
  #   including 3d & 4d nifti volumes applying (optional) mask
  # 
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
    target_img = RNifti::asNifti(target_img[,,,1])
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
  
  #############################################
  resampled = RNiftyReg::niftyreg.linear(source_img, 
                                         target_img,
                                         init=diag(4), nLevels=0)
  #############################################
  return(resampled$image)
} #############################################

#############################################
flatten_img <- function(img, mask=NA){
  #############################################
  # Flattens 4d nifti vol., loaded w/ readNifti()
  #   applying optional mask,
  #     returns matrix as [voxels x time]
  #
  # Reversed by:
  #   unflatten_img()
  # Requires: 
  #   resample_img() if dims. of input differ
  #
  
  if (is.character(img) && file.exists(img)){
    filename = img
    img = RNifti::readNifti(img)
  }
  
  if (!all(is.na(mask))){
    if (is.character(mask) && file.exists(mask)){
      mask_filename = mask
      mask = RNifti::readNifti(mask)
    }
    if (!all(dim(img)[1:3] == dim(mask))){
      mask = resample_img(mask, img)
    }
    stopifnot(all(dim(img)[1:3] == dim(mask)))
    stopifnot(all(RNifti::pixdim(img)[1:3] == RNifti::pixdim(mask)))
    stopifnot(all(RNifti::origin(img) == RNifti::origin(mask)))
    
    if (any(is.na(mask))){  mask[is.na(mask)] = 0   }
    mask = mask != 0
  }else{
    mask = rep(T, prod(dim(img)[1:3], na.rm=T))
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
get_data_dims <- function(Y, var=NA, mask=NA,
                          unstacked=T, return.data=F, verbose=F){
  #############################################
  # Find dimensions of data, 
  #   loading 1st subj. if needed
  #   finding saved 'var' or applying 'mask' if indicated
  #
  # OUTPUT:
  #   v : length of 1st dim. ~ number of voxels
  #   S : number of subjs. (2nd dim. is S*p)
  #   p : number of time points/ICs (~2nd dim.)
  #   datafiles : paths to saved files
  #   Y : (optional) if indicated by return.data=T,
  #         loaded data for all subjs.
  #
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
      Y = load.data(datafiles, var, mask, verbose=verbose)
      spacetime = get.spacetime(Y[[1]])
    }else{
      Y = load.data(datafiles[1], var, mask, verbose=verbose)
      spacetime = get.spacetime(Y)
    }
    
  }else if (is.character(Y) && all(file.exists(Y))){
    if (stacked){  # load all subjs. into memory at once
      datafiles = paste(basename(Y), collapse=', ')
      if (verbose){cat(paste0('Loading & stacking: ', datafiles), '...\n')}
      Y = load.data(Y, var, mask, verbose=verbose)
      S = length(datafiles)
      spacetime = get.spacetime(Y[[1]])
      
    }else{            # load subjs. one-by-one
      datafiles = Y
      S = length(datafiles)
      Y = load.data(datafiles[1], var, mask)
      spacetime = get.spacetime(Y)
      
    }
  }else{
    if (is.character(Y) && !all(file.exists(Y))){
      cat(paste0('WARNING: could not find file as input:'))
      print(Y)
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
  return(list('v'=v, 'S'=S, 'p'=p, 
              'datafiles'=datafiles,
              'Y'=Y))
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
}

#############################################
zeroMean_Yi <- function(Y){
  #############################################
  # Removes column mean of subj.-level data Y
  #     ~centering every vol. at every timepoint
  #   where Y is either list of Yi   [v by p] matrices
  #   or Y is array with dim.   [S by v by p]
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
get_eigenValsVecs_fromPCs <- function(Y, X,
                                      var=NA, mask=NA,
                                      unstacked=F, verbose=F){
  #############################################
  # Find eigenvalues & eigenvectors from data Y & (optional) variable var,
  #   given input Principal Components X.
  # Used to get eigenvalues
  #   following PCA methods that only calc. PCs
  #
  # Theory: 
  #   If Y is the [v x Mp] data matrix,
  #     & X is the [v x k] group PCs,
  #     solve for the [Mp x k] matrix U of eigenvectors,
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
  
  data.dims = get_data_dims(Y, var, mask,
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
                           
                           Yi = load.data(sub, var, mask,
                                          verbose=verbose)
                           Yi = zeroMean_Yi(Yi)
                           return( t(Yi) %*% X )
                           
                         }
  }else{
    Y = load.data(Y, var, mask)
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
find_PCAtype_best <- function(data, k=NA, 
                              var=NA, mask=NA,
                              verbose=F){
  #############################################
  # Finds best PCA algorithm,
  #   based on recommendations in Figure 8, Rachakonda et al. (2016)
  #
  
  ### Get data dims ###
  data.dims = get_data_dims(data, var, mask)
  S = data.dims$S
  v = data.dims$v
  p = data.dims$p
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
    if (v*S*p > memFree){
      PCAtype = 'PCA.MultiPowerIteration'
      unstacked = T
    }else{
      PCAtype = 'PCA.MultiPowerIteration'
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
                                     var=var, mask=mask)
  if (length(memInfo) > 1){
    memNeeded = min(memInfo[[1]]$mem, memInfo[[2]]$mem, na.rm=T)
  }else{
    memNeeded = memInfo[[1]]$mem
  }
  
  if (verbose){
    cat(paste0('\n      requiring  ',signif(memNeeded/1024),' Mb  of Memory'))
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
calc_PCA_memRequirements <- function(data, k, pca_type='all', 
                                     var=NA, mask=NA,
                                     verbose=F){
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
  #   var  : if Y is file(s), name of saved var. w/ data
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
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
    pca_types = c('PCA.EVDtime', 
                  'PCA.EVDvoxel',
                  'PCA.SubsampledTime',
                  'PCA.MultiPowerIteration')
  }else{  pca_types = pca_type  }
  bytes = 8 # double precision, numeric data
  if (!exists('l')){  l = 5   }  # default multiplier for MPOWIT
  if (!exists('g')){  g = 20  }  # default number of subjs. / group in STP
  if (!exists('k_prime')){
    k_prime = max(k*5, 500) # dim. of PCA during group update, defaults suggested number in paper
  }
  
  ### Get data dims ###
  data.dims = get_data_dims(data, var, mask)
  S = data.dims$S
  v = data.dims$v
  p = data.dims$p

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
      g1 = min(g, S)
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
      mem_stacked = v * S * p       # total data size stacked...
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
PCA.EVDtime <- function(Y, k=NA, 
                        var=NA, mask=NA,
                        whitenPCs=T,
                        verbose=F){
  #############################################
  # Principal Components Analysis via Eigenvalue Decomposition (well, actually ~SVD)
  #   to find subj. or group-level PCA space
  #
  # INPUT:
  #   Y : subjs. or group-level data,
  #         as filename(s)
  #         as temp.-concatenated data,    [v by Mp] matrix
  #         as list of [v by p] matrices
  #         or as array w/ subjs. indices in 1st dim.
  #   k : number of group-level top components, after group-level PCA
  #   var  : if Y is file(s), name of saved var. w/ data
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  # OUTPUT:
  #   X : group-level PCA space,            [v by k] matrix
  #   U : group-level eigenvector basis     [Mp by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # REQUIRES:
  #   load.data()         : data munging
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #
  
  if (verbose){cat("\nCalculating PCA space using EVD on time*subjs. dimension...")}

  if (is.character(Y) && all(file.exists(Y))){
    Y = load.data(Y, var, mask) #loads all subjs. into mem. at once
  }else if (is.list(Y)){
    stopifnot(is.numeric(Y[[1]]))
  }else{  stopifnot(is.numeric(Y))  }
  
  Y = zeroMean_Yi(Y)
  Y = concatTemporal_Yi(Y)
  v = nrow(Y)
  Mp = ncol(Y)
  k = min(k, v, Mp, na.rm=T)
  
  C = (t(Y) %*% Y) / (v-1)  # eq. (2)
  
  EVD = svd(C, k)
  
  if (whitenPCs){
    if (verbose){cat('\n...projecting data onto eigenvectors & whitening by scaling to eigenvalues...')}
    X = Y %*% EVD$u %*% diag( 1 / sqrt(EVD$d[1:k])) # eq. (4)
  }else{
    if (verbose){cat('\n...projecting data onto eigenvectors...')}
    X = Y %*% EVD$u
  }
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'U'      = EVD$u,       # Eigenvectors
              'Lambda' = EVD$d[1:k],  # Eigenvalues
              'whitenedPCs' = whitenPCs))
} #############################################


#############################################
PCA.EVDvoxel <- function(Y, k=NA, 
                         var=NA, mask=NA,
                         whitenPCs=T,
                         unstacked=F, 
                         return.eigenvectors=F,
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
  #   var  : if Y is file(s), name of saved var. w/ data
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #                       (NOTE: requires parallelization cluster set up)
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # OUTPUT:
  #   X   : group-level PCA space,          [v by k] matrix
  #   U  : group-level eigenvector basis   [Mp by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # REQUIRES:
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.data()         :  "     "
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
  
  data.dims = get_data_dims(Y, var, mask)
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
                   
            Yi = load.data(sub, var, mask, verbose=verbose)
            Yi = zeroMean_Yi(Yi)
            return((Yi %*% t(Yi)) / (v-1))  # eq. (3)
            
          }
  }else{
    Y = load.data(Y, var, mask)
    Y = zeroMean_Yi(Y)
    Y = concatTemporal_Yi(Y)
    Cv = (Y %*% t(Y)) / (v-1)  # eq. (3)
  }
  
  ### Principal Components ###
  EVD = svd(Cv, k)
  X = EVD$u * sqrt(v-1)  # eq. (4) left
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
                                  var=var, mask=mask,
                                  unstacked=unstacked)$U
    
  }else{  'U' = NA  }
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'U'      = U,           # Eigenvectors
              'Lambda' = Lambda,      # Eigenvalues
              'whitenedPCs' = whitenPCs))
} #############################################



#############################################
PCA.SubsampledTime <- function(Y, k=NA, g=20, parFactor=4,
                               var=NA, mask=NA,
                               whitenPCs=T,
                               unstacked=F, 
                               return.eigenvectors=F,
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
  #   var  : if Y is file(s), name of saved var. w/ data
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # OUTPUT:
  #   X : group-level PCA space,            [v by k] matrix
  #   U : group-level eigenvector basis    [Mp by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #
  # REQUIRES:
  #   PCA.EVDtime()
  #   check_stacked()     : data munging
  #   get_data_dims()     :  "     "
  #   load.data()         :  "     "
  #   zeroMean_Yi()       :  "     "
  #   concatTemporal_Yi() :  "     "
  #  
  # SPEED:   time complexity  O(v(Mp)^2)
  #
  
  if (verbose){cat("\nCalculating PCA space by subsampling time*subjs. dimension...")}
  
  ### Required fns. ###
  partition_group <- function(subj.count, group.size){ 
    # fn. to randomly partitions subjs. into groups, with set group size
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
  
  calc_GroupSubSampleTimeSVD <- function(data, partition_inds, k_prime,  
                                         var=NA, mask=NA, Xg_0=NA){
    ### Calc. for C_g ###
    if (is.character(data)){
      Yg = load.data(data[unlist(partition_inds)], var, mask)
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
    EVD.Yg = PCA.EVDtime(Yg, k_prime)               # eq. (11) & (12) in Rachakonda et al. 2016
    Xg = Yg %*% EVD.Yg$U                            # eq. (13) in Rachakonda et al. 2016
    
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
      
      EVD.g0g1 = svd(Cg0g1, k_prime)                 # eq. (14) right in Rachakonda et al. 2016
      Xg = cbind(Xg_0, Xg) %*% EVD.g0g1$u            # eq. (15)  in Rachakonda et al. 2016
    }
    return(Xg)
  } ###################
  
  
  ### Defaults & Inputs ###
  stacked = check_stacked(Y) || !as.logical(unstacked) # Check requested arg. inputs vs. reality of data input Y
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  
  data.dims = get_data_dims(Y, var, mask, 
                            unstacked, return.data=stacked)
  S = data.dims$S
  v = data.dims$v
  p = data.dims$p
  datafiles = data.dims$datafiles
  if (stacked){  Y = data.dims$Y  }
  stopifnot(exists('S')); stopifnot(exists('v')); stopifnot(exists('p'))
  k = min(v, S*p, k, na.rm=T)
  k_prime = min(max(k*5, 500), # dim. of PCA during group update, defaults suggested number in paper
                v, S*p)        #  limited by rank of matrix
  
  
  ### Subdivide by groups & EVD on each ###
  if (unstacked){
    g_node = ceiling(S / min(S, parFactor))      # number of subjs. on each parallel node
    node.part.list = partition_group(S, g_node)
    numNodes = length(node.part.list)
    Xg = foreach(n.pt.i=icount(numNodes),
                 .export=c('load.data', 'zeroMean_Yi', 'concatTemporal_Yi',
                           'PCA.EVDtime'), 
                 .combine='cbind') %dopar% {
                   dataSubset = datafiles[unlist(node.part.list[[n.pt.i]])]
                   grp.part.list = partition_group(length(dataSubset), g)
                   Xg_0 = Xg = NA
                   for (pt.i in 1:(length(grp.part.list))){
                     Xg_0 = Xg #re-name prev. iteration
                     Xg = calc_GroupSubSampleTimeSVD(dataSubset, 
                                                     grp.part.list[[pt.i]], 
                                                     k_prime, var, mask, Xg_0)
                   }
                   return(Xg)
                 }
    EVD.g_final = svd((t(Xg) %*% Xg) / (v-1)) # combine individual subgroup results
    Xg = Xg %*% EVD.g_final$u
  }else{
    grp.part.list = partition_group(S, g)
    Xg_0 = Xg = NA
    for (pt.i in 1:(length(grp.part.list))){
      Xg_0 = Xg  #re-name prev. iteration
      Xg = calc_GroupSubSampleTimeSVD(Y, 
                                      grp.part.list[[pt.i]],
                                      k_prime, var, mask, Xg_0)
    }
    EVD.g_final = svd(t(Xg) %*% Xg / (v-1))
  }
  
  ### Group PCA space ###
  if (whitenPCs){
    if (verbose){cat('projecting data onto eigenvectors & whitening by scaling to eigenvalues...')}
    X = Xg[,1:k] %*% diag( 1 / sqrt(EVD.g_final$d[1:k])) # Rachakonda et al. (2016) final paragraph for STP section
  }else{
    X = Xg[,1:k]
  }
  
  ### Eigenvectors & Eigenvalues ###
  Lambda = EVD.g_final$d[1:k]
  if (return.eigenvectors){
    if (verbose){cat('finding eigenvectors by projecting PCs onto data...')}
    if (all(is.na(datafiles))){  datafiles = Y  }
    
    U = get_eigenValsVecs_fromPCs(datafiles, X, 
                                  var=var, mask=mask,
                                  unstacked=unstacked)$U
    
  }else{  'U' = NA  }
  
  #############################################
  return(list('X'      = X,           # Principal Comps.
              'U'      = U,           # Eigenvectors
              'Lambda' = Lambda,      # Eigenvalues
              'whitenedPCs' = whitenPCs))
} #############################################



#############################################
PCA.MultiPowerIteration <- function(Y, k=NA, 
                                    X0=NA, l=5, delta=0.001,
                                    var=NA, mask=NA,
                                    whitenPCs=T,
                                    unstacked=F,
                                    return.eigenvectors=F,
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
  #   var  : if Y is file(s), name of saved var. w/ data
  #   mask : if Y is a Nifti file(s), optional binary mask used to limit analysis to contained voxels
  #   whitenPCs : return whitened PCs, by scaling eigenvectors to sqrt. of eigenvalues
  #   unstacked : if Y is a vector of file paths, load & calc. cov. matrix individually
  #                 in order to minimize memory use, despite very large number of subjs.
  #               If unstacked=T, will run in parallel using R doParallel library
  #   return.eigenvectors : calc. & return eigenvectors from group PCA & original data
  #
  # OUTPUT:
  #   X : group-level PCA space,            [v by k] matrix
  #   U : group-level eigenvector basis     [Mp by k] matrix
  #   Lambda : group-level eigenvalues      [k by k] diagonal matrix
  #   delta  : tolerence level
  #
  # REQUIRES:
  #   zeroMean_Yi()       : data munging
  #   concatTemporal_Yi() : data munging for arrays & lists
  #   get_data_dims()     : find number of subjs., voxels, time points, etc.
  #   spatialSubsample()  : sampling, taking into account spatial structure of data
  #   PCA.SubsampledTime() : used to initialize PCs to speed convergence
  #
  # EQUATIONS:
  #   X_0 = G_v,lk,  Chi_0 = (YY^T)X_0,  Lambda_0 = 0          eq (28)
  #   X_j = ortho(Chi_j_1) = Chi_j_1 %*% U %*% L^-1 >= 1       eq (29)
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
    U = svd(t(X) %*% X)$u
    L = apply(X %*% U, 2, L2.norm) # get L2 norms of cols.
    return(X %*% U %*% diag( 1 / L))
  }
  
  ### Defaults & Input ###
  stacked = check_stacked(Y) || !as.logical(unstacked)
  unstacked = !stacked    # if data is already loaded in memory, ignore irrelevant input
  if (unstacked){  library(foreach); library(iterators)  }
  
  l = as.integer(l)
  stopifnot(l > 0)
  
  data.dims = get_data_dims(Y, var, mask, 
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
    X_j = matrix(rnorm(v*lk), v, lk)            # eq (28) 
    Lambda_j = matrix(0, k, k)                  # eq (28) 
  }else if (is.character(X0)){
    stopifnot(X0 %in% c('STP'))
    if (X0 == 'STP'){
      if (verbose){cat('\n...initializing MPOWIT w/ Subsampled Time PCA...')}
      pca.y = PCA.SubsampledTime(Y, mask=mask, g=4, unstacked=unstacked)
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
                    .export=c('load.data',
                              'zeroMean_Yi',
                              'concatTemporal_Yi'),
                    .combine='+') %dopar% {
                      Yi = load.data(sub, var, mask, verbose=verbose)
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
                      .export=c('load.data',
                                'zeroMean_Yi'), 
                      .combine='+') %dopar% {
                        Yi = load.data(sub, var, mask)
                        Yi = zeroMean_Yi(Yi)
                        return(Yi %*% t(t(X_j) %*% Yi))   # eq (30) right hand side
                      }
    }else{
      Y = zeroMean_Yi(Y)
      Y = concatTemporal_Yi(Y)
      Chi_j = Y %*% t(Y) %*% X_j                # eq (30)
    }
    
    
    EVD.XChi = svd((t(X_j) %*% Chi_j) / (v-1))  # eq (31)
    Lambda_j = EVD.XChi$d[1:k]
    
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
  W_j = EVD.XChi$u[,1:k]
  if (whitenPCs){   # whitening PCs recommended for comparison w/ other PCA methods
    if (verbose){cat('projecting data onto eigenvectors & whitening by scaling to eigenvalues and s.d....')}
    X = X_j %*% W_j %*% diag( 1 / sqrt(Lambda_j))   # eq (33)
    X = X %*% diag( 1 / apply(X,2,sd))              # scale principal components
  }else{
    X = X_j %*% W_j
  }
  

  ### Eigenvectors & Eigenvalues ###
  if (return.eigenvectors){
    if (verbose){cat('finding eigenvectors by projecting PCs onto data...')}
    if (all(is.na(datafiles))){  datafiles = Y  }
    U = get_eigenValsVecs_fromPCs(datafiles, X, 
                                  var=var, mask=mask,
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



