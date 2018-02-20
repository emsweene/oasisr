#' @title OASIS Training Data Frame
#' @description This function creates the training vectors from a single MRI
#' study that has FLAIR, T1, T2, and PD volumes
#' as well as binary masks of lesions. The function can create a brain mask
#' for the data (or the user can supply a brain mask),
#' can preprocess the data, and the user may supply already normalized data
#' if they wish to use an alternative normalization method.
#' @param flair FLAIR volume of class \code{\link{nifti}}
#' @param t1 T1 volume of class \code{\link{nifti}}
#' @param t2 T2 volume of class \code{\link{nifti}}
#' @param pd PD volume of class \code{\link{nifti}}
#' @param gold_standard gold standard lesion segmentation mask of class
#' \code{\link{nifti}}
#' @param brain_mask brain mask of class \code{\link{nifti}}, if NULL a
#' brain mask will be created using \code{\link{fslbet}}
#' @param preproc is a logical value that determines whether to call the
#' \code{\link{oasis_preproc}} function
#' and performs the necessary preprocessing steps for OASIS
#' @param normalize is a logical value that determines whether
#' to perform z-score normalization of the image over the brain mask,
#' should be \code{TRUE} unless you train model
#' using an alternative normalization
#' @param slices vector of desired slices to train on, if \code{NULL}
#' then train over the entire brain mask
#' @param orientation string value telling which orientation the
#' training slices are specified in, can take the values of "axial",
#' "sagittal", or "coronal"
#' @param return_preproc is a logical value that indicates whether
#' the preprocessed images should be returned
#' @param cores numeric indicating the number of cores to be used
#' (no more than 4 is useful for this software implementation)
#' @param sigma Sigmas used to smooth the data, default is 10,20
#' @param verbose print diagnostic output
#' @param eroder Should \code{\link{fslerode}} or \code{\link{oasis_erode}} 
#' be used
#'
#' @return If \code{return_preproc = FALSE} the function returns a
#' \code{data.frame} for use with the \code{\link{oasis_training}} function.
#' Otherwise, the function returns a list containing:
#' a \code{data.frame} for use with the \code{\link{oasis_training}} 
#' function,
#' the FLAIR volume, the T1 volume, the T2 volume,
#' the PD volume, the brain mask for the subject, and the voxel selection 
#' mask.
#' @seealso \code{\link{oasis_training}}
#' @export
#' @importFrom neurobase zscore_img
#' @importFrom fslr fslbet fslerode
#' @examples 
#' library(neurobase)
#' dl_file = function(url) {
#'    tfile = tempfile(fileext = ".nii.gz")
#'    req <- httr::GET(url,
#'    httr::write_disk(path = tfile))
#'    httr::stop_for_status(req)
#'    tfile
#' }
#' in_ci <- function() {
#'  nzchar(Sys.getenv("CI"))
#' }
#' on_cran = function() {
#'  identical(Sys.getenv("NOT_CRAN"), "false")
#' } 
#' if (in_ci() || on_cran()) {
#'   if (fslr::have.fsl() && require(httr)) {
#'     mods = c("FLAIR", "T1W", "T2W", "consensus_gt", "brainmask")
#'     base_url = file.path(
#'       "https://raw.githubusercontent.com/muschellij2/open_ms_data", 
#'       "master/cross_sectional/coregistered/patient01/")
#'     files = paste0(base_url, mods, ".nii.gz")
#'     files = sapply(files, dl_file)
#'     names(files) = mods
#' 
#'     flair <- readnii(files["FLAIR"])
#'     t1 <- readnii(files["T1W"])
#'     t2 <- readnii(files["T2W"])
#'     brain_mask <- readnii(files["brainmask"])
#'     gold_standard = readnii(files["consensus_gt"])
#'     oasis_preprocessed_data <- oasis_train_dataframe(flair, t1, t2, 
#'       brain_mask = brain_mask, gold_standard = gold_standard)
#'   } 
#' }
oasis_train_dataframe <- function(
  flair, ##flair volume of class nifti
  t1, ##t1 volume of class nifti
  t2, ##t2 volume of class nifti
  pd = NULL, ##pd volume of class nifti
  gold_standard = NULL, ##gold standard mask of class nifti
  brain_mask = NULL, ##brain mask of class nifti
  preproc = FALSE, ##option to preprocess the data
  normalize = TRUE, ##option to normalize
  slices = NULL, #slice vector
  orientation = c("axial", "coronal", "sagittal"),
  return_preproc = FALSE,
  cores = 1,
  sigma = c(10, 20),
  verbose = TRUE,
  eroder = c("fsl", "oasis")
)
{
  if (verbose) {
    message("Checking File inputs")
  }
  check_nifti2 = function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      return(check_nifti(x))
    }
  }
  flair = check_nifti2(flair)
  t1 = check_nifti2(t1)
  t2 = check_nifti2(t2)
  pd = check_nifti2(pd)
  
  correct_image_dim2 = function(x) {
    if (is.null(x)) {
      return(NULL)
    } else {
      return(correct_image_dim(x))
    }
  }
  ##correct image dimmension
  flair <- correct_image_dim2(flair)
  t1 <- correct_image_dim2(t1)
  t2 <- correct_image_dim2(t2)
  pd <- correct_image_dim2(pd)
  
  ##image preproceesing
  if (preproc == TRUE) {
    if (verbose) {
      message("OASIS Preprocessing")
    }
    ## the image preproceesing
    preprocess <- oasis_preproc(flair = flair,
                                t1 = t1,
                                t2 = t2,
                                pd = pd,
                                cores = cores,
                                brain_mask = brain_mask,
                                verbose = verbose)
    oasis_study <- preprocess[c("flair","t1", "t2", "pd")]
    brain_mask <- preprocess$brain_mask
  } else {
    ## no preprocessing
    oasis_study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  }
  # REMOVE NULL
  nulls = sapply(oasis_study, is.null)
  oasis_study = oasis_study[!nulls]
  
  ###############################
  # Making brain mask if one not needed
  ###############################
  if (is.null(brain_mask) & !preproc) {
    if (verbose) {
      message("Getting Brain Mask")
    }
    ## create a brain mask if not supplied
    brain_mask <- fslbet(infile = oasis_study$t1, retimg = TRUE)
  }
  
  brain_mask = check_nifti(brain_mask)
  brain_mask <- brain_mask > 0
  brain_mask <- datatyper(brain_mask, trybyte = TRUE)
  
  
  ##adjust brain mask for OASIS
  brain_mask <- correct_image_dim(brain_mask)
  if (verbose) {
    message("Eroding Brain Mask")
  }
  eroder = match.arg(eroder)
  if (eroder == "fsl") {
    ero_brain_mask <- fslerode(brain_mask,
                               kopts = "-kernel box 5x5x5",
                               retimg = TRUE)
  }
  if (eroder == "oasis") {
    ero_brain_mask <- oasis_erode(mask = brain_mask,
                                  mm = c(5,5,5))
  }
  
  ##removing voxels below 15th quantile
  brain_mask <- voxel_selection(flair = oasis_study$flair,
                                brain_mask = ero_brain_mask,
                                cutoff = .15)
  
  # cutpoint <- quantile(oasis_study$flair[brain_mask == 1], probs = .15)
  # brain_mask[oasis_study$flair <= cutpoint] <- 0
  
  
  
  ## the image normalization
  if (normalize == TRUE) {
    if (verbose) {
      message("Normalizing Images using Z-score")
    }
    oasis_study <- lapply(oasis_study, zscore_img,
                          mask = brain_mask,
                          margin = NULL)
  }
  
  
  if (verbose) {
    message("Voxel Selection Procedure")
  }
  ##create and apply the voxel selection mask
  top_voxels <- voxel_selection(flair = oasis_study$flair,
                                brain_mask = brain_mask,
                                cutoff = .85)
  
  
  orig_study = oasis_study
  cnames = toupper(names(orig_study))
  names(oasis_study) = cnames
  
  if (length(sigma) > 0) {
    for (isigma in sigma) {
      if (verbose) {
        message(paste0("Smoothing Images: Sigma = ", isigma))
      }
      
      
      # New 2017May03
      smoothed_mask = fslr::fsl_smooth(
        file = brain_mask,
        sigma = isigma,
        smooth_mask = FALSE)
      
      ## smooth the images using fslsmooth from the fslr package
      smooth <- mclapply(orig_study, fslsmooth,
                         sigma = isigma,
                         mask = brain_mask,
                         retimg = TRUE,
                         smooth_mask = TRUE,
                         smoothed_mask = smoothed_mask,
                         mc.cores = cores)
      names(smooth) = paste0(cnames, "_", isigma)
      oasis_study = c(oasis_study, smooth)
      rm(list = c("smooth"))
    }
  }
  
  gold_standard = check_nifti2(gold_standard)
  gold_standard = correct_image_dim2(gold_standard)
  
  oasis_study$GoldStandard <- gold_standard
  
  oasis_study$top_voxels <- top_voxels
  
  #######################################
  # Make data.frame
  #######################################
  oasis_dataframe = lapply(oasis_study, c)
  oasis_dataframe = as.data.frame(oasis_dataframe)
  rownames(oasis_dataframe) = NULL
  
  ######################
  # Adding in slice indicators
  ######################
  inds = niftiarr(brain_mask, 1)
  inds = which(inds == 1, arr.ind = TRUE)
  orientations = c("axial", "coronal", "sagittal")
  colnames(inds) = orientations
  oasis_dataframe = cbind(oasis_dataframe, inds)
  
  ######################
  # Keeping Voxel Selection
  ######################
  oasis_dataframe = oasis_dataframe[ oasis_dataframe$top_voxels == 1, ]
  oasis_dataframe$top_voxels = NULL
  
  ######################
  # Keep only the slices indicated if slices aren't NULL
  ######################
  if (!is.null(slices)) {
    orientation = match.arg(orientation)
    
    oasis_dataframe = oasis_dataframe[
      oasis_dataframe[, orientation] %in% slices, ]
  }
  cn = colnames(oasis_dataframe)
  cn = setdiff(cn, orientations)
  oasis_dataframe = oasis_dataframe[, cn]
  
  ######################
  # Return Preproc
  ######################
  if (!return_preproc) {
    orig_study = NULL
  }
  L = list(oasis_dataframe = oasis_dataframe,
           brain_mask = brain_mask,
           voxel_selection = top_voxels,
           ero_brain_mask = ero_brain_mask)
  L$preproc = orig_study
  
  return(L)
}