#' @title OASIS Prediction
#' @description This function creates the OASIS probability map from a single 
#' MRI study with FLAIR, T1, T2, and PD volumes.
#' @param flair flair volume of class \code{\link{nifti}}
#' @param t1 t1 volume of class \code{\link{nifti}}
#' @param t2 t2 volume of class \code{\link{nifti}}
#' @param pd pd volume of class \code{\link{nifti}}
#' @param brain_mask brain mask of class \code{\link{nifti}},
#' if \code{NULL} a brain mask will be created using \code{\link{fslbet}}.
#' Note that provided brain masks should be in the same space as the T1 volume
#' if \code{preproc = TRUE}, as all volumes will be registered to this space
#' @param model an object of class \code{\link{glm}} used to make the 
#' OASIS predictions
#' @param return_preproc is a logical value that indicates whether the
#' preprocessed images should be returned, if \code{NULL}
#' then the model from the OASIS paper will be used
#' @param binary logical indicating whether a binary map
#' should be returned by thresholding the probability map
#' @param threshold numeric indicating the threshold value
#' for the probability map, with default of 0.16 for the OASIS paper
#' @param verbose print diagnostic messages
#' @param oasis_dataframe if \code{\link{oasis_train_dataframe}}
#' was already run, specify the \code{data.frame} and 
#' \code{voxel_selection} and \code{brain_mask} to make prediction
#' @param voxel_selection image of selected voxels.  
#' If \code{\link{oasis_train_dataframe}}
#' was already run, specify the \code{data.frame} and 
#' \code{voxel_selection} and \code{brain_mask} to make prediction
#' 
#' @param ... options passed to \code{\link{oasis_train_dataframe}}
#'
#' @return A list of volumes:
#' the OASIS probability map, the preprocessed volumes (if 
#' \code{return_preproc = TRUE}),
#' the brain mask for the subject, the voxel selection mask, and a thresholded,
#' binary mask (if \code{binary = TRUE}) .
#' @examples 
#' library(ROCR)
#' p = predict( oasis::oasis_model,
#'     newdata = example_oasis_df,
#'     type = 'response')
#' nopd_p = predict( oasis::nopd_oasis_model,
#'     newdata = example_oasis_df,
#'     type = 'response')    
#' y =  example_oasis_df$GOLD_Lesions
#' pred = ROCR::prediction(p, y)
#' perf = ROCR::performance(pred, "tpr", "fpr")
#' plot(perf)
#' 
#' 
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
#'     oasis_preprocessed_data <- oasis_predict(flair, t1, t2, 
#'       brain_mask = brain_mask, preproc = TRUE)
#'   } 
#' }
#' @export
#' @importFrom stats predict
#' @importFrom fslr fslsmooth fsl_smooth
#' @importFrom neurobase datatyper
#' @importFrom oro.nifti convert.bitpix convert.datatype

oasis_predict <- function(flair, ##flair volume of class nifti
                          t1, ##t1 volume of class nifti
                          t2, ##t2 volume of class nifti
                          pd = NULL, ##pd volume of class nifti
                          brain_mask = NULL, ##brain mask of class nifti
                          model = NULL, ##an OASIS model of class glm
                          return_preproc = FALSE, 
                          ##option to return the preprocessed data
                          binary = FALSE,
                          threshold = 0.16,
                          verbose = TRUE,
                          oasis_dataframe = NULL,
                          voxel_selection = NULL,
                          ...
) {
  if (!is.null(oasis_dataframe) && !is.null(voxel_selection)
      && !is.null(brain_mask)) {
    if (return_preproc) {
      warning(paste0("No preprocessing was done as the ", 
                     "voxel_selection, brain mask, and data frame were", 
                     " specified"))
      return_preproc = FALSE
    }
    preproc = list()
  } else {
    L = oasis_train_dataframe(flair = flair,
                              t1 = t1,
                              t2 = t2,
                              pd = pd,
                              brain_mask = brain_mask,
                              verbose = verbose,
                              return_preproc = TRUE,
                              ...)
    
    oasis_dataframe = L$oasis_dataframe
    brain_mask = L$brain_mask
    voxel_selection = L$voxel_selection
    preproc = L$preproc
    rm(list = "L")
  }
  
  
  if (verbose) {
    message("Model Prediction")
  }
  ## make the model predictions
  if (is.null(model)) {
    predictions <- predict( oasis::oasis_model,
                            newdata = oasis_dataframe,
                            type = 'response')
  } else {
    predictions <- predict(model,
                           newdata = oasis_dataframe,
                           type = 'response')
  }
  
  
  ##put the predictions onto the brain
  predictions_nifti <- niftiarr(brain_mask, 0)
  predictions_nifti[voxel_selection == 1] <- predictions
  predictions_nifti = datatyper(predictions_nifti, 
                                datatype = convert.datatype()$FLOAT32,
                                bitpix = convert.bitpix()$FLOAT32
  )
  if (verbose) {
    message("Smoothing Prediction")
  }
  ##smooth the probability map
  prob_map <- fslsmooth(predictions_nifti, sigma = 1.25,
                        mask = brain_mask, retimg = TRUE,
                        smooth_mask = TRUE)
  
  binary_map = NULL
  if (binary == TRUE) {
    if (verbose) {
      message("Thresholding Smoothed Prediction")
    }
    binary_map <- prob_map
    binary_map[prob_map > threshold] <- 1
    binary_map[prob_map <= threshold] <- 0
  }
  L = list(oasis_map = prob_map,
           binary_map = binary_map,
           unsmoothed_map = predictions_nifti,
           flair = preproc$flair,
           t1 = preproc$t1, t2 = preproc$t2,
           pd = preproc$pd,
           brain_mask = brain_mask,
           voxel_selection = voxel_selection)
  L$threshold = threshold
  
  if (!return_preproc) {
    L$flair = L$t1 = L$t2 = L$pd = NULL
  }
  
  if (!binary) {
    L$binary_map = NULL
    L$threshold = NULL
  }
  
  return(L)
  
}


