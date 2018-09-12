#' @title OASIS Image Preprocessing
#' @description This function does the required preprocessing for OASIS 
#' for the FLAIR, T2,
#' T1, and PD volumes using FSL through \code{fslr}.
#' The preprocessing steps are
#' (1) inhomogeneity correct using \code{\link{fsl_biascorrect}}
#' and (2) rigid registration using \code{\link{flirt}} to the T1 space.
#' @param flair FLAIR volume of class \code{\link{nifti}}
#' @param t1 T1 volume of class \code{\link{nifti}}
#' @param t2 T2 volume of class \code{\link{nifti}}
#' @param pd PD volume of class \code{\link{nifti}}
#' @param brain_mask binary mask volume of class \code{\link{nifti}}
#' @param verbose a logical value for printing diagnostic output
#' @param bias_correct should the images be bias-field corrected?
#' @param cores numeric indicating the number of cores to be used 
#' (no more than 4 is useful for this software implementation)
#' 
#' @return Returns a list of objects of class \code{\link{nifti}},
#' namely the inhomogeneity corrected FLAIR, T1, T2, and PD registered to the
#' space of the T1 volume.
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
#'     oasis_preprocessed_data <- oasis_preproc(flair, t1, t2, 
#'       brain_mask = brain_mask)
#'   } 
#' }
#' @export
#' @importFrom fslr flirt fslbet fsl_biascorrect
#' @importFrom parallel mclapply
#' @importFrom neurobase check_nifti datatyper mask_img
oasis_preproc <- function(flair, #flair volume of class nifti
                          t1, # t1 volume of class nifti
                          t2 = NULL, # t2 volume of class nifti
                          pd = NULL, # pd volume of class nifti
                          brain_mask = NULL,
                          bias_correct = TRUE,
                          verbose = TRUE,
                          cores = 1
){
  
  
  study <- list(flair = flair, t1 = t1, t2 = t2, pd = pd)
  
  # REMOVE NULL
  nulls = sapply(study, is.null)
  study = study[!nulls]
  
  if (verbose) {
    message("Rigidly Registering Data to T1 Space\n")
  }
  
  seqs = c("flair","t2", "pd")
  seqs = intersect(names(study), seqs)
  
  ##rigidly register to the flair, t2, and pd to the t1 using fsl flirt
  study[seqs] <- mclapply(
    study[seqs], function(x)  {
      tfile = tempfile(fileext = ".mat")
      flirt(infile = x, omat = tfile,
            reffile = study$t1,
            retimg = TRUE,  dof = 6)
    }, mc.cores = cores)
  
  if (verbose) {
    message("Running Brain Extraction Tool\n")
  }
  if (is.null(brain_mask)) {
    brain_mask <- fslbet(infile = study$t1, retimg = TRUE)
  }
  brain_mask = check_nifti(brain_mask)
  brain_mask <- brain_mask > 0
  brain_mask <- datatyper(brain_mask, trybyte = TRUE)
  
  study = check_nifti(study)
  study <- mclapply(study, mask_img, mask = brain_mask, 
                    mc.cores = cores)
  
  ## inhomogeneity correction for all four modalities using fsl bias correct
  if (bias_correct) {
    if (verbose) {
      message("# Running Inhomogeneity Correction\n")
    }
    study <- mclapply(study, fsl_biascorrect,
                      retimg = TRUE,
                      verbose = verbose,
                      mc.cores = cores)
  }
  
  study$brain_mask <- brain_mask
  
  ##return a list with the preprocessed images and a brain mask
  return(study)
  
}