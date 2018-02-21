#' @title OASIS Erode Mask
#' @description An alternative to using \code{fslerode} for mask erosion of a 
#' brain mask by a box kernel defined by millimeter
#' @param mask object of class \code{\link{nifti}}
#' @param mm Number of erosion (in millimeters)
#'
#' @return Object of class \code{\link{nifti}}
#' @export
#' @importFrom mmand erode shapeKernel
#' @importFrom neurobase niftiarr zero_pad
#' @importFrom oro.nifti voxdim
#' @examples
#' library(neurobase)
#' library(fslr)
#' library(oasis)
#' niis = tempfile(fileext = ".nii.gz")
#' if (require(httr)) {
#'    url = paste0("https://s3.us-east-2.amazonaws.com/brainder/software/",
#'    "flair/templates/GG-853-FLAIR-2.0mm.nii.gz")
#'    req <- httr::GET(url,
#'    httr::write_disk(path = niis))
#'    httr::stop_for_status(req)
#'
#'  flair <- fast_readnii(niis)
#'  res = oasis_erode(flair > 50, mm = c(2, 2, 2))
#' }
oasis_erode = function(mask, mm = c(5,5,5)) {

  mask = check_nifti(mask)
  vd = voxdim(mask)
  
  if ( length(mm) < 3 ) {
    mm = c(mm, rep(mm, length = 3 - length(mm)))
  }
  
  nvoxels = round(mm/vd)
  # nvoxels = c(5,5,5)
  ### check voxel spec
  stopifnot(length(nvoxels) %in% c(1, 3))
  
  # is.wholenumber <- function(x, tol = .Machine$double.eps ^ 0.5) {
  #   abs(x - round(x)) < tol  
  # }
  
  box = shapeKernel(width = nvoxels, 
                    type = "box")
  img = zero_pad(mask, kdim = nvoxels)
  img = erode(img, kernel = box)
  img = zero_pad(img, kdim = nvoxels, invert = TRUE)
  ero = niftiarr(mask, img)
  return(ero)
}