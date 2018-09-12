#' @title OASIS glm Object 
#'
#' @description A \code{glm} object containing the OASIS model from the OASIS paper. 
#' @name original_oasis_model
#' @docType data
#' @format A \code{glm} object 
"original_oasis_model"


#' @title Updated OASIS glm Object 
#'
#' @description A \code{glm} object containing the OASIS model from the updated, corrected
#' data.  The original model from the OASIS paper is located at 
#' \code{\link{original_oasis_model}}.
#' @name oasis_model
#' @docType data
#' @format A \code{glm} object 
"oasis_model"

#' @title OASIS model without PD
#'
#' @description A \code{glm} object containing the OASIS model 
#' from the updated, corrected
#' data, fit without the PD data.
#' @name nopd_oasis_model
#' @docType data
#' @format A \code{glm} object 
"nopd_oasis_model"

#' @title OASIS model without PD or T2
#'
#' @description A \code{glm} object containing the OASIS model 
#' from the updated, corrected
#' data, fit without the PD or T2 data.
#' @name not2_nopd_oasis_model
#' @docType data
#' @format A \code{glm} object 
"not2_nopd_oasis_model"


#' @title Example \code{data.frame} of MS voxels
#'
#' @description A \code{data.frame} of MS voxels for OASIS prediction
#' @name example_oasis_df
#' @docType data
#' @format A \code{data.frame} with 100 rows and 13 columns, corresponding to the
#' predictors of the OASIS model.
"example_oasis_df"
