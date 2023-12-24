#' Coordinates of bears in the most recent survey of each survey region
#'
#'
#' @format A data frame with 2763 rows and 9 variables:
#' \describe{
#'   \item{id}{ID of the bear from Rovbase 3.0}
#'   \item{survey_region}{Region where the bear was observed}
#'   \item{n_samples}{Number of scat samples collected from the bear}
#'   \item{mean_east}{Mean easting of the scat samples (SWEREF 99 TM)}
#'   \item{mean_north}{Mean northing of the scat samples (SWEREF 99 TM)}
#'   \item{sd_east}{Standard deviation of the easting of the scat samples (SWEREF 99 TM). NA if n_samples = 1}
#'   \item{sd_north}{Standard deviation of the northing of the scat samples (SWEREF 99 TM). NA if n_samples = 1}
#'   \item{dist_to_border}{Distance to the border of the survey region (m)}
#' }
"bear_coordinates"

#' Polygons of survey regions
#'
#'
#' @format An sf object with four polygons corresponding to survey region borders
#' \describe{
#'   \item{survey_region}{Name of the survey region}
#'   \item{geometry}{Polygon of the survey region (SWEREF 99 TM)}
#' }
"survey_regions"
