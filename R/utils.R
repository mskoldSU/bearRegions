dist_to_border <- function(east, north, survey_region){
  #'  Calculate the distance to the border of a survey region
  #'
  #' @param east, @param north coordinates of the point
  #' @param survey_region name of the survey region
  #'
  #' @return distance to the border of the survey region
  #'
  #' @export
  region <- survey_regions$geometry[survey_regions$survey_region == survey_region] |>
    sf::st_cast("LINESTRING")
  point <- sf::st_point(c(east, north)) |>
    sf::st_sfc(crs = sf::st_crs(region))
  sf::st_distance(point, region) |> as.numeric()
}

in_region <- function(east, north, region){
  #'  Check if a point is in a survey region
  #'
  #' @param east, @param north coordinates of the point
  #' @param region name of the survey region
  #'
  #' @return TRUE if the point is in the survey region, FALSE otherwise
  #'
  #' @export

  region <- survey_regions$geometry[which(survey_regions$survey_region == region)]
  point <- sf::st_point(c(east, north)) |> sf::st_sfc() |> sf::st_set_crs(sf::st_crs(region))
  sf::st_contains(point, region, sparse = TRUE)
}

region <- function(east, north){
  #'  Find the survey region of a point.
  #'  If the point is not in a survey region, it is assigned to either Norway (west of survey region) or
  #'  Bothnian bay (east of survey region)
  #'
  #' @param east, @param north coordinates of the point
  #'
  #' @return name of the survey region
  #'
  #' @export

  data.frame(east = east, north = north) |>
    sf::st_as_sf(coords = c("east", "north"), crs = sf::st_crs(survey_regions)) |>
    sf::st_join(survey_regions) |>
    dplyr::mutate(survey_region = dplyr::case_when(is.na(survey_region) & (east < 600000) ~ "Norway",
                                     is.na(survey_region) & (east > 600000) ~ "Bothnian bay",
                                     TRUE ~ survey_region)) |>
    dplyr::pull(survey_region)
}

threshold_NB <- function(k, sigma, lambda, theta){
  #'  Calculate the threshold for the Negative Binomial distribution
  #'
  #' @param k number of individuals
  #' @param sigma standard deviation of home range distribution
  #' @param lambda mean of the Negative Binomial distribution
  #' @param theta dispersion parameter of the Negative Binomial distribution
  #'
  #' @return threshold for the Negative Binomial distribution
  #'
  #' @export

  (theta + k) * sigma * lambda / (sqrt(2 * pi) * k * (theta + lambda / 2))
}

threshold_Po <- function(k, sigma, lambda){
  #'  Calculate the threshold for the Poisson distribution
  #'
  #' @param k number of individuals
  #' @param sigma standard deviation of home range distribution
  #' @param lambda mean of the Poisson distribution
  #'
  #' @return threshold for the Poisson distribution
  #'
  #' @export

  lambda * sigma / (sqrt(2 * pi) * k)
}
sigma_table <- function(data, trim = 20000){
  #'  Calculate three versions of sigma by sex
  #'
  #' @param trim threshold for the trimmed mean
  #'
  #' @return data frame with three versions of sigma (ac_sd) and the method used to calculate it (sd_method)
  #'
  #' @export
  bind_rows(
    summarise(data, ac_sd = median((sd_east + sd_north)/ 2, na.rm = TRUE),
              sd_method = "median", .by = "sex"),
    summarise(data |> filter(sqrt((sd_east^2 + sd_north^2) / 2) < trim),
              ac_sd = sqrt(sum((n_samples - 1) * (sd_east^2 + sd_north^2) / 2, na.rm = TRUE) / sum(n_samples - 1)),
              sd_method = "trimmed", .by = "sex"),
    summarise(data |> filter(n_samples > 1),
              ac_sd = sqrt(sum((n_samples -1) * (sd_east^2 + sd_north^2) / 2, na.rm = TRUE) / sum(n_samples - 1)),
              sd_method = "pooled", .by = "sex")
  )
}



