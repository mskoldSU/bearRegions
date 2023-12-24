tab1 <- function(){
  #' Table 1
  #' Estimates of brown bear home range size sigma in meters.
  #'
  #' @return A data frame.
  #' @export
  #'
  #' @examples
  #' tab1()
  bind_rows(
    summarise(bear_coordinates |> filter(dist_to_border > 20000,
                                         n_samples > 1),
              ac_sd = median((sd_east + sd_north)/ 2),
              sd_method = "median", .by = "sex"),
    summarise(bear_coordinates |> filter(dist_to_border > 20000,
                                         sqrt((sd_east^2 + sd_north^2) / 2) < 20000,
                                         n_samples > 1),
              ac_sd = sqrt(sum((n_samples -1) * (sd_east^2 + sd_north^2)/ 2) / sum(n_samples - 1)),
              sd_method = "trimmed", .by = "sex"),
    summarise(bear_coordinates |> filter(dist_to_border > 20000,
                                         n_samples > 1),
              ac_sd = sqrt(sum((n_samples -1) * (sd_east^2 + sd_north^2)/ 2) / sum(n_samples - 1)),
              sd_method = "pooled", .by = "sex")
  )
}


tab2 <- function(){
  #' Table 2
  #' Each row in the table lists the number (and percentages) of individuals
  #' located to either region during that particular survey. From the Region A
  #' survey, two individuals living near the coast were placed in the sea (Bothnian bay)
  #' and Region D has some exchange with the Norwegian population.
  #'
  #' @return A data frame.
  #' @export
  #'
  #' @examples
  #' tab2()

  fitted_values_NB |> ungroup() |>
    filter(sd_method == "trimmed") |>
    count(region, survey_region, name = "result") |>
    mutate(result = paste0(result, " (", scales::percent(result / sum(result), accuracy = .1), ")"), .by = "survey_region") |>
    pivot_wider(names_from = "region", values_from = "result", values_fill = "-") |>
    arrange(survey_region) |>
    select(all_of(c("survey_region", "Region A", "Region B", "Region C", "Region D", "Bothnian bay", "Norway")))
}




tabS1 <- function(){
  #' Table S1
  #' Placement of individual activity centers with Negative binomial P(n)
  #'
  #' @return A data frame.
  #' @export
  #'
  #' @examples
  #' tabS1()
  fitted_values_NB |>
    ungroup() |>
    summarise(result = n(), .by = c("sd_method", "region", "survey_region")) |>
    mutate(result = paste0(result, " (", scales::percent(result / sum(result), accuracy = .1), ")"),
           .by = c("sd_method", "survey_region")) |>
    pivot_wider(names_from = "region", values_from = "result", values_fill = "-") |>
    arrange(sd_method, survey_region) |>
    select(all_of(c("sd_method", "survey_region", "Region A", "Region B", "Region C", "Region D", "Bothnian bay", "Norway")))
}

tabS2 <- function(){
  #' Table S2
  #' Placement of individual activity centers with Poisson P(n)
  #'
  #' @return A data frame.
  #' @export
  #'
  #' @examples
  #' tabS2()
  #'
  fitted_values_Po |>
    ungroup() |>
    summarise(result = n(), .by = c("sd_method", "region", "survey_region")) |>
    mutate(result = paste0(result, " (", scales::percent(result / sum(result), accuracy = .1), ")"),
           .by = c("sd_method", "survey_region")) |>
    pivot_wider(names_from = "region", values_from = "result", values_fill = "-") |>
    arrange(sd_method, survey_region) |>
    select(all_of(c("sd_method", "survey_region", "Region A", "Region B", "Region C", "Region D", "Bothnian bay", "Norway")))
}

tabS3 <- function(){
  #' Table S3
  #' Number of individuals placed in neighbouring region using Likelihood based estimate of μ and Threshold.
  #' Negative binomial P(n)
  #'
  #' @return A data frame.
  #' @export
  #'
  #' @examples
  #' tabS3()
  #'
  fitted_values_NB |>
    ungroup() |>
    summarise(n_likelihood = sum(survey_region != region),
              n_threshold = sum(dist_to_border < threshold_NB(k = n_samples,
                                                              sigma = ac_sd,
                                                              theta = NB_pars$theta,
                                                              lambda = NB_pars$lambda
              )
              ), .by = c("survey_region", "sd_method")) |>
    arrange(sd_method, survey_region)
}

tabS4 <- function(){
  #' Table S4
  #' Number of individuals placed in neighbouring region using Likelihood based estimate of μ and Threshold.
  #' Poisson P(n)
  #'
  #' @return A data frame.
  #' @export
  #'
  #' @examples
  #' tabS4()
  #'
  fitted_values_Po |>
    ungroup() |>
    summarise(n_likelihood = sum(survey_region != region),
              n_threshold = sum(dist_to_border < threshold_Po(k = n_samples,
                                                              sigma = ac_sd,
                                                              lambda = Po_par$lambda

              )
              ), .by = c("survey_region", "sd_method")) |>
    arrange(sd_method, survey_region)
}
