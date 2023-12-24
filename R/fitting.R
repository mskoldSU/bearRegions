LpolycubNB <- function(mu, x, y, k, sigma, region, lambda, theta){
  #' Calculate the log-likelihood of the Negative Binomial model
  P_mu <- polyCub::polyCub.exact.Gauss(region, mean = mu, Sigma = diag(c(sigma^2, sigma^2)))
  - (theta + k)*log(theta/lambda + P_mu) - k *((x - mu[1])^2 + (y - mu[2])^2)/(2*sigma^2)
}

fit_center_NB <- function(east, north, k, sigma, region_name, lambda, theta){
  #' Numerically optimize the log-likelihood of the Negative Binomial model
  #'
  #' @param east, @param north average sample coordinates
  #' @param k number of samples
  #' @param sigma standard deviation of home range distribution
  #' @param region_name name of the region in which the samples were collected (as in survey_regions table)
  #' @param lambda mean of the NB distribution
  #' @param theta dispersion parameter of the NB distribution
  #'
  #' @return: estimated center of the home range and proportion of the region occupied by the home range
  #'
  #' @export

  # Area over which to integrate is a 10 sigma box around the point intersected with the region
  clipbox <- sf::st_bbox(c(xmin = east - 10*sigma, xmax = east + 10*sigma,  ymax = north + 10*sigma, ymin = north - 10*sigma),
                         crs = sf::st_crs(survey_regions)) |>
    sf::st_as_sfc()
  region <- sf::st_intersection(clipbox, survey_regions$geometry[survey_regions$survey_region == region_name])
  # Optimize the log-likelihood
  mu_hat <- optim(c(east, north), function(mu) -LpolycubNB(mu, east, north, k, sigma, region[[1]], lambda, theta),
                  lower = (c(east, north) - 4 * sigma), upper = (c(east, north) + 4 * sigma), method = "L-BFGS-B",
                  control = list(parscale = c(1000, 1000)))
  P_mu <- polyCub::polyCub.exact.Gauss(region[[1]], mean = mu_hat$par, Sigma = diag(c(sigma^2, sigma^2)))
  list(mu_east = mu_hat$par[1], mu_north = mu_hat$par[2], p = P_mu)
}

LpolycubPo <- function(mu, x, y, k, sigma, region, lambda){
  #' Calculate the log-likelihood of the Poisson model
  P_mu <- polyCub::polyCub.exact.Gauss(region, mean = mu, Sigma = diag(c(sigma^2, sigma^2)))
  - P_mu * lambda - k * ((x - mu[1])^2 + (y - mu[2])^2) / (2*sigma^2)
}

fit_center_Po <- function(east, north, k, sigma, region_name, lambda){
  #' Numerically optimize the log-likelihood of the Poisson model
  #'
  #' @param east, @param north average sample coordinates
  #' @param k number of samples
  #' @param sigma standard deviation of home range distribution
  #' @param region_name name of the region in which the samples were collected (as in survey_regions table)
  #' @param lambda mean of the Poisson distribution
  #'
  #' @return: estimated center of the home range and proportion of the region occupied by the home range
  #'
  #' @export

  # Area over which to integrate is a 10 sigma box around the point intersected with the region
  clipbox <- sf::st_bbox(c(xmin = east - 10*sigma, xmax = east + 10*sigma,  ymax = north + 10*sigma, ymin = north - 10*sigma),
                         crs = sf::st_crs(survey_regions)) |>
    sf::st_as_sfc()
  region <- sf::st_intersection(clipbox, survey_regions$geometry[survey_regions$survey_region == region_name])
  # Optimize the log-likelihood
  mu_hat <- optim(c(east, north), function(mu) -LpolycubPo(mu, east, north, k, sigma, region[[1]], lambda),
                  lower = (c(east, north) - 4*sigma), upper = (c(east, north) + 4*sigma), method = "L-BFGS-B",
                  control = list(parscale = c(1000, 1000)))
  P_mu <- polyCub::polyCub.exact.Gauss(region[[1]], mean = mu_hat$par, Sigma = diag(c(sigma^2, sigma^2)))
  list(mu_east = mu_hat$par[1], mu_north = mu_hat$par[2], p = P_mu)
}


fit_ztNB <- function(n){
  #' Numerically optimize the log-likelihood of the zero-truncated Negative Binomial distribution
  #'
  #' @param n vector of counts
  #' @return estimated parameters of the zero-truncated Negative Binomial distribution
  #'
  #' @export
  loglik <- function(pars, n){
    -sum(dnbinom(n, pars[1], mu = pars[2], log = TRUE) - pnbinom(0, pars[1], mu =  pars[2], log.p = TRUE, lower.tail = FALSE))
  }
  fit <- nlm(loglik, c(5, mean(n)), n = n)
  theta <- fit$estimate[1]
  lambda <- fit$estimate[2]
  list(theta = theta, lambda = lambda)
}

fit_ztPo <- function(n){
  #' Numerically optimize the log-likelihood of the zero-truncated Poisson distribution
  #'
  #' @param n vector of counts
  #'
  #' @return estimated parameter of the zero-truncated Poisson distribution
  #'
  #' @export
  loglik <- function(lambda, n){
    -sum(dpois(n, lambda, log = TRUE) - ppois(0, lambda, log.p = TRUE, lower.tail = FALSE))
  }
  fit <- nlm(loglik, mean(n), n = n)
  list(lambda = fit$estimate)
}

