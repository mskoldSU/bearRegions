fig1A <- function(){
  #' Figure 1A
  #'
  #' Original and transformed coordinate axes.
  #' The y′ axis is aligned with the border (dashed line) to the neighbouring region (shaded area)
  #' while the x′ axis passes through the mean observed coordinate (filled point).
  #' In the new coordinate system, observed coordinates (circles) are represented by their distance
  #' d = x′ to the border (represented by dotted lines). The circle shows a 90% contour line of the
  #' home region with center μ (cross).
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig1A()
  polygon <- tibble(x = c(-2, 3, 3, -2), y = c(3, 3, -2, 3))
  points <- tibble(x = c(-1, -0.5, .5), y = c(1, -1, -.3), in_region = (x + y) < 1,
                   x0 = x/2 + 1/2 -y/2, y0 = (1 - x)/2 + y/2) |>
    mutate(x0 = if_else(in_region == TRUE, x0, NA),
           y0 = if_else(in_region == TRUE, y0, NA))

  center <- points |>
    filter(in_region == TRUE) |>
    summarise(x = mean(x),
              y = mean(y),
              x0 = mean(x0),
              y0 = mean(y0))

  fontsize <- 5

  p <- ggplot(points) + geom_point(aes(x = x, y = y), shape = 1, size = 5) +
    geom_abline(slope = -1, intercept = 1, linetype = "dashed") +
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = 2.17), color = "grey") +
    geom_segment(aes(x = x, y = y, xend = x0, yend = y0), linetype = "dotted") +
    geom_point(data = center, aes(x = x, y = y), size = 5) +
    geom_segment(data = center, aes(x = x, y = y, xend = x0, yend = y0)) +
    annotate("point", x = 0, y = 0, shape = 3) +
    geom_polygon(data = polygon, aes(x = x, y = y), fill = "black", alpha = .2, color = NA) +
    coord_equal(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
    annotate("segment", x = -2.5, xend = 2, y = -2, yend = -2, arrow = arrow()) +
    annotate("segment", y = -2.5, yend = 2, x = -2, xend = -2, arrow = arrow()) +
    annotate("segment", y = center$y0- 2*(center$y - center$y0), yend = center$y0 + 2*(center$y - center$y0),
             x = center$x0 - 2*(center$x - center$x0), xend = center$x0 + 2*(center$x - center$x0), arrow = arrow()) +
    annotate("segment", y = center$y0+ 2.5*(center$y - center$y0), yend = center$y0 - 1.5*(center$y - center$y0),
             x = center$x0 - 2.5*(center$x - center$x0), xend = center$x0 + 1.5*(center$x - center$x0), arrow = arrow()) +
    annotate("text", label = "x", x = 2.2, y = -2, size = fontsize) +
    annotate("text", label = "y", y = 2.2, x = -2, size = fontsize) +
    annotate("text", label = "x'", x = center$x0 + 2.2*(center$x - center$x0), y = center$y0 + 2.2*(center$y - center$y0), angle = 45, size = fontsize) +
    annotate("text", label = "y'", x = center$x0 + 1.7*(center$x - center$x0) -.1, y = center$y0 - 1.7*(center$y - center$y0)-.1, angle = 45, size = fontsize) +
    annotate("text", label = latex2exp::TeX("$\\delta$"), x = -0.25, y = 0.25, angle = 45, size = fontsize) +
    annotate("segment", x = -.07, y = .07, xend = -.16, yend = .16) +
    annotate("text", label = "d", x = -0.15 + center$x, y = 0.15+center$y, angle = 45, size = fontsize) +
    theme_void()
  p
}

fig1B <- function(){
  #' Figure 1B
  #'
  #' MLE of delta as a function of observed average distance to border d in units of sigma,
  #' for selected ratios of lambda and k, where n is Poisson(lambda).
  #' The absolute value |delta| correspond to the estimated distance of the activity center to the
  #' border while the sign denotes its side with negative values in the neighbouring region
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig1B()
  log_lik <- function(delta, ratio, d) {ratio * pnorm(0, delta) - (d - delta)^2/2}

  delta_hat <- expand_grid(ratio = c(.5, 1, 2, 4),
                         d = seq(0, 2.5, length.out = 500)) |>
    rowwise() |>
    mutate(delta = optimise(interval = c(-10, 8), function(delta) -log_lik(delta, ratio, d))$minimum)

  p <- delta_hat |>
    ggplot(aes(x = d, y = delta, color = as.factor(ratio))) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 2.5), breaks = c(0, 1, 2), expand = expansion(mult = c(0, .05)),
                       labels = latex2exp::TeX(c("0", "$\\sigma$", "$2\\sigma$"))) +
    scale_y_continuous(breaks = c(-1, 0, 1, 2),
                       labels = latex2exp::TeX(c( "$-\\sigma$", "0", "$\\sigma$", "$2\\sigma$")),
                       limits = c(-1, 2.5)) +
    scale_color_discrete(labels = c(latex2exp::TeX("$\\lambda/k = 1/2$"),
                                    latex2exp::TeX("$\\lambda/k = 1$"),
                                    latex2exp::TeX("$\\lambda/k = 2$"),
                                    latex2exp::TeX("$\\lambda/k = 4$"))) +
    geom_abline(intercept = 0, slope = 0) +
    labs(x = latex2exp::TeX("Observed average distance to border $\\bar{d}$"),
         y = latex2exp::TeX("Estimated $\\delta$"), color = "")
  p
}

fig2 <- function(){
  #' Figure 2
  #'
  #' The four regions of Swedish national monitoring of the brown bear. Dots represent average scat sample
  #' locations of unique individuals encountered in each of the regions most recent survey.
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig2()
  p <- bear_coordinates |>  ggplot() +
    geom_sf(data = SWEcounties, color = NA, fill = "grey") +
    geom_sf(data = sf::st_as_sf(survey_regions), aes(fill = survey_region), alpha = .5) +
    geom_point(aes(x = mean_east, y = mean_north), size = .1, alpha = .5) +
    theme_void() + labs(fill = "") +
    ggspatial::annotation_scale(location = "br", style = "tick") +
    ggspatial::annotation_north_arrow(location = "tl",
                                      height = unit(1, "cm"),
                                      width = unit(1, "cm"))

  p
}

fig3A <- function(){
  #' Figure 3A
  #'
  #' Number of individuals with a given number of traces recorded (bars). Filled circles represent a
  #' zero-truncated Negative Binomial fit and empty a zero-truncated Poisson.
  #' The x-axis is truncated, excluding a small number of individuals with k > 20
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig3A()
  k <- bear_coordinates |>
    filter(dist_to_border > 20000) |>
    pull(n_samples)
  n_k <- length(k)
  parsNB <- fit_ztNB(k)
  parsPo <- fit_ztPo(k)
  p0NB <- pnbinom(0, parsNB$theta, mu = parsNB$lambda, lower.tail = FALSE)
  p0Po <- ppois(0, parsPo$lambda, lower.tail = FALSE)

  p <- tibble(k = k) |> ggplot(aes(x = k)) + geom_bar(fill = "darkgrey") +
    geom_line(data = tibble(k = 1:30,
                            p =  n_k * dnbinom(k, parsNB$theta, mu = parsNB$lambda) / p0NB),
              aes(y = p), alpha = .5) +
    geom_point(data = tibble(k = 1:30,
                             p =  n_k * dnbinom(k, parsNB$theta, mu = parsNB$lambda) / p0NB),
               aes(y = p)) +
    geom_line(data = tibble(k = 1:30,
                            p =  n_k * dpois(k, parsPo$lambda) / p0Po),
              aes(y = p), alpha = .5) +
    geom_point(data = tibble(k = 1:30,
                             p =  n_k * dpois(k, parsPo$lambda) / p0Po),
               pch = 21, fill = "white", aes(y = p)) +
    theme_bw() +
    labs(x = "Number of traces per individual", y = "Number of individuals") +
    scale_x_continuous(limits = c(0, 20)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, .05))) +
    theme(panel.grid = element_blank())
  p
}

fig3B <- function(){
  #' Figure 3B
  #'
  #' Histograms of observed individual coordinate sample standard deviations (grey). Dashed
  #' lines show level of truncation. Outlined (black lines) are corresponding histograms
  #' based on simulated coordinates under model assumptions with σ = spool.
  #' The x-axis is truncated, excluding a small number of individuals with s > 50000
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig3B()
  sim_data <- bear_coordinates |>
    mutate(ac_sd = sqrt(sum((n_samples - 1) * (sd_east^2 + sd_north^2), na.rm = TRUE) / 2 / sum(n_samples - 1)), .by = "sex") |>
    rowwise() |>
    mutate(sd_east = sd(rnorm(n_samples, 0, ac_sd)), sd_north = sd(rnorm(n_samples, 0, ac_sd)),
           sd = sqrt((sd_east^2 + sd_north^2)/2),
           sd = 1000 * floor(sd / 1000)) |>
    count(sd, sex) |>
    filter(!is.na(sd))
  p <- bear_coordinates |>
    mutate(sd = sqrt((sd_east^2 + sd_north^2) / 2)) |>
    ggplot(aes(x = sd)) +
    stat_bin(binwidth = 1000, fill = "darkgrey", center = 500) + geom_rug(alpha = .5) +
    geom_step(data = sim_data, aes(x = sd, y = n), color = "black", linewidth = .2) +
    geom_vline(xintercept = 20000, linetype = "dashed", color = "darkgrey") +
    theme_bw() +
    scale_x_continuous(limits = c(0, 50000), expand = expansion(c(0, .05))) +
    facet_wrap(~sex, ncol = 1, scales = "free_y") +
    labs(x = "Coordinate sample standard deviation (m)", y = "Number of individuals") +
    theme(panel.grid = element_blank())
  p
}

fig4A <- function(){
  #' Figure 4A
  #'
  #' Illustration of the relocation process in a high density area along the border between Region C and D.
  #' Arrows, colored according to survey, start at the observed sample mean coordinate and point at the
  #' estimated activity center.
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig4A()
  p <- fitted_values_NB |>
    filter(sd_method == "trimmed", survey_region %in% c("Region C", "Region D")) |>
    mutate(n_samples = pmin(n_samples, 5)) |>
    ggplot() +
    geom_sf(data = survey_regions, fill = NA, linetype = "dashed", linewidth = 1, color = "grey") +
    geom_segment(aes(x = mean_east, y = mean_north, xend = mu_east, yend = mu_north, color = survey_region),
                 arrow = arrow(type = "open", length = unit(0.2, "cm")), show.legend = FALSE) +
    geom_point(aes(x = mean_east, y = mean_north, size = n_samples, color = survey_region), shape = 1) +
    coord_sf(xlim = c(4.8e5, 5.2e5), ylim = c(6858000, 6900000), ndiscr = 0) +
    theme_bw() + scale_x_continuous(breaks = NULL, minor_breaks = NULL) +
    scale_y_continuous(breaks = NULL, minor_breaks = NULL) +
    theme(panel.grid.major = element_blank()) + labs(x = "", y = "", color = "", size = "Samples") +
    scale_size(labels = c(1, 2, 3, 4, "> 4"), breaks = 1:5) +
    ggspatial::annotation_scale(location = "br", style = "tick") +
    guides(color = "none") +
    theme(legend.position = "top")

  p
}

fig4B <- function(){
  #' Figure 4B
  #'
  #' Number of individuals placed outside the survey region in relation to estimated home range size σ and
  #' sampling distribution. Dots show estimated numbers based on the home range sizes reported in Table 1, lines show
  #' fitted proportional relations (straight lines with zero intercept)
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig4B()
  fitted <- bind_rows(
    fitted_values_Po |> mutate(model = "Poisson", lambda = Po_par$lambda,),
    fitted_values_NB |> mutate(model = "Negative binomial", lambda = NB_pars$lambda)
  )
  not_inside <- fitted |>
    filter(survey_region != region) |>
    count(sex, ac_sd, sd_method, model, lambda)
  p <- not_inside |>
    ggplot(aes(x = ac_sd, y = n, color = sex)) +
    geom_point() +
    facet_wrap(~model) +
    theme_bw() +
    geom_smooth(aes(group = sex), method = "lm", formula = y ~ -1 + x, se = FALSE, linewidth = .5) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, .05))) +
    scale_x_continuous(limits = c(2000, 11000)) +
    labs(x = latex2exp::TeX("Estimated home range size ($\\sigma$)"), y = "Estimated number of individuals", color = "") +
    theme(legend.position = "top")
  p
}

fig5 <- function(){
  #' Figure 5
  #'
  #' Kernel density estimates of simulated values for delta. Means of simulated values are represented by vertical
  #' dashed lines. Negative values corresponds to individuals in the neighbouring region.
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig5()
  range <-   c(-1, -0.5, 0, 0.5, 1, 2)
  p <- simulation_results |> filter(delta %in% range) |>
    ggplot(aes(x = delta_hat, color = as.factor(delta))) +
    stat_density(adjust = 1.2, geom = "line",position = "identity") +
    facet_wrap(~lambda, labeller = as_labeller(function(string) latex2exp::TeX(paste0("$\\lambda =$", string)), default = label_parsed)) +
    theme_bw() +
    scale_x_continuous(breaks = range,
                       labels = latex2exp::TeX(c( "$-\\sigma$", "$-0.5\\sigma$", "$0$", "$0.5\\sigma$", "$\\sigma$", "$2\\sigma$")), minor_breaks = NULL) +
    coord_cartesian(xlim = c(-1, 3)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, .05)), minor_breaks = NULL, breaks = NULL) +
    geom_vline(data = simulation_results |>
                 ungroup() |>
                 filter(delta %in% range) |>
                 summarise(deltahat_mean = mean(delta_hat), .by = c("lambda", "delta")),
               aes(xintercept = deltahat_mean, color = as.factor(delta)), linetype = "dashed", show.legend = FALSE) +
    labs(y = "", x = latex2exp::TeX("Estimated $\\delta$"), color = "") +
    scale_color_discrete(labels = c(latex2exp::TeX(paste("$\\delta=-\\sigma$")),
                                    latex2exp::TeX(paste("$\\delta=-0.5\\sigma$")),
                                    latex2exp::TeX(paste("$\\delta=0$")),
                                    latex2exp::TeX(paste("$\\delta=0.5\\sigma$")),
                                    latex2exp::TeX(paste("$\\delta=\\sigma$")),
                                    latex2exp::TeX(paste("$\\delta=2\\sigma$")))) +
    geom_vline(xintercept = 0, linewidth = .1) + theme(legend.position = "top") +
    guides(color = guide_legend(nrow = 1))
  p
}

fig6 <- function(){
  #' Figure 6
  #'
  #' Probability of a discovered individual being placed within the survey region given the distance|delta| of its
  #' activity center to the border. Negative delta corresponds to individuals located in the neighbouring region
  #'
  #' @return A ggplot object.
  #' @export
  #'
  #' @examples
  #' fig6()
  p <- simulation_results |> ungroup() |>
    summarise(p = mean(delta_hat > 0), .by = c("delta", "lambda")) |>
    ggplot(aes(x = delta, y = p, color = as.factor(lambda))) +
    geom_point() + geom_line(alpha = .5) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(0)) +
    theme_bw() + labs(x = latex2exp::TeX("Distance from border $\\delta$"), y = "Probability of being placed in survey region", color = "") +
    scale_color_discrete(labels = c(latex2exp::TeX(paste("$\\lambda=1$")),
                                    latex2exp::TeX(paste("$\\lambda=2$")),
                                    latex2exp::TeX(paste("$\\lambda=3$")),
                                    latex2exp::TeX(paste("$\\lambda=4$")))) +
    geom_vline(xintercept = 0, linewidth = .1) + theme(legend.position = "top") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2), labels = latex2exp::TeX(c("$-\\sigma$", "0", "\\sigma", "2\\sigma")))
  p
}








