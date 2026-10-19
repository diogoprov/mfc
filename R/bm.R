library(dplyr)
library(tibble)
library(purrr)

# Simulações de evolução de um atributo contínuo ao longo do tempo (sem
# árvore), para as figuras conceituais. Todas fixam a própria semente.

#' Movimento browniano: `n` trajetórias independentes
#' @return tibble com rep, t e x
sim_bm_tempo <- function(n = 1, passos = 100, sigma2 = 1, z0 = 0, seed = 1) {
  set.seed(seed)
  purrr::map_dfr(seq_len(n), \(r) tibble(
    rep = r, t = 0:passos,
    x = z0 + c(0, cumsum(rnorm(passos, 0, sqrt(sigma2))))))
}

#' Ornstein-Uhlenbeck em tempo discreto (Euler): puxão α(θ − x) a cada passo
sim_ou_tempo <- function(n = 1, passos = 100, sigma2 = 1, alfa = .1,
                         theta = 0, z0 = 0, seed = 1) {
  set.seed(seed)
  purrr::map_dfr(seq_len(n), \(r) {
    x <- numeric(passos + 1); x[1] <- z0
    for (k in seq_len(passos)) x[k + 1] <- x[k] + alfa * (theta - x[k]) + rnorm(1, 0, sqrt(sigma2))
    tibble(rep = r, t = 0:passos, x = x)
  })
}

#' Early burst: a variância de cada passo cai como σ²₀·e^(r·t), com r < 0
sim_eb_tempo <- function(n = 1, passos = 100, sigma2 = 4, r = -.04, z0 = 0, seed = 1) {
  set.seed(seed)
  purrr::map_dfr(seq_len(n), \(k) tibble(
    rep = k, t = 0:passos,
    x = z0 + c(0, cumsum(rnorm(passos, 0, sqrt(sigma2 * exp(r * (0:(passos - 1)))))))))
}

#' BM com tendência: cada passo tem média μ
sim_trend_tempo <- function(n = 1, passos = 100, sigma2 = 1, mu = .3, z0 = 0, seed = 1) {
  set.seed(seed)
  purrr::map_dfr(seq_len(n), \(r) tibble(
    rep = r, t = 0:passos,
    x = z0 + c(0, cumsum(rnorm(passos, mu, sqrt(sigma2))))))
}

#' Linhagens que se separam: BM numa "árvore" desenhada no tempo
#'
#' A linhagem 1 começa em t = 0; a linhagem k nasce da linhagem `mae[k]` no
#' tempo `t_div[k]` e segue independente dali em diante.
#' @return tibble com linhagem, t, x
sim_linhagens <- function(t_div = c(0, 40, 70), mae = c(NA, 1, 2),
                          passos = 100, sigma2 = 1, seed = 1) {
  set.seed(seed)
  inc <- matrix(rnorm(length(t_div) * passos, 0, sqrt(sigma2)), ncol = length(t_div))
  traj <- matrix(NA_real_, passos + 1, length(t_div))
  traj[, 1] <- c(0, cumsum(inc[, 1]))
  for (k in seq_along(t_div)[-1]) {
    t0 <- t_div[k]
    traj[1:(t0 + 1), k] <- traj[1:(t0 + 1), mae[k]]
    traj[(t0 + 2):(passos + 1), k] <- traj[t0 + 1, mae[k]] + cumsum(inc[(t0 + 1):passos, k])
  }
  purrr::map_dfr(seq_along(t_div), \(k) tibble(
    linhagem = k, t = 0:passos, x = traj[, k], nasce = t_div[k], mae = mae[k]))
}

#' Soma de `n_dados` dados, jogados `jogadas` vezes, em `n_series` séries
sim_dados <- function(n_dados = 6, jogadas = 10, n_series = 1, seed = 1) {
  set.seed(seed)
  purrr::map_dfr(seq_len(n_series), \(s) tibble(
    serie = s, jogada = seq_len(jogadas),
    soma = replicate(jogadas, sum(sample(1:6, n_dados, replace = TRUE)))))
}
