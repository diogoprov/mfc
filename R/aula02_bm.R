# Aula 2, bloco 2: o movimento browniano. Simulações no lugar das figuras
# do livro de L. Harmon (CC BY-SA), com a mesma sequência do original.
library(patchwork)
source(here::here("R/tema.R"))
source(here::here("R/bm.R"))

#' Partículas de pólen na água: passeios aleatórios em duas dimensões
fig_brown <- function(n = 5, passos = 300, seed = 3) {
  set.seed(seed)
  d <- purrr::map_dfr(seq_len(n), \(k) tibble(
    p = k, x = cumsum(c(0, rnorm(passos))), y = cumsum(c(0, rnorm(passos)))))
  ini <- d |> group_by(p) |> slice_head(n = 1)
  ggplot(d, aes(x, y, colour = factor(p))) +
    geom_path(linewidth = .5, alpha = .9) +
    geom_point(data = ini, size = 3) +
    scale_colour_manual(values = cores_linhagem, guide = "none") +
    coord_equal() + theme_void() +
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}

#' Um passeio aleatório em uma dimensão (o atributo ao longo do tempo)
fig_passeio <- function(seed = 5) {
  ggplot(sim_bm_tempo(1, 200, seed = seed), aes(t, x)) +
    geom_line(colour = cores[["accent"]], linewidth = .9) +
    labs(x = "tempo", y = "atributo") +
    tema_aula(16, conceitual = TRUE)
}

#' Soma de seis dados: uma série (etapa 1) e muitas, com o histograma (2)
fig_dados <- function(etapa = 1, n_series = 60) {
  um <- sim_dados(6, 10, 1, seed = 7)
  base <- ggplot(mapping = aes(jogada, soma)) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = seq(6, 36, 6)) + coord_cartesian(ylim = c(6, 36)) +
    labs(x = "jogada", y = "soma de 6 dados") + tema_aula(16)
  if (etapa == 1)
    return(base + geom_line(data = um, colour = cores[["accent"]], linewidth = 1) +
             geom_point(data = um, colour = cores[["accent"]], size = 3) +
             geom_text(data = um, aes(label = soma), vjust = -1, size = 5))
  muitos <- sim_dados(6, 10, n_series, seed = 8)
  linhas <- base + geom_line(data = muitos, aes(group = serie), colour = cores[["accent"]],
                             alpha = .35, linewidth = .5)
  hist <- ggplot(muitos, aes(y = soma)) +
    geom_histogram(binwidth = 2, fill = cores[["ocre"]], colour = cores[["papel"]]) +
    coord_cartesian(ylim = c(6, 36)) +
    labs(x = "frequência", y = NULL) + tema_aula(16) +
    theme(axis.text.y = element_blank())
  linhas + hist + plot_layout(widths = c(3, 1))
}

#' Uma trajetória de BM crescendo no tempo (etapas 1–3) e muitas, com o
#' histograma do valor final (etapa 4)
fig_bm_traj <- function(etapa = 1, passos = 100) {
  um <- sim_bm_tempo(1, passos, seed = 21)
  lim <- c(-30, 30)
  base <- ggplot(mapping = aes(t, x)) +
    scale_x_continuous(limits = c(0, passos)) + coord_cartesian(ylim = lim) +
    labs(x = "tempo", y = "atributo") + tema_aula(16)
  if (etapa <= 3) {
    ate <- c(33, 66, 100)[etapa]
    return(base + geom_line(data = filter(um, t <= ate), colour = cores[["accent"]], linewidth = 1))
  }
  muitos <- sim_bm_tempo(200, passos, seed = 22)
  linhas <- base + geom_line(data = muitos, aes(group = rep), colour = cores[["accent"]],
                             alpha = .15, linewidth = .4) +
    geom_line(data = um, colour = cores[["accent_dark"]], linewidth = 1)
  hist <- ggplot(filter(muitos, t == passos), aes(y = x)) +
    geom_histogram(binwidth = 3, fill = cores[["ocre"]], colour = cores[["papel"]]) +
    coord_cartesian(ylim = lim) + labs(x = "frequência", y = NULL) + tema_aula(16) +
    theme(axis.text.y = element_blank())
  linhas + hist + plot_layout(widths = c(3, 1))
}

#' Linhagens que se separam (BM com especiação no tempo)
#'
#' @param k 2 ou 3 linhagens
#' @param muitas FALSE = uma réplica; TRUE = 150 réplicas e o histograma do
#'   valor final de cada linhagem
fig_linhagens <- function(k = 2, muitas = FALSE, reps = 150) {
  t_div <- c(0, 50, 75)[1:k]; mae <- c(NA, 1, 2)[1:k]
  um <- sim_linhagens(t_div, mae, seed = 31)
  lim <- c(-32, 32)
  base <- ggplot(mapping = aes(t, x, colour = factor(linhagem))) +
    scale_colour_manual(values = cores_linhagem, guide = "none") +
    coord_cartesian(ylim = lim) + labs(x = "tempo", y = "atributo") + tema_aula(16)
  desenha <- function(d, alfa, larg)
    geom_line(data = filter(d, t >= nasce), aes(group = interaction(rep, linhagem)),
              alpha = alfa, linewidth = larg)
  um$rep <- 0
  if (!muitas) {
    nos <- tibble(t = t_div[-1],
                  x = purrr::map2_dbl(t_div[-1], mae[-1], \(tt, m) um$x[um$linhagem == m & um$t == tt]))
    return(base + desenha(um, 1, 1) +
             geom_point(data = nos, aes(t, x), inherit.aes = FALSE, size = 3.5,
                        colour = cores[["texto"]]))
  }
  muitas_d <- purrr::map_dfr(seq_len(reps), \(r) sim_linhagens(t_div, mae, seed = 1000 + r) |> mutate(rep = r))
  linhas <- base + desenha(muitas_d, .12, .35)
  fim <- filter(muitas_d, t == max(t))
  hist <- ggplot(fim, aes(y = x, fill = factor(linhagem))) +
    geom_histogram(binwidth = 3, colour = cores[["papel"]], position = "identity", alpha = .6) +
    scale_fill_manual(values = cores_linhagem, guide = "none") +
    coord_cartesian(ylim = lim) + labs(x = "frequência", y = NULL) + tema_aula(16) +
    theme(axis.text.y = element_blank())
  linhas + hist + plot_layout(widths = c(3, 1))
}

#' Teorema do limite central: a soma de muitos passos ±1 é normal
fig_tlc <- function(passos = 100, reps = 10000, seed = 41) {
  set.seed(seed)
  s <- replicate(reps, sum(sample(c(-1, 1), passos, replace = TRUE)))
  ggplot(tibble(s), aes(s)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 2,
                   fill = cores[["ocre"]], colour = cores[["papel"]]) +
    stat_function(fun = dnorm, args = list(0, sqrt(passos)), colour = cores[["accent_dark"]],
                  linewidth = 1.2) +
    labs(x = paste("soma de", passos, "passos de +1 ou −1"), y = "densidade") + tema_aula(16)
}

#' Duas taxas σ², mesmas sementes
fig_sigmas <- function(s2 = c(.01, .1), reps = 100, passos = 100) {
  d <- purrr::map_dfr(s2, \(s) sim_bm_tempo(reps, passos, sigma2 = s, seed = 51) |>
                        mutate(painel = paste0("sigma^2 == \"", format(s, decimal.mark = ","), "\"")))
  ggplot(d, aes(t, x, group = rep)) +
    geom_line(colour = cores[["accent"]], alpha = .3, linewidth = .4) +
    facet_wrap(~painel, labeller = label_parsed) +
    labs(x = "tempo", y = "atributo") + tema_aula(16)
}

#' Variância entre réplicas ao longo do tempo e a reta σ²t
var_tempo <- function(reps = 2000, passos = 100, s2 = 1) {
  sim_bm_tempo(reps, passos, sigma2 = s2, seed = 61) |>
    group_by(t) |> summarise(v = var(x), .groups = "drop")
}
fig_var_tempo <- function(s2 = 1) {
  d <- var_tempo(s2 = s2)
  ggplot(d, aes(t, v)) +
    geom_abline(slope = s2, intercept = 0, colour = cores[["ferrugem"]], linewidth = 1.2,
                linetype = "22") +
    geom_point(colour = cores[["accent"]], size = 1.8) +
    annotate("text", x = 70, y = 40, label = "Var(X) == sigma^2 * t", parse = TRUE,
             colour = cores[["ferrugem"]], size = 7) +
    labs(x = "tempo (t)", y = "variância entre réplicas") + tema_aula(16)
}

#' Quatro processos que produzem um padrão de BM
fig_quando_bm <- function(n = 12) {
  set.seed(71)
  passos <- 200
  pontuada <- purrr::map_dfr(seq_len(n), \(r) {
    salto <- rbinom(passos, 1, .05) * rnorm(passos, 0, 3)
    tibble(rep = r, t = 0:passos, x = c(0, cumsum(salto)))
  })
  flutuante <- purrr::map_dfr(seq_len(n), \(r) {
    theta <- cumsum(rnorm(passos + 1, 0, .8)); x <- numeric(passos + 1)
    for (k in seq_len(passos)) x[k + 1] <- x[k] + .3 * (theta[k] - x[k]) + rnorm(1, 0, .3)
    tibble(rep = r, t = 0:passos, x = x)
  })
  fraca <- sim_ou_tempo(n, passos, sigma2 = 1, alfa = .002, theta = 0, seed = 72)
  d <- bind_rows(
    sim_bm_tempo(n, passos, seed = 73) |> mutate(p = "Deriva genética"),
    pontuada |> mutate(p = "Mudança pontuada\naleatória"),
    fraca |> mutate(p = "Seleção fraca\nna escala de tempo"),
    flutuante |> mutate(p = "Seleção que muda\nde direção ao acaso")) |>
    mutate(p = factor(p, levels = unique(p)))
  ggplot(d, aes(t, x, group = rep)) +
    geom_line(colour = cores[["accent"]], alpha = .5, linewidth = .5) +
    facet_wrap(~p, nrow = 1, scales = "free_y") +
    labs(x = "tempo", y = "atributo") + tema_aula(15, conceitual = TRUE)
}
