# Aula 5: modelos de evolução contínua (verossimilhança do BM, OU, EB,
# tendência, seleção e adequação de modelos, BBM/FPK/Lévy, vários regimes).
library(patchwork)
source(here::here("R/aula03_base.R"))
source(here::here("R/bm.R"))

# ---- Bloco 1: verossimilhança do BM ---------------------------------------

#' log L do BM para o atributo A da árvore-guia, em função de σ² (z₀ no MLE)
bm_verossimilhanca <- function() {
  tr <- arvore_guia(); a <- atributos_guia(); x <- setNames(a$A, a$sp)
  C <- ape::vcv(tr)[names(x), names(x)]; Ci <- solve(C); n <- length(x); um <- rep(1, n)
  z0 <- as.numeric((um %*% Ci %*% x) / (um %*% Ci %*% um))
  s2 <- as.numeric(t(x - z0) %*% Ci %*% (x - z0)) / n
  logL <- function(sig2) -.5 * (n * log(2 * pi) + n * log(sig2) + log(det(C)) +
                                  as.numeric(t(x - z0) %*% Ci %*% (x - z0)) / sig2)
  list(z0 = z0, s2 = s2, logL = logL, logLmax = logL(s2))
}
fig_bm_logL <- function(v = bm_verossimilhanca()) {
  g <- tibble(s2 = seq(.05, 1.2, length.out = 200)) |> mutate(L = purrr::map_dbl(s2, v$logL))
  ggplot(g, aes(s2, L)) +
    geom_line(colour = cores[["accent"]], linewidth = 1.3) +
    geom_vline(xintercept = v$s2, linetype = "22", colour = cores[["ferrugem"]]) +
    annotate("point", x = v$s2, y = v$logLmax, size = 4, colour = cores[["ferrugem"]]) +
    annotate("text", x = v$s2 + .04, y = v$logLmax - 1.2, hjust = 0, size = 6,
             colour = cores[["ferrugem"]], label = paste0("σ² = ", num(v$s2, 3))) +
    labs(x = expression(sigma^2), y = "log-verossimilhança") + tema_aula(17)
}

# ---- Bloco 2: Ornstein-Uhlenbeck ------------------------------------------

#' Simulações de OU variando um parâmetro (θ, α ou σ²)
fig_ou_param <- function(qual = "theta", n = 15, passos = 100) {
  base <- list(alfa = .1, theta = 8, sigma2 = .5)
  vals <- switch(qual, theta = c(4, 8), alfa = c(.03, .3), sigma2 = c(.1, 1))
  rot  <- switch(qual, theta = "θ = %s", alfa = "α = %s", sigma2 = "σ² = %s")
  d <- purrr::map_dfr(vals, \(v) {
    p <- base; p[[qual]] <- v
    sim_ou_tempo(n, passos, sigma2 = p$sigma2, alfa = p$alfa, theta = p$theta, seed = 501) |>
      mutate(painel = sprintf(rot, format(v, decimal.mark = ",")), theta = p$theta)
  })
  ggplot(d, aes(t, x, group = rep)) +
    geom_hline(aes(yintercept = theta), linetype = "22", colour = cores[["ferrugem"]]) +
    geom_line(colour = cores[["accent"]], alpha = .5, linewidth = .5) +
    facet_wrap(~painel) + labs(x = "tempo", y = "atributo") + tema_aula(16)
}

#' BM diversificante × OU de equilíbrio: variância entre réplicas no tempo
fig_bm_vs_ou <- function(reps = 1500, passos = 100, alfa = .1, s2 = 1) {
  v <- bind_rows(
    sim_bm_tempo(reps, passos, sigma2 = s2, seed = 502) |> mutate(modelo = "BM"),
    sim_ou_tempo(reps, passos, sigma2 = s2, alfa = alfa, theta = 0, seed = 502) |> mutate(modelo = "OU")) |>
    group_by(modelo, t) |> summarise(v = var(x), .groups = "drop")
  # variância estacionária do OU discreto: σ² / (1 − (1 − α)²) ≈ σ²/2α
  est <- s2 / (1 - (1 - alfa)^2)
  ggplot(v, aes(t, v, colour = modelo)) +
    geom_line(linewidth = 1.3) +
    geom_hline(yintercept = est, linetype = "22", colour = cores[["ferrugem"]]) +
    annotate("text", x = 5, y = est, vjust = -.6, hjust = 0, size = 5.5, colour = cores[["ferrugem"]],
             label = "OU: estabiliza em ≈ σ²/2α") +
    annotate("text", x = 70, y = 70, hjust = 1, size = 5.5, colour = cores[["accent"]],
             label = "BM: cresce como σ²t") +
    scale_colour_manual(values = c(BM = cores[["accent"]], OU = cores[["ferrugem"]]), name = NULL) +
    labs(x = "tempo", y = "variância entre réplicas") + tema_aula(16)
}

# ---- Bloco 3: EB e tendência ------------------------------------------------

fig_acdc <- function(n = 20) {
  d <- bind_rows(sim_eb_tempo(n, r = -.04, sigma2 = 4, seed = 503) |> mutate(p = "early burst (r < 0)"),
                 sim_eb_tempo(n, r = .02, sigma2 = .3, seed = 503) |> mutate(p = "late burst (r > 0)"))
  ggplot(d, aes(t, x, group = rep)) + geom_line(colour = cores[["accent"]], alpha = .55, linewidth = .5) +
    facet_wrap(~p) + labs(x = "tempo", y = "atributo") + tema_aula(16, conceitual = TRUE)
}

#' Variação dentro × entre subclados sob BM, OU e EB (ideia de Harmon et al. 2010)
fig_dentro_entre <- function(nrep = 200, seed = 504) {
  an <- anolis(); tr <- an$arvore
  cl <- cutree(as.hclust(phytools::force.ultrametric(tr, "extend", message = FALSE)), k = 8)[tr$tip.label]
  V <- list("BM" = ape::vcv(tr), "OU, α = 10" = vcv_ou_arvore(tr, 10), "EB, r = −5" = ape::vcv(eb_arvore(tr, -5)))
  set.seed(seed)
  d <- purrr::imap_dfr(V, \(v, nome) {
    X <- sim_mvn(v[tr$tip.label, tr$tip.label], nrep)
    tibble(modelo = nome, dentro = apply(X, 2, \(x) sum(tapply(x, cl, \(z) sum((z - mean(z))^2))) / sum((x - mean(x))^2)))
  }) |> mutate(modelo = factor(modelo, levels = names(V)))
  ggplot(d, aes(modelo, dentro, fill = modelo)) +
    geom_boxplot(width = .55, outlier.size = .7) +
    scale_fill_manual(values = cores_linhagem[c(1, 2, 3)], guide = "none") +
    scale_y_continuous(labels = \(v) paste0(100 * v, " %")) +
    labs(x = NULL, y = "variação dentro dos 8 maiores clados") + tema_aula(17)
}

fig_tendencia <- function(n = 10) {
  d <- bind_rows(sim_trend_tempo(n, mu = .25, seed = 505) |> mutate(p = "μ = 0,25"),
                 sim_trend_tempo(n, mu = -.25, seed = 506) |> mutate(p = "μ = −0,25"))
  ggplot(d, aes(t, x, group = interaction(rep, p), colour = p)) +
    geom_line(linewidth = .6, alpha = .8) +
    scale_colour_manual(values = c(cores[["accent"]], cores[["ferrugem"]]), name = NULL) +
    labs(x = "tempo", y = "atributo") + tema_aula(16, conceitual = TRUE)
}

# ---- Bloco 4: comparação de modelos nos Anolis ------------------------------

ajustes_anolis <- function(seed = 509) {
  an <- anolis(); tr <- an$arvore
  set.seed(seed)   # o fitContinuous sorteia pontos de partida
  fits <- lapply(c(BM = "BM", OU = "OU", EB = "EB"), \(m) suppressWarnings(geiger::fitContinuous(tr, an$svl, model = m)))
  tab <- purrr::imap_dfr(fits, \(f, m) tibble(modelo = m, k = f$opt$k, logL = f$opt$lnL, AICc = f$opt$aicc))
  tab$dAICc <- tab$AICc - min(tab$AICc)
  tab$peso <- exp(-.5 * tab$dAICc) / sum(exp(-.5 * tab$dAICc))
  list(tab = tab, fits = fits)
}

# ---- Bloco 5: adequação (ciclo do Arbutus) ----------------------------------

fig_arbutus <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo,
    0, 2, 3.4, 1, "1. Ajustar o modelo\naos dados",
    4, 2, 3.4, 1, "2. Transformar a árvore\npelo modelo ajustado",
    8, 2, 3.4, 1, "3. Estatísticas\nnos dados observados",
    8, 0, 3.4, 1, "4. Simular muitos conjuntos\nna árvore transformada",
    4, 0, 3.4, 1, "5. Estatísticas\nnos dados simulados",
    0, 0, 3.4, 1, "6. Comparar: os dados\nparecem simulados pelo modelo?")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend,
                        1.75, 2, 2.25, 2, 5.75, 2, 6.25, 2, 8, 1.45, 8, .55, 6.25, 0, 5.75, 0, 2.25, 0, 1.75, 0)
  diagrama(cx, st, base = 4.6)
}

# ---- Bloco 6: BBM, FPK e Lévy -------------------------------------------------

#' BM refletido entre limites (BBM)
fig_bbm <- function(n = 8, passos = 200, lim = c(-6, 6), seed = 507) {
  set.seed(seed)
  d <- purrr::map_dfr(seq_len(n), \(r) {
    x <- numeric(passos + 1)
    for (k in seq_len(passos)) {
      y <- x[k] + rnorm(1, 0, .6)
      if (y > lim[2]) y <- 2 * lim[2] - y
      if (y < lim[1]) y <- 2 * lim[1] - y
      x[k + 1] <- y
    }
    tibble(rep = r, t = 0:passos, x = x)
  })
  ggplot(d, aes(t, x, group = rep)) +
    geom_hline(yintercept = lim, colour = cores[["ferrugem"]], linewidth = 1.1) +
    geom_line(colour = cores[["accent"]], alpha = .5, linewidth = .5) +
    labs(x = "tempo", y = "atributo") + tema_aula(16, conceitual = TRUE)
}

#' Potenciais macroevolutivos: BM (plano), OU (poço), BBM (caixa), FPK (vários picos)
fig_potenciais <- function() {
  x <- seq(-1, 1, length.out = 400)
  d <- bind_rows(
    tibble(x, V = 0, m = "BM: plano"),
    tibble(x, V = 3 * x^2, m = "OU: um poço"),
    tibble(x, V = ifelse(abs(x) > .8, 3, 0), m = "BBM: uma caixa"),
    tibble(x, V = 8 * (x^4 - .5 * x^2) + .3 * x, m = "FPK: forma livre")) |>
    mutate(m = factor(m, levels = unique(m)))
  ggplot(d, aes(x, V)) +
    geom_line(colour = cores[["accent"]], linewidth = 1.4) +
    facet_wrap(~m, nrow = 1) +
    labs(x = "valor do atributo", y = "potencial") +
    tema_aula(15, conceitual = TRUE)
}

#' Evolução pulsada: BM com saltos raros (processo de Lévy)
fig_levy <- function(n = 12, passos = 200, seed = 508) {
  set.seed(seed)
  d <- purrr::map_dfr(seq_len(n), \(r) {
    inc <- rnorm(passos, 0, .3) + rbinom(passos, 1, .02) * rnorm(passos, 0, 5)
    tibble(rep = r, t = 0:passos, x = c(0, cumsum(inc)))
  })
  ggplot(d, aes(t, x, group = rep)) + geom_line(colour = cores[["accent"]], alpha = .6, linewidth = .5) +
    labs(x = "tempo", y = "atributo") + tema_aula(16, conceitual = TRUE)
}

# ---- Bloco 7: paisagem adaptativa -------------------------------------------

# ---- Bloco 8: várias taxas e vários regimes ----------------------------------

#' Árvore dos Anolis com os ecomorfos mapeados, reescalada para altura 1
anolis_simmap <- function() {
  data(anoletree, package = "phytools", envir = environment())
  data(anole.data, package = "phytools", envir = environment())
  tr <- anoletree; h <- max(phytools::nodeHeights(tr))
  tr$edge.length <- tr$edge.length / h
  tr$maps <- lapply(tr$maps, \(m) m / h); tr$mapped.edge <- tr$mapped.edge / h
  list(arvore = tr, svl = setNames(anole.data$SVL, rownames(anole.data))[tr$tip.label])
}
ajustes_regimes <- function() {
  a <- anolis_simmap()
  r <- list(BM1 = mvMORPH::mvBM(a$arvore, a$svl, model = "BM1", echo = FALSE, diagnostic = FALSE),
            BMM = mvMORPH::mvBM(a$arvore, a$svl, model = "BMM", echo = FALSE, diagnostic = FALSE),
            OU1 = mvMORPH::mvOU(a$arvore, a$svl, model = "OU1", echo = FALSE, diagnostic = FALSE, param = list(root = FALSE)),
            OUM = mvMORPH::mvOU(a$arvore, a$svl, model = "OUM", echo = FALSE, diagnostic = FALSE, param = list(root = FALSE)))
  tab <- purrr::imap_dfr(r, \(f, m) tibble(modelo = m, logL = f$LogLik, AICc = f$AICc))
  tab$peso <- exp(-.5 * (tab$AICc - min(tab$AICc))); tab$peso <- tab$peso / sum(tab$peso)
  list(tab = tab, oum = r$OUM, dados = a)
}
fig_otimos <- function(a = ajustes_regimes()) {
  th <- a$oum$theta; d <- tibble(eco = rownames(th), theta = th[, 1])
  pts <- tibble(sp = names(a$dados$svl), svl = a$dados$svl,
                eco = phytools::getStates(a$dados$arvore, "tips")[names(a$dados$svl)])
  ggplot(pts, aes(eco, svl, colour = eco)) +
    geom_jitter(width = .15, size = 2.8, alpha = .7) +
    geom_point(data = d, aes(eco, theta), shape = 23, size = 6, fill = cores[["texto"]], colour = "white") +
    scale_colour_manual(values = cores_linhagem[1:6], guide = "none") +
    labs(x = "ecomorfo", y = "log SVL") + tema_aula(16)
}
brownie_sunfish <- function() {
  data(sunfish.tree, package = "phytools", envir = environment())
  data(sunfish.data, package = "phytools", envir = environment())
  phytools::brownie.lite(sunfish.tree, setNames(sunfish.data$gape.width, rownames(sunfish.data)))
}
