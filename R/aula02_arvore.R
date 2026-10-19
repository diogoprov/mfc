# Aula 2, bloco 4: o BM numa árvore.
library(patchwork)
library(nlme)
source(here::here("R/tema.R"))
source(here::here("R/arvore.R"))
source(here::here("R/matriz.R"))
source(here::here("R/aula02_intro.R"))   # arvore_dois_clados()

#' Filogenia em estrela com a mesma altura e as mesmas pontas da árvore-guia
arvore_estrela <- function() {
  tr <- arvore_guia()
  ape::stree(7, tip.label = tr$tip.label) |> (\(s) { s$edge.length <- rep(1.44, 7); s })()
}

#' Estrela × árvore-guia: a árvore em cima, o fenograma do BM embaixo
fig_estrela_vs_arvore <- function(seed = 8) {
  st <- arvore_estrela(); tr <- arvore_guia()
  cor_ponta <- function(t) {
    cr <- rep(cores[["texto"]], nrow(t$edge))
    term <- t$edge[, 2] <= 7
    cr[term] <- cores_linhagem[t$edge[term, 2]]
    cr
  }
  a1 <- plot_arvore(st, cor_ramo = cor_ponta(st), base = 5) + labs(title = "Filogenia em estrela") +
    theme(plot.title = element_text(hjust = .5, size = 18, face = "bold", colour = cores[["texto"]]))
  a2 <- plot_arvore(tr, cor_ramo = cor_ponta(tr), base = 5) + labs(title = "Filogenia “real”") +
    theme(plot.title = element_text(hjust = .5, size = 18, face = "bold", colour = cores[["texto"]]))
  lim <- NULL
  f1 <- fenograma(sim_bm_arvore(st, seed = seed), st, cor_ponta(st), rotulos = FALSE)
  f2 <- fenograma(sim_bm_arvore(tr, seed = seed), tr, cor_ponta(tr), rotulos = FALSE)
  lim <- range(c(layer_scales(f1)$y$range$range, layer_scales(f2)$y$range$range))
  (a1 | a2) / ((f1 + ylim(lim)) | (f2 + ylim(lim))) + plot_layout(heights = c(1, 1.3))
}

#' 1000 réplicas de BM na árvore-guia (σ² = 1, raiz = 0)
sims_guia <- function(nsim = 1000, seed = 81) {
  set.seed(seed)
  tr <- arvore_guia()
  t(replicate(nsim, ape::rTraitCont(tr, "BM", sigma = 1)[tr$tip.label]))
}

#' Pares de espécies nas 1000 réplicas: o parentesco cria covariação
fig_pares <- function() {
  X <- sims_guia()
  pares <- list(c(3, 4), c(1, 2), c(1, 3), c(1, 5))
  d <- purrr::map_dfr(pares, \(p) tibble(
    x = X[, p[1]], y = X[, p[2]],
    par = sprintf("sp%d × sp%d\nC = %s · r = %s", p[1], p[2],
                  format(ape::vcv(arvore_guia())[p[1], p[2]], nsmall = 2, decimal.mark = ","),
                  format(round(cor(X[, p[1]], X[, p[2]]), 2), nsmall = 2, decimal.mark = ","))))
  d$par <- factor(d$par, levels = unique(d$par))
  ggplot(d, aes(x, y)) +
    geom_point(colour = cores[["accent"]], alpha = .25, size = 1.2) +
    facet_wrap(~par, nrow = 1) + coord_equal() +
    labs(x = "valor da primeira espécie", y = "valor da segunda") + tema_aula(14)
}

#' Covariância das réplicas × matriz C
fig_cov_sim <- function() {
  tr <- arvore_guia()
  S <- cov(sims_guia()); dimnames(S) <- dimnames(ape::vcv(tr))
  m1 <- plot_matriz(ape::vcv(tr), base = 4.3) + labs(title = "C (da árvore)")
  m2 <- plot_matriz(S, base = 4.3) + labs(title = "covariância de 1000 réplicas")
  tt <- theme(plot.title = element_text(hjust = .5, size = 18, face = "bold", colour = cores[["texto"]]))
  (m1 + tt) + (m2 + tt)
}

#' Erro tipo I: dois atributos independentes, OLS × PGLS
#'
#' @return tibble com o método e a proporção de p < 0,05
erro_tipo1 <- function(nsim = 1000, seed = 91) {
  tr <- arvore_dois_clados()
  set.seed(seed)
  p <- t(replicate(nsim, {
    d <- data.frame(sp = tr$tip.label,
                    x = ape::rTraitCont(tr, "BM")[tr$tip.label],
                    y = ape::rTraitCont(tr, "BM")[tr$tip.label])
    c(OLS  = summary(lm(y ~ x, d))$coefficients[2, 4],
      PGLS = summary(nlme::gls(y ~ x, d, correlation = ape::corBrownian(1, tr, form = ~sp)))$tTable[2, 4])
  }))
  tibble(metodo = c("OLS (espécies\nindependentes)", "PGLS (com C)"),
         taxa = colMeans(p < .05))
}

fig_erro_tipo1 <- function(res = erro_tipo1()) {
  arv <- plot_arvore(arvore_dois_clados(), rotulos = FALSE, larg = 1.1) +
    labs(title = "20 espécies em dois clados") +
    theme(plot.title = element_text(hjust = .5, size = 17, colour = cores[["texto"]]))
  barras <- ggplot(res, aes(metodo, taxa)) +
    geom_col(fill = c(cores[["ferrugem"]], cores[["accent"]]), width = .6) +
    geom_hline(yintercept = .05, linetype = "22", linewidth = 1, colour = cores[["texto"]]) +
    annotate("text", x = 2.45, y = .05, label = "5 %", vjust = -.6, size = 5.5) +
    geom_text(aes(label = paste0(format(round(100 * taxa, 1), decimal.mark = ","), " %")),
              vjust = -.5, size = 6.5, fontface = "bold") +
    scale_y_continuous(labels = \(x) paste0(100 * x, " %"), limits = c(0, max(res$taxa) * 1.15)) +
    labs(x = NULL, y = "simulações com p < 0,05") + tema_aula(17)
  arv + barras + plot_layout(widths = c(1, 1.5))
}
