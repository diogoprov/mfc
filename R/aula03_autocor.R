# Aula 3, blocos 5–6: autocorrelação filogenética (Moran, Abouheif), Mantel
# e distogramas.
source(here::here("R/aula03_base.R"))

#' A matriz W (inverso da distância patrística) da árvore-guia
fig_W <- function() {
  tr <- arvore_guia()
  D <- ape::cophenetic.phylo(tr)[paste0("sp", 1:7), paste0("sp", 1:7)]
  W <- 1 / D; diag(W) <- 0
  plot_arvore(tr, base = 5.5) + plot_matriz(W, base = 4.6) + plot_layout(widths = c(1, 1.3))
}

#' Correlogramas médios sob BM e sob OU (α = 10), 50 réplicas na árvore dos Anolis
fig_correlograma_modelos <- function(nrep = 50, seed = 311) {
  an <- anolis(); tr <- an$arvore
  set.seed(seed)
  V <- list("BM" = ape::vcv(tr), "OU, α = 10" = vcv_ou_arvore(tr, 10))
  d <- purrr::imap_dfr(V, \(v, nome) {
    X <- sim_mvn(v[tr$tip.label, tr$tip.label], nrep)
    purrr::map_dfr(seq_len(nrep), \(r) correlograma(tr, X[, r], 4) |> mutate(rep = r, modelo = nome))
  })
  med <- d |> group_by(modelo, classe) |> summarise(d_media = mean(d_media), I = mean(I), .groups = "drop")
  ggplot(d, aes(d_media, I, colour = modelo)) +
    geom_hline(yintercept = 0, colour = cores[["neutro"]]) +
    geom_line(aes(group = interaction(rep, modelo)), alpha = .12) +
    geom_line(data = med, linewidth = 1.6) + geom_point(data = med, size = 4) +
    scale_colour_manual(values = cores_linhagem[c(1, 2)], name = NULL) +
    labs(x = "distância filogenética (classe)", y = "I de Moran") + tema_aula(17)
}

#' Correlograma do SVL dos Anolis
#' @param etapa 1 = histograma das distâncias com as 4 classes; 2 = correlograma
fig_correlograma_anolis <- function(etapa = 2) {
  an <- anolis(); tr <- an$arvore
  D <- ape::cophenetic.phylo(tr)
  cortes <- seq(0, max(D), length.out = 5)
  if (etapa == 1) {
    return(ggplot(tibble(d = D[upper.tri(D)]), aes(d)) +
             geom_histogram(bins = 50, fill = cores[["neutro"]], colour = cores[["papel"]]) +
             geom_vline(xintercept = cortes[2:4], linetype = "22", colour = cores[["ferrugem"]], linewidth = 1) +
             annotate("text", x = (cortes[-5] + cortes[-1]) / 2, y = Inf, vjust = 1.4,
                      label = paste("classe", 1:4), size = 5.5, colour = cores[["ferrugem"]]) +
             scale_y_sqrt() +
             labs(x = "distância patrística entre pares de espécies", y = "pares (escala raiz)") + tema_aula(16))
  }
  cg <- correlograma(tr, an$svl, 4)
  ggplot(cg, aes(d_media, I)) +
    geom_hline(yintercept = 0, colour = cores[["neutro"]]) +
    geom_line(colour = cores[["accent"]], linewidth = 1.4) +
    geom_point(size = 5, colour = cores[["accent"]]) +
    geom_text(aes(label = paste(pares, "pares")), vjust = -1.2, size = 5) +
    scale_y_continuous(expand = expansion(mult = c(.1, .2))) +
    scale_x_continuous(expand = expansion(mult = .1)) +
    labs(x = "distância filogenética (centro da classe)", y = "I de Moran do SVL") + tema_aula(17)
}

#' Exemplo de Abouheif (1999): massa corporal na ordem das pontas
abouheif_dados <- function() {
  y <- c(53.6, 51.2, 49.2, 47.5, 10.2, 11.1, 13.2, 12.6)
  list(y = y, d2 = diff(y)^2, sq = sum((y - mean(y))^2),
       eta = sum(diff(y)^2) / sum((y - mean(y))^2),
       C = 1 - sum(diff(y)^2) / (2 * sum((y - mean(y))^2)))
}
fig_abouheif <- function(mostrar_d2 = TRUE) {
  a <- abouheif_dados()
  d <- tibble(i = seq_along(a$y), y = a$y)
  g <- ggplot(d, aes(i, y)) +
    geom_col(fill = cores[["accent"]], width = .65) +
    geom_text(aes(label = format(y, decimal.mark = ",")), vjust = -.5, size = 5.5) +
    scale_x_continuous(breaks = 1:8, labels = paste("ponta", 1:8)) +
    scale_y_continuous(expand = expansion(mult = c(0, .25))) +
    labs(x = "pontas, na ordem da árvore", y = "massa corporal (kg)") + tema_aula(16)
  if (mostrar_d2) {
    dd <- tibble(x = 1:7 + .5, d2 = a$d2, y = pmax(a$y[-1], a$y[-8]) + 9)
    g <- g + geom_label(data = dd, aes(x, y, label = format(round(d2, 2), decimal.mark = ",")),
                        size = 4.6, fill = cores[["papel"]], colour = cores[["ferrugem"]], linewidth = .3)
  }
  g
}

#' Teste de Mantel na árvore-guia: distância patrística × distância no
#' atributo A, com a distribuição nula por permutação
mantel_guia <- function(nperm = 999, seed = 321) {
  tr <- arvore_guia(); a <- atributos_guia(); x <- setNames(a$A, a$sp)
  Df <- ape::cophenetic.phylo(tr)[names(x), names(x)]
  Da <- as.matrix(dist(x))
  r <- function(P, A) cor(P[upper.tri(P)], A[upper.tri(A)])
  set.seed(seed)
  nulo <- replicate(nperm, { o <- sample(names(x)); A2 <- Da[o, o]; dimnames(A2) <- dimnames(Da); r(Df, A2) })
  list(Df = Df, Da = Da, r = r(Df, Da), nulo = nulo, p = (sum(nulo >= r(Df, Da)) + 1) / (nperm + 1))
}
fig_mantel <- function(m = mantel_guia()) {
  m1 <- plot_matriz(m$Df, base = 4.4, digitos = 1) + labs(title = "distância patrística")
  m2 <- plot_matriz(m$Da, base = 4.4, digitos = 1) + labs(title = "distância em A")
  tt <- theme(plot.title = element_text(hjust = .5, size = 16, face = "bold", colour = cores[["texto"]]))
  h <- ggplot(tibble(r = m$nulo), aes(r)) +
    geom_histogram(bins = 30, fill = cores[["neutro"]], colour = cores[["papel"]]) +
    geom_vline(xintercept = m$r, colour = cores[["ferrugem"]], linewidth = 1.3) +
    labs(x = "r de Mantel (999 permutações)", y = NULL) + tema_aula(14)
  (m1 + tt) + (m2 + tt) + h + plot_layout(widths = c(1, 1, .8))
}

#' Distograma: semivariância do atributo por classe de distância patrística,
#' com a reta esperada sob BM (γ = σ² d / 2)
fig_distograma <- function(n_classes = 8) {
  an <- anolis(); tr <- an$arvore; x <- an$svl
  D <- ape::cophenetic.phylo(tr)[names(x), names(x)]
  G <- outer(x, x, \(a, b) (a - b)^2 / 2)
  d <- tibble(dist = D[upper.tri(D)], gama = G[upper.tri(G)]) |>
    mutate(classe = cut(dist, n_classes)) |>
    group_by(classe) |> summarise(dist = mean(dist), gama = mean(gama), n = n(), .groups = "drop")
  C <- ape::vcv(tr)[names(x), names(x)]; Ci <- solve(C); um <- rep(1, length(x))
  a <- as.numeric((um %*% Ci %*% x) / (um %*% Ci %*% um))
  sigma2 <- as.numeric(t(x - a) %*% Ci %*% (x - a)) / length(x)   # MLE de σ² sob BM
  ggplot(d, aes(dist, gama)) +
    geom_abline(slope = sigma2 / 2, intercept = 0, linetype = "22", colour = cores[["ferrugem"]], linewidth = 1.1) +
    geom_point(aes(size = n), colour = cores[["accent"]]) +
    annotate("text", x = .25, y = sigma2 / 2 * 1.6, hjust = 0, colour = cores[["ferrugem"]], size = 5.5,
             label = "esperado sob BM") +
    scale_size_area(max_size = 9, name = "pares") +
    labs(x = "distância patrística", y = "semivariância do SVL") + tema_aula(16) +
    theme(legend.position = "right")
}
