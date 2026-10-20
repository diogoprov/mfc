# Aula 3, blocos 1–4: o que é sinal, K de Blomberg e λ de Pagel.
source(here::here("R/aula03_base.R"))

#' Barras do atributo ao lado das pontas de uma árvore
arvore_barras <- function(tr, valores, titulo = NULL, cor = cores[["accent"]], base = 5) {
  cd <- coord_arvore(tr)
  pt <- dplyr::filter(cd$nos, ponta)
  pt$v <- valores[pt$rotulo]
  xmax <- max(cd$nos$x); esc <- .6 * xmax / max(abs(valores))
  x0 <- xmax * 1.35
  plot_arvore(tr, base = base, folga = .9) +
    geom_segment(data = pt, aes(x = x0, xend = x0 + v * esc, y = y, yend = y),
                 linewidth = 7, colour = cor) +
    annotate("segment", x = x0, xend = x0, y = .5, yend = max(pt$y) + .5,
             colour = cores[["neutro"]]) +
    labs(title = titulo) +
    theme(plot.title = element_text(hjust = .5, size = 18, face = "bold", colour = cores[["texto"]]))
}

#' A definição na árvore-guia: um atributo com sinal (BM) e o mesmo atributo
#' embaralhado entre as pontas (sem sinal)
fig_definicao <- function(seed = 4) {
  tr <- arvore_guia()
  set.seed(seed)
  x <- ape::rTraitCont(tr, "BM", sigma = 1.5)[tr$tip.label]
  x <- x - min(x) + .3
  y <- setNames(sample(x), names(x))
  arvore_barras(tr, x, "Com sinal (BM)") + arvore_barras(tr, y, "Sem sinal (embaralhado)", cor = cores[["ferrugem"]])
}

#' Os dois sentidos de "sinal" (slide 26 do original)
fig_dois_sentidos <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    0, 1, 2.6, 1, "\u201cSinal\u201d\nfilogenético", cores[["accent"]], cores[["accent"]], "white",
    4.6, 2, 4.4, 1.4, "Inferência filogenética:\nhabilidade de um caráter de\nreconstruir a história das espécies", cores[["papel"]], cores[["accent_dark"]], cores[["accent_dark"]],
    4.6, 0, 4.4, 1.4, "MFC, ecofilogenética e macroevolução:\nespécies aparentadas são\nmais parecidas", cores[["papel"]], cores[["accent_dark"]], cores[["accent_dark"]])
  st <- tibble::tribble(~x, ~y, ~xend, ~yend, 1.35, 1.25, 2.35, 1.9, 1.35, .75, 2.35, .1)
  diagrama(cx, st, base = 5.6)
}

#' K e λ simulados sob BM, OU e EB na árvore dos Anolis (altura 1)
sim_sinal_modelos <- function(nrep = 100, seed = 301) {
  an <- anolis(); tr <- an$arvore
  modelos <- list(
    "BM" = ape::vcv(tr),
    "OU, α = 2" = vcv_ou_arvore(tr, 2),
    "OU, α = 10" = vcv_ou_arvore(tr, 10),
    "EB, r = −5" = ape::vcv(eb_arvore(tr, -5)))
  set.seed(seed)
  purrr::imap_dfr(modelos, \(V, nome) {
    X <- sim_mvn(V[tr$tip.label, tr$tip.label], nrep)
    purrr::map_dfr(seq_len(nrep), \(r) {
      x <- X[, r]
      tibble(modelo = nome,
             K = as.numeric(phytools::phylosig(tr, x, method = "K")),
             lambda = phytools::phylosig(tr, x, method = "lambda")$lambda)
    })
  }) |> mutate(modelo = factor(modelo, levels = names(modelos)))
}

#' Distribuições de K e de λ por modelo
#' @param destaque modelos a destacar (os demais esmaecidos); NULL = todos
fig_sinal_modelos <- function(d = sim_sinal_modelos(), destaque = NULL) {
  dl <- tidyr::pivot_longer(d, c(K, lambda), names_to = "medida", values_to = "valor") |>
    mutate(medida = factor(medida, levels = c("K", "lambda"), labels = c("K de Blomberg", "λ de Pagel")),
           on = if (is.null(destaque)) TRUE else modelo %in% destaque)
  ggplot(dl, aes(modelo, valor, fill = modelo, alpha = on)) +
    geom_hline(data = tibble(medida = factor(c("K de Blomberg", "λ de Pagel")), y = 1),
               aes(yintercept = y), linetype = "22", colour = cores[["texto"]]) +
    geom_boxplot(outlier.size = .8, width = .6) +
    facet_wrap(~medida, scales = "free_y") +
    scale_fill_manual(values = cores_linhagem[c(1, 4, 2, 3)], guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = .9, `FALSE` = .2), guide = "none") +
    labs(x = NULL, y = NULL) + tema_aula(16) +
    theme(panel.grid.major.x = element_blank())
}

#' Régua de leitura do K
fig_regua_k <- function() {
  marcas <- tibble(x = c(0, 1, 2), r = c("K = 0", "K = 1", "K > 1"),
                   t = c("sem sinal:\natributo independente\nda filogenia",
                         "o esperado\nsob BM",
                         "parentes mais parecidos\ndo que o esperado sob BM"))
  ggplot() +
    annotate("rect", xmin = -.2, xmax = 1, ymin = -.15, ymax = .15, fill = cores[["ferrugem"]], alpha = .25) +
    annotate("rect", xmin = 1, xmax = 2.3, ymin = -.15, ymax = .15, fill = cores[["accent"]], alpha = .25) +
    annotate("segment", x = -.2, xend = 2.35, y = 0, yend = 0, linewidth = 1.2, colour = cores[["texto"]],
             arrow = arrow(length = unit(.3, "cm"), type = "closed")) +
    geom_point(data = marcas, aes(x, 0), size = 5, colour = cores[["texto"]]) +
    geom_text(data = marcas, aes(x, .35, label = r), size = 7, fontface = "bold", colour = cores[["texto"]]) +
    geom_text(data = marcas, aes(x, -.5, label = t), size = 5, lineheight = .95, colour = cores[["texto"]]) +
    annotate("text", x = .4, y = .95, label = "K < 1: menos parecidos que o esperado\n(homoplasia, adaptação não correlacionada\ncom a filogenia ou erro de medida)",
             size = 4.8, colour = cores[["ferrugem"]], lineheight = .95) +
    coord_cartesian(xlim = c(-.35, 2.5), ylim = c(-.95, 1.25)) + theme_void() +
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}

#' K observado do SVL dos Anolis e a distribuição nula por permutação
k_nulo <- function(nperm = 999, seed = 302) {
  an <- anolis(); tr <- an$arvore; x <- an$svl
  set.seed(seed)
  list(obs = as.numeric(phytools::phylosig(tr, x, method = "K")),
       nulo = replicate(nperm, as.numeric(phytools::phylosig(tr, setNames(sample(x), names(x)), method = "K"))))
}
fig_k_nulo <- function(kn = k_nulo()) {
  ggplot(tibble(K = kn$nulo), aes(K)) +
    geom_histogram(bins = 40, fill = cores[["neutro"]], colour = cores[["papel"]]) +
    geom_vline(xintercept = kn$obs, colour = cores[["ferrugem"]], linewidth = 1.4) +
    annotate("text", x = kn$obs, y = Inf, vjust = 1.5, hjust = 1.05, size = 6, colour = cores[["ferrugem"]],
             label = paste0("K observado = ", num(kn$obs))) +
    labs(x = "K com as pontas embaralhadas (999 permutações)", y = "frequência") + tema_aula(16)
}

#' Poder dos testes de K e de λ sob BM, por número de espécies
poder_sinal <- function(ns = c(10, 20, 50, 100), nrep = 100, seed = 303) {
  set.seed(seed)
  purrr::map_dfr(ns, \(n) purrr::map_dfr(seq_len(nrep), \(r) {
    tr <- phytools::pbtree(n = n, scale = 1)
    x <- phytools::fastBM(tr)
    tibble(n = n,
           K = phytools::phylosig(tr, x, method = "K", test = TRUE, nsim = 199)$P < .05,
           lambda = phytools::phylosig(tr, x, method = "lambda", test = TRUE)$P < .05)
  })) |>
    group_by(n) |> summarise(K = mean(K), lambda = mean(lambda), .groups = "drop")
}
fig_poder <- function(p = poder_sinal()) {
  dl <- tidyr::pivot_longer(p, c(K, lambda), names_to = "teste", values_to = "poder") |>
    mutate(teste = ifelse(teste == "K", "K (permutação)", "λ (razão de verossimilhança)"))
  ggplot(dl, aes(n, poder, colour = teste)) +
    geom_line(linewidth = 1.2) + geom_point(size = 4) +
    scale_x_continuous(breaks = p$n) +
    scale_y_continuous(labels = \(v) paste0(100 * v, " %"), limits = c(0, 1)) +
    scale_colour_manual(values = cores_linhagem[c(1, 2)], name = NULL) +
    labs(x = "número de espécies", y = "testes significativos (BM verdadeiro)") + tema_aula(17)
}

#' O exemplo do λ: árvore de 5 espécies e a matriz, com λ = 1, 0,5 ou 0,05
#' @param marcar destacar as células multiplicadas por λ (com "λ ·")
fig_lambda <- function(lam = 1, marcar = FALSE) {
  tr <- arvore_lambda5()
  ord <- paste0("sp", 1:5)
  C <- ape::vcv(tr)[ord, ord]
  M <- C; M[upper.tri(M) | lower.tri(M)] <- lam * C[upper.tri(C) | lower.tri(C)]
  dest <- (M != 0) & !diag(5)
  if (marcar) {
    M2 <- matrix(format(M, decimal.mark = ","), 5, 5, dimnames = dimnames(M))
    M2[dest] <- paste0("λ·", C[dest])
    mat <- plot_matriz(M2, destaque = dest, base = 6.5)
  } else {
    mat <- plot_matriz(M, destaque = if (lam < 1) dest else NULL,
                       digitos = if (lam == .05) 2 else if (lam == .5) 1 else 0, base = 6.5)
  }
  arv <- plot_arvore(lambda_arvore(tr, lam), comprimentos = lam == 1, base = 6.5, dig = 0) +
    labs(title = paste0("λ = ", format(lam, decimal.mark = ","))) +
    theme(plot.title = element_text(hjust = .5, size = 20, face = "bold", colour = cores[["texto"]]))
  arv + mat + plot_layout(widths = c(1, 1.1))
}

#' Perfil da verossimilhança de λ para o SVL dos Anolis
fig_perfil_lambda <- function() {
  an <- anolis()
  pl <- phytools::phylosig(an$arvore, an$svl, method = "lambda", test = TRUE)
  g <- tibble(lambda = seq(0, 1, by = .01)) |> mutate(logL = purrr::map_dbl(lambda, pl$lik))
  ggplot(g, aes(lambda, logL)) +
    geom_line(colour = cores[["accent"]], linewidth = 1.3) +
    geom_vline(xintercept = pl$lambda, linetype = "22", colour = cores[["ferrugem"]]) +
    annotate("point", x = 0, y = pl$lik(0), size = 4, colour = cores[["texto"]]) +
    annotate("text", x = .03, y = pl$lik(0), hjust = 0, size = 5.5, label = "λ = 0 (estrela)") +
    annotate("text", x = pl$lambda - .02, y = max(g$logL) - 6, hjust = 1, size = 5.5,
             colour = cores[["ferrugem"]], label = paste0("λ estimado = ", num(pl$lambda))) +
    labs(x = "λ", y = "log-verossimilhança") + tema_aula(17)
}
