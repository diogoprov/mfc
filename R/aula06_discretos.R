# Aula 6: atributos discretos (modelos de Markov, Mk, multiestado, taxas
# ocultas, correlação de Pagel, modelo de limiar).
library(patchwork)
source(here::here("R/aula03_base.R"))   # tema, árvore, matriz, diagrama
source(here::here("R/bm.R"))            # sim_bm_tempo (modelo de limiar)

# ---- Dados ------------------------------------------------------------------

#' Squamata (Brandley et al. 2008, dados da prática): árvore e estados
squamata <- function() {
  d <- read.csv(here::here("aula-02-bm-matriz-c/brandley_table.csv"))
  tr <- ape::read.nexus(here::here("aula-02-bm-matriz-c/squamate.tre"))
  sp <- sub(" ", "_", d$Species)
  hind <- d$HLL != 0; fore <- d$FLL != 0
  bin <- setNames(ifelse(!hind & !fore, "sem", "com"), sp)
  q4 <- setNames(ifelse(!hind & !fore, "sem membros", ifelse(fore & !hind, "só anteriores",
                 ifelse(!fore & hind, "só posteriores", "os quatro"))), sp)
  m <- suppressWarnings(geiger::treedata(tr, bin))
  tr <- m$phy
  list(arvore = tr, bin = bin[tr$tip.label], q4 = q4[tr$tip.label])
}

ajustes_squamata <- function() {
  s <- squamata()
  er  <- geiger::fitDiscrete(s$arvore, s$bin, model = "ER")
  ard <- geiger::fitDiscrete(s$arvore, s$bin, model = "ARD")
  m4 <- lapply(c(ER = "ER", SYM = "SYM", ARD = "ARD"), \(m) geiger::fitDiscrete(s$arvore, s$q4, model = m))
  tab4 <- purrr::imap_dfr(m4, \(f, m) tibble(modelo = m, k = f$opt$k, logL = f$opt$lnL, AICc = f$opt$aicc))
  tab4$peso <- exp(-.5 * (tab4$AICc - min(tab4$AICc))); tab4$peso <- tab4$peso / sum(tab4$peso)
  list(s = s, er = er$opt, ard = ard$opt, tab4 = tab4)
}

pagel_pratica <- function() {
  X <- read.csv(here::here("aula-05-modelos/fitPagel-data.csv"), row.names = 1)
  tr <- ape::read.tree(here::here("aula-05-modelos/fitPagel-tree.tre"))
  phytools::fitPagel(tr, setNames(X$x, rownames(X)), setNames(X$y, rownames(X)))
}

# ---- Bloco 1: do DNA ao Mk ----------------------------------------------------

#' Transições (purina ↔ purina, pirimidina ↔ pirimidina) e transversões
fig_dna <- function() {
  b <- tibble(x = c(0, 4, 0, 4), y = c(3, 3, 0, 0), l = c("A", "G", "C", "T"),
              tipo = c("purina", "purina", "pirimidina", "pirimidina"))
  seta <- arrow(length = unit(.25, "cm"), type = "closed", ends = "both")
  ggplot() +
    annotate("segment", x = .5, xend = 3.5, y = 3, yend = 3, arrow = seta, linewidth = 1.3, colour = cores[["accent"]]) +
    annotate("segment", x = .5, xend = 3.5, y = 0, yend = 0, arrow = seta, linewidth = 1.3, colour = cores[["accent"]]) +
    annotate("segment", x = 0, xend = 0, y = .5, yend = 2.5, arrow = seta, linewidth = 1.1, colour = cores[["ferrugem"]]) +
    annotate("segment", x = 4, xend = 4, y = .5, yend = 2.5, arrow = seta, linewidth = 1.1, colour = cores[["ferrugem"]]) +
    annotate("segment", x = .4, xend = 3.6, y = .4, yend = 2.6, arrow = seta, linewidth = 1.1, colour = cores[["ferrugem"]]) +
    annotate("segment", x = .4, xend = 3.6, y = 2.6, yend = .4, arrow = seta, linewidth = 1.1, colour = cores[["ferrugem"]]) +
    geom_point(data = b, aes(x, y, fill = tipo), shape = 21, size = 20, colour = "white") +
    geom_text(data = b, aes(x, y, label = l), size = 10, fontface = "bold", colour = "white") +
    annotate("text", x = 2, y = 3.45, label = "transição", colour = cores[["accent"]], size = 6) +
    annotate("text", x = 2, y = -.45, label = "transição", colour = cores[["accent"]], size = 6) +
    annotate("label", x = 2, y = 1.5, label = "transversões", colour = cores[["ferrugem"]], size = 6, fill = cores[["papel"]], border.colour = NA) +
    scale_fill_manual(values = c(purina = cores[["verde"]], pirimidina = cores[["ocre"]]), name = NULL) +
    coord_equal(xlim = c(-.8, 4.8), ylim = c(-.8, 3.8)) + theme_void() +
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA), legend.position = "right",
          legend.text = element_text(size = 15))
}

#' As três matrizes Q de O'Meara (2012): GTR, correlação binária, ordenada
fig_tres_Q <- function() {
  nm <- c("A", "C", "G", "T")
  gtr <- matrix(c("–", "rAC", "rAG", "rAT", "rAC", "–", "rCG", "rCT", "rAG", "rCG", "–", "rGT", "rAT", "rCT", "rGT", "–"), 4, byrow = TRUE, dimnames = list(nm, nm))
  nb <- c("00", "01", "10", "11")
  bin <- matrix(c("–", "q1", "q2", "0", "q3", "–", "0", "q4", "q5", "0", "–", "q6", "0", "q7", "q8", "–"), 4, byrow = TRUE, dimnames = list(nb, nb))
  no <- c("0", "1", "2", "3")
  ord <- matrix(c("–", "q01", "0", "0", "q10", "–", "q12", "0", "0", "q21", "–", "q23", "0", "0", "q32", "–"), 4, byrow = TRUE, dimnames = list(no, no))
  pinta <- function(M) matrix(ifelse(M == "0", "#e6e0d3", ifelse(M == "–", "#f3eee3", cores[["accent"]])), 4, 4)
  tt <- theme(plot.title = element_text(hjust = .5, size = 17, face = "bold", colour = cores[["texto"]]))
  (plot_matriz(gtr, fill = pinta(gtr), base = 5.2) + labs(title = "GTR (nucleotídeos)") + tt) +
    (plot_matriz(bin, fill = pinta(bin), base = 5.2) + labs(title = "dois caracteres binários") + tt) +
    (plot_matriz(ord, fill = pinta(ord), base = 5.2) + labs(title = "estados ordenados") + tt)
}

# ---- Bloco 2: Mk binário ------------------------------------------------------

#' Os dois estados e as setas (em camadas)
#' @param etapa 1 = só os estados; 2 = setas; 3 = taxas q01 e q10
fig_estados <- function(etapa = 3) {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    0, 0, 1.6, 1.4, "0\nsem membros", cores[["texto"]], cores[["texto"]], "white",
    5, 0, 1.6, 1.4, "1\ncom membros", cores[["papel"]], cores[["texto"]], cores[["texto"]])
  st <- if (etapa >= 2) tibble::tribble(~x, ~y, ~xend, ~yend, 1, .25, 4, .25, 4, -.25, 1, -.25) else NULL
  g <- diagrama(cx, st, base = 7.5)
  if (etapa >= 3) g <- g + annotate("text", x = 2.5, y = .7, label = 'q["01"]', parse = TRUE, size = 9, colour = cores[["accent"]]) +
    annotate("text", x = 2.5, y = -.7, label = 'q["10"]', parse = TRUE, size = 9, colour = cores[["ferrugem"]])
  g
}

#' Um ramo: o estado pode mudar ao longo dele (mudancas = TRUE: 0 | 1 | 0)
#' ou só sabemos o estado no início (1) e perguntamos o do fim (?)
fig_ramo <- function(mudancas = FALSE) {
  if (mudancas) {
    d <- tibble(x0 = c(0, 3, 7), x1 = c(3, 7, 10), e = c("0", "1", "0"))
    return(ggplot(d) +
      geom_rect(aes(xmin = x0, xmax = x1, ymin = 0, ymax = 1, fill = e), colour = cores[["texto"]]) +
      geom_text(aes((x0 + x1) / 2, .5, label = e, colour = e), size = 9) +
      scale_fill_manual(values = c("0" = cores[["texto"]], "1" = cores[["papel"]]), guide = "none") +
      scale_colour_manual(values = c("0" = "white", "1" = cores[["texto"]]), guide = "none") +
      annotate("text", x = 5, y = -.35, label = "tempo ao longo do ramo", size = 6) +
      coord_cartesian(ylim = c(-.6, 1.2)) + theme_void() +
      theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA)))
  }
  ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1, fill = "#d9d2c3", colour = cores[["texto"]]) +
    annotate("text", x = 0, y = 1.3, label = "1", size = 10, fontface = "bold") +
    annotate("text", x = 10, y = 1.3, label = "?", size = 10, fontface = "bold", colour = cores[["ferrugem"]]) +
    annotate("text", x = 5, y = -.35, label = "italic(t)", parse = TRUE, size = 9) +
    coord_cartesian(ylim = c(-.6, 1.6)) + theme_void() +
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}

#' Taxas iguais (setas da mesma espessura) × taxas desiguais (espessuras
#' diferentes), como em A2 115–116
fig_er_ard <- function(modelo = "ER") {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    0, 0, 1.6, 1.4, "0\nsem membros", cores[["texto"]], cores[["texto"]], "white",
    5, 0, 1.6, 1.4, "1\ncom membros", cores[["papel"]], cores[["texto"]], cores[["texto"]])
  lw <- if (modelo == "ER") c(1.6, 1.6) else c(2.8, .8)
  rot <- if (modelo == "ER") c("q", "q") else c('q["01"]', 'q["10"]')
  seta <- arrow(length = unit(.3, "cm"), type = "closed")
  diagrama(cx, NULL, base = 7.5) +
    annotate("segment", x = 1, xend = 4, y = .25, yend = .25, linewidth = lw[1], colour = cores[["accent"]], arrow = seta) +
    annotate("segment", x = 4, xend = 1, y = -.25, yend = -.25, linewidth = lw[2], colour = cores[["ferrugem"]], arrow = seta) +
    annotate("text", x = 2.5, y = c(.7, -.7), label = rot, parse = TRUE, size = 9,
             colour = c(cores[["accent"]], cores[["ferrugem"]])) +
    annotate("text", x = 2.5, y = 1.25, size = 8, fontface = "bold", colour = cores[["texto"]],
             label = if (modelo == "ER") "taxas iguais (ER)" else "taxas desiguais (ARD)")
}

#' P(t) para q01 = q10 = 0,1: probabilidade de continuar no estado ou mudar
fig_Pt <- function(q = .1, tmark = 5) {
  t <- seq(0, 30, length.out = 300)
  d <- bind_rows(tibble(t, p = .5 + .5 * exp(-2 * q * t), qual = "P11(t): continuar em 1"),
                 tibble(t, p = .5 - .5 * exp(-2 * q * t), qual = "P10(t): mudar para 0"))
  pm <- .5 + .5 * exp(-2 * q * tmark)
  ggplot(d, aes(t, p, colour = qual)) +
    geom_hline(yintercept = .5, linetype = "22", colour = cores[["neutro"]]) +
    geom_line(linewidth = 1.4) +
    geom_vline(xintercept = tmark, linetype = "22") +
    annotate("point", x = tmark, y = c(pm, 1 - pm), size = 4) +
    annotate("text", x = tmark + .6, y = c(pm, 1 - pm), hjust = 0, size = 5.5,
             label = c(num(pm), num(1 - pm))) +
    scale_colour_manual(values = c(cores[["ferrugem"]], cores[["accent"]]), name = NULL) +
    labs(x = "tempo (milhões de anos)", y = "probabilidade") + tema_aula(16)
}

#' Pontos coloridos nas pontas de uma árvore grande, sem nomes
arvore_estados <- function(tr, estados, cores_e, titulo = NULL, base = 4) {
  cd <- coord_arvore(tr); pt <- dplyr::filter(cd$nos, ponta)
  pt$e <- estados[pt$rotulo]
  xmax <- max(cd$nos$x)
  plot_arvore(tr, rotulos = FALSE, larg = .35, folga = .06) +
    geom_point(data = pt, aes(x + xmax * .02, y, colour = e), size = 1.2, shape = 15) +
    scale_colour_manual(values = cores_e, name = NULL) +
    labs(title = titulo) +
    theme(plot.title = element_text(hjust = .5, size = 15, face = "bold", colour = cores[["texto"]]),
          legend.position = "bottom", legend.text = element_text(size = 13))
}
cores_bin <- c(sem = cores[["ferrugem"]], com = cores[["accent"]])

fig_squamata <- function(s = squamata()) {
  arvore_estados(s$arvore, s$bin, cores_bin) + guides(colour = guide_legend(override.aes = list(size = 5)))
}

#' ER lento, ER rápido e taxas desiguais simulados na árvore dos Squamata
fig_mk_sim <- function(s = squamata(), seed = 601) {
  tr <- s$arvore; h <- max(phytools::nodeHeights(tr))
  Qm <- function(a, b) matrix(c(-a, a, b, -b), 2, byrow = TRUE, dimnames = list(c("com", "sem"), c("com", "sem")))
  set.seed(seed)
  sims <- list("ER lento" = Qm(.2 / h, .2 / h), "ER rápido" = Qm(5 / h, 5 / h), "taxas desiguais" = Qm(1 / h, 4 / h))
  gs <- purrr::imap(sims, \(Q, nome) {
    x <- phytools::sim.Mk(tr, Q, anc = "com")
    arvore_estados(tr, setNames(as.character(x), names(x)), cores_bin, nome) + theme(legend.position = "none")
  })
  gs[[1]] + gs[[2]] + gs[[3]]
}

# ---- Bloco 4: multiestado -------------------------------------------------------

fig_multi <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill,
    0, 2, 2.2, 1, "0\nsem membros", cores_linhagem[1],
    4, 2, 2.2, 1, "1\nsó anteriores", cores_linhagem[2],
    0, 0, 2.2, 1, "2\nsó posteriores", cores_linhagem[3],
    4, 0, 2.2, 1, "3\nos quatro", cores_linhagem[4]) |>
    mutate(cor = fill, txt = "white")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend,
                        1.2, 2.15, 2.8, 2.15, 2.8, 1.85, 1.2, 1.85, 1.2, .15, 2.8, .15, 2.8, -.15, 1.2, -.15,
                        -.15, 1.45, -.15, .55, .15, .55, .15, 1.45, 3.85, 1.45, 3.85, .55, 4.15, .55, 4.15, 1.45,
                        1.15, 1.5, 2.9, .5, 2.85, .5, 1.1, 1.5, 1.15, .5, 2.9, 1.5, 2.85, 1.5, 1.1, .5)
  diagrama(cx, st, base = 6.5)
}

#' As matrizes Q de 4 estados: ER, SYM, ARD (uma cor por parâmetro)
fig_Q4 <- function(modelo = "ER") {
  e <- as.character(0:3)
  M <- matrix(outer(0:3, 0:3, \(i, j) ifelse(i == j, "–", paste0("q", i, j))), 4, dimnames = list(e, e))
  par <- matrix(0, 4, 4)
  if (modelo == "ER") par[row(par) != col(par)] <- 1
  if (modelo == "SYM") { k <- 0; for (i in 1:3) for (j in (i + 1):4) { k <- k + 1; par[i, j] <- par[j, i] <- k } }
  if (modelo == "ARD") par[row(par) != col(par)] <- 1:12
  pal <- grDevices::hcl.colors(12, "Dark 3")
  fill <- matrix(ifelse(par == 0, "#f3eee3", pal[pmax(par, 1)]), 4, 4)
  if (modelo == "ER") fill[par == 1] <- cores[["accent"]]
  nparam <- c(ER = 1, SYM = 6, ARD = 12)[modelo]
  plot_matriz(M, fill = fill, base = 6.5) +
    labs(title = paste0(modelo, ": ", nparam, ifelse(nparam == 1, " parâmetro", " parâmetros"))) +
    theme(plot.title = element_text(hjust = .5, size = 20, face = "bold", colour = cores[["texto"]]))
}

# ---- Bloco 5: taxas ocultas -------------------------------------------------------

fig_ocultas <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    0, 2, 1.8, 1, "0A", cores[["texto"]], cores[["texto"]], "white",
    4, 2, 1.8, 1, "1A", cores[["papel"]], cores[["texto"]], cores[["texto"]],
    0, 0, 1.8, 1, "0B", cores[["texto"]], cores[["texto"]], "white",
    4, 0, 1.8, 1, "1B", cores[["papel"]], cores[["texto"]], cores[["texto"]])
  st <- tibble::tribble(~x, ~y, ~xend, ~yend,
                        1, 2.15, 3, 2.15, 3, 1.85, 1, 1.85, 1, .15, 3, .15, 3, -.15, 1, -.15,
                        -.15, 1.45, -.15, .55, .15, .55, .15, 1.45, 3.85, 1.45, 3.85, .55, 4.15, .55, 4.15, 1.45)
  diagrama(cx, st, base = 7) +
    annotate("text", x = 2, y = 2.55, label = "categoria A: taxas lentas", size = 5.5, colour = cores[["accent"]]) +
    annotate("text", x = 2, y = -.6, label = "categoria B: taxas rápidas", size = 5.5, colour = cores[["ferrugem"]]) +
    annotate("text", x = -1.4, y = 1, label = "troca de\ncategoria", size = 5, colour = cores[["neutro"]])
}

# ---- Bloco 6: Pagel ---------------------------------------------------------------

#' Os 4 estados combinados de dois caracteres binários
#' @param dependente FALSE = taxas de um caráter não dependem do outro (setas
#'   paralelas com a mesma cor); TRUE = todas as 8 taxas diferentes
fig_pagel <- function(dependente = FALSE) {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo,
    0, 2, 1.6, 1, "0,0", 4, 2, 1.6, 1, "0,1", 0, 0, 1.6, 1, "1,0", 4, 0, 1.6, 1, "1,1")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend, ~c,
                        .9, 2.15, 3.1, 2.15, "y01", 3.1, 1.85, .9, 1.85, "y10",
                        .9, .15, 3.1, .15, "y01", 3.1, -.15, .9, -.15, "y10",
                        -.15, 1.45, -.15, .55, "x01", .15, .55, .15, 1.45, "x10",
                        3.85, 1.45, 3.85, .55, "x01", 4.15, .55, 4.15, 1.45, "x10")
  pal <- if (dependente) grDevices::hcl.colors(8, "Dark 3") else
    c(y01 = cores_linhagem[1], y10 = cores_linhagem[2], x01 = cores_linhagem[3], x10 = cores_linhagem[4])[st$c]
  g <- diagrama(cx, NULL, base = 6.5)
  for (i in seq_len(nrow(st))) g <- g + annotate("segment", x = st$x[i], y = st$y[i], xend = st$xend[i], yend = st$yend[i],
                                                linewidth = 1.6, colour = pal[i], arrow = arrow(length = unit(.25, "cm"), type = "closed"))
  g + labs(title = if (dependente) "dependente: 8 taxas" else "independente: 4 taxas") +
    theme(plot.title = element_text(hjust = .5, size = 18, face = "bold", colour = cores[["texto"]]))
}

#' Pseudorreplicação: uma única mudança conjunta, no ramo que leva a um clado,
#' deixa X = 1 e Y = 1 em todas as espécies do clado
fig_pseudorep <- function(n = 40, seed = 602) {
  set.seed(seed)
  tr <- phytools::pbtree(n = n, scale = 1)
  nos <- (n + 2):(n + tr$Nnode)
  ndesc <- sapply(nos, \(k) sum(phytools::getDescendants(tr, k) <= n))
  no <- nos[which.min(abs(ndesc - n / 3))]
  desc <- phytools::getDescendants(tr, no)
  cor <- ifelse(tr$edge[, 2] %in% desc, cores[["ferrugem"]], cores[["neutro"]])
  cd <- coord_arvore(tr); pt <- dplyr::filter(cd$nos, ponta)
  pt$um <- pt$no %in% desc
  xmax <- max(cd$nos$x)
  tronco <- dplyr::filter(cd$ramos, filho == no)
  col <- tibble(xc = xmax * c(1.06, 1.14), l = c("X", "Y"))
  pc <- tidyr::crossing(pt, col)
  plot_arvore(tr, cor_ramo = cor, larg = .7, rotulos = FALSE, folga = .25) +
    geom_point(data = pc, aes(xc, y), shape = 15, size = 3.2,
               colour = ifelse(pc$um, cores[["ferrugem"]], "#cfc8b8")) +
    annotate("text", x = col$xc, y = max(pt$y) + 1.6, label = col$l, size = 6, fontface = "bold") +
    annotate("point", x = (tronco$x0 + tronco$x1) / 2, y = tronco$y, shape = 8, size = 6, stroke = 1.5,
             colour = cores[["ferrugem"]]) +
    expand_limits(y = max(pt$y) + 2.5) +
    annotate("label", x = tronco$x0, y = tronco$y - 2.2, hjust = .3, size = 5.5, fill = cores[["papel"]],
             border.colour = NA, colour = cores[["ferrugem"]], label = "X e Y mudam juntos, uma vez")
}

# ---- Bloco 7: limiar ----------------------------------------------------------------

#' A liability (contínua, não observada) e o limiar que define o estado
fig_limiar <- function(n = 6, passos = 400, seed = 603) {
  d <- sim_bm_tempo(n, passos, sigma2 = .02, seed = seed)
  d$estado <- ifelse(d$x > 0, "1", "0")
  ggplot(d, aes(t, x, group = rep)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf, fill = cores[["accent"]], alpha = .12) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = cores[["ferrugem"]], alpha = .12) +
    geom_hline(yintercept = 0, linewidth = 1.2, colour = cores[["texto"]]) +
    geom_line(colour = cores[["texto"]], alpha = .75, linewidth = .6) +
    annotate("label", x = passos, y = max(d$x), hjust = 1, vjust = 1, size = 6, fill = "white", border.colour = NA, colour = cores[["accent"]], label = "estado 1") +
    annotate("label", x = passos, y = min(d$x), hjust = 1, vjust = 0, size = 6, fill = "white", border.colour = NA, colour = cores[["ferrugem"]], label = "estado 0") +
    annotate("label", x = 0, y = 0, hjust = 0, size = 5.5, label = "limiar", fill = "white", border.colour = NA) +
    labs(x = "tempo", y = "liability") + tema_aula(16, conceitual = TRUE)
}

#' Teste de correlação com taxas ocultas (Boyko & Beaulieu 2023) nos dados
#' da prática: modelos independente e correlacionado, com e sem uma segunda
#' classe de taxas (corHMM::fitCorrelationTest; ~3 min)
correlacao_corhmm <- function() {
  X <- read.csv(here::here("aula-05-modelos/fitPagel-data.csv"), row.names = 1)
  tr <- ape::read.tree(here::here("aula-05-modelos/fitPagel-tree.tre"))
  r <- corHMM::fitCorrelationTest(tr, data.frame(sp = rownames(X), x = X$x, y = X$y))
  nomes <- c(independent_model_fit = "independente",
             hidden_Markov_independent_model_fit = "independente + taxas ocultas",
             correlated_model_fit = "correlacionado",
             hidden_Markov_correlated_model_fit = "correlacionado + taxas ocultas")
  tab <- purrr::imap_dfr(r, \(f, m) tibble(modelo = nomes[[m]], logL = f$loglik, AICc = f$AICc))
  tab$peso <- exp(-.5 * (tab$AICc - min(tab$AICc))); tab$peso <- tab$peso / sum(tab$peso)
  tab
}
