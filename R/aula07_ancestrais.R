# Aula 7: estados ancestrais (contínuos e discretos) e visualização.
library(patchwork)
source(here::here("R/aula03_base.R"))   # tema, árvore, matriz, diagrama, anolis()
source(here::here("R/bm.R"))

fundo <- function() par(bg = cores[["papel"]], fg = cores[["texto"]], col.axis = cores[["texto"]],
                        col.lab = cores[["texto"]], family = "sans")

# ---- Dados ------------------------------------------------------------------

#' Árvore-guia e atributo A: estimativas por ML (BM), com e sem comprimentos
guia_ancestrais <- function() {
  tr <- arvore_guia(); a <- atributos_guia(); A <- setNames(a$A, a$sp)
  ml <- ape::ace(A, tr, method = "ML")
  un <- ape::ace(A, ape::compute.brlen(tr, 1), method = "ML")
  list(arvore = tr, A = A, ml = ml$ace, ic = ml$CI95, un = un$ace,
       no13 = mean(A[c("sp6", "sp7")]))
}

#' Enguias (Collar et al. 2014, dados da prática): modo de alimentação
enguias <- function() {
  X <- read.csv(here::here("aula-03-sinal/elopomorph.csv"), row.names = 1)
  tr <- ape::read.tree(here::here("aula-03-sinal/elopomorph.tre"))
  x <- setNames(factor(X$feed_mode, levels = c("bite", "suction"),
                       labels = c("mordida", "sucção")), rownames(X))[tr$tip.label]
  list(arvore = tr, x = x)
}
cores_eng <- c(mordida = cores[["ferrugem"]], sucção = cores[["accent"]])

ajustes_enguias <- function(e = enguias(), nsim = 100, seed = 701) {
  er  <- phytools::fitMk(e$arvore, e$x, model = "ER")
  ard <- phytools::fitMk(e$arvore, e$x, model = "ARD")
  mg <- phytools::ancr(er); jt <- phytools::ancr(er, type = "joint")
  st_m <- colnames(mg$ace)[apply(mg$ace, 1, which.max)]
  set.seed(seed)
  sm <- phytools::make.simmap(e$arvore, e$x, model = "ER", nsim = nsim, message = FALSE)
  list(e = e, er = er, ard = ard, marg = mg$ace, joint = as.character(jt$ace),
       discorda = names(jt$ace)[st_m != as.character(jt$ace)],
       incertos = sum(apply(mg$ace, 1, max) < .8), mapas = sm,
       contagens = summary(sm)$count)
}

#' Squamata (Aula 6): perdas e reganhos dos membros por mapeamento estocástico
mapas_squamata <- function(nsim = 100, seed = 702) {
  source(here::here("R/aula06_discretos.R"), local = TRUE)
  s <- squamata(); set.seed(seed)
  sm <- phytools::make.simmap(s$arvore, s$bin, model = "ARD", nsim = nsim, message = FALSE)
  summary(sm)$count
}

#' Tendência com fósseis: valores verdadeiros nos nós × estimativas sob BM e
#' sob BM com tendência (phytools::anc.trend)
tendencia_fosseis <- function(n = 60, mu = 2, seed = 703) {
  set.seed(seed)
  tr <- phytools::pbtree(b = 1, d = .4, n = n, extant.only = FALSE, scale = 1)
  tr <- ape::compute.brlen(tr, tr$edge.length * 10)
  x <- phytools::fastBM(tr, a = 0, mu = mu, sig2 = 1, internal = TRUE)
  pt <- x[seq_len(ape::Ntip(tr))]; nos <- x[-seq_len(ape::Ntip(tr))]
  bm <- phytools::fastAnc(tr, pt)
  td <- phytools::anc.trend(tr, pt)
  tibble(verdadeiro = rep(nos, 2), estimado = c(bm, td$ace),
         modelo = rep(c("BM (sem tendência)", "BM com tendência"), each = length(nos)))
}

# ---- Bloco 1: por que estimar ---------------------------------------------------

#' Árvore-guia com A nas pontas; etapa 2 = "?" nos nós; 3 = estimativas
fig_guia_nos <- function(g = guia_ancestrais(), etapa = 3) {
  tr <- g$arvore
  tr$tip.label <- sprintf("%s  %s", tr$tip.label, num(g$A[tr$tip.label], 2))
  nos <- as.character(8:13)
  rot <- if (etapa == 2) setNames(rep("?", 6), nos) else if (etapa == 3) setNames(num(g$ml, 2), nos) else NULL
  plot_arvore(tr, rot_nos = rot, base = 6.5, folga = .35, larg = 1.5)
}

# ---- Bloco 2: contínuos -----------------------------------------------------------

fig_equivalencias <- function(etapa = 2) {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    0, 3, 3.6, 1.1, "parcimônia de mudança\nquadrada não ponderada", cores[["accent"]], cores[["accent"]], "white",
    0, 0, 3.6, 1.1, "parcimônia de mudança\nquadrada ponderada", cores[["texto"]], cores[["texto"]], "white",
    5.5, 1.5, 3.4, 1.1, "máxima verossimilhança\n(BM)", cores[["ferrugem"]], cores[["ferrugem"]], "white")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend, 0, 2.4, 0, .6, 1.85, 2.6, 3.75, 1.95)
  if (etapa >= 2) st <- bind_rows(st, tibble(x = 1.85, y = .4, xend = 3.75, yend = 1.05))
  g <- diagrama(cx, st, base = 5.8) +
    annotate("text", x = -.2, y = 1.5, hjust = 1, size = 5.2, label = "comprimentos\nde ramo iguais") +
    annotate("text", x = 3.3, y = 2.65, hjust = 0, size = 5.2, label = "comprimentos\nde ramo iguais")
  if (etapa >= 2) g <- g + annotate("text", x = 3.3, y = .3, hjust = 0, size = 5.2,
                                    label = "mesmas estimativas", colour = cores[["ferrugem"]])
  g
}

#' Nó 13: só com os descendentes (etapa 1) × com a árvore toda (etapa 2)
fig_no13 <- function(g = guia_ancestrais(), etapa = 1) {
  tr <- g$arvore
  tr$tip.label <- sprintf("%s  %s", tr$tip.label, num(g$A[tr$tip.label], 2))
  sub <- tr$edge[, 2] %in% c(6, 7)
  cor <- if (etapa == 1) ifelse(sub, cores[["ferrugem"]], "#cfc8b8") else rep(cores[["ferrugem"]], nrow(tr$edge))
  rot <- if (etapa == 1) c("13" = num(g$no13, 2)) else c("13" = num(g$ml[["13"]], 2))
  plot_arvore(tr, cor_ramo = cor, rot_nos = rot, base = 6.5, folga = .35, larg = 1.5)
}

#' Estimativas e IC 95% contra a distância da raiz (árvore-guia, fenograma)
fig_ic_guia <- function(g = guia_ancestrais()) {
  tr <- g$arvore; n <- ape::Ntip(tr)
  prof <- ape::node.depth.edgelength(tr)
  val <- c(g$A[tr$tip.label], g$ml)
  ed <- tibble(x0 = prof[tr$edge[, 1]], x1 = prof[tr$edge[, 2]],
               y0 = val[tr$edge[, 1]], y1 = val[tr$edge[, 2]])
  nd <- tibble(no = names(g$ml), x = prof[(n + 1):(n + tr$Nnode)], y = g$ml,
               lo = g$ic[, 1], hi = g$ic[, 2])
  pt <- tibble(x = prof[1:n], y = g$A[tr$tip.label], l = tr$tip.label)
  ggplot() +
    geom_segment(data = ed, aes(x0, y0, xend = x1, yend = y1), colour = cores[["neutro"]], linewidth = .9) +
    geom_errorbar(data = nd, aes(x, ymin = lo, ymax = hi), width = .04, colour = cores[["ferrugem"]], linewidth = 1.1) +
    geom_point(data = nd, aes(x, y), size = 3.5, colour = cores[["ferrugem"]]) +
    geom_text(data = dplyr::filter(nd, no != "11"), aes(x, hi, label = no), vjust = -.6, size = 5, colour = cores[["ferrugem"]]) +
    geom_text(data = dplyr::filter(nd, no == "11"), aes(x - .03, y, label = no), hjust = 1, size = 5, colour = cores[["ferrugem"]]) +
    geom_point(data = pt, aes(x, y), size = 3, colour = cores[["texto"]]) +
    geom_text(data = pt, aes(x + .04, y, label = l), hjust = 0, size = 5) +
    scale_x_continuous(expand = expansion(mult = c(.03, .1))) +
    scale_y_continuous(expand = expansion(mult = c(.05, .12))) +
    labs(x = "distância da raiz", y = "atributo A") + tema_aula(16)
}

#' Anolis: variância das estimativas × distância do nó à raiz
fig_var_prof <- function() {
  an <- anolis(); tr <- an$arvore; x <- an$svl; n <- ape::Ntip(tr)   # anole.data já está em log natural
  fa <- phytools::fastAnc(tr, x, vars = TRUE)
  d <- tibble(prof = ape::node.depth.edgelength(tr)[(n + 1):(n + tr$Nnode)], v = fa$var)
  ggplot(d, aes(prof, v)) +
    geom_point(colour = cores[["accent"]], size = 2.8, alpha = .75) +
    labs(x = "distância do nó à raiz (altura da árvore = 1)", y = "variância da estimativa") +
    tema_aula(16)
}

fig_contmap <- function() {
  an <- anolis(); x <- an$svl   # anole.data já está em log natural
  obj <- phytools::contMap(an$arvore, x, plot = FALSE)
  obj <- phytools::setMap(obj, c("#1c6e8c", "#e9dfc9", "#a3472e"))
  fundo(); par(mar = c(3, 0, 0, 0), xpd = NA)   # a barra fica abaixo da árvore
  plot(obj, fsize = c(.45, .9), lwd = c(3, 5), outline = FALSE, legend = FALSE, mar = c(3, 0, 0, 0))
  phytools::add.color.bar(.4, obj$cols, title = "log SVL", lims = obj$lims, digits = 2, prompt = FALSE,
                          x = 0, y = -4, subtitle = "", fsize = 1)
}

fig_fenograma95 <- function() {
  an <- anolis(); x <- an$svl   # anole.data já está em log natural
  fundo(); phytools::fancyTree(an$arvore, type = "phenogram95", x = x, spread.labels = FALSE,
                               ftype = "off", ylab = "log SVL", xlab = "tempo (altura = 1)")
}

fig_tendencia_fosseis <- function(d = tendencia_fosseis()) {
  ggplot(d, aes(verdadeiro, estimado)) +
    geom_abline(linetype = "22", colour = cores[["neutro"]]) +
    geom_point(colour = cores[["accent"]], size = 2.6, alpha = .8) +
    facet_wrap(~modelo) +
    labs(x = "valor verdadeiro no nó", y = "estimativa") + tema_aula(16)
}

# ---- Bloco 3: discretos -----------------------------------------------------------

fig_enguias_pontas <- function(e = enguias()) {
  fundo(); phytools::plotTree(e$arvore, type = "fan", ftype = "off", lwd = 2, color = cores[["texto"]])
  ape::tiplabels(pie = phytools::to.matrix(e$x, levels(e$x)), piecol = cores_eng, cex = .45)
  legend("bottomleft", legend = names(cores_eng), pch = 21, pt.bg = cores_eng, pt.cex = 2.2, bty = "n", cex = 1.2)
}

#' Pizzas das probabilidades marginais; anéis nos nós em que a conjunta discorda
fig_enguias_marginal <- function(aj, aneis = FALSE) {
  e <- aj$e; tr <- e$arvore
  fundo(); phytools::plotTree(tr, type = "fan", ftype = "off", lwd = 2, color = cores[["texto"]])
  par(fg = "transparent")
  ape::nodelabels(node = as.integer(rownames(aj$marg)), pie = aj$marg, piecol = cores_eng, cex = .7)
  ape::tiplabels(pie = phytools::to.matrix(e$x, levels(e$x)), piecol = cores_eng, cex = .4)
  par(fg = cores[["texto"]])
  if (aneis) {
    pp <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
    k <- as.integer(aj$discorda)
    points(pp$xx[k], pp$yy[k], cex = 5, lwd = 3.5, col = cores[["texto"]])
  }
  legend("bottomleft", legend = names(cores_eng), pch = 21, pt.bg = cores_eng, pt.cex = 2.2, bty = "n", cex = 1.2)
}

fig_um_mapa <- function(aj, i = 1) {
  fundo(); phytools::plotSimmap(aj$mapas[[i]], colors = cores_eng, type = "fan", ftype = "off", lwd = 3.5)
}

#' Prancha com os primeiros mapas estocásticos e o número de mudanças de cada um
fig_prancha_mapas <- function(aj, n = 8) {
  fundo(); par(mfrow = c(2, n / 2))
  for (i in seq_len(n)) {
    m <- aj$mapas[[i]]
    phytools::plotSimmap(m, colors = cores_eng, type = "fan", ftype = "off", lwd = 1.6,
                         mar = c(0, 0, 1.6, 0))
    title(sprintf("mapa %d: %d mudanças", i, phytools::countSimmap(m)$N),
          cex.main = 1.3, font.main = 1, col.main = cores[["texto"]])
  }
}

fig_densitymap <- function(aj) {
  dm <- phytools::densityMap(aj$mapas, plot = FALSE, states = levels(aj$e$x))
  dm <- phytools::setMap(dm, c(cores[["ferrugem"]], "#e9dfc9", cores[["accent"]]))
  fundo(); par(xpd = NA)   # a barra fica abaixo da árvore
  plot(dm, type = "fan", fsize = c(.55, .9), lwd = c(4, 5), outline = FALSE, ftype = "off", legend = FALSE,
       mar = c(4, 0, 0, 0))
  u <- par()$usr
  phytools::add.color.bar(.35 * diff(u[1:2]), dm$cols, title = "P(sucção)", lims = c(0, 1), digits = 1, prompt = FALSE,
                          x = u[1] + .325 * diff(u[1:2]), y = u[3] - .04 * diff(u[3:4]), subtitle = "", fsize = 1)
}

fig_contagens <- function(cont, rotulos) {
  d <- as_tibble(cont[, -1, drop = FALSE]) |>
    tidyr::pivot_longer(everything(), names_to = "tipo", values_to = "n") |>
    mutate(tipo = factor(rotulos[tipo], levels = rotulos))
  ggplot(d, aes(n, fill = tipo)) +
    geom_histogram(binwidth = 1, colour = "white", alpha = .9) +
    facet_wrap(~tipo) +
    scale_fill_manual(values = c(cores[["ferrugem"]], cores[["accent"]]), name = NULL) +
    labs(x = "número de mudanças por mapa", y = "mapas") + tema_aula(16) +
    theme(legend.position = "none")
}

# ---- Bloco 4: visualização ---------------------------------------------------------

#' Filomorfoespaço: a árvore projetada no espaço de dois atributos, com os
#' nós nas estimativas por ML
fig_filomorfo <- function() {
  an <- anolis(); data(anole.data, package = "phytools", envir = environment())
  X <- as.matrix(anole.data[an$arvore$tip.label, c("SVL", "TL")])   # já em log
  colnames(X) <- c("log SVL", "log cauda")
  fundo(); phytools::phylomorphospace(an$arvore, X, label = "off", node.size = c(.6, 1.1),
                                      control = list(col.node = setNames(
                                        c(rep(cores[["accent"]], ape::Ntip(an$arvore)), rep(cores[["ferrugem"]], an$arvore$Nnode)),
                                        1:(ape::Ntip(an$arvore) + an$arvore$Nnode))),
                                      xlab = "log SVL", ylab = "log cauda")
}

# ---- Padrões no filomorfoespaço (ideia da Fig. 2 de Adams & Collyer 2019) ----------

#' Dois clados de 16 espécies e quatro cenários simulados em dois atributos:
#' sobreposição, divergência, disparidade diferente e tendência no clado 2
padroes_filomorfo <- function(semente = 707) {
  set.seed(semente)
  c1 <- phytools::pbtree(n = 16, scale = .7, tip.label = paste0("a", 1:16))
  c2 <- phytools::pbtree(n = 16, scale = .7, tip.label = paste0("b", 1:16))
  nw <- function(t) sub(";$", "", ape::write.tree(t))
  tr <- ape::read.tree(text = sprintf("(%s:0.3,%s:0.3);", nw(c1), nw(c2)))
  # clado: valor inicial (raiz do clado), taxa e tendência por atributo
  sim <- function(t, a, s2 = 1, mu = c(0, 0))
    sapply(1:2, \(k) a[k] + phytools::fastBM(t, sig2 = s2, mu = mu[k]))
  X1 <- sim(c1, c(0, 0))                    # o clado 1 é o mesmo nos quatro painéis
  c2_esp <- ape::compute.brlen(c2, 1)       # tendência por evento de especiação
  X2 <- list(
    "Sobreposição" = sim(c2, c(0, 0)),
    "Divergência"  = sim(c2, c(3, 0)),
    "Disparidade diferente" = sim(c2, c(2, 0), .02),
    "Tendência"    = sim(c2_esp, c(.3, -.3), .004, c(.4, .25)))
  purrr::imap(X2, \(x2, nome) {
    X <- rbind(X1, x2)
    X <- X[tr$tip.label, ]
    anc <- apply(X, 2, \(x) phytools::fastAnc(tr, x))
    list(arvore = tr, X = X, anc = anc, nome = nome)
  })
}

fig_padrao_filomorfo <- function(p) {
  tr <- p$arvore; n <- ape::Ntip(tr)
  xy <- rbind(p$X, p$anc)
  ed <- tibble(x0 = xy[tr$edge[, 1], 1], y0 = xy[tr$edge[, 1], 2],
               x1 = xy[tr$edge[, 2], 1], y1 = xy[tr$edge[, 2], 2])
  pt <- tibble(x = p$X[, 1], y = p$X[, 2],
               clado = ifelse(startsWith(rownames(p$X), "a"), "clado 1", "clado 2"))
  ggplot() +
    geom_segment(data = ed, aes(x0, y0, xend = x1, yend = y1), colour = "#b9b2a5", linewidth = .7) +
    geom_point(aes(p$anc[, 1], p$anc[, 2]), colour = cores[["neutro"]], size = 1.6) +
    geom_point(data = pt, aes(x, y, colour = clado), size = 2.6) +
    scale_colour_manual(values = c("clado 1" = cores[["accent"]], "clado 2" = cores[["verde"]]), name = NULL) +
    labs(title = p$nome, x = "PC 1", y = "PC 2") +
    tema_aula(15) + theme(legend.position = "none", plot.title = element_text(size = 17))
}

#' A árvore dos dois clados, com as cores dos painéis
fig_arvore_padroes <- function(p) {
  tr <- p$arvore
  cl <- ifelse(startsWith(tr$tip.label, "a"), cores[["accent"]], cores[["verde"]])
  cor_ramo <- rep(cores[["neutro"]], nrow(tr$edge))
  fundo(); par(mar = c(1, 1, 1, 1))
  ape::plot.phylo(tr, show.tip.label = FALSE, edge.color = cor_ramo, edge.width = 1.6)
  ape::tiplabels(pch = 19, col = cl, cex = 1.1)
  legend("topleft", legend = c("clado 1", "clado 2", "estados ancestrais"), pch = 19,
         col = c(cores[["accent"]], cores[["verde"]], cores[["neutro"]]), bty = "n", cex = 1.3)
}

fig_padroes_filomorfo <- function(pp = padroes_filomorfo()) {
  patchwork::wrap_plots(lapply(pp, fig_padrao_filomorfo), ncol = 2) &
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}
