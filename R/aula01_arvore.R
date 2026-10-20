# Aula 1: a árvore e os dados (não independência, ler e escrever árvores,
# de onde vem a árvore, forma da árvore).
library(patchwork)
source(here::here("R/aula03_base.R"))   # tema, árvore, matriz, diagrama

fundo <- function() par(bg = cores[["papel"]], fg = cores[["texto"]], col.axis = cores[["texto"]],
                        col.lab = cores[["texto"]], family = "sans")

# ---- Bloco 0: o problema do parentesco comum ---------------------------------

#' Dois clados de 10 espécies que se separaram cedo; X e Y evoluem por BM,
#' independentes entre si
dois_clados <- function(seed = 101) {
  set.seed(seed)
  a <- phytools::pbtree(n = 10, scale = .3, tip.label = paste0("a", 1:10))
  b <- phytools::pbtree(n = 10, scale = .3, tip.label = paste0("b", 1:10))
  op <- options(OutDec = "."); on.exit(options(op))
  nw <- function(t) sub(";$", "", ape::write.tree(t))
  tr <- ape::read.tree(text = sprintf("(%s:0.7,%s:0.7);", nw(a), nw(b)))
  options(op)
  repeat {
    x <- phytools::fastBM(tr); y <- phytools::fastBM(tr)
    d <- tibble(sp = names(x), x = x, y = y[names(x)], clado = ifelse(grepl("^a", names(x)), "A", "B"))
    p_tot <- summary(lm(y ~ x, d))$coef[2, 4]
    p_den <- sapply(c("A", "B"), \(k) summary(lm(y ~ x, d[d$clado == k, ]))$coef[2, 4])
    if (p_tot < .01 && all(p_den > .2)) break
  }
  list(arvore = tr, dados = d, p = p_tot, r = cor(d$x, d$y))
}
cores_clado <- c(A = cores[["accent"]], B = cores[["ferrugem"]])

fig_parentesco <- function(dc = dois_clados(), etapa = 1) {
  d <- dc$dados
  if (etapa == 1)
    return(ggplot(d, aes(x, y)) + geom_smooth(method = "lm", se = FALSE, colour = cores[["texto"]], linewidth = 1.2) +
             geom_point(size = 4, colour = cores[["texto"]]) + labs(x = "X", y = "Y") + tema_aula(17))
  g1 <- ggplot(d, aes(x, y, colour = clado)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
    geom_point(size = 4) + scale_colour_manual(values = cores_clado, guide = "none") +
    labs(x = "X", y = "Y") + tema_aula(17)
  cor <- ifelse(dc$arvore$edge[, 2] %in% c(which(grepl("^a", dc$arvore$tip.label)),
                                         phytools::getDescendants(dc$arvore, ape::getMRCA(dc$arvore, paste0("a", 1:10)))),
                cores_clado[["A"]], cores_clado[["B"]])
  cor[dc$arvore$edge[, 1] == ape::Ntip(dc$arvore) + 1] <- cores[["texto"]]
  g2 <- plot_arvore(dc$arvore, cor_ramo = cor, rotulos = FALSE, larg = 1.2, folga = .05)
  g1 + g2 + plot_layout(widths = c(1.4, 1))
}

#' A árvore "real" (dois clados) × a filogenia em estrela
fig_estrela <- function(dc = dois_clados()) {
  tr <- dc$arvore
  est <- ape::stree(ape::Ntip(tr), "star", tip.label = tr$tip.label); est$edge.length <- rep(1, ape::Ntip(tr))
  t1 <- plot_arvore(tr, rotulos = FALSE, larg = 1.2, folga = .05) + labs(title = "filogenia \"real\"")
  t2 <- plot_arvore(est, rotulos = FALSE, larg = 1.2, folga = .05) + labs(title = "filogenia em estrela")
  tt <- theme(plot.title = element_text(hjust = .5, size = 20, face = "bold", colour = cores[["texto"]]))
  (t1 + tt) + (t2 + tt)
}

#' Linha do tempo das obras que fundaram e consolidaram os MFC: vertical,
#' uma obra por linha (fora de escala), para os rótulos não se sobreporem
fig_linha_tempo <- function() {
  d <- tibble::tribble(
    ~ano, ~autor, ~obra,
    1985, "Cheverud, Dow & Leutenegger", "autocorrelação filogenética",
    1985, "Felsenstein", "contrastes independentes",
    1989, "Grafen", "regressão filogenética",
    1991, "Harvey & Pagel", "The Comparative Method in Evolutionary Biology",
    2011, "Nunn", "The Comparative Approach in Evolutionary Anthropology and Biology",
    2014, "Garamszegi (org.)", "Modern Phylogenetic Comparative Methods",
    2022, "Revell & Harmon", "Phylogenetic Comparative Methods in R")
  d$y <- rev(seq_len(nrow(d)))
  ggplot(d) +
    annotate("segment", x = 0, xend = 0, y = .5, yend = nrow(d) + .5, linewidth = 1.2, colour = cores[["texto"]]) +
    geom_point(aes(0, y), size = 4, colour = cores[["ferrugem"]]) +
    geom_text(aes(-.15, y, label = ano), hjust = 1, size = 6, fontface = "bold", colour = cores[["ferrugem"]]) +
    geom_text(aes(.15, y + .1, label = autor), hjust = 0, size = 6, fontface = "bold", colour = cores[["texto"]]) +
    geom_text(aes(.15, y - .3, label = obra), hjust = 0, size = 4.8, colour = "#4a5568") +
    scale_x_continuous(limits = c(-1.2, 8)) + scale_y_continuous(limits = c(0.3, nrow(d) + .6)) +
    theme_void() + theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}

# ---- Bloco 1: o que é uma filogenia ----------------------------------------------

fig_partes <- function() {
  tr <- arvore_guia(); cd <- coord_arvore(tr); n <- cd$nos
  xmax <- max(n$x); no10 <- n[n$no == 10, ]; raiz <- n[n$no == 8, ]
  r6 <- dplyr::filter(cd$ramos, filho == 10)
  no10 <- n[n$no == 11, ]
  cl <- dplyr::filter(n, no %in% 5:7)
  cor <- ifelse(tr$edge[, 2] %in% c(5:7, 12, 13), cores[["accent"]], cores[["texto"]])
  plot_arvore(tr, cor_ramo = cor, base = 6, folga = .45, larg = 1.4) +
    annotate("point", x = raiz$x, y = raiz$y, size = 5, colour = cores[["ferrugem"]]) +
    annotate("label", x = raiz$x + .05, y = raiz$y + 1.1, label = "raiz", size = 6, colour = cores[["ferrugem"]],
             fill = cores[["papel"]], border.colour = NA) +
    annotate("point", x = no10$x, y = no10$y, size = 5, colour = cores[["ferrugem"]]) +
    annotate("label", x = no10$x - .05, y = no10$y - .8, hjust = 1, label = "nó (ancestral comum de sp3 e sp4)", size = 5.5,
             colour = cores[["ferrugem"]], fill = cores[["papel"]], border.colour = NA) +
    annotate("label", x = (r6$x0 + r6$x1) / 2, y = r6$y + .45, label = "ramo", size = 6, fill = cores[["papel"]],
             border.colour = NA) +
    annotate("segment", x = r6$x0, xend = r6$x1, y = r6$y - .4, yend = r6$y - .4, colour = cores[["neutro"]],
             arrow = arrow(length = unit(.18, "cm"), ends = "both")) +
    annotate("text", x = (r6$x0 + r6$x1) / 2, y = r6$y - .75, label = "comprimento do ramo", size = 5,
             colour = cores[["neutro"]]) +
    annotate("text", x = xmax * 1.16, y = 4, hjust = 0, label = "← pontas (espécies)", size = 6) +
    annotate("segment", x = xmax * 1.42, xend = xmax * 1.42, y = min(cl$y) - .2, yend = max(cl$y) + .2,
             linewidth = 1.4, colour = cores[["accent"]]) +
    annotate("text", x = xmax * 1.45, y = mean(cl$y), hjust = 0, label = "clado", size = 6, colour = cores[["accent"]])
}

#' A mesma árvore não enraizada e enraizada em dois lugares
fig_enraizar <- function() {
  u <- ape::unroot(ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);"))
  r1 <- ape::root(u, "A", resolve.root = TRUE); r1 <- ape::compute.brlen(r1, 1)
  r2 <- ape::read.tree(text = "((A:1,B:1):1,(C:1,D:1):1);")
  fundo(); par(mfrow = c(1, 3), mar = c(1, 2, 3, 2), cex.main = 1.8, xpd = NA)
  ape::plot.phylo(u, type = "unrooted", lab4ut = "horizontal", cex = 2.4, label.offset = .15, edge.width = 3, edge.color = cores[["texto"]],
                  main = "não enraizada", no.margin = FALSE)
  ape::plot.phylo(r1, cex = 2.2, edge.width = 3, edge.color = cores[["texto"]], main = "raiz no ramo de A", label.offset = .1)
  ape::plot.phylo(r2, cex = 2.2, edge.width = 3, edge.color = cores[["texto"]], main = "raiz no ramo central", label.offset = .1)
}

fig_politomia <- function() {
  d <- ape::compute.brlen(ape::read.tree(text = "(A,(B,(C,D)));"), method = "Grafen")
  p <- ape::compute.brlen(ape::read.tree(text = "(A,(B,C,D));"), method = "Grafen")
  tt <- theme(plot.title = element_text(hjust = .5, size = 20, face = "bold", colour = cores[["texto"]]))
  (plot_arvore(d, base = 7, larg = 1.6) + labs(title = "dicotomias") + tt) +
    (plot_arvore(p, base = 7, larg = 1.6) + labs(title = "politomia") + tt)
}

#' Cladograma, filograma e cronograma da árvore-guia
fig_tipos <- function(seed = 102) {
  tr <- arvore_guia()
  cl <- ape::compute.brlen(tr, method = "Grafen")
  set.seed(seed); fi <- tr; fi$edge.length <- tr$edge.length * exp(rnorm(nrow(tr$edge), 0, .6))
  tt <- theme(plot.title = element_text(hjust = .5, size = 19, face = "bold", colour = cores[["texto"]]))
  (plot_arvore(cl, base = 5.5, larg = 1.3) + labs(title = "cladograma") + tt) +
    (plot_arvore(fi, base = 5.5, larg = 1.3) + labs(title = "filograma") + tt) +
    (plot_arvore(tr, base = 5.5, larg = 1.3) + labs(title = "cronograma") + tt)
}

#' Bootstrap de verdade: NJ (distância K80) no woodmouse do ape, 100 réplicas,
#' árvores enraizadas na primeira sequência
suporte_woodmouse <- function(B = 100, seed = 103) {
  data(woodmouse, package = "ape", envir = environment())
  f <- function(x) ape::root(ape::nj(ape::dist.dna(x, model = "K80")), 1, resolve.root = TRUE)
  tr <- f(woodmouse)
  set.seed(seed)
  bt <- ape::boot.phylo(tr, woodmouse, f, B = B, trees = TRUE, quiet = TRUE, rooted = TRUE)
  tr$edge.length[tr$edge.length < 0] <- 0
  list(arvore = tr, suporte = bt$BP)
}

fig_suporte <- function(s = suporte_woodmouse()) {
  tr <- s$arvore; n <- ape::Ntip(tr)
  rot <- setNames(as.character(s$suporte), n + seq_len(tr$Nnode))
  rot <- rot[!is.na(s$suporte) & s$suporte >= 50]
  rot <- rot[names(rot) != as.character(n + 1)]
  plot_arvore(tr, rot_nos = rot, base = 4.6, folga = .2, larg = 1.1)
}

# ---- Bloco 2: como escrever uma árvore -------------------------------------------

newick_passos <- c("(A,B);", "(A:2.75,B:3.76);", "((A:2.75,B:3.76):2.0,C:2.6);",
                   "((A:2.75,B:3.76)N1:2.0,C:2.6);", "((A:2.75,B:3.76)N1:2.0,C:2.6)N2;")

fig_newick <- function(etapa = 1) {
  tr <- ape::read.tree(text = newick_passos[etapa])
  if (etapa == 1) tr <- ape::compute.brlen(tr, 1)
  rot <- if (etapa >= 4) c(setNames(tr$node.label, ape::Ntip(tr) + seq_len(tr$Nnode))) else NULL
  if (!is.null(rot)) rot <- rot[rot != ""]
  plot_arvore(tr, comprimentos = etapa >= 2, rot_nos = rot, base = 8, larg = 1.8, folga = .12, dig = 2) +
    labs(title = newick_passos[etapa]) +
    theme(plot.title = element_text(family = "mono", size = 26, face = "bold", colour = cores[["accent_dark"]],
                                    hjust = 0, margin = margin(b = 10)))
}

# ---- Bloco 3: de onde vem a árvore --------------------------------------------------

fig_fluxo <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    -.2, 2, 3.0, .9, "sequências homólogas\n(GenBank ou bancada)", cores[["texto"]], cores[["texto"]], "white",
    3.3, 2, 2.6, .9, "alinhamento\n(MAFFT, MUSCLE)", cores[["texto"]], cores[["texto"]], "white",
    7, 3.4, 3.2, .9, "distâncias corrigidas por\num modelo (ex.: K80)", cores[["papel"]], cores[["accent"]], cores[["accent"]],
    10.8, 3.4, 2.6, .9, "UPGMA, NJ", cores[["accent"]], cores[["accent"]], "white",
    7, 2, 3.2, .9, "escolha do modelo\n(por partição)", cores[["papel"]], cores[["ferrugem"]], cores[["ferrugem"]],
    10.8, 2, 2.6, .9, "ML, inferência\nbayesiana", cores[["ferrugem"]], cores[["ferrugem"]], "white",
    7, .6, 3.2, .9, "menor número\nde mudanças", cores[["papel"]], cores[["verde"]], cores[["verde"]],
    10.8, .6, 2.6, .9, "máxima\nparcimônia", cores[["verde"]], cores[["verde"]], "white")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend,
                        1.3, 2, 2, 2, 4.6, 2.2, 5.4, 3.3, 4.6, 2, 5.4, 2, 4.6, 1.8, 5.4, .7,
                        8.6, 3.4, 9.5, 3.4, 8.6, 2, 9.5, 2, 8.6, .6, 9.5, .6)
  diagrama(cx, st, base = 4.6)
}

# ---- Bloco 6: desenho ----------------------------------------------------------------

fig_ignorancias <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    0, 2.4, 3.4, .9, "topologia", cores[["papel"]], cores[["accent"]], cores[["accent"]],
    0, 1.2, 3.4, .9, "comprimentos de ramo\ne tempos de divergência", cores[["papel"]], cores[["accent"]], cores[["accent"]],
    0, 0, 3.4, .9, "modelo de evolução\ndos atributos", cores[["papel"]], cores[["accent"]], cores[["accent"]],
    5.5, 1.2, 3.2, 1.2, "inferência\ncomparativa", cores[["ferrugem"]], cores[["ferrugem"]], "white")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend, 1.75, 2.4, 3.85, 1.5, 1.75, 1.2, 3.85, 1.2, 1.75, 0, 3.85, .9)
  diagrama(cx, st, base = 6)
}

fig_phylotargeting <- function() {
  passos <- c("problema: espécies sem dados da variável resposta",
              "requisitos: variáveis resposta e preditoras, hipótese, filogenia",
              "PhyloTargeting: um escore para cada par de espécies",
              "coletar dados das espécies com escores altos",
              "testar a hipótese (PGLS, ANOVA, PGLMM...)")
  cx <- tibble(x = 0, y = rev(seq_along(passos)) * 1.05, w = 8.6, h = .78, rotulo = passos,
               fill = c(cores[["papel"]], cores[["papel"]], cores[["ferrugem"]], cores[["papel"]], cores[["papel"]]),
               cor = c(cores[["texto"]], cores[["texto"]], cores[["ferrugem"]], cores[["texto"]], cores[["texto"]]),
               txt = c(cores[["texto"]], cores[["texto"]], "white", cores[["texto"]], cores[["texto"]]))
  st <- tibble(x = 0, y = cx$y[-5] - .39, xend = 0, yend = cx$y[-1] + .39)
  diagrama(cx, st, base = 5.8)
}

# ---- Bloco 7: a forma da árvore -------------------------------------------------------

#' Índice de Colless: soma, nos nós internos, de |pontas à esquerda − à direita|
colless <- function(tr) {
  n <- ape::Ntip(tr)
  nd <- sapply((n + 1):(n + tr$Nnode), \(k) {
    f <- tr$edge[tr$edge[, 1] == k, 2]
    if (length(f) != 2) return(NA)
    s <- sapply(f, \(x) if (x <= n) 1 else length(phytools::getDescendants(tr, x)[phytools::getDescendants(tr, x) <= n]))
    abs(s[1] - s[2])
  })
  sum(nd, na.rm = TRUE)
}

formas <- function() {
  bal <- ape::compute.brlen(ape::stree(8, "balanced"), method = "Grafen")
  pen <- ape::compute.brlen(ape::stree(8, "left"), method = "Grafen")
  bal$tip.label <- pen$tip.label <- paste0("sp", 1:8)
  list(bal = bal, pen = pen,
       ic = c(bal = colless(bal), pen = colless(pen)),
       gama = c(bal = ape::gammaStat(bal), pen = ape::gammaStat(pen)))
}

fig_formas <- function(f = formas()) {
  tt <- theme(plot.title = element_text(hjust = .5, size = 19, face = "bold", colour = cores[["texto"]]))
  (plot_arvore(f$bal, base = 5.5) + labs(title = paste0("balanceada: I = ", f$ic[["bal"]])) + tt) +
    (plot_arvore(f$pen, base = 5.5) + labs(title = paste0("pente: I = ", f$ic[["pen"]])) + tt)
}

fig_formas_C <- function(f = formas()) {
  tt <- theme(plot.title = element_text(hjust = .5, size = 18, face = "bold", colour = cores[["texto"]]))
  m <- function(tr) { C <- ape::vcv(tr); C / max(C) }
  pinta <- function(C) matrix(grDevices::colorRampPalette(c("#f3eee3", cores[["accent"]]))(101)[round(C * 100) + 1], nrow(C))
  (plot_matriz(m(f$bal), fill = pinta(m(f$bal)), digitos = 1, base = 3.6) + labs(title = "C da balanceada") + tt) +
    (plot_matriz(m(f$pen), fill = pinta(m(f$pen)), digitos = 1, base = 3.6) + labs(title = "C do pente") + tt)
}
