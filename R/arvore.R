library(ape)
library(dplyr)
library(tibble)
library(ggplot2)
source(here::here("R/tema.R"))

#' A árvore-guia do curso: as 7 espécies de Freckleton et al. (2002)
#'
#' Comprimentos de ramo lidos da figura usada no original (Aula 2, slides
#' 36–44). Ultramétrica, altura 1,44.
arvore_guia <- function() {
  ape::read.tree(text = paste0(
    "(((sp1:0.18,sp2:0.18):1.06,(sp3:0.02,sp4:0.02):1.22):0.20,",
    "(sp5:1.37,(sp6:0.23,sp7:0.23):1.14):0.07);"))
}

#' Atributos A e B das 7 espécies, do mesmo exemplo de Freckleton et al. (2002)
atributos_guia <- function() {
  tibble::tibble(
    sp = paste0("sp", 1:7),
    A  = c(2.51, 2.80, 2.23, 2.32, 1.64, 1.18, 0.58),
    B  = c(0.38, 0.89, 0.88, 1.60, 0.33, 1.42, 0.31))
}

#' Coordenadas de uma árvore para desenhar em ggplot2
#'
#' x = distância da raiz; y = posição vertical das pontas (n a 1, de cima
#' para baixo: a primeira ponta fica no alto, como na figura do original) e
#' dos nós (média dos filhos).
#' @return lista com `ramos` (um por aresta: x0, x1, y, pai, filho, aresta),
#'   `verticais` (um por nó interno) e `nos` (todos os nós)
coord_arvore <- function(tr) {
  n  <- ape::Ntip(tr)
  x  <- ape::node.depth.edgelength(tr)
  y  <- n + 1 - ape::node.height(tr)
  nos <- tibble::tibble(no = seq_along(x), x = x, y = y, ponta = no <= n,
                        rotulo = c(tr$tip.label, rep(NA, length(x) - n)))
  ramos <- tibble::tibble(aresta = seq_len(nrow(tr$edge)),
                          pai = tr$edge[, 1], filho = tr$edge[, 2],
                          x0 = x[pai], x1 = x[filho], y = y[filho],
                          comp = tr$edge.length)
  verticais <- ramos |>
    group_by(pai) |>
    summarise(x = first(x0), y0 = min(y), y1 = max(y), .groups = "drop")
  list(ramos = ramos, verticais = verticais, nos = nos)
}

#' Arestas no caminho da raiz até um nó (ou ponta)
caminho_raiz <- function(tr, no) {
  if (is.character(no)) no <- match(no, tr$tip.label)
  raiz <- ape::Ntip(tr) + 1
  out <- integer(0)
  while (no != raiz) {
    e <- which(tr$edge[, 2] == no)
    out <- c(e, out)
    no <- tr$edge[e, 1]
  }
  out
}

#' Arestas compartilhadas pelos caminhos raiz → i e raiz → j
caminho_comum <- function(tr, i, j) intersect(caminho_raiz(tr, i), caminho_raiz(tr, j))

#' Desenha uma árvore em ggplot2, no estilo da árvore-guia
#'
#' @param cor_ramo cor de cada aresta (vetor do tamanho de `nrow(tr$edge)`);
#'   NULL = cor do texto
#' @param larg espessura de cada aresta (vetor ou escalar)
#' @param rotulos escrever o nome das pontas
#' @param comprimentos escrever o comprimento sobre cada ramo
#' @param nos numerar os nós internos (como no ape)
#' @param rot_nos rótulos livres para os nós internos (vetor nomeado pelo
#'   número do nó)
#' @param base tamanho do texto
plot_arvore <- function(tr, cor_ramo = NULL, larg = 1.3, rotulos = TRUE,
                        comprimentos = FALSE, nos = FALSE, rot_nos = NULL,
                        base = 5.2, folga = .18, dig = 2) {
  cd <- coord_arvore(tr)
  r  <- cd$ramos
  r$cor  <- if (is.null(cor_ramo)) cores[["texto"]] else cor_ramo
  r$larg <- rep_len(larg, nrow(r))
  # a vertical de cada nó leva a cor do ramo que chega nele (preto na raiz)
  v <- cd$verticais |>
    left_join(select(r, filho, cor_pai = cor, larg_pai = larg), by = c("pai" = "filho")) |>
    mutate(cor_pai = ifelse(is.na(cor_pai), cores[["texto"]], cor_pai),
           larg_pai = ifelse(is.na(larg_pai), larg[1], larg_pai))
  xmax <- max(cd$nos$x)
  g <- ggplot() +
    geom_segment(data = v, aes(x = x, xend = x, y = y0, yend = y1),
                 colour = v$cor_pai, linewidth = v$larg_pai, lineend = "round") +
    geom_segment(data = r, aes(x = x0, xend = x1, y = y, yend = y),
                 colour = r$cor, linewidth = r$larg, lineend = "round")
  if (rotulos) {
    pt <- filter(cd$nos, ponta)
    g <- g + geom_text(data = pt, aes(x + xmax * .03, y, label = rotulo),
                       hjust = 0, size = base, colour = cores[["texto"]])
  }
  if (comprimentos) {
    # ramos curtos demais para o rótulo centrado: rótulo terminando no fim
    # do ramo, acima dele
    # (abaixo do ramo quando a ponta fica abaixo da irmã, para não colidir)
    r$curto <- r$comp < .1 * xmax
    y_pai <- cd$nos$y[r$pai]
    r$dy <- ifelse(r$curto & r$y < y_pai, -.28, .26)
    g <- g + geom_text(data = r, aes(ifelse(curto, x1, (x0 + x1) / 2), y + dy,
                                     label = format(round(comp, dig), nsmall = dig, decimal.mark = ",")),
                       hjust = ifelse(r$curto, 1, .5),
                       size = base * .78, colour = "#4a5568")
  }
  if (nos) {
    ni <- filter(cd$nos, !ponta)
    g <- g + geom_label(data = ni, aes(x, y, label = no), size = base * .7,
                        fill = cores[["papel"]], colour = cores[["accent_dark"]],
                        linewidth = .3)
  }
  if (!is.null(rot_nos)) {
    ni <- filter(cd$nos, no %in% as.integer(names(rot_nos))) |>
      mutate(l = rot_nos[as.character(no)])
    g <- g + geom_label(data = ni, aes(x, y, label = l), size = base * .78,
                        fill = cores[["papel"]], colour = cores[["accent_dark"]],
                        linewidth = .3)
  }
  g +
    scale_x_continuous(expand = expansion(mult = c(.06, folga))) +
    scale_y_continuous(expand = expansion(add = .5)) +
    theme_void() +
    theme(plot.background  = element_rect(fill = cores[["papel"]], colour = NA),
          panel.background = element_rect(fill = cores[["papel"]], colour = NA))
}

#' Simula um BM ao longo da árvore, em passos pequenos, para o fenograma
#'
#' @return tibble com aresta, t (tempo desde a raiz), x (valor) e ponta
#'   (rótulo da ponta para as arestas terminais)
sim_bm_arvore <- function(tr, sigma2 = 1, dt = .005, z0 = 0, seed = 1) {
  set.seed(seed)
  cd <- coord_arvore(tr)
  n  <- ape::Ntip(tr)
  valor <- numeric(max(tr$edge)); valor[n + 1] <- z0
  ordem <- ape::reorder.phylo(tr, "cladewise")$edge
  out <- vector("list", nrow(ordem))
  for (k in seq_len(nrow(ordem))) {
    e  <- which(tr$edge[, 1] == ordem[k, 1] & tr$edge[, 2] == ordem[k, 2])
    pai <- tr$edge[e, 1]; filho <- tr$edge[e, 2]
    t0 <- cd$nos$x[pai]; t1 <- cd$nos$x[filho]
    passos <- max(2, ceiling((t1 - t0) / dt))
    tt <- seq(t0, t1, length.out = passos + 1)
    xx <- valor[pai] + c(0, cumsum(rnorm(passos, 0, sqrt(sigma2 * diff(tt)))))
    valor[filho] <- xx[length(xx)]
    out[[k]] <- tibble::tibble(aresta = e, t = tt, x = xx,
                               ponta = if (filho <= n) tr$tip.label[filho] else NA)
  }
  dplyr::bind_rows(out)
}

#' Fenograma: valor do atributo × tempo, um traço por ramo
fenograma <- function(sim, tr, cor_ramo = NULL, base = 5, rotulos = TRUE) {
  cor_ramo <- if (is.null(cor_ramo)) rep(cores[["texto"]], nrow(tr$edge)) else cor_ramo
  sim$cor <- cor_ramo[sim$aresta]
  fim <- sim |> filter(!is.na(ponta)) |> group_by(ponta) |> slice_tail(n = 1)
  g <- ggplot(sim, aes(t, x, group = aresta)) +
    geom_path(colour = sim$cor, linewidth = .7)
  if (rotulos) g <- g + geom_text(data = fim, aes(t, x, label = ponta), inherit.aes = FALSE,
                                  hjust = -.15, size = base, colour = cores[["texto"]])
  g +
    scale_x_continuous(expand = expansion(mult = c(.02, .1))) +
    labs(x = "tempo desde a raiz", y = "valor do atributo") +
    tema_aula()
}
