# Aula 3, bloco 7: PVR (Diniz-Filho et al. 1998), curva PSR (Diniz-Filho et
# al. 2012) e partição da variação. Ideia dos esquemas: slides de J. A. F.
# Diniz-Filho, redesenhados.
source(here::here("R/aula03_base.R"))

#' T = P + S
fig_tps <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo, ~fill, ~cor, ~txt,
    0, 1, 2.8, 1.1, "Variação total\nT", cores[["accent"]], cores[["accent"]], "white",
    4.2, 2, 3.4, 1.1, "Componente filogenético\nP", cores[["ocre"]], cores[["ocre"]], "white",
    4.2, 0, 3.4, 1.1, "Componente específico\nS", cores[["verde"]], cores[["verde"]], "white")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend, 1.45, 1.3, 2.45, 1.85, 1.45, .7, 2.45, .15)
  diagrama(cx, st, base = 6) + annotate("text", x = 2.1, y = -1, label = "T = P + S", size = 9,
                                        fontface = "bold", colour = cores[["texto"]])
}

#' PVR passo a passo na árvore-guia (atributo A)
#' @param etapa 1 = distâncias; 2 = matriz duplamente centralizada; 3 = os
#'   três primeiros autovetores nas pontas; 4 = P e S
fig_pvr <- function(etapa = 1) {
  tr <- arvore_guia(); a <- atributos_guia(); x <- setNames(a$A, a$sp)
  p <- pvr_passos(tr, x, k = 2)
  arv <- plot_arvore(tr, base = 5.2)
  if (etapa == 1) return(arv + plot_matriz(p$D, base = 4.6, digitos = 2) + plot_layout(widths = c(1, 1.3)))
  if (etapa == 2) {
    G <- p$G; dimnames(G) <- dimnames(p$D)
    return(arv + plot_matriz(G, base = 4.4, digitos = 2) + plot_layout(widths = c(1, 1.3)))
  }
  if (etapa == 3) {
    cd <- coord_arvore(tr); pt <- dplyr::filter(cd$nos, ponta)
    d <- purrr::map_dfr(1:3, \(j) tibble(rotulo = rownames(p$V), v = p$V[, j],
                                         ev = paste0("autovetor ", j, " (", num(100 * p$valores[j] / sum(p$valores), 0), " %)"))) |>
      left_join(dplyr::select(pt, rotulo, y), by = "rotulo")
    barras <- ggplot(d, aes(v, y)) +
      geom_vline(xintercept = 0, colour = cores[["neutro"]]) +
      geom_segment(aes(x = 0, xend = v, yend = y, colour = v > 0), linewidth = 7) +
      facet_wrap(~ev, nrow = 1) +
      scale_colour_manual(values = c(`TRUE` = cores[["accent"]], `FALSE` = cores[["ferrugem"]]), guide = "none") +
      scale_y_continuous(expand = expansion(add = .5)) +
      theme_void() + theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA),
                           strip.text = element_text(size = 15, face = "bold", colour = cores[["texto"]]))
    return(arv + barras + plot_layout(widths = c(1, 2.2)))
  }
  d <- tibble(sp = names(x), P = p$P, S = p$S) |>
    tidyr::pivot_longer(c(P, S), names_to = "comp", values_to = "v") |>
    mutate(comp = ifelse(comp == "P", "P: valores estimados pelos autovetores", "S: resíduos"),
           sp = factor(sp, levels = rev(names(x))))
  ggplot(d, aes(v, sp, fill = comp)) +
    geom_vline(xintercept = 0, colour = cores[["neutro"]]) +
    geom_col(width = .7) +
    facet_wrap(~comp, scales = "free_x") +
    scale_fill_manual(values = c(cores[["ocre"]], cores[["verde"]]), guide = "none") +
    labs(x = NULL, y = NULL) + tema_aula(16)
}

#' Autovalores de árvores balanceada e pectinada (ideia de Martins &
#' Housworth 2002)
fig_autovalores_formas <- function() {
  bal <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  pec <- ape::read.tree(text = "(a:7,(b:6,(c:5,(d:4,(e:3,(f:2,(g:1,h:1):1):1):1):1):1):1);")
  ev <- function(tr, nome) {
    x <- setNames(rep(0, 8), tr$tip.label)
    p <- pvr_passos(tr, x + rnorm(8), 1)
    tibble(forma = nome, k = seq_along(p$valores), pct = 100 * p$valores / sum(p$valores))
  }
  set.seed(1)
  d <- bind_rows(ev(bal, "balanceada"), ev(pec, "pectinada"))
  a1 <- plot_arvore(bal, rotulos = FALSE) + labs(title = "balanceada")
  a2 <- plot_arvore(pec, rotulos = FALSE) + labs(title = "pectinada")
  tt <- theme(plot.title = element_text(hjust = .5, size = 17, face = "bold", colour = cores[["texto"]]))
  b <- ggplot(d, aes(factor(k), pct, fill = forma)) +
    geom_col(position = "dodge", width = .7) +
    scale_fill_manual(values = cores_linhagem[c(1, 2)], name = NULL) +
    labs(x = "autovetor", y = "% da variação filogenética") + tema_aula(15)
  ((a1 + tt) / (a2 + tt)) | b
}

#' Curva PSR do SVL dos Anolis
fig_psr <- function() {
  an <- anolis()
  d <- psr_curva(an$arvore, an$svl)
  ggplot(d, aes(repr, R2)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = cores[["texto"]]) +
    geom_ribbon(aes(ymin = pmin(R2, repr), ymax = pmax(R2, repr)), fill = cores[["neutro_cl"]]) +
    geom_line(colour = cores[["accent"]], linewidth = 1.2) +
    geom_point(colour = cores[["accent"]], size = 1.8) +
    annotate("text", x = .7, y = .62, label = "esperado sob BM", angle = 38, size = 5.5) +
    scale_x_continuous(labels = \(v) paste0(100 * v, " %")) +
    labs(x = "representação filogenética (autovalores acumulados)", y = "R² do PVR (SVL)") +
    coord_equal() + tema_aula(15)
}

#' Partição da variação entre filogenia e ambiente (Desdevises et al. 2003)
fig_particao <- function() {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 1, fill = cores[["papel"]], colour = cores[["texto"]]) +
    annotate("rect", xmin = 0, xmax = 6.5, ymin = 1.4, ymax = 1.9, fill = cores[["ocre"]], alpha = .8) +
    annotate("rect", xmin = 2.5, xmax = 8, ymin = 2.2, ymax = 2.7, fill = cores[["verde"]], alpha = .8) +
    annotate("text", x = 3.25, y = 1.65, label = "explicada pela filogenia", colour = "white", size = 6, fontface = "bold") +
    annotate("text", x = 5.25, y = 2.45, label = "explicada pelo ambiente", colour = "white", size = 6, fontface = "bold") +
    annotate("segment", x = c(2.5, 6.5, 8), xend = c(2.5, 6.5, 8), y = 0, yend = 1, linetype = "22") +
    annotate("text", x = c(1.25, 4.5, 7.25, 9), y = .5, label = c("[a]", "[b]", "[c]", "[d]"), size = 8, fontface = "bold") +
    annotate("text", x = c(1.25, 4.5, 7.25, 9), y = -.35,
             label = c("só filogenia", "compartilhada", "só ambiente", "não explicada"), size = 5) +
    coord_cartesian(xlim = c(-.2, 10.2), ylim = c(-.6, 2.9)) + theme_void() +
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}

#' Bocaina (Garey et al. 2018): PVR da frequência dominante (log)
pvr_bocaina <- function(k = 3) {
  tb <- ape::read.tree(here::here("aula-03-sinal/filogenia_bocaina.txt"))
  cn <- read.csv(here::here("aula-03-sinal/canto.csv"), sep = ";", row.names = 1)
  x <- setNames(log(cn$F_dom), rownames(cn))[tb$tip.label]
  list(arvore = tb, x = x, pvr = pvr_passos(tb, x, k = k))
}
fig_pvr_bocaina <- function(b = pvr_bocaina()) {
  tr <- b$arvore
  tr$tip.label <- gsub("_", " ", tr$tip.label)
  d <- tibble(sp = gsub("_", " ", names(b$x)), P = b$pvr$P - mean(b$x), S = b$pvr$S)
  cd <- coord_arvore(tr); pt <- dplyr::filter(cd$nos, ponta) |> left_join(d, by = c("rotulo" = "sp"))
  arv <- plot_arvore(tr, base = 3.6, folga = .55) + theme(plot.margin = margin(0, 0, 0, 0))
  barras <- pt |> tidyr::pivot_longer(c(P, S), names_to = "comp", values_to = "v") |>
    mutate(comp = ifelse(comp == "P", "P (filogenético)", "S (específico)")) |>
    ggplot(aes(v, y, fill = comp)) +
    geom_vline(xintercept = 0, colour = cores[["neutro"]]) +
    geom_col(orientation = "y", width = .7) +
    facet_wrap(~comp, scales = "free_x") +
    scale_fill_manual(values = c(cores[["ocre"]], cores[["verde"]]), guide = "none") +
    scale_y_continuous(expand = expansion(add = .5)) +
    labs(x = "log da frequência dominante, centrado", y = NULL) + tema_aula(13) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
  arv + barras + plot_layout(widths = c(1.3, 1))
}
