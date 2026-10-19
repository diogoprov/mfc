# Aula 2, abertura e bloco 1 (Do processo ao padrão).
library(patchwork)
source(here::here("R/tema.R"))
source(here::here("R/diagrama.R"))
source(here::here("R/arvore.R"))
source(here::here("R/bm.R"))

#' Árvore de dois clados para a recapitulação e o erro tipo I
#'
#' Dois clados de `n` espécies cada, com ramos internos longos e pontas
#' curtas: o "pior caso" de Felsenstein (1985).
arvore_dois_clados <- function(n = 10, tronco = .9, pontas = .1) {
  # sprintf, e não paste0: com options(OutDec = ","), paste0 escreveria 0,1
  clado <- function(pref) paste0("(", paste(sprintf("%s%d:%g", pref, seq_len(n), pontas), collapse = ","), ")")
  ape::read.tree(text = sprintf("(%s:%g,%s:%g);", clado("a"), tronco, clado("b"), tronco))
}

#' Dois atributos independentes simulados por BM na árvore de dois clados
dados_dois_clados <- function(seed = 11) {
  tr <- arvore_dois_clados()
  set.seed(seed)
  tibble::tibble(sp = tr$tip.label,
                 x = ape::rTraitCont(tr, "BM", sigma = 1)[tr$tip.label],
                 y = ape::rTraitCont(tr, "BM", sigma = 1)[tr$tip.label],
                 clado = ifelse(substr(sp, 1, 1) == "a", "clado A", "clado B"))
}

#' Recapitulação da Aula 1: a correlação que vem só do parentesco
#'
#' @param etapa 1 = dispersão com a reta de OLS; 2 = os pontos coloridos
#'   pelos dois clados, com a árvore ao lado
fig_recap_clados <- function(etapa = 1) {
  d  <- dados_dois_clados()
  cc <- c("clado A" = cores_linhagem[1], "clado B" = cores_linhagem[2])
  g <- ggplot(d, aes(x, y)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = cores[["neutro"]], linewidth = 1) +
    labs(x = "atributo X", y = "atributo Y") +
    tema_aula(16) + theme(axis.text = element_blank())
  if (etapa == 1) return(g + geom_point(size = 4, shape = 21, fill = cores[["accent"]],
                                         colour = "white", stroke = .8))
  g <- g + geom_point(aes(fill = clado), size = 4, shape = 21, colour = "white", stroke = .8) +
    scale_fill_manual(values = cc, name = NULL)
  tr <- arvore_dois_clados()
  cor_r <- ifelse(tr$edge[, 2] %in% c(which(substr(tr$tip.label, 1, 1) == "a"),
                                      ape::getMRCA(tr, tr$tip.label[1:10])),
                  cc[1], cc[2])
  arv <- plot_arvore(tr, cor_ramo = cor_r, rotulos = FALSE, larg = 1.1)
  g + arv + plot_layout(widths = c(1.4, 1))
}

#' Hansen & Martins (1996): dos mecanismos microevolutivos aos dados
#'
#' @param etapa 1 = o caminho conceitual (para baixo); 2 = acrescenta o
#'   caminho inverso, com "?"
fig_hansen_martins <- function(etapa = 1) {
  sim <- sim_linhagens(t_div = c(0, 35, 60), mae = c(NA, 1, 1), sigma2 = 1, seed = 4)
  fim <- sim |> group_by(linhagem) |> slice_tail(n = 1) |> ungroup()
  fim$x <- separa_rotulos(fim$t, fim$x, dx = 10, dy = 3)
  traj <- ggplot(sim, aes(t, x, colour = factor(linhagem))) +
    geom_line(linewidth = .9) +
    geom_text(data = fim, aes(label = paste0("X[", linhagem, "](t)")), parse = TRUE,
              hjust = -.2, size = 5.5, show.legend = FALSE) +
    scale_colour_manual(values = cores_linhagem, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(.02, .15))) +
    labs(x = "tempo", y = "fenótipo") +
    tema_aula(16, conceitual = TRUE)
  caixas <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo,
    0, 4, 4.6, 1.3, "Mecanismos microevolutivos\n(deriva, seleção, mutação, fluxo gênico)",
    0, 2, 4.6, .9,  "Modelos evolutivos",
    0, 0, 4.6, .9,  "Dados interespecíficos")
  setas <- tibble::tribble(~x, ~y, ~xend, ~yend, -.6, 3.3, -.6, 2.5, -.6, 1.5, -.6, .5)
  dg <- diagrama(caixas, setas, base = 5)
  if (etapa == 2) {
    dg <- dg +
      annotate("segment", x = 2.75, y = .2, xend = 2.75, yend = 3.6, linewidth = 2.2,
               colour = cores[["ferrugem"]],
               arrow = arrow(length = unit(.35, "cm"), type = "closed")) +
      annotate("text", x = 3.25, y = 1.9, label = "?", size = 16, fontface = "bold",
               colour = cores[["ferrugem"]])
  }
  traj + dg + plot_layout(widths = c(1.15, 1))
}

#' Os quatro modelos de evolução fenotípica, lado a lado (prévia da Aula 5)
fig_quatro_modelos <- function(n = 25) {
  d <- bind_rows(
    sim_bm_tempo(n, seed = 2) |> mutate(modelo = "BM"),
    sim_eb_tempo(n, seed = 2) |> mutate(modelo = "EB"),
    sim_trend_tempo(n, seed = 2) |> mutate(modelo = "Tendência"),
    sim_ou_tempo(n, alfa = .15, theta = 6, seed = 2) |> mutate(modelo = "OU")) |>
    mutate(modelo = factor(modelo, levels = c("BM", "EB", "Tendência", "OU")))
  ggplot(d, aes(t, x, group = rep)) +
    geom_line(colour = cores[["accent"]], alpha = .45, linewidth = .5) +
    facet_wrap(~modelo, nrow = 1) +
    labs(x = "tempo", y = "atributo") +
    tema_aula(16, conceitual = TRUE)
}
