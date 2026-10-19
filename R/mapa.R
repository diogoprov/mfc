library(purrr)
library(tibble)
library(ggplot2)
library(dplyr)
source(here::here("R/tema.R"))   # cores, cores_familia

# Retângulo de cantos arredondados como polígono, para geom_polygon().
ret_arred <- function(x, y, w, h, r = .16, n = 12) {
  t <- seq(0, pi / 2, length.out = n)
  cx <- c(x + w / 2 - r, x - w / 2 + r, x - w / 2 + r, x + w / 2 - r)
  cy <- c(y + h / 2 - r, y + h / 2 - r, y - h / 2 + r, y - h / 2 + r)
  purrr::map_dfr(1:4, \(k) tibble::tibble(
    px = cx[k] + r * cos(t + (k - 1) * pi / 2),
    py = cy[k] + r * sin(t + (k - 1) * pi / 2)))
}

# Células do mapa: pergunta (coluna) × tipo de atributo (linha), com os
# métodos (do mapa mental do Coggle do Diogo) e a(s) aula(s) em que entram.
mapa_celulas <- tibble::tribble(
  ~col, ~lin, ~metodos,                                            ~aulas,
  1, 3, "K de Blomberg · λ de Pagel\nI de Moran · Mantel\nPVR",    list(3),
  1, 2, "D de Fritz · estatística δ",                              list(3),
  1, 1, "K multivariado",                                          list(8),
  2, 3, "BM · OU · EB · tendência\nmulti-taxas\nmulti-regime",     list(5),
  2, 2, "Mk: ER · SYM · ARD\ntaxas ocultas",                       list(6),
  2, 1, "mvMORPH\ntaxas por módulo\nPhylogeneticEM",               list(8),
  3, 3, "máxima verossimilhança\ncontMap",                         list(7),
  3, 2, "marginal · conjunta\nmapeamento\nestocástico",             list(7),
  3, 1, "filomorfoespaço",                                         list(c(7, 8)),
  4, 3, "PIC · PGLS\nANOVA filogenética\nGLM filogenético",        list(c(2, 4)),
  4, 2, "Pagel (1994) · threshold\nregressão logística\nfilogenética", list(6),
  4, 1, "PGLS multivariado\nPCA filogenética",                     list(8)
)

#' Mapa dos métodos: o fio condutor da disciplina
#'
#' Colunas = as perguntas (cada uma com a sua cor, `cores_familia`);
#' linhas = o tipo de atributo; embaixo, as duas faixas de base (Aula 1, a
#' árvore e os dados; Aula 2, o BM e a matriz C). As células da aula em
#' `aula` ficam com a cor cheia; as demais, esmaecidas.
#'
#' @param aula número(s) da(s) aula(s) a destacar; NULL destaca tudo
#' @param base tamanho base do texto
mapa_plot <- function(aula = NULL, base = 4.8) {
  perg <- names(cores_familia)[-1]
  sub  <- c("padrão", "processo", "história", "relação")
  tipos <- c("Multivariado", "Discreto", "Contínuo")   # lin = 1, 2, 3
  neutro <- "#ece6da"; neutro_txt <- "#6b6558"

  on <- function(a) if (is.null(aula)) TRUE else any(unlist(a) %in% aula)
  cel <- mapa_celulas |>
    mutate(x = col * 3, y = lin * 1.6 + 1.2,
           ligado = map_lgl(aulas, on),
           fill = ifelse(ligado, cores_familia[perg[col]], neutro),
           txt  = ifelse(ligado, "white", neutro_txt),
           rot_aula = map_chr(aulas, \(a) paste("Aula", paste(unlist(a), collapse = " e "))))
  caixas <- purrr::pmap_dfr(cel, \(x, y, fill, ...) ret_arred(x, y, 2.85, 1.3, r = .1) |>
                              mutate(id = paste(x, y), fill = fill))

  faixas <- tibble(y = c(1.78, 1.12), aula = c(2, 1),
                   rot = c("Aula 2 · O movimento browniano e a matriz C: a base de todas as colunas",
                           "Aula 1 · A árvore e os dados: obter a árvore, descompasso árvore–atributos, ler e interpretar")) |>
    mutate(ligado = if (is.null(aula)) TRUE else aula %in% !!aula,
           fill = ifelse(ligado, cores_familia[["Base"]], neutro),
           txt  = ifelse(ligado, "white", neutro_txt))
  pol_faixas <- purrr::pmap_dfr(faixas, \(y, fill, ...) ret_arred(7.5, y, 11.85, .55, r = .1) |>
                                  mutate(id = y, fill = fill))
  cab <- tibble(x = 1:4 * 3, y = 7.2, l = paste0(perg, "\n(", sub, ")"),
                fill = cores_familia[perg])
  pol_cab <- purrr::pmap_dfr(cab, \(x, y, fill, ...) ret_arred(x, y, 2.85, .8, r = .1) |>
                               mutate(id = x, fill = fill))

  ggplot() +
    geom_polygon(data = pol_cab, aes(px, py, group = id), fill = pol_cab$fill) +
    geom_text(data = cab, aes(x, y, label = l), colour = "white", size = base * .8,
              fontface = "bold", lineheight = .9) +
    geom_polygon(data = caixas, aes(px, py, group = id), fill = caixas$fill) +
    geom_text(data = cel, aes(x, y - .1, label = metodos), colour = cel$txt,
              size = base * .8, lineheight = .95) +
    geom_text(data = cel, aes(x + 1.32, y + .5, label = rot_aula), colour = cel$txt,
              size = base * .62, hjust = 1, fontface = "bold") +
    geom_text(data = tibble(y = 1:3 * 1.6 + 1.2, l = tipos), aes(x = 1.35, y, label = l),
              size = base * 1.02, fontface = "bold", hjust = 1, colour = cores[["texto"]]) +
    geom_polygon(data = pol_faixas, aes(px, py, group = id), fill = pol_faixas$fill) +
    geom_text(data = faixas, aes(7.5, y, label = rot), colour = faixas$txt,
              size = base * .85, fontface = "bold") +
    coord_cartesian(xlim = c(-0.9, 13.5), ylim = c(.78, 7.65), expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.background  = element_rect(fill = cores[["papel"]], colour = NA),
          panel.background = element_rect(fill = cores[["papel"]], colour = NA),
          plot.margin = margin(2, 4, 2, 4))
}
