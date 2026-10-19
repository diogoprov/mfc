library(dplyr)
library(tibble)
library(ggplot2)
source(here::here("R/tema.R"))

#' Desenha uma matriz como tabela, com células destacadas
#'
#' Para montar matrizes passo a passo (C, distâncias patrísticas, Q, P(t)).
#'
#' @param M matriz numérica ou de texto (os nomes de linhas e colunas viram
#'   os rótulos)
#' @param destaque matriz lógica do tamanho de M (ou NULL): células acesas
#' @param cor_destaque cor das células acesas
#' @param fill matriz de cores de fundo (opcional; sobrepõe-se ao padrão)
#' @param ocultar matriz lógica: células ainda vazias (para revelar aos poucos)
#' @param digitos casas decimais dos números
#' @param base tamanho do texto
#' @param colchetes desenhar colchetes, como numa matriz matemática
plot_matriz <- function(M, destaque = NULL, cor_destaque = cores[["ocre"]],
                        fill = NULL, ocultar = NULL, digitos = 2, base = 5.5,
                        colchetes = FALSE, rot_linhas = rownames(M),
                        rot_colunas = colnames(M), cor_texto = NULL) {
  n <- nrow(M); m <- ncol(M)
  txt <- if (is.numeric(M)) format(round(M, digitos), nsmall = digitos, decimal.mark = ",") else M
  d <- tibble(i = rep(seq_len(n), m), j = rep(seq_len(m), each = n),
              txt = as.vector(txt))
  d$fill <- if (is.null(fill)) "#f3eee3" else as.vector(fill)
  if (!is.null(destaque)) d$fill[as.vector(destaque)] <- cor_destaque
  d$cor <- ifelse(d$fill %in% c("#f3eee3", cores[["papel"]], "white"), cores[["texto"]], "white")
  if (!is.null(cor_texto)) d$cor <- as.vector(cor_texto)
  if (!is.null(ocultar)) d$txt[as.vector(ocultar)] <- ""
  g <- ggplot(d, aes(j, -i)) +
    geom_tile(fill = d$fill, colour = cores[["papel"]], linewidth = 1.2) +
    geom_text(aes(label = txt), colour = d$cor, size = base)
  if (!is.null(rot_linhas))
    g <- g + annotate("text", x = .4, y = -seq_len(n), label = rot_linhas,
                      hjust = 1, size = base * .9, colour = cores[["texto"]], fontface = "bold")
  if (!is.null(rot_colunas))
    g <- g + annotate("text", x = seq_len(m), y = -.35, label = rot_colunas,
                      vjust = 0, size = base * .9, colour = cores[["texto"]], fontface = "bold")
  if (colchetes) {
    for (s in c(-1, 1)) {
      x0 <- if (s < 0) .45 else m + .55
      g <- g + annotate("path", x = c(x0 - s * .15, x0, x0, x0 - s * .15),
                        y = c(-.45, -.45, -n - .55, -n - .55),
                        linewidth = 1, colour = cores[["texto"]])
    }
  }
  g + coord_equal(clip = "off") + theme_void() +
    theme(plot.background  = element_rect(fill = cores[["papel"]], colour = NA),
          panel.background = element_rect(fill = cores[["papel"]], colour = NA),
          plot.margin = margin(10, 10, 10, 30))
}
