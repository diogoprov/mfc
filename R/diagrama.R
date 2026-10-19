library(tibble)
library(ggplot2)
source(here::here("R/mapa.R"))   # ret_arred(), cores, cores_familia

#' Diagrama de caixas ligadas por setas
#'
#' Base dos fluxogramas das aulas (árvores de decisão, passos de um
#' método). As coordenadas são livres; `coord_equal()` mantém as caixas
#' com a proporção desenhada.
#'
#' @param caixas tibble com x, y (centro), w, h, rotulo e, opcionais, fill
#'   (fundo), cor (borda e texto quando o fundo é claro), txt (cor do
#'   texto), tam (multiplicador do texto) e face ("plain"/"bold")
#' @param setas tibble com x, y, xend, yend e, opcional, tipo ("solid",
#'   "22" para tracejada) e curv (curvatura; 0 = reta)
#' @param base tamanho base do texto
diagrama <- function(caixas, setas = NULL, base = 4.6) {
  cx <- caixas
  if (!"fill" %in% names(cx)) cx$fill <- cores[["papel"]]
  if (!"cor" %in% names(cx))  cx$cor  <- cores[["accent_dark"]]
  if (!"txt" %in% names(cx))  cx$txt  <- cx$cor
  if (!"tam" %in% names(cx))  cx$tam  <- 1
  if (!"face" %in% names(cx)) cx$face <- "plain"
  cx$id <- seq_len(nrow(cx))
  poli <- do.call(rbind, lapply(cx$id, \(i)
    transform(ret_arred(cx$x[i], cx$y[i], cx$w[i], cx$h[i],
                        r = min(.12, cx$h[i] / 3)), id = i)))
  poli$fill <- cx$fill[poli$id]
  poli$cor  <- cx$cor[poli$id]

  g <- ggplot()
  if (!is.null(setas) && nrow(setas) > 0) {
    if (!"tipo" %in% names(setas)) setas$tipo <- "solid"
    if (!"curv" %in% names(setas)) setas$curv <- 0
    for (i in seq_len(nrow(setas))) {
      g <- g + annotate("curve", x = setas$x[i], y = setas$y[i],
                        xend = setas$xend[i], yend = setas$yend[i],
                        curvature = setas$curv[i], linetype = setas$tipo[i],
                        linewidth = .7, colour = cores[["accent_dark"]],
                        arrow = arrow(length = unit(.18, "cm"), type = "closed"))
    }
  }
  g +
    geom_polygon(data = poli, aes(px, py, group = id),
                 fill = poli$fill, colour = poli$cor, linewidth = .6) +
    geom_text(data = cx, aes(x, y, label = rotulo), colour = cx$txt,
              size = base * cx$tam, fontface = cx$face, lineheight = .95) +
    coord_equal(clip = "off") +
    theme_void() +
    theme(plot.background  = element_rect(fill = cores[["papel"]], colour = NA),
          panel.background = element_rect(fill = cores[["papel"]], colour = NA),
          plot.margin = margin(4, 4, 4, 4))
}
