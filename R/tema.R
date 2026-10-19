library(ggplot2)

# Vírgula decimal nos eixos e em format(): os slides são em português.
options(OutDec = ",")

# Paleta do Biodiversity Synthesis Lab, a mesma do site e do tema dos slides
# (copiada de ~/GitHub/multidimensionais/R/tema.R).
cores <- c(
  accent      = "#1c6e8c",
  accent_dark = "#144f64",
  verde       = "#527238",
  ocre        = "#96661f",
  ferrugem    = "#a15e3f",
  neutro      = "#8a8172",
  neutro_cl   = "#e6e0d3",
  texto       = "#2d3748",
  papel       = "#faf7f0"
)

# Identidade visual do mapa dos métodos: uma cor por pergunta, usada no
# mapa (R/mapa.R) e nos selos dos slides ($familia-* em slides/theme.scss;
# mantenha as duas listas iguais). Todas dão contraste de pelo menos 5:1
# com texto branco. "Base" é a faixa das Aulas 1 e 2.
cores_familia <- c(
  "Base"                  = "#6f5f45",
  "Há sinal?"             = "#3d7560",
  "Como evolui?"          = "#a3472e",
  "Como era o ancestral?" = "#6b4c8a",
  "Estão associados?"     = "#1c6e8c"
)

# Cores fixas das linhagens nas simulações (a mesma linhagem, a mesma cor,
# em todas as figuras).
cores_linhagem <- c("#1c6e8c", "#a3472e", "#3d7560", "#96661f", "#6b4c8a",
                    "#8c3b5e", "#4a5568")

#' Tema base dos gráficos das aulas
#'
#' @param base tamanho base do texto (pt); 15 fica legível projetado
#' @param conceitual TRUE apaga números dos eixos e desenha setas nas
#'   pontas, para gráficos que ilustram uma ideia e não um conjunto de dados
tema_aula <- function(base = 15, conceitual = FALSE) {
  t <- theme_minimal(base_size = base) +
    theme(
      plot.background  = element_rect(fill = cores[["papel"]], colour = NA),
      panel.background = element_rect(fill = cores[["papel"]], colour = NA),
      panel.grid.minor = element_blank(),
      text             = element_text(colour = cores[["texto"]]),
      axis.text        = element_text(colour = "#4a5568"),
      strip.text       = element_text(face = "bold", size = rel(1),
                                      colour = cores[["texto"]]),
      plot.title       = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle    = element_text(colour = "#4a5568"),
      legend.position  = "bottom"
    )
  if (conceitual) {
    t <- t + theme(
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.line  = element_line(
        colour = cores[["texto"]], linewidth = .7,
        arrow = grid::arrow(length = grid::unit(.22, "cm"), type = "closed"))
    )
  }
  t
}

#' Afasta rótulos que se sobrepõem, empurrando-os na vertical
#'
#' Para os nomes nas pontas das setas de biplots, que se amontoam quando as
#' setas são quase paralelas. Sem pacote extra (ggrepel).
#'
#' @param x,y posições dos rótulos
#' @param dx,dy largura e altura aproximadas de um rótulo, nas unidades
#'   dos dados
#' @return y ajustado
separa_rotulos <- function(x, y, dx, dy, iter = 100) {
  for (k in seq_len(iter)) {
    mexeu <- FALSE
    for (i in seq_along(y)) for (j in seq_along(y)) {
      if (i < j && abs(x[i] - x[j]) < dx && abs(y[i] - y[j]) < dy) {
        d <- (dy - abs(y[i] - y[j])) / 2 + 1e-9
        s <- if (y[i] >= y[j]) 1 else -1
        y[i] <- y[i] + s * d; y[j] <- y[j] - s * d
        mexeu <- TRUE
      }
    }
    if (!mexeu) break
  }
  y
}
