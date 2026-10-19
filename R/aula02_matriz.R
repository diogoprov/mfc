# Aula 2, bloco 3: a matriz de covariância filogenética C, montada passo a
# passo na árvore-guia (exemplo de Freckleton et al. 2002).
library(patchwork)
source(here::here("R/tema.R"))
source(here::here("R/arvore.R"))
source(here::here("R/matriz.R"))

#' Árvore de 3 espécies do slide de distância patrística (A: 5; B e C: 3 + 2)
arvore_tres <- function() ape::read.tree(text = "(A:5,(B:2,C:2):3);")

#' Distância patrística, revelada célula a célula
#'
#' @param etapa 0 = matriz vazia; 1 = diagonal; 2 = + d(A,B) e d(A,C);
#'   3 = completa (+ d(B,C))
fig_patristica <- function(etapa = 3) {
  tr <- arvore_tres()
  D <- ape::cophenetic.phylo(tr)[c("A", "B", "C"), c("A", "B", "C")]
  mostra <- matrix(FALSE, 3, 3)
  if (etapa >= 1) diag(mostra) <- TRUE
  if (etapa >= 2) mostra[1, 2:3] <- mostra[2:3, 1] <- TRUE
  if (etapa >= 3) mostra[2, 3] <- mostra[3, 2] <- TRUE
  novo <- matrix(FALSE, 3, 3)
  if (etapa == 2) novo[1, 2:3] <- novo[2:3, 1] <- TRUE
  if (etapa == 3) novo[2, 3] <- novo[3, 2] <- TRUE
  cr <- rep(cores[["texto"]], nrow(tr$edge))
  if (etapa == 2) cr[c(caminho_raiz(tr, "A"), caminho_raiz(tr, "B"))] <- cores[["ocre"]]
  if (etapa == 3) cr[setdiff(c(caminho_raiz(tr, "B"), caminho_raiz(tr, "C")),
                             caminho_comum(tr, "B", "C"))] <- cores[["ocre"]]
  arv <- plot_arvore(tr, cor_ramo = cr, larg = ifelse(cr == cores[["ocre"]], 2.6, 1.4),
                     comprimentos = TRUE, base = 7, dig = 0)
  mat <- plot_matriz(D, destaque = novo, ocultar = !mostra, digitos = 0, base = 8)
  arv + mat + plot_layout(widths = c(1.1, 1))
}

#' Rótulos "var"/"cov" da matriz C (slide 37 do original)
matriz_simbolica <- function(n = 7, nomes = paste0("sp", 1:n)) {
  M <- matrix("cov", n, n, dimnames = list(nomes, nomes)); diag(M) <- "var"; M
}

#' A matriz C montada na árvore-guia
#'
#' @param etapa "simbolica" (var/cov), "blocos" (os dois grandes clados),
#'   "var1", "cov12", "cov14", "cov34" (caminho pintado e célula acesa) ou
#'   "completa"
fig_C <- function(etapa = "completa", base_mat = 5.2) {
  tr <- arvore_guia()
  C  <- ape::vcv(tr)
  cr <- rep(cores[["texto"]], nrow(tr$edge)); lg <- rep(1.3, nrow(tr$edge))
  dest <- matrix(FALSE, 7, 7); fill <- NULL; M <- C; ocultar <- NULL
  pares <- list(var1 = c(1, 1), cov12 = c(1, 2), cov14 = c(1, 4), cov34 = c(3, 4))
  if (etapa == "simbolica") {
    M <- matriz_simbolica()
    fill <- matrix(ifelse(diag(7) == 1, cores[["ferrugem"]], cores[["accent"]]), 7, 7)
  } else if (etapa == "blocos") {
    cl1 <- 1:4; cl2 <- 5:7
    f <- matrix("#f3eee3", 7, 7)
    f[cl1, cl1] <- cores_linhagem[1]; f[cl2, cl2] <- cores_linhagem[3]
    fill <- f
    cr[tr$edge[, 2] %in% c(1:4, 9:11)] <- cores_linhagem[1]
    cr[tr$edge[, 2] %in% c(5:7, 12:13)] <- cores_linhagem[3]
  } else if (etapa %in% names(pares)) {
    p <- pares[[etapa]]
    dest[p[1], p[2]] <- dest[p[2], p[1]] <- TRUE
    e <- if (p[1] == p[2]) caminho_raiz(tr, p[1]) else caminho_comum(tr, p[1], p[2])
    cr[e] <- cores[["ocre"]]; lg[e] <- 2.8
    # nas etapas anteriores à completa, só as células já calculadas aparecem
    ordem <- names(pares)
    feitas <- pares[seq_len(match(etapa, ordem))]
    vis <- matrix(FALSE, 7, 7)
    for (q in feitas) vis[q[1], q[2]] <- vis[q[2], q[1]] <- TRUE
    ocultar <- !vis
  }
  arv <- plot_arvore(tr, cor_ramo = cr, larg = lg, comprimentos = TRUE, base = 5.5)
  mat <- plot_matriz(M, destaque = dest, fill = fill, ocultar = ocultar, base = base_mat)
  arv + mat + plot_layout(widths = c(1, 1.15))
}

#' Covariância sob OU numa árvore ultramétrica de altura T
#'
#' Cᵢⱼ = σ²/(2α) · e^(−2α(T − sᵢⱼ)) · (1 − e^(−2α sᵢⱼ)), com sᵢⱼ o caminho
#' compartilhado (Hansen 1997, raiz em equilíbrio omitida).
vcv_ou <- function(tr, alfa = 2, sigma2 = 1) {
  S <- ape::vcv(tr); Tt <- max(diag(S))
  sigma2 / (2 * alfa) * exp(-2 * alfa * (Tt - S)) * (1 - exp(-2 * alfa * S))
}

#' A mesma árvore, três matrizes: BM, λ = 0,5 e OU
fig_tres_matrizes <- function() {
  tr <- arvore_guia(); C <- ape::vcv(tr)
  lam <- C * .5; diag(lam) <- diag(C)
  ou <- vcv_ou(tr, alfa = 2); ou <- ou / max(ou) * max(C)   # mesma escala
  para_long <- function(M, nome) tibble(i = rep(1:7, 7), j = rep(1:7, each = 7),
                                       v = as.vector(M) / max(M), modelo = nome)
  d <- bind_rows(para_long(C, "BM"),
                 para_long(lam, "λ = 0,5 (Aula 3)"),
                 para_long(ou, "OU, α = 2 (Aula 5)")) |>
    mutate(modelo = factor(modelo, levels = unique(modelo)))
  ggplot(d, aes(j, -i, fill = v)) +
    geom_tile(colour = cores[["papel"]], linewidth = .8) +
    scale_fill_gradient(low = "#f3eee3", high = cores[["accent_dark"]],
                        name = "covariância\n(relativa)", limits = c(0, 1)) +
    facet_wrap(~modelo) + coord_equal() +
    scale_x_continuous(breaks = 1:7, labels = paste0("sp", 1:7), position = "top") +
    scale_y_continuous(breaks = -(1:7), labels = paste0("sp", 1:7)) +
    labs(x = NULL, y = NULL) + tema_aula(15) +
    theme(panel.grid = element_blank(), legend.position = "right",
          axis.text.x = element_text(angle = 90, vjust = .5))
}
