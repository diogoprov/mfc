# Aula 2, bloco 5: contrastes independentes (Felsenstein 1985), passo a
# passo na árvore-guia, com o atributo A do exemplo de Freckleton et al.
library(patchwork)
library(nlme)
source(here::here("R/tema.R"))
source(here::here("R/arvore.R"))

#' O algoritmo dos contrastes, nó a nó (em pós-ordem, como o `ape::pic`)
#'
#' @return tibble com um nó por linha: filhos, valores e ramos dos filhos,
#'   contraste bruto e padronizado, valor estimado do nó e acréscimo ao
#'   ramo que leva ao nó
pic_passos <- function(tr = arvore_guia(), x = setNames(atributos_guia()$A, atributos_guia()$sp)) {
  n <- ape::Ntip(tr)
  val <- numeric(max(tr$edge)); val[1:n] <- x[tr$tip.label]
  v   <- numeric(max(tr$edge))           # comprimento (corrigido) do ramo acima de cada nó
  v[tr$edge[, 2]] <- tr$edge.length
  po <- ape::reorder.phylo(tr, "postorder")$edge
  nos <- unique(po[, 1])
  out <- purrr::map_dfr(nos, \(k) {
    f <- tr$edge[tr$edge[, 1] == k, 2]
    i <- f[1]; j <- f[2]
    bruto <- val[i] - val[j]
    pad   <- bruto / sqrt(v[i] + v[j])
    val[k] <<- (val[i] / v[i] + val[j] / v[j]) / (1 / v[i] + 1 / v[j])
    extra <- v[i] * v[j] / (v[i] + v[j])
    if (k != n + 1) v[k] <<- v[k] + extra
    tibble(no = k, i = i, j = j, xi = val[i], xj = val[j], vi = v[i], vj = v[j],
           bruto = bruto, contraste = pad, valor_no = val[k], extra = extra,
           dp = sqrt(v[i] + v[j]))
  })
  out
}

rot_no <- function(k, tr = arvore_guia()) {
  n <- ape::Ntip(tr); ifelse(k <= n, tr$tip.label[k], paste0("nó ", k))
}

#' A árvore-guia com o atributo A e os passos dos contrastes
#'
#' @param etapa 0 = árvore com os valores; 1 = contraste sp1–sp2;
#'   2 = valor do nó 10; 3 = ramo do nó 10 alongado; 4 = todos os nós
fig_pic <- function(etapa = 0) {
  tr <- arvore_guia(); a <- atributos_guia()
  ps <- pic_passos()
  tr2 <- tr; tr2$tip.label <- paste0(a$sp, "  (", format(a$A, nsmall = 2, decimal.mark = ","), ")")
  cr <- rep(cores[["texto"]], nrow(tr$edge)); lg <- rep(1.3, nrow(tr$edge))
  rot <- NULL
  f <- function(z) format(round(z, 2), nsmall = 2, decimal.mark = ",")
  if (etapa >= 1) {
    e <- which(tr$edge[, 2] %in% c(1, 2)); cr[e] <- cores[["ocre"]]; lg[e] <- 2.8
  }
  if (etapa >= 2) rot <- c(`10` = f(ps$valor_no[ps$no == 10]))
  if (etapa >= 3) { e <- which(tr$edge[, 2] == 10); cr[e] <- cores[["ferrugem"]]; lg[e] <- 2.8 }
  if (etapa >= 4) {
    cr[] <- cores[["texto"]]; lg[] <- 1.3
    rot <- setNames(f(ps$valor_no), ps$no)
  }
  g <- plot_arvore(tr2, cor_ramo = cr, larg = lg, comprimentos = etapa %in% 1:3,
                   rot_nos = rot, base = 5.6, folga = .3)
  if (etapa == 3) {
    cd <- coord_arvore(tr)
    y10 <- cd$nos$y[10]; x10 <- cd$nos$x[10]
    g <- g + annotate("text", x = x10 - .45, y = y10 - .5, colour = cores[["ferrugem"]], size = 5.4,
                      label = paste0("1,06 + ", f(ps$extra[ps$no == 10]), " = ",
                                     f(1.06 + ps$extra[ps$no == 10])))
  }
  g
}

#' Tabela dos contrastes de A e de B (com rótulos dos nós)
tabela_pic <- function() {
  tr <- arvore_guia(); a <- atributos_guia()
  pa <- pic_passos(tr, setNames(a$A, a$sp)); pb <- pic_passos(tr, setNames(a$B, a$sp))
  tibble(`Nó` = pa$no,
         `Contraste entre` = paste(rot_no(pa$i), "e", rot_no(pa$j)),
         `Contraste de A` = pa$contraste, `Contraste de B` = pb$contraste)
}

#' Regressão dos contrastes pela origem
fig_pic_reg <- function() {
  d <- tabela_pic()
  b <- coef(lm(`Contraste de B` ~ `Contraste de A` - 1, d))
  ggplot(d, aes(`Contraste de A`, `Contraste de B`)) +
    geom_hline(yintercept = 0, colour = cores[["neutro_cl"]]) +
    geom_vline(xintercept = 0, colour = cores[["neutro_cl"]]) +
    geom_abline(slope = b, intercept = 0, colour = cores[["accent"]], linewidth = 1.2) +
    geom_point(size = 4.5, shape = 21, fill = cores[["ocre"]], colour = "white", stroke = .8) +
    geom_text(aes(label = paste0("nó ", `Nó`)), vjust = -1.1, size = 4.8) +
    scale_y_continuous(expand = expansion(mult = c(.08, .15))) +
    scale_x_continuous(expand = expansion(mult = .12)) +
    coord_cartesian(clip = "off") +
    labs(x = "contraste de A", y = "contraste de B") + tema_aula(17)
}

#' Dados brutos × contrastes, com as inclinações
fig_bruto_vs_pic <- function() {
  a <- atributos_guia() |> mutate(clado = ifelse(sp %in% paste0("sp", 1:4), "sp1–sp4", "sp5–sp7"))
  b_ols <- coef(lm(B ~ A, a))[2]
  d <- tabela_pic(); b_pic <- coef(lm(`Contraste de B` ~ `Contraste de A` - 1, d))
  f <- function(z) format(round(z, 2), nsmall = 2, decimal.mark = ",")
  g1 <- ggplot(a, aes(A, B)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = cores[["ferrugem"]], linewidth = 1.1) +
    geom_point(aes(fill = clado), size = 4.5, shape = 21, colour = "white", stroke = .8) +
    geom_text(aes(label = sp), vjust = -1.1, size = 4.6) +
    scale_fill_manual(values = cores_linhagem[c(1, 3)], name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(.08, .15))) +
    labs(title = paste0("Espécies (OLS): inclinação = ", f(b_ols)), x = "A", y = "B") + tema_aula(16)
  g2 <- fig_pic_reg() + labs(title = paste0("Contrastes: inclinação = ", f(b_pic))) + tema_aula(16)
  g1 + g2
}

#' Diagnóstico: |contraste| × desvio padrão (√ da soma dos ramos)
fig_pic_diag <- function() {
  ps <- pic_passos()
  ggplot(ps, aes(dp, abs(contraste))) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = cores[["neutro"]], linetype = "22") +
    geom_point(size = 4.5, shape = 21, fill = cores[["accent"]], colour = "white", stroke = .8) +
    labs(x = "desvio padrão do contraste (√(vᵢ + vⱼ))", y = "|contraste| de A") + tema_aula(17)
}
