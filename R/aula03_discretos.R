# Aula 3, bloco 8: sinal em atributos categóricos (decomposição de Rao,
# D de Fritz & Purvis, δ).
source(here::here("R/aula03_base.R"))

#' Os quatro passos da decomposição da entropia quadrática (Pavoine et al. 2010)
fig_rao_passos <- function() {
  cx <- tibble::tribble(
    ~x, ~y, ~w, ~h, ~rotulo,
    0,  0, 2.6, 1.4, "Entropia quadrática\nde Rao do conjunto\nde espécies",
    3.2, 0, 2.6, 1.4, "Decompõe a entropia\nnos nós da filogenia\n(só a topologia)",
    6.4, 0, 2.6, 1.4, "Testa se a diversidade\nse concentra perto da\nraiz ou das pontas",
    9.6, 0, 2.6, 1.4, "Aleatorização\npara a\nsignificância")
  st <- tibble::tribble(~x, ~y, ~xend, ~yend, 1.35, 0, 1.85, 0, 4.55, 0, 5.05, 0, 7.75, 0, 8.25, 0)
  diagrama(cx, st, base = 4.6)
}

#' Duas árvores de 8 espécies: diversidade concentrada perto da raiz (A) ou
#' das pontas (B)
fig_rao_arvores <- function() {
  tr <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
  estados <- list(A = c("1", "1", "1", "1", "2", "2", "2", "2"),
                  B = c("1", "2", "1", "2", "1", "2", "1", "2"))
  cor_e <- c("1" = cores_linhagem[1], "2" = cores_linhagem[2])
  um <- function(nome, lado) {
    cd <- coord_arvore(tr); pt <- dplyr::filter(cd$nos, ponta)
    pt$e <- estados[[nome]][match(pt$rotulo, tr$tip.label)]
    plot_arvore(tr, rotulos = FALSE, base = 5) +
      geom_point(data = pt, aes(x + .25, y, colour = e), size = 7) +
      scale_colour_manual(values = cor_e, guide = "none") +
      labs(title = paste0(nome, ": diversidade perto ", lado)) +
      theme(plot.title = element_text(hjust = .5, size = 17, face = "bold", colour = cores[["texto"]]))
  }
  um("A", "da raiz") + um("B", "das pontas")
}

#' D de Fritz & Purvis (2010) para um atributo binário dos Anolis:
#' ecomorfo trunk-ground (TG) ou não
d_anolis <- function(nperm = 1000, seed = 331) {
  an <- anolis()
  dados <- data.frame(sp = names(an$ecomorfo), tg = as.integer(an$ecomorfo == "TG"))
  cd <- caper::comparative.data(an$arvore, dados, names.col = "sp")
  set.seed(seed)
  r <- caper::phylo.d(cd, binvar = tg, permut = nperm)
  list(D = r$DEstimate, p_aleatorio = r$Pval1, p_bm = r$Pval0, n_tg = sum(dados$tg), n = nrow(dados))
}

#' Régua de leitura do D
fig_regua_d <- function(d = NULL) {
  marcas <- tibble(x = c(0, 1), r = c("D = 0", "D = 1"),
                   t = c("como o esperado sob BM\ncom limiar (agrupado)", "distribuição aleatória\nnas pontas"))
  g <- ggplot() +
    annotate("segment", x = -1.1, xend = 1.6, y = 0, yend = 0, linewidth = 1.2, colour = cores[["texto"]],
             arrow = arrow(length = unit(.3, "cm"), ends = "both", type = "closed")) +
    geom_point(data = marcas, aes(x, 0), size = 5, colour = cores[["texto"]]) +
    geom_text(data = marcas, aes(x, .35, label = r), size = 7, fontface = "bold", colour = cores[["texto"]]) +
    geom_text(data = marcas, aes(x, -.45, label = t), size = 5, lineheight = .95, colour = cores[["texto"]]) +
    annotate("text", x = -.75, y = .35, label = "D < 0:\nmais agrupado", size = 4.8, colour = cores[["accent"]], lineheight = .9) +
    annotate("text", x = 1.45, y = .35, label = "D > 1:\nsobredisperso", size = 4.8, colour = cores[["ferrugem"]], lineheight = .9)
  if (!is.null(d))
    g <- g + annotate("point", x = d$D, y = 0, size = 7, shape = 25, fill = cores[["ocre"]], colour = cores[["ocre"]]) +
      annotate("text", x = d$D, y = -.95, hjust = 0, size = 5.5, colour = cores[["ocre"]], fontface = "bold",
               label = paste0("Anolis, trunk-ground: D = ", num(d$D)))
  g + coord_cartesian(xlim = c(-1.2, 1.75), ylim = c(-1.15, .75)) + theme_void() +
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}
