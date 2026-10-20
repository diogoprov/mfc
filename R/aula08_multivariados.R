# Aula 8: atributos multivariados (morfometria geométrica, K e PGLS
# multivariados, mvMORPH, modularidade, PCA filogenética).
library(patchwork)
source(here::here("R/aula03_base.R"))   # tema, árvore, matriz, diagrama, anolis()

# ---- Dados ------------------------------------------------------------------

#' Plethodon (geomorph): 40 exemplares de duas espécies, 12 marcos, alinhados
pleth_ind <- function() {
  e <- new.env(); utils::data("plethodon", package = "geomorph", envir = e)
  p <- e$plethodon
  gp <- geomorph::gpagen(p$land, print.progress = FALSE)
  list(bruto = p$land, gpa = gp, links = p$links, especie = p$species, local = p$site)
}

#' Nove espécies de Plethodon com árvore (geomorph::plethspecies)
pleth_sp <- function() {
  e <- new.env(); utils::data("plethspecies", package = "geomorph", envir = e)
  gp <- geomorph::gpagen(e$plethspecies$land, print.progress = FALSE)
  list(gpa = gp, arvore = e$plethspecies$phy)
}

#' As análises multivariadas testadas no inventário (sementes fixas)
multivar_resultados <- function(pin = pleth_ind(), ps = pleth_sp()) {
  tr <- ps$arvore
  set.seed(801); k <- geomorph::physignal(ps$gpa$coords, tr, iter = 999, print.progress = FALSE)
  gdf <- geomorph::geomorph.data.frame(ps$gpa, phy = tr)
  set.seed(804); pg <- geomorph::procD.pgls(coords ~ log(Csize), phy = phy, data = gdf, iter = 999, print.progress = FALSE)
  Y <- geomorph::two.d.array(ps$gpa$coords)[tr$tip.label, ]
  gic <- sapply(c(BM = "BM", OU = "OU", EB = "EB"), \(m)
    mvMORPH::GIC(mvMORPH::mvgls(Y ~ 1, tree = tr, model = m, method = "PL-LOOCV"))$GIC)
  mods12 <- rep(c("A", "B"), c(5, 7)); mods11 <- rep(c("A", "B"), c(5, 6))
  set.seed(802); mt <- geomorph::modularity.test(pin$gpa$coords, mods12, CI = FALSE, iter = 999, print.progress = FALSE)
  set.seed(805); it <- geomorph::integration.test(pin$gpa$coords, partition.gp = mods12, iter = 999, print.progress = FALSE)
  set.seed(803); mr <- geomorph::compare.multi.evol.rates(ps$gpa$coords, gp = mods11, phy = tr, iter = 999, print.progress = FALSE)
  pt <- summary(pg)$table
  list(k = k$phy.signal, k_p = k$pvalue, k_nula = k$random.K,
       pgls_r2 = pt[1, "Rsq"], pgls_p = pt[1, "Pr(>F)"], gic = gic,
       cr = mt$CR, cr_p = mt$P.value, pls = it$r.pls, pls_p = it$P.value,
       taxa_mod = mr$sigma.d.gp, razao_mod = mr$sigma.d.ratio, razao_p = mr$P.value)
}

#' Anolis: log SVL e log cauda, BM com uma matriz de taxas × uma por ecomorfo
anolis_bmm <- function() {
  source(here::here("R/aula05_modelos.R"), local = TRUE)
  tr <- anolis_simmap()$arvore
  e <- new.env(); utils::data("anole.data", package = "phytools", envir = e)
  X <- as.matrix(e$anole.data[tr$tip.label, c("SVL", "TL")])   # já em log
  b1 <- mvMORPH::mvBM(tr, X, model = "BM1", diagnostic = FALSE, echo = FALSE)
  bm <- mvMORPH::mvBM(tr, X, model = "BMM", diagnostic = FALSE, echo = FALSE)
  list(aicc = c(BM1 = b1$AICc, BMM = bm$AICc), sigma = bm$sigma)
}

# ---- Bloco 1: morfometria geométrica ----------------------------------------------

desenha_marcos <- function(M, links, cor, rotulo = NULL) {
  pts <- tibble(x = M[, 1], y = M[, 2], conj = rotulo)
  seg <- tibble(x = M[links[, 1], 1], y = M[links[, 1], 2], xend = M[links[, 2], 1], yend = M[links[, 2], 2], conj = rotulo)
  list(geom_segment(data = seg, aes(x, y, xend = xend, yend = yend), colour = cor, linewidth = .9),
       geom_point(data = pts, aes(x, y), colour = cor, size = 3))
}

#' Superposição de Procrustes passo a passo: exemplar 2 sobre o exemplar 1
#' @param etapa 1 = brutos; 2 = centrados; 3 = tamanho do centroide = 1; 4 = rotação
fig_procrustes_passos <- function(pin = pleth_ind(), etapa = 1) {
  A <- pin$bruto[, , 1]; B <- pin$bruto[, , 30]
  # para a demonstração, o exemplar 30 é girado 40° e ampliado 1,5 vez
  a <- 40 * pi / 180; B <- 1.5 * B %*% matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2) + 1
  cen <- function(M) sweep(M, 2, colMeans(M)); tam <- function(M) sqrt(sum(cen(M)^2))
  if (etapa >= 2) { A <- cen(A); B <- cen(B) }
  if (etapa >= 3) { A <- A / tam(A); B <- B / tam(B) }
  if (etapa >= 4) { s <- svd(t(A) %*% B); R <- s$v %*% t(s$u); B <- B %*% R }
  tit <- c("coordenadas brutas", "1. translação: centrar", "2. escala: tamanho do centroide = 1", "3. rotação: mínimos quadrados")[etapa]
  ggplot() + desenha_marcos(A, pin$links, cores[["accent"]]) + desenha_marcos(B, pin$links, cores[["ferrugem"]]) +
    coord_equal() + labs(title = tit, x = NULL, y = NULL) + tema_aula(15, conceitual = TRUE) +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = .5), axis.text = element_blank())
}

#' Os 40 exemplares antes e depois da GPA, com o consenso
fig_gpa <- function(pin = pleth_ind()) {
  bruto <- purrr::map_dfr(1:40, \(i) tibble(x = pin$bruto[, 1, i], y = pin$bruto[, 2, i], ex = i))
  ali <- purrr::map_dfr(1:40, \(i) tibble(x = pin$gpa$coords[, 1, i], y = pin$gpa$coords[, 2, i], ex = i))
  cons <- pin$gpa$consensus
  g1 <- ggplot(bruto, aes(x, y)) + geom_point(colour = cores[["neutro"]], size = 1.5, alpha = .7) +
    labs(title = "antes", x = NULL, y = NULL) + tema_aula(15, conceitual = TRUE)
  g2 <- ggplot(ali, aes(x, y)) + geom_point(colour = cores[["neutro"]], size = 1.5, alpha = .6) +
    desenha_marcos(cons, pin$links, cores[["ferrugem"]]) +
    coord_equal() + labs(title = "depois da GPA (consenso em ferrugem)", x = NULL, y = NULL) + tema_aula(15, conceitual = TRUE)
  tt <- theme(plot.title = element_text(size = 18, face = "bold", hjust = .5), axis.text = element_blank())
  (g1 + tt) + (g2 + tt) + plot_layout(widths = c(1, 1.6))
}

fig_morfoespaco <- function(pin = pleth_ind()) {
  pc <- geomorph::gm.prcomp(pin$gpa$coords)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  d <- tibble(PC1 = pc$x[, 1], PC2 = pc$x[, 2], especie = pin$especie, local = pin$local) |>
    mutate(grupo = paste(ifelse(especie == "Jord", "P. jordani", "P. teyahalee"),
                         ifelse(local == "Allo", "(alopatria)", "(simpatria)")))
  ggplot(d, aes(PC1, PC2, colour = grupo, shape = grupo)) + geom_point(size = 4) +
    scale_colour_manual(values = cores_linhagem[1:4], name = NULL) + scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
    labs(x = sprintf("PC1 (%s%%)", num(100 * ve[1], 0)), y = sprintf("PC2 (%s%%)", num(100 * ve[2], 0))) +
    tema_aula(16) + theme(legend.position = "right")
}

fig_etapas_gm <- function() {
  passos <- c("digitalizar os marcos\n(e semimarcos)", "superposição de\nProcrustes (GPA)",
              "estatística\n(PCA, MANOVA, PGLS...)", "visualização\n(grades de deformação)")
  cx <- tibble(x = 0:3 * 4, y = 0, w = 3.6, h = 1.3, rotulo = passos,
               fill = c(cores[["papel"]], cores[["ferrugem"]], cores[["papel"]], cores[["papel"]]),
               cor = c(cores[["texto"]], cores[["ferrugem"]], cores[["texto"]], cores[["texto"]]),
               txt = c(cores[["texto"]], "white", cores[["texto"]], cores[["texto"]]))
  st <- tibble(x = cx$x[-4] + 1.8, y = 0, xend = cx$x[-1] - 1.8, yend = 0)
  diagrama(cx, st, base = 4.5)
}

# ---- Bloco 2: sinal, taxa e PGLS ----------------------------------------------------

fig_kmult <- function(r) {
  d <- tibble(K = r$k_nula)
  ggplot(d, aes(K)) + geom_histogram(bins = 40, fill = "#cfc8b8", colour = "white") +
    geom_vline(xintercept = r$k, colour = cores[["ferrugem"]], linewidth = 1.4) +
    annotate("text", x = r$k, y = Inf, vjust = 1.5, hjust = -.08, size = 6, colour = cores[["ferrugem"]],
             label = paste0("K[mult] == \"", num(r$k, 2), "\""), parse = TRUE) +
    labs(x = "K multivariado com as pontas permutadas", y = "permutações") + tema_aula(16)
}

#' Filomorfoespaço da forma: PC1 × PC2 com a árvore (nós estimados por ML)
fig_filomorfo_forma <- function(ps = pleth_sp()) {
  tr <- ps$arvore
  pc <- geomorph::gm.prcomp(ps$gpa$coords)
  S <- pc$x[tr$tip.label, 1:2]
  an <- apply(S, 2, \(v) phytools::fastAnc(tr, v))
  todos <- rbind(S, an); n <- ape::Ntip(tr)
  rownames(todos) <- c(tr$tip.label, rownames(an))
  ed <- tibble(x = todos[tr$edge[, 1], 1], y = todos[tr$edge[, 1], 2], xend = todos[tr$edge[, 2], 1], yend = todos[tr$edge[, 2], 2])
  pts <- tibble(x = S[, 1], y = S[, 2], l = sub("P_", "P. ", tr$tip.label))
  ggplot() + geom_segment(data = ed, aes(x, y, xend = xend, yend = yend), colour = cores[["neutro"]], linewidth = .9) +
    geom_point(data = tibble(x = an[, 1], y = an[, 2]), aes(x, y), colour = cores[["ferrugem"]], size = 2.5) +
    geom_point(data = pts, aes(x, y), colour = cores[["accent"]], size = 4) +
    geom_text(data = pts, aes(x, y, label = l), vjust = -1, size = 4.6, fontface = "italic") +
    scale_y_continuous(expand = expansion(mult = .12)) +
    labs(x = "PC1 da forma", y = "PC2 da forma") + tema_aula(16)
}

# ---- Bloco 3: modelos -----------------------------------------------------------------

elipse <- function(S, n = 100, nivel = .95) {
  t <- seq(0, 2 * pi, length.out = n)
  r <- sqrt(stats::qchisq(nivel, 2))
  e <- eigen(S, symmetric = TRUE)
  pts <- cbind(cos(t), sin(t)) %*% diag(sqrt(pmax(e$values, 0)) * r) %*% t(e$vectors)
  tibble(x = pts[, 1], y = pts[, 2])
}

#' Uma elipse por ecomorfo, a partir da matriz de taxas do BMM
fig_elipses <- function(b = anolis_bmm()) {
  nomes <- c(CG = "crown-giant", GB = "grass-bush", TC = "trunk-crown", TG = "trunk-ground", Tr = "trunk", Tw = "twig")
  d <- purrr::map_dfr(dimnames(b$sigma)[[3]], \(k) elipse(b$sigma[, , k]) |> mutate(eco = nomes[[k]]))
  ggplot(d, aes(x, y, colour = eco)) + geom_path(linewidth = 1.3) +
    scale_colour_manual(values = cores_linhagem[1:6], name = NULL) +
    coord_equal() + labs(x = "taxa em log SVL", y = "taxa em log cauda") + tema_aula(16) +
    theme(legend.position = "right")
}

#' Esquema das restrições sobre as matrizes de taxas de dois regimes
fig_restricoes <- function() {
  R <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)
  m <- function(v, a) R(a) %*% diag(v) %*% t(R(a))
  casos <- list(
    "sem\nrestrição" = list(m(c(1, .25), .5), m(c(.6, .5), -.6)),
    "proporcional" = list(m(c(1, .25), .5), 0.45 * m(c(1, .25), .5)),
    "autovetores\ncompartilhados" = list(m(c(1, .25), .5), m(c(.45, .35), .5)),
    "variâncias\niguais" = list(m(c(1, .25), .5), m(c(1, .25), -.5)),
    "matrizes iguais" = list(m(c(1, .25), .5), m(c(1, .25), .5)))
  d <- purrr::imap_dfr(casos, \(par, nome) bind_rows(elipse(par[[1]]) |> mutate(reg = "regime 1"),
                                                     elipse(par[[2]]) |> mutate(reg = "regime 2")) |> mutate(caso = nome))
  d$caso <- factor(d$caso, names(casos))
  ggplot(d, aes(x, y, colour = reg, linetype = reg)) + geom_path(linewidth = 1.2) +
    facet_wrap(~caso, nrow = 1) + coord_equal() +
    scale_colour_manual(values = c(cores[["accent"]], cores[["ferrugem"]]), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    labs(x = "atributo 1", y = "atributo 2") + tema_aula(14, conceitual = TRUE) +
    theme(legend.position = "bottom", axis.text = element_blank(), strip.text = element_text(size = 13))
}

# ---- Bloco 5: PCA filogenética -----------------------------------------------------------

pca_filo_anolis <- function() {
  an <- anolis(); tr <- an$arvore
  e <- new.env(); utils::data("anole.data", package = "phytools", envir = e)
  X <- as.matrix(e$anole.data[tr$tip.label, ])   # já em log
  pp <- phytools::phyl.pca(tr, X, method = "BM", mode = "cov")
  pc <- stats::prcomp(X)
  ve_f <- diag(pp$Eval) / sum(diag(pp$Eval)); ve_c <- pc$sdev^2 / sum(pc$sdev^2)
  list(L_filo = pp$L[, 1:2], L_comum = cor(X, pc$x[, 1:2]), ve_f = ve_f, ve_c = ve_c)
}

fig_cargas <- function(p = pca_filo_anolis()) {
  nomes <- c(SVL = "SVL", HL = "cabeça", HLL = "membro post.", FLL = "membro ant.", LAM = "lamelas", TL = "cauda")
  d <- bind_rows(tibble(var = rownames(p$L_comum), carga = abs(p$L_comum[, 1]), tipo = "PCA comum"),
                 tibble(var = rownames(p$L_filo), carga = abs(p$L_filo[, 1]), tipo = "PCA filogenética")) |>
    mutate(var = factor(nomes[var], nomes))
  ggplot(d, aes(var, carga, fill = tipo)) + geom_col(position = "dodge", width = .7) +
    scale_fill_manual(values = c(cores[["neutro"]], cores[["accent"]]), name = NULL) +
    labs(x = NULL, y = "|correlação com o PC1|") + tema_aula(16) + theme(legend.position = "top")
}

# ---- PACA (Collyer & Adams 2021) --------------------------------------------------

#' PCA comum × PACA na forma da cabeça das 9 espécies de Plethodon
paca_pleth <- function(ps = pleth_sp()) {
  tr <- ps$arvore
  pc <- geomorph::gm.prcomp(ps$gpa$coords, phy = tr)
  pa <- geomorph::gm.prcomp(ps$gpa$coords, phy = tr, align.to.phy = TRUE)
  k1 <- function(x) phytools::phylosig(tr, setNames(x[, 1], rownames(x)), method = "K")
  invisible(capture.output(s_pa <- summary(pa)))   # o summary() imprime a tabela
  list(arvore = tr, pc = pc, pa = pa,
       var_pc1 = pc$d[1] / sum(pc$d),                      # variância no PC1 da PCA
       cov_pa1 = s_pa$PC.summary[2, 1],              # covariância com a filogenia no eixo 1
       var_pa1 = var(pa$x[, 1]) / sum(apply(pa$x, 2, var)), # variância no eixo 1 da PACA
       k_pc1 = k1(pc$x), k_pa1 = k1(pa$x))
}

#' Filomorfoespaços lado a lado: PCA comum e PACA
fig_paca <- function(p = paca_pleth()) {
  tr <- p$arvore; n <- ape::Ntip(tr)
  um <- function(obj, titulo, eixo) {
    xy <- rbind(obj$x[tr$tip.label, 1:2], obj$anc.x[, 1:2])
    ed <- tibble(x0 = xy[tr$edge[, 1], 1], y0 = xy[tr$edge[, 1], 2],
                 x1 = xy[tr$edge[, 2], 1], y1 = xy[tr$edge[, 2], 2])
    pt <- tibble(x = xy[1:n, 1], y = xy[1:n, 2], l = sub("^P_", "P. ", tr$tip.label))
    ggplot() +
      geom_segment(data = ed, aes(x0, y0, xend = x1, yend = y1), colour = "#b9b2a5", linewidth = .8) +
      geom_point(aes(xy[-(1:n), 1], xy[-(1:n), 2]), colour = cores[["neutro"]], size = 2) +
      geom_point(data = pt, aes(x, y), colour = cores[["accent"]], size = 3.2) +
      geom_text(data = pt, aes(x, y, label = l), size = 4, vjust = -.8, fontface = "italic", colour = cores[["texto"]]) +
      scale_x_continuous(expand = expansion(mult = .15)) + scale_y_continuous(expand = expansion(mult = .15)) +
      labs(title = titulo, x = eixo[1], y = eixo[2]) + tema_aula(15)
  }
  patchwork::wrap_plots(
    um(p$pc, "PCA comum: máxima variância", c("PC 1", "PC 2")),
    um(p$pa, "PACA: máxima covariância com a filogenia", c("PAC 1", "PAC 2")), ncol = 2) &
    theme(plot.background = element_rect(fill = cores[["papel"]], colour = NA))
}
