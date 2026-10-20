# Aula 4: regressão filogenética (PGLS), diagnóstico, ANOVA/ANCOVA
# filogenéticas, GLM filogenético, matrizes por regime e incerteza.
library(nlme)
source(here::here("R/aula03_base.R"))   # anolis(), lambda_arvore(), vcv_ou_arvore(), sim_mvn()

#' Dados dos Anolis para os modelos (log SVL, log membro posterior, lamelas)
anolis_dados <- function() {
  an <- anolis(); tr <- an$arvore
  data(anole.data, package = "phytools", envir = environment())
  d <- data.frame(sp = tr$tip.label,
                  svl = anole.data[tr$tip.label, "SVL"],
                  hll = anole.data[tr$tip.label, "HLL"],
                  lam = round(exp(anole.data[tr$tip.label, "LAM"])),
                  eco = factor(an$ecomorfo, levels = c("TG", "Tr", "TC", "GB", "CG", "Tw")))
  rownames(d) <- d$sp
  list(arvore = tr, dados = d)
}

#' Os modelos dos Anolis usados nos slides
modelos_anolis <- function() {
  a <- anolis_dados(); d <- a$dados; tr <- a$arvore
  ols  <- lm(hll ~ svl, d)
  pgls <- gls(hll ~ svl, d, correlation = corPagel(1, tr, form = ~sp), method = "ML")
  bm   <- gls(hll ~ svl, d, correlation = corBrownian(1, tr, form = ~sp))
  d$res <- resid(bm)
  list(arvore = tr, dados = d, ols = ols, pgls = pgls, bm = bm,
       lambda = as.numeric(pgls$modelStruct$corStruct))
}

#' Relembrando o OLS: dispersão com a reta e os resíduos
fig_ols <- function(seed = 41) {
  set.seed(seed)
  d <- tibble(x = runif(14, 1, 10)) |> mutate(y = 2 + .8 * x + rnorm(14, 0, 1.2))
  m <- lm(y ~ x, d); d$yh <- fitted(m)
  i <- which.max(abs(resid(m)))
  ggplot(d, aes(x, y)) +
    geom_segment(aes(xend = x, yend = yh), colour = cores[["ferrugem"]], linewidth = .7) +
    geom_abline(intercept = coef(m)[1], slope = coef(m)[2], colour = cores[["accent"]], linewidth = 1.2) +
    geom_point(size = 4, shape = 21, fill = cores[["verde"]], colour = "white", stroke = .8) +
    annotate("text", x = d$x[i] + .25, y = (d$y[i] + d$yh[i]) / 2, label = "e[i]", parse = TRUE,
             colour = cores[["ferrugem"]], size = 7, hjust = 0) +
    annotate("text", x = 1.3, y = max(d$y), hjust = 0, size = 7, colour = cores[["accent"]],
             label = "y[i] == b[0] + b[1] * x[i] + e[i]", parse = TRUE) +
    labs(x = "x", y = "y") + tema_aula(17, conceitual = TRUE)
}

#' GLS passo a passo na árvore-guia (B ~ A)
gls_passos <- function() {
  tr <- arvore_guia(); a <- atributos_guia()
  X <- cbind(1, a$A); y <- a$B
  C <- ape::vcv(tr)[a$sp, a$sp]; Ci <- solve(C)
  b_gls <- solve(t(X) %*% Ci %*% X, t(X) %*% Ci %*% y)
  b_ols <- solve(t(X) %*% X, t(X) %*% y)
  list(b_gls = as.numeric(b_gls), b_ols = as.numeric(b_ols))
}
fig_ols_gls <- function() {
  a <- atributos_guia() |> mutate(clado = ifelse(sp %in% paste0("sp", 1:4), "sp1–sp4", "sp5–sp7"))
  g <- gls_passos()
  ggplot(a, aes(A, B)) +
    geom_abline(intercept = g$b_ols[1], slope = g$b_ols[2], colour = cores[["ferrugem"]], linewidth = 1.2) +
    geom_abline(intercept = g$b_gls[1], slope = g$b_gls[2], colour = cores[["accent"]], linewidth = 1.2) +
    geom_point(aes(fill = clado), size = 5, shape = 21, colour = "white", stroke = .8) +
    geom_text(aes(label = sp), vjust = -1.1, size = 5) +
    scale_fill_manual(values = cores_linhagem[c(1, 3)], name = NULL) +
    scale_y_continuous(expand = expansion(mult = .15)) +
    labs(x = "atributo A", y = "atributo B") + tema_aula(17)
}

#' A inclinação depende do modelo evolutivo em C (ideia de Mazel et al. 2016)
fig_modelos_C <- function() {
  m <- modelos_anolis(); d <- m$dados; tr <- m$arvore
  fit <- function(V, nome) {
    Vi <- solve(V[d$sp, d$sp]); X <- cbind(1, d$svl)
    b <- solve(t(X) %*% Vi %*% X, t(X) %*% Vi %*% d$hll)
    tibble(modelo = nome, b0 = b[1], b1 = b[2])
  }
  r <- bind_rows(
    fit(diag(nrow(d)) |> `dimnames<-`(list(d$sp, d$sp)), "OLS (estrela)"),
    fit(ape::vcv(tr), "BM"),
    fit(ape::vcv(lambda_arvore(tr, .5)), "λ = 0,5"),
    fit(vcv_ou_arvore(tr, 5), "OU, α = 5")) |>
    mutate(modelo = factor(modelo, levels = modelo))
  ggplot(d, aes(svl, hll)) +
    geom_point(colour = cores[["neutro"]], size = 2.5) +
    geom_abline(data = r, aes(intercept = b0, slope = b1, colour = modelo), linewidth = 1.1) +
    scale_colour_manual(values = cores_linhagem[c(2, 1, 3, 4)],
                        labels = paste0(r$modelo, ": ", num(r$b1, 3)), name = "inclinação") +
    labs(x = "log SVL", y = "log membro posterior") + tema_aula(16) +
    theme(legend.position = "right")
}

#' Sinal nas variáveis × nos resíduos (ideia de Revell 2010)
sinal_residuos <- function(seed = 42) {
  an <- anolis(); tr <- an$arvore
  set.seed(seed)
  x <- phytools::fastBM(tr, sig2 = 1)[tr$tip.label]
  e <- setNames(rnorm(length(x), 0, .3), names(x))       # resíduos sem sinal
  y <- 1 + 1.2 * x + e
  r <- resid(lm(y ~ x))
  list(x = x, y = y, e = r, tr = tr,
       Kx = as.numeric(phytools::phylosig(tr, x)), Ky = as.numeric(phytools::phylosig(tr, y)),
       Ke = as.numeric(phytools::phylosig(tr, resid(lm(y ~ x)))),
       lambda_res = phytools::phylosig(tr, r, method = "lambda")$lambda)
}
fig_sinal_residuos <- function(s = sinal_residuos()) {
  d <- bind_rows(tibble(qual = "X", K = s$Kx), tibble(qual = "Y", K = s$Ky),
                 tibble(qual = "resíduos de Y ~ X", K = s$Ke)) |>
    mutate(qual = factor(qual, levels = qual))
  ggplot(d, aes(qual, K, fill = qual)) +
    geom_col(width = .6) +
    geom_hline(yintercept = 1, linetype = "22") +
    geom_text(aes(label = num(K)), vjust = -.5, size = 6.5, fontface = "bold") +
    scale_fill_manual(values = c(cores[["accent"]], cores[["accent"]], cores[["ferrugem"]]), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .15))) +
    labs(x = NULL, y = "K de Blomberg") + tema_aula(17)
}

#' Anolis: OLS × PGLS com λ
fig_anolis_pgls <- function(m = modelos_anolis()) {
  d <- m$dados
  ggplot(d, aes(svl, hll)) +
    geom_abline(intercept = coef(m$ols)[1], slope = coef(m$ols)[2], colour = cores[["ferrugem"]], linewidth = 1.1, linetype = "22") +
    geom_abline(intercept = coef(m$pgls)[1], slope = coef(m$pgls)[2], colour = cores[["accent"]], linewidth = 1.3) +
    geom_point(aes(fill = eco), size = 3.5, shape = 21, colour = "white", stroke = .5) +
    scale_fill_manual(values = cores_linhagem[1:6], name = "ecomorfo") +
    labs(x = "log SVL", y = "log membro posterior") + tema_aula(16) +
    theme(legend.position = "right")
}

#' Diagnóstico do PGLS com resíduos normalizados
fig_diagnostico <- function(m = modelos_anolis()) {
  d <- tibble(aj = fitted(m$pgls), rn = resid(m$pgls, type = "normalized"), svl = m$dados$svl)
  base <- tema_aula(14)
  g1 <- ggplot(d, aes(aj, rn)) + geom_hline(yintercept = 0, colour = cores[["neutro"]]) +
    geom_point(colour = cores[["accent"]], size = 2.5) + labs(x = "valores ajustados", y = "resíduos normalizados") + base
  g2 <- ggplot(d, aes(sample = rn)) + stat_qq(colour = cores[["accent"]], size = 2.5) + stat_qq_line(colour = cores[["ferrugem"]]) +
    labs(x = "quantis teóricos", y = "quantis dos resíduos") + base
  g3 <- ggplot(d, aes(rn)) + geom_histogram(bins = 15, fill = cores[["accent"]], colour = cores[["papel"]]) +
    labs(x = "resíduos normalizados", y = "frequência") + base
  g4 <- ggplot(d, aes(svl, rn)) + geom_hline(yintercept = 0, colour = cores[["neutro"]]) +
    geom_point(colour = cores[["accent"]], size = 2.5) + labs(x = "log SVL", y = "resíduos normalizados") + base
  (g1 | g2 | g3 | g4)
}

#' Membro posterior relativo (resíduo do PGLS sob BM) por ecomorfo
fig_ecomorfos <- function(m = modelos_anolis()) {
  d <- m$dados
  ggplot(d, aes(eco, res, colour = eco)) +
    geom_hline(yintercept = 0, colour = cores[["neutro"]]) +
    geom_jitter(width = .15, size = 3, alpha = .8) +
    stat_summary(fun = mean, geom = "crossbar", width = .5, colour = cores[["texto"]]) +
    scale_colour_manual(values = cores_linhagem[1:6], guide = "none") +
    scale_x_discrete(labels = c(TG = "TG\ntrunk-ground", Tr = "Tr\ntrunk", TC = "TC\ntrunk-crown",
                                GB = "GB\ngrass-bush", CG = "CG\ncrown-giant", Tw = "Tw\ntwig")) +
    labs(x = NULL, y = "membro posterior relativo") + tema_aula(15)
}

#' ANOVA filogenética por simulação (Garland et al. 1993)
anova_filo <- function(m = modelos_anolis(), nsim = 999, seed = 43) {
  d <- m$dados
  set.seed(seed)
  r <- phytools::phylANOVA(m$arvore, setNames(d$eco, d$sp), setNames(d$res, d$sp), nsim = nsim)
  # F nulo, simulado sob BM com a taxa estimada (como no phylANOVA)
  s2 <- mean(ape::pic(setNames(d$res, d$sp), m$arvore)^2)
  Fnulo <- replicate(nsim, {
    y <- phytools::fastBM(m$arvore, sig2 = s2)[d$sp]
    summary(aov(y ~ d$eco))[[1]][1, "F value"]
  })
  list(F = r$F, p = r$Pf, Fnulo = Fnulo, Fcomum = summary(aov(d$res ~ d$eco))[[1]][1, "Pr(>F)"])
}
fig_anova_filo <- function(a = anova_filo()) {
  ggplot(tibble(F = a$Fnulo), aes(F)) +
    geom_histogram(bins = 50, fill = cores[["neutro"]], colour = cores[["papel"]]) +
    geom_vline(xintercept = a$F, colour = cores[["ferrugem"]], linewidth = 1.3) +
    geom_vline(xintercept = qf(.95, 5, 76), colour = cores[["texto"]], linetype = "22") +
    annotate("text", x = a$F, y = Inf, vjust = 1.5, hjust = 1.05, size = 5.5, colour = cores[["ferrugem"]],
             label = paste0("F observado = ", num(a$F, 1))) +
    annotate("text", x = qf(.95, 5, 76), y = Inf, vjust = 3.5, hjust = -.05, size = 5,
             label = "F crítico da ANOVA comum") +
    labs(x = "F simulado sob BM (999 simulações)", y = "frequência") + tema_aula(16)
}

#' GLM filogenético: logística e Poisson (phylolm::phyloglm)
glm_filo <- function(m = modelos_anolis()) {
  d <- m$dados; d$tg <- as.integer(d$eco == "TG"); d$res_z <- as.numeric(scale(d$res))
  lg <- phylolm::phyloglm(tg ~ res_z, data = d, phy = m$arvore, method = "logistic_MPLE")
  po <- phylolm::phyloglm(lam ~ svl, data = d, phy = m$arvore, method = "poisson_GEE")
  list(d = d, lg = lg, po = po)
}
fig_glm_filo <- function(g = glm_filo()) {
  d <- g$d
  cl <- coef(g$lg); cp <- coef(g$po)
  xs <- seq(min(d$res_z), max(d$res_z), length.out = 100)
  sv <- seq(min(d$svl), max(d$svl), length.out = 100)
  g1 <- ggplot(d, aes(res_z, tg)) +
    geom_point(colour = cores[["accent"]], size = 3, alpha = .7, position = position_jitter(height = .03, seed = 1)) +
    geom_line(data = tibble(res_z = xs, tg = plogis(cl[1] + cl[2] * xs)), colour = cores[["ferrugem"]], linewidth = 1.3) +
    scale_y_continuous(breaks = c(0, 1), labels = c("outro", "trunk-ground")) +
    labs(title = "Logística", x = "membro posterior relativo (desvios padrão)", y = NULL) + tema_aula(15)
  g2 <- ggplot(d, aes(svl, lam)) +
    geom_point(colour = cores[["accent"]], size = 3, alpha = .7) +
    geom_line(data = tibble(svl = sv, lam = exp(cp[1] + cp[2] * sv)), colour = cores[["ferrugem"]], linewidth = 1.3) +
    labs(title = "Poisson", x = "log SVL", y = "número de lamelas") + tema_aula(15)
  g1 + g2
}

#' Centrarchidae: árvore com os regimes e dispersão por regime
fig_sunfish <- function() {
  data(sunfish.tree, package = "phytools", envir = environment())
  data(sunfish.data, package = "phytools", envir = environment())
  tr <- sunfish.tree
  reg <- sapply(tr$maps, \(mp) names(mp)[which.max(mp)])          # regime predominante de cada ramo
  cor_r <- ifelse(reg == "pisc", cores[["ferrugem"]], cores[["accent"]])
  tr2 <- tr; class(tr2) <- "phylo"; tr2$maps <- NULL; tr2$mapped.edge <- NULL
  tr2$tip.label <- gsub("_", " ", tr2$tip.label)
  arv <- plot_arvore(tr2, cor_ramo = cor_r, base = 3.8, folga = .6)
  d <- sunfish.data |> tibble::rownames_to_column("sp") |>
    mutate(regime = ifelse(feeding.mode == "pisc", "piscívoros", "não piscívoros"))
  disp <- ggplot(d, aes(buccal.length, gape.width, fill = regime)) +
    geom_point(size = 4, shape = 21, colour = "white", stroke = .6) +
    scale_fill_manual(values = c(cores[["accent"]], cores[["ferrugem"]]), name = NULL) +
    labs(x = "comprimento bucal (corrigido pelo tamanho)", y = "abertura da boca (corrigida)") + tema_aula(15)
  arv + disp + plot_layout(widths = c(1.1, 1))
}
evol_vcv_sunfish <- function() {
  data(sunfish.tree, package = "phytools", envir = environment())
  data(sunfish.data, package = "phytools", envir = environment())
  phytools::evol.vcv(sunfish.tree, as.matrix(sunfish.data[, c("gape.width", "buccal.length")]))
}

#' Incerteza filogenética: a inclinação do PGLS em 100 árvores com
#' comprimentos de ramo perturbados (ilustração)
incerteza <- function(n = 100, sd_log = .4, seed = 44) {
  m <- modelos_anolis(); d <- m$dados; tr <- m$arvore
  set.seed(seed)
  b <- replicate(n, {
    t2 <- tr; t2$edge.length <- tr$edge.length * exp(rnorm(length(tr$edge.length), 0, sd_log))
    C <- ape::vcv(t2)[d$sp, d$sp]; Ci <- solve(C); X <- cbind(1, d$svl)
    solve(t(X) %*% Ci %*% X, t(X) %*% Ci %*% d$hll)[2]
  })
  list(b = b, b_orig = unname(coef(m$bm)[2]))
}
fig_incerteza <- function(i = incerteza()) {
  ggplot(tibble(b = i$b), aes(b)) +
    geom_histogram(bins = 25, fill = cores[["accent"]], colour = cores[["papel"]]) +
    geom_vline(xintercept = i$b_orig, colour = cores[["ferrugem"]], linewidth = 1.3) +
    annotate("text", x = i$b_orig, y = Inf, vjust = 1.5, hjust = -.05, size = 5.5, colour = cores[["ferrugem"]],
             label = "árvore original") +
    labs(x = "inclinação de membro posterior × SVL (PGLS, BM)", y = "árvores") + tema_aula(16)
}
