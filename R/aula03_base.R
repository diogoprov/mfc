# Aula 3: funções comuns (dados, transformações de árvore, medidas de sinal
# calculadas à mão e conferidas com ape/phytools).
library(ape)
library(phytools)
library(dplyr)
library(tibble)
library(purrr)
library(patchwork)
source(here::here("R/tema.R"))
source(here::here("R/arvore.R"))
source(here::here("R/matriz.R"))
source(here::here("R/diagrama.R"))

#' Anolis do phytools: árvore com altura 1 e o SVL (log) na ordem das pontas
anolis <- function() {
  data(anoletree, package = "phytools", envir = environment())
  data(anole.data, package = "phytools", envir = environment())
  tr <- anoletree
  tr$edge.length <- tr$edge.length / max(phytools::nodeHeights(tr))
  tr <- ape::compute.brlen(tr, tr$edge.length)   # tira o mapa de estados (simmap)
  class(tr) <- "phylo"; tr$maps <- NULL; tr$mapped.edge <- NULL
  list(arvore = tr,
       svl = setNames(anole.data$SVL, rownames(anole.data))[tr$tip.label],
       ecomorfo = phytools::getStates(anoletree, "tips")[tr$tip.label])
}

#' Árvore de 5 espécies do exemplo do λ (lida da figura do original)
arvore_lambda5 <- function() ape::read.tree(text = "((sp1:2,sp2:2):8,(sp3:7,(sp4:1,sp5:1):6):3);")

#' Transformação λ de Pagel: ramos internos × λ, pontas alongadas para
#' manter a altura de cada ponta (λ = 0: estrela)
lambda_arvore <- function(tr, lam) {
  h <- ape::node.depth.edgelength(tr)
  n <- ape::Ntip(tr)
  e <- tr$edge
  nova <- tr$edge.length * lam
  term <- e[, 2] <= n
  nova[term] <- h[e[term, 2]] - lam * h[e[term, 1]]
  tr$edge.length <- nova
  tr
}

#' Árvore transformada pelo early burst (taxa r < 0): cada ramo de t1 a t2
#' vira (e^(r t2) − e^(r t1)) / r
eb_arvore <- function(tr, r) {
  h <- ape::node.depth.edgelength(tr)
  t1 <- h[tr$edge[, 1]]; t2 <- h[tr$edge[, 2]]
  tr$edge.length <- (exp(r * t2) - exp(r * t1)) / r
  tr
}

#' Covariância sob OU (raiz fixa), numa árvore ultramétrica
vcv_ou_arvore <- function(tr, alfa, sigma2 = 1) {
  S <- ape::vcv(tr); Tt <- max(diag(S))
  sigma2 / (2 * alfa) * exp(-2 * alfa * (Tt - S)) * (1 - exp(-2 * alfa * S))
}

#' Simula n réplicas de um atributo com covariância V (normal multivariada)
sim_mvn <- function(V, n = 1) {
  L <- t(chol(V))
  X <- L %*% matrix(rnorm(nrow(V) * n), nrow(V))
  rownames(X) <- rownames(V)
  X
}

#' K de Blomberg passo a passo (Blomberg et al. 2003)
k_passos <- function(tr, x) {
  C <- ape::vcv(tr)[names(x), names(x)]
  Ci <- solve(C); n <- length(x); um <- rep(1, n)
  a <- as.numeric((um %*% Ci %*% x) / (um %*% Ci %*% um))      # média filogenética
  mse0 <- sum((x - a)^2) / (n - 1)
  mse  <- as.numeric(t(x - a) %*% Ci %*% (x - a)) / (n - 1)
  esp  <- (sum(diag(C)) - n / sum(Ci)) / (n - 1)
  list(a = a, mse0 = mse0, mse = mse, obs = mse0 / mse, esp = esp, K = (mse0 / mse) / esp)
}

#' I de Moran com uma matriz de pesos W
#' (`padroniza = TRUE` divide cada linha de W pela sua soma, como o
#' `ape::Moran.I`)
moran_I <- function(x, W, padroniza = FALSE) {
  z <- x - mean(x); diag(W) <- 0
  if (padroniza) { rs <- rowSums(W); W[rs > 0, ] <- W[rs > 0, ] / rs[rs > 0] }
  length(x) / sum(W) * sum(W * outer(z, z)) / sum(z^2)
}

#' Correlograma filogenético: I de Moran por classe de distância patrística
#' (classes de largura igual, de 0 à maior distância)
correlograma <- function(tr, x, n_classes = 5) {
  D <- ape::cophenetic.phylo(tr)[names(x), names(x)]
  cortes <- seq(0, max(D), length.out = n_classes + 1)
  purrr::map_dfr(seq_len(n_classes), \(k) {
    W <- (D > cortes[k] - 1e-9 & D <= cortes[k + 1] + 1e-9) * 1; diag(W) <- 0
    if (sum(W) == 0) return(NULL)
    tibble(classe = k, d_media = mean(D[W == 1]), pares = sum(W) / 2, I = moran_I(x, W, padroniza = TRUE))
  })
}

#' PVR passo a passo: distâncias → dupla centralização → autovetores →
#' regressão do atributo nos k primeiros autovetores
pvr_passos <- function(tr, x, k = 2) {
  D <- ape::cophenetic.phylo(tr)[names(x), names(x)]
  n <- nrow(D); J <- diag(n) - 1 / n
  G <- -0.5 * J %*% D^2 %*% J
  ev <- eigen(G, symmetric = TRUE)
  pos <- ev$values > 1e-8
  V <- ev$vectors[, pos, drop = FALSE]; rownames(V) <- names(x)
  m <- lm(x ~ V[, seq_len(k), drop = FALSE])
  list(D = D, G = G, valores = ev$values[pos], V = V, modelo = m,
       P = fitted(m), S = resid(m), R2 = summary(m)$r.squared)
}

#' Curva PSR (Diniz-Filho et al. 2012): R² do PVR ao acrescentar autovetores
#' em sequência, contra a fração acumulada dos autovalores
psr_curva <- function(tr, x, max_k = NULL) {
  p <- pvr_passos(tr, x, k = 1)
  m <- if (is.null(max_k)) length(p$valores) - 1 else max_k
  tibble(k = seq_len(m),
         repr = cumsum(p$valores)[seq_len(m)] / sum(p$valores),
         R2 = purrr::map_dbl(seq_len(m), \(j) summary(lm(x ~ p$V[, seq_len(j)]))$r.squared))
}

# Formatação de números nos slides
num <- function(x, d = 2) format(round(x, d), nsmall = d, decimal.mark = ",")
