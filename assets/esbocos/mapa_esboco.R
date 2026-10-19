library(ggplot2); library(dplyr); library(tibble)
perg <- c("Há sinal?\n(padrão)", "Como evolui?\n(processo)", "Como era o ancestral?\n(história)", "Estão associados?\n(relação)")
cor_p <- c("#5b8e7d", "#c8553d", "#8e6c8a", "#3d6a99")
tipos <- c("Contínuo", "Discreto", "Multivariado")
celulas <- tribble(
  ~col, ~lin, ~txt, ~aula,
  1, 3, "K de Blomberg · λ de Pagel\nI de Moran · Mantel · PVR", "3",
  1, 2, "D de Fritz · estatística δ", "3",
  1, 1, "K multivariado", "8",
  2, 3, "BM · OU · EB · tendência\nmulti-taxas · multi-regime", "5",
  2, 2, "Mk: ER · SYM · ARD\ntaxas ocultas", "6",
  2, 1, "mvMORPH · taxas por módulo\nPhylogeneticEM", "8",
  3, 3, "máxima verossimilhança\ncontMap", "7",
  3, 2, "marginal · conjunta\nmapeamento estocástico", "7",
  3, 1, "filomorfoespaço", "7–8",
  4, 3, "PIC · PGLS · ANOVA filogenética\nmodelos mistos (brms)", "2, 4",
  4, 2, "Pagel (1994) · threshold\nregressão logística filogenética", "6",
  4, 1, "PGLS multivariado\nPCA filogenética", "8") |>
  mutate(x = col * 3, y = lin * 1.6 + 1.2)
destaque <- function(aula) grepl(paste0("\\b", aula, "\\b"), celulas$aula)
p <- ggplot() +
  annotate("text", x = 7.5, y = 7.75, label = "Que atributo eu tenho? Que pergunta eu faço?", size = 7, fontface = "bold", colour = "grey20") +
  geom_tile(data = tibble(x = 1:4 * 3, y = 7.05), aes(x, y), width = 2.85, height = 0.8, fill = cor_p) +
  geom_text(data = tibble(x = 1:4 * 3, y = 7.05, l = perg), aes(x, y, label = l), colour = "white", size = 4.3, fontface = "bold", lineheight = 0.9) +
  geom_tile(data = celulas, aes(x, y, fill = factor(col), alpha = destaque(2) | destaque(4)), width = 2.85, height = 1.45, colour = "white") +
  geom_text(data = celulas, aes(x, y, label = txt), size = 3.6, lineheight = 0.95, colour = "grey15") +
  geom_label(data = celulas, aes(x + 1.2, y + 0.55, label = paste("A", aula)), size = 3, linewidth = 0, fill = "white", colour = "grey30") +
  geom_text(data = tibble(y = 1:3 * 1.6 + 1.2, l = rev(tipos)), aes(x = 0.9, y, label = l), size = 4.5, fontface = "bold", hjust = 1, colour = "grey25") +
  annotate("tile", x = 7.5, y = 1.95, width = 11.85, height = 0.55, fill = "#d9a441") +
  annotate("text", x = 7.5, y = 1.95, label = "A 2 · O movimento browniano e a matriz C: a base de todas as colunas", size = 4.3, fontface = "bold", colour = "grey15") +
  annotate("tile", x = 7.5, y = 1.3, width = 11.85, height = 0.55, fill = "grey75") +
  annotate("text", x = 7.5, y = 1.3, label = "A 1 · A árvore e os dados: obter a árvore, descompasso árvore–atributos, árvore sem comprimento de ramo", size = 4, colour = "grey15") +
  scale_fill_manual(values = cor_p, guide = "none") + scale_alpha_manual(values = c(`TRUE` = 0.55, `FALSE` = 0.15), guide = "none") +
  coord_cartesian(xlim = c(-0.6, 13.5), ylim = c(0.9, 8.1), expand = FALSE) + theme_void() +
  theme(plot.background = element_rect(fill = "#faf7f0", colour = NA))
ggsave(file.path(Sys.getenv("OUT"), "mapa_esboco.png"), p, width = 13, height = 7.2, dpi = 110)
