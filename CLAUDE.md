# Métodos Filogenéticos Comparativos: reformulação em Quarto/revealjs

Repositório da disciplina de pós-graduação **Métodos Filogenéticos
Comparativos** (201521439, 45 h), Mestrado em Biologia Animal, UFMS. Site
publicado em <https://provetelab.org/mfc/>. Responda em português.

O desenho da reformulação (diagnóstico, nova sequência de aulas, fio
condutor, mapeamento original → nova, decisões pendentes) está em
`PLANO-conversao-qmd.md`. Este arquivo guarda as regras e a situação atual.
As regras herdadas valem quando este arquivo não disser outra coisa:
`~/GitHub/multidimensionais/CLAUDE.md` e `WORKFLOW-conversao-qmd.md` (a
referência mais próxima), depois `~/GitHub/bioestatistica/`.

## Objetivo

**Reformular**, não só converter: reorganizar as 6 aulas de 2020 (472
slides em `.pptx`) em 8 aulas em torno do mapa *tipo de atributo ×
pergunta* (do Coggle do Diogo), em slides Quarto/revealjs
(`slides/aNN.qmd`), com figuras geradas em R a partir de dados de pacotes,
simulações ou da árvore-guia.

## Situação

| Nova aula | Título | Origem principal | Inventário | Deck |
|---:|---|---|---|---|
| 1 | A árvore e os dados | A1 8–158, A5 1–7 | **`inventarios/a01.md`** | **`slides/a01.qmd` ✔** (revisada e publicada) |
| 2 | O movimento browniano e a matriz C | A2 1–49 | **`inventarios/a02.md`** | **`slides/a02.qmd` ✔** (aula-modelo; revisada e publicada) |
| 3 | Sinal filogenético | A3 24–61, A5 7–10 | **`inventarios/a03.md`** | **`slides/a03.qmd` ✔** (revisada e publicada) |
| 4 | Regressão filogenética (PGLS) | A2 126–141, A4 50–58 | **`inventarios/a04.md`** | **`slides/a04.qmd` ✔** (revisada e publicada) |
| 5 | Modelos de evolução contínua | A2 50–88, A4 18–49 | **`inventarios/a05.md`** | **`slides/a05.qmd` ✔** (revisada e publicada) |
| 6 | Atributos discretos | A2 89–126, A4 1–17, A1 65–77 | **`inventarios/a06.md`** | **`slides/a06.qmd` ✔** (revisada e publicada) |
| 7 | Estados ancestrais e visualização | A3 1–23, A5 11 | **`inventarios/a07.md`** | **`slides/a07.qmd` ✔** (revisada e publicada) |
| 8 | Atributos multivariados | A6 1–39 | **`inventarios/a08.md`** | **`slides/a08.qmd` ✔** (revisada e publicada) |

**Fase 0 (plano) feita em 19/out/2026.** Aprovados pelo Diogo: as 8 aulas,
o mapa dos métodos e a Aula 2 como aula-modelo; as decisões 4–11 do plano
aguardam resposta. Feito: plano, este
arquivo, inventário da aula-modelo e o esboço do mapa
(`assets/esbocos/`). Próximo, depois das respostas: fase 1
(infraestrutura) e a aula-modelo. **Não escrever o deck antes das respostas
às decisões 4–11.**

**Site reorganizado (20/out/2026)**: oito pastas `aula-0N-*` (`arvore-dados`,
`bm-matriz-c`, `sinal`, `pgls`, `modelos`, `discretos`, `ancestrais`,
`multivariados`), cada uma com `index.qmd` (slides, em 3 minutos, conteúdo,
roteiro, dados, desafio, leituras). Os walkthroughs de 2020 ficaram inteiros
na pasta da aula com mais conteúdo deles e têm `aliases` para as URLs
antigas; as páginas antigas das 6 aulas redirecionam para as novas. Páginas
gerais: `comece-aqui.qmd`, `programa.qmd` (ementa literal + distribuição nas
aulas), `praticas.qmd`, `glossario.qmd`; cronograma em cartões
(`_cronograma.ejs`, campos `gancho` e `familias`: base, sinal, modelos,
ancestral, relacao). Os dados das práticas mudaram de pasta: os caminhos em
`R/` já apontam para as novas.

**Roteiros executados (20/out/2026)**: um `roteiro.qmd` por pasta de aula,
com `eval: true` e freeze, alinhado aos slides e ao desafio (as dicas do
desafio apontam as seções). Os dados são lidos relativos à pasta: os
arquivos de que as Aulas 6 e 7 precisam foram copiados para as pastas
delas. `praticas.qmd` lista os oito roteiros e, à parte, os walkthroughs de
2020. Os roteiros 1 e 3 usam `rotl`, `adephylo` e `phylobase` (fora das
figuras). Desafios e roteiros conferidos e aprovados pelo Diogo
(20/out/2026).

**Formato**: a disciplina é ofertada condensada em uma semana; cada aula
ocupa um período (manhã ou tarde) com teoria e prática, 8 aulas em 10
períodos. Decks de ~50–65 slides.

## Como trabalhar com o Diogo

- Português na conversa; slides em português, termos técnicos podem ficar
  em inglês. Texto destinado a publicação (artigos) em inglês.
- **Commit e push só com aprovação.** Na multidimensionais, depois da
  aula-modelo, passou a "converta sem esperar, eu reviso no site"; aqui,
  esperar ele dizer que vale o mesmo.
- **Inventário slide a slide antes de cada deck** (`inventarios/aNN.md`),
  com a coluna do que entra de novo; conferir o deck pronto contra ele.
- A revisão chega slide a slide, por URL; trate na hora e volte ao que
  estava fazendo.
- Aponte os erros do original e relate as correções.
- Aplique a skill **karpathy-guidelines**.
- Antes de reescrever muito, cortar conteúdo ou mudar escopo, combine.

## Regras dos slides (herdadas, por ordem de importância)

1. **Preserve os esquemas do original e alterne visual e texto.** Nada de
   sequências só de texto. Sequências montadas aos poucos no original (a
   **C** em 9 passos, as linhagens de BM, a **Q** do Mk) viram vários slides
   com a mesma figura crescendo.
2. **Fidelidade ao original** (Diogo, 19/out/2026): o deck reformulado deve
   ficar o mais fiel possível ao conteúdo do original, corrigindo erros,
   duplicatas e redundâncias. Antes de publicar uma aula, rode
   `python3 -I ferramentas/fidelidade.py <texto-do-pptx> <slide-ini> <slide-fim> slides/aNN.qmd`
   (o texto vem de `ferramentas/extrai.py <pptx> <saida.txt>`): ele lista os
   slides do original cujas palavras-chave não aparecem no deck. Cuidado com
   leituras erradas dos métodos (ex.: os círculos de Pavoine são contribuições
   à diversidade, não estados ancestrais).
3. **Não corte conteúdo.** Todo bloco do original tem destino na tabela da
   seção 4 do plano. Encolher só o que o Diogo aprovar (plano, seção 6).
4. **Passo a passo das contas**, um passo por fragmento (montar **C**, os
   contrastes, *P(t) = e^Qt*).
5. **Texto grande e espaço aproveitado** (pedido do Diogo na revisão da
   Aula 2): corpo em 32 px, nada abaixo de ~26 px fora das `.ref`. Slide de
   texto com espaço sobrando leva `{.grande}` (1,15em) ou `{.muito-grande}`
   (1,3em); `{.medio}` só quando não couber. Figuras na largura toda
   (`out.width: "100%"`). Meça a altura usada: mire em 600–745 px.
6. **Referências com link** para o DOI, conferidas no Crossref.
7. **Reaproveite** cursos irmãos: Procrustes e PCA da multidimensionais
   (Aulas 6 e 4) na Aula 8; OLS e pressupostos da univariada/Bioestatística
   na Aula 4.
8. **Pergunte antes de adicionar um pacote.**

## Decisões da disciplina

Fio condutor (mapa) e estrutura de 8 aulas **aprovados**; o resto são
propostas **a confirmar** com o Diogo (plano, seção 12). Atualizar aqui
quando ele decidir.

- **Fonte**: os `.pptx` de out/2021 em
  `~/Library/CloudStorage/OneDrive-Pessoal/Aulas/2017 - PCM UFSCAR/`
  (`Aula 1.pptx` a `Aula 5.pptx`; a Aula 5 é de 2016) e
  `2017 - PCM UFMS/Aula sexta.pptx` (a Aula 6). As pastas `2016 - Rio
  Preto PCM` e `2016 - PCM unifesp` têm versões mais antigas. Os PDFs de
  `slides/` não batem com os `.pptx` (salvo `aula-01.pdf` e `aula-05.pdf`):
  exportar o `.pptx` para PDF pelo PowerPoint via AppleScript
  (`save active presentation in ... as save as PDF`) para as folhas de
  contato.
- **Fio condutor**: o mapa dos métodos (`mapa_plot(destaque = ...)`, em
  `R/mapa.R`, a partir de `assets/esbocos/mapa_esboco.R`). Colunas = as
  perguntas (*Há sinal?*, *Como evolui?*, *Como era o ancestral?*, *Estão
  associados?*); linhas = tipo de atributo (contínuo, discreto,
  multivariado); faixas de base = Aula 1 (a árvore e os dados) e Aula 2 (BM
  e **C**).
- **Selos**: `[Há sinal?]{.familia .sinal}`; classes `.base`, `.sinal`,
  `.modelos`, `.ancestral`, `.relacao`. Cores em `cores_familia`
  (`R/tema.R`) e `$familia-*` (`slides/theme.scss`): mude nos dois lugares.
- **Árvore-guia** (aprovada pelo Diogo): a de Freckleton et al. (2002), 7 espécies,
  `(((sp1:0.18,sp2:0.18):1.06,(sp3:0.02,sp4:0.02):1.22):0.20,(sp5:1.37,(sp6:0.23,sp7:0.23):1.14):0.07);`
  com os atributos A e B do mesmo slide. Sempre desenhada pela mesma função.
- **Dados que atravessam o curso** (aprovado): `phytools::anoletree` e
  `anole.data`, no curso todo.
  Exemplo brasileiro: cantos de anuros da Bocaina (`canto.csv`,
  `filogenia_bocaina.txt`), de Garey et al. (2018)
  `10.1093/biolinnean/bly129`, uso autorizado pelo Diogo.
- **Livro-texto**: Revell & Harmon (2022), *Phylogenetic Comparative
  Methods in R* (Princeton; sem DOI, link da editora). Toda aula cita os
  capítulos/seções nas leituras; mapa capítulos → aulas no plano, seção 7.3.
- **Leituras prévias**: Rezende & Diniz-Filho (2012), O'Meara (2012) e
  Cornwallis & Griffin (2024) (no lugar de Garland, Bennett & Rezende 2005,
  troca do Diogo em 20/out/2026), lidas antes da disciplina; a Aula 1 as
  supõe lidas (abertura: plano, seção 3.7).
- **Aula 1**: bases, megafilogenias, pipelines e softwares resumidos em
  tabelas com links; NGS em um slide (aprovado).
- **Inferência filogenética na Aula 1**: resumida (aprovado). Só o que o
  aluno precisa: de onde vem a árvore, como é construída, como ler e
  interpretar, e as limitações. Modelos de DNA → Aula 6 (ponte para o Mk).
- **Pacotes das figuras**: `ggplot2`, `dplyr`, `tidyr`,
  `purrr`, `tibble`, `patchwork`, `here`, `knitr`, `ape`, `phytools`,
  `geiger`, `nlme` (**aprovados**); `caper` (aprovado para o D de Fritz &
  Purvis, Aula 3); `phylolm` (Aula 4); `mvMORPH` (OU multi-regime, Aula 5;
  previsto no plano para a Aula 8); `corHMM` (taxas ocultas e
  `fitCorrelationTest`, Aula 6) e `rr2` (R² de PGLS, Aula 4), ambos sugeridos
  pelo Diogo; `geomorph` (Aula 8, aprovado em 20/out/2026). PhylogeneticEM
  só citado, sem instalar. Sem `ggtree` (árvores por função
  própria) e sem `brms`. `rphylopic` só para baixar silhuetas, fora do render.
- **Não entram** (como na multidimensionais): Kahoot, webR, "Teste-se",
  atividades extras nos slides.
- **Atividade de cada aula = o walkthrough** (prática no mesmo período);
  além dele, **um desafio por aula**: um problema para o aluno resolver
  **sozinho, com código** (dados de pacotes; no fim da página da aula). As
  respostas não vão para o repositório antes do prazo.
- **Avaliação**: documento **Quarto** (pergunta, análise, resultados,
  discussão curta), descrito na página do programa.
- **Memes**: pinturas em domínio público, de preferência brasileiras
  (Almeida Júnior, Belmiro de Almeida, Amoedo, Pedro Américo, Visconti, os
  Timótheo da Costa, Ismael Nery; Portinari, Tarsila, Di Cavalcanti e Anita
  Malfatti **ainda não**). Nada de símbolos patrióticos ou militares dos
  EUA ou de outras potências. Todo meme com título `##` próprio. Lista dos
  já usados na multidimensionais em `~/GitHub/multidimensionais/CLAUDE.md`:
  não repetir sem motivo.
- **Regra geral: use o máximo possível de imagens CC BY e CC0** (Diogo,
  19/out/2026): figura de artigo ou do Commons com essas licenças vem antes
  do redesenho. CC BY-NC só quando não houver alternativa.
- **Aula 6, correlação de discretos** (Diogo, 19/out/2026): o modelo de
  Pagel (1994) fica como base conceitual; depois o problema da
  pseudorreplicação (Maddison & FitzJohn 2015; Uyeda et al. 2018) e, como
  ferramenta principal, o teste com taxas ocultas do `corHMM` (Boyko &
  Beaulieu 2023). BayesTraits/`btw` só como alternativa citada.
- **Aula 4**: além do PGLS, só `phylolm::phyloglm` (logística e Poisson).
  Nada de `brms` nem de modelos mistos (exigem background demais).
- **Abertura da Aula 1** (aprovada): 3–4 slides que retomam as leituras
  prévias (perguntas, mapa como síntese, não independência, linha do tempo).
- **Créditos de slides de colegas**: o original não marca tudo. Já
  identificados pelo Diogo: o esquema de Hansen & Martins e o caminho inverso
  (J. A. F. Diniz-Filho, Aula 2), os dados como processo estocástico (José
  Hidasi-Neto, Aula 2). Os slides de **fundo azul** do original (ex.: o formato
  Newick, Aula 1 original, slides 42–46) são de J. A. F. Diniz-Filho:
  refazer/adaptar, com crédito. `.ref`: "Ideia de um slide de ...". Na dúvida,
  pergunte de quem é um slide sem crédito.
- **Figuras de terceiros**: livro de Harmon (2019) é CC BY-SA 4.0 →
  redesenhar. Slides de J. A. F. Diniz-Filho (Aula 3 original) → aproveitar a
  ideia e **redesenhar pesadamente** (aprovado); `.ref`: "Ideia de um slide
  de J. A. F. Diniz-Filho". Licença conferida no Crossref/Unpaywall
  antes de usar qualquer figura de artigo; registro em
  `assets/terceiros/FONTES.md`.

## Arquivos de apoio

- `R/tema.R`: `cores`, `cores_familia` (perguntas do mapa), `cores_linhagem`,
  `tema_aula()`, `separa_rotulos()`.
- `R/mapa.R`: `mapa_plot(aula = ...)`; células em `mapa_celulas`.
- `R/arvore.R`: `arvore_guia()`, `atributos_guia()`, `coord_arvore()`,
  `caminho_raiz()`, `caminho_comum()`, `plot_arvore()`, `sim_bm_arvore()`,
  `fenograma()`.
- `R/matriz.R`: `plot_matriz()` (destaque, células ocultas, colchetes).
- `R/bm.R`: simulações no tempo (`sim_bm_tempo`, `sim_ou_tempo`,
  `sim_eb_tempo`, `sim_trend_tempo`, `sim_linhagens`, `sim_dados`).
- `R/diagrama.R`: `diagrama(caixas, setas)` (da multidimensionais).
- `R/aula03_base.R`: `anolis()`, `lambda_arvore()`, `eb_arvore()`,
  `vcv_ou_arvore()`, `sim_mvn()`, `k_passos()`, `moran_I()`, `correlograma()`,
  `pvr_passos()`, `psr_curva()` (todas conferidas com ape/phytools);
  `R/aula03_*.R`: figuras da Aula 3.
- `R/aula05_modelos.R`: verossimilhança do BM, simulações de OU/EB/tendência,
  `ajustes_anolis()` (`fitContinuous`, com semente), BBM/FPK/Lévy,
  `ajustes_regimes()` (`mvMORPH`, OUM com os ecomorfos), `brownie_sunfish()`.
- `R/aula06_discretos.R`: `squamata()` (Brandley et al. 2008, binário e 4
  estados), `ajustes_squamata()` (`fitDiscrete`; reproduz o log L do
  original), `pagel_pratica()`, `correlacao_corhmm()` (~3 min) e as figuras
  da Aula 6.
- `R/aula08_multivariados.R`: `pleth_ind()`/`pleth_sp()` (`geomorph`),
  `multivar_resultados()` (K mult., PGLS, `mvgls`, modularidade, taxas por
  módulo; sementes 801–805), `anolis_bmm()`, `pca_filo_anolis()` e as figuras
  (Procrustes passo a passo, GPA, elipses, restrições). Não use `pi` como
  nome de variável: mascara a constante usada em `ret_arred()`.
- `R/aula01_arvore.R`: `dois_clados()` (X e Y independentes que parecem
  associados), `suporte_woodmouse()` (NJ + bootstrap reais), `formas()`
  (Colless e γ) e as figuras da Aula 1 (Newick em 5 passos, partes da árvore,
  enraizamento, tipos de árvore, fluxos).
- `R/aula07_ancestrais.R`: `guia_ancestrais()` (ML e parcimônia na
  árvore-guia), `enguias()`/`ajustes_enguias()` (Collar et al. 2014; marginal
  × conjunta com `phytools::ancr`, 100 mapas estocásticos),
  `mapas_squamata()`, `tendencia_fosseis()` e as figuras (algumas em gráficos
  base do `phytools`: `contMap`, `fancyTree`, `densityMap`,
  `phylomorphospace`).
- `R/aula04_pgls.R`: `anolis_dados()`, `modelos_anolis()`, `gls_passos()`,
  `anova_filo()`, `glm_filo()`, `evol_vcv_sunfish()`, `incerteza()` e as
  figuras da Aula 4.
- `R/aula02_*.R`: figuras da Aula 2 (`intro`, `bm`, `matriz`, `arvore`,
  `pic`). `pic_passos()` faz os contrastes nó a nó.
- `.claude/launch.json`: servidor `docs` (porta 4321) para o teste de
  transbordamento.

## Estrutura de um deck

Igual à multidimensionais: chunk `setup` com `source(here::here("R/..."))`
→ `## Hoje` → `## Onde estamos` (mapa destacado) → blocos
`# N. Título {.secao}` com o selo da pergunta → mapa de
fechamento → leituras com o aviso do Coggle → meme e prática.

- Não ponha `---` antes de um `# ` (slide em branco).
- `fig-alt` descritivo em português em toda figura.
- Números citados no texto saem de `` `r ...` ``, nunca digitados.
- Sementes fixas em toda simulação. Um arquivo `R/aulaNN_bloco.R` por
  bloco, só com funções `fig_*()`.
- Números com vírgula em LaTeX: `1{,}44`.

## Build e publicação

- O GitHub Actions roda só o Quarto: `freeze: auto` com `_freeze/`
  versionado (a configurar na fase 1). Mudou função em `R/`? Renderize o
  `.qmd` explicitamente. Acrescentou ou removeu chunk? Apague
  `_freeze/slides/aNN` antes.
- `README.md`, `CLAUDE.md`, `PLANO-conversao-qmd.md` e `inventarios/` ficam
  fora do render (lista `render` do `_quarto.yml`).
- Teste de transbordamento no navegador embutido, viewport 1250 × 750, em
  unidades de slide (ver o CLAUDE.md da multidimensionais).

## Armadilhas (previstas e encontradas)

- `geiger::fitContinuous` sorteia pontos de partida: fixe a semente, ou o log L
  muda entre renderizações. `mvMORPH::mvOU` não converge na árvore dos
  *Anolis* com altura 6: reescale para altura 1 (ramos e `maps`).

- Com `options(OutDec = ",")`, `paste0()` de números escreve vírgula: um
  Newick montado assim quebra (comprimentos negativos). Use `sprintf("%g")`.
- Rótulos `plotmath` com número decimal: ponha o número entre aspas
  (`sigma^2 == "0,01"`), senão a vírgula vira erro de sintaxe.
- Figuras no revealjs: `out.width: "100%"` no YAML do deck; não combine com
  `out-height` (distorce). Ajuste a altura pelo `fig-height`.
- `R/arvore.R`: as pontas são desenhadas de cima para baixo (sp1 no alto,
  como no original); os nós internos seguem a numeração do `ape` (8 = raiz
  da árvore-guia; 9 = sp1–sp4; 10 = sp1+sp2; 11 = sp3+sp4; 12 = sp5–sp7;
  13 = sp6+sp7).

- `phytools::anole.data` **já está em log natural** (SVL 3,46–5,11 = 32–166
  mm): não aplique `log()` de novo. `geiger::treedata()` não põe as linhas na
  ordem das pontas sem `sort = TRUE`.
- Árvore não ultramétrica ou sem comprimentos de ramo; nomes das pontas que
  não batem com os dados (`geiger::name.check`); ordem das linhas em
  `nlme::gls` com `corBrownian`/`corPagel` (use `form = ~sp`).
- λ e K com poucas espécies; α e σ² do OU correlacionados (Ho & Ané 2014).
- O original usa θ tanto para o valor inicial do BM quanto para o ótimo do
  OU: θ só para o ótimo; valor inicial = *z*₀.
- `ape::plot.phylo(plot = FALSE)` abre um dispositivo e gera imagem em
  branco: envolva em `pdf(NULL)`/`dev.off()`.
- O `convert` do ImageMagick instalado é binário Intel e não roda: folhas
  de contato com `gs -sDEVICE=png16m` e o pacote `magick` do R
  (`image_convert(m, "png")` antes de `image_write` do montage).
- O Coggle não exporta a hierarquia sem login; o texto dos nós sai da
  página, e a estrutura, pelas cores dos ramos.
- Python local sem certificados SSL: use `curl` para o Crossref.

## Cursos irmãos

- `~/GitHub/multidimensionais`: Análises Multidimensionais (o modelo mais
  próximo)
- `~/GitHub/bioestatistica`: Bioestatística (graduação)
- `~/GitHub/modelos-lineares-2026-2`: Análise de Dados Univariados
- `~/GitHub/intro-r`: Introdução ao R
