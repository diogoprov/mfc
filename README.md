# Métodos Filogenéticos Comparativos

Materiais da disciplina de pós-graduação **Métodos Filogenéticos
Comparativos** (201521439, 45 h) do
[Mestrado em Biologia Animal](https://ppgbioanimal.ufms.br/) da Universidade
Federal de Mato Grosso do Sul.

**Site publicado: <https://provetelab.org/mfc/>**

## O que há aqui

Seis aulas, cada uma em sua própria pasta, contendo a página da aula
(`index.qmd`), os scripts em R usados em sala, os conjuntos de dados e os
*walkthroughs* — versões HTML comentadas do código, para o aluno acompanhar
linha a linha.

| Aula | Tópico |
|-----:|--------|
| 1 | Introdução e desenho de estudos comparativos |
| 2 | Modelos de evolução de atributos contínuos e PGLS |
| 3 | Estimativa de estado ancestral e sinal filogenético |
| 4 | Correlação de atributos e modelo Ornstein-Uhlenbeck multi-regime |
| 5 | Métodos baseados em autovetores e visualização |
| 6 | Métodos para atributos multivariados |

Os slides em PDF ficam em `slides/`.

## Como o site é construído

Site [Quarto](https://quarto.org) publicado no GitHub Pages por GitHub
Actions (`.github/workflows/publish.yml`). **A pasta `docs/` não é
versionada** — ela é gerada a cada push. Para publicar uma alteração, basta
editar o `.qmd`, commitar e dar push.

Para trabalhar localmente:

```bash
quarto preview      # pré-visualiza com recarga automática
quarto render       # gera docs/
```

O código R das aulas é exibido mas **não executado** durante o render (todos
os walkthroughs trazem `eval: false`), então não é preciso ter os pacotes
instalados para construir o site.

## Origem

Estes materiais viviam em `teaching/mfc/` no repositório do site do
laboratório ([diogoprov.github.io](https://github.com/diogoprov/diogoprov.github.io)),
de onde foram copiados como ponto de partida deste repositório. O histórico
anterior continua disponível lá.

---

[Biodiversity Synthesis Lab](https://provetelab.org/) · Diogo B. Provete ·
Instituto de Biociências, UFMS
