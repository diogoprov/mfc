#Aula 1 - Segunda

#Nesta aula vamos aprender a importar e manipular dados no R.

#carregando pacotes
library(ape)
library(phytools)
library(geiger)

#Básico do R

N<-2
n<-3
3+2
N+n
x<-1:5
mode(x)
length(x)
z<-c("ordem","superfamília","família","gênero","espécies")
z
mode(z)
class(z)
z[2]
z[c(1,3)]
z[-2]
f<-c("Male","Male","Female","Female","Female")
f<-factor(f)
f
X<-matrix(1:9,3,3)
X
X[3,2]
X[,3]
X[3,]
Y<-data.frame(z,y=1:5,x=5:1)
Y
class(Y)
#Importanto filogenias em diferentes formatos

##Arquivo em formato nexus

barb<-read.nexus("BarbetTree.nex")

##Arquivo em formato newick

cent<-read.tree("Centrarchidae.tre")

#Essa extensão .tre é a extensão que muitos programas trabalham, p.ex., o BEAST. 
#É o formato nativo do FigTree, um programa de visualização de filogenias

#Vamos explorar os formatos e entender a classe de dados
cent
str(cent)
barb
str(barb)

barb$edge

plot.phylo(barb)
nodelabels()

#newick simples em txt
gomes<-read.tree("Gomes_phylo.txt")
gomes_dat<-read.tree("Gomes_tree_datada.txt")

gomes_dat
gomes

plot(gomes_dat);axisPhylo()
plot(gomes);axisPhylo()

is.ultrametric(gomes)
is.ultrametric(gomes_dat)

branch<-compute.brlen(gomes, 1)
plot(branch);axisPhylo()
is.ultrametric(branch)
branch1<-compute.brlen(gomes, method = "Grafen", 0.8)
plot(branch1);axisPhylo()

branch4<-compute.brlen(gomes, method = "Grafen", 2)
plot(branch4);axisPhylo()

#É muito commum que nossos dados de atributos estejam em uma planilha. Quando 
#importada no R essa planilha de dados torna-se um data frame

obj<-read.csv("Centrarchidae.csv",header=TRUE,row.names=1)
head(obj)
#Vejam que os nomes das especies estao na 1a coluna, que esta nomeada
#agora precisamos extrair os dados para buccal length e guardar num vetor, mas 
#conservando o nome das especies em cada elemento do vetor. Uma forma de fazer isso 
#é

buccal.length<-setNames(obj[,"buccal.length"],rownames(obj))

gape<-obj$gape.width
names(gape)<-rownames(obj)
gape

#outra forma é usando o as.matrix
X<-read.csv("anole.data.csv",header=TRUE,row.names=1)
X
head(X)
Y<-as.matrix(X)
Y[,1]#pega os nomes das linhas e atribui aos elementos do vetor

#Tudo vai depender do formato dos dados

#---Demonstrando efeito de correlação filogenética nos resultados
set.seed(21) ## 10 ok, 2 ok, 3 ok, 6 ok, 7 ok, 
## simulate a coalescent shaped tree
tree<-rcoal(n=100)
plotTree(tree,ftype="off")
## simulate uncorrelated Brownian evolution
x<-fastBM(tree)
y<-fastBM(tree)
par(mar=c(5.1,4.1,2.1,2.1))
plot(x,y)
#test relationship without taking into account phylogeny 
fit<-lm(y~x)
abline(fit,lwd=2,lty="dashed",col="red")
anova(fit)#very significant relationship
phylomorphospace(tree,cbind(x,y),label="off",node.size=c(0.5,0.7))
abline(fit,lwd=2,lty="dashed",col="red")
ix<-pic(x,tree)
iy<-pic(y,tree)
fit1<-lm(iy~ix-1)
anova(fit1)#no significance
summary(fit1)

##---Demonstrando o phyndr
#use some external information,such as a topological hypothesis or taxonomic resource, 
#to swap out species in the tree that don't have trait data with "phylogenetically 
#equivalent" species that do.

library(ape)
library(phyndr)
library(taxonlookup)
#https://github.com/traitecoevo/phyndr-ms/blob/master/vignette/vignette.md
#Filogenia de angiospermas datada com 798 taxa
angio_phy <- read.nexus("Magallon_etal_PL_ML.nex")  # Magallón et al. 2015, Dryad doi:10.5061/dryad.k4227

#dados para plant growth form (i.e., woody v. herbaceous) para 238 spp
wood_dat<-read.csv("GlobalWoodinessDatabase.csv", row.names = 1)  # Zanne et al. 2014, Dryad doi:10.5061/dryad.63q27
head(wood_dat)
dim(wood_dat)

#We want to build a taxonomic table that includes all species in the phylogenetic tree
#AND all species in the trait data.

head(plant_lookup())
angio_spp <- unique(c(angio_phy$tip.label, rownames(wood_dat)))#Create a vector of taxa for the union of names in the tree and the data
angio_tax <- lookup_table(angio_spp, by_species=TRUE)#Get a taxonomic lookup table for this vector of names
head(angio_tax)
angio_tax_alt <- lookup_table(angio_spp)
head(angio_tax_alt)
angio_tax <- angio_tax[,c("genus", "family", "order")]#subset the taxonomy so we only consider the genus, family and order columns

#we need to supply the chronogram, the trait data with rownames set to species names, and the taxonomic table.
angio_phyndr <- phyndr_taxonomy(angio_phy, rownames(wood_dat), angio_tax)
angio_phyndr
#Encontraram dados para mais 531 spp. Ficando a filogenia fully-matched para 778 sp
#com dados de história de vida

class(angio_phyndr)
plot.phylo(angio_phy, type = "fan", cex=0.3)

#Demonstrando Rphylopars
library(Rphylopars)
library(phytools) # for simulating pure-birth phylogenies
set.seed(21) # Set the seed for reproducible results
trait_cov <- matrix(c(4,2,2.2,2,3,1.5,2.2,1.5,2.5),nrow = 3,ncol = 3)
tree <- pbtree(n = 20)
sim_data <- simtraits(v = trait_cov,tree = tree,nmissing = 10)
original<-sim_data$trait_data[,2:4] # Original data with missing values
which(is.na(original))#quais espécies tem dados faltantes
original[rowSums(is.na(original)) > 0,]#mostra só as linhas com NA
phylo.heatmap(tree,as.matrix(original),standardize=TRUE,lwd=3,
              pts=FALSE)

p_BM <- phylopars(trait_data = sim_data$trait_data,tree = sim_data$tree, model = "BM")

imput<-p_BM$anc_recon[1:20,] # Data with imputed species means
phylo.heatmap(tree,as.matrix(imput),standardize=TRUE,lwd=3,
              pts=FALSE)

p_BM$anc_var[1:20,] # Variances for each estimate
p_BM$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
p_BM$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI

#Demonstrando o rtol

library(rotl)
apes <- c("Pongo", "Pan", "Gorilla", "Hoolock", "Homo")
resolved_names <- tnrs_match_names(apes)
resolved_names
tr <- tol_induced_subtree(ott_ids=ott_id(resolved_names))
plot(tr)#;axisPhylo()

treefrogs<-c("Hypsiboas", "Bokermannohyla", "Xenohyla", "Hyla", "Trachycephalus")
resolved_names <- tnrs_match_names(treefrogs)
resolved_names
tr <- tol_induced_subtree(ott_ids=ott_id(resolved_names))
plot(tr);axisPhylo()
