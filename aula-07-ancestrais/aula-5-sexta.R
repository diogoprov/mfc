#Aula 5 sexta manhã

#---Demonstração plots phytools
library(phytools)

tree<-pbtree(n=26,tip.label=LETTERS)
plot.phylo(tree)
plotTree(tree)
x<-fastBM(tree)## simulate data under Brownian motion
x
phenogram(tree,x,spread.labels=TRUE)
dotTree(tree,x)
dotTree(tree,x,standardize=TRUE)
X<-fastBM(tree,nsim=10)
dotTree(tree,X)
phylo.heatmap(tree,X)
fancyTree(tree,type="phenogram95",x=x,spread.cost=c(1,0))
obj<-contMap(tree,x)
plot(obj,type="fan",lwd=5)
obj<-setMap(obj,invert=TRUE)
plot(obj)
obj<-setMap(obj,colors=c("white","black"))
plot(obj)

plotTree.wBars(tree,exp(x),scale=0.2)
plotTree.wBars(tree,exp(x),scale=0.2, tip.labels = TRUE)

par(mar = c(4, 4, 2, 2))
Y<-fastBM(tree,nsim=3) ## simulate 3 characters
phylomorphospace(tree,Y[,c(1,2)],xlab="trait 1",ylab="trait 2")
phylomorphospace3d(tree,Y)
fancyTree(tree,type="scattergram",X=X)

#Mapas
set.seed(1)
tree<-pbtree(n=26,scale=100)
tree$tip.label<-LETTERS[26:1]
lat<-fastBM(tree,sig2=10,bounds=c(-90,90))
long<-fastBM(tree,sig2=80,bounds=c(-180,180))
# now plot projection
xx<-phylo.to.map(tree,cbind(lat,long),plot=FALSE)
plot(xx,type="phylogram",asp=1.3,mar=c(0.1,0.5,3.1,0.1))

#---Demonstração plots phylosignal
#http://www.francoiskeck.fr/phylosignal/demo_plots.html

library(phylosignal)
data(carni19)
tre <- read.tree(text = carni19$tre)
dat <- data.frame(carni19$bm)
dat$bm <- rTraitCont(tre)
p4d <- phylo4d(tre, dat)
barplot(p4d)
dotplot(p4d)
gridplot(p4d)
gridplot(p4d, tree.type = "fan", tip.cex = 0.6, show.trait = FALSE)

mat.e <- matrix(abs(rnorm(19 * 3, 0, 0.5)), ncol = 3,
                dimnames = list(tipLabels(p4d), names(tdata(p4d))))
barplot(p4d, error.bar.sup = mat.e, error.bar.inf = mat.e)

mat.col <- ifelse(tdata(p4d, "tip") < 0, "red", "grey35")
barplot(p4d, center = FALSE, bar.col = mat.col)

#---Phylobabase
library(phylobase)
filo<-read.tree("filogenia_bocaina.txt")
trait<-read.csv2("canto.csv", h=TRUE, dec = ".")
trait<-trait[,-c(1,2,4,8)]
head(trait)
habitat<-trait[,c(1,2)]
head(habitat)

#---phylo PCA e PVR
obj<-phylo4d(filo, trait[,-7])
table.phylo4d(obj, box = FALSE, cex.label = 0.7)
phylo.pca<-ppca(obj, method = "patristic",scannf=FALSE, nfposi=2)
plot(phylo.pca, useLag = TRUE)
scatter(phylo.pca, useLag = TRUE)
vectors<-phylo.pca$li#eigenvectors of PCA
espectral=as.matrix(vectors)[,1]
names(pulses)=rownames(vectors)
temporal=as.matrix(vectors)[,2]
names(temporal)=rownames(vectors)
new.data2<-treedata(filo, habitat)

x<-PVRdecomp(filo)
psrcurve<-PSR(x, trait = vectors, null.model = TRUE, Brownian.model = TRUE, times = 999)
PSRplot(psrcurve, info = "both")#Top plot refers to Null model; bottom plot refers to neutral model (Brownian Motion)
anal_pvr<-PVR(psrcurve, trait = espectral, envVar = new.data2$data, method = "psr")
anal_pvr@PVR$R2
anal_pvr@PVR$p
VarPartplot(anal_pvr)#plot of variation partitioning
str(anal_pvr@VarPart)#R2 of each component of variation partitioning
