#Aula 2 - Terça 

set.seed(21)#reprodutibilidade

t<-0:100 # time
sig2<-0.01

#simulando 100 traits sem uma filogenia
nsim<-100
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2)),nsim,length(t)-1)
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum)))
plot(t,X[1,],xlab="time",ylab="phenotype",ylim=c(-2,2),type="l")
apply(X[2:nsim,],1,function(x,t) lines(t,x),t=t)

#mude sig2, divida por 10 -> a variância irá diminuir
X<-matrix(rnorm(n=nsim*(length(t)-1),sd=sqrt(sig2/10)),nsim,length(t)-1)
X<-cbind(rep(0,nsim),t(apply(X,1,cumsum)))
plot(t,X[1,],xlab="time",ylab="phenotype",ylim=c(-2,2),type="l")
apply(X[2:nsim,],1,function(x,t) lines(t,x),t=t)

#Simulando Movimento Browniano sem considerar a filogenia
BM<- function(t = 100, sig2 = 0.01, sp = 100, ...) { 
       #written by Bruno Vilela Jun 2016 PCM workshop       
       par(mar = c(4, 4, 2, 2))
       t2 <- 0:t # time
       x<- rnorm(n = t, sd=sqrt(sig2))
       ## now compute their cumulative sum
       x<- c(0, cumsum(x))
       plot(t2, x, type="l", ...)
       
       for(i in 1:sp) {
              ## first, simulate a set of random deviates
              x<- rnorm(n = t, sd=sqrt(sig2))
              ## now compute their cumulative sum
              x<- c(0, cumsum(x))
              lines(t2, x, type="l", col = rainbow(100)[sample(100, 1)], ...)
       }
}

BM(100, 0.01, 1000, ylim = c(-10, 10), ylab = "Phenotype", xlab = "time",
   main = "Brownian motion")

#Simulando BM numa filogenia

set.seed(1003)
library(phytools)
tree<-make.era.map(pbtree(n=10,tip.label=LETTERS[1:10],scale=100),0:200/2)
plot(tree)
obj<-map.to.singleton(tree)
x<-fastBM(obj,internal=TRUE)
par(mar=c(5.1,4.1,2.1,2.1))
phenogram(obj,x,fsize=0.8,spread.cost=c(1,0))

#Simulando um OU
tree<-pbtree(n=26,tip.label=LETTERS,scale=10)
x<-fastBM(tree,a=0,theta=3,alpha=0.2,sig2=0.1, internal=TRUE)
phenogram(tree,x,spread.labels=TRUE,spread.cost= c(1,0.01))

#mudando o theta
y<-fastBM(tree,a=0,theta=6,alpha=0.2,sig2=0.1, internal=TRUE)
phenogram(tree,y,spread.labels=TRUE,spread.cost= c(1,0.01))

#mudando o alpha
b<-fastBM(tree,a=0,theta=3,alpha=0.8,sig2=0.1, internal=TRUE)
phenogram(tree,b,spread.labels=TRUE,spread.cost= c(1,0.01))

#Simulando Trend 
tree2<-pbtree(n=100)
w<-fastBM(tree2,mu=3,a=1,sig2=0.1, internal=TRUE) # with a trend
phenogram(tree2,w,spread.labels=TRUE,spread.cost= c(1,0.01))

#--- Comparando ajuste de modelos evolutivos
library(geiger)
library(picante)
library(phytools)

anoleData <- read.csv("anolisDataAppended.csv", row.names = 1)
head(anoleData)
dim(anoleData)
anoleTree <- read.tree("anolis.phy")
anoleSize <- anoleData[, 1]
names(anoleSize) <- rownames(anoleData)
cn<-match.phylo.data(anoleTree, anoleSize)
cn$data
(brownianModel <- fitContinuous(cn$phy, cn$data))
(OUModel <- fitContinuous(cn$phy, cn$data, model = "OU"))
(EBModel <- fitContinuous(cn$phy, cn$data, model = "EB"))

bmAICC <- brownianModel$opt$aicc
ouAICC <- OUModel$opt$aicc
ebAICC <- EBModel$opt$aicc

aicc <- c(bmAICC, ouAICC, ebAICC)
aiccD <- aicc - min(aicc)
aw <- exp(-0.5 * aiccD)
aiccW <- aw/sum(aw)
aiccW

##---Testando model adequacy
#install.packages(devtools)
#library(devtools)
#devtools::install_github("mwpennell/arbutus")
library(arbutus)
modadeq<-arbutus(brownianModel)
plot(modadeq)
modadeq1<-arbutus(OUModel)
plot(modadeq1)
modadeq2<-arbutus(EBModel)
plot(modadeq2)

#--- Ajustando modelos de evolução para caracteres discretos

sqData<-read.csv("brandley_table.csv")
head(sqData)
sqTree<-read.nexus("squamate.tre")
par(mfrow=c(1,1))
plotTree(sqTree,type="fan",lwd=1,fsize=0.5)

hindLimbs<-sqData[,"HLL"]!=0
foreLimbs<-sqData[,"FLL"]!=0
limbless<-!hindLimbs & !foreLimbs
sum(limbless)

speciesNames<-sqData[,1]
species<-sub(" ", "_", speciesNames)
names(limbless)<-species

limbless<-as.character(limbless)
names(limbless)<-species
table(foreLimbs, hindLimbs)

#ER model - simpler one = forward same as backward - same probability
#of gaining and losing limbs
name.check(sqTree, limbless)

compfilo<-match.phylo.data(sqTree, limbless)

erModel<-fitDiscrete(sqTree, compfilo$data, model="ER")
plot(erModel)

ardModel<-fitDiscrete(sqTree, limbless, model="ARD")
plot(ardModel)
#irreversible - impossible to regain limbs
irrMat<-cbind(c(NA, 1),c(0,NA))
irrModel<-fitDiscrete(sqTree, limbless, model = irrMat)

aic.vals <- setNames(c(erModel$opt$aicc, ardModel$opt$aicc, irrModel$opt$aicc),
                     c("ER","ARD","iRR"))
aic.vals
aic.w<-function(aic){
       d.aic<-aic-min(aic)
       exp(-1/2*d.aic)/sum(exp(-1/2*d.aic))
}
aic.w(aic.vals)

limbState<-paste(foreLimbs, hindLimbs)
names(limbState)<-species
erModel<-fitDiscrete(sqTree, limbState, model="ER")
plot(erModel)
symModel<-fitDiscrete(sqTree, limbState, model="SYM")
plot(symModel)
ardModel<-fitDiscrete(sqTree, limbState, model="ARD")
plot(ardModel)

#Convergence diagnostics
##frequency of best fit => proporção dos cadeias que chegam encontram a mesma solução
#quanto menor esse número pior

aicc <- c(erModel$opt$aicc, symModel$opt$aicc, ardModel$opt$aicc)
aiccD <- aicc - min(aicc)
aw <- exp(-0.5 * aiccD)
aiccW <- aw/sum(aw)
aiccW

##Exercício:
#montar a matriz de transição para limbState de forma que 
#membros podem ser perdidos, um de cada vez ou os dois ao mesmo tempo, 
#mas não ganhados novamente

#---PGLS
library(nlme)
library(geiger)
library(picante)
datos<-read.csv("Barbetdata.csv",header=TRUE,row.names=1)
head(datos)
arbol<-read.nexus("BarbetTree.nex")
obj<-name.check(arbol,datos)
obj
arbol.cortado<-drop.tip(arbol, obj$tree_not_data)
name.check(arbol.cortado,datos)
aa<-match.phylo.data(arbol.cortado, datos)
modelo1<-gls(Lnote~Lnalt, data=aa$data, correlation=corBrownian(1, arbol.cortado))
summary(modelo1)

modelo2<-gls(Lnote~Lnalt, data=aa$data, correlation=corPagel(1,arbol.cortado))
plot(modelo2)#diganóstico
summary(modelo2)

modelo3<-gls(Lnote~Lnalt+Length, data=aa$data)
modelo2<-gls(Lnote~Length, data=aa$data, correlation=corPagel(1,arbol.cortado))
summary(modelo2)

#---Testando adequabilidade dos modelos
adeq<-arbutus(modelo1)
plot(adeq)
adeq
plot(adeq)
adeq2<-arbutus(modelo2)
plot(adeq2)
adeq2
plot(Lnote~Lnalt, data=aa$data)
abline(modelo2,lwd=2,lty="dashed",col="red")
abline(modelo1,lwd=2,lty="dashed",col="blue")
abline(modelo3,lwd=2,lty="dashed",col="green")

#---ANOVA filogenética
library(phytools)
trait<-read.table("testes_matrix_r.txt", h=TRUE, dec = ",")
head(trait)
filogenia<-read.tree("filogenia_gabi.txt")

saz<-trait$sazonalidade
names(saz)<-rownames(trait)
diam<-trait$diametro
names(diam)<-rownames(trait)
repmode<-as.factor(trait$modos)
names(repmode)<-rownames(trait)
ret<-trait$reticular
names(ret)<-rownames(trait)

##---Testing the influence of reproductive modes on testicular characters
phylANOVA(filogenia, repmode, log(diam), nsim = 2000, p.adj = "BY")
phylANOVA(filogenia, repmode, log(ret))

str(trait)

summary(gls(diametro~as.factor(modos), data=dados_novos$data, correlation = corPagel(1, dados_novos$phy)))
pairs(lsmeans(a,  ~ modos))
summary(gls(reticular~as.factor(modos), data=dados_novos$data, correlation = corPagel(1, dados_novos$phy)))

ggplot(trait, aes(modos,diametro))+
       geom_col()


