set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/ICBP2011/Vs1/ICBP2011.RSS0.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/ICBP2011/Vs1/ICBP2011.dtlesssignif.vs1.annot.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0,]

~~~
> dim(gl)
[1] 796  20
> gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
>
>
> dim(gl)
[1] 783  20
~~~

Z = cbind(gl$Z.SBP,gl$Z.DBP,gl$Z.PP,gl$Z.MAP)
n = cbind(gl$n_SBP, gl$n_DBP, gl$n_PP, gl$n_MAP)
n = apply(n,1,min)
gl$nmin=n
sigmaa=c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)
lbf=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa)
lbf.bigmat=do.call(rbind,lbf$lbf) #lbf$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows 

#collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
collapse = function(x,nsigmaa){return(apply(matrix(x,ncol=nsigmaa),1,sum))}
  
#compute empirical bayes priors from global lipids hit SNPs
gl.glhits = gl[gl$annot==1,] # subset of SNPs reported in Teslovich 2010 paper
lbf.glhits = lbf.bigmat[,gl$annot==1]
ebprior.glhits = em.priorprobs(lbf.glhits,lbf$prior,100)
ebprior.glhits2 = em.priorprobs(lbf.glhits,lbf$prior*runif(length(lbf$prior)),100)
ebprior.glhits3 = em.priorprobs(lbf.glhits,lbf$prior*runif(length(lbf$prior)),100)
ebprior.glhits4 = em.priorprobs(lbf.glhits,lbf$prior*runif(length(lbf$prior)),100)

ebprior.glhits.collapse =collapse(ebprior.glhits,length(sigmaa))
ebprior.glhits.collapse2 =collapse(ebprior.glhits2,length(sigmaa))
ebprior.glhits.collapse3 =collapse(ebprior.glhits3,length(sigmaa))
ebprior.glhits.collapse4 =collapse(ebprior.glhits4,length(sigmaa))

pp.glhits = posteriorprob(lbf.glhits,ebprior.glhits) #posterior prob on models for gl hits
pp.glhits.collapse =  apply(pp.glhits,2,collapse, nsigmaa=length(sigmaa))
  
# this returns a list with elements  pU pD and pI
marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

#looking at which models are favored by the prior
cumsum(sort(ebprior.glhits.collapse,decreasing=TRUE))
lbf$gamma[order(ebprior.glhits.collapse,decreasing=TRUE),]
modelmatrix = cbind(lbf$gamma,ebprior.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
modelmatrix = data.frame(cbind(modelmatrix,cumsum(modelmatrix[,5])))
colnames(modelmatrix)= c("SBP","DBP","PP","MAP","p","cump")

allassoc=(apply((modelmatrix[,1:4]>0),1,sum)==4) #vector of which represent situations in which all 4 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:4]>0),1,sum)==3) #vector of which represent situations in which 3 phenotypes are associated

sum(modelmatrix[allassoc,5]) #0.9756873
sum(modelmatrix[allbut1assoc,5]) #0.02431274

#look at weight on each of SBP, DBP, PP, MAP being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5]) 
# 1.711418e-54, 6.298546e-39, 0.02431274, 2.91612e-55 

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notMAP, or notPP, or notDBP or notSBP)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"gl.bestmodel.vs2.txt",quote=F,sep=" & ",row.names=F)

#     [,1]         [,2] [,3]
# [1,] "rs17477177" "4"  "0.6"

#Previous results -- kept here as a reminder to try and get gene annotations at some point? This information (gene annotations per SNP) was provided for Matthew from someone else previously
#      [,1]              [,2] [,3]               
#  [1,] "CILP2_TC:TG:LDL" "3"  "0.545021355642582"
#  [2,] "GCKR_TG:TC"      "3"  "0.801152441063756"
#  [3,] "LPL_TG:HDL"      "2"  "0.991234434984749"
#  [4,] "LIPC_HDL:TC:TG"  "2"  "0.912277207888376"
#  [5,] "MLXIPL_TG:HDL"   "2"  "0.974368366454864"
#  [6,] "HNF4A_HDL:TC"    "4"  "0.83326782439584" 
#  [7,] "HPR_TC:LDL"      "3"  "0.583327091921807"
#  [8,] "CAPN3_TG"        "2"  "0.628280193423631"
#  [9,] "SORT1_LDL:TC"    "4"  "0.548904956312579"
# [10,] "LDLR_LDL:TC"     "4"  "0.614955457455026"
# [11,] "LIPG_HDL:TC"     "4"  "0.780055882099995"


#Compare multivariate results with univariate results for the GL SNPs.
lbf.av.glhits= lbf.av(lbf.glhits,ebprior.glhits)
nsigma=length(sigmaa)
origprior = rep(c(0,lbf$prior[-1]),nsigma)
origprior = normalize(origprior)
lbf.av.origprior.glhits = lbf.av(lbf.glhits,origprior)
lbf.uni.glhits = lbf.uni(lbf.glhits,lbf$gamma)
lbf.all.glhits = lbf.all(lbf.glhits,lbf$gamma)


#Comparisons with lbf.av.glhits are problematic because the lbf.av.glhits
#prior are fitted to data (particularly for estimating sigmaa)
#whereas the lbf.uni are not (but could be if wanted?)
#So stick to comparisons of lbfs computed without learning prior from data
pdf("nana1")
layout(rbind(c(1,2)))
plot(lbf.uni.glhits,lbf.av.origprior.glhits,xlim=c(4,15),ylim=c(4,15))
abline(a=0,b=1)
plot(lbf.uni.glhits,lbf.all.glhits,xlim=c(4,15),ylim=c(4,15))
abline(a=0,b=1)
dev.off()

#this establishes a threshold for subsequent analysis
min(lbf.av.glhits)
#[1] 4.058776


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]


sub = gl[gl$nmin>20000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.058776 & sub$nmin>20000,c(4,2:3,7,22:25,30)]
newhits = sub[sub$annot==0 & sub$lbfav>4.058776 & sub$nmin>20000,]

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>20000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>4.058776 & sub$nmin>20000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#10 in new results
#bestclass                    
#[1,] "rs10505863" "1"       "0.999992288214027"
#[2,] "rs11066320" "1"       "0.999999999978562"
#[3,] "rs12509057" "1"       "0.999975416673578"
#[4,] "rs1486936"  "1"       "0.998596495289966"
#[5,] "rs17677264" "1"       "0.99999952070509"
#[6,] "rs199205"   "1"       "0.999986683453266"
#[7,] "rs2052313"  "1"       "0.999999645173759"
#[8,] "rs2969070"  "1"       "0.999999995877814"
#[9,] "rs6544619"  "1"       "0.999999978119228"
#[10,] "rs783622"   "1"       "0.999709482253066"

#Old results
#                   bestclass                    
#  [1,] "rs12739698" "1"       "0.811711793061224"
#  [2,] "rs267733"   "1"       "0.783023729441067"
#  [3,] "rs10490632" "1"       "0.99998311557302" 
#  [4,] "rs13326165" "1"       "0.97851833794812" 
#  [5,] "rs762861"   "1"       "0.787090817832785"
#  [6,] "rs2746150"  "1"       "0.999366545708106"
#  [7,] "rs998584"   "1"       "0.835529167903279"
#  [8,] "rs6951245"  "1"       "0.931741017952404"
#  [9,] "rs4722551"  "1"       "0.69763687764716" 
# [10,] "rs17134533" "2"       "0.753742777779664"
# [11,] "rs10904908" "1"       "0.792803343102496"
# [12,] "rs970548"   "1"       "0.776315520754646"
# [13,] "rs1408579"  "1"       "0.718123030686882"
# [14,] "rs11246602" "1"       "0.959441987095743"
# [15,] "rs11229252" "1"       "0.973330722568743"
# [16,] "rs11227638" "1"       "0.997884622570898"
# [17,] "rs499974"   "1"       "0.99640117575064" 
# [18,] "rs4942505"  "1"       "0.631116110836567"
# [19,] "rs10422101" "1"       "0.898137410416297"

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/ICBP2011/Vs1/GLC.MTedits.ForICBP2011.vs1.ForCSHLPoster.ModelsVslbfavNewOld.vs1.jpeg", res=300, width=5000, height=4000)
par(mfrow=c(2,2))

head(modelmatrix)
modelmatrix.paste <- apply(modelmatrix[,1:4], 1, function(x) paste(x, collapse="_"))
table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),5]

snpmodels <- cbind(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)], apply(ppmatrix.newhits, 2, max), apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] }), apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-2)[length(x)-2] }))
topmodels <- c()

for (i in unique(snpmodels[,1])) {
	snpmodels.subset <- matrix(snpmodels[snpmodels[,1]==i,], nrow=length(which(snpmodels[,1]==i)))
	topmodels <- rbind(topmodels, cbind(i, nrow(snpmodels.subset), mean(as.numeric(snpmodels.subset[,2])), sd(as.numeric(snpmodels.subset[,2])), mean(as.numeric(snpmodels.subset[,2]) - as.numeric(snpmodels.subset[,3])), sd(as.numeric(snpmodels.subset[,2]) - as.numeric(snpmodels.subset[,3])), modelmatrix[modelmatrix.paste == i,5]))
}
topmodels <- topmodels[order(as.numeric(topmodels[,2]), decreasing=TRUE),]

#topmodels <- cbind(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])), table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]), modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),5])
#topmodels <- topmodels[order(as.numeric(topmodels[,2]), decreasing=FALSE),]

#plot difference between first best model and second best model against lbfav 
plot(newhits$lbfav, apply(ppmatrix.newhits, 2, max), main="New Hits", xlab="lbfav", ylab="max posteriorprob")
abline(a=0, b=.1, col="BLACK")
abline(lm(apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav), col="RED")
plot(newhits$lbfav, (apply(ppmatrix.newhits, 2, max) - apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })), main="New Hits", xlab="lbfav", ylab="max posteriorprob - 2nd max posteriorprob")
abline(a=0, b=.1, col="BLACK")
abline(lm((apply(ppmatrix.newhits, 2, max) - apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ newhits$lbfav), col="RED")
summary(lm(apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav))
summary(lm((apply(ppmatrix.newhits, 2, max) - apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ newhits$lbfav))

oldsnpmodels <- cbind(modelmatrix.paste[apply(ppmatrix, 2, which.max)], apply(ppmatrix, 2, max), apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] }), apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-2)[length(x)-2] }))
oldtopmodels <- c()

for (i in unique(oldsnpmodels[,1])) {
        oldsnpmodels.subset <- matrix(oldsnpmodels[oldsnpmodels[,1]==i,], nrow=length(which(oldsnpmodels[,1]==i)))
	oldtopmodels <- rbind(oldtopmodels, cbind(i, nrow(oldsnpmodels.subset), mean(as.numeric(oldsnpmodels.subset[,2])), sd(as.numeric(oldsnpmodels.subset[,2])), mean(as.numeric(oldsnpmodels.subset[,2]) - as.numeric(oldsnpmodels.subset[,3])), sd(as.numeric(oldsnpmodels.subset[,2]) - as.numeric(oldsnpmodels.subset[,3])), modelmatrix[modelmatrix.paste == i,5])) 
}
oldtopmodels <- oldtopmodels[order(as.numeric(oldtopmodels[,2]), decreasing=TRUE),]

#uniformprior thing

#plot difference between first best model and second best model against lbfav 
plot(gl[gl$annot==1,]$lbfav, apply(ppmatrix, 2, max), main="Old Hits", xlab="lbfav", ylab="max posteriorprob")
abline(a=0, b=.1, col="BLACK")
abline(lm(apply(ppmatrix, 2, max) ~ gl[gl$annot==1,]$lbfav), col="RED")
plot(gl[gl$annot==1,]$lbfav, (apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })), main="Old Hits", xlab="lbfav", ylab="max posteriorprob - 2nd max posteriorprob")
abline(a=0, b=.1, col="BLACK")
abline(lm((apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ gl[gl$annot==1,]$lbfav), col="RED")
summary(lm(apply(ppmatrix, 2, max) ~ gl[gl$annot==1,]$lbfav))
summary(lm((apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ gl[gl$annot==1,]$lbfav))

dev.off()


#~~~
#> head(modelmatrix)
#SBP DBP PP MAP          p      cump
#1   1   2  1   1 0.41777495 0.4177749
#2   2   1  2   1 0.24279598 0.6605709
#3   1   1  1   1 0.16404167 0.8246126
#4   1   1  2   1 0.08613134 0.9107439
#5   1   1  1   2 0.06494331 0.9756873
#6   1   0  1   2 0.02431274 1.0000000
#> modelmatrix.paste <- apply(modelmatrix[,1:4], 1, function(x) paste(x, collapse="_"))
#> table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
#
#1_2_1_1 2_1_2_1
#5       5
#> modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),5]
#[1] 0.4177749 0.2427960
#>
#> topmodels
#i                                                                        
#[1,] "1_2_1_1" "5" "0.705067923682195" "0.176614190985742" "0.550918323631555"
#[2,] "2_1_2_1" "5" "0.612730537308071" "0.157083608119753" "0.403194173973646"
#
#[1,] "0.223748658522453" "0.417774947763434"
#[2,] "0.267058373429467" "0.242795983159925"
#> oldtopmodels
#i                                                      
#[1,] "1_2_1_1" "18" "0.75364615631689"  "0.165377340486151" 
#[2,] "2_1_2_1" "11" "0.67351484602083"  "0.159513178044761" 
#[3,] "1_1_1_1" "3"  "0.417842923458669" "0.0920116402108855"
#[4,] "1_1_2_1" "1"  "0.562144728471013" NA                  
#[5,] "1_0_1_2" "1"  "0.598201621019533" NA                  
#[6,] "1_1_1_2" "1"  "0.397649151405326" NA                  
#
#[1,] "0.599576534434112" "0.247506951504211" "0.417774947763434" 
#[2,] "0.500463227484356" "0.228729509928718" "0.242795983159925" 
#[3,] "0.123130551127911" "0.127067242559379" "0.164041670998494" 
#[4,] "0.330290214408235" NA                  "0.0861313446973542"
#[5,] "0.355628684471552" NA                  "0.0243127383256028"
#[6,] "0.144281758144174" NA                  "0.0649433119638657"
#> 
#> summary(lm(apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav)
#
#Residuals:
#Min      1Q  Median      3Q     Max
#-0.2450 -0.1411  0.0208  0.1411  0.1923
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.78817    0.20226   3.897  0.00456 **
#newhits$lbfav -0.02213    0.03338  -0.663  0.52591
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.1703 on 8 degrees of freedom
#Multiple R-squared:  0.0521,    Adjusted R-squared:  -0.06639
#F-statistic: 0.4397 on 1 and 8 DF,  p-value: 0.5259
#
#> summary(lm((apply(ppmatrix.newhits, 2, max) - apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ newhits$lbfav))
#
#Call:
#lm(formula = (apply(ppmatrix.newhits, 2, max) - apply(ppmatrix.newhits,
#2, function(x) {
#sort(x, partial = length(x) - 1)[length(x) - 1]
#})) ~ newhits$lbfav)
#
#Residuals:
#Min      1Q  Median      3Q     Max
#-0.3645 -0.1923  0.0181  0.2043  0.2730
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.77714    0.28828   2.696   0.0273 *
#newhits$lbfav -0.05138    0.04758  -1.080   0.3117
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2427 on 8 degrees of freedom
#Multiple R-squared:  0.1272,    Adjusted R-squared:  0.01814
#F-statistic: 1.166 on 1 and 8 DF,  p-value: 0.3117
#
#> summary(lm(apply(ppmatrix, 2, max) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix, 2, max) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.33126 -0.13863  0.05503  0.14195  0.25567
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.647911   0.088672   7.307 2.19e-08 ***
#gl[gl$annot == 1, ]$lbfav 0.003591   0.009400   0.382    0.705
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.1855 on 33 degrees of freedom
#Multiple R-squared:  0.004402,  Adjusted R-squared:  -0.02577
#F-statistic: 0.1459 on 1 and 33 DF,  p-value: 0.7049
#
#> summary(lm((apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = (apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) {
#sort(x, partial = length(x) - 1)[length(x) - 1]
#})) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min      1Q  Median      3Q     Max
#-0.4798 -0.2103  0.0783  0.1968  0.3966
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.477161   0.128221   3.721 0.000736 ***
#gl[gl$annot == 1, ]$lbfav 0.002579   0.013593   0.190 0.850685
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2683 on 33 degrees of freedom
#Multiple R-squared:  0.00109,   Adjusted R-squared:  -0.02918
#F-statistic: 0.036 on 1 and 33 DF,  p-value: 0.8507
#~~~



pdf("plots.bychr.vs2.pdf")
for(i in 1:22){
 if(sum(sub$chr==i)>0){
plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>50000))
abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>4.058776],col=3)
}
}
dev.off()

#try to plot all hits to see relationships?
tophits =  sub[sub$lbfav>4.058776 & sub$nmin>20000,]
betahat = tophits[,c(8,11,14,17)]
betahat.scaled = betahat/apply(betahat,1,sd)
betahat.scaled.pr = prcomp(betahat.scaled,center=F)
pdf("plots.betahats.vs2.pdf")
plot(abs(betahat.scaled.pr$x[,1]),abs(betahat.scaled.pr$x[,2]),type="n")
text(abs(betahat.scaled.pr$x[,1]),abs(betahat.scaled.pr$x[,2]),as.character(tophits$chr))
dev.off()

hla = gl[gl$chr==6 & gl$pos<40e6 & gl$pos>28e6,]
hla.z = hla[,22:25]

# > plot(hla.pr$x[,1])
# > plot(hla.pr$x[,2])
# > plot(hla.pr$x[,2]^2)
# > plot(hla.pr$x[,2]^2,hla.pr$x[,1]^2)
# > image(hla.z)
# Error in image.default(hla.z) : 'z' must be a matrix
# > image(as.matrix(hla.z))
# > image(as.matrix(hla.z^2))
# > image(as.matrix(abs(hla.z)))
# 

zhla = head(gl[gl$chr==6 & gl$pos<40e6 & gl$pos>28e6,],600)[,22:25]
pdf("plots.cor.zhla.vs2.pdf")
image(cor(t(zhla))^2)
dev.off()

pdf("plots.chr1lbfav.vs2.pdf")
plot(gl$pos[gl$chr==1],gl$lbfav[gl$chr==1],ylim=c(0,10),xlim=c(27102620-10^5,27102620+10^5))
dev.off()

write.table(file="newtophits.vs2.txt",cbind(sub[sub$lbfav>4.058776 & sub$nmin>20000 & sub$annot==0,c(1:3,13)],round(sub[sub$lbfav>4.058776 & sub$nmin>20000 & sub$annot==0,c(4,15:22)],digits=5)),quote=FALSE,sep= " ", row.names=FALSE)

c(4,2:3,7,22:25,30)

#newhits = gl[gl$annot==0 & gl$lbfav>4.058776 & gl$nmin>20000,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)




