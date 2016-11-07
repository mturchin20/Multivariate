library(ggplot2)
library(reshape2)
set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.4Phenos_2.RSS0.vs3.SignCrrct.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.4Phenos_2.dtlesslesssignif.vs3.SignCrrct.vs1.annot.vs1.MAF.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0,]

Z = cbind(gl$Z.dHb,gl$Z.OxHb,gl$Z.Pulse,gl$Z.LvBrth)

#~~~
> dim(gl)
[1] 2202   20
> gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
> dim(gl)
[1] 2202   20
> VYY
           Z.dHb      Z.OxHb     Z.Pulse    Z.LvBrth
[1,]  1.00000000  0.01181027  0.14204221 -0.02747314
[2,]  0.01181027  1.00000000  0.07792366 -0.04976397
[3,]  0.14204221  0.07792366  1.00000000 -0.10999706
[4,] -0.02747314 -0.04976397 -0.10999706  1.00000000
> head(gl)
         snp chr       pos   maf        p_OxHb n_OxHb      p_dHb n_dHb
1 rs10006686   4 177993314 0.213  7.461186e-05    915 -0.6300364   915
2 rs10011618   4 177999800 0.211  9.727020e-05    916 -0.5423901   916
3 rs10012477   4 177994703 0.213 -7.461186e-05    915  0.6300364   915
4 rs10012578   4 177994820 0.213  7.461186e-05    915 -0.6300364   915
5 rs10012654   4 177994748 0.213  7.461186e-05    915 -0.6300364   915
6 rs10016312   4 177996423 0.212 -8.639529e-05    917  0.6484305   917
     p_Pulse n_Pulse   p_LvBrth n_LvBrth annot    Z.OxHb      Z.dHb    Z.Pulse
1  0.3517754     914  0.3150894      975     0  3.961079 -0.4816756  0.9311511
2  0.3756171     915  0.3832814      976     0  3.897302 -0.6092027  0.8860008
3 -0.3517754     914 -0.3150894      975     0 -3.961079  0.4816756 -0.9311511
4  0.3517754     914  0.3150894      975     0  3.961079 -0.4816756  0.9311511
5  0.3517754     914  0.3150894      975     0  3.961079 -0.4816756  0.9311511
6 -0.3522391     916 -0.2797214      977     0 -3.925930  0.4559436 -0.9302549
    Z.LvBrth   mvstat      mvp     unip
1  1.0046002 18.08029 3.373256 4.127192
2  0.8718662 17.35590 3.223973 4.012020
3 -1.0046002 18.08029 3.373256 4.127192
4  1.0046002 18.08029 3.373256 4.127192
5  1.0046002 18.08029 3.373256 4.127192
6 -1.0809454 17.97125 3.350766 4.063510
~~~

n = cbind(gl$n_dHb, gl$n_OxHb, gl$n_Pulse, gl$n_LvBrth)
n = apply(n,1,min)
gl$nmin=n
sigmaa=c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)
lbf=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa)
lbf.bigmat=do.call(rbind,lbf$lbf) #lbf$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows 

#collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
collapse = function(x,nsigmaa){return(apply(matrix(x,ncol=nsigmaa),1,sum))}
  
#compute empirical bayes priors from GWAS hit SNPs
gl.glhits = gl[gl$annot==1,] # subset of GWAS SNPs
#Expecting way more than this...need to specifically do something to forcefully insert top SNPs?
#~~~
> gl.glhits
             snp chr       pos   maf        p_OxHb n_OxHb      p_dHb n_dHb
1178 rs372272284   2  46584859 0.248  5.713546e-09    914 -0.8825847   914
1595   rs6711319   2 179717217 0.247 -4.038030e-01    921  0.6952335   921
       p_Pulse n_Pulse      p_LvBrth n_LvBrth annot    Z.OxHb      Z.dHb
1178 0.7307703     913  7.775514e-01      974     1 5.8249332 -0.1476935
1595 0.2904827     920 -2.892655e-09      981     1 0.8348485 -0.3917628
       Z.Pulse  Z.LvBrth   mvstat      mvp     unip nmin
1178  0.344101 0.2825114 34.29920 6.766246 8.243094  913
1595 -1.057063 5.9375768 36.80627 7.296126 8.538703  920
#~~~
lbf.glhits = as.matrix(lbf.bigmat[,gl$annot==1])
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
colnames(modelmatrix)= c("dHb","OxHb","Pulse","LvBrth","p","cump")

allassoc=(apply((modelmatrix[,1:4]>0),1,sum)==4) #vector of which represent situations in which all 4 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:4]>0),1,sum)==3) #vector of which represent situations in which 3 phenotypes are associated

sum(modelmatrix[allassoc,5]) #9.147548e-08
sum(modelmatrix[allbut1assoc,5]) #0.4771804

#~~~
> sum(modelmatrix[allassoc,5]) #1.537052e-11
[1] 9.147548e-08
> sum(modelmatrix[allbut1assoc,5]) #0.06239647
[1] 0.4771804
#~~~

#look at weight on each of LvBrth, Pulse, OxHb, dHb being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5]) 
# 0.02148162, 1.998338e-09, 0.455698, 7.03089e-07

#~~~
> sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5])
[1] 0.02148162
> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5])
[1] 1.998338e-09
> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5])
[1] 0.455698
> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5])
[1] 7.03089e-07
#~~~

#20160810 NOTE -- Only one SNP so hard to do these calculations?

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLvBrth, or notPulse, or notOxHb or notdHb)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(matrix(ppmatrix)[allassoc,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,4]==0,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,3]==0,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,2]==0,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"Choongwon2016.RanHalfFlip.4Phenos_2.gl.bestmodel.vs3.SignCrrct.vs1.txt",quote=F,sep=" & ",row.names=F)

#~~~

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
#~~~

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
min(lbf.av.origprior.glhits)
#~~~
> min(lbf.av.glhits)
[1] 5.924648
> min(lbf.av.origprior.glhits)
[1] 4.503723
#~~~




lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]

#lbf.av.all.flatprior = lbf.av(lbf.bigmat, rep(lbf$prior,nsigma))
lbf.av.all.unifprior = lbf.av(lbf.bigmat, normalize(rep(c(0,lbf$prior[-1]),nsigma))) 
gl$lbfavunif = lbf.av.all.unifprior




> lbf.av
function(lbf,prior){
  max = apply(lbf,2,max) #remove the max, and store it, to prevent overflow
  lbf=centered.lbf(lbf)
  bfav = prior * 10^lbf
  return(log10(apply(bfav,2,sum))+max)
}

matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE)

#This one goes across the 14 sigmaa rows per gamma to take the mean, producing a single per gamma value across each SNP -- final matrix of 81 x 84 (gammas x SNPs)
MeanAcrossSigmaas <- function(lbf, ngamma, nsigmaa) {
        lbfMeanOverSigmas <- matrix(0, ncol=ncol(lbf), nrow=ngamma)
        for (i in 1:ngamma) {
                coords <- seq.int(from=i, by=ngamma, length.out=nsigmaa)
                max <- apply(lbf[coords,], 2, max)
                lbf[coords,] <- lbf[coords,] - matrix(max, nrow=nrow(lbf[coords,]), ncol=ncol(lbf[coords,]), byrow=TRUE)
                lbfMeanOverSigmas[i,] <- log10(apply(10^lbf[coords,], 2, mean)) + max
        }
        return(lbfMeanOverSigmas)
}
SumAcrossSigmaas.pp.NoMax <- function(pp, ngamma, nsigmaa) {
        ppSumOverSigmas <- matrix(0, ncol=ncol(pp), nrow=ngamma)
        for (i in 1:ngamma) {
                coords <- seq.int(from=i, by=ngamma, length.out=nsigmaa)
                ppSumOverSigmas[i,] <- apply(pp[coords,], 2, sum)
        }
        return(ppSumOverSigmas)
}
MeanAcrossSigmaas.wPriorAvg <- function(lbf, prior, ngamma, nsigmaa) {
        lbfMeanOverSigmas <- matrix(0, ncol=ncol(lbf), nrow=ngamma)
        for (i in 1:ngamma) {
                coords <- seq.int(from=i, by=ngamma, length.out=nsigmaa)
                max <- apply(lbf[coords,], 2, max)
                lbf[coords,] <- lbf[coords,] - matrix(max, nrow=nrow(lbf[coords,]), ncol=ncol(lbf[coords,]), byrow=TRUE)
		lbfMeanOverSigmas[i,] <- log10(apply(prior[coords,] * 10^lbf[coords,], 2, mean)) + max
        }
        return(lbfMeanOverSigmas)
}
SumAcrossSigmaas.wPriorAvg <- function(lbf, prior, ngamma, nsigmaa) {
        lbfSumOverSigmas <- matrix(0, ncol=ncol(lbf), nrow=ngamma)
        for (i in 1:ngamma) {
                coords <- seq.int(from=i, by=ngamma, length.out=nsigmaa)
                max <- apply(lbf[coords,], 2, max)
                lbf[coords,] <- lbf[coords,] - matrix(max, nrow=nrow(lbf[coords,]), ncol=ncol(lbf[coords,]), byrow=TRUE)
                lbfSumOverSigmas[i,] <- log10(apply(prior[coords,] * 10^lbf[coords,], 2, sum)) + c(max - log10(1/sum(prior[coords,1])))
        }
        return(lbfSumOverSigmas)
}
#SumAcrossSigmaas.wPriorAvg.Adjusted <- function(lbf, prior, ngamma, nsigmaa) {
#        lbfSumOverSigmas <- matrix(0, ncol=ncol(lbf), nrow=ngamma)
#        for (i in 1:ngamma) {
#                coords <- seq.int(from=i, by=ngamma, length.out=nsigmaa)
#                max <- apply(lbf[coords,], 2, max)
#                lbf[coords,] <- lbf[coords,] - matrix(max, nrow=nrow(lbf[coords,]), ncol=ncol(lbf[coords,]), byrow=TRUE)
#                lbfSumOverSigmas[i,] <- log10(apply(prior[coords,] * 10^lbf[coords,], 2, sum)) + c(max - log10(1/sum(prior[coords,1]))) + log10(apply(10^lbf[coords,], 2, mean)
#        }
#        return(lbfSumOverSigmas)
#}



#> seq.int(from = 1, by=27, length.out=14)
# [1]    1   82  163  244  325  406  487  568  649  730  811  892  973 1054

prior is column of all values across gammas (with multiple entries per sigmaa) so per row it's the same priors for each SNP





#sub = gl[gl$nmin>50000,]
sub = gl
#l=indephits(sub$lbfav,sub$chr,sub$pos)
l.lbfav=indephits(sub$lbfav,sub$chr,sub$pos)
l.lbfavunif=indephits(sub$lbfavunif,sub$chr,sub$pos)
sub.lbfav=sub[l.lbfav==1,]
sub.lbfavunif=sub[l.lbfavunif==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>5.924648 & sub$nmin>100,c(1:4,11:14,19:20)]
sub.lbfav[sub.lbfav$annot==0 & sub.lbfav$lbfav>5.924648,]
sub.lbfavunif[sub.lbfavunif$annot==0 & sub.lbfavunif$lbfavunif>4.503723,]

#~~~
> dim(sub)
[1] 2202   23
> l.lbfav=indephits(sub$lbfav,sub$chr,sub$pos)
> l.lbfavunif=indephits(sub$lbfavunif,sub$chr,sub$pos)
> sub.lbfav=sub[l.lbfav==1,]
> sub.lbfavunif=sub[l.lbfavunif==1,]
> dim(sub.lbfav)
[1] 311  23
> dim(sub.lbfavunif)
[1] 298  23
> sub.lbfav[sub.lbfav$annot==0 & sub.lbfav$lbfav>5.924648,]
 [1] snp       chr       pos       maf       p_OxHb    n_OxHb    p_dHb    
 [8] n_dHb     p_Pulse   n_Pulse   p_LvBrth  n_LvBrth  annot     Z.OxHb   
[15] Z.dHb     Z.Pulse   Z.LvBrth  mvstat    mvp       unip      nmin     
[22] lbfav     lbfavunif
<0 rows> (or 0-length row.names)
> sub.lbfavunif[sub.lbfavunif$annot==0 & sub.lbfavunif$lbfavunif>4.503723,]
 [1] snp       chr       pos       maf       p_OxHb    n_OxHb    p_dHb    
 [8] n_dHb     p_Pulse   n_Pulse   p_LvBrth  n_LvBrth  annot     Z.OxHb   
[15] Z.dHb     Z.Pulse   Z.LvBrth  mvstat    mvp       unip      nmin     
[22] lbfav     lbfavunif
<0 rows> (or 0-length row.names)
#~~~

#extract lbfs for all the new hits
#lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.bigmat
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>5.924648 & sub$nmin>100]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notPulse, or notOxHb or notdHb)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.newtophits.vs3.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)







#### in results






lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 81, 14)
lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, function(x) { sum1 <- sum(x); if (sum1 == length(x)) { sum1 <- 1; }; return(sum1); })), lbf.gl)
##lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 81, 14)
lbf.gl.prior <- SumAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 81, 14)
lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, function(x) { sum1 <- sum(x); if (sum1 == length(x)) { sum1 <- 1; }; return(sum1); })), lbf.gl.prior)
lbf.gl.prior.adjusted <- SumAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(rep(c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); }))/nsigma,nsigma), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 81, 14)
lbf.gl.prior.adjusted.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior.adjusted, 1, function(x) { sum1 <- sum(x); if (sum1 == length(x)) { sum1 <- 1; }; return(sum1); })), lbf.gl.prior.adjusted) 
rep(c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); }))/10,10)[1:5] 
posterior.gl.NoSum <- posteriorprob(lbf.bigmat, normalize(rep(c(0,lbf$prior[-1]),nsigma)))
posterior.gl.Sum <- SumAcrossSigmaas.pp.NoMax(posterior.gl.NoSum, 81, 14)
posterior.gl.Sum.format <- cbind(lbf$gamma, normalize(10^lbf.gl.prior.format[,5]), posterior.gl.Sum)  

#sub = gl
#l=indephits(sub$lbfavunif,sub$chr,sub$pos)
#sub=sub[l==1,]
~~~
> dim(lbf.gl.prior.format)
[1]   81 2207
> dim(posterior.gl.Sum.format)
[1]   81 2207
~~~

##lbf.sub <- MeanAcrossSigmaas(lbf.bigmat[,l.lbfavunif==1], 81, 14)
##lbf.sub.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub, 1, function(x) { sum1 <- sum(x); if (sum1 == length(x)) { sum1 <- 1; }; return(sum1); })), lbf.sub)
##lbf.sub.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat[,l.lbfavunif==1], matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat[,l.lbfavunif==1]), ncol=ncol(lbf.bigmat[,l.lbfavunif==1]), byrow=FALSE), 81, 14)
lbf.sub.prior <- SumAcrossSigmaas.wPriorAvg(lbf.bigmat[,l.lbfavunif==1], matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat[,l.lbfavunif==1]), ncol=ncol(lbf.bigmat[,l.lbfavunif==1]), byrow=FALSE), 81, 14)
lbf.sub.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub.prior, 1, function(x) { sum1 <- sum(x); if (sum1 == length(x)) { sum1 <- 1; }; return(sum1); })), lbf.sub.prior)
lbf.sub.prior.adjusted <- SumAcrossSigmaas.wPriorAvg(lbf.bigmat[,l.lbfavunif==1], matrix(rep(c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); }))/nsigma,nsigma), nrow = nrow(lbf.bigmat[,l.lbfavunif==1]), ncol=ncol(lbf.bigmat[,l.lbfavunif==1]), byrow=FALSE), 81, 14)
lbf.sub.prior.adjusted.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub.prior.adjusted, 1, function(x) { sum1 <- sum(x); if (sum1 == length(x)) { sum1 <- 1; }; return(sum1); })), lbf.sub.prior.adjusted)
posterior.sub.NoSum <- posteriorprob(lbf.bigmat[,l.lbfavunif==1], normalize(rep(c(0,lbf$prior[-1]),nsigma)))
posterior.sub.Sum <- SumAcrossSigmaas.pp.NoMax(posterior.gl.NoSum, 81, 14)
posterior.sub.Sum.format <- cbind(lbf$gamma, normalize(10^lbf.sub.prior.format[,5]), posterior.sub.Sum)  


##lbf.gl.format.ChoongwonWrite <- as.data.frame(lbf.gl.format)
##names(lbf.gl.format.ChoongwonWrite) <- c("dHb", "OxHb", "Pulse", "LvBrth", "lbfTotal", as.character(gl[,1]))
##row.names(lbf.gl.format.ChoongwonWrite) <- apply(lbf.gl.format.ChoongwonWrite[,1:4], 1, function(x) { return(paste(x, collapse="_"))})
##lbf.gl.format.ChoongwonWrite <- lbf.gl.format.ChoongwonWrite[,5:ncol(lbf.gl.format.ChoongwonWrite)]
lbf.gl.prior.format.ChoongwonWrite <- as.data.frame(lbf.gl.prior.format)
names(lbf.gl.prior.format.ChoongwonWrite) <- c("dHb", "OxHb", "Pulse", "LvBrth", "lbfTotal", as.character(gl[,1]))
row.names(lbf.gl.prior.format.ChoongwonWrite) <- apply(lbf.gl.prior.format.ChoongwonWrite[,1:4], 1, function(x) { return(paste(x, collapse="_"))})
lbf.gl.prior.format.ChoongwonWrite <- lbf.gl.prior.format.ChoongwonWrite[,5:ncol(lbf.gl.prior.format.ChoongwonWrite)]
#lbf.gl.prior.format.ChoongwonWrite.Posterior <- cbind(lbf.gl.prior.format.ChoongwonWrite[,1:4], apply(10^lbf.gl.prior.format.ChoongwonWrite[,5:ncol(lbf.gl.prior.format.ChoongwonWrite)], 2, normalize))
lbf.gl.prior.adjusted.format.ChoongwonWrite <- as.data.frame(lbf.gl.prior.adjusted.format)
names(lbf.gl.prior.adjusted.format.ChoongwonWrite) <- c("dHb", "OxHb", "Pulse", "LvBrth", "lbfTotal", as.character(gl[,1]))
row.names(lbf.gl.prior.adjusted.format.ChoongwonWrite) <- apply(lbf.gl.prior.adjusted.format.ChoongwonWrite[,1:4], 1, function(x) { return(paste(x, collapse="_"))})
lbf.gl.prior.adjusted.format.ChoongwonWrite <- lbf.gl.prior.adjusted.format.ChoongwonWrite[,5:ncol(lbf.gl.prior.adjusted.format.ChoongwonWrite)]
#posterior.gl.prior <- SumAcrossSigmaas.pp.NoMax(posteriorprob(lbf.bigmat, normalize(rep(c(0,lbf$prior[-1]),nsigma))), 81, 14)
posterior.gl.Sum.format.ChoongwonWrite <- as.data.frame(posterior.gl.Sum.format)
names(posterior.gl.Sum.format.ChoongwonWrite) <- c("dHb", "OxHb", "Pulse", "LvBrth", "lbfTotalPosterior", as.character(gl[,1]))
row.names(posterior.gl.Sum.format.ChoongwonWrite) <- apply(posterior.gl.Sum.format.ChoongwonWrite[,1:4], 1, function(x) { return(paste(x, collapse="_"))})
posterior.gl.Sum.format.ChoongwonWrite <- posterior.gl.Sum.format.ChoongwonWrite[,5:ncol(posterior.gl.Sum.format.ChoongwonWrite)]


~~~
#> posterior.gl.prior <- SumAcrossSigmaas.pp.NoMax(posteriorprob(lbf.bigmat, normalize(rep(c(0,lbf$prior[-1]),nsigma))), 81, 14)
#> apply(posterior.gl.prior, 2, sum)[1:10]
# [1] 1 1 1 1 1 1 1 1 1 1
# > dim(posterior.gl.prior)
# [1]   81 2087
# > sum(apply(posterior.gl.prior, 2, sum))     
# [1] 2087
#> sum(apply(posterior.gl.NoSum, 2, sum))
#[1] 2087
#> sum(apply(posterior.gl.Sum, 2, sum))
#[1] 2087
~~~ 

##lbf.sub.format.ChoongwonWrite <- as.data.frame(lbf.sub.format)
##names(lbf.sub.format.ChoongwonWrite) <- c("dHb", "OxHb", "Pulse", "LvBrth", "lbfTotal", as.character(sub[,1]))
##row.names(lbf.sub.format.ChoongwonWrite) <- apply(lbf.sub.format.ChoongwonWrite[,1:4], 1, function(x) { return(paste(x, collapse="_"))})
##lbf.sub.format.ChoongwonWrite <- lbf.sub.format.ChoongwonWrite[,5:ncol(lbf.sub.format.ChoongwonWrite)]
#lbf.sub.prior.format.ChoongwonWrite <- as.data.frame(lbf.sub.prior.format)
#names(lbf.sub.prior.format.ChoongwonWrite) <- c("dHb", "OxHb", "Pulse", "LvBrth", "lbfTotal", as.character(sub[,1]))
#row.names(lbf.sub.prior.format.ChoongwonWrite) <- apply(lbf.sub.prior.format.ChoongwonWrite[,1:4], 1, function(x) { return(paste(x, collapse="_"))})
#lbf.sub.prior.format.ChoongwonWrite <- lbf.sub.prior.format.ChoongwonWrite[,5:ncol(lbf.sub.prior.format.ChoongwonWrite)]
#lbf.sub.prior.format.ChoongwonWrite.Posterior <- apply(10^lbf.sub.prior.format.ChoongwonWrite, 2, normalize)



##write.table(t(rbind(names(lbf.gl.format.ChoongwonWrite), apply(lbf.gl.format.ChoongwonWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.NotPruned.logBFs.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.gl.format.ChoongwonWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.gl.prior.format.ChoongwonWrite), apply(lbf.gl.prior.format.ChoongwonWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.NotPruned.logBFs.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.gl.prior.format.ChoongwonWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.gl.prior.adjusted.format.ChoongwonWrite), apply(lbf.gl.prior.adjusted.format.ChoongwonWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.NotPruned.logBFsAdjusted.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.gl.prior.adjusted.format.ChoongwonWrite))), quote=FALSE)
##write.table(t(rbind(names(lbf.sub.format.ChoongwonWrite), apply(lbf.sub.format.ChoongwonWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.Pruned.logBFs.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.sub.format.ChoongwonWrite))), quote=FALSE)
#write.table(t(rbind(names(lbf.sub.prior.format.ChoongwonWrite), apply(lbf.sub.prior.format.ChoongwonWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.Pruned.logBFs.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.sub.prior.format.ChoongwonWrite))), quote=FALSE)
write.table(t(rbind(names(posterior.gl.Sum.format.ChoongwonWrite), apply(posterior.gl.Sum.format.ChoongwonWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.NotPruned.Posteriors.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(posterior.gl.Sum.format.ChoongwonWrite))), quote=FALSE)

write.table(gl[,c(1:4,13,15,14,16:23)], "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalSNPsAnalyzed.txt", row.names=FALSE, quote=FALSE)


#20161010 CHECK_0 -- Prob: Unresolved issue here, why are stronger lbf.gl entries much lower in lbf.gl.prior? Shouldn't be so probably/likely, right?
~~~
> 10^lbf.gl.prior[1:5,1:5]
            [,1]       [,2]       [,3]       [,4]       [,5]
[1,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000
[2,] 0.020442982 0.04748552 0.04748552 0.05230575 0.04873808
[3,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000
[4,] 0.048015130 0.04633202 0.04633202 0.05457150 0.01998008
[5,] 0.008032528 0.01263842 0.01263842 0.01540204 0.00673513
> 10^lbf.gl[1:5,1:5]
           [,1]      [,2]      [,3]      [,4]        [,5]
[1,]  1.0000000 1.0000000 1.0000000 1.0000000   1.0000000
[2,] 19.6119359 0.7570558 0.7570558 0.8356388   0.7776777
[3,]  1.0000000 1.0000000 1.0000000 1.0000000   1.0000000
[4,]  0.7656448 0.7381430 0.7381430 0.8726094 357.7542061
[5,] 11.3981202 0.6015976 0.6015976 0.7373183 154.5960049
> dim(10^lbf.gl.prior + matrix(apply(10^lbf.gl, 2, mean) * (1-mean(normalize(c(0,lbf$prior[-1])))), ncol=ncol(lbf.gl), nrow=nrow(lbf.gl), byrow=TRUE))
> normalize(c(0,lbf$prior[-1])))[1:5]
Error: unexpected ')' in "normalize(c(0,lbf$prior[-1])))"
> normalize(c(0,lbf$prior[-1]))[1:5]
[1] 0.00000000 0.06250000 0.00000000 0.06250000 0.02083333
> 19.6119359 * .0625
[1] 1.225746
> normalize(rep(c(0,lbf$prior[-1]),nsigma))[1:5]
[1] 0.000000000 0.004464286 0.000000000 0.004464286 0.001488095
> .0625/14
[1] 0.004464286
.
.
.
> log10(50) - log10(25)
[1] 0.30103
> log(50)
[1] 3.912023
> log10(50)
[1] 1.69897
> log10(25)
[1] 1.39794
> (10^log10(50) - log10(25))
[1] 48.60206
> 10^(log10(50) - log10(25))
[1] 2
> log10(.5)
[1] -0.30103
> log10(2)
[1] 0.30103
> log10(50) * log10(2)
[1] 0.5114409
> 10^(log10(50) * log10(2))
[1] 3.246691
> 10^(log10(50) - log10(2))
[1] 25
.
.
.
> mean(c(0,lbf$prior[-1]))
[1] 0.00617284
> .03125/.00617284
[1] 5.0625
> centered.lbf
centered.lbf
> centered.lbf
function(lbf){
  return(t(t(lbf) - apply(lbf,2,max)))
}
> c(c(0,lbf$prior[-1])-max(c(0,lbf$prior[-1]))+1)[1:5]
[1] 0.9687500 1.0000000 0.9687500 1.0000000 0.9791667
> c(0,lbf$prior[-1])[1:5]
[1] 0.00000000 0.03125000 0.00000000 0.03125000 0.01041667
> c(c(0,lbf$prior[-1])+(1-0.00617284))[1:5]
[1] 0.9938272 1.0250772 0.9938272 1.0250772 1.0042438
> c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); })[1:5]
+ 
> c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); }))[1:5]
[1] 0.000000 1.025077 0.000000 1.025077 1.004244
> quantile(c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); })))    
       0%       25%       50%       75%      100% 
0.0000000 0.9972994 0.9990355 1.0016397 1.0250772 
> mean(c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); })))
[1] 0.8036885
> normalize(rep(c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); })), 10))[1:5]
[1] 0.000000000 0.001574649 0.000000000 0.001574649 0.001542647
> rep(c(c(0,lbf$prior[-1])+(1-0.00617284)*sapply(c(0,lbf$prior[-1]), function(x) { if (x > 0) { x <- 1; } ; return(x); }))/10,10)[1:5]
[1] 0.0000000 0.1025077 0.0000000 0.1025077 0.1004244
~~~

posterior.sub.NoSum.marginals <- marginal.postprobs(posterior.sub.NoSum, lbf$gamma, length(sigmaa))
posterior.sub.NoSum.marginals.DirAssoc <- t(posterior.sub.NoSum.marginals[[2]])
posterior.sub.NoSum.marginals.DirAssoc <- data.frame(cbind(as.character(sub.lbfavunif$snp), posterior.sub.NoSum.marginals.DirAssoc))
posterior.sub.NoSum.marginals.DirAssoc[,2] <- as.numeric(as.character(posterior.sub.NoSum.marginals.DirAssoc[,2]))
posterior.sub.NoSum.marginals.DirAssoc[,3] <- as.numeric(as.character(posterior.sub.NoSum.marginals.DirAssoc[,3]))
posterior.sub.NoSum.marginals.DirAssoc[,4] <- as.numeric(as.character(posterior.sub.NoSum.marginals.DirAssoc[,4]))
posterior.sub.NoSum.marginals.DirAssoc[,5] <- as.numeric(as.character(posterior.sub.NoSum.marginals.DirAssoc[,5]))
colnames(posterior.sub.NoSum.marginals.DirAssoc) <- c("snp", "dHb", "OxHb", "Pulse", "LvBrth")
posterior.sub.NoSum.marginals.DirAssoc <- posterior.sub.NoSum.marginals.DirAssoc[order(sub.lbfavunif$lbfavunif, decreasing=TRUE),]
posterior.sub.NoSum.marginals.DirAssoc$snp <- factor(posterior.sub.NoSum.marginals.DirAssoc$snp, levels=sub.lbfavunif[order(sub.lbfavunif$lbfavunif),][,1])










#jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MultiPhenoGroups_Vs_lBFs.jpeg", width=2000, height=2000, res=300)
#
#plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(4,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n")
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#
#dev.off()
#
#jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MultiPhenoGroups_Vs_lBFs.EachPheno.jpeg", width=4000, height=4000, res=300)
#par(mfrow=c(2,2))
#
##"dHb", "OxHb", "Pulse", "LvBrth"
##sapply(lbf.gl.format[1:10,2], function(x) { val1 <- "YELLOW"; if (x > 0) { val1 <- "RED"; } else { val1 <- "BLACK" }; return(val1);})
#
#plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,1], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("topleft", c("dHb DirAssoc", "dHb IndirAssoc", "dHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,2], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("topleft", c("OxHb DirAssoc", "OxHb IndirAssoc", "OxHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,3], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("topleft", c("Pulse DirAssoc", "Pulse IndirAssoc", "Pulse NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,4], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("topleft", c("LvBrth DirAssoc", "LvBrth IndirAssoc", "LvBrth NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#dev.off()
#
#
##plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,1], function(x) { val1 <- "YELLOW"; if (x > 0) { val1 <- "RED"; } else { val1 <- "BLACK" }; return(val1);}))
##axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
##legend("topleft", c("dHb Assoc", "dHb NotAssoc"), col=c("RED","BLACK"), pch=c(1,1))
##plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,2], function(x) { val1 <- "YELLOW"; if (x > 0) { val1 <- "RED"; } else { val1 <- "BLACK" }; return(val1);}))
##axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
##legend("topleft", c("OxHb Assoc", "OxHb NotAssoc"), col=c("RED","BLACK"), pch=c(1,1))
##plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,3], function(x) { val1 <- "YELLOW"; if (x > 0) { val1 <- "RED"; } else { val1 <- "BLACK" }; return(val1);}))
##axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
##legend("topleft", c("Pulse Assoc", "Pulse NotAssoc"), col=c("RED","BLACK"), pch=c(1,1))
##plot(apply(lbf.gl.format[,1:4]>0, 1, sum), lbf.gl.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.format[lbf.gl.format[,5]>0,5])-.1,max(lbf.gl.format[,5])+.1), main="Distribution of lbfTotals in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "lbfTotal", xaxt="n", col=sapply(lbf.gl.format[,4], function(x) { val1 <- "YELLOW"; if (x > 0) { val1 <- "RED"; } else { val1 <- "BLACK" }; return(val1);}))
##axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
##
#
#
#jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MultiPhenoGroups_Vs_lBFs.wPrior.EachPheno.jpeg", width=4000, height=4000, res=300)
#par(mfrow=c(2,2))
#
#plot(apply(lbf.gl.prior.format[,1:4]>0, 1, sum), lbf.gl.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.prior.format[lbf.gl.prior.format[,5]>0,5])-.1,max(lbf.gl.prior.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.gl.prior.format[,1], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("bottomleft", c("dHb DirAssoc", "dHb IndirAssoc", "dHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#plot(apply(lbf.gl.prior.format[,1:4]>0, 1, sum), lbf.gl.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.prior.format[lbf.gl.prior.format[,5]>0,5])-.1,max(lbf.gl.prior.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.gl.prior.format[,2], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("bottomleft", c("OxHb DirAssoc", "OxHb IndirAssoc", "OxHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#plot(apply(lbf.gl.prior.format[,1:4]>0, 1, sum), lbf.gl.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.prior.format[lbf.gl.prior.format[,5]>0,5])-.1,max(lbf.gl.prior.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.gl.prior.format[,3], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("bottomleft", c("Pulse DirAssoc", "Pulse IndirAssoc", "Pulse NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#plot(apply(lbf.gl.prior.format[,1:4]>0, 1, sum), lbf.gl.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.gl.prior.format[lbf.gl.prior.format[,5]>0,5])-.1,max(lbf.gl.prior.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.gl.prior.format[,4], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
#axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
#legend("bottomleft", c("LvBrth DirAssoc", "LvBrth IndirAssoc", "LvBrth NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16))
#
#dev.off()


jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MultiPhenoGroups_Vs_logBFs.wPrior.Pruned.EachPheno.jpeg", width=4000, height=4000, res=300)
par(mfrow=c(2,2))

plot(apply(lbf.sub.prior.format[,1:4]>0, 1, sum), lbf.sub.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.format[lbf.sub.prior.format[,5]>0,5])-.1,max(lbf.sub.prior.format[,5])+.1), main="Distribution of per model prior*log10(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model prior*log10(BF)", xaxt="n", col=sapply(lbf.sub.prior.format[,1], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("dHb DirAssoc", "dHb IndirAssoc", "dHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(lbf.sub.prior.format[,1:4]>0, 1, sum), lbf.sub.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.format[lbf.sub.prior.format[,5]>0,5])-.1,max(lbf.sub.prior.format[,5])+.1), main="Distribution of per model prior*log10(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model prior*log10(BF)", xaxt="n", col=sapply(lbf.sub.prior.format[,2], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("OxHb DirAssoc", "OxHb IndirAssoc", "OxHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(lbf.sub.prior.format[,1:4]>0, 1, sum), lbf.sub.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.format[lbf.sub.prior.format[,5]>0,5])-.1,max(lbf.sub.prior.format[,5])+.1), main="Distribution of per model prior*log10(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model prior*log10(BF)", xaxt="n", col=sapply(lbf.sub.prior.format[,3], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("Pulse DirAssoc", "Pulse IndirAssoc", "Pulse NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(lbf.sub.prior.format[,1:4]>0, 1, sum), lbf.sub.prior.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.format[lbf.sub.prior.format[,5]>0,5])-.1,max(lbf.sub.prior.format[,5])+.1), main="Distribution of per model prior*log10(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model prior*log10(BF)", xaxt="n", col=sapply(lbf.sub.prior.format[,4], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("LvBrth DirAssoc", "LvBrth IndirAssoc", "LvBrth NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

dev.off()

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MultiPhenoGroups_Vs_logBFs.wPrior.Pruned.Adjusted.EachPheno.jpeg", width=4000, height=4000, res=300)
par(mfrow=c(2,2))

plot(apply(lbf.sub.prior.adjusted.format[,1:4]>0, 1, sum), lbf.sub.prior.adjusted.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.adjusted.format[lbf.sub.prior.adjusted.format[,5]>0,5])-.1,max(lbf.sub.prior.adjusted.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.sub.prior.adjusted.format[,1], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("dHb DirAssoc", "dHb IndirAssoc", "dHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(lbf.sub.prior.adjusted.format[,1:4]>0, 1, sum), lbf.sub.prior.adjusted.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.adjusted.format[lbf.sub.prior.adjusted.format[,5]>0,5])-.1,max(lbf.sub.prior.adjusted.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.sub.prior.adjusted.format[,2], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("OxHb DirAssoc", "OxHb IndirAssoc", "OxHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(lbf.sub.prior.adjusted.format[,1:4]>0, 1, sum), lbf.sub.prior.adjusted.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.adjusted.format[lbf.sub.prior.adjusted.format[,5]>0,5])-.1,max(lbf.sub.prior.adjusted.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.sub.prior.adjusted.format[,3], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("Pulse DirAssoc", "Pulse IndirAssoc", "Pulse NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(lbf.sub.prior.adjusted.format[,1:4]>0, 1, sum), lbf.sub.prior.adjusted.format[,5], xlim=c(1,4), ylim=c(min(lbf.sub.prior.adjusted.format[lbf.sub.prior.adjusted.format[,5]>0,5])-.1,max(lbf.sub.prior.adjusted.format[,5])+.1), main="Distribution of per model log(BF) in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model log(BF)", xaxt="n", col=sapply(lbf.sub.prior.adjusted.format[,4], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("LvBrth DirAssoc", "LvBrth IndirAssoc", "LvBrth NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

dev.off()

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MultiPhenoGroups_Vs_Posteriors.wPrior.Pruned.EachPheno.jpeg", width=4000, height=4000, res=300)
par(mfrow=c(2,2))

plot(apply(posterior.sub.Sum.format[,1:4]>0, 1, sum), posterior.sub.Sum.format[,5], xlim=c(1,4), ylim=c(min(posterior.sub.Sum.format[posterior.sub.Sum.format[,5]>0,5])-.1,max(posterior.sub.Sum.format[,5])+.1), main="Distribution of per model avg Posterior in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model avg Posterior", xaxt="n", col=sapply(posterior.sub.Sum.format[,1], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("dHb DirAssoc", "dHb IndirAssoc", "dHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(posterior.sub.Sum.format[,1:4]>0, 1, sum), posterior.sub.Sum.format[,5], xlim=c(1,4), ylim=c(min(posterior.sub.Sum.format[posterior.sub.Sum.format[,5]>0,5])-.1,max(posterior.sub.Sum.format[,5])+.1), main="Distribution of per model avg Posterior in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model avg Posterior", xaxt="n", col=sapply(posterior.sub.Sum.format[,2], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("OxHb DirAssoc", "OxHb IndirAssoc", "OxHb NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(posterior.sub.Sum.format[,1:4]>0, 1, sum), posterior.sub.Sum.format[,5], xlim=c(1,4), ylim=c(min(posterior.sub.Sum.format[posterior.sub.Sum.format[,5]>0,5])-.1,max(posterior.sub.Sum.format[,5])+.1), main="Distribution of per model avg Posterior in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model avg Posterior", xaxt="n", col=sapply(posterior.sub.Sum.format[,3], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("Pulse DirAssoc", "Pulse IndirAssoc", "Pulse NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

plot(apply(posterior.sub.Sum.format[,1:4]>0, 1, sum), posterior.sub.Sum.format[,5], xlim=c(1,4), ylim=c(min(posterior.sub.Sum.format[posterior.sub.Sum.format[,5]>0,5])-.1,max(posterior.sub.Sum.format[,5])+.1), main="Distribution of per model avg Posterior in each tier of multivariate assocation", xlab = "Number of phenotypes associated in model", ylab= "per model avg Posterior", xaxt="n", col=sapply(posterior.sub.Sum.format[,4], function(x) { val1 <- "YELLOW"; if (x == 1) { val1 <- "RED"; } else if (x == 2) { val1 <- "BLUE" } else { val1 <- "BLACK" }; return(val1);}), pch=16)
axis(1, at=c(1,2,3,4), xlab = "Number of phenotypes associated in model")
legend("bottomleft", c("LvBrth DirAssoc", "LvBrth IndirAssoc", "LvBrth NotAssoc"), col=c("RED","BLUE","BLACK"), pch=c(16,16,16), bg="transparent")

dev.off()




jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.DirAssocMarginalPosteriors.HeatPlot.hclust.vs1.jpeg", width=6000, height=2000, res=300)

ggPlotData <- melt(posterior.sub.NoSum.marginals.DirAssoc[hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,])
ggPlotData$snp <- factor(ggPlotData$snp, levels=posterior.sub.NoSum.marginals.DirAssoc[,1][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order])

ggplot(ggPlotData, aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue") + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.DirAssocMarginalPosteriors.HeatPlot.lbfavunif.vs1.jpeg", width=6000, height=2000, res=300)

ggPlotData <- melt(posterior.sub.NoSum.marginals.DirAssoc)
ggPlotData$snp <- factor(ggPlotData$snp, levels=sub.lbfavunif[order(sub.lbfavunif$lbfavunif),][,1]) 

ggplot(ggPlotData, aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue") + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()








#20161028
#Anna followup requests -- blue/red heatmaps showing direction of effect (w/ posterior probs); include indirect marginal posteriorprobs; give list of rsIDs in same order as heatmaps (do all three for each phenotype)

Z.Directions <- apply(Z, c(1,2), function(x) { y <- 0; if (x > 0) { y <- 1; } else if (x < 0) { y <- -1; } else if (x == 0) { y <- 0; } else { y <- NA; }; return(y);})

#Z.Directions.sub.lbfavunif <- Z.Directions[l.lbfavunif==1,]
Z.Directions.sub.lbfavunif <- Z.Directions[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),]

posterior.sub.NoSum.marginals.DirAssoc.wDir <- posterior.sub.NoSum.marginals.DirAssoc
posterior.sub.NoSum.marginals.DirAssoc.wDir[,2:5] <- posterior.sub.NoSum.marginals.DirAssoc[,2:5] * Z.Directions.sub.lbfavunif

#jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.DirAssocMarginalPosteriors.HeatPlot.wDir.hclust.vs1.jpeg", width=6000, height=2000, res=300)
jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.DirAssocMarginalPosteriors.HeatPlot.wDir.hclust_onDirVals.vs1.jpeg", width=6000, height=2000, res=300)

ggPlotData <- melt(posterior.sub.NoSum.marginals.DirAssoc.wDir[hclust(dist(posterior.sub.NoSum.marginals.DirAssoc.wDir[,2:5]))$order,])
#ggPlotData$snp <- factor(ggPlotData$snp, levels=posterior.sub.NoSum.marginals.DirAssoc[,1][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order])
ggPlotData$snp <- factor(ggPlotData$snp, levels=posterior.sub.NoSum.marginals.DirAssoc.wDir[,1][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc.wDir[,2:5]))$order])

#ggplot(ggPlotData, aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "darkred", high="steelblue") + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(ggPlotData, aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

posterior.sub.NoSum.marginals.IndirAssoc <- t(posterior.sub.NoSum.marginals[[3]])
posterior.sub.NoSum.marginals.IndirAssoc <- data.frame(cbind(as.character(sub.lbfavunif$snp), posterior.sub.NoSum.marginals.IndirAssoc))
posterior.sub.NoSum.marginals.IndirAssoc[,2] <- as.numeric(as.character(posterior.sub.NoSum.marginals.IndirAssoc[,2]))
posterior.sub.NoSum.marginals.IndirAssoc[,3] <- as.numeric(as.character(posterior.sub.NoSum.marginals.IndirAssoc[,3]))
posterior.sub.NoSum.marginals.IndirAssoc[,4] <- as.numeric(as.character(posterior.sub.NoSum.marginals.IndirAssoc[,4]))
posterior.sub.NoSum.marginals.IndirAssoc[,5] <- as.numeric(as.character(posterior.sub.NoSum.marginals.IndirAssoc[,5]))
colnames(posterior.sub.NoSum.marginals.DirAssoc) <- c("snp", "dHb", "OxHb", "Pulse", "LvBrth")
posterior.sub.NoSum.marginals.IndirAssoc <- posterior.sub.NoSum.marginals.IndirAssoc[order(sub.lbfavunif$lbfavunif, decreasing=TRUE),]
posterior.sub.NoSum.marginals.IndirAssoc$snp <- factor(posterior.sub.NoSum.marginals.IndirAssoc$snp, levels=sub.lbfavunif[order(sub.lbfavunif$lbfavunif),][,1])

~~~
> head(cbind(as.character(sub[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,]$snp), Z[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,]))
     [,1]         [,2]                 [,3]                 [,4]               
[1,] "rs73521166" "-0.344158088083837" "-0.450229177170152" "4.27465312595573" 
[2,] "rs28666472" "0.40062342015514"   "0.811619148991553"  "4.60637458376269" 
[3,] "rs7075577"  "0.125018000251463"  "0.301858149833601"  "4.4401879635215"  
[4,] "rs73003243" "0.573446815280319"  "0.30531611964425"   "-3.98566728056487"
[5,] "rs77110410" "-0.576465622303058" "-0.108922197920176" "3.97834701511612" 
[6,] "rs17162012" "0.179994609604668"  "-0.445409289773187" "3.90005592712398" 
     [,5]                 
[1,] "0.623550237361774"  
[2,] "-0.0222462807662346"
[3,] "0.271694279753069"  
[4,] "1.19297720646042"   
[5,] "-0.95646016621093"  
[6,] "-0.108089622460911" 
> head(posterior.sub.NoSum.marginals.DirAssoc.wDir.wIndir[hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,])
           snp      dHb_D     OxHb_D    Pulse_D   LvBrth_D   LvBrth_I
242 rs73521166 -0.2075156 -0.2044559  0.9963811  0.2207132  0.1838201
152 rs28666472  0.1943010  0.2126776  0.9978930 -0.1933157 -0.2578040
229  rs7075577  0.1899353  0.1876191  0.9969269  0.1977447  0.2241711
239 rs73003243  0.2711983  0.2449013 -0.9822726  0.2944129  0.2621156
252 rs77110410 -0.2826538 -0.2502683  0.9743237 -0.2804821 -0.2573470
99  rs17162012  0.2456704 -0.2582733  0.9660319 -0.2442024 -0.2396297
> head(sub[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,])
            snp chr       pos   maf     p_OxHb n_OxHb      p_dHb n_dHb
1772 rs73521166  19   3832972 0.452 -0.6525452    900 -0.7307274   900
1083 rs28666472  18  66176774 0.397  0.4170102    887  0.6886974   887
1678  rs7075577  10  14781053 0.380 -0.7627602    902 -0.9005093   902
1742 rs73003243   3 150052369 0.237 -0.7601254    915 -0.5663422   915
1895 rs77110410   6  37161228 0.153  0.9132642    913  0.5643005   913
788  rs17162012   7 140946310 0.147  0.6560240    904 -0.8571568   904
           p_Pulse n_Pulse   p_LvBrth n_LvBrth annot     Z.OxHb      Z.dHb
1772  1.914351e-05     899  0.5329230      959     0 -0.4502292 -0.3441581
1083  4.097500e-06     886 -0.9822515      944     0  0.8116191  0.4006234
1678 -8.988034e-06     901 -0.7858571      961     0  0.3018581  0.1250180
1742  6.729075e-05     914 -0.2328783      975     0  0.3053161  0.5734468
1895 -6.939604e-05     912  0.3388398      972     0 -0.1089222 -0.5764656
788  -9.617047e-05     903  0.9139246      963     0 -0.4454093  0.1799946
       Z.Pulse    Z.LvBrth   mvstat      mvp     unip nmin      lbfav lbfavunif
1772  4.274653  0.62355024 20.93026 3.963188 4.717978  899 -0.4666497  2.076461
1083  4.606375 -0.02224628 21.74438 4.132371 5.387481  886 -0.4434408  2.271945
1678  4.440188  0.27169428 20.55194 3.884661 5.046335  901 -0.4939709  2.107829
1742 -3.985667  1.19297721 18.23058 3.404266 4.172045  914 -0.2906451  1.413396
1895  3.978347 -0.95646017 17.63984 3.282452 4.158665  912 -0.3046780  1.255730
788   3.900056 -0.10808962 16.00023 2.945440 4.016958  903 -0.3709700  1.108375
~~~

posterior.sub.NoSum.marginals.IndirAssoc.wDir <- posterior.sub.NoSum.marginals.IndirAssoc
posterior.sub.NoSum.marginals.IndirAssoc.wDir[,2:5] <- posterior.sub.NoSum.marginals.IndirAssoc[,2:5] * Z.Directions.sub.lbfavunif

posterior.sub.NoSum.marginals.DirAssoc.wIndir <- cbind(posterior.sub.NoSum.marginals.DirAssoc, posterior.sub.NoSum.marginals.IndirAssoc[,5])
colnames(posterior.sub.NoSum.marginals.DirAssoc.wIndir) <- c("snp", "dHb_D", "OxHb_D", "Pulse_D", "LvBrth_D", "LvBrth_I")
posterior.sub.NoSum.marginals.DirAssoc.wDir.wIndir <- cbind(posterior.sub.NoSum.marginals.DirAssoc.wDir, posterior.sub.NoSum.marginals.IndirAssoc.wDir[,5])
colnames(posterior.sub.NoSum.marginals.DirAssoc.wDir.wIndir) <- c("snp", "dHb_D", "OxHb_D", "Pulse_D", "LvBrth_D", "LvBrth_I")


#jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalPosteriors.HeatPlot.hclust.vs1.jpeg", width=6000, height=2000, res=300)
#jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalPosteriors.HeatPlot.wDir.hclust.vs1.jpeg", width=6000, height=2000, res=300)
#jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalPosteriors.HeatPlot.wDir.hclust_onDirVals.vs1.jpeg", width=6000, height=2000, res=300)
jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalPosteriors.HeatPlot.hclust_onDirVals.vs1.jpeg", width=6000, height=2000, res=300)

ggPlotData <- melt(posterior.sub.NoSum.marginals.DirAssoc.wIndir[hclust(dist(posterior.sub.NoSum.marginals.DirAssoc.wIndir[,2:6]))$order,])
#ggPlotData <- melt(posterior.sub.NoSum.marginals.DirAssoc.wDir.wIndir[hclust(dist(posterior.sub.NoSum.marginals.DirAssoc.wDir.wIndir[,2:6]))$order,])
#ggPlotData$snp <- factor(ggPlotData$snp, levels=posterior.sub.NoSum.marginals.DirAssoc[,1][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order])
ggPlotData$snp <- factor(ggPlotData$snp, levels=posterior.sub.NoSum.marginals.DirAssoc.wDir[,1][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc.wDir[,2:5]))$order])

ggplot(ggPlotData, aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue") + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggplot(ggPlotData, aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradientn(colours=c("darkred", "white", "steelblue")) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()


write.table(file="/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalSNPs.Pruned.Posteriors.vs1.txt", posterior.sub.NoSum.marginals.DirAssoc.wDir.wIndir[hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,], col.names=TRUE, quote=FALSE,sep= " ", row.names=FALSE)

#write.table(file="/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalSNPs.Pruned.ZScores.vs1.txt", cbind(as.character(sub[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),]$snp), Z[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),]), col.names=c("snp", "Sat", "Hb", "Pulse", "LvBrth"), quote=FALSE,sep= " ", row.names=FALSE)
write.table(file="/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs3/Choongwon2016.RanHalfFlip.PhenoGroup2.MarginalSNPs.Pruned.ZScores.vs1.txt", cbind(as.character(sub[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,]$snp), Z[l.lbfavunif==1,][order(sub.lbfavunif$lbfavunif, decreasing=TRUE),][hclust(dist(posterior.sub.NoSum.marginals.DirAssoc[,2:5]))$order,]), col.names=c("snp", "Sat", "Hb", "Pulse", "LvBrth"), quote=FALSE,sep= " ", row.names=FALSE)







posterior.gl.NoSum.marginals <- marginal.postprobs(posterior.gl.NoSum, lbf$gamma, length(sigmaa))
posterior.gl.NoSum.marginals.IndirAssoc <- t(posterior.gl.NoSum.marginals[[3]])
posterior.gl.NoSum.marginals.IndirAssoc <- data.frame(cbind(as.character(gl$snp), posterior.gl.NoSum.marginals.IndirAssoc))
posterior.gl.NoSum.marginals.IndirAssoc[,2] <- as.numeric(as.character(posterior.gl.NoSum.marginals.IndirAssoc[,2]))
posterior.gl.NoSum.marginals.IndirAssoc[,3] <- as.numeric(as.character(posterior.gl.NoSum.marginals.IndirAssoc[,3]))
posterior.gl.NoSum.marginals.IndirAssoc[,4] <- as.numeric(as.character(posterior.gl.NoSum.marginals.IndirAssoc[,4]))
posterior.gl.NoSum.marginals.IndirAssoc[,5] <- as.numeric(as.character(posterior.gl.NoSum.marginals.IndirAssoc[,5]))
colnames(posterior.gl.NoSum.marginals.IndirAssoc) <- c("snp", "Sat", "Hb", "Pulse", "LvBrth")
posterior.gl.NoSum.marginals.IndirAssoc <- posterior.gl.NoSum.marginals.IndirAssoc[order(gl.lbfavunif$lbfavunif, decreasing=TRUE),]
posterior.gl.NoSum.marginals.IndirAssoc$snp <- factor(posterior.gl.NoSum.marginals.IndirAssoc$snp, levels=gl.lbfavunif[order(gl.lbfavunif$lbfavunif),][,1])














#20160810 -- All of the below is old from the first Vs1, just keeping here for now
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.newtophits.vs3.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)
#write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.newtophits.vs3.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>5 & sub$nmin>100 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>5 & sub$nmin>100 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)

#newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>100,c(1:3,7,18:20,25)]

#newhits = gl[gl$annot==0 & gl$lbfav>3.323736 & gl$nmin>100,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)










