set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.3Phenos_2.RSS0.vs1.SignCrrct.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.3Phenos_2.dtlesslesssignif.vs1.SignCrrct.vs1.annot.vs1.MAF.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0,]

Z = cbind(gl$Z.dHb,gl$Z.OxHb,gl$Z.Pulse)

#~~~
> dim(gl)
[1] 3736   17
> gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
> dim(gl)
[1] 3736   17
> VYY
          Z.dHb       Z.OxHb     Z.Pulse
[1,]  1.0000000 -0.380453596 0.117220010
[2,] -0.3804536  1.000000000 0.008529997
[3,]  0.1172200  0.008529997 1.000000000
> head(gl)
         snp chr       pos   maf       p_OxHb n_OxHb         p_dHb n_dHb
1 rs10006686   4 177993314 0.213 7.461186e-05    915 -0.6300364000   915
2 rs10007496   4  17406405 0.255 2.964302e-03    920  0.0258232700   920
3 rs10009224   4 133257015 0.330 2.163151e-01    889  0.0002333501   889
4 rs10011618   4 177999800 0.211 9.727020e-05    916 -0.5423901000   916
5 rs10012073   4 147015936 0.130 1.818857e-01    921  0.0003085042   921
6 rs10012313   4  17405513 0.256 3.289038e-03    921  0.0267067900   921
     p_Pulse n_Pulse annot    Z.OxHb      Z.dHb    Z.Pulse   mvstat      mvp
1 0.35177540     914     0  3.961079 -0.4816756  0.9311511 17.49612 3.252847
2 0.72615460     919     0 -2.971416 -2.2288593 -0.3502453 22.04681 4.195284
3 0.02956022     888     0 -1.236386 -3.6798730 -2.1759337 24.21073 4.646410
4 0.37561710     915     0  3.897302 -0.6092027  0.8860008 16.62290 3.073227
5 0.03942898     920     0 -1.334971 -3.6080524 -2.0596817 23.78433 4.557389
6 0.73987210     920     0 -2.939348 -2.2157768 -0.3320227 21.66394 4.115643
      unip
1 4.127192
2 2.528078
3 3.631992
4 4.012020
5 3.510739
6 2.482931
~~~

n = cbind(gl$n_dHb, gl$n_OxHb, gl$n_Pulse)
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
> dim(gl.glhits)
[1]  1 18
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

~~~
> order(ebprior.glhits.collapse,decreasing=TRUE)
 [1] 22  4 13 23  5 14  6 24 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
[26] 26 27
> order(ebprior.glhits.collapse2,decreasing=TRUE)
 [1] 22  4 13 23  5 14 24  6 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
[26] 26 27
> order(ebprior.glhits.collapse3,decreasing=TRUE)
 [1] 22  4 13 23  5 14  6 24 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
[26] 26 27
> order(ebprior.glhits.collapse4,decreasing=TRUE)
 [1] 22  4 13 23  5 14 24  6 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
[26] 26 27
~~~

pp.glhits = posteriorprob(lbf.glhits,ebprior.glhits) #posterior prob on models for gl hits
pp.glhits.collapse =  apply(pp.glhits,2,collapse, nsigmaa=length(sigmaa))

# this returns a list with elements  pU pD and pI
marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

#looking at which models are favored by the prior
cumsum(sort(ebprior.glhits.collapse,decreasing=TRUE))
lbf$gamma[order(ebprior.glhits.collapse,decreasing=TRUE),]
modelmatrix = cbind(lbf$gamma,ebprior.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
modelmatrix = data.frame(cbind(modelmatrix,cumsum(modelmatrix[,4])))
colnames(modelmatrix)= c("dHb","OxHb","Pulse","p","cump")

allassoc=(apply((modelmatrix[,1:3]>0),1,sum)==3) #vector of which represent situations in which all 3 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:3]>0),1,sum)==2) #vector of which represent situations in which 2 phenotypes are associated

sum(modelmatrix[allassoc,4]) #0.9574025
sum(modelmatrix[allbut1assoc,4]) #0.03937951

#~~~
> sum(modelmatrix[allassoc,4]) #0.9574025
[1] 4.266743e-47
> sum(modelmatrix[allbut1assoc,4]) #0.03937951
[1] 0.9279949
#~~~

#look at weight on each of Pulse, OxHb, dHb being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4]) 
# 2.327247e-48, 0, 0.9279949 

#~~~
> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4])
[1] 2.327247e-48
> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4])
[1] 0
> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4])
[1] 0.9279949
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notPulse, or notOxHb or notdHb)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(matrix(ppmatrix)[allassoc,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,3]==0,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,2]==0,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"Choongwon2016.3Phenos_2.gl.bestmodel.vs2.SignCrrct.vs1.txt",quote=F,sep=" & ",row.names=F)

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
#> min(lbf.av.glhits)
#[1] 4.346141
#> min(lbf.av.origprior.glhits)
#[1] 3.006386
> min(lbf.av.glhits)
[1] 7.15689
> min(lbf.av.origprior.glhits)
[1] 5.446706
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

#> seq.int(from = 1, by=27, length.out=14)
# [1]    1   82  163  244  325  406  487  568  649  730  811  892  973 1054

prior is column of all values across gammas (with multiple entries per sigmaa) so per row it's the same priors for each SNP

> dim(lbf.newhits)
[1] 378 193
> dim(lbf.bigmat)
[1]  378 1484


lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 27, 14) 
lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)[order(log10(apply(10^lbf.gl, 1, sum))),]
lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 27, 14) 
lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)[order(log10(apply(10^lbf.gl.prior, 1, sum))),]

lbf.sub <- MeanAcrossSigmaas(lbf.bigmat[,l==1], 27, 14)
lbf.sub.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub, 1, sum)), lbf.sub)[order(log10(apply(10^lbf.sub, 1, sum))),]
lbf.sub.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat[,l==1], matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat[,l==1]), ncol=ncol(lbf.bigmat[,l==1]), byrow=FALSE), 27, 14) 
lbf.sub.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub.prior, 1, sum)), lbf.sub.prior)[order(log10(apply(10^lbf.sub.prior, 1, sum))),]

gl[l==1,1]

sub = gl
l=indephits(sub$lbfavunif,sub$chr,sub$pos)
sub=sub[l==1,]
#lbf.newhits= lbf.bigmat
#lbf.newhits= lbf.newhits[,l==1]
newhits.unif = sub[sub$annot==0 & sub$lbfavunif>3.006386 & sub$nmin>100,c(1:4,11:14,19:20)]

#lbf.newhits.sigmaaMean <- MeanAcrossSigmaas(lbf.newhits, 81, 14)
#lbf.newhits.sigmaaMean.plusGammaDescrips <- cbind(lbf$gamma, lbf.newhits.sigmaaMean)

names(lbf.gl.prior.format) <- c("dHb", "OxHb", "Pulse", "lbfTotal_perModel", as.character(gl[,1]))
blah <- as.data.frame(lbf.gl.prior.format)
names(blah) <- c("dHb", "OxHb", "Pulse", "lbfTotal", as.character(gl[,1]))
row.names(blah) <- apply(blah[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
blah <- rbind(blah, log10(apply(10^blah, 2, sum)))
row.names(blah) <- c(row.names(blah), "lbfTotal_perSNP")

apply(blah[1:3,1:3], 1, function(x) { return(paste(x, collapse="_"))})

write.table(t(blah)[4:nrow(t(blah)),], "ouasdouhasd", row.names=TRUE, quote=FALSE)

write.table(lbf.gl.prior.format, "ouasdouhasd")
write.table(blah, "ouasdouhasd", row.names=FALSE, quote=FALSE)

~~~
> sub[order(sub$lbfavunif, decreasing=TRUE),][1:10,]
> sub[order(sub$lbfav, decreasing=TRUE),][1:10,]
> VYY
> apply(10^lbf.gl.prior, 2, normalize)[1:10,1:10]
> apply(apply(10^lbf.gl.prior, 2, normalize), 2, sum)[1:10]
> lbf.gl.prior.format.ChoongownWrite.Posterior[1:10,1:10]
~~~




newhits.unif = sub[sub$annot==0 & sub$lbfavnif>3.006386 & sub$nmin>100,c(1:4,11:14,19:20)]

lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 27, 14)
#lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)[order(log10(apply(10^lbf.gl, 1, sum))),]
lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)
lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 27, 14)
#lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)[order(log10(apply(10^lbf.gl.prior, 1, sum))),]
lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)

sub = gl
l=indephits(sub$lbfavunif,sub$chr,sub$pos)
sub=sub[l==1,]

lbf.sub <- MeanAcrossSigmaas(lbf.bigmat[,l==1], 27, 14)
#lbf.sub.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub, 1, sum)), lbf.sub)[order(log10(apply(10^lbf.sub, 1, sum))),]
lbf.sub.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub, 1, sum)), lbf.sub)
lbf.sub.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat[,l==1], matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat[,l==1]), ncol=ncol(lbf.bigmat[,l==1]), byrow=FALSE), 27, 14)
#lbf.sub.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub.prior, 1, sum)), lbf.sub.prior)[order(log10(apply(10^lbf.sub.prior, 1, sum))),]
lbf.sub.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub.prior, 1, sum)), lbf.sub.prior)

#blah <- as.data.frame(lbf.gl.prior.format)
#names(blah) <- c("dHb", "OxHb", "Pulse", "lbfTotal", as.character(gl[,1]))
#row.names(blah) <- apply(blah[,1:3], 1, function(x) { return(paste(x, collapse="_"))})

lbf.gl.format.ChoongownWrite <- as.data.frame(lbf.gl.format)
names(lbf.gl.format.ChoongownWrite) <- c("dHb", "OxHb", "Pulse", "lbfTotal", as.character(gl[,1]))
row.names(lbf.gl.format.ChoongownWrite) <- apply(lbf.gl.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.gl.format.ChoongownWrite <- lbf.gl.format.ChoongownWrite[,5:ncol(lbf.gl.format.ChoongownWrite)]
lbf.gl.prior.format.ChoongownWrite <- as.data.frame(lbf.gl.prior.format)
names(lbf.gl.prior.format.ChoongownWrite) <- c("dHb", "OxHb", "Pulse", "lbfTotal", as.character(gl[,1]))
row.names(lbf.gl.prior.format.ChoongownWrite) <- apply(lbf.gl.prior.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.gl.prior.format.ChoongownWrite <- lbf.gl.prior.format.ChoongownWrite[,5:ncol(lbf.gl.prior.format.ChoongownWrite)]
lbf.gl.prior.format.ChoongownWrite.Posterior <- apply(10^lbf.gl.prior.format.ChoongownWrite, 2, normalize)

lbf.sub.format.ChoongownWrite <- as.data.frame(lbf.sub.format)
names(lbf.sub.format.ChoongownWrite) <- c("dHb", "OxHb", "Pulse", "lbfTotal", as.character(sub[,1]))
row.names(lbf.sub.format.ChoongownWrite) <- apply(lbf.sub.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.sub.format.ChoongownWrite <- lbf.sub.format.ChoongownWrite[,5:ncol(lbf.sub.format.ChoongownWrite)]
lbf.sub.prior.format.ChoongownWrite <- as.data.frame(lbf.sub.prior.format)
names(lbf.sub.prior.format.ChoongownWrite) <- c("dHb", "OxHb", "Pulse", "lbfTotal", as.character(sub[,1]))
row.names(lbf.sub.prior.format.ChoongownWrite) <- apply(lbf.sub.prior.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.sub.prior.format.ChoongownWrite <- lbf.sub.prior.format.ChoongownWrite[,5:ncol(lbf.sub.prior.format.ChoongownWrite)]
lbf.sub.prior.format.ChoongownWrite.Posterior <- apply(10^lbf.sub.prior.format.ChoongownWrite, 2, normalize)



#write.table(t(lbf.gl.format.ChoongownWrite), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup2.logBFs.txt", row.names=TRUE, quote=FALSE)
#write.table(t(apply(lbf.gl.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) })), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup2.NotPruned.logBFs.txt", row.names=TRUE, col.names=as.character(c("rsID", row.names(lbf.gl.format.ChoongownWrite))), quote=FALSE)
#write.table(t(apply(lbf.gl.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) })), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup2.NotPruned.logBFs.txt", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(t(rbind(names(lbf.gl.format.ChoongownWrite), apply(lbf.gl.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup2.NotPruned.logBFs.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.gl.format.ChoongownWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.gl.prior.format.ChoongownWrite), apply(lbf.gl.prior.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup2.NotPruned.logBFs.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.gl.prior.format.ChoongownWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.sub.format.ChoongownWrite), apply(lbf.sub.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup2.Pruned.logBFs.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.sub.format.ChoongownWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.sub.prior.format.ChoongownWrite), apply(lbf.sub.prior.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup2.Pruned.logBFs.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.sub.prior.format.ChoongownWrite))), quote=FALSE)





#write.table(t(blah)[4:nrow(t(blah)),], "ouasdouhasd", row.names=TRUE, quote=FALSE)


~~~
> blah[1:10,1:10]
      dHb OxHb Pulse lbfTotal rs10042098 rs10069542  rs1008264 rs10086302
0_0_0   0  0     0     -Inf       -Inf       -Inf       -Inf       -Inf
2_0_0   2  0     0     -Inf       -Inf       -Inf       -Inf       -Inf
0_2_0   0  2     0     -Inf       -Inf       -Inf       -Inf       -Inf
2_2_0   2  2     0     -Inf       -Inf       -Inf       -Inf       -Inf
0_0_2   0  0     2     -Inf       -Inf       -Inf       -Inf       -Inf
2_0_2   2  0     2     -Inf       -Inf       -Inf       -Inf       -Inf
0_2_2   0  2     2     -Inf       -Inf       -Inf       -Inf       -Inf
2_2_2   2  2     2     -Inf       -Inf       -Inf       -Inf       -Inf
0_2_1   0  2     1 1.753260 -0.4781951 -0.4781951 -1.0258958  -2.845679
2_0_1   2  0     1 1.771887 -0.5709561 -0.5709561 -0.8767494  -2.844661
      rs10086401 rs10086462
0_0_0       -Inf       -Inf
2_0_0       -Inf       -Inf
0_2_0       -Inf       -Inf
2_2_0       -Inf       -Inf
0_0_2       -Inf       -Inf
2_0_2       -Inf       -Inf
0_2_2       -Inf       -Inf
2_2_2       -Inf       -Inf
0_2_1  -2.808332  -2.779030
2_0_1  -2.762436  -2.715298
> as.character(gl[,1])[1:10]
 [1] "rs10042098" "rs10069542" "rs1008264"  "rs10086302" "rs10086401"
 [6] "rs10086462" "rs10086940" "rs10088178" "rs10089136" "rs10089402"
> gl[1:10,1]
 [1] rs10042098 rs10069542 rs1008264  rs10086302 rs10086401 rs10086462
 [7] rs10086940 rs10088178 rs10089136 rs10089402
1484 Levels: rs10042098 rs10069542 rs1008264 rs10086302 rs10086401 ... rs9988539
> as.character(sub[,1])[1:10]
 [1] "rs10086940" "rs10105487" "rs10112319" "rs10162740" "rs10184862"
 [6] "rs10246790" "rs10271491" "rs10271931" "rs10462914" "rs10473935"
> sub[1:10,1]
 [1] rs10086940 rs10105487 rs10112319 rs10162740 rs10184862 rs10246790
 [7] rs10271491 rs10271931 rs10462914 rs10473935
1484 Levels: rs10042098 rs10069542 rs1008264 rs10086302 rs10086401 ... rs9988539
> lbf.gl.format.ChoongownWrite[1:10,1:10]
      rs10042098 rs10069542   rs1008264 rs10086302  rs10086401  rs10086462
0_0_0  0.0000000  0.0000000  0.00000000  0.0000000  0.00000000  0.00000000
1_0_0 -0.1227195 -0.1227195 -0.08106369 -0.1286427  2.18653718  1.75208927
2_0_0  0.0000000  0.0000000  0.00000000  0.0000000  0.00000000  0.00000000
0_1_0 -0.1175143 -0.1175143 -0.08680725  2.5042353 -0.09677778 -0.10011317
1_1_0 -0.2108375 -0.2108375 -0.15342430  2.1416281  1.87684932  1.46039084
2_1_0 -0.1197784 -0.1197784 -0.09059113  2.4695742 -0.11639016 -0.11671216
0_2_0  0.0000000  0.0000000  0.00000000  0.0000000  0.00000000  0.00000000
1_2_0 -0.1249627 -0.1249627 -0.08488790 -0.1433005  2.13913968  1.71303431
2_2_0  0.0000000  0.0000000  0.00000000  0.0000000  0.00000000  0.00000000
0_0_1  2.1703030  2.1703030  1.71919198 -0.1430408 -0.10316689 -0.07242698
      rs10086940  rs10088178 rs10089136 rs10089402
0_0_0  0.0000000  0.00000000  0.0000000  0.0000000
1_0_0  1.7164710 -0.07910429 -0.1459305 -0.1286427
2_0_0  0.0000000  0.00000000  0.0000000  0.0000000
0_1_0 -0.1425457 -0.12958663  2.0655419  2.5042353
1_1_0  1.3836130 -0.18483234  1.7212912  2.1416281
2_1_0 -0.1453264 -0.12984813  2.0708063  2.4695742
0_2_0  0.0000000  0.00000000  0.0000000  0.0000000
1_2_0  1.7101185 -0.07939044 -0.1436614 -0.1433005
2_2_0  0.0000000  0.00000000  0.0000000  0.0000000
0_0_1  0.0565179  1.69516201 -0.1460475 -0.1430408
> apply(lbf.sub.prior.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) })[1:10,1:10]
       rs10086940  rs10105487 rs10112319  rs10162740  rs10184862 rs10246790
0_0_0  0.00000000  0.00000000  0.0000000  0.00000000  0.00000000  0.0000000
1_0_0 -2.24321919  0.01444402 -1.7825777 -2.14808633 -1.23875018 -2.1904665
2_0_0  0.00000000  0.00000000  0.0000000  0.00000000  0.00000000  0.0000000
0_1_0  0.01923957 -2.13468798 -2.2032743  0.09956698 -0.66644329 -0.1159674
1_1_0 -0.60660326 -0.29892854 -2.2348683 -0.44041858 -0.05224558 -0.6993633
2_1_0 -0.56050830 -2.67772050 -2.8049141 -0.58248603 -1.49064222 -0.7367191
0_2_0  0.00000000  0.00000000  0.0000000  0.00000000  0.00000000  0.0000000
1_2_0 -2.83564108 -0.47547724 -2.3839242 -2.78369733 -2.03893278 -2.7995253
2_2_0  0.00000000  0.00000000  0.0000000  0.00000000  0.00000000  0.0000000
0_0_1 -1.76430132 -2.23536644 -0.0498365 -2.18707599 -2.14616706 -2.0666245
      rs10271491 rs10271931 rs10462914 rs10473935
0_0_0  0.0000000  0.0000000   0.000000  0.0000000
1_0_0 -0.6525289 -0.6525289  -2.116168 -2.2249534
2_0_0  0.0000000  0.0000000   0.000000  0.0000000
0_1_0 -2.1581901 -2.1581901  -2.234301 -1.1852284
1_1_0 -1.0871908 -1.0871908  -2.536995 -1.7361437
2_1_0 -2.7427615 -2.7427615  -2.833527 -1.7875045
0_2_0  0.0000000  0.0000000   0.000000  0.0000000
1_2_0 -1.2085229 -1.2085229  -2.714809 -2.8271139
2_2_0  0.0000000  0.0000000   0.000000  0.0000000
0_0_1 -2.1057183 -2.1057183   0.318316  0.1149702
> dim(lbf.gl.prior.format.ChoongownWrite)
[1]   27 1484
> dim(lbf.sub.prior.format.ChoongownWrite)
[1]  27 196
> head(subunif[order(subunif$lbfavunif, decreasing=TRUE),])
             snp chr       pos   maf       p_OxHb n_OxHb         p_dHb n_dHb
2653   rs6539167  12 105172975 0.375 4.925449e-03    908  9.752924e-06   908
2012 rs372272284   2  46584859 0.248 5.713546e-09    914 -8.825847e-01   914
1159   rs1499702   6  66628474 0.394 1.202528e-02    920  4.052976e-05   920
1391   rs1938121   6  66626824 0.394 1.202528e-02    920  4.052976e-05   920
1392   rs1938122   6  66626689 0.394 1.202528e-02    920  4.052976e-05   920
3549   rs9363491   6  66625418 0.394 1.202528e-02    920  4.052976e-05   920
       p_Pulse n_Pulse annot    Z.OxHb      Z.dHb     Z.Pulse   mvstat      mvp
2653 0.9413777     907     0 -2.811869 -4.4225800 -0.07353838 43.71141 8.759960
2012 0.7307703     913     1  5.824933 -0.1476935  0.34410104 38.93231 7.746195
1159 0.6366713     919     0  2.511402  4.1044388  0.47235797 36.31259 7.191709
1391 0.6366713     919     0  2.511402  4.1044388  0.47235797 36.31259 7.191709
1392 0.6366713     919     0  2.511402  4.1044388  0.47235797 36.31259 7.191709
3549 0.6366713     919     0  2.511402  4.1044388  0.47235797 36.31259 7.191709
         unip nmin    lbfav lbfavunif
2653 5.010865  907 4.297146  5.927270
2012 8.243094  913 7.156890  5.446706
1159 4.392226  919 3.418779  4.559091
1391 4.392226  919 3.418779  4.559091
1392 4.392226  919 3.418779  4.559091
3549 4.392226  919 3.418779  4.559091
> head(subtrain[order(subtrain$lbfav, decreasing=TRUE),])
             snp chr       pos   maf       p_OxHb n_OxHb       p_dHb n_dHb
2012 rs372272284   2  46584859 0.248 5.713546e-09    914 -0.88258470   914
1871   rs2897724  21  21674336 0.382 2.795057e-06    907  0.17461770   907
2784    rs688874  11 130084507 0.358 2.638686e-05    915  0.06055327   915
2293   rs4944119  11  76428055 0.383 1.008179e-05    911  0.26821750   911
3203    rs768095   8 113450779 0.485 1.678790e-05    893  0.21505100   893
2665   rs6567111  18  57323790 0.177 2.288221e-05    921  0.09592008   921
        p_Pulse n_Pulse annot    Z.OxHb      Z.dHb     Z.Pulse   mvstat
2012  0.7307703     913     1  5.824933 -0.1476935  0.34410104 38.93231
1871  0.1511684     906     0 -4.685333 -1.3575148 -1.43541666 34.40536
2784  0.2674815     914     0  4.202598  1.8767434  1.10888105 32.13269
2293 -0.5887631     910     0  4.415412  1.1071768 -0.54062933 29.53693
3203  0.8144962     893     0 -4.303819 -1.2397956 -0.23462964 28.24590
2665 -0.9834672     920     0 -4.234734 -1.6649633  0.02072227 30.75540
          mvp     unip nmin    lbfav lbfavunif
2012 7.746195 8.243094  913 7.156890  5.446706
1871 6.788661 5.553609  906 5.822011  4.307006
2784 6.309176 4.578612  914 5.141330  3.812437
2293 5.762720 4.996462  910 4.965660  3.592176
3203 5.491469 4.775004  893 4.832804  3.411828
2665 6.019061 4.640502  920 4.829406  3.424386
~~~




> min(lbf.av.glhits)
[1] 7.15689 
> min(lbf.av.origprior.glhits)
[1] 5.446706    


sub = gl
lunif=indephits(sub$lbfavunif,sub$chr,sub$pos)
subunif=sub[lunif==1,]
#lbf.newhits= lbf.bigmat
#lbf.newhits= lbf.newhits[,l==1]
newhits.unif = subunif[subunif$annot==0 & subunif$lbfavunif>5.446706 & subunif$nmin>100,c(1:4,11:14,19:20)]

~~~
> dim(subunif)
[1] 443  20
> dim(subtrain)
[1] 417  20
> newhits.unif
           snp chr       pos   maf annot    Z.OxHb    Z.dHb     Z.Pulse
2653 rs6539167  12 105172975 0.375     0 -2.811869 -4.42258 -0.07353838
        lbfav lbfavunif
2653 4.297146   5.92727
> newhits
 [1] snp       chr       pos       maf       annot     Z.OxHb    Z.dHb    
 [8] Z.Pulse   lbfav     lbfavunif
<0 rows> (or 0-length row.names)
~~~

#sub = gl[gl$nmin>50000,]
sub = gl
ltrain=indephits(sub$lbfav,sub$chr,sub$pos)
subtrain=sub[ltrain==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
#newhits = sub[sub$annot==0 & sub$lbfav>4.346141 & sub$nmin>100,c(1:3,7,18:20,25)]
newhits = subtrain[subtrain$annot==0 & subtrain$lbfav>7.15689 & subtrain$nmin>100,c(1:4,11:14,19:20)]












#> newhits
#           snp chr       pos   maf annot      Z.OxHb    Z.dHb     Z.Pulse   lbfav
#1045 rs6539167  12 105172975 0.375     0 -4.811546 3.317785 -0.07353838 4.80128
> newhits
           snp chr       pos   maf annot      Z.OxHb    Z.dHb     Z.Pulse
1045 rs6539167  12 105172975 0.375     0 -4.811546 3.317785 -0.07353838
        lbfav lbfavflat
1045 5.391825  3.736339











#~~~
> dim(sub)
[1] 1484   19i
> l=indephits(sub$lbfav,sub$chr,sub$pos)
> sub=sub[l==1,]
> dim(sub)
[1] 193  19
#~~~

#extract lbfs for all the new hits
#lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.bigmat
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>4.346141 & sub$nmin>100]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notPulse, or notOxHb or notdHb)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))


#### in results
#~~~
> cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))
                        bestclass                    
 [1,] "chr10:124064829" "1"       "0.967765408967158"
 [2,] "chr11:46723937"  "4"       "0.549865354645701"
 [3,] "rs10196674"      "1"       "0.999469526408363"
 [4,] "rs10783573"      "1"       "0.95424438170892" 
 [5,] "rs10946458"      "1"       "0.865187718309875"
 [6,] "rs11024028"      "1"       "0.914947935771833"
 [7,] "rs13046645"      "1"       "0.883951671613376"
 [8,] "rs13204965"      "1"       "0.999202725076594"
 [9,] "rs35481967"      "1"       "0.521502423223891"
[10,] "rs4671925"       "1"       "0.517051826742596"
[11,] "rs4776341"       "1"       "0.998602764452439"
[12,] "rs6500348"       "1"       "0.919130784452271"
[13,] "rs7209460"       "1"       "0.880655565505357"
[14,] "rs7560852"       "1"       "0.999970495459465"

#~~~

pdf("plots.bychr.vs2.pdf")
for(i in 1:22){
 if(sum(sub$chr==i)>0){
plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>1000))
abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>4.608818],col=3)
}
}
dev.off()

#try to plot all hits to see relationships?
tophits =  sub[sub$lbfav>3.323736 & sub$nmin>100,]
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
write.table(file="Choongwon2016.3Phenos_2.newtophits.vs1.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)
#write.table(file="Choongwon2016.3Phenos_2.newtophits.vs1.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>5 & sub$nmin>100 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>5 & sub$nmin>100 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)

#newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>100,c(1:3,7,18:20,25)]

#newhits = gl[gl$annot==0 & gl$lbfav>3.323736 & gl$nmin>100,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)










