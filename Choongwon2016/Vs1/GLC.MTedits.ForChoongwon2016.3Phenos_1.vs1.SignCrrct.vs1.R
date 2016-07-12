set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.3Phenos_1.RSS0.vs1.SignCrrct.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.3Phenos_1.dtlesslesssignif.vs1.SignCrrct.vs1.annot.vs1.MAF.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0,]

Z = cbind(gl$Z.Sat,gl$Z.Hb,gl$Z.Pulse)

#~~~
> dim(gl)
[1] 1484   17
> gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
> 
> dim(gl)
[1] 1484   17
> head(gl)
> head(gl)
         snp chr       pos   maf         p_Hb n_Hb         p_Sat n_Sat
1 rs10042098   5 176143265 0.342 4.324555e-01  912 -9.952110e-01   912
2 rs10069542   5 176143314 0.342 4.324555e-01  912 -9.952110e-01   912
3  rs1008264   7  25618485 0.192 7.521300e-01  919  2.877942e-01   919
4 rs10086302   8 113467607 0.486 1.760304e-05  891  8.806490e-01   891
5 rs10086401   8  91572818 0.240 7.862087e-01  901  2.314939e-05   901
6 rs10086462   8  91558089 0.240 8.203675e-01  891  8.043717e-05   891
        p_Pulse n_Pulse annot       Z.Hb        Z.Sat    Z.Pulse   mvstat
1  3.635102e-05     911     0 -0.7849967  0.006002157 -4.1295271 17.51117
2  3.635102e-05     911     0 -0.7849967  0.006002157 -4.1295271 17.51117
3 -8.035984e-05     918     0  0.3158320  1.062972991 -3.9433244 16.38799
4  8.787754e-01     891     0 -4.2933082 -0.150146543 -0.1525218 18.65599
5 -3.524451e-01     900     0 -0.2712371 -4.232124789  0.9298571 18.52878
6 -2.268448e-01     890     0 -0.2270723 -3.943093819  1.2085265 16.51506
       mvp     unip
1 3.255945 4.439483
2 3.255945 4.439483
3 3.024988 4.094961
4 3.492103 4.754412
5 3.465827 4.635460
6 3.051077 4.094543
~~~

n = cbind(gl$n_Sat, gl$n_Hb, gl$n_Pulse)
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
 [1] 23  5 22  4 24  6 14 13 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
[26] 26 27
> order(ebprior.glhits.collapse2,decreasing=TRUE)
 [1] 23  5 22  4 24  6 14 13 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
[26] 26 27
> order(ebprior.glhits.collapse3,decreasing=TRUE)
 [1] 23  5 22  4 24  6 14 13 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
[26] 26 27
> order(ebprior.glhits.collapse4,decreasing=TRUE)
 [1] 23  5 22  4 24  6 14 13 15  1  2  3  7  8  9 10 11 12 16 17 18 19 20 21 25
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
colnames(modelmatrix)= c("Sat","Hb","Pulse","p","cump")

allassoc=(apply((modelmatrix[,1:3]>0),1,sum)==3) #vector of which represent situations in which all 3 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:3]>0),1,sum)==2) #vector of which represent situations in which 2 phenotypes are associated

sum(modelmatrix[allassoc,4]) #0.9574025
sum(modelmatrix[allbut1assoc,4]) #0.03937951

#~~~
> sum(modelmatrix[allassoc,4]) #0.9574025
[1] 0.8467565
> sum(modelmatrix[allbut1assoc,4]) #0.03937951
[1] 0.1532435
#~~~

#look at weight on each of Pulse, Hb, Sat being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4]) 
# 0.07841946, 0, 0.7579819 

#~~~
> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4])
[1] 0.1532408
> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4])
[1] 0
> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4])
[1] 2.693174e-06
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notPulse, or notHb or notSat)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(matrix(ppmatrix)[allassoc,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,3]==0,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,2]==0,]),colSums(matrix(ppmatrix)[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"Choongwon2016.3Phenos_1.gl.bestmodel.vs2.SignCrrct.vs1.txt",quote=F,sep=" & ",row.names=F)

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
[1] 4.990048
> min(lbf.av.origprior.glhits)
[1] 3.568456
#~~~


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]

#lbf.av.all.flatprior = lbf.av(lbf.bigmat, rep(lbf$prior,nsigma))
lbf.av.all.flatprior = lbf.av(lbf.bigmat, normalize(rep(c(0,lbf$prior[-1]),nsigma))) 
gl$lbfavflat = lbf.av.all.flatprior






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
l=indephits(sub$lbfavflat,sub$chr,sub$pos)
sub=sub[l==1,]
#lbf.newhits= lbf.bigmat
#lbf.newhits= lbf.newhits[,l==1]
newhits.flat = sub[sub$annot==0 & sub$lbfavflat>3.006386 & sub$nmin>100,c(1:4,11:14,19:20)]

#lbf.newhits.sigmaaMean <- MeanAcrossSigmaas(lbf.newhits, 81, 14)
#lbf.newhits.sigmaaMean.plusGammaDescrips <- cbind(lbf$gamma, lbf.newhits.sigmaaMean)

names(lbf.gl.prior.format) <- c("Sat", "Hb", "Pulse", "lbfTotal_perModel", as.character(gl[,1]))
blah <- as.data.frame(lbf.gl.prior.format)
names(blah) <- c("Sat", "Hb", "Pulse", "lbfTotal", as.character(gl[,1]))
row.names(blah) <- apply(blah[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
blah <- rbind(blah, log10(apply(10^blah, 2, sum)))
row.names(blah) <- c(row.names(blah), "lbfTotal_perSNP")

apply(blah[1:3,1:3], 1, function(x) { return(paste(x, collapse="_"))})

write.table(t(blah)[4:nrow(t(blah)),], "ouasdouhasd", row.names=TRUE, quote=FALSE)

write.table(lbf.gl.prior.format, "ouasdouhasd")
write.table(blah, "ouasdouhasd", row.names=FALSE, quote=FALSE)

~~~
> sub[order(sub$lbfavflat, decreasing=TRUE),][1:10,]
             snp chr       pos   maf         p_Hb n_Hb         p_Sat n_Sat
1045   rs6539167  12 105172975 0.375 1.497669e-06  908 -9.073437e-04   908
775  rs372272284   2  46584859 0.248 2.952695e-07  914  9.929316e-02   914
724    rs2897724  21  21674336 0.382 1.383738e-06  907  9.569834e-01   907
1097    rs688874  11 130084507 0.358 2.929710e-06  915 -6.531881e-01   915
711   rs28666472  18  66176774 0.397 4.021709e-01  887  9.851052e-01   887
744   rs34055428   4 125978367 0.101 9.983088e-01  909  1.281524e-06   909
758   rs35179589   5  13069765 0.437 8.087318e-01  885 -4.738930e-06   885
515    rs1840812  19  28774567 0.314 1.754192e-01  914 -5.077551e-06   914
622    rs2404529   2 127936624 0.484 8.060646e-01  921 -1.309655e-01   921
1160   rs7275795  21  40663670 0.326 7.639523e-01  889 -1.194650e-05   889
           p_Pulse n_Pulse annot         Z.Hb       Z.Sat     Z.Pulse   mvstat
1045  9.413777e-01     907     0 -4.811546472  3.31778477 -0.07353838 32.32984
775   7.307703e-01     913     1  5.126441656  1.64829008  0.34410104 30.41378
724   1.511684e-01     906     0 -4.827332767 -0.05393946 -1.43541666 24.78446
1097  2.674815e-01     914     0  4.675687118 -0.44933765  1.10888105 22.48053
711   4.097500e-06     886     0  0.837750401  0.01866895  4.60637458 21.77198
744   3.946114e-01     908     0 -0.002119606 -4.84259799 -0.85128438 25.37145
758   7.491955e-01     884     0 -0.242062416  4.57602541 -0.31970034 20.97023
515   4.687430e-01     913     0  1.354994843 -4.56155774  0.72452579 21.90038
622  -2.178049e-05     920     0 -0.245506081  1.51030546  4.24581025 22.04049
1160 -2.633668e-01     888     0  0.300294799 -4.37856192 -1.11846886 21.67488
          mvp     unip nmin      lbfav lbfavflat
1045 6.350735 5.824584  907  4.8012799  3.736339
775  5.947159 6.529781  913  4.9497605  3.568456
724  4.766279 5.858946  906  4.0813368  2.821698
1097 4.285573 5.533175  914  3.7213755  2.497416
711  4.138112 5.387481  886 -0.4498697  2.436257
744  4.889024 5.892273  908  2.8495889  2.424619
758  3.971487 5.324320  884  2.5084894  2.379899
515  4.164819 5.294346  913  2.6634983  2.272652
622  4.193971 4.661932  920 -0.5018111  2.265494
1160 4.117917 4.922759  888  2.2620380  2.252230
> sub[order(sub$lbfav, decreasing=TRUE),][1:10,]
             snp chr       pos   maf         p_Hb n_Hb         p_Sat n_Sat
775  rs372272284   2  46584859 0.248 2.952695e-07  914  0.0992931600   914
1045   rs6539167  12 105172975 0.375 1.497669e-06  908 -0.0009073437   908
724    rs2897724  21  21674336 0.382 1.383738e-06  907  0.9569834000   907
1097    rs688874  11 130084507 0.358 2.929710e-06  915 -0.6531881000   915
1051   rs6567111  18  57323790 0.177 6.872462e-06  921 -0.8580935000   921
893    rs4944119  11  76428055 0.383 1.104393e-05  911  0.8528308000   911
989   rs60685370   5 171910564 0.179 7.527685e-06  903 -0.0356554600   903
624     rs240569  11  64733678 0.169 8.070854e-06  919 -0.7224649000   919
1270    rs768095   8 113450779 0.485 1.590864e-05  893 -0.9667614000   893
503   rs17245616   5  34997038 0.166 9.575504e-05  911 -0.0008820181   911
         p_Pulse n_Pulse annot      Z.Hb       Z.Sat      Z.Pulse   mvstat
775   0.73077030     913     1  5.126442  1.64829008  0.344101041 30.41378
1045  0.94137770     907     0 -4.811546  3.31778477 -0.073538385 32.32984
724   0.15116840     906     0 -4.827333 -0.05393946 -1.435416665 24.78446
1097  0.26748150     914     0  4.675687 -0.44933765  1.108881049 22.48053
1051 -0.98346720     920     0 -4.497601  0.17880159  0.020722275 20.35319
893  -0.58876310     910     0  4.395655  0.18550775 -0.540629326 20.22017
989  -0.99379840     902     0  4.478195 -2.10083492 -0.007772631 23.47723
624   0.06335743     918     0 -4.463296  0.35516645 -1.856674834 22.33860
1270  0.81449620     893     0 -4.315716  0.04167046 -0.234629643 18.70042
503   0.32145070     910     0  3.901104 -3.32568417  0.991481134 24.70901
          mvp     unip nmin    lbfav lbfavflat
775  5.947159 6.529781  913 4.949761  3.568456
1045 6.350735 5.824584  907 4.801280  3.736339
724  4.766279 5.858946  906 4.081337  2.821698
1097 4.285573 5.533175  914 3.721376  2.497416
1051 3.843432 5.162888  920 3.316069  2.043388
893  3.815849 4.956876  910 3.301370  2.189904
989  4.493312 5.123339  902 3.269011  2.144649
624  4.256019 5.093081  918 3.217221  2.019685
1270 3.501283 4.798367  893 3.133022  1.997087
503  4.750511 4.018838  910 3.121285  2.057169
> VYY
           Z.Sat        Z.Hb     Z.Pulse
[1,]  1.00000000 -0.07230977 -0.10345962
[2,] -0.07230977  1.00000000  0.06952607
[3,] -0.10345962  0.06952607  1.00000000
> apply(10^lbf.gl.prior, 2, normalize)[1:10,1:10]
              [,1]         [,2]         [,3]         [,4]        [,5]
 [1,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
 [2,] 0.0015347725 0.0015347725 0.0031866404 0.3823175889 0.001116057
 [3,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
 [4,] 0.0011869729 0.0011869729 0.0042946639 0.0007763802 0.349312872
 [5,] 0.0005835984 0.0005835984 0.0016641515 0.0822858707 0.081531122
 [6,] 0.0002986758 0.0002986758 0.0010749804 0.0001880807 0.080670387
 [7,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
 [8,] 0.0003866812 0.0003866812 0.0007974443 0.0886939911 0.000270248
 [9,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
[10,] 0.3682141630 0.3682141630 0.3093336457 0.0007967158 0.001366796
              [,6]         [,7]         [,8]         [,9]        [,10]
 [1,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
 [2,] 0.0024930609 0.0012126885 0.0024680383 0.3656340984 0.3823175889
 [3,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
 [4,] 0.3183861503 0.2219259158 0.0034676634 0.0018712516 0.0007763802
 [5,] 0.0767860659 0.0525248948 0.0013868273 0.0828386603 0.0822858707
 [6,] 0.0753379510 0.0584063609 0.0008889474 0.0004721977 0.0001880807
 [7,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
 [8,] 0.0006093872 0.0003099755 0.0006296020 0.0934159994 0.0886939911
 [9,] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
[10,] 0.0035317634 0.0036531469 0.3247848345 0.0020363324 0.0007967158
> apply(apply(10^lbf.gl.prior, 2, normalize), 2, sum)[1:10]
 [1] 1 1 1 1 1 1 1 1 1 1
> lbf.gl.prior.format.ChoongownWrite.Posterior[1:10,1:10]
        rs10042098   rs10069542    rs1008264   rs10086302  rs10086401
0_0_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
1_0_0 0.0015347725 0.0015347725 0.0031866404 0.3823175889 0.001116057
2_0_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
0_1_0 0.0011869729 0.0011869729 0.0042946639 0.0007763802 0.349312872
1_1_0 0.0005835984 0.0005835984 0.0016641515 0.0822858707 0.081531122
2_1_0 0.0002986758 0.0002986758 0.0010749804 0.0001880807 0.080670387
0_2_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
1_2_0 0.0003866812 0.0003866812 0.0007974443 0.0886939911 0.000270248
2_2_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.000000000
0_0_1 0.3682141630 0.3682141630 0.3093336457 0.0007967158 0.001366796
        rs10086462   rs10086940   rs10088178   rs10089136   rs10089402
0_0_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
1_0_0 0.0024930609 0.0012126885 0.0024680383 0.3656340984 0.3823175889
2_0_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
0_1_0 0.3183861503 0.2219259158 0.0034676634 0.0018712516 0.0007763802
1_1_0 0.0767860659 0.0525248948 0.0013868273 0.0828386603 0.0822858707
2_1_0 0.0753379510 0.0584063609 0.0008889474 0.0004721977 0.0001880807
0_2_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
1_2_0 0.0006093872 0.0003099755 0.0006296020 0.0934159994 0.0886939911
2_2_0 0.0000000000 0.0000000000 0.0000000000 0.0000000000 0.0000000000
0_0_1 0.0035317634 0.0036531469 0.3247848345 0.0020363324 0.0007967158
~~~




newhits.flat = sub[sub$annot==0 & sub$lbfavflat>3.006386 & sub$nmin>100,c(1:4,11:14,19:20)]

lbf.gl <- MeanAcrossSigmaas(lbf.bigmat, 27, 14)
#lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)[order(log10(apply(10^lbf.gl, 1, sum))),]
lbf.gl.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl, 1, sum)), lbf.gl)
lbf.gl.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat, matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat), ncol=ncol(lbf.bigmat), byrow=FALSE), 27, 14)
#lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)[order(log10(apply(10^lbf.gl.prior, 1, sum))),]
lbf.gl.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.gl.prior, 1, sum)), lbf.gl.prior)

sub = gl
l=indephits(sub$lbfavflat,sub$chr,sub$pos)
sub=sub[l==1,]

lbf.sub <- MeanAcrossSigmaas(lbf.bigmat[,l==1], 27, 14)
#lbf.sub.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub, 1, sum)), lbf.sub)[order(log10(apply(10^lbf.sub, 1, sum))),]
lbf.sub.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub, 1, sum)), lbf.sub)
lbf.sub.prior <- MeanAcrossSigmaas.wPriorAvg(lbf.bigmat[,l==1], matrix(normalize(rep(c(0,lbf$prior[-1]),nsigma)), nrow = nrow(lbf.bigmat[,l==1]), ncol=ncol(lbf.bigmat[,l==1]), byrow=FALSE), 27, 14)
#lbf.sub.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub.prior, 1, sum)), lbf.sub.prior)[order(log10(apply(10^lbf.sub.prior, 1, sum))),]
lbf.sub.prior.format <- cbind(lbf$gamma, log10(apply(10^lbf.sub.prior, 1, sum)), lbf.sub.prior)

#blah <- as.data.frame(lbf.gl.prior.format)
#names(blah) <- c("Sat", "Hb", "Pulse", "lbfTotal", as.character(gl[,1]))
#row.names(blah) <- apply(blah[,1:3], 1, function(x) { return(paste(x, collapse="_"))})

lbf.gl.format.ChoongownWrite <- as.data.frame(lbf.gl.format)
names(lbf.gl.format.ChoongownWrite) <- c("Sat", "Hb", "Pulse", "lbfTotal", as.character(gl[,1]))
row.names(lbf.gl.format.ChoongownWrite) <- apply(lbf.gl.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.gl.format.ChoongownWrite <- lbf.gl.format.ChoongownWrite[,5:ncol(lbf.gl.format.ChoongownWrite)]
lbf.gl.prior.format.ChoongownWrite <- as.data.frame(lbf.gl.prior.format)
names(lbf.gl.prior.format.ChoongownWrite) <- c("Sat", "Hb", "Pulse", "lbfTotal", as.character(gl[,1]))
row.names(lbf.gl.prior.format.ChoongownWrite) <- apply(lbf.gl.prior.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.gl.prior.format.ChoongownWrite <- lbf.gl.prior.format.ChoongownWrite[,5:ncol(lbf.gl.prior.format.ChoongownWrite)]
lbf.gl.prior.format.ChoongownWrite.Posterior <- apply(10^lbf.gl.prior.format.ChoongownWrite, 2, normalize)

lbf.sub.format.ChoongownWrite <- as.data.frame(lbf.sub.format)
names(lbf.sub.format.ChoongownWrite) <- c("Sat", "Hb", "Pulse", "lbfTotal", as.character(sub[,1]))
row.names(lbf.sub.format.ChoongownWrite) <- apply(lbf.sub.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.sub.format.ChoongownWrite <- lbf.sub.format.ChoongownWrite[,5:ncol(lbf.sub.format.ChoongownWrite)]
lbf.sub.prior.format.ChoongownWrite <- as.data.frame(lbf.sub.prior.format)
names(lbf.sub.prior.format.ChoongownWrite) <- c("Sat", "Hb", "Pulse", "lbfTotal", as.character(sub[,1]))
row.names(lbf.sub.prior.format.ChoongownWrite) <- apply(lbf.sub.prior.format.ChoongownWrite[,1:3], 1, function(x) { return(paste(x, collapse="_"))})
lbf.sub.prior.format.ChoongownWrite <- lbf.sub.prior.format.ChoongownWrite[,5:ncol(lbf.sub.prior.format.ChoongownWrite)]
lbf.sub.prior.format.ChoongownWrite.Posterior <- apply(10^lbf.sub.prior.format.ChoongownWrite, 2, normalize)



#write.table(t(lbf.gl.format.ChoongownWrite), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup1.logBFs.txt", row.names=TRUE, quote=FALSE)
#write.table(t(apply(lbf.gl.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) })), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup1.NotPruned.logBFs.txt", row.names=TRUE, col.names=as.character(c("rsID", row.names(lbf.gl.format.ChoongownWrite))), quote=FALSE)
#write.table(t(apply(lbf.gl.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) })), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup1.NotPruned.logBFs.txt", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(t(rbind(names(lbf.gl.format.ChoongownWrite), apply(lbf.gl.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup1.NotPruned.logBFs.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.gl.format.ChoongownWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.gl.prior.format.ChoongownWrite), apply(lbf.gl.prior.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup1.NotPruned.logBFs.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.gl.prior.format.ChoongownWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.sub.format.ChoongownWrite), apply(lbf.sub.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup1.Pruned.logBFs.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.sub.format.ChoongownWrite))), quote=FALSE)
write.table(t(rbind(names(lbf.sub.prior.format.ChoongownWrite), apply(lbf.sub.prior.format.ChoongownWrite, c(1,2), function(x) { if (!is.finite(x)) { x <- 0; }; return(x) }))), "/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Choongwon2016/Vs1/Choongwon2016.PhenoGroup1.Pruned.logBFs.unifpriors.txt", row.names=FALSE, col.names=as.character(c("rsID", row.names(lbf.sub.prior.format.ChoongownWrite))), quote=FALSE)





#write.table(t(blah)[4:nrow(t(blah)),], "ouasdouhasd", row.names=TRUE, quote=FALSE)


~~~
> blah[1:10,1:10]
      Sat Hb Pulse lbfTotal rs10042098 rs10069542  rs1008264 rs10086302
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
~~~






























> min(lbf.av.glhits)
[1] 4.990048
> min(lbf.av.origprior.glhits)
[1] 3.568456    



sub = gl
l=indephits(sub$lbfavflat,sub$chr,sub$pos)
sub=sub[l==1,]
#lbf.newhits= lbf.bigmat
#lbf.newhits= lbf.newhits[,l==1]
newhits.flat = sub[sub$annot==0 & sub$lbfavflat>3.568456 & sub$nmin>100,c(1:4,11:14,19:20)]

~~~
> dim(sub)
[1] 197  20
> newhits.flat
           snp chr       pos   maf annot      Z.Hb    Z.Sat     Z.Pulse
1045 rs6539167  12 105172975 0.375     0 -4.811546 3.317785 -0.07353838
        lbfav lbfavflat
1045 5.391825  3.736339
~~~



#sub = gl[gl$nmin>50000,]
sub = gl
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
#newhits = sub[sub$annot==0 & sub$lbfav>4.346141 & sub$nmin>100,c(1:3,7,18:20,25)]
newhits = sub[sub$annot==0 & sub$lbfav>4.990048 & sub$nmin>100,c(1:4,11:14,19:20)]

#> newhits
#           snp chr       pos   maf annot      Z.Hb    Z.Sat     Z.Pulse   lbfav
#1045 rs6539167  12 105172975 0.375     0 -4.811546 3.317785 -0.07353838 4.80128
> newhits
           snp chr       pos   maf annot      Z.Hb    Z.Sat     Z.Pulse
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
#(all assoc, or notPulse, or notHb or notSat)
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
write.table(file="Choongwon2016.3Phenos_1.newtophits.vs1.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>4.346141 & sub$nmin>100 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)
#write.table(file="Choongwon2016.3Phenos_1.newtophits.vs1.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>5 & sub$nmin>100 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>5 & sub$nmin>100 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)

#newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>100,c(1:3,7,18:20,25)]

#newhits = gl[gl$annot==0 & gl$lbfav>3.323736 & gl$nmin>100,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)










