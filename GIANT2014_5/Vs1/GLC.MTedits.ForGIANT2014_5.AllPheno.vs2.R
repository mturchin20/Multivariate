set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.RSS0.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$chrbp) & gl$chrbp != "None",]
Z = cbind(gl$Z.height,gl$Z.BMI,gl$Z.WHRadjBMI,gl$Z.WHR,gl$Z.HIPadjBMI,gl$Z.HIP,gl$Z.WCadjBMI,gl$Z.WC)

#~~~
#> dim(gl)
#[1] 56322    43
#> gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$chrbp) & gl$chrbp != "None",]
#> dim(gl)
#[1] 56301    43
#> head(gl)
#snp chr       pos a1 a2       chrbp   maf annot beta_height se_height
#1  rs6531098   2  18457622  T  C  2_18457622 0.183     2      -0.019    0.0037
#2  rs6437061   2 232893295  A  C 2_232893295 0.433     2       0.034    0.0031
#3  rs6437062   2 232926991  A  C 2_232926991 0.333     2      -0.029    0.0035
#4 rs11978390   7  19615684  T  G  7_19615684 0.017     2      -0.058    0.0100
#5 rs16968231  15  74519634  T  C 15_74519634 0.075     2      -0.034    0.0059
#6 rs16968233  15  74520101  A  G 15_74520101 0.075     2       0.034    0.0059
#n_height beta_BMI se_BMI  n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI
#1   247989  -0.0001 0.0048 230302         0.0034       0.0055      140543
#2   250263  -0.0079 0.0031 321964         0.0037       0.0035      209385
#3   245891  -0.0116 0.0044 232538         0.0000       0.0051      141130
#4   231926   0.0181 0.0133 217020         0.0038       0.0160      124659
#5   253240  -0.0140 0.0075 234015        -0.0016       0.0085      142759
#6   253220  -0.0135 0.0075 234006        -0.0020       0.0085      142746
#beta_WHR se_WHR  n_WHR beta_HIPadjBMI se_HIPadjBMI n_HIPadjBMI beta_HIP
#1   0.0010 0.0054 142387         -0.011       0.0057      141597  -0.0100
#2  -0.0009 0.0034 212067          0.021       0.0036      210421   0.0072
#3  -0.0059 0.0051 144049          0.015       0.0054      142152  -0.0016
#4   0.0066 0.0150 128036         -0.058       0.0160      125273  -0.0280
#5  -0.0047 0.0085 144597          0.030       0.0090      143812   0.0110
#6  -0.0050 0.0085 144584          0.030       0.0090      143802   0.0120
#se_HIP  n_HIP beta_WCadjBMI se_WCadjBMI n_WCadjBMI beta_WC  se_WC   n_WC
#1 0.0056 143246        0.1833      0.0054     151263  0.1833 0.0054 151739
#2 0.0037 212906        0.4333      0.0035     229808  0.4333 0.0035 231963
#3 0.0053 144908        0.3333      0.0050     151849  0.3333 0.0052 153401
#4 0.0160 128457        0.0167      0.0150     134453  0.0167 0.0160 137855
#5 0.0090 145457        0.0750      0.0085     152083  0.0750 0.0086 153949
#6 0.0090 145439        0.0750      0.0085     153466  0.0750 0.0086 153936
#Z.height      Z.BMI Z.WHRadjBMI     Z.WHR Z.HIPadjBMI     Z.HIP Z.WCadjBMI
#1  5.099807 0.02080652   0.6128130 0.1891184    1.873495 1.8734955  0.6588377
#2 11.053855 2.53514746   1.0581216 0.2793190    5.843876 1.9773684  4.7746724
#3  8.186023 2.63636313   0.0000000 1.1503494    2.753288 0.3054808  2.1834865
#4  5.630201 1.36104332   0.2404260 0.4261480    3.533562 1.6953977  2.2414027
#5  5.831292 1.86665344   0.1891184 0.5533847    3.290527 1.2535654  1.9685917
#6  5.840642 1.80000404   0.2404260 0.5828415    3.332725 1.3105791  1.9514798
#Z.WC    mvstat       mvp      unip
#1 1.4395315  32.87966  5.898185  6.468521
#2 0.8064212 172.56951 35.532054 27.677781
#3 1.4757910  84.65738 16.746342 15.568636
#4 1.1030626  47.94740  9.014163  7.744727
#5 0.1763742  50.87584  9.625321  8.259637
#6 0.2018935  50.79996  9.609467  8.283997
#~~~

gl.priorHits <- gl[gl$annot==1,]
Z.priorHits <- Z[gl$annot==1,]

n = cbind(gl$n_height, gl$n_BMI, gl$n_WHRadjBMI, gl$n_WHR, gl$n_HIPadjBMI, gl$n_HIP, gl$n_WCadjBMI, gl$n_WC)
n = apply(n,1,min)
gl$nmin=n
sigmaa=c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)
#lbf=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa)
#lbf.bigmat=do.call(rbind,lbf$lbf) #lbf$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows 

gl.priorHits$nmin=n[gl$annot==1]
lbf.priorHits=compute.allBFs.fromZscores(Z.priorHits,VYY,gl.priorHits$nmin,gl.priorHits$maf,sigmaa)
lbf.priorHits.bigmat=do.call(rbind,lbf.priorHits$lbf)

lbf_sig1=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[1])
lbf=list(gamma=lbf_sig1$gamma,prior=lbf_sig1$prior)

#collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
collapse = function(x,nsigmaa){return(apply(matrix(x,ncol=nsigmaa),1,sum))}

#compute empirical bayes priors from GWAS hit SNPs
#gl.glhits = gl[gl$annot==1,] # subset of GWAS SNPs
gl.glhits = gl.priorHits # subset of GWAS SNPs
#CHECK_0: Note -- some top hits dropped because they did not have 'maf'; expecting ~881 below, not 843
#~~~
#> dim(gl.glhits)
#[1] 843  44
#~~~
lbf.glhits = lbf.bigmat[,gl$annot==1]
#lbf.glhits = lbf.priorHits.bigmat[,gl$annot==1]
lbf.glhits = lbf.priorHits.bigmat
ebprior.glhits = em.priorprobs(lbf.glhits,lbf$prior,100)
#ebprior.glhits2 = em.priorprobs(lbf.glhits,lbf$prior*runif(length(lbf$prior)),100)
#ebprior.glhits3 = em.priorprobs(lbf.glhits,lbf$prior*runif(length(lbf$prior)),100)
#ebprior.glhits4 = em.priorprobs(lbf.glhits,lbf$prior*runif(length(lbf$prior)),100)

ebprior.glhits.collapse =collapse(ebprior.glhits,length(sigmaa))
#ebprior.glhits.collapse2 =collapse(ebprior.glhits2,length(sigmaa))
#ebprior.glhits.collapse3 =collapse(ebprior.glhits3,length(sigmaa))
#ebprior.glhits.collapse4 =collapse(ebprior.glhits4,length(sigmaa))

pp.glhits = posteriorprob(lbf.glhits,ebprior.glhits) #posterior prob on models for gl hits
pp.glhits.collapse =  apply(pp.glhits,2,collapse, nsigmaa=length(sigmaa))

# this returns a list with elements  pU pD and pI
marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

#looking at which models are favored by the prior
cumsum(sort(ebprior.glhits.collapse,decreasing=TRUE))
lbf$gamma[order(ebprior.glhits.collapse,decreasing=TRUE),]
modelmatrix = cbind(lbf$gamma,ebprior.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
modelmatrix = data.frame(cbind(modelmatrix,cumsum(modelmatrix[,9])))
colnames(modelmatrix)= c("height","BMI","WHRadjBMI","WHR", "HIPadjBMI", "HIP", "WCadjBMI", "WC", "p","cump")

allassoc=(apply((modelmatrix[,1:8]>0),1,sum)==8) #vector of which represent situations in which all 8 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:8]>0),1,sum)==7) #vector of which represent situations in which 7 phenotypes are associated

sum(modelmatrix[allassoc,9]) #0.1902515
sum(modelmatrix[allbut1assoc,9]) #0.1538153

#~~~
#> sum(modelmatrix[allassoc,9])
#[1] 0.1902515
#> sum(modelmatrix[allbut1assoc,9])
#[1] 0.1538153
#~~~

#look at weight on each of WC, WCadjBMI, HIP, HIPadjBMi, WHR, WHRadjBMI, BMI, height being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,8]==0,9]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,7]==0,9]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,6]==0,9]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,5]==0,9]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,9]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,9]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,9]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,9]) 
#
# 0.04653936, 0.004498162, 0.01375071, 2.516849e-09, 0.02369146, 0.04216148, 1.45596e-09, 0.02317413

#~~~
> sum(modelmatrix[allbut1assoc & modelmatrix[,8]==0,9])
[1] 0.04653936
> sum(modelmatrix[allbut1assoc & modelmatrix[,7]==0,9])
[1] 0.004498162
> sum(modelmatrix[allbut1assoc & modelmatrix[,6]==0,9])
[1] 0.01375071
> sum(modelmatrix[allbut1assoc & modelmatrix[,5]==0,9])
[1] 2.516849e-09
> sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,9])
[1] 0.02369146
> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,9])
[1] 0.04216148
> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,9])
[1] 1.45596e-09
> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,9])
[1] 0.02317413
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or not WC, not WCadjBMI, not HIP, not HIPadjBMI, not WHR, notWHRadjBMI, or notBMI or notheight)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,8]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,7]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,6]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,5]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"GLC.MTedits.ForGIANT2014_5.AllPheno.vs2.gl.bestmodel.vs2.txt",quote=F,sep=" & ",row.names=F)

#~~~
> dim(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)))
[1] 261   3
#~~~

#~~~
#> head(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)))
#[,1]         [,2] [,3]
#[1,] "rs10929925" "7"  "0.52"
#[2,] "rs17806888" "2"  "0.02"
#[3,] "rs9914578"  "9"  "0.38"
#[4,] "rs780094"   "6"  "0.41"
#[5,] "rs1528435"  "7"  "0.63"
#[6,] "rs998584"   "2"  "1"
#.
#.
#.
#~~~

#Compare multivariate results with univariate results for the GL SNPs.
lbf.av.glhits= lbf.av(lbf.glhits,ebprior.glhits)
nsigma=length(sigmaa)
origprior = rep(c(0,lbf$prior[-1]),nsigma)
origprior = normalize(origprior)
lbf.av.origprior.glhits = lbf.av(lbf.glhits,origprior)
lbf.uni.glhits = lbf.uni(lbf.glhits,lbf$gamma)
lbf.all.glhits = lbf.all(lbf.glhits,lbf$gamma)

###Comparisons with lbf.av.glhits are problematic because the lbf.av.glhits
###prior are fitted to data (particularly for estimating sigmaa)
###whereas the lbf.uni are not (but could be if wanted?)
###So stick to comparisons of lbfs computed without learning prior from data
##pdf("nana1")
##layout(rbind(c(1,2)))
##plot(lbf.uni.glhits,lbf.av.origprior.glhits,xlim=c(4,15),ylim=c(4,15))
##abline(a=0,b=1)
##plot(lbf.uni.glhits,lbf.all.glhits,xlim=c(4,15),ylim=c(4,15))
##abline(a=0,b=1)
##dev.off()

#this establishes a threshold for subsequent analysis
min(lbf.av.glhits)
#~~~
#> min(lbf.av.glhits)
#[1] 2.895304
#~~~

#lbf_sig1=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[1])
#lbf_sig2=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[2])
#lbf.bigmat1 <- do.call(rbind, c(lbf_sig1$lbf, lbf_sig2$lbf))
#lbf=list(gamma=lbf_sig1$gamma,prior=lbf_sig1$prior)
#rm(lbf_sig1, lbf_sig2)
lbf_sig1=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[1])
lbf_sig2=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[2])
lbf.bigmat1 <- do.call(rbind, c(lbf_sig1$lbf, lbf_sig2$lbf))
lbf=list(gamma=lbf_sig1$gamma,prior=lbf_sig1$prior)
rm(lbf_sig1, lbf_sig2)
lbf_sig3=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[3])
lbf_sig4=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[4])
lbf.bigmat2 <- do.call(rbind, c(lbf_sig3$lbf, lbf_sig4$lbf))
rm(lbf_sig3, lbf_sig4)
lbf_sig5=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[5])
lbf_sig6=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[6])
lbf.bigmat3 <- do.call(rbind, c(lbf_sig5$lbf, lbf_sig6$lbf))
rm(lbf_sig5, lbf_sig6)
lbf_sig7=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[7])
lbf_sig8=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[8])
lbf.bigmat4 <- do.call(rbind, c(lbf_sig7$lbf, lbf_sig8$lbf))
rm(lbf_sig7, lbf_sig8)
lbf_sig9=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[9])
lbf_sig10=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[10])
lbf.bigmat5 <- do.call(rbind, c(lbf_sig9$lbf, lbf_sig10$lbf))
rm(lbf_sig9, lbf_sig10)
lbf_sig11=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[11])
lbf_sig12=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[12])
lbf.bigmat6 <- do.call(rbind, c(lbf_sig11$lbf, lbf_sig12$lbf))
rm(lbf_sig11, lbf_sig12)
lbf_sig13=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[13])
lbf_sig14=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[14])
lbf.bigmat7 <- do.call(rbind, c(lbf_sig13$lbf, lbf_sig14$lbf))
rm(lbf_sig13, lbf_sig14)

lbf.av.all.sub1 = lbf.av(rbind(lbf.bigmat1[,1:10000],lbf.bigmat2[,1:10000],lbf.bigmat3[,1:10000],lbf.bigmat4[,1:10000],lbf.bigmat5[,1:10000],lbf.bigmat6[,1:10000],lbf.bigmat7[,1:10000]),ebprior.glhits)
lbf.av.all.sub2 = lbf.av(rbind(lbf.bigmat1[,10001:20000],lbf.bigmat2[,10001:20000],lbf.bigmat3[,10001:20000],lbf.bigmat4[,10001:20000],lbf.bigmat5[,10001:20000],lbf.bigmat6[,10001:20000],lbf.bigmat7[,10001:20000]),ebprior.glhits)
lbf.av.all.sub3 = lbf.av(rbind(lbf.bigmat1[,20001:30000],lbf.bigmat2[,20001:30000],lbf.bigmat3[,20001:30000],lbf.bigmat4[,20001:30000],lbf.bigmat5[,20001:30000],lbf.bigmat6[,20001:30000],lbf.bigmat7[,20001:30000]),ebprior.glhits)
lbf.av.all.sub4 = lbf.av(rbind(lbf.bigmat1[,30001:40000],lbf.bigmat2[,30001:40000],lbf.bigmat3[,30001:40000],lbf.bigmat4[,30001:40000],lbf.bigmat5[,30001:40000],lbf.bigmat6[,30001:40000],lbf.bigmat7[,30001:40000]),ebprior.glhits)
lbf.av.all.sub5 = lbf.av(rbind(lbf.bigmat1[,40001:50000],lbf.bigmat2[,40001:50000],lbf.bigmat3[,40001:50000],lbf.bigmat4[,40001:50000],lbf.bigmat5[,40001:50000],lbf.bigmat6[,40001:50000],lbf.bigmat7[,40001:50000]),ebprior.glhits)
lbf.av.all.sub6 = lbf.av(rbind(lbf.bigmat1[,50001:56301],lbf.bigmat2[,50001:56301],lbf.bigmat3[,50001:56301],lbf.bigmat4[,50001:56301],lbf.bigmat5[,50001:56301],lbf.bigmat6[,50001:56301],lbf.bigmat7[,50001:56301]),ebprior.glhits)

length(c(lbf.av.all.sub1, lbf.av.all.sub2, lbf.av.all.sub3, lbf.av.all.sub4, lbf.av.all.sub5, lbf.av.all.sub6))
dim(lbf.bigmat1)
dim(lbf.bigmat2)
dim(lbf.bigmat3)
dim(lbf.bigmat4)
dim(lbf.bigmat5)
dim(lbf.bigmat6)
dim(lbf.bigmat7)

#~~~
#> dim(lbf.bigmat1)
#[1] 13122 56301
#> dim(lbf.bigmat2)
#[1] 13122 56301
#> dim(lbf.bigmat3)
#[1] 13122 56301
#> dim(lbf.bigmat4)
#[1] 13122 56301
#> dim(lbf.bigmat5)
#[1] 13122 56301
#> dim(lbf.bigmat6)
#[1] 13122 56301
#> dim(lbf.bigmat7)
#[1] 13122 56301
#> length(c(lbf.av.all.sub1, lbf.av.all.sub2, lbf.av.all.sub3, lbf.av.all.sub4, lbf.av.all.sub5, lbf.av.all.sub6))
#[1] 56301
#~~~

#lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = c(lbf.av.all.sub1, lbf.av.all.sub2, lbf.av.all.sub3, lbf.av.all.sub4, lbf.av.all.sub5, lbf.av.all.sub6)
#o = order(gl$chr, gl$pos)
#gl = gl[o,]

sub = gl[gl$nmin>50000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
newhits = sub[sub$annot==0 & sub$lbfav>2.895304 & sub$nmin>50000,c(1:5,7,33:40,45)]

#~~~
#> sub = gl[gl$nmin>50000,]
#> dim(sub)
#[1] 56301    45
#> sub=sub[l==1,]
#> dim(sub)
#[1] 1041   45
#~~~

#extract lbfs for all the new hits
#lbf.newhits= lbf.bigmat[,gl$nmin>50000]
#lbf.newhits= lbf.newhits[,l==1]
lbf.newhits <- rbind(lbf.bigmat1[,l==1],lbf.bigmat2[,l==1],lbf.bigmat3[,l==1],lbf.bigmat4[,l==1],lbf.bigmat5[,l==1],lbf.bigmat6[,l==1],lbf.bigmat7[,l==1]) 
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>2.895304 & sub$nmin>50000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#~~~
#> dim(lbf.newhits)
#[1] 91854  1041
#> lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>2.895304 & sub$nmin>50000]
#> dim(lbf.newhits)
#[1] 91854   411
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notWC, or notWCadjBMI, or notHIP, or notHIPadjBMI, or notWHR, or notWHRadjBMI, or notBMI or notheight)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,8]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,7]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,6]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,5]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#411 in results
#~~~
#> as.character(newhits$snp)
#.
#.
#.
#> length(bestclass)
#[1] 411
#> bestclass
#[1] 1 1 7 7 6 7 1 7 1 4 1 1 1 1 1 7 1 2 1 1 1 1 1 1 1 1 9 1 1 9 1 1 1 1 1 1 9
#[38] 1 1 9 1 1 1 7 1 1 1 7 1 1 1 1 1 7 1 1 1 1 1 2 1 1 1 1 1 1 9 1 1 1 1 6 1 1
#[75] 1 1 1 3 1 1 1 1 1 7 1 1 1 1 1 1 1 1 1 6 1 1 9 1 1 1 7 2 1 1 1 1 1 1 1 1 1
#[112] 1 7 1 1 1 1 1 1 7 1 7 2 1 1 1 1 1 9 1 1 1 7 1 1 1 1 6 2 1 7 1 1 1 1 1 9 1
#[149] 1 1 1 1 6 1 1 1 1 1 1 1 1 1 1 7 1 1 1 1 1 7 1 1 1 1 1 1 1 1 1 1 9 2 1 9 1
#[186] 1 1 1 7 7 1 1 1 1 1 1 7 4 1 1 1 1 1 1 1 7 9 2 1 6 6 1 7 1 1 1 1 1 1 1 1 9
#[223] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 7 1 1 6 1 1 1 7 9 1 7 1 1 1 1 1
#[260] 2 1 1 1 7 9 1 1 6 1 1 1 1 1 1 1 1 1 7 2 1 1 1 7 1 1 1 1 1 1 1 7 1 1 2 1 1
#[297] 1 1 1 1 9 1 1 7 1 9 1 1 1 1 1 7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 7 4 1
#[334] 1 7 1 1 7 1 7 1 1 1 1 1 2 7 2 1 1 1 1 1 2 1 1 7 1 1 1 1 7 1 1 1 1 1 2 1 1
#[371] 6 1 1 1 4 1 1 1 1 2 1 9 1 6 2 1 1 1 1 1 9 1 1 4 1 1 1 1 1 1 1 6 9 1 1 9 7
#[408] 1 1 1 7
#> apply(ppmatrix.newhits, 2,which.max)
#[1]  1  1  9 21 31  9  5  9  5 18  1  5 13  1  4  9  1  8  1 43  2  1  5  3  1
#[26] 20 10  4 13 10  5  1 11  5 45  1 10  1 13 10  1  1  5  9 43 15  3  9 20  3
#[51]  5  5  1 22  4  1  1 15  1 33  5 20  5  4  1  5 10  5  1  5 14 46  4  1  1
#[76]  4  3 27  1 11 14  4  1  9  5  5  5 14  1  1  1  4  5 26 17 43 10  1  1  2
#[101]  9  1  5  5  1  5 14 21  1  1  5  5  9 14 14  1  5  5  1  9  5  9  1  2 13
#[126]  1  1  3 10  1 11 43  9 43  4  5  5  2  8 13  9  1  2  1  1  1 10 11  5 15
#[151]  5  3  2  5  1  5  5 11  1  1  1  5  1 21  1  3  5  5  5 21  5 13  5  1  1
#[176]  5  5 31  2  5 10 33 46 10  2  4  5  1  9  9  1 43  1  2 17  6  9 18  5  1
#[201]  5  5  1  1  5  9 10 25  1  1 31  4 22 65  1  4 21  2 20  5 14 10  5  5  2
#[226]  5  5  5  5  5  1  1 14 15  5 13  1  1  4  5  5  2  1  9 20  1 31  4 15 11
#[251]  9 10  1  9  4 43  5 14  2 11  4 13  1  1 10  5  1 46 11  5  1 13  2  1  1
#[276]  1  1  9  8  1 14  1  9  5  1  3  1  1  1 43  9  5  1  8 13 13  5 70  1 13
#[301] 46  5  5 21  1 46  3  5  4  1  3  9 17 13  1  1  9 17  5  5 14  1  2  6  1
#[326]  2 13  1  1  5 22 18 43 14  9  1  5  9 14 16  1  1 20  5  3 33  9 10 13  5
#[351]  5  5  1 33  1  1 22  1  1  5  3  9  1  4  5  1  1  1  5  5  9 13  4 15 18
#[376] 11  1  3  1  1  1 10  3 26 33  1  1 14  4  1 10  5  9 52  1  5 13 11  2  1
#[401]  4 31 10  1 43 10 22 44 14  1  9
#> apply(ppmatrix.newhits, 2,max)[apply(ppmatrix.newhits, 2,max)>.75]
#[1] 0.9189713 0.8542716 0.7506763 0.7524897 0.8471900 0.8205852 0.8691441
#[8] 0.7783979 0.8433007 0.9315604 0.8283082 0.8704295 0.8048510 0.8831607
#[15] 0.8000743 0.8068011 0.7844906 0.7948894 0.8510579 0.8153530 0.8283186
#[22] 0.8002370 0.8247344 0.7738733 0.7940703 0.9625288 0.8434042 0.9133117
#[29] 0.9304080 0.8767345 0.7721371 0.8960261 0.7521920 0.8488992 0.7950071
#[36] 0.8437116 0.8821659 0.7832316 0.8369998 0.8757812
#> apply(ppmatrix.newhits, 2,which.max)[apply(ppmatrix.newhits, 2,max)>.5]
#[1]  1  9 21 31  9  5  5  1  5 13  1  5 10  5  1  1 10 10  1 15  9  5  1  1 15
#[26]  1  5  5  1  5  5  1  1  1  5  5 14  1  1  5 17 43 10  1  1  5  5  1  5  9
#[51]  1  5  5  1 43 43  5  5  1 10 11  5 15  5  5  1  5  1  5  1  1  5  5 21  5
#[76]  5  5  5  2  5 10  5  1  9  1  1 17 18  5  5  5  1  5  9 10  2  5  5  5  5
#[101]  5  1  5  5  1 11 10  1 43  5  1  5  1 11  1  1  5  1  1  1 43  5 13  5  1
#[126] 13 21  5  1  5  5 14  1 13  5 18  1  5  5  9  5  5  1  1  5  5  1  5 18 11
#[151]  1 10 33  1  5 52  5 13 10 10  1
#> apply(ppmatrix.newhits, 2,max)[apply(ppmatrix.newhits, 2,max)>.5]
#[1] 0.5004890 0.5784610 0.5205169 0.5221882 0.6781449 0.9189713 0.5872205
#[8] 0.6865888 0.8542716 0.6027869 0.6075081 0.7506763 0.5536017 0.7524897
#[15] 0.5158207 0.7060345 0.8471900 0.5194476 0.6016343 0.8205852 0.5045100
#[22] 0.6521827 0.6275769 0.7427017 0.6673476 0.6612960 0.7069084 0.8691441
#[29] 0.6661230 0.6213632 0.6510104 0.6547819 0.7783979 0.7213823 0.6866550
#[36] 0.8433007 0.5892078 0.7171824 0.6903062 0.9315604 0.6920957 0.5271843
#[43] 0.5235760 0.5145779 0.8283082 0.8704295 0.6997309 0.6314845 0.8048510
#[50] 0.5188278 0.7266389 0.8831607 0.6471322 0.7116318 0.6797641 0.8000743
#[57] 0.8068011 0.5037725 0.7072781 0.7445771 0.7324901 0.6410747 0.6664809
#[64] 0.6967524 0.7844906 0.5366190 0.6283673 0.6756784 0.5870642 0.7948894
#[71] 0.7152171 0.5161147 0.7474826 0.5420676 0.8510579 0.8153530 0.5333000
#[78] 0.8283186 0.8002370 0.8247344 0.6732977 0.6138360 0.5911746 0.5358608
#[85] 0.5119000 0.5928136 0.7738733 0.5164170 0.6687901 0.5078716 0.6353523
#[92] 0.7940703 0.6728617 0.5437806 0.5019119 0.9625288 0.8434042 0.9133117
#[99] 0.9304080 0.8767345 0.7721371 0.7130767 0.8960261 0.5937418 0.6362862
#[106] 0.5937072 0.7019343 0.5189174 0.6618059 0.7190836 0.7099095 0.5636841
#[113] 0.5843118 0.6317536 0.5220670 0.7396838 0.6699234 0.5430604 0.7521920
#[120] 0.5765165 0.8488992 0.6375972 0.5388773 0.6534821 0.5443532 0.5420972
#[127] 0.5081312 0.7950071 0.6702207 0.6236940 0.5217426 0.5482685 0.5949710
#[134] 0.5552022 0.6243978 0.8437116 0.6269425 0.8821659 0.6167578 0.6031187
#[141] 0.5465188 0.7272837 0.7832316 0.5035239 0.7416224 0.5377178 0.7481916
#[148] 0.8369998 0.8757812 0.7097081 0.7129656 0.7310120 0.5623154 0.6282847
#[155] 0.7101294 0.6675403 0.5769318 0.6301486 0.7150226 0.7224416 0.5216311
#~~~

##pdf("plots.bychr.vs2.pdf")
##for(i in 1:22){
## if(sum(sub$chr==i)>0){
##plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>50000))
##abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
##abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>4.608818],col=3)
}
}
#dev.off()

###try to plot all hits to see relationships?
##tophits =  sub[sub$lbfav>3.323736 & sub$nmin>50000,]
##betahat = tophits[,c(8,11,14,17)]
##betahat.scaled = betahat/apply(betahat,1,sd)
##betahat.scaled.pr = prcomp(betahat.scaled,center=F)
##pdf("plots.betahats.vs2.pdf")
##plot(abs(betahat.scaled.pr$x[,1]),abs(betahat.scaled.pr$x[,2]),type="n")
##text(abs(betahat.scaled.pr$x[,1]),abs(betahat.scaled.pr$x[,2]),as.character(tophits$chr))
##dev.off()

##hla = gl[gl$chr==6 & gl$pos<40e6 & gl$pos>28e6,]
##hla.z = hla[,22:25]

### > plot(hla.pr$x[,1])
### > plot(hla.pr$x[,2])
### > plot(hla.pr$x[,2]^2)
### > plot(hla.pr$x[,2]^2,hla.pr$x[,1]^2)
### > image(hla.z)
### Error in image.default(hla.z) : 'z' must be a matrix
### > image(as.matrix(hla.z))
### > image(as.matrix(hla.z^2))
### > image(as.matrix(abs(hla.z)))
### 

##zhla = head(gl[gl$chr==6 & gl$pos<40e6 & gl$pos>28e6,],600)[,22:25]
##pdf("plots.cor.zhla.vs2.pdf")
##image(cor(t(zhla))^2)
##dev.off()

##pdf("plots.chr1lbfav.vs2.pdf")
##plot(gl$pos[gl$chr==1],gl$lbfav[gl$chr==1],ylim=c(0,10),xlim=c(27102620-10^5,27102620+10^5))
##dev.off()

write.table(file="GLC.MTedits.ForGIANT2014_5.AllPheno.vs2.newtophits.vs1.txt",cbind(sub[sub$lbfav>2.895304 & sub$nmin>50000 & sub$annot==0,c(1:5,7)],round(sub[sub$lbfav>2.895304 & sub$nmin>50000 & sub$annot==0,c(33:40,45)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)










