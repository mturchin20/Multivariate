set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/HaemgenRBC2012/Vs1/HaemgenRBC2012.RSS0.vs1.SignCrrct.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/HaemgenRBC2012/Vs1/HaemgenRBC2012.dtlesssignif.vs1.annot.vs1.SignCrrct.vs1.MAF.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$pos) & gl$pos != "None",]
Z = cbind(gl$Z.RBC,gl$Z.PCV,gl$Z.MCV,gl$Z.MCH,gl$Z.MCHC,gl$Z.Hb)

#~~~
#> dim(gl)
#[1] 4875   26
#> gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$pos) & gl$pos != "None",]
#> dim(gl)
#[1] 4875   26
#> 
#> head(gl)
#snp chr       pos   maf     p_RBC n_RBC     p_PCV n_PCV     p_MCV
#1 rs10015639   4  88110589 0.167 6.918e-01 53394 1.275e-02 52821 3.893e-06
#2 rs10023020   4 122944054 0.633 1.069e-03 53370 7.410e-01 52798 2.474e-08
#3 rs10023050   4  88283455 0.650 2.702e-04 49161 6.699e-07 52816 2.782e-02
#4 rs10023056   4  88283490 0.758 4.335e-01 53397 8.750e-05 52824 1.611e-06
#5  rs1002488   1  26738509 0.754 1.799e-06 53380 1.877e-03 52806 5.919e-05
#6 rs10029915   4  88278641 0.788 6.214e-01 53297 4.243e-04 52724 1.791e-06
#n_MCV     p_MCH n_MCH p_MCHC n_MCHC      p_Hb  n_Hb annot     Z.RBC    Z.PCV
#1 57848 8.436e-07 51445 0.2132  56208 4.890e-03 60888     0 0.3964135 4.617015
#2 57823 1.664e-05 51421 0.4811  56185 6.573e-01 60864     2 3.2717078 5.575091
#3 51961 2.709e-02 47287 0.9659  50394 1.249e-07 54966     0 3.6423123 2.199815
#4 57850 7.186e-07 51447 0.2793  56211 8.145e-06 60891     0 0.7832165 4.796950
#5 57832 4.393e-04 51430 0.4254  56194 8.674e-04 60874     0 4.7747843 4.016017
#6 57750 1.115e-06 51347 0.3544  56111 4.508e-05 60791     0 0.4938672 4.775681
#Z.MCV    Z.MCH   Z.MCHC       Z.Hb   mvstat       mvp     unip
#1 2.4906769 4.925001 2.814192 1.24481435 44.15364  7.161705 6.073863
#2 0.3305294 4.305778 0.443644 0.70453446 42.72551  6.878838 7.606600
#3 4.9698926 2.210219 5.286176 0.04275103 43.60996  7.053918 6.903438
#4 3.9228710 4.956268 4.461337 1.08189317 60.06907 10.360721 6.143513
#5 3.1090342 3.515274 3.330340 0.79708787 44.73717  7.277523 5.744969
#6 3.5244893 4.870174 4.079771 0.92608848 55.49899  9.434698 5.952725
> dim(gl)
[1] 3995   26
> gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$pos) & gl$pos != "None",]
> dim(gl)
[1] 3995   26
> head(gl)
         snp chr       pos   maf     p_RBC n_RBC      p_PCV n_PCV      p_MCV
1 rs10015639   4  88110589 0.167 6.918e-01 53394 -1.275e-02 52821 -3.893e-06
2 rs10023020   4 122944054 0.367 1.069e-03 53370  7.410e-01 52798 -2.474e-08
3 rs10023050   4  88283455 0.350 2.702e-04 49161  6.699e-07 52816  2.782e-02
4 rs10023056   4  88283490 0.242 4.335e-01 53397  8.750e-05 52824  1.611e-06
5 rs10032244   4  88136955 0.258 5.507e-01 53541 -9.583e-02 52969 -2.020e-04
6 rs10073340   5   1374873 0.133 5.983e-06 48888  1.743e-01 48314 -1.283e-05
  n_MCV      p_MCH n_MCH    p_MCHC n_MCHC       p_Hb  n_Hb annot      Z.RBC
1 57848 -8.436e-07 51445 -0.213200  56208 -4.890e-03 60888     0  0.3964135
2 57823 -1.664e-05 51421  0.481100  56185  6.573e-01 60864     2 -3.2717078
3 51961  2.709e-02 47287  0.965900  50394  1.249e-07 54966     0  3.6423123
4 57850  7.186e-07 51447  0.279300  56211  8.145e-06 60891     0 -0.7832165
5 57995 -2.986e-07 51593 -0.008596  56356 -1.443e-02 61035     0  0.5967115
6 53353 -9.924e-06 46955  0.909100  51311  1.078e-01 56386     0  4.5269892
      Z.PCV      Z.MCV     Z.MCH    Z.MCHC        Z.Hb   mvstat      mvp
1 -4.617015 -2.4906769 -4.925001 -2.814192 -1.24481435 28.84407 4.186335
2  5.575091 -0.3305294  4.305778 -0.443644 -0.70453446 40.87195 6.513031
3  2.199815  4.9698926  2.210219  5.286176  0.04275103 33.87158 5.147362
4 -4.796950 -3.9228710 -4.956268 -4.461337 -1.08189317 35.83439 5.527444
5 -3.716502 -1.6654149 -5.124328 -2.446377 -2.62771698 30.02832 4.410898
6 -4.362985  1.3585160 -4.418823  1.608161  0.11417383 40.51522 6.442809
      unip
1 6.073863
2 7.606600
3 6.903438
4 6.143513
5 6.524910
6 5.223081
#~~~

n = cbind(gl$n_RBC, gl$n_PCV, gl$n_MCV, gl$n_MCH, gl$n_MCHC, gl$n_Hb)
n = apply(n,1,min)
gl$nmin=n
sigmaa=c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)
lbf=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa)
lbf.bigmat=do.call(rbind,lbf$lbf) #lbf$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows 

#collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
collapse = function(x,nsigmaa){return(apply(matrix(x,ncol=nsigmaa),1,sum))}

#compute empirical bayes priors from GWAS hit SNPs
gl.glhits = gl[gl$annot==1,] # subset of GWAS SNPs
#20151006 CHECK_0: This CHECK_0 is from GIANT2014_5 not sure if it applies here though obviously 4 top hits are missing: CHECK_0: Note -- some top hits dropped because they did not have 'maf'; expecting ~881 below, not 843
#~~~
#> dim(gl.glhits)
#[1] 71 27
> dim(gl.glhits)
[1] 66 27
#~~~
lbf.glhits = lbf.bigmat[,gl$annot==1]
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
modelmatrix = data.frame(cbind(modelmatrix,cumsum(modelmatrix[,7])))
colnames(modelmatrix)= c("RBC","PCV","MCV","MCH", "MCHC", "Hb", "p","cump")

allassoc=(apply((modelmatrix[,1:6]>0),1,sum)==6) #vector of which represent situations in which all 6 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:6]>0),1,sum)==5) #vector of which represent situations in which 5 phenotypes are associated

sum(modelmatrix[allassoc,7]) # 0.5976083
sum(modelmatrix[allbut1assoc,7]) # 0.1716537

#~~~
#> sum(modelmatrix[allassoc,7])
#[1] 0.5976083
#> sum(modelmatrix[allbut1assoc,7])
#[1] 0.1716537
> sum(modelmatrix[allassoc,7]) # 0.5976083
[1] 0.2287069
> sum(modelmatrix[allbut1assoc,7]) # 0.1716537
[1] 0.4213086
#~~~

#look at weight on each of Hb, MCHC, MCH, MCV, PCV, RBC being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,6]==0,7]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,5]==0,7]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,7]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,7]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,7]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,7]) 
# 7.918026e-05, 4.605305e-17, 4.58361e-28, 0.07212402, 0.05711201, 0.04233848 
# 

#~~~
#> sum(modelmatrix[allbut1assoc & modelmatrix[,6]==0,7])
#[1] 7.918026e-05
#> sum(modelmatrix[allbut1assoc & modelmatrix[,5]==0,7])
#[1] 4.605305e-17
#> sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,7])
#[1] 4.58361e-28
#> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,7])
#[1] 0.07212402
#> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,7])
#[1] 0.05711201
#> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,7])
#[1] 0.04233848
> sum(modelmatrix[allbut1assoc & modelmatrix[,6]==0,7])
[1] 0.1549341
> sum(modelmatrix[allbut1assoc & modelmatrix[,5]==0,7])
[1] 0.1233723
> sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,7])
[1] 0.03926759
> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,7])
[1] 2.615689e-09
> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,7])
[1] 0.05094429
> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,7])
[1] 0.05279029
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or not Hb, not MCHC, not MCH, notMCV, or notPCV or notRBC)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,6]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,5]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"GLC.MTedits.ForHaemgenRBC2012.AllPheno.vs1.gl.bestmodel.vs2.SignCrrct.vs1.txt",quote=F,sep=" & ",row.names=F)

#~~~
#> dim(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)))
#[1] 14  3
#> cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2))
#[,1]         [,2] [,3]  
#[1,] "rs10445033" "5"  "0.27"
#[2,] "rs11042125" "6"  "0.88"
#[3,] "rs12150672" "6"  "0.06"
#[4,] "rs13061823" "5"  "0.45"
#[5,] "rs13152701" "5"  "0.43"
#[6,] "rs13219787" "7"  "0.98"
#[7,] "rs1532085"  "6"  "0.11"
#[8,] "rs2097775"  "7"  "0.96"
#[9,] "rs2159213"  "6"  "0.01"
#[10,] "rs2867932"  "6"  "0.97"
#[11,] "rs3184504"  "6"  "0.11"
#[12,] "rs4969184"  "6"  "0.88"
#[13,] "rs7936461"  "6"  "0.09"
#[14,] "rs855791"   "7"  "1"   

dim(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)))
[1] 51  3

> table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2))[,2])

 2  3  4  6  7 
24 13  3  7  4 

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
#[1] 6.695078
> min(lbf.av.glhits)
[1] 5.083439
#~~~

lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]

sub = gl[gl$nmin>20000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
newhits = sub[sub$annot==0 & sub$lbfav>5.083439 & sub$nmin>20000,c(1:4,17:28)]

#quantile(gl[gl$annot==1,][apply(gl[gl$annot==1,18:23], 1, which.max)==3,]$Z.MCV)

#~~~
#> sub = gl[gl$nmin>20000,]
#> dim(sub)
#[1] 4872    28
#> sub=sub[l==1,]
#> dim(sub)
#[1] 142   28
> l=indephits(sub$lbfav,sub$chr,sub$pos)
> sub=sub[l==1,]
> dim(sub)
[1] 149  28
#~~~

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>20000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>5.083439 & sub$nmin>20000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#~~~
#> dim(lbf.newhits)
#[1] 10206   142
#> lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>5.083439 & sub$nmin>20000]
#> dim(lbf.newhits)
#[1] 10206    26
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notHb, or notMCHC, or notMCH, or notMCV, or notPCV or notRBC)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,6]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,5]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#411 in results
#~~~

> dim(newhits)
[1] 70 16

#> cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))
#bestclass                     
#[1,] "rs11672923" "6"       "0.576713718675208" 
#[2,] "rs12422977" "6"       "0.177162709036857" 
#[3,] "rs12513744" "1"       "0.779014061654097" 
#[4,] "rs1256061"  "6"       "0.325786969939427" 
#[5,] "rs12878795" "1"       "0.869373519623004" 
#[6,] "rs13064931" "6"       "0.84707175242557"  
#[7,] "rs13194491" "7"       "0.819668170046279" 
#[8,] "rs13254494" "1"       "0.837082879491064" 
#[9,] "rs13708"    "1"       "0.987923470673524" 
#[10,] "rs17630235" "1"       "0.118715380660887" 
#[11,] "rs17701999" "1"       "0.907914504295556" 
#[12,] "rs1883349"  "6"       "0.177712778597281" 
#[13,] "rs1934661"  "1"       "0.0730559273380986"
#[14,] "rs2242652"  "1"       "0.960807995180489" 
#[15,] "rs2427350"  "6"       "0.0615688217144139"
#[16,] "rs2793678"  "1"       "0.212683462191235" 
#[17,] "rs2853925"  "1"       "0.989177760195318" 
#[18,] "rs3118228"  "1"       "0.813757957371494" 
#[19,] "rs3118362"  "7"       "0.924857508663405" 
#[20,] "rs3856444"  "1"       "0.643585444840029" 
#[21,] "rs6060399"  "1"       "0.970629751174068" 
#[22,] "rs6656196"  "1"       "0.985052140241061" 
#[23,] "rs6770091"  "1"       "0.0807155825316071"
#[24,] "rs7909074"  "1"       "0.982088215261971" 
#[25,] "rs891463"   "1"       "0.976520214799613" 
#[26,] "rs916888"   "1"       "0.604894156909518" 
#~~~


 newhits_check = sub[sub$annot==0 & sub$nmin>20000,c(1:4,17:28)]
 newhits_check <- cbind(newhits_check, apply(newhits_check[,6:11], 1, function(x) { val1 <- sum(abs(x)); return(val1) } ))
 cor(newhits_check[,16], newhits_check[,17])
 > cor(newhits_check[,16], newhits_check[,17])
 [1] 0.7064361
> cor(newhits_check$mvstat, newhits_check$lbfav)
[1] 0.7389623
> cor(newhits_check$mvp, newhits_check$lbfav)
[1] 0.7405769
> cor(newhits_check$unip, newhits_check$lbfav)
[1] 0.557489

> dim(newhits[apply(newhits[,6:11], 1, function(x) { val1 <- max(abs(x)); return(val1) }) > 6,])
[1]  3 16
> dim(gl.glhits[apply(gl.glhits[,18:23], 1, function(x) { val1 <- max(abs(x)); return(val1) } ) > 6,])
[1] 47 27




##pdf("plots.bychr.vs2.pdf")
##for(i in 1:22){
## if(sum(sub$chr==i)>0){
##plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>20000))
##abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
##abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>4.608818],col=3)
}
}
#dev.off()

###try to plot all hits to see relationships?
##tophits =  sub[sub$lbfav>3.323736 & sub$nmin>20000,]
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

write.table(file="GLC.MTedits.ForHaemgenRBC2012.AllPheno.vs1.newtophits.vs1.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>5.083439 & sub$nmin>20000 & sub$annot==0,c(1:5,17)],round(sub[sub$lbfav>5.083439 & sub$nmin>20000 & sub$annot==0,c(18:28)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)










