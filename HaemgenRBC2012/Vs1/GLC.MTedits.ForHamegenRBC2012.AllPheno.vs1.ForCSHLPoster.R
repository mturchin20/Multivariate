set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/HaemgenRBC2012/Vs1/HaemgenRBC2012.RSS0.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/HaemgenRBC2012/Vs1/HaemgenRBC2012.dtlesssignif.vs1.annot.txt.gz", header=T)
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
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or not Hb, not MCHC, not MCH, notMCV, or notPCV or notRBC)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,6]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,5]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"GLC.MTedits.ForHaemgenRBC2012.AllPheno.vs1.gl.bestmodel.vs2.txt",quote=F,sep=" & ",row.names=F)

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
#~~~

lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]

sub = gl[gl$nmin>20000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
newhits = sub[sub$annot==0 & sub$lbfav>6.695078 & sub$nmin>20000,c(1:4,17:28)]

#quantile(gl[gl$annot==1,][apply(gl[gl$annot==1,18:23], 1, which.max)==3,]$Z.MCV)

#~~~
#> sub = gl[gl$nmin>20000,]
#> dim(sub)
#[1] 4872    28
#> sub=sub[l==1,]
#> dim(sub)
#[1] 142   28
#~~~

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>20000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>6.695078 & sub$nmin>20000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#~~~
#> dim(lbf.newhits)
#[1] 10206   142
#> lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>6.695078 & sub$nmin>20000]
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

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/HaemgenRBC2012/Vs1/GLC.MTedits.ForHaemgenRBC2012.AllPheno.vs1.ForCSHLPoster.ModelsVslbfavNewOld.vs1.jpeg", res=300, width=5000, height=4000)
par(mfrow=c(2,2))

head(modelmatrix)
modelmatrix.paste <- apply(modelmatrix[,1:6], 1, function(x) paste(x, collapse="_"))
table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),7]

snpmodels <- cbind(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)], apply(ppmatrix.newhits, 2, max), apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] }), apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-2)[length(x)-2] }))
topmodels <- c()

for (i in unique(snpmodels[,1])) {
	snpmodels.subset <- matrix(snpmodels[snpmodels[,1]==i,], nrow=length(which(snpmodels[,1]==i)))
	topmodels <- rbind(topmodels, cbind(i, nrow(snpmodels.subset), mean(as.numeric(snpmodels.subset[,2])), sd(as.numeric(snpmodels.subset[,2])), mean(as.numeric(snpmodels.subset[,2]) - as.numeric(snpmodels.subset[,3])), sd(as.numeric(snpmodels.subset[,2]) - as.numeric(snpmodels.subset[,3])), modelmatrix[modelmatrix.paste == i,7]))
}
topmodels <- topmodels[order(as.numeric(topmodels[,2]), decreasing=TRUE),]

#topmodels <- cbind(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])), table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]), modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),7])
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
	oldtopmodels <- rbind(oldtopmodels, cbind(i, nrow(oldsnpmodels.subset), mean(as.numeric(oldsnpmodels.subset[,2])), sd(as.numeric(oldsnpmodels.subset[,2])), mean(as.numeric(oldsnpmodels.subset[,2]) - as.numeric(oldsnpmodels.subset[,3])), sd(as.numeric(oldsnpmodels.subset[,2]) - as.numeric(oldsnpmodels.subset[,3])), modelmatrix[modelmatrix.paste == i,7])) 
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
#RBC PCV MCV MCH MCHC Hb          p      cump
#1   2   1   2   1    2  1 0.17555634 0.1755563
#2   1   0   1   0    1  0 0.10805000 0.2836063
#3   1   1   2   1    2  2 0.09919851 0.3828049
#4   1   2   1   2    1  2 0.09264908 0.4754539
#5   2   1   2   1    2  2 0.08352962 0.5589835
#6   1   0   1   2    1  1 0.05679712 0.6157807
#> modelmatrix.paste <- apply(modelmatrix[,1:6], 1, function(x) paste(x, collapse="_"))
#> table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
#
#0_1_2_1_1_1 1_0_1_0_1_0 1_0_1_2_1_1 1_1_2_1_2_2 1_2_1_2_1_2 2_1_2_1_2_1
#2           8           2           6           3           5
#> modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),7]
#[1] 0.04233848 0.10805000 0.05679712 0.09919851 0.09264908 0.17555634
#>
#> topmodels
#i                                                         
#[1,] "1_0_1_0_1_0" "8" "0.746450507195395" "0.165286680790238" 
#[2,] "1_1_2_1_2_2" "6" "0.696986382037346" "0.211359842867695" 
#[3,] "2_1_2_1_2_1" "5" "0.567017758247753" "0.111652627100599" 
#[4,] "1_2_1_2_1_2" "3" "0.832357270421061" "0.205534912783212" 
#[5,] "1_0_1_2_1_1" "2" "0.70998562682716"  "0.193851709747913" 
#[6,] "0_1_2_1_1_1" "2" "0.872262839354842" "0.0743800946446978"
#
#[1,] "0.594933612087672" "0.250271007423985" "0.108049997636366" 
#[2,] "0.569887532905965" "0.242754334717206" "0.0991985105860755"
#[3,] "0.408117779055712" "0.129364264364991" "0.175556341828659" 
#[4,] "0.745028301570999" "0.327696686722976" "0.092649077183807" 
#[5,] "0.489027162017573" "0.384698811986865" "0.05679711530471"  
#[6,] "0.785602061687708" "0.12722128214335"  "0.0423384821635489"
#> oldtopmodels
#i                                                          
#[1,] "2_1_2_1_2_1" "16" "0.652963546591385" "0.159357382367687" 
#[2,] "1_0_1_0_1_0" "8"  "0.893394658392044" "0.0973146804835773"
#[3,] "1_2_1_2_1_2" "7"  "0.866201314021426" "0.179384042175176" 
#[4,] "2_1_0_1_2_2" "7"  "0.341653331948777" "0.0510115808906829"
#[5,] "1_1_2_1_2_2" "7"  "0.655060876551673" "0.154737113843472" 
#[6,] "1_1_0_0_1_1" "5"  "0.764370035680309" "0.180268215135302" 
#[7,] "2_1_2_1_2_2" "5"  "0.606673830055071" "0.246577735061375" 
#[8,] "1_0_1_2_1_1" "4"  "0.784690176871887" "0.233841830755469" 
#[9,] "0_1_2_1_1_1" "3"  "0.978832820180443" "0.016882845448066" 
#[10,] "1_1_0_1_0_2" "3"  "0.698140347782883" "0.251065009298668" 
#[11,] "1_1_1_1_1_1" "2"  "0.919257699101056" "0.11418677848921"  
#[12,] "1_1_1_1_2_1" "2"  "0.7616747467773"   "0.239487947329969" 
#[13,] "1_1_2_1_2_1" "1"  "0.93886551921131"  NA                  
#[14,] "1_1_2_2_1_1" "1"  "0.861616284582913" NA                  
#
#[1,] "0.486447046424441" "0.225233031216826"  "0.175556341828659" 
#[2,] "0.820638380792041" "0.15986277715676"   "0.108049997636366" 
#[3,] "0.788140507848196" "0.276444001237045"  "0.092649077183807" 
#[4,] "0.106028971821827" "0.0438295383356529" "0.0468006491072463"
#[5,] "0.476864880209375" "0.228328918253945"  "0.0991985105860755"
#[6,] "0.64000216022323"  "0.288245167127103"  "0.0541002932997882"
#[7,] "0.443932277448672" "0.332351369422735"  "0.0835296171404135"
#[8,] "0.662091906324268" "0.373306418163191"  "0.05679711530471"  
#[9,] "0.957809492540878" "0.0337535296254892" "0.0423384821635489"
#[10,] "0.560184154872811" "0.367199871782779"  "0.0415046168905986"
#[11,] "0.838773275745415" "0.228008865857161"  "0.040242466831677" 
#[12,] "0.586473779255963" "0.389704848566987"  "0.0278442374380658"
#[13,] "0.90397654343789"  NA                   "0.0558601906828733"
#[14,] "0.815258372785272" NA                   "0.0227171546876558"
#> summary(lm(apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.38526 -0.09511  0.02435  0.12064  0.27211
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.591452   0.083457   7.087 2.51e-07 ***
#newhits$lbfav 0.012987   0.007865   1.651    0.112
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.1732 on 24 degrees of freedom
#Multiple R-squared:  0.102,     Adjusted R-squared:  0.06459
#F-statistic: 2.726 on 1 and 24 DF,  p-value: 0.1117
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
#Min       1Q   Median       3Q      Max
#-0.50970 -0.09422  0.01996  0.16502  0.41145
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.41223    0.11515   3.580  0.00151 **
#newhits$lbfav  0.01701    0.01085   1.567  0.13016   
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.239 on 24 degrees of freedom
#Multiple R-squared:  0.09284,   Adjusted R-squared:  0.05504 
#F-statistic: 2.456 on 1 and 24 DF,  p-value: 0.1302
#
#> 
#> summary(lm(apply(ppmatrix, 2, max) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix, 2, max) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.43042 -0.16763  0.03325  0.19552  0.28849
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.6671288  0.0315896  21.119   <2e-16 ***
#gl[gl$annot == 1, ]$lbfav 0.0020129  0.0007763   2.593   0.0116 *
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2141 on 69 degrees of freedom
#Multiple R-squared:  0.08878,   Adjusted R-squared:  0.07558
#F-statistic: 6.723 on 1 and 69 DF,  p-value: 0.01161
#
#> summary(lm((apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = (apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) {
#sort(x, partial = length(x) - 1)[length(x) - 1]
#})) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.49921 -0.27842 -0.01093  0.28315  0.42383
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)                0.51338    0.04437  11.570   <2e-16 ***
#gl[gl$annot == 1, ]$lbfav  0.00273    0.00109   2.503   0.0147 *
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.3007 on 69 degrees of freedom
#Multiple R-squared:  0.08326,   Adjusted R-squared:  0.06997
#F-statistic: 6.267 on 1 and 69 DF,  p-value: 0.01467
#
#> 
#> dev.off()
#~~~




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

write.table(file="GLC.MTedits.ForHaemgenRBC2012.AllPheno.vs1.newtophits.vs1.txt",cbind(sub[sub$lbfav>6.695078 & sub$nmin>20000 & sub$annot==0,c(1:5,17)],round(sub[sub$lbfav>6.695078 & sub$nmin>20000 & sub$annot==0,c(18:28)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)





#20160428
~~~
#Prior values
> modelmatrix[1:5,]
  RBC PCV MCV MCH MCHC Hb          p      cump
  1   2   1   2   1    2  1 0.17555634 0.1755563
  2   1   0   1   0    1  0 0.10805000 0.2836063
  3   1   1   2   1    2  2 0.09919851 0.3828049
  4   1   2   1   2    1  2 0.09264908 0.4754539
  5   2   1   2   1    2  2 0.08352962 0.5589835


#> modelmatrix[1:5,]
#  TC TG HDL LDL         p      cump
#  1  1  1   1   1 0.2384938 0.2384938
#  2  1  2   1   1 0.1538137 0.3923075
#  3  2  1   1   2 0.1359187 0.5282262
#  4  1  1   2   1 0.1169160 0.6451422
#  5  2  1   1   1 0.1061912 0.7513333

gl.glhits$lbfav = lbf.av.glhits
gl.glhits$lbfall = lbf.all.glhits
gl.glhits$lbfuni = lbf.uni.glhits

c(1:4,17:28)
#newhits = sub[sub$annot==0 & sub$lbfav>6.695078 & sub$nmin>50000,c(1:4,17:28)]
failhits = sub[sub$annot==0 & sub$lbfav<=6.695078,c(1:4,17:28)]
prevhits = gl.glhits[,c(1:4,17:28)]

GLC.MTedits.ForHamegenRBC2012.AllPheno.vs1.ForCSHLPoster
write.table(file="failtophits.vs3.ForCSHLPoster.txt",cbind(sub[sub$lbfav<=6.695078 & sub$nmin>50000 & sub$annot==0,c(4,2:3,7)],round(sub[sub$lbfav<=6.695078 & sub$nmin>50000 & sub$annot==0,c(22:25,30)],digits=1)),quote=FALSE,sep= " ", row.names=FALSE)


write.table(file="glhits.vs3.ForCSHLPoster.txt",gl.glhits[,c(4,2,3,7,22,23,24,25,30,31,32)],quote=FALSE,row.names=FALSE)

jpeg("HaemgenRBC2012.vs3.ForCSHLPoster.NHSSlides.jpg", width=4500, height=2250, res=300)
#jpeg("HaemgenRBC2012.vs3.ForCSHLPoster.NHSSlides.SDvsMax.jpg", width=2000, height=2000, res=300)
#jpeg("HaemgenRBC2012.vs3.ForCSHLPoster.NHSSlides.ScaledSDvsMax.jpg", width=2250, height=2250, res=300)
#jpeg("HaemgenRBC2012.vs3.ForCSHLPoster.NHSSlides.ScaledSDvsLBF.jpg", width=2000, height=2000, res=300)

#par(xpd=FALSE, mar=c(6,6,6,6), oma=c(2,2,2,2))
par(mfrow=c(1,2), mar=c(4,5,5,4), oma=c(2.25,2.25,2.25,2.25))

#plot(c(apply(prevhits[,6:11], 1, sd), apply(newhits[,6:11], 1, sd)), c(apply(prevhits[,6:11], 1, max), apply(newhits[,6:11], 1, max)), main="SD vs Max Z-Scores of Old and New Hits", xlab="SD of Zscores", ylab="Max Zscore", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])))
#plot(c(apply(prevhits[,6:11], 1, sd), apply(newhits[,6:11], 1, sd), apply(failhits[,6:11], 1, sd)), c(apply(prevhits[,6:11], 1, max), apply(newhits[,6:11], 1, max), apply(failhits[,6:11], 1, max)), main="SD vs Max Z-Scores of Old and New Hits", xlab="SD of Zscores", ylab="Max Zscore", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])))

plot(c(apply(prevhits[,6:11], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(newhits[,6:11], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(failhits[,6:11], 1, function(x) { val1 <- sd(x/max(x)); return(val1);})), c(apply(prevhits[,6:11], 1, max), apply(newhits[,6:11], 1, max), apply(failhits[,6:11], 1, max)), main="SD(Zscore/MaxZ) vs MaxZ of Prev, New, & Failed Hits", xlab="SD(Zscore/MaxZ)", ylab="Max Zscore", pch=c(16,16,16), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5) 

legend(.225, 22.25, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))

#plot(c(prevhits$lbfav, newhits$lbfav), c(apply(prevhits[,6:11], 1, sd), apply(newhits[,6:11], 1, sd)), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])))
#plot(c(prevhits$lbfav, newhits$lbfav, failhits$lbfav), c(apply(prevhits[,6:11], 1, sd), apply(newhits[,6:11], 1, sd), apply(failhits[,6:11], 1, sd)), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])))

plot(c(apply(prevhits[,6:11], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(newhits[,6:11], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(failhits[,6:11], 1, function(x) { val1 <- sd(x/max(x)); return(val1);})), c(prevhits$lbfav, newhits$lbfav, failhits$lbfav), main="SD(Zscore/MaxZ) vs LBF of Prev, New, & Failed Hits", xlab="SD(Zscore/MaxZ)", ylab="LBF", ylim=c(2.5,100), pch=c(16,16,16), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5)

legend(.225, 100, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))

dev.off()








