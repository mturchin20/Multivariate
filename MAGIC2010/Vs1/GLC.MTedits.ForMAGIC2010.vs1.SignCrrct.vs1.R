set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/MAGIC2010/Vs1/MAGIC2010.RSS0.vs1.SignCrrct.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/MAGIC2010/Vs1/MAGIC2010.dtlesssignif.vs1.SignCrrct.vs1.annot.vs1.MAF.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0,]

~~~
> dim(gl)
[1] 796  20
> gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
> dim(gl)
[1] 783  20

> dim(gl)
[1] 475  20
> gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
> 
> dim(gl)
[1] 475  20
~~~

Z = cbind(gl$Z.FastIns,gl$Z.FastGlu,gl$Z.HOMA_B,gl$Z.HOMA_IR)
n = cbind(gl$n_FastIns, gl$n_FastGlu, gl$n_HOMA_B, gl$n_HOMA_IR)
n = apply(n,1,min)
gl$nmin=n
#NOTE 20160225 -- running code with nmin=20001 produces almost essentially identical results
#gl$nmin=20001
#
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
colnames(modelmatrix)= c("FastIns","FastGlu","HOMA_B","HOMA_IR","p","cump")

allassoc=(apply((modelmatrix[,1:4]>0),1,sum)==4) #vector of which represent situations in which all 4 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:4]>0),1,sum)==3) #vector of which represent situations in which 3 phenotypes are associated

sum(modelmatrix[allassoc,5]) #0.9756873
sum(modelmatrix[allbut1assoc,5]) #0.02431274

~~~
> sum(modelmatrix[allassoc,5]) #0.9756873
[1] 4.269205e-11
> sum(modelmatrix[allbut1assoc,5]) #0.02431274
[1] 0.133333
~~~

~~~
> modelmatrix[1:20,]
   FastIns FastGlu HOMA_B HOMA_IR             p      cump
1        1       0      0       1  8.666670e-01 0.8666670
2        0       1      1       1  6.666663e-02 0.9333336
3        1       2      2       0  6.664448e-02 0.9999781
4        1       2      1       0  2.191608e-05 1.0000000
5        1       2      2       2  4.269205e-11 1.0000000
6        1       0      1       0  9.308232e-12 1.0000000
7        1       2      1       2  1.361671e-42 1.0000000
8        1       2      1       1  6.477864e-44 1.0000000
9        1       0      0       0  1.073483e-46 1.0000000
10       1       1      1       0  5.796651e-49 1.0000000
11       1       2      2       1  1.826810e-51 1.0000000
12       1       0      1       1  7.162499e-65 1.0000000
13       1       1      0       1  9.859988e-68 1.0000000
14       1       1      2       2  9.902631e-69 1.0000000
15       1       1      1       2  6.034307e-69 1.0000000
16       1       1      1       1  1.599572e-75 1.0000000
17       1       1      2       0  9.040324e-79 1.0000000
18       1       1      0       0 3.129001e-101 1.0000000
19       1       1      2       1 3.429365e-117 1.0000000
20       1       0      2       2 2.050915e-313 1.0000000
~~~

#look at weight on each of FastIns, FastGlu, HOMA_B, HOMA_IR being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5]) 
# 1.711418e-54, 6.298546e-39, 0.02431274, 2.91612e-55 

~~~
> sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5])
[1] 0.0666664
> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5])
[1] 9.859988e-68
> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5])
[1] 7.162499e-65
> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5])
[1] 0.06666663
~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notHOMA_IR, or notHOMA_B, or notFastGlu or notFastIns)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"gl.bestmodel.vs2.SignCrrct.vs1.txt",quote=F,sep=" & ",row.names=F)

#     [,1]         [,2] [,3]
# [1,] "rs17477177" "4"  "0.6"
> cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2))
     [,1] [,2] [,3]

#Below isn't expected to work since best model is generally only 2 phenotypes being associated, not 'max-1' which the above and below code assumes
> cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],apply(pp.classmatrix,2,max)[bestclass!=1])
      [,1]         [,2] [,3]                  
 [1,] "rs10830963" "4"  "7.33542677012644e-65"
 [2,] "rs10885122" "5"  "1.19209756074573e-11"
 [3,] "rs11071657" "5"  "1.17445183626875e-07"
 [4,] "rs11558471" "5"  "3.54034024688273e-12"
 [5,] "rs11605924" "5"  "1.31361397938308e-10"
 [6,] "rs11708067" "5"  "8.28608127670209e-10"
 [7,] "rs174550"   "5"  "2.54906379185451e-09"
 [8,] "rs2191349"  "5"  "1.06328695660413e-17"
 [9,] "rs340874"   "5"  "3.25871072682321e-08"
[10,] "rs35767"    "5"  "0.999999289928613"   
[11,] "rs4506565"  "5"  "2.23256660762174e-08"
[12,] "rs4607517"  "5"  "1.29807347756377e-35"
[13,] "rs560887"   "4"  "1.68805142271953e-65"
[14,] "rs780094"   "2"  "0.999995948998168"   
[15,] "rs7944584"  "5"  "1.12480466179796e-10"



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
#[1] 8.51211


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]


sub = gl[gl$nmin>20000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>8.51211 & sub$nmin>20000,c(4,2:3,7,22:25,30)]
newhits = sub[sub$annot==0 & sub$lbfav>8.51211 & sub$nmin>20000,]

~~~
> dim(sub)
[1] 475  22
> l=indephits(sub$lbfav,sub$chr,sub$pos)
> sub=sub[l==1,]
> dim(sub)
[1] 27 22
~~~

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>20000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>8.51211 & sub$nmin>20000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#9 in new results
#bestclass                    
> cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))
                   bestclass                       
 [1,] "rs10012946" "2"       "2.83440185124899e-06"
 [2,] "rs1015497"  "5"       "0.997224337672807"   
 [3,] "rs10768118" "5"       "0.984508385813887"   
 [4,] "rs10997908" "5"       "0.00521897281190426" 
 [5,] "rs12759838" "5"       "0.999991275765589"   
 [6,] "rs17390909" "5"       "2.9727699646322e-09" 
 [7,] "rs187678"   "5"       "0.938855063939861"   
 [8,] "rs4912494"  "5"       "0.999915595321191"   
 [9,] "rs7708079"  "5"       "0.998925790054594"   


pdf("plots.bychr.vs2.pdf")
for(i in 1:22){
 if(sum(sub$chr==i)>0){
plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>50000))
abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>8.51211],col=3)
}
}
dev.off()

#try to plot all hits to see relationships?
tophits =  sub[sub$lbfav>8.51211 & sub$nmin>20000,]
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

write.table(file="newtophits.vs2.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>8.51211 & sub$nmin>20000 & sub$annot==0,c(1:3,13)],round(sub[sub$lbfav>8.51211 & sub$nmin>20000 & sub$annot==0,c(4,15:22)],digits=5)),quote=FALSE,sep= " ", row.names=FALSE)

c(4,2:3,7,22:25,30)

#newhits = gl[gl$annot==0 & gl$lbfav>8.51211 & gl$nmin>20000,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)




