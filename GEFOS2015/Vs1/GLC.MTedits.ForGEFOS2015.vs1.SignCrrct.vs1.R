set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GEFOS2015/Vs1/GEFOS2015.RSS0.vs1.SignCrrct.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GEFOS2015/Vs1/GEFOS2015.dtlesssignif.vs1.SignCrrct.vs1.annot.vs1.MAF.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0,]

Z = cbind(gl$Z.fa,gl$Z.fn,gl$Z.ls)

#~~~
> dim(gl)
[1] 2866   17
> gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
> dim(gl)
[1] 2866   17
> head(gl)
              snp chr       pos      maf     p_fa  n_fa         p_fn  n_fn
1 chr10:124064829  10 124064829 0.411456 0.049199 32965  0.000000914 32965
2  chr10:79459841  10  79459841 0.221480 0.185891 32965 -0.965944000 32965
3  chr11:46723937  11  46723937 0.394590 0.864396 32965 -0.000000318 32965
4  chr11:68177510  11  68177510 0.144710 0.007741 32965  0.000011000 32965
5  chr11:68181621  11  68181621 0.144630 0.006701 32965  0.000031700 32965
6  chr11:68243038  11  68243038 0.232776 0.005431 32965  0.001453000 32965
       p_ls  n_ls annot      Z.fa        Z.fn      Z.ls   mvstat       mvp
1  7.67e-07 32965     0  1.966863  4.90930443  4.943582 35.89473  7.103357
2  5.55e-09 32965     2  1.322833 -0.04269583  5.829781 41.94161  8.384235
3 -1.87e-04 32965     2 -0.170781  5.11245426  3.735961 29.98309  5.856545
4  8.47e-11 32965     2 -2.663163 -4.39652045 -6.492008 49.09713  9.905249
5  2.53e-11 32965     2 -2.711349 -4.16089703 -6.671617 50.47752 10.199204
6  2.05e-10 32965     2 -2.780292 -3.18391036 -6.357548 45.01141  9.036157
       unip
1  6.115205
2  8.255707
3  6.497573
4 10.072117
5 10.596879
6  9.688246
~~~

n = cbind(gl$n_fa, gl$n_fn, gl$n_ls)
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
[1] 30 18
#~~~
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
modelmatrix = data.frame(cbind(modelmatrix,cumsum(modelmatrix[,4])))
colnames(modelmatrix)= c("fa","fn","ls","p","cump")

allassoc=(apply((modelmatrix[,1:3]>0),1,sum)==3) #vector of which represent situations in which all 3 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:3]>0),1,sum)==2) #vector of which represent situations in which 2 phenotypes are associated

sum(modelmatrix[allassoc,4]) #0.9574025
sum(modelmatrix[allbut1assoc,4]) #0.03937951

#~~~
> sum(modelmatrix[allassoc,4]) #0.9574025
[1] 0.8665836
> sum(modelmatrix[allbut1assoc,4]) #0.03937951
[1] 0.1334164
#~~~

#look at weight on each of ls, fn, fa being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4]) 
# 8.189632e-20, 1.872029e-58, 0.1334164

#~~~
> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4])
[1] 8.189632e-20
> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4])
[1] 1.872029e-58
> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4])
[1] 0.1334164
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notls, or notfn or notfa)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"GEFOS2015.gl.bestmodel.vs2.SignCrrct.vs1.txt",quote=F,sep=" & ",row.names=F)

#~~~
> cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2))
     [,1]        [,2] [,3]  
[1,] "rs1357651" "4"  "0.63"
[2,] "rs2235811" "4"  "0.66"
[3,] "rs2566752" "4"  "0.77"

> modelmatrix
   fa fn ls             p      cump
1   2  1  1  5.855122e-01 0.5855122
2   1  1  1  2.810715e-01 0.8665836
3   0  1  1  1.334164e-01 1.0000000
4   1  1  0  8.189632e-20 1.0000000


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
#~~~
> min(lbf.av.glhits)
[1] 6.269104
#~~~


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]


#sub = gl[gl$nmin>50000,]
sub = gl
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
#newhits = sub[sub$annot==0 & sub$lbfav>6.269104 & sub$nmin>1000,c(1:3,7,18:20,25)]
newhits = sub[sub$annot==0 & sub$lbfav>6.269104 & sub$nmin>1000,c(1:4,11:14,19)]

#~~~
> dim(sub)
[1] 2866   19
> l=indephits(sub$lbfav,sub$chr,sub$pos)
> sub=sub[l==1,]
> dim(sub)
[1] 57 19
#~~~

#extract lbfs for all the new hits
#lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.bigmat
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>6.269104 & sub$nmin>1000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notls, or notfn or notfa)
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
tophits =  sub[sub$lbfav>3.323736 & sub$nmin>1000,]
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
write.table(file="GEFOS2015.newtophits.vs1.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>6.269104 & sub$nmin>1000 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>6.269104 & sub$nmin>1000 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)
#write.table(file="GEFOS2015.newtophits.vs1.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>5 & sub$nmin>1000 & sub$annot==0,c(1:4,11)],round(sub[sub$lbfav>5 & sub$nmin>1000 & sub$annot==0,c(12:14,19)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)

#newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>1000,c(1:3,7,18:20,25)]

#newhits = gl[gl$annot==0 & gl$lbfav>3.323736 & gl$nmin>1000,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)










