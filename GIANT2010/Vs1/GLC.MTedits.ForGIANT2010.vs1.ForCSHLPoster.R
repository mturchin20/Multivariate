set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.RSS0.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
#gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$chrbp),]
Z = cbind(gl$Z.height,gl$Z.BMI,gl$Z.WHRadjBMI)

#~~~
#> gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz", header=T)
#> dim(gl)
#[1] 8841   23
#> gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$chrbp),]
#> dim(gl)
#[1] 7988   23
#~~~

n = cbind(gl$n_height, gl$n_BMI, gl$n_WHRadjBMI)
n = apply(n,1,min)
gl$nmin=n
sigmaa=c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)
lbf=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa)
lbf.bigmat=do.call(rbind,lbf$lbf) #lbf$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows 

#collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
collapse = function(x,nsigmaa){return(apply(matrix(x,ncol=nsigmaa),1,sum))}
  
#compute empirical bayes priors from GWAS hit SNPs
gl.glhits = gl[gl$annot==1,] # subset of GWAS SNPs
#CHECK_0: Note -- some top hits dropped because they did not have 'maf'; expecting ~226 below, not 176
#~~~
#> dim(gl.glhits)
#[1] 176  24
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
colnames(modelmatrix)= c("height","BMI","WHRadjBMI","p","cump")

allassoc=(apply((modelmatrix[,1:3]>0),1,sum)==3) #vector of which represent situations in which all 3 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:3]>0),1,sum)==2) #vector of which represent situations in which 2 phenotypes are associated

sum(modelmatrix[allassoc,4]) #0.7596256
sum(modelmatrix[allbut1assoc,4]) #0.2315347

#~~~
#> sum(modelmatrix[allassoc,4])
#[1] 0.7596256
#> sum(modelmatrix[allbut1assoc,4])
#[1] 0.2315347
#~~~

#look at weight on each of WHRadjBMI, BMI, height being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4]) 
# 0.2259918, 0.003215401, 0.002327482 

#~~~
#> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4])
#[1] 0.2259918
#> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4])
#[1] 0.003215401
#> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4])
#[1] 0.002327482
#~~~

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notWHRadjBMI, or notBMI or notheight)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"gl.bestmodel.vs2.txt",quote=F,sep=" & ",row.names=F)

#       Empty for height -- all SNPs have 'bestclass' as 1

#Below are results from Lipids just as a reminder
#     [,1]         [,2] [,3]
#     [1,] "rs838880"   "2"  "0.75"
#     [2,] "rs386000"   "2"  "0.64"
#     [3,] "rs1800961"  "4"  "0.87"
#     [4,] "rs1532085"  "2"  "0.54"
#     [5,] "rs16942887" "2"  "0.91"
#     [6,] "rs737337"   "2"  "0.55"

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
#~~~
#> min(lbf.av.glhits)
#[1] 3.323736
#~~~


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]


sub = gl[gl$nmin>50000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>50000,c(1:3,7,18:20,25)]

#~~~
#> sub = gl[gl$nmin>50000,]
#> dim(sub)
#[1] 7977   25
#> sub=sub[l==1,]
#> dim(sub)
#[1] 189  25
#~~~

#Investigating 'indephits'
#~~~
#> sub$lbfav[sub$chr==2 & sub$pos<900000 & sub$pos>100000]
#[1]  4.100670 17.314986  5.286055 17.312360  5.565013 17.737481 17.679755
#[8] 17.760487 17.212332 17.249917 17.518145  5.390687 17.218923 17.571293
#[15]  5.138245 18.073186 17.368890  4.921732 17.476747 17.009358 13.228503
#[22] 17.632784 17.501315 17.105200  5.131919 17.194006 17.677015 17.751375
#[29] 17.445604 17.734091 17.144550 17.831652 17.498369 17.948064 17.759142
#[36]  7.862099 17.289904 17.514460 17.683073 15.033423  5.230502  5.159357
#[43] 17.505760  5.091114  5.309113  5.744831 15.204181 17.379275 16.634037
#[50] 17.875640 17.542579  5.105418 17.346049 17.101164 17.661485  5.423729
#[57] 17.450871 17.504994 17.737586 17.760496 17.650788  5.090595 17.246934
#[64] 17.404437 17.455963 18.043795 17.647656 17.164451 17.690611 17.502048
#[71] 17.242195 17.534686  5.489030 17.264382 16.731337 17.912680 17.259981
#[78]  6.099039 17.190026 17.396661 17.388246 17.415921 17.631646 17.408549
#[85] 17.866815 17.575844  5.289236 17.336140 18.113258  5.210825  3.121887
#[92]  3.609672 17.084093 17.612929 17.837278  5.110351 17.323426 17.404353
#[99]  5.089691  5.550603 17.970325  5.087274 17.453137  5.225956  5.189750
#[106] 17.461156 17.934960 17.408513 17.220760  5.090597 14.250492 17.682808
#[113] 17.706136
#> max(sub$lbfav[sub$chr==2 & sub$pos<900000 & sub$pos>100000])
#[1] 18.11326
#~~~

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>50000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notWHRadjBMI, or notBMI or notheight)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#19 in results
#> cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))
#bestclass
#[1,] "rs17016663" "1"       "0.943259069758143"
#[2,] "rs10805383" "1"       "0.937650792648966"
#[3,] "rs7601531"  "1"       "0.733867477052145"
#[4,] "rs11835818" "1"       "0.907304645214397"
#[5,] "rs1809889"  "1"       "0.637590556831822"
#[6,] "rs12204421" "1"       "0.644883927737445"
#[7,] "rs2025151"  "1"       "0.659267026848644"
#[8,] "rs34651"    "1"       "0.656838708886582"
#[9,] "rs6824258"  "1"       "0.691982088687886"
#[10,] "rs4735692"  "1"       "0.923288493973102"
#[11,] "rs7614120"  "1"       "0.905348670213562"
#[12,] "rs7081678"  "1"       "0.954339521237619"
#[13,] "rs12534698" "1"       "0.645066623519918"
#[14,] "rs648831"   "1"       "0.635682512834355"
#[15,] "rs2390312"  "1"       "0.911776480593067"
#[16,] "rs389883"   "1"       "0.990111195777526"
#[17,] "rs17783015" "1"       "0.674581804925153"
#[18,] "rs10040888" "1"       "0.642088830668284"
#[19,] "rs10822129" "1"       "0.937243836077944"

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GLC.MTedits.ForGIANT2010.vs1.ForCSHLPoster.ModelsVslbfavNewOld.vs1.jpeg", res=300, width=5000, height=4000)
par(mfrow=c(2,2))

head(modelmatrix)
modelmatrix.paste <- apply(modelmatrix[,1:3], 1, function(x) paste(x, collapse="_"))
table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),4]

snpmodels <- cbind(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)], apply(ppmatrix.newhits, 2, max), apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] }), apply(ppmatrix.newhits, 2, function(x) { sort(x,partial=length(x)-2)[length(x)-2] }))
topmodels <- c()

for (i in unique(snpmodels[,1])) {
	snpmodels.subset <- matrix(snpmodels[snpmodels[,1]==i,], nrow=length(which(snpmodels[,1]==i)))
	topmodels <- rbind(topmodels, cbind(i, nrow(snpmodels.subset), mean(as.numeric(snpmodels.subset[,2])), sd(as.numeric(snpmodels.subset[,2])), mean(as.numeric(snpmodels.subset[,2]) - as.numeric(snpmodels.subset[,3])), sd(as.numeric(snpmodels.subset[,2]) - as.numeric(snpmodels.subset[,3])), modelmatrix[modelmatrix.paste == i,4]))
}
topmodels <- topmodels[order(as.numeric(topmodels[,2]), decreasing=TRUE),]

#topmodels <- cbind(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])), table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]), modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),4])
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
	oldtopmodels <- rbind(oldtopmodels, cbind(i, nrow(oldsnpmodels.subset), mean(as.numeric(oldsnpmodels.subset[,2])), sd(as.numeric(oldsnpmodels.subset[,2])), mean(as.numeric(oldsnpmodels.subset[,2]) - as.numeric(oldsnpmodels.subset[,3])), sd(as.numeric(oldsnpmodels.subset[,2]) - as.numeric(oldsnpmodels.subset[,3])), modelmatrix[modelmatrix.paste == i,4]))
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
#height BMI WHRadjBMI          p      cump
#1      1   2         2 0.31153912 0.3115391
#2      1   1         2 0.22311231 0.5346514
#3      1   2         0 0.20450584 0.7391573
#4      1   2         1 0.15293796 0.8920952
#5      2   1         2 0.04854781 0.9406430
#6      1   1         0 0.02142094 0.9620640
#> modelmatrix.paste <- apply(modelmatrix[,1:3], 1, function(x) paste(x, collapse="_"))
#> table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
#
#1_1_1 1_1_2 1_2_1 1_2_2 
#2     5     2    10 
#> modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),4]
#[1] 0.01703055 0.22311231 0.15293796 0.31153912
#> topmodels
#i                                                                        
#[1,] "1_2_2" "10" "0.475660686616409" "0.0450199468671559" "0.158555046638225"
#[2,] "1_1_2" "5"  "0.735493411817986" "0.149817380855992"  "0.545307527610242"
#[3,] "1_2_1" "2"  "0.694340612082986" "0.0785915217494392" "0.549443823371466"
#[4,] "1_1_1" "2"  "0.750354939914833" "0.153260740201386"  "0.577754171709174"
#
#[1,] "0.0297205116348298" "0.311539120527113" 
#[2,] "0.291206347124426"  "0.223112308215958" 
#[3,] "0.114740453845823"  "0.152937962820109" 
#[4,] "0.234455677190261"  "0.0170305459846473"
#> oldtopmodels
#i                                                                         
#[1,] "1_2_2" "103" "0.472674743310942" "0.0544671658388599" "0.157693512439468"
#[2,] "1_1_2" "39"  "0.632941577213098" "0.180299470254791"  "0.39901833878412" 
#[3,] "1_2_1" "27"  "0.685578426667608" "0.223675256723207"  "0.532613930616903"
#[4,] "2_1_2" "6"   "0.546919083579652" "0.0658991740796973" "0.177176118808093"
#[5,] "1_1_1" "1"   "0.631953581464673" NA                   "0.440453072824029"
#
#[1,] "0.0333288298072369" "0.311539120527113" 
#[2,] "0.303781911901942"  "0.223112308215958" 
#[3,] "0.337531278805877"  "0.152937962820109" 
#[4,] "0.146081668617743"  "0.048547814103867" 
#[5,] NA                   "0.0170305459846473"
#> summary(lm(apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.21174 -0.10609 -0.08813  0.11194  0.28773
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.54000    0.11776   4.586 0.000263 ***
#newhits$lbfav  0.01047    0.02091   0.501 0.623096
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.1609 on 17 degrees of freedom
#Multiple R-squared:  0.01453,   Adjusted R-squared:  -0.04344
#F-statistic: 0.2506 on 1 and 17 DF,  p-value: 0.6231
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
#-0.2560 -0.1845 -0.1665  0.2188  0.4551
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.27752    0.18997   1.461    0.162
#newhits$lbfav  0.01273    0.03373   0.377    0.711
#
#Residual standard error: 0.2596 on 17 degrees of freedom
#Multiple R-squared:  0.008309,  Adjusted R-squared:  -0.05003
#F-statistic: 0.1424 on 1 and 17 DF,  p-value: 0.7105
#
#> summary(lm(apply(ppmatrix, 2, max) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix, 2, max) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.26973 -0.08500 -0.04026  0.01555  0.41374
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.521457   0.018881   27.62   <2e-16 ***
#gl[gl$annot == 1, ]$lbfav 0.002508   0.001629    1.54    0.125
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.1551 on 174 degrees of freedom
#Multiple R-squared:  0.01344,   Adjusted R-squared:  0.007772
#F-statistic: 2.371 on 1 and 174 DF,  p-value: 0.1254
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
#-0.27731 -0.12260 -0.09724  0.01159  0.66229
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.246308   0.029926   8.230 4.24e-14 ***
#gl[gl$annot == 1, ]$lbfav 0.002707   0.002581   1.049    0.296
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2458 on 174 degrees of freedom
#Multiple R-squared:  0.006282,  Adjusted R-squared:  0.0005709
#F-statistic:   1.1 on 1 and 174 DF,  p-value: 0.2957
#
#>
#~~~





pdf("plots.bychr.vs2.pdf")
for(i in 1:22){
 if(sum(sub$chr==i)>0){
plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>50000))
abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>4.608818],col=3)
}
}
dev.off()

#try to plot all hits to see relationships?
tophits =  sub[sub$lbfav>3.323736 & sub$nmin>50000,]
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

write.table(file="GIANT2010.newtophits.vs1.txt",cbind(sub[sub$lbfav>3.323736 & sub$nmin>50000 & sub$annot==0,c(1:3,7)],round(sub[sub$lbfav>3.323736 & sub$nmin>50000 & sub$annot==0,c(18:20,25)],digits=1)),quote=FALSE,sep= " ", row.names=FALSE)

#newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>50000,c(1:3,7,18:20,25)]

#newhits = gl[gl$annot==0 & gl$lbfav>3.323736 & gl$nmin>50000,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)










