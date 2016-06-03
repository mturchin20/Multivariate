set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/RSS0.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/dtlesssignif.annot.txt", header=T)
Z = cbind(gl$Z.tc,gl$Z.tg,gl$Z.hdl,gl$Z.ldl)

n = cbind(gl$n_TC, gl$n_TG, gl$n_HDL, gl$n_LDL)
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
colnames(modelmatrix)= c("TC","TG","HDL","LDL","p","cump")

allassoc=(apply((modelmatrix[,1:4]>0),1,sum)==4) #vector of which represent situations in which all 4 phenotypes are associated
allbut1assoc=(apply((modelmatrix[,1:4]>0),1,sum)==3) #vector of which represent situations in which 3 phenotypes are associated

sum(modelmatrix[allassoc,5]) #0.77
sum(modelmatrix[allbut1assoc,5]) #0.21

#look at weight on each of LDL, HDL, TG, TC being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5]) 
# 0.08, 0.07, 0.05, 0.00

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
gg=gl.glhits$gene
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"gl.bestmodel.txt",quote=F,sep=" & ",row.names=F)


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
layout(rbind(c(1,2)))
plot(lbf.uni.glhits,lbf.av.origprior.glhits,xlim=c(4,15),ylim=c(4,15))
abline(a=0,b=1)
plot(lbf.uni.glhits,lbf.all.glhits,xlim=c(4,15),ylim=c(4,15))
abline(a=0,b=1)


#this establishes a threshold for subsequent analysis
min(lbf.av.glhits)
#[1] 4.349679


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
o = order(gl$chr, gl$pos)
gl = gl[o,]


sub = gl[gl$nmin>50000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))
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

head(modelmatrix)
modelmatrix.paste <- apply(modelmatrix[,1:4], 1, function(x) paste(x, collapse="_"))
table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2010/Vs1/GLC.MTedits.ForGlobalLipids2010.vs1.ForCSHLPoster.ModelsVslbfavNewOld.vs1.jpeg", res=300, width=5000, height=4000)
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
#TC TG HDL LDL          p      cump
#1  1  1   1   2 0.38057200 0.3805720
#2  1  2   1   2 0.07952934 0.4601013
#3  2  1   1   2 0.06249012 0.5225915
#4  2  1   1   0 0.04642033 0.5690118
#5  1  1   2   2 0.04622803 0.6152398
#6  2  2   2   1 0.04577695 0.6610168
#> modelmatrix.paste <- apply(modelmatrix[,1:4], 1, function(x) paste(x, collapse="_"))
#> table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
#
#1_1_1_1 1_1_1_2 1_2_0_2 1_2_2_2 2_1_1_0 2_1_1_2 2_2_1_0 2_2_2_1 
#1      10       1       2       1       1       1       2 
#> 
#> topmodels
#i                                                      
#[1,] "1_1_1_2" "10" "0.648146156223119" "0.208644055894459" 
#[2,] "2_2_2_1" "2"  "0.672982625286163" "0.0729814274809761"
#[3,] "1_2_2_2" "2"  "0.293631952449074" "0.0277544774992512"
#[4,] "1_1_1_1" "1"  "0.999979880081639" NA                  
#[5,] "2_1_1_2" "1"  "0.497955528594509" NA                  
#[6,] "2_2_1_0" "1"  "0.305787348808541" NA                  
#[7,] "2_1_1_0" "1"  "0.796504451183269" NA                  
#[8,] "1_2_0_2" "1"  "0.380854602394163" NA                  
#
#[1,] "0.48043327407024"  "0.265161498069431"  "0.380571997601203" 
#[2,] "0.47298748864434"  "0.117000572799571"  "0.0457769477606886"
#[3,] "0.092198030260899" "0.0270997987484532" "0.0359827597613109"
#[4,] "0.999963268766235" NA                   "0.0323155579074516"
#[5,] "0.026714581058676" NA                   "0.0624901203276318"
#[6,] "0.079570479688507" NA                   "0.0239889845317868"
#[7,] "0.697509065664131" NA                   "0.0464203307564649"
#[8,] "0.141905208654991" NA                   "0.0341502727063856"
#> oldtopmodels
#i                                                       
#[1,] "1_1_1_2" "56" "0.594916963550666" "0.218761856241571"  
#[2,] "2_2_2_1" "7"  "0.57629648196499"  "0.130271117017175"  
#[3,] "1_2_2_2" "6"  "0.323555221970788" "0.0333562926025526" 
#[4,] "1_2_0_2" "5"  "0.372293353897778" "0.122804033760489"  
#[5,] "2_1_1_0" "5"  "0.653788922243709" "0.312111918832214"  
#[6,] "2_2_1_1" "3"  "0.459426311978285" "0.242474358430201"  
#[7,] "1_2_1_1" "3"  "0.734567634682521" "0.238856238009038"  
#[8,] "2_1_1_2" "3"  "0.299405435927484" "0.161772040415935"  
#[9,] "1_0_1_2" "2"  "0.792445237351275" "0.0255341303056004" 
#[10,] "2_2_1_0" "2"  "0.321296889334506" "0.0359320469804227" 
#[11,] "1_1_1_1" "2"  "0.99791539573249"  "0.00294517167466403"
#[12,] "1_0_1_1" "2"  "0.637291355213706" "0.0626930594276656" 
#[13,] "1_1_0_2" "1"  "0.91631477029998"  NA                   
#[14,] "2_1_1_1" "1"  "0.945410424475398" NA                   
#[15,] "2_1_0_2" "1"  "0.834210124574721" NA                   
#[16,] "2_1_2_0" "1"  "0.763202932726069" NA                   
#[17,] "1_2_1_2" "1"  "0.517542986147095" NA                   
#[18,] "2_0_1_0" "1"  "0.308336541599954" NA                   
#
#[1,] "0.382516030096653" "0.309663529764632"   "0.380571997601203"  
#[2,] "0.317465686901165" "0.152657307545512"   "0.0457769477606886" 
#[3,] "0.104021775649593" "0.0660443074047328"  "0.0359827597613109" 
#[4,] "0.18093631170802"  "0.13988171072146"    "0.0341502727063856" 
#[5,] "0.498515028803513" "0.454165003831325"   "0.0464203307564649" 
#[6,] "0.224036798194631" "0.321585390235337"   "0.0440613699607882" 
#[7,] "0.596960528892902" "0.258304628897546"   "0.0249532653600498" 
#[8,] "0.005605479043644" "0.00356398169596546" "0.0624901203276318" 
#[9,] "0.671208842599156" "0.0146094757861486"  "0.0283173783453159" 
#[10,] "0.103419202899761" "0.0865284885543862"  "0.0239889845317868" 
#[11,] "0.995831218079251" "0.00589093440479308" "0.0323155579074516" 
#[12,] "0.274921720704251" "0.125864941123202"   "0.0147643986093055" 
#[13,] "0.87673522796214"  NA                    "0.0173738961700482" 
#[14,] "0.901705149104183" NA                    "0.0199074856276527" 
#[15,] "0.746414330019498" NA                    "0.00857738046841892"
#[16,] "0.543731464591985" NA                    "0.0125643058985003" 
#[17,] "0.131306484791022" NA                    "0.0795293417741806" 
#[18,] "0.031936972762708" NA                    "0.0135778437323181" 
#> summary(lm(apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav)
#
#Residuals:
#Min      1Q  Median      3Q     Max
#-0.3735 -0.1852  0.0153  0.1715  0.3884
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.88886    0.35610   2.496   0.0231 *
#newhits$lbfav -0.05408    0.06586  -0.821   0.4229
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2334 on 17 degrees of freedom
#Multiple R-squared:  0.03815,   Adjusted R-squared:  -0.01843
#F-statistic: 0.6742 on 1 and 17 DF,  p-value: 0.4229
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
#-0.41880 -0.26112 -0.03008  0.22763  0.57409
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)    0.68791    0.47270   1.455    0.164
#newhits$lbfav -0.05110    0.08743  -0.584    0.567
#
#Residual standard error: 0.3099 on 17 degrees of freedom
#Multiple R-squared:  0.0197,    Adjusted R-squared:  -0.03796 
#F-statistic: 0.3416 on 1 and 17 DF,  p-value: 0.5666
#
#> summary(lm(apply(ppmatrix, 2, max) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix, 2, max) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.38469 -0.21113 -0.01924  0.19129  0.42517
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               5.730e-01  2.633e-02  21.761   <2e-16 ***
#gl[gl$annot == 1, ]$lbfav 7.634e-05  4.152e-04   0.184    0.854
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2341 on 100 degrees of freedom
#Multiple R-squared:  0.000338,  Adjusted R-squared:  -0.009659
#F-statistic: 0.03381 on 1 and 100 DF,  p-value: 0.8545
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
#-0.36285 -0.28038 -0.07511  0.25217  0.63398
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.3613734  0.0356241  10.144   <2e-16 ***
#gl[gl$annot == 1, ]$lbfav 0.0001976  0.0005617   0.352    0.726
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.3167 on 100 degrees of freedom
#Multiple R-squared:  0.001236,  Adjusted R-squared:  -0.008751
#F-statistic: 0.1238 on 1 and 100 DF,  p-value: 0.7257
#
#>
#~~~



pdf("plots.bychr.pdf")
for(i in 1:22){
 if(sum(sub$chr==i)>0){
plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>50000))
abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>4.3],col=3)
}
}
dev.off()

#try to plot all hits to see relationships?
tophits =  sub[sub$lbfav>4.3 & sub$nmin>50000,]
betahat = tophits[,c(6,9,12,15)]
betahat.scaled = betahat/apply(betahat,1,sd)
betahat.scaled.pr = prcomp(betahat.scaled,center=F)
plot(abs(betahat.scaled.pr$x[,1]),abs(betahat.scaled.pr$x[,2]),type="n")
text(abs(betahat.scaled.pr$x[,1]),abs(betahat.scaled.pr$x[,2]),as.character(tophits$chr))

hla = gl[gl$chr==6 & gl$pos<40e6 & gl$pos>28e6,]
hla.z = hla[,20:23]

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

zhla = head(gl[gl$chr==6 & gl$pos<40e6 & gl$pos>28e6,],600)[,20:23]
image(cor(t(zhla))^2)

plot(gl$pos[gl$chr==1],gl$lbfav[gl$chr==1],ylim=c(0,10),xlim=c(27102620-10^5,27102620+10^5))

write.table(file="newtophits.txt",cbind(sub[sub$lbfav>4.3 & sub$nmin>50000 & sub$annot==0,1:4],round(sub[sub$lbfav>4.3 & sub$nmin>50000 & sub$annot==0,c(20:23,28)],digits=1)),quote=FALSE,sep= " & ", row.names=FALSE)

newhits = gl[gl$annot==0 & gl$lbfav>4.3 & gl$nmin>50000,c(1:3,28)]
write.table(file="newhits.wr",newhits,quote=FALSE,row.names=FALSE)

#             snp chr       pos   maf   Z.tg   Z.tc  Z.ldl  Z.hdl    lbfav
# 1612 rs12739698   1  27102620 0.083  3.315  3.394  4.626 -5.096 7.186728 #8kb 5' upstream NROB2
# 3683   rs267733   1 149225460 0.141 -0.970  3.469  5.404 -1.901 4.821367 #ANXA9, coding non-synon
# 
# 245  rs10490632   2 118295555 0.082  0.845  4.718  5.384 -1.987 5.111464 #intronic, DDX18
# 
# 1837 rs13326165   3  52507158 0.207 -4.375 -1.390 -2.053  5.004 4.740907 #STAB1 intronic
# 
# 7061   rs762861   4   3411809 0.261  4.543  4.909  4.380 -2.102 5.241455 #HGFAC 5'upstream, RGS12 3' downstream 
# 
# 8062   rs998584   6  43865874 0.491  5.066  0.334 -0.032 -4.299 4.458487 #VEGFA 3'downstream (4kb)
# 
# 6476  rs6951245   7   1024719 0.156  1.980  5.417  3.561  3.254 5.410647 #C7ORF50 intronic
# 5112  rs4722551   7  25958351 0.177  3.921 -3.792 -4.891 -1.330 7.509293 #miR-148a
# 
# 2302 rs17134533  10   5237098 0.146 -4.919 -3.235 -1.163 -1.817 6.011030 #AKR1C4 intronic
# 572  rs10904908  10  17300296 0.432 -1.567 -4.917 -3.434 -3.640 4.649076 # 10kb 5' upstream of VIM
# 7949   rs970548  10  45333283 0.246 -0.026 -3.750 -2.097 -5.245 5.448010 #MARCH8 intronic
# 1947  rs1408579  10 101902184 0.470  2.599  3.242  0.760  3.930 4.829068 #ERLIN1 intronic
# 
# 
# 853  rs11246602  11  51368666 0.126 -0.168 -2.536 -0.708 -5.395 5.289490
# 839  rs11229252  11  54886216 0.091  0.538  3.373  1.540  4.966 4.932485
# 836  rs11227638  11  55776161 0.118 -0.637 -2.285 -0.445 -5.026 4.772437 #OR5T3 #these 3 are in an area of high LD; lots of ORs
# 5493   rs499974  11  75132669 0.175 -2.098 -2.383  0.266 -4.110 4.439553 #16kb 5' upstream DGAT2
# 
# 5435  rs4942505  13  31861707 0.476 -1.346 -4.533 -5.330  3.092 5.829394 #BRCA2 intronic
# 
# 203  rs10422101  19  57011927 0.265  0.950 -3.608 -2.437 -5.411 5.602997 #intronic, FPR3
# 
# 3792  rs2746150   6  29550680 0.089 -3.764 -2.943 -1.044 -2.915 5.249762 #MAS1LP1; looks like it is due to HLA

#NROB2 (small heterodimer partner) regulates metabolic pathways, including hepatic bile acid, lipid, and glucose homeostasis. \cite{huang:2007}

#ANXA9 codes for the Annexin-A9 protein. The annexins are a family of calcium-dependent phospholipid-binding proteins, and studies in cows have associated variation in or near ANXA9 with milk-fat yield \cite{martinez-royo:2010}.

#STAB1 (also known as FEEL-1, and CLEVER-1) codes for the protein stabilin-1 which acts as a scavenger receptor for acetylated low density lipoprotein and oxidized LDL \cite{li:2011,kzhyshkowska:2006, adachi:2002}.

#VEGFA codes for the protein "Vascular endothelial growth factor A", and variants near VEGFA have been implicated in a range of clinical conditions, including diabetic retinopathy \cite{al-kateb:2007,abhary:2009,yang:2011} and age-related macular degeneration \cite{yu:2011}.

#AKR1C4 codes for the enzyme Aldo-keto reductase family 1 member C4, and plays a major role in bile acid biosynthesis \cite{russell:2003}, which is a major pathway of cholesterol catabolism in mammals.


#VIM codes for Vimentin, which assists in the transport of LDL cholesterol from a lysosome to the site of esterification \cite{sarria:1992}.

#DGAT2 encodes one of two enzymes which catalyzes the final reaction in the synthesis of triglycerides, has been implicated as a major target for the action of niacin in regulating lipids \cite{ganji:2004,hu:2012}.


#E3 ubiquitin-protein ligase MARCH8. snps in it rs11239550
#have been associated with mean corpuscular volume in a large meta-analysis \cite{ganesh:2009}


#ERLIN1 codes for the protein Erlin-1, which is a member of the prohibitin family of proteins that define lipid-raft-like domains of the endoplasmic reticulum \cite{browman:2006} and SNPs near ERLIN1 have been associated with  plasma levels of alanine-aminotransferase \cite{yuan:2008}, an important liver enzyme. 


