library(ggplot2)
library(reshape2)
set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GlobalLipids2013.RSS0.vs1.SignCrrct.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GlobalLipids2013.dtlesssignif.vs1.SignCrrct.vs1.annot.MAF.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0,]
#Z = cbind(gl$Z.tc,gl$Z.tg,gl$Z.hdl,gl$Z.ldl)
Z = cbind(gl$Z.ldl, gl$Z.hdl, gl$Z.tc, gl$Z.tg)

#> VYY
#           Z.LDL       Z.HDL       Z.TC       Z.TG
#[1,]  1.00000000 -0.05101127 0.83617727  0.1446093
#[2,] -0.05101127  1.00000000 0.09840344 -0.2298445
#[3,]  0.83617727  0.09840344 1.00000000  0.3389824
#[4,]  0.14460929 -0.22984452 0.33898241  1.0000000



#n = cbind(gl$n_TC, gl$n_TG, gl$n_HDL, gl$n_LDL)
n = cbind(gl$n_LDL, gl$n_HDL, gl$n_TC, gl$n_TG)
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

sum(modelmatrix[allassoc,5]) #0.9300127
sum(modelmatrix[allbut1assoc,5]) #0.05885692

#> sum(modelmatrix[allassoc,5]) #0.9300127
#[1] 0.3659239
#> sum(modelmatrix[allbut1assoc,5]) #0.05885692
#[1] 0.3117781
 
#> sum(modelmatrix[allassoc,5]) #0.9300127
#[1] 0.6147474
#> sum(modelmatrix[allbut1assoc,5]) #0.05885692
#[1] 0.2430032


#look at weight on each of LDL, HDL, TG, TC being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5]) 
# 0.04843585, 1.852496e-05, 0.01040254, 1.343718e-12

#> sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5])
#[1] 1.571201e-18
#> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5])
#[1] 0.1866107
#> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5])
#[1] 0.09731063
#> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5])
#[1] 0.02785675

#> sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5])
#[1] 0.09227745
#> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5])
#[1] 0.0587368
#> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5])
#[1] 0.06043325
#> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5])
#[1] 0.03155571


#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
#Best class is '1' for all SNPs?
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"GlobalLipids2013.gl.bestmodel.vs3.SignCrrct.vs1.txt",quote=F,sep=" & ",row.names=F)

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
#[1] 4.608818
#[1] 7.312795
[1] 3.989499



lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]


sub = gl[gl$nmin>50000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
#newhits = sub[sub$annot==0 & sub$lbfav>4.608818 & sub$nmin>50000,c(4,2:3,7,22:25,30)]
newhits = sub[sub$annot==0 & sub$lbfav>3.989499 & sub$nmin>50000,c(4,2:3,7,22:25,30)]

~~~
> dim(newhits)
[1] 92  9
~~~

#~~~
#> l=indephits(sub$lbfav,sub$chr,sub$pos)
#> length(l)
#[1] 11417
#> dim(sub)
#[1] 11417    30
#> sub=sub[l==1,]
#> #newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
#> newhits = sub[sub$annot==0 & sub$lbfav>4.608818 & sub$nmin>50000,c(4,2:3,7,22:25,30)]
#> dim(sub)
#[1] 268  30
#~~~

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>3.989499 & sub$nmin>50000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#84 in new results

#74 in new results
#bestclass                    
#[1,] "rs9856765"  "1"       "0.999981661378945"
#[2,] "rs7033354"  "1"       "0.999999651388728"
#[3,] "rs4808802"  "1"       "0.998392155719864"
#[4,] "rs2201637"  "1"       "0.798969363078283"
#[5,] "rs9646133"  "1"       "0.999918450115978"
#[6,] "rs10733608" "1"       "0.976641194419346"
#[7,] "rs11040329" "1"       "0.664682677261752"
#[8,] "rs211617"   "1"       "0.999985276934298"
#[9,] "rs6968554"  "1"       "0.995464675137669"
#[10,] "rs7076938"  "1"       "0.97713506628236" 
#[11,] "rs10757056" "1"       "0.995494084318009"
#[12,] "rs5752792"  "1"       "0.99413629689843" 
#[13,] "rs1688030"  "1"       "0.999189887866103"
#[14,] "rs17630235" "1"       "0.997715650680412"
#[15,] "rs550136"   "1"       "0.999275070485492"
#[16,] "rs721772"   "1"       "0.992654628495879"
#[17,] "rs4976033"  "1"       "0.993951799460061"
#[18,] "rs1482852"  "1"       "0.993277552876587"
#[19,] "rs2275774"  "1"       "0.999967704377268"
#[20,] "rs661171"   "1"       "0.982851269201709"
#[21,] "rs10496123" "1"       "0.993990119565431"
#[22,] "rs2278093"  "1"       "0.999920112253773"
#[23,] "rs11987974" "1"       "0.992036366229579"
#[24,] "rs4683438"  "1"       "0.999967165871738"
#[25,] "rs8069974"  "1"       "0.999985463967641"
#[26,] "rs11556341" "1"       "0.995595440355304"
#[27,] "rs10861661" "1"       "0.993455628244922"
#[28,] "rs2521921"  "1"       "0.995076925013423"
#[29,] "rs13379043" "1"       "0.999834188347365"
#[30,] "rs17309825" "1"       "0.975918361725404"
#[31,] "rs6822892"  "1"       "0.992019292195826"
#[32,] "rs4728614"  "1"       "0.999988384696052"
#[33,] "rs10832027" "1"       "0.997915468524301"
#[34,] "rs1872167"  "2"       "0.564635523180757"
#[35,] "rs2820426"  "1"       "0.999785060604486"
#[36,] "rs4988235"  "1"       "0.998665643961838"
#[37,] "rs10513688" "1"       "0.999915735396396"
#[38,] "rs12602912" "1"       "0.999947424372267"
#[39,] "rs895954"   "1"       "0.991227660321327"
#[40,] "rs2327277"  "1"       "0.999911309627049"
#[41,] "rs4871137"  "1"       "0.853489124696034"
#[42,] "rs7255743"  "1"       "0.999327714298719"
#[43,] "rs2845885"  "1"       "0.988970446492519"
#[44,] "rs1132990"  "1"       "0.999922339498545"
#[45,] "rs4897361"  "1"       "0.99980837739916" 
#[46,] "rs74458891" "1"       "0.889594973919558"
#[47,] "rs2247056"  "1"       "0.999999991989186"
#[48,] "rs4850047"  "1"       "0.997102669241985"
#[49,] "rs1174604"  "1"       "0.99999732835551" 
#[50,] "rs176813"   "1"       "0.999803882681035"
#[51,] "rs583484"   "1"       "0.999930295118799"
#[52,] "rs1062219"  "1"       "0.996858085577632"
#[53,] "rs6059932"  "1"       "0.996555190111254"
#[54,] "rs12133576" "1"       "0.995463958692015"
#[55,] "rs2746150"  "1"       "0.993118839632978"
#[56,] "rs13194504" "1"       "0.981868259064509"
#[57,] "rs2454722"  "1"       "0.913791896472197"
#[58,] "rs884366"   "1"       "0.856702579683981"
#[59,] "rs10408163" "1"       "0.995031140417091"
#[60,] "rs2274517"  "1"       "0.999936465724214"
#[61,] "rs17715343" "1"       "0.871397958891365"
#[62,] "rs12423664" "1"       "0.999939387280693"
#[63,] "rs13205911" "1"       "0.981016328167447"
#[64,] "rs2180314"  "1"       "0.99730646658639" 
#[65,] "rs1473886"  "1"       "0.999972842887614"
#[66,] "rs10501321" "1"       "0.983504674485564"
#[67,] "rs4075205"  "1"       "0.997391138591597"
#[68,] "rs11229606" "1"       "0.749500163117371"
#[69,] "rs6066141"  "1"       "0.997091041310893"
#[70,] "rs2976940"  "1"       "0.998001128169822"
#[71,] "rs1045241"  "1"       "0.9999031067678"  
#[72,] "rs9425592"  "1"       "0.999916356750539"
#[73,] "rs2215169"  "1"       "0.995456756366827"
#[74,] "rs2862954"  "1"       "0.733283872231295"


#Old results
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


pdf("plots.bychr.vs2.pdf")
for(i in 1:22){
 if(sum(sub$chr==i)>0){
plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>50000))
abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>3.989499],col=3)
}
}
dev.off()

#try to plot all hits to see relationships?
tophits =  sub[sub$lbfav>3.989499 & sub$nmin>50000,]
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

write.table(file="GlobalLipids2013.newtophits.vs3.SignCrrct.vs1.txt",cbind(sub[sub$lbfav>3.989499 & sub$nmin>50000 & sub$annot==0,c(4,2:3,7)],round(sub[sub$lbfav>3.989499 & sub$nmin>50000 & sub$annot==0,c(22:25,30)],digits=1)),quote=FALSE,sep= " ", row.names=FALSE)

c(4,2:3,7,22:25,30)

#newhits = gl[gl$annot==0 & gl$lbfav>3.989499 & gl$nmin>50000,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)

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





#NOTE -- everything below is 100% code I have added
~~~~~~~~~~~~~~

#lbf.uni.newhits = lbf.uni(lbf.newhits,lbf$gamma)
#lbf.all.newhits = lbf.all(lbf.newhits,lbf$gamma)
#
#lbf.uni.newhits < lbf.all.newhits
#
#newhits.BFall <- newhits[lbf.uni.newhits < lbf.all.newhits,]
#newhits.BFuni <- newhits[lbf.uni.newhits > lbf.all.newhits,]
#
#pdf("plots.newhits.BFcompare.vs2.pdf")
#
#plot(lbf.uni.newhits, lbf.all.newhits, xlab= "Univariate BFs", ylab="Multivariate BFs", main="Newhits: Univariate vs. Multivariate BFs")
#abline(a=0, b=1)
#abline(lm(lbf.all.newhits~lbf.uni.newhits), col="RED")
#
#dev.off()
#
#
#lbf.all.newhits[lbf.uni.newhits < lbf.all.newhits]
# [1]   9.805931   8.815734   7.231060   4.894935   5.033745  10.401903
#  [7]  10.983507   4.659718  20.541163   4.102946  29.222851  62.303929
#  [13]   3.900526   5.221707  19.880882  11.820379  11.132608  19.212606
#  [19]   8.972023 140.081233   5.558127   8.683993   4.307358  10.276061
#  [25]   8.120480  16.679203   5.092154  11.343710   4.585494   5.977829
#  [31]   3.978291  21.572096  48.048088   7.436775   9.149239   4.054293
#  [37]   8.025160  20.350255   8.864311  14.709923  86.440319  36.891124
#  [43]   6.539004   4.139504  29.890960  21.913561   3.892162   4.339384
#  [49] 139.058244   5.179194   7.629123   4.398708 135.404293  26.783135
#  [55]   8.002221  15.617592  12.197956  26.419773
#
#
#
#lbf.all.newhits[lbf.uni.newhits > lbf.all.newhits]
# [1]  6.225993  2.753010 11.857290 14.678976  7.044515  4.527221  3.265205
#  [8]  5.899027 11.843284  4.207585  4.994016  6.309116  2.705491  3.631810
#  [15]  3.699085  3.102687
#
#
#
#newhits[lbf.all.newhits > 113,]
#snp chr       pos     maf     Z.tg       Z.tc      Z.ldl     Z.hdl
#6866 rs4897361   6 130364328 0.43400 4.733115 3.36820251  2.7053620 2.8532835
#2753  rs721772  15  41829230 0.48150 3.629746 0.01040269  0.8219594 5.1686429
#6219 rs7255743  19  46018119 0.02639 2.624603 7.13597998 10.4537794 0.6782708
#lbfav
#6866  4.757524
#2753  5.355168
#6219 20.611843
#
#
#
#newhits$lbfall = lbf.all.newhits
#newhits$lbfuni = lbf.uni.newhits
#
#write.table(file="newhits.vs2.txt",newhits,quote=FALSE,row.names=FALSE)
#
#gl.glhits$lbfav = lbf.av.glhits
#gl.glhits$lbfall = lbf.all.glhits
#gl.glhits$lbfuni = lbf.uni.glhits
#
#write.table(file="glhits.vs2.txt",gl.glhits[,c(4,2,7,22,23,24,25,30,31,32)],quote=FALSE,row.names=FALSE)
#
#
#NOTE -- all the below was to help troubleshoot the odd lbf* results where it seemed like the wrong lbf* values were coming up for the top new variants. This was because of the o = order(gl$chr, gl$pos).... command and rearranging the gl matrix /after/ getting the lbf* so that the wrong rows were being called when using the previous order from lbf* on the new ordered gl matrix. Keeping all the old results for records. 
#gl2 = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz", header=T)
#gl2[row.names(newhits),]$annot <- 3
#gl2 <- gl2[!is.na(gl2$maf) & gl2$maf > 0,]
#gl2$lbfav = lbf.av.all
#
#lbf.newhits2 = lbf.bigmat[,gl2$annot==3]
#
#lbf.uni.newhits2 = lbf.uni(lbf.newhits2,lbf$gamma)
#lbf.all.newhits2 = lbf.all(lbf.newhits2,lbf$gamma)
#
#lbf.all.newhits2[lbf.uni.newhits2 < lbf.all.newhits2]
# [1]  5.113600  9.238360  6.887086  4.339257  4.501904  4.378643  5.930985
#  [8]  4.371434  5.184725  4.125246  5.819409 17.316921  4.514223  4.770906
#  [15]  5.259777  5.139357  6.664481  5.397071  4.275684  4.560132  4.154638
#  [22]  4.857058  5.761283  5.226976  6.581458  4.150895  5.575327  4.252227
#  [29]  4.292621  4.981670  6.722394 25.754510  5.658189 11.867786  4.854878
#  [36]  7.653760  6.827153  5.810097  5.007101  5.123809  5.980351  4.851327
#  [43] 31.807344  4.781310  5.977829  4.754488  6.047135  5.258339  9.998777
#  [50]  6.890773  5.988893 11.425993  4.190938  6.475299  4.686178  4.280777
#  [57]  6.004999  5.289372  4.777772 11.061500 36.511011  8.021977  5.313751
#  [64]  4.773029  5.407697  5.001026  5.908662
#
#
#lbf.all.newhits2[lbf.uni.newhits2 > lbf.all.newhits2]
#[1]  4.767202  8.055220 20.992954  5.048607  5.396590  6.225892  4.796701
#
#
#newhits2 = gl2[gl2$annot==3,]
#
#newhits2$lbfall = lbf.all.newhits2
#newhits2$lbfuni = lbf.uni.newhits2
#
#write.table(file="newhits2.vs1.txt",newhits2[,c(4,2,7,22,23,24,25,29,30,31)],quote=FALSE,row.names=FALSE)
#
#
#
#lbf.glhits = lbf.bigmat[,gl$annot==1]
#
#
##newhits = sub[sub$annot==0 & sub$lbfav>3.989499 & sub$nmin>50000,c(4,2:3,7,22:25,30)]
#
#
#lbf.uni.allhits = lbf.uni(lbf.bigmat, lbf$gamma)
#lbf.all.allhits = lbf.all(lbf.bigmat, lbf$gamma)
#
#lbf.all.allhits == 139.058244 
#


lbf.uni.newhits = lbf.uni(lbf.newhits,lbf$gamma)
lbf.all.newhits = lbf.all(lbf.newhits,lbf$gamma)

lbf.uni.newhits < lbf.all.newhits

newhits.BFall <- newhits[lbf.uni.newhits < lbf.all.newhits,]
newhits.BFuni <- newhits[lbf.uni.newhits > lbf.all.newhits,]

pdf("plots.newhits.BFcompare.vs2.pdf")

plot(lbf.uni.newhits, lbf.all.newhits, xlab= "Univariate BFs", ylab="Multivariate BFs", main="Newhits: Univariate vs. Multivariate BFs")
abline(a=0, b=1)
abline(lm(lbf.all.newhits~lbf.uni.newhits), col="RED")

dev.off()


lbf.all.newhits[lbf.uni.newhits < lbf.all.newhits]
 [1]  5.113600  9.238360  6.887086  4.339257  4.501904  4.378643  5.930985
  [8]  4.371434  5.184725  4.125246  5.819409 17.316921  4.514223  4.770906
  [15]  5.259777  5.139357  6.664481  5.397071  4.275684  4.560132  4.154638
  [22]  4.857058  5.761283  5.226976  6.581458  4.150895  5.575327  4.252227
  [29]  4.292621  4.981670  6.722394 25.754510  5.658189 11.867786  4.854878
  [36]  7.653760  6.827153  5.810097  5.007101  5.123809  5.980351  4.851327
  [43] 31.807344  4.781310  5.977829  4.754488  6.047135  5.258339  9.998777
  [50]  6.890773  5.988893 11.425993  4.190938  6.475299  4.686178  4.280777
  [57]  6.004999  5.289372  4.777772 11.061500 36.511011  8.021977  5.313751
  [64]  4.773029  5.407697  5.001026  5.908662

lbf.all.newhits[lbf.uni.newhits > lbf.all.newhits]
[1]  4.767202  8.055220 20.992954  5.048607  5.396590  6.225892  4.796701

newhits$lbfall = lbf.all.newhits
newhits$lbfuni = lbf.uni.newhits

write.table(file="newhits.vs2.txt",newhits,quote=FALSE,row.names=FALSE)

#newhits[lbf.all.newhits > 31.8 & lbf.all.newhits < 31.9,]
newhits[lbf.all.newhits > 8,]

snp chr       pos     maf      Z.tg      Z.tc      Z.ldl
130    rs7033354   9  16904846 0.35360 5.0490124 3.6199148  4.8218961
1433  rs11040329  11  49339071 0.14250 1.3626273 2.3391050  0.3564548
2291  rs17630235  12 112591686 0.42220 1.2233435 7.9045052  6.5536699
5270   rs1872167  11  47901269 0.14780 2.2769368 4.4005035  1.1452636
5419   rs4988235   2 136608646 0.47630 0.6931746 7.5618096  6.6362335
6219   rs7255743  19  46018119 0.02639 2.6246032 7.1359800 10.4537794
6948   rs2247056   6  31265490 0.21770 9.4363205 9.2520149  5.6710801
7876  rs12133576   1  93816400 0.35490 2.5142081 4.4282290  2.8918342
8442   rs2454722  12 123171218 0.14510 1.7660349 3.6741955  1.8017802
11051  rs1473886   2  20368519 0.46170 4.2783309 6.2173255  3.5235542
11056 rs10501321  11  47294626 0.31400 5.6719272 1.5121331  2.0296709
11460  rs4075205   8 144284709 0.45120 1.9467059 0.1502087  3.3455147
Z.hdl     lbfav    lbfall    lbfuni
130    3.9923200  9.561825  9.238360  3.852086
1433   6.8930726  8.966206  8.055220  8.398072
2291   5.8145698 18.025654 17.316921 11.572406
5270  10.8221024 26.424143 25.754510 23.323655
5419   3.2509996 12.435684 11.867786 10.439136
6219   0.6782708 20.611843 20.992954 21.417230
6948   2.8949657 32.067766 31.807344 17.353809
7876   6.5400458 10.583552  9.998777  7.368020
8442   7.5856595 12.331865 11.425993 10.532322
11051  4.3893308 11.244442 11.061500  6.492045
11056 12.9184728 37.346393 36.511011 34.062316
11460  5.9042754  8.057769  8.021977  5.688962


#rs10501321 falls in MADD but is just downstream of NR1H3 which comes up as an 'other enriched gene' in supplement of 2013 paper, specifically pathway analysis where 'other enriched gene' means it's a gene that does not reach genome-wide significance but is enriched in the analysis
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs10501321
'chr11:47251202 chr11:47294626  rs10501321      c       t       0.0483  0.0036  186984.00       3.541e-38       0.314
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs10501321
chr11:47251202  chr11:47294626  rs10501321      t       c       0.0087  0.0039  172916.00       0.04239 0.686
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs10501321
chr11:47251202  chr11:47294626  rs10501321      t       c       0.0216  0.0035  177680.00       1.412e-08       0.686
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs10501321
chr11:47251202  chr11:47294626  rs10501321      c       t       0.006   0.0038  187172.00       0.1305  0.314

#rs2247056 falls in HLA-B and there's already one SNP representing the HLA-DRA in the 2013 paper (rs3177928) -- so maybe just masked any other HLA SNP? HLA-B and HLA-DRA are ~1Mb apart. Also one of the few SNPs where lbfuni > lbfall. The other SNPs that fall in this category have values between 4-8
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs2247056
chr6:31373469   chr6:31265490   rs2247056       c       t       0.0121  0.0040  183368.00       0.003792        0.7823
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs2247056
chr6:31373469   chr6:31265490   rs2247056       c       t       0.0248  0.0043  169288.00       1.419e-08       0.7823
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs2247056
chr6:31373469   chr6:31265490   rs2247056       c       t       0.0378  0.0039  174062.00       3.861e-21       0.7823
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs2247056
chr6:31373469   chr6:31265490   rs2247056       c       t       0.0391  0.0041  183557.00       2.203e-20       0.7823
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs3177928
chr6:32520413   chr6:32412435   rs3177928       a       g       0.0148  0.0048  179804.90       0.002925        0.1807
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs3177928
chr6:32520413   chr6:32412435   rs3177928       a       g       0.0452  0.0052  165751.00       3.096e-17       0.1807
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs3177928
chr6:32520413   chr6:32412435   rs3177928       a       g       0.0135  0.0048  170516.00       0.003179        0.1807
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs3177928
chr6:32520413   chr6:32412435   rs3177928       a       g       0.0482  0.0050  179995.90       9.778e-22       0.1807

#rs1872167 falls downstream of NUP160 and is in a region enriched for H3K4Me1 in Gm12878 (Female LCL) based on ENCODE track. NUP160 is not in 2013 supplement. Does not show up in GTEx portal. Maybe other nearby SNPs are associated somewhere?
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs1872167
chr11:47857845  chr11:47901269  rs1872167       t       c       0.0548  0.0048  178710.01       2.705e-27       0.1478
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs1872167
chr11:47857845  chr11:47901269  rs1872167       t       c       0.0049  0.0051  164709.01       0.2521  0.1478
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs1872167
chr11:47857845  chr11:47901269  rs1872167       c       t       0.0124  0.0046  169360.02       0.02279 0.8522
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs1872167
chr11:47857845  chr11:47901269  rs1872167       t       c       0.0222  0.0049  178854.32       1.08e-05        0.1478

#rs7255743 falls in VASP which is also an 'other enriched gene' in the supplement of the 2013 paper. Falls in region with a ton of H3K4Me1 enrichment in multiple tissues
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs7255743
chr19:50709959  chr19:46018119  rs7255743       a       g       0.0027  0.0146  137200.49       0.4976  0.02639
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs7255743
chr19:50709959  chr19:46018119  rs7255743       g       a       0.1647  0.0156  123958.77       1.408e-25       0.97361
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs7255743
chr19:50709959  chr19:46018119  rs7255743       a       g       0.0323  0.0138  127634.90       0.008675        0.02639
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs7255743
chr19:50709959  chr19:46018119  rs7255743       g       a       0.1066  0.0150  135895.27       9.61e-13        0.97361

#rs17630235 downstream of TRAFD1, which is not in supplement of 2013 paper. Nothing in GTEx or dbSNP
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs17630235
chr12:111076069 chr12:112591686 rs17630235      g       a       0.0213  0.0035  187081.90       6.079e-09       0.5778
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs17630235
chr12:111076069 chr12:112591686 rs17630235      g       a       0.026   0.0038  173014.10       5.614e-11       0.5778
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs17630235
chr12:111076069 chr12:112591686 rs17630235      a       g       0.0067  0.0034  177781.90       0.2212  0.4222
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs17630235
chr12:111076069 chr12:112591686 rs17630235      g       a       0.0298  0.0036  187267.00       2.69e-15        0.5778

#rs4988235 falls in MCM6 which is not in supplement of 2013 paper (downstream of LCT which is also not in 2013 paper). Nothing in GTEx, gets pulled in with LCT for 'Lactase Persistence' as a clinical phenotype in dbSNP
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs4988235
chr2:136325116  chr2:136608646  rs4988235       g       a       0.0136  0.0039  183569.61       0.00115 0.4763
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs4988235
chr2:136325116  chr2:136608646  rs4988235       g       a       0.0278  0.0042  169531.02       3.218e-11       0.4763
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs4988235
chr2:136325116  chr2:136608646  rs4988235       a       g       0.0031  0.0038  174266.97       0.4882  0.5237
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs4988235
chr2:136325116  chr2:136608646  rs4988235       g       a       0.0308  0.0040  183760.79       3.975e-14       0.4763

#rs2454722 falls in GPR81 which is not in supplement of 2013 paper (upstream of NIACR1 and GPR146 is in supplement though on a different chromosome). Nothing in GTEx or dbSNP
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs2454722
chr12:121737171 chr12:123171218 rs2454722       g       a       0.0351  0.0044  186355.01       3.308e-14       0.1451
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs2454722
chr12:121737171 chr12:123171218 rs2454722       g       a       0.0089  0.0047  172323.04       0.07158 0.1451
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs2454722
chr12:121737171 chr12:123171218 rs2454722       a       g       0.0067  0.0043  177056.04       0.07739 0.8549
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs2454722
chr12:121737171 chr12:123171218 rs2454722       g       a       0.0172  0.0045  186538.48       0.0002386       0.1451

#rs1473886 upstream of SDC1 which is also an 'other enriched gene' in the suppleent of 2013 paper. In a region of ENCODE mark enrichment though right at the SNP nothing particularly large. Nothing in GTEx or dbSNP
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs1473886
chr2:20232000   chr2:20368519   rs1473886       g       t       0.0163  0.0034  186880.00       1.137e-05       0.5383
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs1473886
chr2:20232000   chr2:20368519   rs1473886       g       t       0.0135  0.0036  172815.00       0.0004258       0.5383
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs1473886
chr2:20232000   chr2:20368519   rs1473886       g       t       0.0157  0.0033  177577.10       1.883e-05       0.5383
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs1473886
chr2:20232000   chr2:20368519   rs1473886       g       t       0.0236  0.0035  187070.00       5.057e-10       0.5383

#rs12133576 falls in DR1 which is not in supplement of 2013 paper (downstream of CR609342 and AL832786). GTEx shows a significant association with adipose tissue and DR1 http://www.gtexportal.org/home/eqtls/bySnp?snpId=rs12133576&tissueName=All and nothing in dbSNP
#ENSG00000117505.7	DR1	rs12133576	1.7E-7	Adipose_Subcutaneous	Show eQTL box plot
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs12133576
chr1:93588988   chr1:93816400   rs12133576      a       g       0.0243  0.0035  187123.10       6.15e-11        0.3549
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs12133576
chr1:93588988   chr1:93816400   rs12133576      a       g       0.0104  0.0038  173042.00       0.00383 0.3549
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs12133576
chr1:93588988   chr1:93816400   rs12133576      g       a       0.009   0.0034  177818.00       0.01193 0.6451
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs12133576
chr1:93588988   chr1:93816400   rs12133576      a       g       0.0162  0.0037  187314.00       9.501e-06       0.3549

#rs7033354 falls downstream of BNC2 which is not in supplement of 2013 paper. Nothing in GTEx or dbSNP
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs7033354
chr9:16894846   chr9:16904846   rs7033354       t       c       0.0154  0.0035  187047.00       6.543e-05       0.6464
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs7033354
chr9:16894846   chr9:16904846   rs7033354       c       t       0.0189  0.0038  172976.00       1.422e-06       0.3536
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs7033354
chr9:16894846   chr9:16904846   rs7033354       c       t       0.019   0.0034  177742.00       4.441e-07       0.3536
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs7033354
chr9:16894846   chr9:16904846   rs7033354       c       t       0.0139  0.0036  187232.90       0.0002947       0.3536

#rs11040329 falls downstream of FOLH1, FGCP and PSMAL none of which are in supplement of 2013 paper. Nothing in GTEx or dbSNP
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs11040329
chr11:49295647  chr11:49339071  rs11040329      c       t       0.038   0.0054  147407.93       5.46e-12        0.1425
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs11040329
chr11:49295647  chr11:49339071  rs11040329      c       t       0.0001  0.0058  139598.03       0.7215  0.1425
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs11040329
chr11:49295647  chr11:49339071  rs11040329      t       c       0.0102  0.0051  143686.93       0.173   0.8575
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs11040329
chr11:49295647  chr11:49339071  rs11040329      c       t       0.012   0.0056  151532.13       0.01933 0.1425

#rs4075205 falls upstream of GPIHBP1 which is not in supplement of 2013 paper. Nothing in GTEx or dbSNP 
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs4075205
chr8:144356084  chr8:144284709  rs4075205       t       c       0.0224  0.0035  186000.00       3.542e-09       0.5488
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs4075205
chr8:144356084  chr8:144284709  rs4075205       c       t       0.0119  0.0038  171977.00       0.0008213       0.4512
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs4075205
chr8:144356084  chr8:144284709  rs4075205       c       t       0.009   0.0034  176721.00       0.05157 0.4512
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs4075205
chr8:144356084  chr8:144284709  rs4075205       t       c       0.0003  0.0036  186169.00       0.8806  0.5488



newhits.BFavOverUni <- newhits[newhits$lbfav > newhits$lbfuni,]
newhits.BFuniOverAv <- newhits[newhits$lbfav < newhits$lbfuni,]
newhits.BFavOverAll <- newhits[newhits$lbfav > newhits$lbfall,]
newhits.BFavOverBoth <- newhits[newhits$lbfav > newhits$lbfuni & newhits$lbfav > newhits$lbfall,]

> newhits.BFuniOverAv
snp chr      pos     maf     Z.tg      Z.tc     Z.ldl     Z.hdl
6219  rs7255743  19 46018119 0.02639 2.624603 7.1359800 10.453779 0.6782708
7424  rs1062219   8 11616410 0.49210 6.025488 0.4162870  1.773173 1.3312706
12171 rs2976940   8  8279462 0.44720 5.708763 0.4038813  1.977546 0.7941639
lbfav    lbfall    lbfuni
6219  20.611843 20.992954 21.417230
7424   5.727299  5.396590  5.994709
12171  5.021955  4.796701  5.209468

#rs1062219 is in GATA4 which is not in the supplement of 2013 paper. Nothing in GTEx or dbSNP 

#rs2976940 is downstream of PRAGMIN and SGK223, neither of which are in the supplment of 2013 paper. eQTL with CTA-398F10.2 with Skin_Sun_Exposed_Lower_leg and nothing in dbSNP
#ENSG00000254153.1	CTA-398F10.2	rs2976940	5.0E-7	Skin_Sun_Exposed_Lower_leg	Show eQTL box plot







#20160630 -- New additions for NHS/committee meeting misc

#newhits = sub[sub$annot==0 & sub$lbfav>3.989499 & sub$nmin>50000,c(4,2:3,7,22:25,30)]
failhits = sub[sub$annot==0 & sub$lbfav<=3.989499 & sub$nmin>50000,c(4,2:3,7,22:25,30)]
prevhits = gl.glhits[,c(2:4,7,22:25)] 


jpeg("GlobalLipids2013.vs3.SignCrrct.vs1.ForCSHLPoster.NHSSlides.ZScorePatterns.1.jpg", width=4000, height=2000, res=300)

par(mfrow=c(1,2), mar=c(4,5,5,4), oma=c(2.25,2.25,2.25,2.25))

#plot(c(apply(prevhits[,5:8], 1, median), apply(newhits[,5:8], 1, median), apply(failhits[,5:8], 1, median)), c(apply(prevhits[,5:8], 1, max), apply(newhits[,5:8], 1, max), apply(failhits[,5:8], 1, max)), main="Median vs Max Z-Scores of Old and New Hits", xlab="Median of Z-scores", ylab="Max of Z-scores^2", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])))
##plot(c(apply(prevhits[,5:8]^2, 1, median), apply(newhits[,5:8]^2, 1, median), apply(failhits[,5:8]^2, 1, median)), c(apply(prevhits[,5:8]^2, 1, max), apply(newhits[,5:8]^2, 1, max), apply(failhits[,5:8]^2, 1, max)), main="Median vs Max Z-Scores of Old and New Hits", xlab="Median of Z-scores", ylab="Max of Z-scores^2", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])))
plot(c(apply(prevhits[,5:8]^2, 1, median), apply(newhits[,5:8]^2, 1, median)), c(apply(prevhits[,5:8]^2, 1, sd), apply(newhits[,5:8]^2, 1, sd)), main="Median vs SD Z-Scores of Old and New Hits", xlab="Median of Zscores", ylab="SD of Zscores", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])))

abline(a=0, b=1)
#legend(.075,38.25, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))
##legend(.075,1450, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))
legend(1,17, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))

#legend(.075,100, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))

plot(c(apply(prevhits[,5:8]^2, 1, max), apply(newhits[,5:8]^2, 1, max)), c(apply(prevhits[,5:8]^2, 1, sum), apply(newhits[,5:8]^2, 1, sum)), main="Max Z-Score^2 vs Sum(Z-Scores^2) of Old and New Hits", xlab="Max Z-Score^2", ylab="Sum(Z-Scores^2)", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])), xlim=c(0,500), ylim=c(0,500))

abline(a=0, b=1)
legend(1,500, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))

dev.off()







#This one goes across the 14 sigmaa rows per gamma to take the mean, producing a single per gamma value across each SNP -- final matrix of 81 x 84 (gammas x SNPs)
SumAcrossSigmaas <- function(lbf, ngamma, nsigmaa) {
	lbfSummedOverSigmas <- matrix(0, ncol=ncol(lbf), nrow=ngamma)
	for (i in 1:ngamma) {
		coords <- seq.int(from=i, by=ngamma, length.out=nsigmaa)
		max <- apply(lbf[coords,], 2, max)
		lbf[coords,] <- lbf[coords,] - matrix(max, nrow=nrow(lbf[coords,]), ncol=ncol(lbf[coords,]), byrow=TRUE)
		lbfSummedOverSigmas[i,] <- log10(apply(10^lbf[coords,], 2, mean)) + max
	}
	return(lbfSummedOverSigmas)
}

#> seq.int(from = 1, by=27, length.out=14)
# [1]    1   82  163  244  325  406  487  568  649  730  811  892  973 1054

lbf.newhits.sigmaaSummed <- SumAcrossSigmaas(lbf.newhits, 81, 14)
lbf.newhits.sigmaaSummed.plusGammaDescrips <- cbind(lbf$gamma, lbf.newhits.sigmaaSummed)


CountGammaValues <- function(x) { count0 <- 0; count1 <- 0; count2 <- 0; for (j in 1:length(x)) { if (x[j] == 0) { count0 <- count0 + 1; }; if (x[j] == 1) { count1 <- count1 + 1; }; if (x[j] == 2) { count2 <- count2 + 1;}; }; return(c(count0, count1, count2));}

UniBFs <- lbf.newhits.sigmaaSummed.plusGammaDescrips[t(apply(lbf.newhits.sigmaaSummed.plusGammaDescrips[,1:4], 1, CountGammaValues))[,2] == 1 & t(apply(lbf.newhits.sigmaaSummed.plusGammaDescrips[,1:4], 1, CountGammaValues))[,1] == 0,]
UniBFs.max <- apply(UniBFs[,5:ncol(UniBFs)], 2, max)
MultiBFs <- lbf.newhits.sigmaaSummed.plusGammaDescrips[t(apply(lbf.newhits.sigmaaSummed.plusGammaDescrips[,1:4], 1, CountGammaValues))[,2] > 1 | ( t(apply(lbf.newhits.sigmaaSummed.plusGammaDescrips[,1:4], 1, CountGammaValues))[,2] == 1 & t(apply(lbf.newhits.sigmaaSummed.plusGammaDescrips[,1:4], 1, CountGammaValues))[,1] > 0),] 
MultiBFs.max <- apply(MultiBFs[,5:ncol(MultiBFs)], 2, max)
AllBFs <- lbf.newhits.sigmaaSummed.plusGammaDescrips[t(apply(lbf.newhits.sigmaaSummed.plusGammaDescrips[,1:4], 1, CountGammaValues))[,2] == 4,]


jpeg("GlobalLipids2013.vs3.SignCrrct.vs1.20160616CommitteeMeetingSlides.pt2.jpg", width=4500, height=4500, res=300)
par(mfrow=c(2,2), mar=c(6,5,4,2))

plot(UniBFs.max, AllBFs[5:length(AllBFs)], main="Global Lipids 2013 -- Max Uni BF vs. All Dir. Assoc. BF per SNP", xlab="Max Univariate BF", ylab="All Directly Associated BF", cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
abline(0,1, col="BLACK")
plot(UniBFs.max, MultiBFs.max, main="Global Lipids 2013 -- Max Uni BF vs. Max Multi BF per SNP", xlab="Max Univariate BF", ylab="Max Multivariate BF", cex=1.5, cex.main=1.5, cex.axis=1.5, cex.lab=1.5)
abline(0,1, col="BLACK")

dev.off()

length(UniBFs.max[UniBFs.max<AllBFs[5:length(AllBFs)]])
length(UniBFs.max[UniBFs.max>AllBFs[5:length(AllBFs)]])
mean(AllBFs[5:length(AllBFs)]-UniBFs.max)
sd(AllBFs[5:length(AllBFs)]-UniBFs.max)
quantile(AllBFs[5:length(AllBFs)]-UniBFs.max)
length(UniBFs.max[UniBFs.max<MultiBFs.max])
length(UniBFs.max[UniBFs.max>MultiBFs.max])
mean(MultiBFs.max-UniBFs.max)
sd(MultiBFs.max-UniBFs.max)
quantile(MultiBFs.max-UniBFs.max)

~~~
> length(UniBFs.max[UniBFs.max<AllBFs[5:length(AllBFs)]])
[1] 84
> length(UniBFs.max[UniBFs.max>AllBFs[5:length(AllBFs)]])
[1] 0
> mean(AllBFs[5:length(AllBFs)]-UniBFs.max)
[1] 30.30219
> sd(AllBFs[5:length(AllBFs)]-UniBFs.max)
[1] 27.89008
> quantile(AllBFs[5:length(AllBFs)]-UniBFs.max)
       0%       25%       50%       75%      100% 
  3.40921  12.21924  20.42154  40.46465 170.52648 
> length(UniBFs.max[UniBFs.max<MultiBFs.max])
[1] 84
> length(UniBFs.max[UniBFs.max>MultiBFs.max])
[1] 0
> mean(MultiBFs.max-UniBFs.max)
[1] 31.36831
> sd(MultiBFs.max-UniBFs.max)
[1] 27.82734
> quantile(MultiBFs.max-UniBFs.max)
        0%        25%        50%        75%       100% 
  3.840307  12.590121  21.477026  41.239183 171.488950 
~~~







#20160830
#Heatmap on per phenotype marginal LBFs, and other pre-ProbGenom ideas
#Likely work to help setup moving onto next phase pre/during September 2016

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

#DirPhenoGroupsPerPheno <- c()
##lbf$gamma
##?
#for (i in 1:ncol(lbf$gamma)) {
#	DirPhenoGroupsPerPheno <- cbind(DirPhenoGroupsPerPheno, lbf$gamma[,i] == 1)
#}

GetDirAssocModels <- function (x) {
	returnValue1 <- FALSE
	if (x == 1) {
		returnValue1 <- TRUE
	}
	return(returnValue1)
}

GetModelsOfSpecificAssocType <- function (x, assocType) {
	returnValue1 <- FALSE
	if (x == assocType) {
		returnValue1 <- TRUE
	}
	return(returnValue1)
}

#DirAssocModels <- apply(lbf$gamma, c(1,2), GetDirAssocModels)
Gamma_DirAssocModels <- apply(lbf$gamma, c(1,2), GetModelsOfSpecificAssocType, assocType=1)
do.call(rbind, replicate(length(sigmaa), Gamma_DirAssocModels, simplify=FALSE))

lbf.bigmat.phenoMarginals <- c()
for (i in seq.int(from=1, by=nrow(lbf$gamma), length.out=length(sigmaa))) {
	lbf.bigmat.phenoMarginals.addOn <- c()
	for (j in 1:ncol(lbf$gamma)) {
		tempVals <- lbf.bigmat[i:(i+nrow(lbf$gamma)-1),][Gamma_DirAssocModels[,j],]
		max <- apply(tempVals, 2, max)
		tempVals <- tempVals - matrix(max, nrow=nrow(tempVals), ncol=ncol(tempVals), byrow=TRUE)
		tempVals <- log10(apply(10^tempVals, 2, sum))
		tempVals <- tempVals + (max * 27)

	
		tempVals <- log10(apply(10^lbf.bigmat[i:(i+nrow(lbf$gamma)-1),][Gamma_DirAssocModels[,j],], 2, sum))
		#Below is for summing multiple 10^0 values produces 1 * whatever, which after log10 does not return 0
		tempVals <- sapply(tempVals, function(x) { val1 <- x; if (val1 == log10(sum(Gamma_DirAssocModels[,j]))) { val1 <- 0; }; return(val1); })
		lbf.bigmat.phenoMarginals.addOn <- rbind(lbf.bigmat.phenoMarginals.addOn, tempVals)
#		lbf.bigmat.phenoMarginals.addOn <- rbind(lbf.bigmat.phenoMarginals.addOn, log10(apply(10^lbf.bigmat[i:(i+nrow(lbf$gamma)-1),][Gamma_DirAssocModels[,j],], 2, sum)))
	}
	lbf.bigmat.phenoMarginals <- rbind(lbf.bigmat.phenoMarginals, lbf.bigmat.phenoMarginals.addOn)
}

lbf.bigmat.phenoMarginals.MeanAcrossSigmaas <- MeanAcrossSigmaas(lbf.bigmat.phenoMarginals, 4, 14)
lbf.bigmat.phenoMarginals.MeanAcrossSigmaas.Normalized <- apply(lbf.bigmat.phenoMarginals.MeanAcrossSigmaas, 2, normalize)

colnames(lbf.bigmat.phenoMarginals.MeanAcrossSigmaas.Normalized) <- gl$snp
row.names(lbf.bigmat.phenoMarginals.MeanAcrossSigmaas.Normalized) <- c("TC", "TG", "HDL", "LDL")





#lbf.bigmat.pp <- posteriorprob(
#sub = gl[gl$nmin>50000,]
#l=indephits(sub$lbfav,sub$chr,sub$pos)
#sub=sub[l==1,]

#20160831 NOTE -- Kept thinking ggplot2 wasn't getting the order right with the factor & level settings, but found out ggplot2 orders thing going from the bottom left up; so things were correct I just needed to reverse all my factor orders in order to get the visual I was looking for
#20160901 NOTE -- Found out about scale_y_reverse(), so going to keep the 'decreasing=TRUE' statements in the order calls, but then use graph + scale_y_reverse() to fix it visually since I normally think/associate the actual data matrices as being in decreasing=TRUE sort orders
#20160901 NOTE -- Then I found out this scale_y_reverse() function doesn't really work when the axes are strings/characters vs. numerical values, so back to original approach of using no 'decreasing=TRUE' calls
#20160901 NOTE -- Also found out that factor level order can be different from the order of the actual matrix, and ggplot2/heatmap/melt will follow order of levels regardless. This happened a bit by accident but as it turns out can keep 'decreasing=TRUE' for matrix but remove it for factor levels, though that could also possibly be less than preferred (probably not expected/intuitive if use matrix for anything after that other than the plot...)

sub.pp <- posteriorprob(lbf.bigmat[,l==1],ebprior.glhits)
sub.pp.marginals <- marginal.postprobs(sub.pp, lbf$gamma, length(sigmaa))
#sub.pp.marginals.DirAssoc <- sub.pp.marginals[[2]]
sub.pp.marginals.DirAssoc <- t(sub.pp.marginals[[2]])
#colnames(sub.pp.marginals.DirAssoc) <- sub$snp
#sub.pp.marginals.DirAssoc <- rbind(as.character(sub$snp), sub.pp.marginals.DirAssoc)
sub.pp.marginals.DirAssoc <- data.frame(cbind(as.character(sub$snp), sub.pp.marginals.DirAssoc))
sub.pp.marginals.DirAssoc[,2] <- as.numeric(as.character(sub.pp.marginals.DirAssoc[,2]))
sub.pp.marginals.DirAssoc[,3] <- as.numeric(as.character(sub.pp.marginals.DirAssoc[,3]))
sub.pp.marginals.DirAssoc[,4] <- as.numeric(as.character(sub.pp.marginals.DirAssoc[,4]))
sub.pp.marginals.DirAssoc[,5] <- as.numeric(as.character(sub.pp.marginals.DirAssoc[,5]))
#rownames(sub.pp.marginals.DirAssoc) <- c("TC", "TG", "HDL", "LDL")
#rownames(sub.pp.marginals.DirAssoc) <- c("LDL", "HDL", "TC", "TG")
#rownames(sub.pp.marginals.DirAssoc) <- c("snp", "LDL", "HDL", "TC", "TG")
colnames(sub.pp.marginals.DirAssoc) <- c("snp", "LDL", "HDL", "TC", "TG")
sub.pp.marginals.DirAssoc <- sub.pp.marginals.DirAssoc[order(sub$lbfav, decreasing=TRUE),]
sub.pp.marginals.DirAssoc$snp <- factor(sub.pp.marginals.DirAssoc$snp, levels=sub[order(sub$lbfav),][,4])

#newhits$snp <- factor(newhits$snp, levels=newhits[order(newhits$lbfav),][,1]) 

sub.pp.marginals.DirAssoc.newhits <- sub.pp.marginals.DirAssoc[sub[order(sub$lbfav, decreasing=TRUE),]$annot==0 & sub[order(sub$lbfav, decreasing=TRUE),]$lbfav>3.989499 & sub[order(sub$lbfav, decreasing=TRUE),]$nmin>50000,]
sub.pp.marginals.DirAssoc.newhits$snp <- factor(sub.pp.marginals.DirAssoc.newhits$snp, levels=sub[sub$annot==0 & sub$lbfav>3.989499 & sub$nmin>50000,][order(sub[sub$annot==0 & sub$lbfav>3.989499 & sub$nmin>50000,]$lbfav),][,4])

sub.pp.marginals.DirAssoc.nonGWASregions <- sub.pp.marginals.DirAssoc[sub[order(sub$lbfav, decreasing=TRUE),]$annot==0,]
sub.pp.marginals.DirAssoc.nonGWASregions$snp <- factor(sub.pp.marginals.DirAssoc.nonGWASregions$snp, levels=sub[sub$annot==0,][order(sub[sub$annot==0,]$lbfav),][,4])

#gl.glhits = gl[gl$annot==1,] # subset of SNPs reported in Teslovich 2010 paper
#marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))
#lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
#gl$lbfav = lbf.av.all
#
#lbf.glhits = lbf.bigmat[,gl$annot==1]
#pp.glhits = posteriorprob(lbf.glhits,ebprior.glhits) #posterior prob on models for gl hits
#pp.glhits.collapse =  apply(pp.glhits,2,collapse, nsigmaa=length(sigmaa))
#marginal.glhits = marginal.postprobs(pp.glhits, lbf$gamma,length(sigmaa))

prevhits <- gl[gl$annot==1,]
prevhits.pp <- posteriorprob(lbf.bigmat[,gl$annot==1],ebprior.glhits)
prevhits.pp.marginals <- marginal.postprobs(prevhits.pp, lbf$gamma, length(sigmaa))
prevhits.pp.marginals.DirAssoc <- t(prevhits.pp.marginals[[2]])
prevhits.pp.marginals.DirAssoc <- data.frame(cbind(as.character(prevhits$snp), prevhits.pp.marginals.DirAssoc))
prevhits.pp.marginals.DirAssoc[,2] <- as.numeric(as.character(prevhits.pp.marginals.DirAssoc[,2]))
prevhits.pp.marginals.DirAssoc[,3] <- as.numeric(as.character(prevhits.pp.marginals.DirAssoc[,3]))
prevhits.pp.marginals.DirAssoc[,4] <- as.numeric(as.character(prevhits.pp.marginals.DirAssoc[,4]))
prevhits.pp.marginals.DirAssoc[,5] <- as.numeric(as.character(prevhits.pp.marginals.DirAssoc[,5]))
colnames(prevhits.pp.marginals.DirAssoc) <- c("snp", "LDL", "HDL", "TC", "TG")
prevhits.pp.marginals.DirAssoc$snp <- factor(prevhits.pp.marginals.DirAssoc$snp, levels=prevhits[order(prevhits$lbfav),][,4])



#newhitsfailhits2[,1] <- factor(newhitsfailhits2$snp, levels=sub[order(sub$lbfav, decreasing=TRUE),][,4])
#
#> VYY
#           Z.LDL       Z.HDL       Z.TC       Z.TG
#Z = cbind(gl$Z.ldl, gl$Z.hdl, gl$Z.tc, gl$Z.tg)

~~~
> head(lbf$gamma == 0)
      [,1]  [,2] [,3] [,4]
[1,]  TRUE  TRUE TRUE TRUE
[2,] FALSE  TRUE TRUE TRUE
[3,] FALSE  TRUE TRUE TRUE
[4,]  TRUE FALSE TRUE TRUE
[5,] FALSE FALSE TRUE TRUE
[6,] FALSE FALSE TRUE TRUE
> t(lbf$gamma == 0)[,1:10]
     [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10]
[1,] TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE FALSE FALSE  TRUE
[2,] TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE
[3,] TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE
[4,] TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
> str(prevhits.pp.marginals.DirAssoc)
'data.frame':   148 obs. of  5 variables:
 $ snp: Factor w/ 148 levels "rs10019888","rs10102164",..: 29 78 75 13 94 24 136 9 143 103 ...
 $ LDL: num  0.402 1 1 1 0.659 ...
 $ HDL: num  1 0.516 0.833 0.837 1 ...
 $ TC : num  0.0455 0.9998 1 1 0.6842 ...
 $ TG : num  1 1 0.999 1 0.999 ...
> prevhits.pp.marginals.DirAssoc$snp <- factor(prevhits.pp.marginals.DirAssoc$snp, levels=prevhits[order(prevhits$lbfav, decreasing=TRUE),][,4])
> 
> str(prevhits.pp.marginals.DirAssoc)
'data.frame':   148 obs. of  5 variables:
 $ snp: Factor w/ 148 levels "rs3764261","rs1532085",..: 6 39 89 38 141 112 106 67 19 85 ...
 $ LDL: num  0.402 1 1 1 0.659 ...
 $ HDL: num  1 0.516 0.833 0.837 1 ...
 $ TC : num  0.0455 0.9998 1 1 0.6842 ...
 $ TG : num  1 1 0.999 1 0.999 ...
> prevhits[order(prevhits$lbfav, decreasing=TRUE),][1:5,]
            posHg18 chr       pos       snp a1 a2     maf beta_LDL se_LDL
3142 chr16:55550825  16  56993324 rs3764261  c  a 0.29420   0.0528 0.0042
7853 chr15:56470658  15  58683366 rs1532085  a  g 0.36680   0.0026 0.0037
9322  chr2:27584444   2  27730940 rs1260326  t  c 0.41290   0.0206 0.0037
4294 chr1:109619829   1 109818306  rs629301  t  g 0.21240   0.1669 0.0049
9903 chr19:11063306  19  11202306 rs6511720  g  t 0.09763   0.2209 0.0061
        n_LDL beta_HDL se_HDL    n_HDL beta_TG  se_TG   n_TG beta_TC  se_TC
3142 164865.0   0.2412 0.0039 177533.0  0.0395 0.0038 169275  0.0503 0.0040
7853 171461.0   0.1068 0.0035 185482.0  0.0310 0.0034 176156  0.0545 0.0036
9322 172995.0   0.0113 0.0035 187062.0  0.1148 0.0034 177765  0.0512 0.0036
4294 142644.0   0.0334 0.0045 155813.0  0.0012 0.0045 146759  0.1340 0.0047
9903 170607.9   0.0249 0.0057 184617.2  0.0084 0.0056 175280  0.1851 0.0059
         n_TC annot gene       Z.tg       Z.tc     Z.ldl      Z.hdl   mvstat
3142 177497.0     1   NA -38.153405 12.2273740 -11.61497  10.442607 2852.227
7853 185639.0     1   NA  29.293909  0.4575163  14.37587   8.740637 1266.253
9322 187251.0     1   NA   3.131865 -5.2514579 -13.61921 -33.038308 1118.931
4294 155873.0     1   NA  -7.168038 33.1514351  27.83363   0.432745 1189.257
9903 184763.8     1   NA  -4.000495 34.5878329  30.32462   1.624355 1253.111
          mvp     unip     nmin     lbfav
3142 616.1989 317.7773 164865.0 1235.6337
7853 272.1611 187.9066 171461.0  703.3767
9322 240.2242 238.6402 172995.0  419.9862
4294 255.4690 240.2676 142644.0  350.7596
9903 269.3120 261.4145 170607.9  329.7005
> sub[order(sub$lbfav, decreasing=TRUE),][1:5,]
             posHg18 chr       pos       snp a1 a2     maf beta_LDL se_LDL
7369  chr19:50103919  19  45412079    rs7412  c  t 0.06596   0.5898 0.0101
11316 chr16:55547091  16  56989590  rs247616  c  t 0.29290   0.0547 0.0041
10290 chr15:56510771  15  58723479 rs1077834  c  t 0.21110   0.0006 0.0045
9322   chr2:27584444   2  27730940 rs1260326  t  c 0.41290   0.0206 0.0037
6800  chr1:109620053   1 109818530  rs646776  t  c 0.21240   0.1602 0.0044
          n_LDL beta_HDL se_HDL     n_HDL beta_TG  se_TG      n_TG beta_TC
7369   82533.09   0.0978 0.0097  92158.17  0.1119 0.0093  86163.71  0.3736
11316 171458.00   0.2430 0.0038 185471.00  0.0393 0.0037 176146.00  0.0499
10290 170313.99   0.1253 0.0041 184349.04  0.0470 0.0041 175012.96  0.0652
9322  172995.00   0.0113 0.0035 187062.00  0.1148 0.0034 177765.00  0.0512
6800  173020.95   0.0340 0.0041 187093.01  0.0034 0.0040 177796.95  0.1272
       se_TC      n_TC annot gene       Z.tg        Z.tc     Z.ldl       Z.hdl
7369  0.0096  92046.16     2   NA  -8.925523 -38.1534046  35.98270 -11.1078442
11316 0.0040 185621.00     2   NA -38.153405  12.7651669 -11.78858  10.4752000
10290 0.0043 184503.96     2   NA  28.594923   0.2076529  14.43141  10.9786694
9322  0.0036 187251.00     1   NA   3.131865  -5.2514579 -13.61921 -33.0383084
6800  0.0042 187287.97     2   NA  -7.903032  35.2710904  29.16920   0.8916137
         mvstat       mvp     unip      nmin     lbfav
7369  25597.065 5554.2248 317.7773  82533.09 4933.5851
11316  2961.166  639.8384 317.7773 171458.00 1277.7186
10290  1288.083  276.8941 179.1096 170313.99  735.4238
9322   1118.931  240.2242 238.6402 172995.00  419.9862
6800   1334.150  286.8821 271.7878 173020.95  402.4187
> sub[sub$annot==0,][order(sub[sub$annot==0,]$lbfav, decreasing=TRUE),][1:5,]
             posHg18 chr      pos        snp a1 a2     maf beta_LDL se_LDL
8902  chr11:47269468  11 47312892 rs10838687  t  g 0.20580   0.0001 0.0045
9948   chr6:31456344   6 31348365  rs4143332  g  a 0.09103   0.0191 0.0060
5780  chr19:50709959  19 46018119  rs7255743  g  a 0.02639   0.1647 0.0156
11103  chr6:30452624   6 30344645  rs3132631  c  t 0.09367   0.0103 0.0061
8183   chr6:28365356   6 28257377 rs13211507  c  t 0.07520   0.0034 0.0069
         n_LDL beta_HDL se_HDL    n_HDL beta_TG  se_TG     n_TG beta_TC  se_TC
8902  163382.0   0.0517 0.0042 177336.0  0.0136 0.0041 167995.0  0.0167 0.0044
9948  169940.2   0.0262 0.0056 182038.8  0.0464 0.0056 174674.0  0.0443 0.0058
5780  123958.8   0.0027 0.0146 137200.5  0.0323 0.0138 127634.9  0.1066 0.0150
11103 161088.4   0.0244 0.0057 173780.5  0.0338 0.0056 165502.5  0.0303 0.0059
8183  167622.0   0.0272 0.0064 180353.0  0.0229 0.0063 172118.0  0.0197 0.0067
          n_TC annot gene        Z.tg       Z.tc     Z.ldl     Z.hdl    mvstat
8902  177467.0     0   NA -11.9551817  0.4565424 -3.961053  2.897708 177.09359
9948  184164.8     0   NA  -4.4538348 -3.1038812 -7.432953 -8.104553 128.17357
5780  135895.3     0   NA  -0.6782708 10.4537794  7.135980 -2.624603 129.01266
11103 173763.6     0   NA  -3.8862611 -1.5280505 -4.781260 -6.020109  72.31703
8183  180382.0     0   NA   4.3977080 -0.6128130  2.909247  3.718764  54.89740
           mvp      unip     nmin    lbfav
8902  36.50333 32.214670 163382.0 67.08808
9948  26.01904 15.276216 169940.2 61.74188
5780  26.19846 24.851397 123958.8 37.46776
11103 14.13339  8.758703 161088.4 34.72637
8183  10.46676  4.960983 167622.0 30.07076
> sub[order(sub$lbfav, decreasing=TRUE),][sub$annot==0,][1:5,]
             posHg18 chr      pos       snp a1 a2     maf beta_LDL se_LDL
7369  chr19:50103919  19 45412079    rs7412  c  t 0.06596   0.5898 0.0101
10290 chr15:56510771  15 58723479 rs1077834  c  t 0.21110   0.0006 0.0045
9322   chr2:27584444   2 27730940 rs1260326  t  c 0.41290   0.0206 0.0037
9903  chr19:11063306  19 11202306 rs6511720  g  t 0.09763   0.2209 0.0061
11152  chr1:62906518   1 63133930 rs4587594  g  a 0.31000   0.0493 0.0039
          n_LDL beta_HDL se_HDL     n_HDL beta_TG  se_TG      n_TG beta_TC
7369   82533.09   0.0978 0.0097  92158.17  0.1119 0.0093  86163.71  0.3736
10290 170313.99   0.1253 0.0041 184349.04  0.0470 0.0041 175012.96  0.0652
9322  172995.00   0.0113 0.0035 187062.00  0.1148 0.0034 177765.00  0.0512
9903  170607.90   0.0249 0.0057 184617.21  0.0084 0.0056 175280.04  0.1851
11152 173007.00   0.0147 0.0036 187073.90  0.0694 0.0035 177772.00  0.0754
       se_TC      n_TC annot gene      Z.tg        Z.tc     Z.ldl      Z.hdl
7369  0.0096  92046.16     2   NA -8.925523 -38.1534046  35.98270 -11.107844
10290 0.0043 184503.96     2   NA 28.594923   0.2076529  14.43141  10.978669
9322  0.0036 187251.00     1   NA  3.131865  -5.2514579 -13.61921 -33.038308
9903  0.0059 184763.78     1   NA -4.000495  34.5878329  30.32462   1.624355
11152 0.0037 187260.00     2   NA -3.872782 -11.8730499 -19.06773 -19.202855
          mvstat       mvp      unip      nmin     lbfav
7369  25597.0648 5554.2248 317.77728  82533.09 4933.5851
10290  1288.0829  276.8941 179.10958 170313.99  735.4238
9322   1118.9309  240.2242 238.64016 172995.00  419.9862
9903   1253.1114  269.3120 261.41454 170607.90  329.7005
11152   587.4199  125.0872  81.45556 173007.00  241.0438
> sub[order(sub$lbfav, decreasing=TRUE),][sub[order(sub$lbfav, decreasing=TRUE),]$annot==0,][1:5,]
             posHg18 chr      pos        snp a1 a2     maf beta_LDL se_LDL
8902  chr11:47269468  11 47312892 rs10838687  t  g 0.20580   0.0001 0.0045
9948   chr6:31456344   6 31348365  rs4143332  g  a 0.09103   0.0191 0.0060
5780  chr19:50709959  19 46018119  rs7255743  g  a 0.02639   0.1647 0.0156
11103  chr6:30452624   6 30344645  rs3132631  c  t 0.09367   0.0103 0.0061
8183   chr6:28365356   6 28257377 rs13211507  c  t 0.07520   0.0034 0.0069
         n_LDL beta_HDL se_HDL    n_HDL beta_TG  se_TG     n_TG beta_TC  se_TC
8902  163382.0   0.0517 0.0042 177336.0  0.0136 0.0041 167995.0  0.0167 0.0044
9948  169940.2   0.0262 0.0056 182038.8  0.0464 0.0056 174674.0  0.0443 0.0058
5780  123958.8   0.0027 0.0146 137200.5  0.0323 0.0138 127634.9  0.1066 0.0150
11103 161088.4   0.0244 0.0057 173780.5  0.0338 0.0056 165502.5  0.0303 0.0059
8183  167622.0   0.0272 0.0064 180353.0  0.0229 0.0063 172118.0  0.0197 0.0067
          n_TC annot gene        Z.tg       Z.tc     Z.ldl     Z.hdl    mvstat
8902  177467.0     0   NA -11.9551817  0.4565424 -3.961053  2.897708 177.09359
9948  184164.8     0   NA  -4.4538348 -3.1038812 -7.432953 -8.104553 128.17357
5780  135895.3     0   NA  -0.6782708 10.4537794  7.135980 -2.624603 129.01266
11103 173763.6     0   NA  -3.8862611 -1.5280505 -4.781260 -6.020109  72.31703
8183  180382.0     0   NA   4.3977080 -0.6128130  2.909247  3.718764  54.89740
           mvp      unip     nmin    lbfav
8902  36.50333 32.214670 163382.0 67.08808
9948  26.01904 15.276216 169940.2 61.74188
5780  26.19846 24.851397 123958.8 37.46776
11103 14.13339  8.758703 161088.4 34.72637
8183  10.46676  4.960983 167622.0 30.07076
> sub[order(sub$lbfav, decreasing=TRUE),][1:5,]
             posHg18 chr       pos       snp a1 a2     maf beta_LDL se_LDL
7369  chr19:50103919  19  45412079    rs7412  c  t 0.06596   0.5898 0.0101
11316 chr16:55547091  16  56989590  rs247616  c  t 0.29290   0.0547 0.0041
10290 chr15:56510771  15  58723479 rs1077834  c  t 0.21110   0.0006 0.0045
9322   chr2:27584444   2  27730940 rs1260326  t  c 0.41290   0.0206 0.0037
6800  chr1:109620053   1 109818530  rs646776  t  c 0.21240   0.1602 0.0044
          n_LDL beta_HDL se_HDL     n_HDL beta_TG  se_TG      n_TG beta_TC
7369   82533.09   0.0978 0.0097  92158.17  0.1119 0.0093  86163.71  0.3736
11316 171458.00   0.2430 0.0038 185471.00  0.0393 0.0037 176146.00  0.0499
10290 170313.99   0.1253 0.0041 184349.04  0.0470 0.0041 175012.96  0.0652
9322  172995.00   0.0113 0.0035 187062.00  0.1148 0.0034 177765.00  0.0512
6800  173020.95   0.0340 0.0041 187093.01  0.0034 0.0040 177796.95  0.1272
       se_TC      n_TC annot gene       Z.tg        Z.tc     Z.ldl       Z.hdl
7369  0.0096  92046.16     2   NA  -8.925523 -38.1534046  35.98270 -11.1078442
11316 0.0040 185621.00     2   NA -38.153405  12.7651669 -11.78858  10.4752000
10290 0.0043 184503.96     2   NA  28.594923   0.2076529  14.43141  10.9786694
9322  0.0036 187251.00     1   NA   3.131865  -5.2514579 -13.61921 -33.0383084
6800  0.0042 187287.97     2   NA  -7.903032  35.2710904  29.16920   0.8916137
         mvstat       mvp     unip      nmin     lbfav
7369  25597.065 5554.2248 317.7773  82533.09 4933.5851
11316  2961.166  639.8384 317.7773 171458.00 1277.7186
10290  1288.083  276.8941 179.1096 170313.99  735.4238
9322   1118.931  240.2242 238.6402 172995.00  419.9862
6800   1334.150  286.8821 271.7878 173020.95  402.4187
> sub.pp.marginals.DirAssoc[1:5,]
        snp       LDL       HDL        TC        TG
1 rs1905053 0.9959086 0.1609070 0.1684652 0.9999986
2 rs2814982 1.0000000 0.5156295 0.9997797 1.0000000
3  rs400058 0.8840995 0.8560100 0.9763497 0.9999967
4 rs4234589 0.9981705 0.9123257 0.8459061 0.9999713
5   rs38855 0.6588561 0.9999940 0.6842426 0.9990371
> str(sub.pp.marginals.DirAssoc)
'data.frame':   260 obs. of  5 variables:
 $ snp: Factor w/ 260 levels "rs10019888","rs1004712",..: 101 138 161 167 158 46 181 159 234 50 ...
 $ LDL: num  0.996 1 0.884 0.998 0.659 ...
 $ HDL: num  0.161 0.516 0.856 0.912 1 ...
 $ TC : num  0.168 1 0.976 0.846 0.684 ...
 $ TG : num  1 1 1 1 0.999 ...
> sub.pp.marginals.DirAssoc[1:5,]
          snp LDL        HDL         TC        TG
165    rs7412   1 1.00000000 1.00000000 1.0000000
252  rs247616   1 1.00000000 1.00000000 1.0000000
231 rs1077834   1 1.00000000 0.07162226 1.0000000
211 rs1260326   1 1.00000000 0.99998974 0.9999927
149  rs646776   1 0.08622262 1.00000000 1.0000000
> sub.pp.marginals.DirAssoc$snp <- factor(sub.pp.marginals.DirAssoc$snp, levels=sub[order(sub$lbfav, decreasing=TRUE),][,4])
> str(sub.pp.marginals.DirAssoc)
'data.frame':   260 obs. of  5 variables:
 $ snp: Factor w/ 260 levels "rs7412","rs247616",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ LDL: num  1 1 1 1 1 ...
 $ HDL: num  1 1 1 1 0.0862 ...
 $ TC : num  1 1 0.0716 1 1 ...
 $ TG : num  1 1 1 1 1 ...
> str(sub.pp.marginals.DirAssoc.nonGWASregions)
'data.frame':   106 obs. of  5 variables:
 $ snp: Factor w/ 260 levels "rs7412","rs247616",..: 101 85 97 210 161 252 132 184 253 65 ...
 $ LDL: num  0.996 0.884 0.998 1 0.833 ...
 $ HDL: num  0.161 0.856 0.912 0.981 0.655 ...
 $ TC : num  0.168 0.976 0.846 0.995 0.765 ...
 $ TG : num  1 1 1 0.989 1 ...
> sub.pp.marginals.DirAssoc.nonGWASregions <- sub.pp.marginals.DirAssoc[sub$annot==0,]
> sub.pp.marginals.DirAssoc.nonGWASregions$snp <- factor(sub.pp.marginals.DirAssoc.nonGWASregions$snp, levels=sub[sub$annot==0,][order(sub[sub$annot==0,]$lbfav, decreasing=TRUE),][,4]) 
> str(sub.pp.marginals.DirAssoc.nonGWASregions)
'data.frame':   106 obs. of  5 variables:
 $ snp: Factor w/ 106 levels "rs10838687","rs4143332",..: NA NA NA NA NA NA NA NA NA NA ...
 $ LDL: num  1 1 1 1 1 ...
 $ HDL: num  1 1 1 0.152 1 ...
 $ TC : num  1 0.0716 1 1 1 ...
 $ TG : num  1 1 1 1 1 ...
> head(sub.pp.marginals.DirAssoc.nonGWASregions)
     snp LDL       HDL         TC        TG
165 <NA>   1 1.0000000 1.00000000 1.0000000
231 <NA>   1 1.0000000 0.07162226 1.0000000
211 <NA>   1 1.0000000 0.99998974 0.9999927
221 <NA>   1 0.1519397 1.00000000 1.0000000
248 <NA>   1 1.0000000 1.00000000 0.9999837
176 <NA>   1 1.0000000 0.65299897 1.0000000
> sub.pp.marginals.DirAssoc.nonGWASregions <- sub.pp.marginals.DirAssoc[sub[order(sub$lbfav, decreasing=TRUE),]$annot==0,]
> sub.pp.marginals.DirAssoc.nonGWASregions$snp <- factor(sub.pp.marginals.DirAssoc.nonGWASregions$snp, levels=sub[sub$annot==0,][order(sub[sub$annot==0,]$lbfav, decreasing=TRUE),][,4])
> head(sub.pp.marginals.DirAssoc.nonGWASregions)
           snp       LDL       HDL        TC        TG
206 rs10838687 0.9998738 0.7160358 0.2597294 1.0000000
222  rs4143332 1.0000000 1.0000000 0.9752928 0.9999925
129  rs7255743 1.0000000 0.8667669 1.0000000 0.9999454
247  rs3132631 0.9999916 0.9999997 0.6598628 0.9998294
186 rs13211507 0.9949236 0.9889688 0.4327028 0.9999861
230  rs1473886 1.0000000 0.9999616 0.9998611 0.9999929
> str(sub.pp.marginals.DirAssoc.nonGWASregions)
'data.frame':   106 obs. of  5 variables:
 $ snp: Factor w/ 106 levels "rs10838687","rs4143332",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ LDL: num  1 1 1 1 0.995 ...
 $ HDL: num  0.716 1 0.867 1 0.989 ...
 $ TC : num  0.26 0.975 1 0.66 0.433 ...
 $ TG : num  1 1 1 1 1 ...
> sub.pp.marginals.DirAssoc[sub.pp.marginals.DirAssoc$snp=="rs3770818",]
          snp       LDL       HDL        TC        TG
257 rs3770818 0.9723532 0.9100756 0.9491748 0.9257847
> sub[sub$snp=="rs3770818",]
            posHg18 chr      pos       snp a1 a2   maf beta_LDL se_LDL n_LDL
11446 chr2:36623001   2 36769497 rs3770818  t  g 0.281   0.0207 0.0059 89888
      beta_HDL se_HDL n_HDL beta_TG  se_TG  n_TG beta_TC  se_TC  n_TC annot
11446   0.0175 0.0054 94311  0.0087 0.0052 91013  0.0218 0.0057 94595     0
      gene      Z.tg      Z.tc     Z.ldl   Z.hdl   mvstat      mvp     unip
11446   NA -2.883246 -4.282506 -4.437381 1.64529 33.94603 6.116666 5.040672
       nmin    lbfav
11446 89888 2.641615
~~~

#20160901 from http://www.sthda.com/english/wiki/ggplot2-rotate-a-graph-reverse-and-flip-the-plot 
~~~
Horizontal plot : coord_flip()

Box plot :
library(ggplot2)
# Basic box plot
bp <- ggplot(PlantGrowth, aes(x=group, y=weight))+
  geom_boxplot()
bp
# Horizontal box plot
bp + coord_flip()
.
.
.
everse y axis

The function scale_y_reverse() can be used as follow :

# Basic histogram
hp
# Y axis reversed
hp + scale_y_reverse()
~~~












#melt(t(sub.pp.marginals.DirAssoc[1:4,1:5]))
#ggplot(t(sub.pp.marginals.DirAssoc[1:4,1:5])) 
#ggplot(melt(prevhits2[,c(1,5:8)]), aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue")


#jpeg("GlobalLipids2013.vs3.ForCSHLPoster.20160830PreProbGenom.vs1.jpg", width=2000, height=4000, res=300)
jpeg("GlobalLipids2013.vs3.ForCSHLPoster.20160830PreProbGenom.MarginalHeatplots.NewHits.vs1.jpg", width=2000, height=4000, res=300)

#ggplot(melt(sub.pp.marginals.DirAssoc), aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue")
#ggplot(melt(newhits[,c(1,5:8)]), aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue")
ggplot(melt(sub.pp.marginals.DirAssoc.newhits), aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue")

dev.off()

jpeg("GlobalLipids2013.vs3.ForCSHLPoster.20160830PreProbGenom.MarginalHeatplots.PrevHits.vs1.jpg", width=2000, height=4000, res=300)

ggplot(melt(prevhits.pp.marginals.DirAssoc), aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue")

dev.off()

jpeg("GlobalLipids2013.vs3.ForCSHLPoster.20160830PreProbGenom.MarginalHeatplots.MargSNPs.vs1.jpg", width=2000, height=4000, res=300)

ggplot(melt(sub.pp.marginals.DirAssoc), aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue")

dev.off()

jpeg("GlobalLipids2013.vs3.ForCSHLPoster.20160830PreProbGenom.MarginalHeatplots.MargSNPs.nonGWASregions.vs1.jpg", width=2000, height=4000, res=300)

ggplot(melt(sub.pp.marginals.DirAssoc.nonGWASregions), aes(variable, snp, value)) + geom_tile(aes(fill=value), colour="white") + scale_fill_gradient(low = "white", high="steelblue")

dev.off()

#Get marginals per alphasigma level first then combine across sigmas












