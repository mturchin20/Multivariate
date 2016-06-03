set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.RSS0.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$chrbp) & gl$chrbp != "None",]
Z = cbind(gl$Z.height,gl$Z.BMI,gl$Z.WHRadjBMI)

#~~~
#> gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz", header=T)
#> dim(gl)
#[1] 42022    23
#> gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$chrbp) & gl$chrbp != "None",]
#> dim(gl)
#[1] 42008    23
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
#Z.height      Z.BMI Z.WHRadjBMI    mvstat       mvp      unip
#1  5.099807 0.02080652   0.6128130  26.31993  4.564237  6.468521
#2 11.053855 2.53514746   1.0581216 129.25344 26.249946 27.677781
#3  8.186023 2.63636313   0.0000000  73.79878 14.446568 15.568636
#4  5.630201 1.36104332   0.2404260  33.52077  6.029490  7.744727
#5  5.831292 1.86665344   0.1891184  37.41533  6.830000  8.259637
#6  5.840642 1.80000404   0.2404260  37.29792  6.805801  8.283997
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
#CHECK_0: Note -- some top hits dropped because they did not have 'maf'; expecting ~861 below, not 821
#~~~
#> dim(gl.glhits)
#[1] 821  24
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

sum(modelmatrix[allassoc,4]) #0.9574025
sum(modelmatrix[allbut1assoc,4]) #0.03937951

#~~~
#> sum(modelmatrix[allassoc,4])
#[1] 0.9574025
#> sum(modelmatrix[allbut1assoc,4])
#[1] 0.03937951
#~~~

#look at weight on each of WHRadjBMI, BMI, height being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4]) 
# 0.01383266, 0.02330547, 0.002241376

#~~~
#> sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,4])
#[1] 0.01383266
#> sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,4])
#[1] 0.02330547
#> sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,4])
#[1] 0.002241376
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
#[1] 2.89926
#~~~


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
#o = order(gl$chr, gl$pos)
#gl = gl[o,]


sub = gl[gl$nmin>50000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
newhits = sub[sub$annot==0 & sub$lbfav>2.89926 & sub$nmin>50000,c(1:3,7,18:20,25)]

#~~~
#> sub = gl[gl$nmin>50000,]
#> dim(sub)
#[1] 42008    25
#> sub=sub[l==1,]
#> dim(sub)
#[1] 775  25
#~~~

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>2.89926 & sub$nmin>50000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notWHRadjBMI, or notBMI or notheight)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#194 in results
#~~~
#> as.character(newhits$snp)
#[1] "rs6719451"  "rs7651617"  "rs2741856"  "rs2347709"  "rs10737541"
#[6] "rs9424148"  "rs11135932" "rs7492628"  "rs4934353"  "rs1344632"
#[11] "rs12875433" "rs2239535"  "rs2481665"  "rs4666483"  "rs11594532"
#[16] "rs751008"   "rs12247907" "rs675144"   "rs948847"   "rs175425"
#[21] "rs11615028" "rs11170468" "rs4738873"  "rs10495563" "rs11924886"
#[26] "rs1856370"  "rs823114"   "rs17469356" "rs3740237"  "rs2727270"
#[31] "rs628751"   "rs4690014"  "rs7398691"  "rs17784714" "rs162965"
#[36] "rs6990042"  "rs10111718" "rs1470374"  "rs4277572"  "rs12328935"
#[41] "rs13097792" "rs10883035" "rs4766230"  "rs38902"    "rs12139158"
#[46] "rs2615075"  "rs5749202"  "rs1248690"  "rs3801427"  "rs11099956"
#[51] "rs2322193"  "rs10882160" "rs10956525" "rs4618485"  "rs12150665"
#[56] "rs6976841"  "rs7729880"  "rs2219993"  "rs504560"   "rs2281514"
#[61] "rs9288807"  "rs2286046"  "rs2270894"  "rs1793249"  "rs6471763"
#[66] "rs1572413"  "rs17375290" "rs6870983"  "rs11104829" "rs2970357"
#[71] "rs2014590"  "rs10739570" "rs1333916"  "rs1372449"  "rs12882850"
#[76] "rs7200543"  "rs811153"   "rs11708067" "rs5017948"  "rs290195"
#[81] "rs365836"   "rs1043595"  "rs4444350"  "rs6708784"  "rs13007001"
#[86] "rs12603582" "rs1869868"  "rs11766345" "rs11654954" "rs1396214"
#[91] "rs6864049"  "rs235292"   "rs2360960"  "rs10792252" "rs17177887"
#[96] "rs13284412" "rs10262697" "rs10760279" "rs11051005" "rs6559588"
#[101] "rs6807897"  "rs2049974"  "rs8031625"  "rs11167042" "rs17630235"
#[106] "rs10510128" "rs6542180"  "rs12282824" "rs1020548"  "rs2310357"
#[111] "rs2303085"  "rs2282611"  "rs12770214" "rs12545970" "rs1816537"
#[116] "rs285575"   "rs432037"   "rs4321967"  "rs17731432" "rs181362"
#[121] "rs4121977"  "rs9908143"  "rs1785702"  "rs2068462"  "rs8123881"
#[126] "rs1074287"  "rs10467367" "rs12339475" "rs6546568"  "rs9893015"
#[131] "rs10950356" "rs10965269" "rs1644634"  "rs12496416" "rs946197"
#[136] "rs3790729"  "rs6903903"  "rs17813106" "rs4447843"  "rs13014679"
#[141] "rs1255504"  "rs4666070"  "rs929641"   "rs2488003"  "rs12601850"
#[146] "rs972540"   "rs7797205"  "rs2405983"  "rs7611238"  "rs6086034"
#[151] "rs11726922" "rs4401606"  "rs1075698"  "rs12153375" "rs11008282"
#[156] "rs11045889" "rs1304463"  "rs4922031"  "rs13361606" "rs4984406"
#[161] "rs1288380"  "rs990211"   "rs11787111" "rs11743851" "rs1857752"
#[166] "rs2222413"  "rs2270204"  "rs288193"   "rs4551650"  "rs7017861"
#[171] "rs7947391"  "rs1370394"  "rs350887"   "rs7163907"  "rs8057807"
#[176] "rs900178"   "rs12191682" "rs9507983"  "rs1465821"  "rs1714383"
#[181] "rs7703022"  "rs10971721" "rs17767294" "rs11183167" "rs17630640"
#[186] "rs17014035" "rs150992"   "rs16893186" "rs17301622" "rs2058051"
#[191] "rs12624843" "rs935728"   "rs12523158" "rs7617596"
#> bestclass
#[1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[38] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[75] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[112] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[149] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#[186] 1 1 1 1 1 1 1 1 1
#> apply(pp.newhits.classmatrix,2,max)
#[1] 0.9469624 0.9395570 0.9421921 0.9451274 0.9577207 0.9412928 0.9397165
#[8] 0.9997772 0.9461989 0.9450618 0.9413021 0.9394811 0.9612050 0.9504490
#[15] 0.9984616 0.9917135 0.9998947 0.9442216 0.9427088 0.9380733 0.9550508
#[22] 0.9801420 0.9943608 0.9960750 0.9503399 0.9554707 0.9946697 0.9496766
#[29] 0.9975098 0.9410904 0.9967279 0.9431784 0.9423323 0.9467280 0.9412098
#[36] 0.9576658 0.9418842 0.9435346 0.9556576 0.9467830 0.9413201 0.9957943
#[43] 0.9530781 0.9979963 0.9445678 0.9766967 0.9419298 0.9422432 0.9696445
#[50] 0.9398143 0.9396037 0.9481592 0.9481916 0.9415437 0.9936665 0.9444344
#[57] 0.9425070 0.9577377 0.9506539 0.9477826 0.9419274 0.9730724 0.9958914
#[64] 0.9396601 0.9973864 0.9500770 0.9400710 0.9961800 0.9392188 0.9528371
#[71] 0.9396433 0.9533637 0.9463143 0.9386868 0.9609857 0.9964937 0.9476491
#[78] 0.9965255 0.9483328 0.9491552 0.9422796 0.9457736 0.9554793 0.9547151
#[85] 0.9434625 0.9444213 0.9441454 0.9976559 0.9610070 0.9999680 0.9597367
#[92] 0.9406636 0.9438881 0.9413222 0.9878613 0.9395708 0.9492613 0.9941207
#[99] 0.9971584 0.9508683 0.9551730 0.9394474 0.9412728 0.9385537 0.9625536
#[106] 0.9460490 0.9619768 0.9398861 0.9632684 0.9399706 0.9484233 0.9412539
#[113] 0.9497841 0.9452230 0.9934516 0.9598152 0.9816674 0.9999494 0.9657346
#[120] 0.9980752 0.9589702 0.9412696 0.9426978 0.9543860 0.9704273 0.9642260
#[127] 0.9417928 0.9986000 0.9462351 0.9848101 0.9524631 0.9444054 0.9623679
#[134] 0.9434728 0.9469063 0.9406878 0.9414425 0.9517826 0.9461573 0.9410797
#[141] 0.9408420 0.9465745 0.9990783 0.9435724 0.9864703 0.9646327 0.9409583
#[148] 0.9462190 0.9841938 0.9413167 0.9417313 0.9637907 0.9393196 0.9406568
#[155] 0.9821392 0.9421412 0.9390136 0.9397773 0.9401868 0.9783352 0.9462376
#[162] 0.9997583 0.9985108 0.9905839 0.9424003 0.9973442 0.9961488 0.9973391
#[169] 0.9444953 0.9521341 0.9509595 0.9862950 0.9535038 0.9656809 0.9401161
#[176] 0.9416789 0.9400099 0.9875251 0.9517646 0.9421261 0.9412000 0.9985334
#[183] 0.9391406 0.9403965 0.9954847 0.9410925 0.9738198 0.9395741 0.9634972
#[190] 0.9407671 0.9641745 0.9455318 0.9505616 0.9427241
#~~~

jpeg("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.Orig3.vs1.ForCSHLPoster.ModelsVslbfavNewOld.vs1.jpeg", res=300, width=5000, height=4000)
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
#1      1   2         2 0.56563014 0.5656301
#2      1   1         1 0.17533016 0.7409603
#3      1   1         2 0.10816493 0.8491252
#4      1   2         1 0.07529206 0.9244173
#5      2   1         2 0.03298522 0.9574025
#6      1   0         2 0.02290875 0.9803113
#> modelmatrix.paste <- apply(modelmatrix[,1:3], 1, function(x) paste(x, collapse="_"))
#> table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)])
#
#1_1_1 1_1_2 1_2_1 1_2_2 2_1_2
#23    15     8   141     7
#> modelmatrix[apply(as.matrix(names(table(modelmatrix.paste[apply(ppmatrix.newhits, 2, which.max)]))), 1, function(x) which(modelmatrix.paste == x)),4]
#[1] 0.17533016 0.10816493 0.07529206 0.56563014 0.03298522
#>
#> topmodels
#i
#[1,] "1_2_2" "141" "0.769816970076792" "0.109224690960025"  "0.678534152428078"
#[2,] "1_1_1" "23"  "0.759929776427074" "0.180319136622165"  "0.565396692070456"
#[3,] "1_1_2" "15"  "0.595495133308273" "0.0912283440372175" "0.300891497929458"
#[4,] "1_2_1" "8"   "0.523346662202924" "0.0832780037794609" "0.124097799210174"
#[5,] "2_1_2" "7"   "0.511305809867544" "0.0584819678380476" "0.226833702383289"
#
#[1,] "0.174652253486766" "0.56563014458915"
#[2,] "0.343429769778505" "0.175330164601129"
#[3,] "0.108744537236256" "0.108164929727224"
#[4,] "0.115974475966936" "0.0752920582366282"
#[5,] "0.068073227470706" "0.0329852200087393"
#> oldtopmodels
#i                                                                         
#[1,] "1_2_2" "568" "0.799635492482729" "0.113343420673262"  "0.715481427138328"
#[2,] "1_1_1" "114" "0.762419752684886" "0.207035936833435"  "0.590655768413298"
#[3,] "1_1_2" "60"  "0.609247546859145" "0.130735416680685"  "0.331389153026427"
#[4,] "1_2_1" "41"  "0.634267970783804" "0.150761937613095"  "0.32881465849437" 
#[5,] "2_1_2" "38"  "0.519379134729573" "0.0598438809005867" "0.212528879937159"
#
#[1,] "0.188175540462659" "0.56563014458915"  
#[2,] "0.347843634338573" "0.175330164601129" 
#[3,] "0.198744382491746" "0.108164929727224" 
#[4,] "0.25853652381315"  "0.0752920582366282"
#[5,] "0.097177545228877" "0.0329852200087393"
#> summary(lm(apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix.newhits, 2, max) ~ newhits$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.36590 -0.09530  0.04199  0.11068  0.25773
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.729980   0.038010  19.205   <2e-16 ***
#newhits$lbfav 0.001130   0.007277   0.155    0.877
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.1399 on 192 degrees of freedom
#Multiple R-squared:  0.0001257, Adjusted R-squared:  -0.005082
#F-statistic: 0.02413 on 1 and 192 DF,  p-value: 0.8767
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
#-0.59282 -0.19166  0.09772  0.20252  0.39238
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.589582   0.068432   8.616 2.53e-15 ***
#newhits$lbfav 0.001425   0.013101   0.109    0.914
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2518 on 192 degrees of freedom
#Multiple R-squared:  6.16e-05,  Adjusted R-squared:  -0.005146
#F-statistic: 0.01183 on 1 and 192 DF,  p-value: 0.9135
#
#> 
#> summary(lm(apply(ppmatrix, 2, max) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = apply(ppmatrix, 2, max) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min       1Q   Median       3Q      Max
#-0.40344 -0.11952  0.06328  0.11520  0.25077
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.7373541  0.0073499  100.32  < 2e-16 ***
#gl[gl$annot == 1, ]$lbfav 0.0016154  0.0003748    4.31 1.83e-05 ***
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.1517 on 819 degrees of freedom
#Multiple R-squared:  0.02218,   Adjusted R-squared:  0.02099
#F-statistic: 18.58 on 1 and 819 DF,  p-value: 1.829e-05
#
#> summary(lm((apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) { sort(x,partial=length(x)-1)[length(x)-1] })) ~ gl[gl$annot==1,]$lbfav))
#
#Call:
#lm(formula = (apply(ppmatrix, 2, max) - apply(ppmatrix, 2, function(x) {
#sort(x, partial = length(x) - 1)[length(x) - 1]
#})) ~ gl[gl$annot == 1, ]$lbfav)
#
#Residuals:
#Min      1Q  Median      3Q     Max
#-0.6391 -0.2203  0.1246  0.2112  0.3837
#
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.5999276  0.0129403  46.361   <2e-16 ***
#gl[gl$annot == 1, ]$lbfav 0.0020265  0.0006598   3.071   0.0022 **
#---
#Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
#
#Residual standard error: 0.2671 on 819 degrees of freedom
#Multiple R-squared:  0.01139,   Adjusted R-squared:  0.01018
#F-statistic: 9.433 on 1 and 819 DF,  p-value: 0.002202
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

write.table(file="GIANT2014_5.Orig3.newtophits.vs1.txt",cbind(sub[sub$lbfav>2.89926 & sub$nmin>50000 & sub$annot==0,c(1:3,7)],round(sub[sub$lbfav>2.89926 & sub$nmin>50000 & sub$annot==0,c(18:20,25)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)

#newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>50000,c(1:3,7,18:20,25)]

#newhits = gl[gl$annot==0 & gl$lbfav>3.323736 & gl$nmin>50000,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)




#20160428
~~~
#Prior values
> modelmatrix[1:5,]
  height BMI WHRadjBMI          p      cump
  1      1   2         2 0.56563014 0.5656301
  2      1   1         1 0.17533016 0.7409603
  3      1   1         2 0.10816493 0.8491252
  4      1   2         1 0.07529206 0.9244173
  5      2   1         2 0.03298522 0.9574025

#> modelmatrix[1:5,]
#  TC TG HDL LDL         p      cump
#  1  1  1   1   1 0.2384938 0.2384938
#  2  1  2   1   1 0.1538137 0.3923075
#  3  2  1   1   2 0.1359187 0.5282262
#  4  1  1   2   1 0.1169160 0.6451422
#  5  2  1   1   1 0.1061912 0.7513333

newhits = sub[sub$annot==0 & sub$lbfav>2.89926 & sub$nmin>50000,c(1:3,7,18:20,25)]

gl.glhits$lbfav = lbf.av.glhits
gl.glhits$lbfall = lbf.all.glhits
gl.glhits$lbfuni = lbf.uni.glhits

#newhits = sub[sub$annot==0 & sub$lbfav>2.89926 & sub$nmin>50000,c(1:3,7,18:20,25)]
failhits = sub[sub$annot==0 & sub$lbfav<=2.89926, c(1:3,7,18:20,25)]
prevhits = gl.glhits[,c(1:3,7,18:20,25)]

GLC.MTedits.ForGIANT2014_5.Orig3.vs1.ForCSHLPoster
#write.table(file="failtophits.vs3.ForCSHLPoster.txt",cbind(sub[sub$lbfav<=2.89926 & sub$nmin>50000 & sub$annot==0,c(4,2:3,7)],round(sub[sub$lbfav<=2.89926 & sub$nmin>50000 & sub$annot==0,c(22:25,30)],digits=1)),quote=FALSE,sep= " ", row.names=FALSE)


write.table(file="glhits.vs3.ForCSHLPoster.txt",gl.glhits[,c(4,2,3,7,22,23,24,25,30,31,32)],quote=FALSE,row.names=FALSE)

jpeg("GIANT2014_5.vs3.ForCSHLPoster.NHSSlides.jpg", width=4500, height=2250, res=300)
#jpeg("GIANT2014_5.vs3.ForCSHLPoster.NHSSlides.SDvsMax.jpg", width=2000, height=2000, res=300)
#jpeg("GIANT2014_5.vs3.ForCSHLPoster.NHSSlides.ScaledSDvsMax.jpg", width=2250, height=2250, res=300)
#jpeg("GIANT2014_5.vs3.ForCSHLPoster.NHSSlides.ScaledSDvsLBF.jpg", width=2000, height=2000, res=300)

#par(xpd=FALSE, mar=c(6,6,6,6), oma=c(2,2,2,2))
par(mfrow=c(1,2), mar=c(4,5,5,4), oma=c(2.25,2.25,2.25,2.25))

#plot(c(apply(prevhits[,5:7], 1, sd), apply(newhits[,5:7], 1, sd)), c(apply(prevhits[,5:7], 1, max), apply(newhits[,5:7], 1, max)), main="SD vs Max Z-Scores of Old and New Hits", xlab="SD of Zscores", ylab="Max Zscore", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])))
#plot(c(apply(prevhits[,5:7], 1, sd), apply(newhits[,5:7], 1, sd), apply(failhits[,5:7], 1, sd)), c(apply(prevhits[,5:7], 1, max), apply(newhits[,5:7], 1, max), apply(failhits[,5:7], 1, max)), main="SD vs Max Z-Scores of Old and New Hits", xlab="SD of Zscores", ylab="Max Zscore", col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])))

#plot(c(apply(prevhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(newhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(failhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);})), c(apply(prevhits[,5:7], 1, max), apply(newhits[,5:7], 1, max), apply(failhits[,5:7], 1, max)), main="SD(Zscore/MaxZ) vs MaxZ of Prev, New, & Failed Hits", xlab="SD(Zscore/MaxZ)", ylab="Max Zscore", pch=c(16,16,16), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5) 
plot(c(apply(prevhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(newhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);})), c(apply(prevhits[,5:7], 1, max), apply(newhits[,5:7], 1, max)), main="SD(Zscore/MaxZ) vs MaxZ of Prev & New Hits", xlab="SD(Zscore/MaxZ)", ylab="Max Zscore", pch=c(16,16), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5) 

#legend(.075,38.25, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))
legend(.084, 26.8, c("Previous", "New"), col=c("BLUE", "RED"), pch=c(16,16))

#plot(c(prevhits$lbfav, newhits$lbfav), c(apply(prevhits[,5:7], 1, sd), apply(newhits[,5:7], 1, sd)), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])))
#plot(c(prevhits$lbfav, newhits$lbfav, failhits$lbfav), c(apply(prevhits[,5:7], 1, sd), apply(newhits[,5:7], 1, sd), apply(failhits[,5:7], 1, sd)), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])))

#plot(c(apply(prevhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(newhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(failhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);})), c(prevhits$lbfav, newhits$lbfav, failhits$lbfav), main="SD(Zscore/MaxZ) vs LBF of Prev, New, & Failed Hits", xlab="SD(Zscore/MaxZ)", ylab="LBF", ylim=c(2.5,100), pch=c(16,16,16), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1]), rep("GREY", dim(failhits)[1])), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5)
plot(c(apply(prevhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);}), apply(newhits[,5:7], 1, function(x) { val1 <- sd(x/max(x)); return(val1);})), c(prevhits$lbfav, newhits$lbfav), main="SD(Zscore/MaxZ) vs LBF of Prev & New Hits", xlab="SD(Zscore/MaxZ)", ylab="LBF", ylim=c(2.5,100), pch=c(16,16), col=c(rep("BLUE", dim(prevhits)[1]), rep("RED", dim(newhits)[1])), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, cex=1.5)

#legend(.075,100, c("Previous", "New", "Failed"), col=c("BLUE", "RED", "GREY"), pch=c(16,16,16))
legend(.084, 100, c("Previous", "New"), col=c("BLUE", "RED"), pch=c(16,16))

dev.off()


cex








