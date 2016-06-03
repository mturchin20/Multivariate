set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.RSS0.vs1.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtsignif.vs1.annot.vs1.txt.gz", header=T)
#gl3 <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
gl <- gl[!is.na(gl$maf) & gl$maf > 0 & !is.na(gl$chrbp) & gl$chrbp != "None",]
Z = cbind(gl$Z.height,gl$Z.BMI,gl$Z.WHRadjBMI,gl$Z.WHR,gl$Z.HIPadjBMI,gl$Z.HIP,gl$Z.WCadjBMI,gl$Z.WC)

#~~~

#~~~


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

n = cbind(gl$n_height, gl$n_BMI, gl$n_WHRadjBMI, gl$n_WHR, gl$n_HIPadjBMI, gl$n_HIP, gl$n_WCadjBMI, gl$n_WC)
n = apply(n,1,min)
gl$nmin=n
sigmaa=c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)
#sigmaa_2=c(0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15)
#sigmaa_2=c(0.005)
#lbf=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa)
#lbf=compute.allBFs.fromZscores.Efficient(Z,VYY,gl$nmin,gl$maf,sigmaa)
#lbf.bigmat=do.call(rbind,lbf$lbf) #lbf$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows 


lbf_sig1=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[1])
lbf_sig2=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[2])
lbf.bigmat <- do.call(rbind, c(lbf_sig1$lbf, lbf_sig2$lbf))
rm(lbf_sig1, lbf_sig2)
lbf_sig3=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[3])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig3$lbf))
rm(lbf_sig3)
lbf_sig4=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[4])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig4$lbf))
rm(lbf_sig4)
lbf_sig5=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[5])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig5$lbf))
rm(lbf_sig5)
lbf_sig6=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[6])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig6$lbf))
rm(lbf_sig6)
lbf_sig7=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[7])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig7$lbf))
rm(lbf_sig7)
lbf_sig8=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[8])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig8$lbf))
rm(lbf_sig8)
lbf_sig9=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[9])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig9$lbf))
rm(lbf_sig9)
lbf_sig10=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[10])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig10$lbf))
rm(lbf_sig10)
lbf_sig11=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[11])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig11$lbf))
rm(lbf_sig11)
lbf_sig12=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[12])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig12$lbf))
rm(lbf_sig12)
lbf_sig13=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[13])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig13$lbf))
rm(lbf_sig13)
lbf_sig14=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa[14])
lbf.bigmat <- do.call(rbind, c(lbf.bigmat, lbf_sig14$lbf))
rm(lbf_sig14)

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

lbf.bigmatSub1<- rbind(lbf.bigmat1, lbf.bigmat2)
rm(lbf.bigmat1, lbf.bigmat2)
lbf.bigmatSub2<- rbind(lbf.bigmat3, lbf.bigmat4)
rm(lbf.bigmat3, lbf.bigmat4)
lbf.bigmatSub3<- rbind(lbf.bigmat5, lbf.bigmat6)
rm(lbf.bigmat5, lbf.bigmat6)


lbf.bigmat <- rbind(lbf.bigmatSub1, lbf.bigmatSub2, lbf.bigmatSub3, lbf.bigmat7)

rm(lbf.bigmatSub1, lbf.bigmatSub2, lbf.bigmatSub3, lbf.bigmat7)

lbf.bigmat<- rbind(lbf.bigmat1, lbf.bigmat2, lbf.bigmat3, lbf.bigmat4, lbf.bigmat5, lbf.bigmat6, lbf.bigmat7)


#lbf.bigmat1=rbind(lbf_sig1$lbf, lbf_sig2$lbf, lbf_sig3$lbf, lbf_sig4$lbf, lbf_sig5$lbf)
#rm(lbf_sig2, lbf_sig3, lbf_sig4, lbf_sig5)
#lbf.bigmat1 <- rbind(lbf.bigmat1, lbf_sig6$lbf, lbf_sig7$lbf, lbf_sig8$lbf, lbf_sig9$lbf, lbf_sig10$lbf)

#lbf.bigmat=rbind(lbf_sig1$lbf, lbf_sig2$lbf, lbf_sig3$lbf, lbf_sig4$lbf, lbf_sig5$lbf, lbf_sig6$lbf, lbf_sig7$lbf, lbf_sig8$lbf, lbf_sig9$lbf, lbf_sig10$lbf, lbf_sig11$lbf, lbf_sig12$lbf, lbf_sig13$lbf, lbf_sig14$lbf)

lbf=list(gamma=lbf_sig1$gamma,prior=lbf__sig1$prior)

#collapse takes a vector that is nsigmmaa stacked m-vectors, and adds them together to produce a single m vector (averages over values of sigmaa)
collapse = function(x,nsigmaa){return(apply(matrix(x,ncol=nsigmaa),1,sum))}

#install.packages("/mnt/lustre/home/mturchin20/Software/bigmemory.sri_0.1.3.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/BH_1.58.0-1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/bigmemory_4.4.6.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/stringi_0.4-1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/magrittr_1.5.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/stringr_1.0.0.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/Rcpp_0.11.6.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/pryr_0.1.1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/plyr_1.8.2.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/reshape2_1.4.1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/chron_2.3-45.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/data.table_1.9.4.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/assertthat_0.1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/R6_2.0.1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/lazyeval_0.1.10.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/DBI_0.3.1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/dplyr_0.4.1.tar.gz", repos = NULL, type="source")
#install.packages("/mnt/lustre/home/mturchin20/Software/", repos = NULL, type="source")

library(bigmemory)
#library(pryr)

compute.allBFs.fromZscores.Efficient = function(Z,VYY,n,f,sigmaa,pi0=0.5,m=0){
	d = dim(Z)[2]
	p = dim(Z)[1]
	if(m==0){m = d-1}
	VXX = 2*f *(1-f)
	VYX = t(sqrt(VXX)*Z/sqrt(n))

	prior = rep(0,3^d)
	gamma=as.big.matrix(matrix(0,nrow=3^d,ncol=d))
	lbf=list()
	for(ss in 1:length(sigmaa)){
		lbf[[ss]] = as.big.matrix(matrix(0,nrow=3^d, ncol=p))
	}

	for(i in 0:(3^d-1)){
		for(j in 1:d){
			gamma[i+1,j]= (i %% 3^j) %/% 3^{j-1}
		}
		prior[i+1] = computeprior(gamma[i+1,],pi0)
		U = (gamma[i+1,]==0)
		D = (gamma[i+1,]==1)
		#print(U); print(D);
		if(prior[i+1]>0){
			for(ss in 1:length(sigmaa)){
				lbf[[ss]][i+1,] = logBF.fromVSummaries(VYX,VYY,VXX,U,D,n,m,d,sigmaa[ss])/log(10)
				#note we just don't bother computing for models with prior = 0
			}
		} else {
			for(ss in 1:length(sigmaa)){
				lbf[[ss]][i+1,] = 0
			}
		}
	}
	prior[1] = pi0
	return(list(lbf=lbf,gamma=gamma,prior=prior))
}



#compute empirical bayes priors from GWAS hit SNPs
gl.glhits = gl[gl$annot==1,] # subset of GWAS SNPs
#CHECK_0: Note -- some top hits dropped because they did not have 'maf'; expecting ~881 below, not 843
#~~~
#> dim(gl.glhits)
#[1] 843  44
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

write.table(file="GIANT2014_5.AllPheno.newtophits.vs1.txt",cbind(sub[sub$lbfav>2.89926 & sub$nmin>50000 & sub$annot==0,c(1:3,7)],round(sub[sub$lbfav>2.89926 & sub$nmin>50000 & sub$annot==0,c(18:20,25)],digits=4)),quote=FALSE,sep= " ", row.names=FALSE)

#newhits = sub[sub$annot==0 & sub$lbfav>3.323736 & sub$nmin>50000,c(1:3,7,18:20,25)]

#newhits = gl[gl$annot==0 & gl$lbfav>3.323736 & gl$nmin>50000,c(4,2:3,30)]
#write.table(file="newhits.vs2.wr",newhits,quote=FALSE,row.names=FALSE)










