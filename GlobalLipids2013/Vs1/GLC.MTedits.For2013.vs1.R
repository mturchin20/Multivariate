set.seed(100)
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/multivariate/test.funcs.R")
source("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/multivariate/globallipids/GLCfuncs.R")
VYY = as.matrix(read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1/RSS0.txt",header=T,sep=","))
gl = read.table("/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1/dtlesssignif.annot.txt.gz", header=T)
gl <- gl[!is.na(gl[,ncol(gl)-1]) & !is.infinite(gl[,ncol(gl)-1]) & !is.na(gl$maf) & gl$maf > 0,]
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

sum(modelmatrix[allassoc,5]) #0.8879285
sum(modelmatrix[allbut1assoc,5]) #0.08868659

#look at weight on each of LDL, HDL, TG, TC being unassociated
sum(modelmatrix[allbut1assoc & modelmatrix[,4]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,3]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,2]==0,5]) 
sum(modelmatrix[allbut1assoc & modelmatrix[,1]==0,5]) 
# 0.05604953, 0.005025145, 0.02758528, 2.663935e-05

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix = cbind(pp.glhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.classmatrix = rbind(colSums(ppmatrix[allassoc,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.classmatrix, 2,which.max)
#gg=gl.glhits$gene
gg=gl.glhits$snp
write.table(cbind(as.character(gg[bestclass !=1]),bestclass[bestclass!=1],round(apply(pp.classmatrix,2,max)[bestclass!=1],2)),"gl.bestmodel.txt",quote=F,sep=" & ",row.names=F)

#	[,1]         [,2] [,3]  
#	[1,] "rs12678919" "2"  "0.01"
#	[2,] "rs2412710"  "2"  "0.34"
#	[3,] "rs629301"   "4"  "0.9" 
#	[4,] "rs17145738" "2"  "0.74"
#	[5,] "rs386000"   "2"  "0.39"
#	[6,] "rs1800961"  "4"  "0.97"
#	[7,] "rs1532085"  "2"  "0.74"
#	[8,] "rs16942887" "2"  "0.91"
#	[9,] "rs7241918"  "4"  "0.75"
#	[10,] "rs6511720"  "4"  "0.58"
#	[11,] "rs737337"   "2"  "0.6" 

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
#[1] 4.078788


lbf.av.all = lbf.av(lbf.bigmat,ebprior.glhits)
gl$lbfav = lbf.av.all
o = order(gl$chr, gl$pos)
gl = gl[o,]


sub = gl[gl$nmin>50000,]
l=indephits(sub$lbfav,sub$chr,sub$pos)
sub=sub[l==1,]
#newhits = sub[sub$annot==0 & sub$lbfav>4.3 & sub$nmin>50000,c(1:4,20:23,28)]
newhits = sub[sub$annot==0 & sub$lbfav>4.07 & sub$nmin>50000,c(4,2:3,7,22:25,30)]

#extract lbfs for all the new hits
lbf.newhits= lbf.bigmat[,gl$nmin>50000]
lbf.newhits= lbf.newhits[,l==1]
lbf.newhits= lbf.newhits[,sub$annot==0 & sub$lbfav>4.07 & sub$nmin>50000]
pp.newhits = posteriorprob(lbf.newhits,ebprior.glhits) #posterior prob on models for new hits
pp.newhits.collapse =  apply(pp.newhits,2,collapse, nsigmaa=length(sigmaa))

#compute for each SNP which class it is most likely assigned to
#(all assoc, or notLDL, or notHDL or notTG)
ppmatrix.newhits = cbind(pp.newhits.collapse)[order(ebprior.glhits.collapse,decreasing=TRUE),]
pp.newhits.classmatrix = rbind(colSums(ppmatrix.newhits[allassoc,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,4]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,3]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,2]==0,]),colSums(ppmatrix.newhits[allbut1assoc & modelmatrix[,1]==0,]))
bestclass= apply(pp.newhits.classmatrix, 2,which.max)
cbind(as.character(newhits$snp),bestclass,apply(pp.newhits.classmatrix,2,max))

#90 in new results
#[1,] "rs6691100"  "1"       "0.982087967129813"
#[2,] "rs12133576" "1"       "0.960902710069874"
#[3,] "rs11556341" "1"       "0.975815483615284"
#[4,] "rs9425592"  "1"       "0.999999710910014"
#[5,] "rs6541229"  "1"       "0.983894432409467"
#[6,] "rs4850047"  "1"       "0.722128287338936"
#[7,] "rs1473886"  "1"       "0.990937467103809"
#[8,] "rs10496123" "1"       "0.972744249814961"
#[9,] "rs10175204" "2"       "0.681290587759992"
#[10,] "rs4988235"  "1"       "0.99999999993547"
#[11,] "rs17715343" "1"       "0.993499191599131"
#[12,] "rs12473460" "1"       "0.996981235062487"
#[13,] "rs1174604"  "1"       "0.999982803534828"
#[14,] "rs6739688"  "1"       "0.999337385040667"
#[15,] "rs9856765"  "1"       "0.913163450558902"
#[16,] "rs2062432"  "1"       "0.811776367061361"
#[17,] "rs4683438"  "1"       "0.943562545902491"
#[18,] "rs900399"   "1"       "0.999260679970406"
#[19,] "rs11720640" "1"       "0.995908685169762"
#[20,] "rs4234589"  "1"       "0.995954732610249"
#[21,] "rs176813"   "1"       "0.993570268477009"
#[22,] "rs4859889"  "1"       "0.895208988753498"
#[23,] "rs974801"   "1"       "0.992412211599868"
#[24,] "rs1425486"  "1"       "0.987186201866784"
#[25,] "rs4976033"  "1"       "0.925285199217479"
#[26,] "rs72812840" "1"       "0.97345588596407"
#[27,] "rs13219354" "1"       "0.682134227285844"
#[28,] "rs13205911" "1"       "0.940490307704971"
#[29,] "rs9257809"  "1"       "0.494172913524211"
#[30,] "rs2247056"  "1"       "0.65577762841494"
#[31,] "rs4714556"  "1"       "0.981315744285572"
#[32,] "rs2274517"  "1"       "0.983584625063006"
#[33,] "rs2239620"  "1"       "0.950062864266521"
#[34,] "rs884366"   "1"       "0.991364755704908"
#[35,] "rs4897361"  "1"       "0.479297159709638"
#[36,] "rs2327277"  "1"       "0.911865138643902"
#[37,] "rs3020339"  "1"       "0.974949164028948"
#[38,] "rs2215169"  "1"       "0.9769687253358"
#[39,] "rs6968554"  "1"       "0.915715820357692"
#[40,] "rs4728614"  "1"       "0.904900351218904"
#[41,] "rs17286838" "1"       "0.527189258267025"
#[42,] "rs2686187"  "1"       "0.721922479516088"
#[43,] "rs4871137"  "1"       "0.967199412634128"
#[44,] "rs10087900" "1"       "0.671638726415354"
#[45,] "rs7033354"  "1"       "0.484313735619213"
#[46,] "rs10757056" "1"       "0.989519623540754"
#[47,] "rs17712651" "1"       "0.993889009473228"
#[48,] "rs10733608" "1"       "0.984851154417819"
#[49,] "rs2275774"  "1"       "0.820196348088633"
#[50,] "rs7090871"  "1"       "0.943341200529069"
#[51,] "rs2862954"  "1"       "0.965624635909413"
#[52,] "rs7076938"  "1"       "0.948874310716251"
#[53,] "rs10832027" "1"       "0.998453917316959"
#[54,] "rs211137"   "2"       "0.386882164035866"
#[55,] "rs17309825" "1"       "0.873634849299708"
#[56,] "rs2121703"  "1"       "0.5496043098202"
#[57,] "rs10501321" "4"       "0.60587121258703"
#[58,] "rs11040329" "1"       "0.999928548257938"
#[59,] "rs2201637"  "1"       "0.971873055173386"
#[60,] "rs11229606" "1"       "0.677351778930427"
#[61,] "rs12787728" "1"       "0.999999987102324"
#[62,] "rs2845885"  "1"       "0.993745454704047"
#[63,] "rs661171"   "1"       "0.703738464101595"
#[64,] "rs4149056"  "1"       "0.954679876985734"
#[65,] "rs2278093"  "1"       "0.975920404003111"
#[66,] "rs10861661" "1"       "0.973371017791773"
#[67,] "rs17630235" "1"       "0.973881838924708"
#[68,] "rs895954"   "1"       "0.98315237292604"
#[69,] "rs2454722"  "1"       "0.976836457054309"
#[70,] "rs12423664" "1"       "0.706225099965516"
#[71,] "rs4393438"  "1"       "0.844408236879002"
#[72,] "rs6573939"  "1"       "0.942433545183986"
#[73,] "rs13379043" "1"       "0.989944653747731"
#[74,] "rs721772"   "1"       "0.630212987666067"
#[75,] "rs580469"   "1"       "0.990160484477328"
#[76,] "rs8042320"  "1"       "0.890936535777157"
#[77,] "rs7498491"  "1"       "0.787691838444672"
#[78,] "rs8069974"  "1"       "0.952625227314757"
#[79,] "rs4791641"  "1"       "0.970443594741761"
#[80,] "rs12602912" "1"       "0.98870595220558"
#[81,] "rs11660101" "1"       "0.99549848080848"
#[82,] "rs4808802"  "1"       "0.977413008888971"
#[83,] "rs1688030"  "1"       "0.939087776627577"
#[84,] "rs400058"   "1"       "0.758860335920724"
#[85,] "rs7255743"  "1"       "0.945436945829168"
#[86,] "rs10408163" "1"       "0.957850401179916"
#[87,] "rs2974224"  "2"       "0.254603623445427"
#[88,] "rs1132990"  "1"       "0.983839574259034"
#[89,] "rs6059932"  "1"       "0.967606676662439"
#[90,] "rs6066141"  "1"       "0.99861768877061"






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



pdf("plots.bychr.pdf")
for(i in 1:22){
 if(sum(sub$chr==i)>0){
plot(10^{-6}*gl$pos[gl$chr==i],gl$lbfav[gl$chr==i],ylim=c(0,100),main=paste("Chromosome ", i),xlab="position (Mb)",ylab="log10(BFav)",col=2-(gl$nmin[gl$chr==i]>50000))
abline(v=10^{-6}*gl$pos[gl$chr==i & gl$annot==1],col=2)
abline(v=10^{-6}*sub$pos[sub$annot==0 & sub$chr==i & sub$lbfav>4.07],col=3)
}
}
dev.off()

#try to plot all hits to see relationships?
tophits =  sub[sub$lbfav>4.07 & sub$nmin>50000,]
betahat = tophits[,c(8,11,14,17)]
betahat.scaled = betahat/apply(betahat,1,sd)
betahat.scaled.pr = prcomp(betahat.scaled,center=F)
pdf("plots.betahats.pdf")
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
pdf("plots.cor.zhla.pdf")
image(cor(t(zhla))^2)
dev.off()

pdf("plots.chr1lbfav.pdf")
plot(gl$pos[gl$chr==1],gl$lbfav[gl$chr==1],ylim=c(0,10),xlim=c(27102620-10^5,27102620+10^5))
dev.off()

write.table(file="newtophits.txt",cbind(sub[sub$lbfav>4.07 & sub$nmin>50000 & sub$annot==0,c(4,2:3,7)],round(sub[sub$lbfav>4.07 & sub$nmin>50000 & sub$annot==0,c(22:25,30)],digits=1)),quote=FALSE,sep= " ", row.names=FALSE)

c(4,2:3,7,22:25,30)

newhits = gl[gl$annot==0 & gl$lbfav>4.07 & gl$nmin>50000,c(4,2:3,30)]
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



#NOTE -- everything below is 100% code I have added
~~~~~~~~~~~~~~

lbf.uni.newhits = lbf.uni(lbf.newhits,lbf$gamma)
lbf.all.newhits = lbf.all(lbf.newhits,lbf$gamma)

lbf.uni.newhits < lbf.all.newhits

newhits.BFall <- newhits[lbf.uni.newhits < lbf.all.newhits,]
newhits.BFuni <- newhits[lbf.uni.newhits > lbf.all.newhits,]

pdf("plots.newhits.BFcompare.pdf")

plot(lbf.uni.newhits, lbf.all.newhits, xlab= "Univariate BFs", ylab="Multivariate BFs", main="Newhits: Univariate vs. Multivariate BFs")
abline(a=0, b=1)
abline(lm(lbf.all.newhits~lbf.uni.newhits), col="RED")

dev.off()


lbf.all.newhits[lbf.uni.newhits < lbf.all.newhits]
 [1]   3.968463  68.990173  11.613760   5.645839   6.094224   4.827436
  [7]  51.505439  31.451394  10.455647 112.019029   7.705638  62.248284
  [13]  18.113024   8.531377  40.611264   8.386052   5.495064  15.196135
  [19]   4.591436   8.477698   5.655390  29.637847   9.394281   5.928839
  [25]  16.090102  17.582049  15.402835   3.849938   7.187856  50.880653
  [31]   5.465450   5.239797   7.993237  46.338930  44.081221   4.046166
  [37]  16.314203  10.583469  25.045259  10.315367 223.354893  17.752708
  [43]  10.412600  14.089767  86.731036  21.800449  24.478120  13.624132
  [49]   3.913438   5.427590   7.086014 197.911143   3.976495  10.055042
  [55]   7.937167  20.102563  42.463351   8.063883   9.615528  17.809583
  [61] 135.822736   7.417306  22.992715   4.674452  12.238796

lbf.all.newhits[lbf.uni.newhits > lbf.all.newhits]
 [1]  3.325144  3.748584  5.160717  5.601728  3.914453  8.574840 27.818895
  [8] 12.248346  9.683589 10.114043  4.064096  3.653187 14.991872  2.574374
  [15]  3.470884  4.308010  5.301452  4.692360  4.303231  3.809404 16.369942
  [22]  6.331422  3.134792  7.248606  5.441453


newhits[lbf.all.newhits > 223,]
snp chr      pos     maf      Z.tg     Z.tc      Z.ldl     Z.hdl
2715 rs2121703  11 45762123 0.09235 -1.828464 1.349427 -0.9415951 -5.595174
lbfav
2715 5.637638



