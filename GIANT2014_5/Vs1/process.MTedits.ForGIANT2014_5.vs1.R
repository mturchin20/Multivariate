#this is an incomplete record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# GIANT 2014_5 data

height=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz",header=T)
BMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz",header=T)
WHRadjBMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz",header=T)
WHR=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz",header=T)
HIPadjBMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz",header=T)
HIP=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz",header=T)
WCadjBMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz",header=T)
WC=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz",header=T)

#~~~
#> dim(height)
#[1] 2550858       9
#> dim(BMI)
#[1] 2554637       9
#> dim(WHRadjBMI)
#[1] 2542431       9
#> dim(WHR)
#[1] 2560781      11
#> dim(HIPadjBMI)
#[1] 2540925      11
#> dim(HIP)
#[1] 2559738       9
#> dim(WCadjBMI)
#[1] 2546073      11
#> dim(WC)
#[1] 2565407       9
#~~~

GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }

height <- cbind(height, apply(as.matrix(height$p), 2, GetZScore))
BMI <- cbind(BMI, apply(as.matrix(BMI$p), 2, GetZScore))
WHRadjBMI <- cbind(WHRadjBMI, apply(as.matrix(WHRadjBMI$p), 2, GetZScore))
WHR <- cbind(WHR, apply(as.matrix(WHR$p), 2, GetZScore))
HIPadjBMI <- cbind(HIPadjBMI, apply(as.matrix(HIPadjBMI$p), 2, GetZScore))
HIP <- cbind(HIP, apply(as.matrix(HIP$p), 2, GetZScore))
WCadjBMI <- cbind(WCadjBMI, apply(as.matrix(WCadjBMI$p), 2, GetZScore))
WC <- cbind(WC, apply(as.matrix(WC$p), 2, GetZScore))


height.colnames <- colnames(height)
height.colnames[10] <- "GC.Zscore"
colnames(height) <- height.colnames
BMI.colnames <- colnames(BMI)
BMI.colnames[1] <- "MarkerName"
BMI.colnames[10] <- "GC.Zscore"
colnames(BMI) <- BMI.colnames
WHRadjBMI.colnames <- colnames(WHRadjBMI)
WHRadjBMI.colnames[10] <- "GC.Zscore"
colnames(WHRadjBMI) <- WHRadjBMI.colnames
WHR.colnames <- colnames(WHR)
WHR.colnames[12] <- "GC.Zscore"
colnames(WHR) <- WHR.colnames
HIPadjBMI.colnames <- colnames(HIPadjBMI)
HIPadjBMI.colnames[12] <- "GC.Zscore"
colnames(HIPadjBMI) <- HIPadjBMI.colnames
HIP.colnames <- colnames(HIP)
HIP.colnames[10] <- "GC.Zscore"
colnames(HIP) <- HIP.colnames
WCadjBMI.colnames <- colnames(WCadjBMI)
WCadjBMI.colnames[12] <- "GC.Zscore"
colnames(WCadjBMI) <- WCadjBMI.colnames
WC.colnames <- colnames(WC)
WC.colnames[10] <- "GC.Zscore"
colnames(WC) <- WC.colnames

#~~~
#> head(height)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU       b     SE       p
#1  rs4747841       A       G                  0.551 -0.0011 0.0029 7.0e-01
#2  rs4749917       T       C                  0.436  0.0011 0.0029 7.0e-01
#3   rs737656       A       G                  0.367 -0.0062 0.0030 4.2e-02
#4   rs737657       A       G                  0.358 -0.0062 0.0030 4.1e-02
#5  rs7086391       T       C                  0.120 -0.0087 0.0038 2.4e-02
#6   rs878177       T       C                  0.300  0.0140 0.0031 8.2e-06
#N        ChrBP GC.Zscore
#1 253213  10_10000134 0.3853205
#2 253213  10_10000264 0.3853205
#3 253116 10_100002728 2.0335201
#4 252156 10_100002879 2.0435300
#5 248425 10_100003552 2.2571292
#6 251271 10_100003804 4.4598945
#> head(BMI)
#MarkerName A1 A2 Freq1.Hapmap       b     se       p      N        ChrBP
#1  rs1000000  G  A       0.6333  0.0001 0.0044 0.98190 231410 12_125456932
#2 rs10000010  T  C       0.5750 -0.0029 0.0030 0.33740 322079   4_21227771
#3 rs10000012  G  C       0.1917 -0.0095 0.0054 0.07853 233933    4_1347324
#4 rs10000013  A  C       0.8333 -0.0095 0.0044 0.03084 233886   4_36901463
#5 rs10000017  C  T       0.7667 -0.0034 0.0046 0.45980 233146   4_84997148
#6 rs10000023  G  T       0.4083  0.0024 0.0038 0.52770 233860   4_95952928
#GC.Zscore
#1 0.02268693
#2 0.95931516
#3 1.75927978
#4 2.15913107
#5 0.73917622
#6 0.63152088
#> head(WHRadjBMI)
#MarkerName Allele1 Allele2 FreqAllele1HapMapCEU     b     se       p      N
#1 rs10011200       C       G               0.5333 0.017 0.0043 0.00011 142475
#2  rs8051831       T       C               0.0917 0.034 0.0089 0.00011 138860
#3 rs17542520       A       G               0.0917 0.028 0.0071 0.00011 142581
#4  rs6954671       G       C               0.4167 0.016 0.0042 0.00011 142637
#5  rs6749617       T       A               0.1083 0.021 0.0053 0.00011 206260
#6 rs16873543       C       T               0.2333 0.019 0.0048 0.00011 142701
#ChrBP GC.Zscore
#1 4_145749152  3.867404
#2  16_4343664  3.867404
#3 15_50622399  3.867404
#4  7_77248249  3.867404
#5  2_43698832  3.867404
#6  6_45685111  3.867404
#> head(WHR)
#MarkerName Chr       Pos Allele1 Allele2 FreqAllele1HapMapCEU     b     se
#1  rs2497311  10  94471975       T       C               0.0833 0.023 0.0059
#2  rs7178130  15  70765255       A       G               0.7833 0.017 0.0044
#3  rs7603236   2  66615915       C       T               0.4500 0.016 0.0042
#4  rs6590684  11 132639204       A       G               0.4583 0.016 0.0042
#5 rs12705980   7 114122079       A       G               0.4833 0.016 0.0042
#6  rs1229758   7 114016375       A       G               0.6000 0.016 0.0042
#p      N        ChrBP GC.Zscore
#1 1e-04 196365  10_94471974  3.890592
#2 1e-04 143826  15_70765254  3.890592
#3 1e-04 143613   2_66615914  3.890592
#4 1e-04 144487 11_132639203  3.890592
#5 1e-04 144572  7_114122078  3.890592
#6 1e-04 144489  7_114016374  3.890592
#> head(HIPadjBMI)
#MarkerName Chr       Pos Allele1 Allele2 FreqAllele1HapMapCEU      b     se
#1 rs17321515   8 126555591       A       G               0.6083 -0.017 0.0043
#2 rs11049399  12  28218413       C       T               0.6500  0.018 0.0047
#3 rs10442655   1 169139023       T       C               0.1525  0.025 0.0065
#4 rs10195055   2  67766101       G       T               0.5417 -0.017 0.0043
#5  rs2476391  13  49516505       C       T               0.0500 -0.037 0.0095
#6  rs6914419   6 130478165       C       T               0.6750 -0.018 0.0047
#p      N       ChrBP GC.Zscore
#1 1e-04 143707 8_126555590  3.890592
#2 1e-04 143788 12_28218412  3.890592
#3 1e-04 141602 1_169139022  3.890592
#4 1e-04 143712  2_67766100  3.890592
#5 1e-04 143810 13_49516504  3.890592
#6 1e-04 132003 6_130478164  3.890592
#> head(HIP)
#MarkerName Allele1 Allele2 FreqAllele1HapMapCEU      b     se     p      N
#1 chr6:43865060       A       T                   NA -0.024 0.0061 1e-04  63525
#2     rs6595438       T       C               0.6000  0.017 0.0044 1e-04 145453
#3    rs17259223       A       C               0.6417  0.018 0.0046 1e-04 145386
#4     rs7641973       A       G               0.3833  0.018 0.0045 1e-04 145366
#5     rs7725897       T       C               0.1750 -0.021 0.0055 1e-04 145447
#6     rs4968100       G       A               0.2000 -0.024 0.0061 1e-04 145380
#ChrBP GC.Zscore
#1        <NA>  3.890592
#2 5_122711894  3.890592
#3 18_19239629  3.890592
#4  3_81993178  3.890592
#5 5_122704269  3.890592
#6   17_547484  3.890592
#> head(WCadjBMI)
#MarkerName Chr      Pos Allele1 Allele2 FreqAllele1HapMapCEU      b     se
#1  rs7638333   3 12781142       T       C               0.3083  0.013 0.0035
#2  rs7947604  11 65198635       G       C               0.6583 -0.017 0.0044
#3 rs12654696   5 64588544       C       T               0.1833  0.020 0.0052
#4  rs2448428   5 55879806       T       C               0.4083  0.013 0.0035
#5  rs2011074   5 88776815       T       C               0.2917  0.018 0.0045
#6   rs679805  13 79645702       T       C               0.5250  0.017 0.0043
#p      N       ChrBP GC.Zscore
#1 1e-04 229089  3_12781141  3.890592
#2 1e-04 151256 11_65198634  3.890592
#3 1e-04 153438  5_64588543  3.890592
#4 1e-04 230984  5_55879805  3.890592
#5 1e-04 152701  5_88776814  3.890592
#6 1e-04 150635 13_79645701  3.890592
#> head(WC)
#MarkerName Allele1 Allele2 FreqAllele1HapMapCEU      b     se     p      N
#1 rs10508465       C       T               0.3583  0.017 0.0044 1e-04 153945
#2  rs2291285       G       A               0.5667  0.017 0.0043 1e-04 153841
#3 rs12967873       G       A               0.2083  0.021 0.0055 1e-04 153802
#4  rs1515098       T       C               0.6417 -0.014 0.0036 1e-04 231812
#5  rs9310067       G       T               0.8000  0.023 0.0060 1e-04 153948
#6  rs2066295       G       A               0.3083 -0.020 0.0051 1e-04 153924
#ChrBP GC.Zscore
#1 10_13725193  3.890592
#2 17_31964474  3.890592
#3 18_38104620  3.890592
#4 2_226782097  3.890592
#5  3_88206787  3.890592
#6  6_26276881  3.890592
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=height$MarkerName,Z.height=height$GC.Zscore,key="rs")
dt2 = data.table(rs=BMI$MarkerName,Z.BMI=BMI$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)

#clear some space 
#rm(height)
#rm(BMI)
rm(dt2)
rm(dt1)

dt1 = data.table(rs=WHRadjBMI$MarkerName,Z.WHRadjBMI=WHRadjBMI$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)

rm(dt1)
#rm(WHRadjBMI)
rm(dt3)

##dt1 = data.table(rs=tg$MarkerName,Z.tg=tg$GC.Zscore,key="rs")
##dt3 = merge(dt1,dt2)
##rm(tg)
##rm(dt1)

#~~~
#> dim(dt2)
#[1] 2505734       4
#~~~

height.maxZ <- max(dt2$Z.height[!is.infinite(dt2$Z.height)])
BMI.maxZ <- max(dt2$Z.BMI[!is.infinite(dt2$Z.BMI)])
WHRadjBMI.maxZ <- max(dt2$Z.WHRadjBMI[!is.infinite(dt2$Z.WHRadjBMI)])
maxZ <- max(c(height.maxZ, BMI.maxZ, WHRadjBMI.maxZ))

#~~~
#> height.maxZ
#[1] 26.80036
#> BMI.maxZ
#[1] 26.33557
#> WHRadjBMI.maxZ
#[1] 12.3743
#> maxZ
#[1] 26.80036
#~~~

#replaceInf <- function(x) { val2 <- x; if (is.infinite(x)) { val2 <- maxZ; } return(val2); }
replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt2$Z.height <- apply(as.matrix(dt2$Z.height), 1, replaceInf)
dt2$Z.BMI <- apply(as.matrix(dt2$Z.BMI), 1, replaceInf)
dt2$Z.WHRadjBMI <- apply(as.matrix(dt2$Z.WHRadjBMI), 1, replaceInf)

#~~~
#> replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }
#>
#> dt2$Z.height <- apply(as.matrix(dt2$Z.height), 1, replaceInf)
#> dt2$Z.BMI <- apply(as.matrix(dt2$Z.BMI), 1, replaceInf)
#> dt2$Z.WHRadjBMI <- apply(as.matrix(dt2$Z.WHRadjBMI), 1, replaceInf)
#>
#~~~

attach(dt2)
nullset = (abs(dt2$Z.BMI)<2) & (abs(dt2$Z.height)<2) & (abs(dt2$Z.WHRadjBMI)<2) #extract null Z values
nullset_lte1 = (abs(dt2$Z.BMI)<1) & (abs(dt2$Z.height)<1) & (abs(dt2$Z.WHRadjBMI)<1) #extract null Z values
nullset_lte3 = (abs(dt2$Z.BMI)<3) & (abs(dt2$Z.height)<3) & (abs(dt2$Z.WHRadjBMI)<3) #extract null Z values
Z = cbind(Z.BMI,Z.height,Z.WHRadjBMI)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

#~~~
#> dim(Z)
#[1] 2505734       3
#> dim(Znull)
#[1] 1848720       3
#> RSS0
#Z.BMI    Z.height Z.WHRadjBMI
#Z.BMI       1.000000000 0.004024732 0.004761944
#Z.height    0.004024732 1.000000000 0.010661975
#Z.WHRadjBMI 0.004761944 0.010661975 1.000000000
#~~~

mvstat =  rowSums(Z * (Z %*% RSS0inv)) # comptues Z RSS0^-1 Z'
dt2$mvstat = mvstat
statchi = -log10(exp(1))*pchisq(mvstat,df=4,log.p=TRUE,lower.tail=FALSE)
dt2$mvp = statchi

maxZ2 = apply(Z^2,1,max)
max.unip = -log10(exp(1))*pchisq(maxZ2,df=1,log.p=TRUE, lower.tail=FALSE)
dt2$unip = max.unip

dtsignif = dt2[dt2$mvp>-log10(5e-8) | dt2$unip>-log10(5e-8),]
dtlesssignif = dt2[dt2$mvp>-log10(1e-6) | dt2$unip>-log10(1e-6),]

write.table(file="GIANT2014_5.Orig3.dtsignif.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.Orig3.dtsignif.rs.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.Orig3.dtlesssignif.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.Orig3.dtlesssignif.rs.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GIANT2014_5.Orig3.dtlesssignif.rs.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GIANT2014_5.Orig3.RSS0.vs1.txt",RSS0,row.names=F,sep=",",quote=F)

#~~~
#> dim(dtsignif)
#[1] 29600     7
#> dim(dtlesssignif)
#[1] 42022     7
#> max(dt2$mvstat)
#[1] 732.2285
#> max(dt2$mvp)
#[1] 156.4366
#> max(dt2$unip)
#[1] 157.4949
#> quantile(dt2$mvp)
#0%          25%          50%          75%         100%
#1.339890e-09 8.452235e-02 2.642825e-01 6.700757e-01 1.564366e+02
#> quantile(dt2$unip)
#0%          25%          50%          75%         100%
#4.364805e-03 5.086383e-01 8.317973e-01 1.377786e+00 1.574949e+02
#~~~

rm(dt2)


#AllPhenos
#create a merged table of Z scores
dt1 = data.table(rs=height$MarkerName,Z.height=height$GC.Zscore,key="rs")
dt2 = data.table(rs=BMI$MarkerName,Z.BMI=BMI$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=WHRadjBMI$MarkerName,Z.WHRadjBMI=WHRadjBMI$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=WHR$MarkerName,Z.WHR=WHR$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)
dt1 = data.table(rs=HIPadjBMI$MarkerName,Z.HIPadjBMI=HIPadjBMI$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=HIP$MarkerName,Z.HIP=HIP$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)
dt1 = data.table(rs=WCadjBMI$MarkerName,Z.WCadjBMI=WCadjBMI$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=WC$MarkerName,Z.WC=WC$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)

#~~~
#> dim(dt3)
#[1] 2502096       9
#> head(dt3)
#rs   Z.height      Z.BMI Z.WHRadjBMI     Z.WHR Z.HIPadjBMI
#1:       rs10 0.35845879 0.57144210  0.22754498 0.2533471  0.89647336
#2:  rs1000000 2.14441062 0.02268693  1.12639113 0.6433454  0.26631061
#3: rs10000010 0.11303854 0.95931516  1.20035886 0.8064212  0.06270678
#4: rs10000012 3.01145376 1.75927978  1.51410189 2.0435300  0.99445788
#5: rs10000013 0.02506891 2.15913107  0.53883603 0.8064212  0.65883769
#6: rs10000017 1.31057911 0.73917622  0.02506891 0.2404260  1.12639113
#Z.HIP Z.WCadjBMI       Z.WC
#1: 0.07526986  0.8416212 0.03760829
#2: 0.26631061  1.1749868 0.37185609
#3: 0.85961736  1.2815516 0.39885507
#4: 1.01522203  0.7388468 1.77438191
#5: 1.90331082  0.2793190 1.97736843
#6: 0.87789630  0.7388468 0.77219321
#~~~

height.maxZ <- max(dt3$Z.height[!is.infinite(dt3$Z.height)])
BMI.maxZ <- max(dt3$Z.BMI[!is.infinite(dt3$Z.BMI)])
WHRadjBMI.maxZ <- max(dt3$Z.WHRadjBMI[!is.infinite(dt3$Z.WHRadjBMI)])
WHR.maxZ <- max(dt3$Z.WHR[!is.infinite(dt3$Z.WHR)])
HIPadjBMI.maxZ <- max(dt3$Z.HIPadjBMI[!is.infinite(dt3$Z.HIPadjBMI)])
HIP.maxZ <- max(dt3$Z.HIP[!is.infinite(dt3$Z.HIP)])
WCadjBMI.maxZ <- max(dt3$Z.WCadjBMI[!is.infinite(dt3$Z.WCadjBMI)])
WC.maxZ <- max(dt3$Z.WC[!is.infinite(dt3$Z.WC)])
maxZ <- max(c(height.maxZ, BMI.maxZ, WHRadjBMI.maxZ, WHR.maxZ, HIPadjBMI.maxZ, HIP.maxZ, WCadjBMI.maxZ, WC.maxZ))

#~~~
#> height.maxZ
#[1] 26.80036
#> BMI.maxZ
#[1] 26.33557
#> WHRadjBMI.maxZ
#[1] 12.3743
#> WHR.maxZ
#[1] 12.99536
#> HIPadjBMI.maxZ
#[1] 13.72197
#> HIP.maxZ
#[1] 19.69673
#> WCadjBMI.maxZ
#[1] 11.12296
#> WC.maxZ
#[1] 21.35245
#> maxZ
#[1] 26.80036
#~~~

dt3$Z.height <- apply(as.matrix(dt3$Z.height), 1, replaceInf)
dt3$Z.BMI <- apply(as.matrix(dt3$Z.BMI), 1, replaceInf)
dt3$Z.WHRadjBMI <- apply(as.matrix(dt3$Z.WHRadjBMI), 1, replaceInf)
dt3$Z.WHR <- apply(as.matrix(dt3$Z.WHR), 1, replaceInf)
dt3$Z.HIPadjBMI <- apply(as.matrix(dt3$Z.HIPadjBMI), 1, replaceInf)
dt3$Z.HIP <- apply(as.matrix(dt3$Z.HIP), 1, replaceInf)
dt3$Z.WCadjBMI <- apply(as.matrix(dt3$Z.WCadjBMI), 1, replaceInf)
dt3$Z.WC <- apply(as.matrix(dt3$Z.WC), 1, replaceInf)

#~~~
#> dt3$Z.height <- apply(as.matrix(dt3$Z.height), 1, replaceInf)
#> dt3$Z.WHRadjBMI <- apply(as.matrix(dt3$Z.WHRadjBMI), 1, replaceInf)
#> dt3$Z.WHR <- apply(as.matrix(dt3$Z.WHR), 1, replaceInf)
#> dt3$Z.HIPadjBMI <- apply(as.matrix(dt3$Z.HIPadjBMI), 1, replaceInf)
#> dt3$Z.HIP <- apply(as.matrix(dt3$Z.HIP), 1, replaceInf)
#> dt3$Z.WCadjBMI <- apply(as.matrix(dt3$Z.WCadjBMI), 1, replaceInf)
#> dt3$Z.WC <- apply(as.matrix(dt3$Z.WC), 1, replaceInf)
#~~~

attach(dt3)
nullset = (abs(dt3$Z.BMI)<2) & (abs(dt3$Z.height)<2) & (abs(dt3$Z.WHRadjBMI)<2) & (abs(dt3$Z.WHR)<2) & (abs(dt3$Z.HIPadjBMI)<2) & (abs(dt3$Z.HIP)<2) & (abs(dt3$Z.WCadjBMI)<2) & (abs(dt3$Z.WC)<2) #extract null Z values
Z = cbind(Z.BMI,Z.height,Z.WHRadjBMI,Z.WHR,Z.HIPadjBMI,Z.HIP,Z.WCadjBMI,Z.WC)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

#~~~
> dim(Z)
[1] 2502096       8
> dim(Znull)
[1] 1662836       8
> RSS0
Z.BMI     Z.height  Z.WHRadjBMI        Z.WHR Z.HIPadjBMI
Z.BMI        1.000000000  0.002750141 -0.021741922  0.070091252 -0.01734839
Z.height     0.002750141  1.000000000  0.001228069 -0.002312525  0.07871526
Z.WHRadjBMI -0.021741922  0.001228069  1.000000000  0.624754356  0.05189382
Z.WHR        0.070091252 -0.002312525  0.624754356  1.000000000  0.04507764
Z.HIPadjBMI -0.017348393  0.078715263  0.051893819  0.045077636  1.00000000
Z.HIP        0.293468426  0.004696159 -0.002337970 -0.023693018  0.17681032
Z.WCadjBMI  -0.020066493  0.038778768  0.358187399  0.211264081  0.05175004
Z.WC         0.359179900 -0.003433090  0.030403583  0.409559364 -0.02311725
Z.HIP   Z.WCadjBMI        Z.WC
Z.BMI        0.293468426 -0.020066493  0.35917990
Z.height     0.004696159  0.038778768 -0.00343309
Z.WHRadjBMI -0.002337970  0.358187399  0.03040358
Z.WHR       -0.023693018  0.211264081  0.40955936
Z.HIPadjBMI  0.176810315  0.051750042 -0.02311725
Z.HIP        1.000000000 -0.009618877  0.48498439
Z.WCadjBMI  -0.009618877  1.000000000  0.13920749
Z.WC         0.484984395  0.139207490  1.00000000
#~~~

mvstat =  rowSums(Z * (Z %*% RSS0inv)) # comptues Z RSS0^-1 Z'
dt3$mvstat = mvstat
statchi = -log10(exp(1))*pchisq(mvstat,df=4,log.p=TRUE,lower.tail=FALSE)
dt3$mvp = statchi

maxZ2 = apply(Z^2,1,max)
max.unip = -log10(exp(1))*pchisq(maxZ2,df=1,log.p=TRUE, lower.tail=FALSE)
dt3$unip = max.unip

dtsignif = dt3[dt3$mvp>-log10(5e-8) | dt3$unip>-log10(5e-8),]
dtlesssignif = dt3[dt3$mvp>-log10(1e-6) | dt3$unip>-log10(1e-6),]

write.table(file="GIANT2014_5.AllPheno.dtsignif.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.AllPheno.dtsignif.rs.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.AllPheno.dtlesssignif.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.AllPheno.dtlesssignif.rs.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GIANT2014_5.AllPheno.dtlesssignif.rs.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GIANT2014_5.AllPheno.RSS0.vs1.txt",RSS0,row.names=F,sep=",",quote=F)

#~~~
#> dim(dtsignif)
#[1] 38977    12
#> dim(dtlesssignif)
#[1] 56322    12
#> max(dt3$mvstat)
#[1] 1055.865
#> max(dt3$mvp)
#[1] 226.5548
#> max(dt3$unip)
#[1] 157.4949
#> quantile(dt3$mvp)
#0%          25%          50%          75%         100%
#5.980426e-05 2.781496e-01 6.297995e-01 1.289772e+00 2.265548e+02
#> quantile(dt3$unip)
#0%          25%          50%          75%         100%
#0.03621217   0.67778071   1.02227639   1.56863624 157.49485002
#~~~

rm(dt3)


#No adj phenos









