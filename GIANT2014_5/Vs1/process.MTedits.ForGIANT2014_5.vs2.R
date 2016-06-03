#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# GIANT2014_5 data

/mnt/gluster/data/external_public_supp/GIANT2014_5/
-rw-rw-r-- 1 mturchin20 mturchin20 50135820 Aug 29 16:01 GIANT_2015_BMI.ForAnalysis.wUCSCGB_snp130.vs3.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 45886394 Aug 29 16:12 GIANT_2014_Height_publicrelease_HapMapCeuFreq.ForAnalysis.wUCSCGB_snp130.vs3.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 48177538 Aug 29 16:20 GIANT_2015_WHRadjBMI_COMBINED_EUR.ForAnalysis.wUCSCGB_snp130.vs3.HeaderFix.txt.gz

Height=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2014_Height_publicrelease_HapMapCeuFreq.ForAnalysis.wUCSCGB_snp130.vs3.txt.gz", header=T)
BMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_BMI.ForAnalysis.wUCSCGB_snp130.vs3.txt.gz", header=T)
WHRadjBMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.ForAnalysis.wUCSCGB_snp130.vs3.HeaderFix.txt.gz", header=T)

dim(Height)
dim(BMI)
dim(WHRadjBMI)

#~~~
#> dim(Height)
#[1] 2550858       9
#> dim(BMI)
#[1] 2554637       9
#> dim(WHRadjBMI)
#[1] 2542431       9
#~~~

GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }

Height <- cbind(Height, apply(as.matrix(Height$p), 2, GetZScore))
BMI <- cbind(BMI, apply(as.matrix(BMI$p), 2, GetZScore))
WHRadjBMI <- cbind(WHRadjBMI, apply(as.matrix(WHRadjBMI$p), 2, GetZScore))

Height.colnames <- colnames(Height)
Height.colnames[10] <- "GC.Zscore"
colnames(Height) <- Height.colnames
BMI.colnames <- colnames(BMI)
BMI.colnames[1] <- "MarkerName"
BMI.colnames[10] <- "GC.Zscore"
colnames(BMI) <- BMI.colnames
WHRadjBMI.colnames <- colnames(WHRadjBMI)
WHRadjBMI.colnames[10] <- "GC.Zscore"
colnames(WHRadjBMI) <- WHRadjBMI.colnames

head(Height)
head(BMI)
head(WHRadjBMI)

#~~~
#> head(Height)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU       b     SE       p
#1  rs4747841       A       G                  0.551 -0.0011 0.0029 7.0e-01
#2  rs4749917       T       C                  0.436  0.0011 0.0029 7.0e-01
#3   rs737656       A       G                  0.367 -0.0062 0.0030 4.2e-02
#4   rs737657       A       G                  0.358 -0.0062 0.0030 4.1e-02
#5  rs7086391       T       C                  0.120 -0.0087 0.0038 2.4e-02
#6   rs878177       T       C                  0.300  0.0140 0.0031 8.2e-06
#N        ChrBP GC.Zscore
#1 253213  10_10000135 0.3853205
#2 253213  10_10000265 0.3853205
#3 253116 10_100002729 2.0335201
#4 252156 10_100002880 2.0435300
#5 248425 10_100003553 2.2571292
#6 251271 10_100003805 4.4598945
#> head(BMI)
#MarkerName A1 A2 Freq1.Hapmap       b     se       p      N        ChrBP
#1  rs1000000  G  A       0.6333  0.0001 0.0044 0.98190 231410 12_125456933
#2 rs10000010  T  C       0.5750 -0.0029 0.0030 0.33740 322079   4_21227772
#3 rs10000012  G  C       0.1917 -0.0095 0.0054 0.07853 233933    4_1347325
#4 rs10000013  A  C       0.8333 -0.0095 0.0044 0.03084 233886   4_36901464
#5 rs10000017  C  T       0.7667 -0.0034 0.0046 0.45980 233146   4_84997149
#6 rs10000023  G  T       0.4083  0.0024 0.0038 0.52770 233860   4_95952929
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
#1 4_145749153  3.867404
#2  16_4343665  3.867404
#3 15_50622400  3.867404
#4  7_77248250  3.867404
#5  2_43698833  3.867404
#6  6_45685112  3.867404
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=Height$MarkerName,Z.Height=Height$GC.Zscore,key="rs")
dt2 = data.table(rs=BMI$MarkerName,Z.BMI=BMI$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=WHRadjBMI$MarkerName,Z.WHRadjBMI=WHRadjBMI$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)

dim(dt2)
head(dt2)

#~~~
#> dim(dt2)
#[1] 2505734       4
#> head(dt2)
#rs   Z.Height      Z.BMI Z.WHRadjBMI
#1:       rs10 0.35845879 0.57144210  0.22754498
#2:  rs1000000 2.14441062 0.02268693  1.12639113
#3: rs10000010 0.11303854 0.95931516  1.20035886
#4: rs10000012 3.01145376 1.75927978  1.51410189
#5: rs10000013 0.02506891 2.15913107  0.53883603
#6: rs10000017 1.31057911 0.73917622  0.02506891
#~~~

Height.maxZ <- max(dt2$Z.Height[!is.infinite(dt2$Z.Height) & !is.na(dt2$Z.Height)])
BMI.maxZ <- max(dt2$Z.BMI[!is.infinite(dt2$Z.BMI) & !is.na(dt2$Z.BMI)])
WHRadjBMI.maxZ <- max(dt2$Z.WHRadjBMI[!is.infinite(dt2$Z.WHRadjBMI) & !is.na(dt2$Z.WHRadjBMI)])
maxZ <- max(c(Height.maxZ, BMI.maxZ, WHRadjBMI.maxZ))

Height.maxZ
BMI.maxZ
WHRadjBMI.maxZ
maxZ

#~~~
#> Height.maxZ
#[1] 26.80036
#> BMI.maxZ
#[1] 26.33557
#> WHRadjBMI.maxZ
#[1] 12.3743
#> maxZ
#[1] 26.80036
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt2$Z.Height <- apply(as.matrix(dt2$Z.Height), 1, replaceInf)
dt2$Z.BMI <- apply(as.matrix(dt2$Z.BMI), 1, replaceInf)
dt2$Z.WHRadjBMI <- apply(as.matrix(dt2$Z.WHRadjBMI), 1, replaceInf)

#~~~
#> dt2$Z.Height <- apply(as.matrix(dt2$Z.Height), 1, replaceInf)
#> dt2$Z.BMI <- apply(as.matrix(dt2$Z.BMI), 1, replaceInf)
#> dt2$Z.WHRadjBMI <- apply(as.matrix(dt2$Z.WHRadjBMI), 1, replaceInf)
#~~~

attach(dt2)
nullset = (abs(dt2$Z.BMI)<2) & (abs(dt2$Z.Height)<2) & (abs(dt2$Z.WHRadjBMI)<2) #extract null Z values
Z = cbind(Z.BMI,Z.Height,Z.WHRadjBMI)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
#> dim(Z)
#[1] 2505734       3
#> dim(Znull)
#[1] 1848720       3
#> RSS0
#Z.BMI    Z.Height Z.WHRadjBMI
#Z.BMI       1.000000000 0.004024732 0.004761944
#Z.Height    0.004024732 1.000000000 0.010661975
#Z.WHRadjBMI 0.004761944 0.010661975 1.000000000
#~~~

mvstat =  rowSums(Z * (Z %*% RSS0inv)) # comptues Z RSS0^-1 Z'
dt2$mvstat = mvstat
statchi = -log10(exp(1))*pchisq(mvstat,df=3,log.p=TRUE,lower.tail=FALSE)
dt2$mvp = statchi

maxZ2 = apply(Z^2,1,max)
max.unip = -log10(exp(1))*pchisq(maxZ2,df=1,log.p=TRUE, lower.tail=FALSE)
dt2$unip = max.unip

dtsignif = dt2[dt2$mvp>-log10(5e-8) | dt2$unip>-log10(5e-8),]
dtlesssignif = dt2[dt2$mvp>-log10(1e-6) | dt2$unip>-log10(1e-6),]

dim(dtsignif)
dim(dtlesssignif)
max(dt2$mvstat)
max(dt2$mvp)
max(dt2$unip)
quantile(dt2$mvp)
quantile(dt2$unip)

#~~~
#> dim(dtsignif)
#[1] 30362     7
#> dim(dtlesssignif)
#[1] 43086     7
#> max(dt2$mvstat)
#[1] 732.2285
#> max(dt2$mvp)
#[1] 157.6665
#> max(dt2$unip)
#[1] 157.4949
#> quantile(dt2$mvp)
#0%          25%          50%          75%         100%
#2.274474e-07 1.689416e-01 4.212353e-01 9.164765e-01 1.576665e+02
#> quantile(dt2$unip)
#0%          25%          50%          75%         100%
#4.364805e-03 5.086383e-01 8.317973e-01 1.377786e+00 1.574949e+02
#~~~

write.table(file="GIANT2014_5.vs2.dtsignif.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.vs2.dtsignif.rs.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.vs2.dtlesssignif.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2014_5.vs2.dtlesssignif.rs.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GIANT2014_5.vs2.dtlesssignif.rs.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GIANT2014_5.vs2.RSS0.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


