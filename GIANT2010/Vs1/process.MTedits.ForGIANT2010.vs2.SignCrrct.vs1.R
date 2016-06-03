#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# GIANT2010 data

/mnt/gluster/data/external_public_supp/GIANT2010/
-rw-rw-r-- 1 mturchin20 mturchin20 33480881 Aug 29 03:18 GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGB_snp130.vs3.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 36138426 Aug 29 03:23 GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGB_snp130.vs3.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 38278896 Aug 29 03:36 GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGB_snp130.vs3.txt.gz

Height=read.table("/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz", header=T)
BMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.MatchedToHeight.txt.gz", header=T)
WHRadjBMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.MatchedToHeight.txt.gz", header=T)

dim(Height)
dim(BMI)
dim(WHRadjBMI)

#~~~
> dim(Height)
[1] 2469635       7
> dim(BMI)
[1] 2471516       7
> dim(WHRadjBMI)
[1] 2483325       7
#~~~

#GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }
GetZScoreAndSign <- function(x) { val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE); if (x < 0) { val1 <- val1 * -1 }; return(val1) }

#HDL <- cbind(HDL, apply(as.matrix(HDL$P.value), 1, GetZScoreAndSign))
#LDL <- cbind(LDL, apply(as.matrix(LDL$P.value), 1, GetZScoreAndSign))
#TC <- cbind(TC, apply(as.matrix(TC$P.value), 1, GetZScoreAndSign))
#TG <- cbind(TG, apply(as.matrix(TG$P.value), 1, GetZScoreAndSign))

Height <- cbind(Height, apply(as.matrix(Height$p), 1, GetZScoreAndSign))
BMI <- cbind(BMI, apply(as.matrix(BMI$p), 1, GetZScoreAndSign))
WHRadjBMI <- cbind(WHRadjBMI, apply(as.matrix(WHRadjBMI$p), 1, GetZScoreAndSign))

Height.colnames <- colnames(Height)
Height.colnames[8] <- "GC.Zscore"
colnames(Height) <- Height.colnames
BMI.colnames <- colnames(BMI)
BMI.colnames[8] <- "GC.Zscore"
colnames(BMI) <- BMI.colnames
WHRadjBMI.colnames <- colnames(WHRadjBMI)
WHRadjBMI.colnames[8] <- "GC.Zscore"
colnames(WHRadjBMI) <- WHRadjBMI.colnames

head(Height)
head(BMI)
head(WHRadjBMI)

#~~~
#> head(Height)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU      p      N        ChrBP
#1       rs10       a       c                 0.0333 0.8826  78380   7_92221824
#2  rs1000000       a       g                 0.3667 0.1858 133822 12_125456933
#3 rs10000010       t       c                  0.575 0.8947 132858   4_21227772
#4 rs10000012       c       g                 0.8083 0.1312 133785    4_1347325
#5 rs10000013       a       c                 0.8333 0.6280 133843   4_36901464
#6 rs10000017       t       c                 0.2333 0.3073 133174   4_84997149
#GC.Zscore
#1 0.1476741
#2 1.3231064
#3 0.1323594
#4 1.5093867
#5 0.4845438
#6 1.0209038
#> head(BMI)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU      p      N        ChrBP
#1       rs10       a       c                 0.0333 0.7080  80566   7_92221824
#2  rs1000000       g       a                 0.6333 0.5060 123865 12_125456933
#3 rs10000010       c       t                  0.425 0.7360 123827   4_21227772
#4 rs10000012       c       g                 0.8083 0.0420 123809    4_1347325
#5 rs10000013       c       a                 0.1667 0.0689 123863   4_36901464
#6 rs10000017       t       c                 0.2333 0.4570 123262   4_84997149
#GC.Zscore
#1 0.3745435
#2 0.6650789
#3 0.3371551
#4 2.0335201
#5 1.8190749
#6 0.7437958
#> head(WHRadjBMI)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU      p     N        ChrBP
#1       rs10       c       a                 0.9667 0.4200 57031   7_92221824
#2  rs1000000       g       a                 0.6333 0.5500 77168 12_125456933
#3 rs10000010       t       c                  0.575 0.0029 77152   4_21227772
#4 rs10000012       g       c                 0.1917 0.9900 77117    4_1347325
#5 rs10000013       a       c                 0.8333 0.8900 77167   4_36901464
#6 rs10000017       c       t                 0.7667 0.6400 77166   4_84997149
#GC.Zscore
#1 0.80642125
#2 0.59776013
#3 2.97814368
#4 0.01253347
#5 0.13830421
#6 0.46769880
> head(Height)
  MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU      p      N        ChrBP
1       rs10       a       c                 0.0333 0.8826  78380   7_92221823
2  rs1000000       a       g                 0.3667 0.1858 133822 12_125456932
3 rs10000010       t       c                  0.575 0.8947 132858   4_21227771
4 rs10000012       c       g                 0.8083 0.1312 133785    4_1347324
5 rs10000013       a       c                 0.8333 0.6280 133843   4_36901463
6 rs10000017       t       c                 0.2333 0.3073 133174   4_84997148
  GC.Zscore
1 0.1476741
2 1.3231064
3 0.1323594
4 1.5093867
5 0.4845438
6 1.0209038
> head(BMI)
  MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU       p      N        ChrBP
1  rs1000000       g       a                 0.6333 -0.5060 123865 12_125456932
2 rs10000010       c       t                  0.425 -0.7360 123827   4_21227771
3 rs10000012       c       g                 0.8083  0.0420 123809    4_1347324
4 rs10000013       c       a                 0.1667 -0.0689 123863   4_36901463
5 rs10000017       t       c                 0.2333  0.4570 123262   4_84997148
6 rs10000023       t       g                 0.5917  0.9390 123756   4_95952928
    GC.Zscore
1 -0.66507895
2 -0.33715508
3  2.03352015
4 -1.81907493
5  0.74379584
6  0.07652679
> head(WHRadjBMI)
  MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU       p     N        ChrBP
1  rs1000000       g       a                 0.6333 -0.5500 77168 12_125456932
2 rs10000010       t       c                  0.575  0.0029 77152   4_21227771
3 rs10000012       g       c                 0.1917 -0.9900 77117    4_1347324
4 rs10000013       a       c                 0.8333  0.8900 77167   4_36901463
5 rs10000017       c       t                 0.7667 -0.6400 77166   4_84997148
6 rs10000023       g       t                 0.4083 -0.0410 77140   4_95952928
    GC.Zscore
1 -0.59776013
2  2.97814368
3 -0.01253347
4  0.13830421
5 -0.46769880
6 -2.04353001
> 
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
#[1] 2436719       4
#> head(dt2)
#rs  Z.Height     Z.BMI Z.WHRadjBMI
#1:       rs10 0.1476741 0.3745435  0.80642125
#2:  rs1000000 1.3231064 0.6650789  0.59776013
#3: rs10000010 0.1323594 0.3371551  2.97814368
#4: rs10000012 1.5093867 2.0335201  0.01253347
#5: rs10000013 0.4845438 1.8190749  0.13830421
#6: rs10000017 1.0209038 0.7437958  0.46769880
> dim(dt2)
[1] 2233314       4
> head(dt2)
           rs  Z.Height       Z.BMI Z.WHRadjBMI
1:  rs1000000 1.3231064 -0.66507895 -0.59776013
2: rs10000010 0.1323594 -0.33715508  2.97814368
3: rs10000012 1.5093867  2.03352015 -0.01253347
4: rs10000013 0.4845438 -1.81907493  0.13830421
5: rs10000017 1.0209038  0.74379584 -0.46769880
6: rs10000023 0.9268585  0.07652679 -2.04353001
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
> Height.maxZ
[1] 15.18464
> BMI.maxZ
[1] 16.67329
> WHRadjBMI.maxZ
[1] 7.773081
> maxZ
[1] 16.67329
> Height.maxZ
[1] 14.82918
> BMI.maxZ
[1] 9.751367
> WHRadjBMI.maxZ
[1] 7.773081
> maxZ
[1] 14.82918

#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt2$Z.Height <- apply(as.matrix(dt2$Z.Height), 1, replaceInf)
dt2$Z.BMI <- apply(as.matrix(dt2$Z.BMI), 1, replaceInf)
dt2$Z.WHRadjBMI <- apply(as.matrix(dt2$Z.WHRadjBMI), 1, replaceInf)

#~~~
#> dt2$Z.Height <- apply(as.matrix(dt2$Z.Height), 1, replaceInf)
#> dt2$Z.BMI <- apply(as.matrix(dt2$Z.BMI), 1, replaceInf)
#> dt2$Z.WHRadjBMI <- apply(as.matrix(dt2$Z.WHRadjBMI), 1, replaceInf)
> dt2$Z.Height <- apply(as.matrix(dt2$Z.Height), 1, replaceInf)
> dt2$Z.BMI <- apply(as.matrix(dt2$Z.BMI), 1, replaceInf)
> dt2$Z.WHRadjBMI <- apply(as.matrix(dt2$Z.WHRadjBMI), 1, replaceInf)

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
#[1] 2436719       3
#> dim(Znull)
#[1] 2041036       3
#> RSS0
#Z.BMI    Z.Height  Z.WHRadjBMI
#Z.BMI       1.0000000000 0.007383229 0.0003445489
#Z.Height    0.0073832291 1.000000000 0.0093910411
#Z.WHRadjBMI 0.0003445489 0.009391041 1.0000000000
> dim(Z)
[1] 2233314       3
> dim(Znull)
[1] 1873246       3
> RSS0
                  Z.BMI     Z.Height  Z.WHRadjBMI
Z.BMI        1.00000000 -0.028708246  0.014570846
Z.Height    -0.02870825  1.000000000 -0.003307104
Z.WHRadjBMI  0.01457085 -0.003307104  1.000000000

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
#[1] 6077    7
#> dim(dtlesssignif)
#[1] 9118    7
#> max(dt2$mvstat)
#[1] 284.9248
#> max(dt2$mvp)
#[1] 60.73981
#> max(dt2$unip)
#[1] 61.68825
#> quantile(dt2$mvp)
#0%          25%          50%          75%         100%
#1.455461e-08 1.272388e-01 3.132055e-01 6.518595e-01 6.073981e+01
#> quantile(dt2$unip)
#0%          25%          50%          75%         100%
#0.001740662  0.438898616  0.703334810  1.099414413 61.688246139
> dim(dtsignif)
[1] 5464    7
> dim(dtlesssignif)
[1] 8277    7
> max(dt2$mvstat)
[1] 282.981
> max(dt2$mvp)
[1] 60.3192
> max(dt2$unip)
[1] 61.68825
> quantile(dt2$mvp)
          0%          25%          50%          75%         100%
1.457602e-08 1.281008e-01 3.148465e-01 6.546450e-01 6.031920e+01
> quantile(dt2$unip)
          0%          25%          50%          75%         100%
 0.001740662  0.437707136  0.700492701  1.096910013 61.688246139

#~~~

write.table(file="GIANT2010.vs2.dtsignif.vs1.SignCrrct.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.vs2.dtsignif.rs.vs1.SignCrrct.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.vs2.dtlesssignif.vs1.SignCrrct.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.vs2.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GIANT2010.vs2.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GIANT2010.vs2.RSS0.vs1.SignCrrct.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


