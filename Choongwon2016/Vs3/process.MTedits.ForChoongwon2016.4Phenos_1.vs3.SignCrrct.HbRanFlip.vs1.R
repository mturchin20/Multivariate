#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# Choongwon2016 data

Hb=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.Hb.assoc.IncAllele.RanHalfFlip.formatted.txt.gz", header=T)
Sat=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.Sat.assoc.IncAllele.MatchedToHbRanFlip.formatted.txt.gz", header=T)
Pulse=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.Pulse.assoc.IncAllele.MatchedToHbRanFlip.formatted.txt.gz", header=T)
LvBrth=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.LvBrth.assoc.IncAllele.MatchedToHbRanFlip.formatted.txt.gz", header=T)

dim(Hb)
dim(Sat)
dim(Pulse)
dim(LvBrth)

#~~~
> dim(Hb)
[1] 3507568       6
> dim(Sat)
[1] 3507568       6
> dim(Pulse)
[1] 3506812       6
> dim(LvBrth)
[1] 3499631       6
#~~~

#GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }
GetZScoreAndSign <- function(x) { val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE); if (x < 0) { val1 <- val1 * -1 }; return(val1) }

Hb <- cbind(Hb, apply(as.matrix(Hb$p), 1, GetZScoreAndSign))
Sat <- cbind(Sat, apply(as.matrix(Sat$p), 1, GetZScoreAndSign))
Pulse <- cbind(Pulse, apply(as.matrix(Pulse$p), 1, GetZScoreAndSign))
LvBrth <- cbind(LvBrth, apply(as.matrix(LvBrth$p), 1, GetZScoreAndSign))

Hb.colnames <- colnames(Hb)
Hb.colnames[7] <- "GC.Zscore"
colnames(Hb) <- Hb.colnames
Sat.colnames <- colnames(Sat)
Sat.colnames[7] <- "GC.Zscore"
colnames(Sat) <- Sat.colnames
Pulse.colnames <- colnames(Pulse)
Pulse.colnames[7] <- "GC.Zscore"
colnames(Pulse) <- Pulse.colnames
LvBrth.colnames <- colnames(LvBrth)
LvBrth.colnames[7] <- "GC.Zscore"
colnames(LvBrth) <- LvBrth.colnames

head(Hb)
head(Sat)
head(Pulse)
head(LvBrth)

#~~~
> head(Hb)
  CHR    POS       SNP A1frq         p   N GC.Zscore
1   1 752566 rs3094315 0.149 0.4168330 921 0.8119279
2   1 752721 rs3131972 0.276 0.7333741 921 0.3406407
3   1 752894 rs3131971 0.149 0.3131140 903 1.0087094
4   1 754334 rs3131967 0.268 0.4013396 897 0.8392312
5   1 754503 rs3115859 0.277 0.8260580 908 0.2197601
6   1 754964 rs3131966 0.277 0.8125771 907 0.2371027
> head(Sat)
  CHR       POS        SNP A1frq          p   N     GC.Zscore
1  10 100000625  rs7899632 0.422 -0.8755691 916 -1.565886e-01
2  10 100000645 rs61875309 0.060 -0.7304266 916 -3.445581e-01
3  10 100002399  rs8181398 0.055 -0.3570072 918 -9.210828e-01
4  10 100003242 rs12258651 0.879 -0.9999230 917 -9.650519e-05
5  10 100003785  rs1359508 0.538 -0.5861611 919 -5.444075e-01
6  10 100004360  rs1048754 0.060 -0.7304266 916 -3.445581e-01
> head(Pulse)
  CHR       POS        SNP A1frq           p   N  GC.Zscore
1  10 100000625  rs7899632 0.422 -0.30098060 915 -1.0343328
2  10 100000645 rs61875309 0.940  0.06235068 915  1.8637937
3  10 100002399  rs8181398 0.055 -0.38966980 917 -0.8602163
4  10 100003242 rs12258651 0.121  0.10113560 916  1.6393730
5  10 100003785  rs1359508 0.538 -0.60860010 918 -0.5120727
6  10 100004360  rs1048754 0.940  0.06235068 915  1.8637937
> head(LvBrth)
  CHR       POS        SNP A1frq           p   N  GC.Zscore
1  10 100000625  rs7899632 0.425 -0.17775840 974 -1.3476891
2  10 100000645 rs61875309 0.060 -0.19837670 975 -1.2861902
3  10 100002399  rs8181398 0.947  0.88547890 976  0.1440273
4  10 100003242 rs12258651 0.120  0.41649100 977  0.8125240
5  10 100003785  rs1359508 0.539 -0.04605376 979 -1.9949003
6  10 100004360  rs1048754 0.060 -0.19837670 975 -1.2861902
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=Hb$SNP,Z.Hb=Hb$GC.Zscore,key="rs")
dt2 = data.table(rs=Sat$SNP,Z.Sat=Sat$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=Pulse$SNP,Z.Pulse=Pulse$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=LvBrth$SNP,Z.LvBrth=LvBrth$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)

dim(dt3)
head(dt3)

#~~~
> dim(dt3)
[1] 3499311       5
> head(dt3)
             rs      Z.Hb      Z.Sat    Z.Pulse   Z.LvBrth
1:  10_10003392 1.5791838  0.3494072  1.4752539 -2.3253383
2: 10_100563437 0.8109792 -0.2097422 -0.7627170  0.6706849
3: 10_100720888 0.1009797  0.8629236  0.0711521  0.1465233
4: 10_102317557 1.4871436 -1.6563970  1.7561462 -0.2953684
5:  10_10236947 0.2015823 -1.4576396 -0.8250900 -0.3812933
6: 10_103886793 1.7126065  0.9075181 -0.3209355  2.2083069
#~~~

Hb.maxZ <- max(dt3$Z.Hb[!is.infinite(dt3$Z.Hb) & !is.na(dt3$Z.Hb)])
Sat.maxZ <- max(dt3$Z.Sat[!is.infinite(dt3$Z.Sat) & !is.na(dt3$Z.Sat)])
Pulse.maxZ <- max(dt3$Z.Pulse[!is.infinite(dt3$Z.Pulse) & !is.na(dt3$Z.Pulse)])
LvBrth.maxZ <- max(dt3$Z.LvBrth[!is.infinite(dt3$Z.LvBrth) & !is.na(dt3$Z.LvBrth)])
maxZ <- max(c(Hb.maxZ, Sat.maxZ, Pulse.maxZ, LvBrth.maxZ))

Hb.maxZ
Sat.maxZ
Pulse.maxZ
LvBrth.maxZ
maxZ

#~~~
#> Hb.maxZ
#[1] 5.126442
#> Sat.maxZ 
#[1] 4.842598
#> Pulse.maxZ
#[1] 4.621363
#> LvBrth.maxZ
#[1] 5.937577
#> maxZ
#[1] 5.937577
> Hb.maxZ
[1] 5.126442
> Sat.maxZ
[1] 4.842598
> Pulse.maxZ
[1] 4.621363
> LvBrth.maxZ
[1] 4.993735
> maxZ
[1] 5.126442
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt3$Z.Hb <- apply(as.matrix(dt3$Z.Hb), 1, replaceInf)
dt3$Z.Sat <- apply(as.matrix(dt3$Z.Sat), 1, replaceInf)
dt3$Z.Pulse <- apply(as.matrix(dt3$Z.Pulse), 1, replaceInf)
dt3$Z.LvBrth <- apply(as.matrix(dt3$Z.LvBrth), 1, replaceInf)

#~~~
> 
> dt3$Z.Hb <- apply(as.matrix(dt3$Z.Hb), 1, replaceInf)
> dt3$Z.Sat <- apply(as.matrix(dt3$Z.Sat), 1, replaceInf)
> dt3$Z.Pulse <- apply(as.matrix(dt3$Z.Pulse), 1, replaceInf)
> dt3$Z.LvBrth <- apply(as.matrix(dt3$Z.LvBrth), 1, replaceInf)
> 
#~~~

attach(dt3)
nullset = (abs(dt3$Z.Sat)<2) & (abs(dt3$Z.Hb)<2) & (abs(dt3$Z.Pulse)<2) & (abs(dt3$Z.LvBrth)<2) #extract null Z values
Z = cbind(Z.Sat,Z.Hb,Z.Pulse,Z.LvBrth)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
#> dim(Z)
#[1] 3499311       4
#> dim(Znull)
#[1] 2882127       4
#> RSS0
#                Z.Sat        Z.Hb     Z.Pulse     Z.LvBrth
#Z.Sat     1.000000000 -0.07228911 -0.10196383  0.003547562
#Z.Hb     -0.072289107  1.00000000  0.06869758 -0.031857655
#Z.Pulse  -0.101963828  0.06869758  1.00000000 -0.107232640
#Z.LvBrth  0.003547562 -0.03185766 -0.10723264  1.000000000
> dim(Z)
[1] 3499311       4
> dim(Znull)
[1] 2882127       4
> RSS0
                Z.Sat        Z.Hb    Z.Pulse     Z.LvBrth
Z.Sat     1.000000000 -0.12856576 -0.1115199  0.007571157
Z.Hb     -0.128565760  1.00000000  0.1209940 -0.049328528
Z.Pulse  -0.111519919  0.12099405  1.0000000 -0.110409289
Z.LvBrth  0.007571157 -0.04932853 -0.1104093  1.000000000
#~~~

mvstat =  rowSums(Z * (Z %*% RSS0inv)) # comptues Z RSS0^-1 Z'
dt3$mvstat = mvstat
statchi = -log10(exp(1))*pchisq(mvstat,df=3,log.p=TRUE,lower.tail=FALSE)
dt3$mvp = statchi

maxZ2 = apply(Z^2,1,max)
max.unip = -log10(exp(1))*pchisq(maxZ2,df=1,log.p=TRUE, lower.tail=FALSE)
dt3$unip = max.unip

dtsignif = dt3[dt3$mvp>-log10(5e-8) | dt3$unip>-log10(5e-8),]
dtlesssignif = dt3[dt3$mvp>-log10(1e-6) | dt3$unip>-log10(1e-6),]
dtlesslesssignif = dt3[dt3$mvp>-log10(1e-4) | dt3$unip>-log10(1e-4),]

dim(dtsignif)
dim(dtlesssignif)
dim(dtlesslesssignif)
max(dt3$mvstat)
max(dt3$mvp)
max(dt3$unip)
quantile(dt3$mvp)
quantile(dt3$unip)

#~~~
#> dim(dtsignif)
#[1] 2 8
#> dim(dtlesssignif)
#[1] 15  8
#> dim(dtlesslesssignif)
#[1] 2247    8
#> max(dt3$mvstat)
#[1] 36.67802
#> max(dt3$mvp)
#[1] 7.268998
#> max(dt3$unip)
#[1] 8.538703
#> quantile(dt3$mvp)
#          0%          25%          50%          75%         100% 
#3.376168e-06 2.306524e-01 4.716323e-01 8.423546e-01 7.268998e+00 
#> quantile(dt3$unip)
#         0%         25%         50%         75%        100% 
#0.008190191 0.533773771 0.802543278 1.170578051 8.538703361 
> dim(dtsignif)
[1] 2 8
> dim(dtlesssignif)
[1] 15  8
> dim(dtlesslesssignif)
[1] 2133    8
> max(dt3$mvstat)
[1] 36.94829
> max(dt3$mvp)
[1] 7.326172
> max(dt3$unip)
[1] 8.538703
> quantile(dt3$mvp)
          0%          25%          50%          75%         100% 
3.147973e-06 2.316149e-01 4.738096e-01 8.452084e-01 7.326172e+00 
> quantile(dt3$unip)
         0%         25%         50%         75%        100% 
0.008190191 0.533773771 0.802543278 1.170578051 8.538703361 

#~~~

write.table(file="Choongwon2016.RanHalfFlip.4Phenos_1.dtsignif.vs3.SignCrrct.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_1.dtsignif.rs.vs3.SignCrrct.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_1.dtlesssignif.vs3.SignCrrct.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_1.dtlesssignif.rs.vs3.SignCrrct.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_1.dtlesslesssignif.vs3.SignCrrct.vs1.txt",dtlesslesssignif,sep=",",row.names=F,quote=F)
#write.table(file="Choongwon2016.4Phenos_1.dtlesssignif.rs.vs2.SignCrrct.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="Choongwon2016.RanHalfFlip.4Phenos_1.RSS0.vs3.SignCrrct.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


