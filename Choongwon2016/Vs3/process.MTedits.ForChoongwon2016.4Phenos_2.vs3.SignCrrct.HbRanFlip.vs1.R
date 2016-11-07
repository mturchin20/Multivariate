#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# Choongwon2016 data

OxHb=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.OxHb.assoc.IncAllele.MatchedToHbRanFlip.formatted.txt.gz", header=T)
dHb=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.dHb.assoc.IncAllele.MatchedToHbRanFlip.formatted.txt.gz", header=T)
Pulse=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.Pulse.assoc.IncAllele.MatchedToHbRanFlip.formatted.txt.gz", header=T)
LvBrth=read.table("/mnt/gluster/data/external_private_supp/Choongwon2016/Tibetan.LvBrth.assoc.IncAllele.MatchedToHbRanFlip.formatted.txt.gz", header=T)

dim(OxHb)
dim(dHb)
dim(Pulse)
dim(LvBrth)

#~~~
> dim(OxHb)
[1] 3507568       6
> dim(dHb)
[1] 3507568       6
> dim(Pulse)
[1] 3506812       6
> dim(LvBrth)
[1] 3499631       6
#~~~

#GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }
GetZScoreAndSign <- function(x) { val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE); if (x < 0) { val1 <- val1 * -1 }; return(val1) }

OxHb <- cbind(OxHb, apply(as.matrix(OxHb$p), 1, GetZScoreAndSign))
dHb <- cbind(dHb, apply(as.matrix(dHb$p), 1, GetZScoreAndSign))
Pulse <- cbind(Pulse, apply(as.matrix(Pulse$p), 1, GetZScoreAndSign))
LvBrth <- cbind(LvBrth, apply(as.matrix(LvBrth$p), 1, GetZScoreAndSign))

OxHb.colnames <- colnames(OxHb)
OxHb.colnames[7] <- "GC.Zscore"
colnames(OxHb) <- OxHb.colnames
dHb.colnames <- colnames(dHb)
dHb.colnames[7] <- "GC.Zscore"
colnames(dHb) <- dHb.colnames
Pulse.colnames <- colnames(Pulse)
Pulse.colnames[7] <- "GC.Zscore"
colnames(Pulse) <- Pulse.colnames
LvBrth.colnames <- colnames(LvBrth)
LvBrth.colnames[7] <- "GC.Zscore"
colnames(LvBrth) <- LvBrth.colnames

head(OxHb)
head(dHb)
head(Pulse)
head(LvBrth)

#~~~
> head(OxHb)
  CHR       POS        SNP A1frq           p   N  GC.Zscore
1  10 100000625  rs7899632 0.578  0.12259150 916  1.5439869
2  10 100000645 rs61875309 0.940 -0.73402050 916 -0.3397823
3  10 100002399  rs8181398 0.945  0.39538810 918  0.8498867
4  10 100003242 rs12258651 0.121  0.19803270 917  1.2871767
5  10 100003785  rs1359508 0.462  0.05650671 919  1.9071072
6  10 100004360  rs1048754 0.940 -0.73402050 916 -0.3397823
> head(dHb)
  CHR       POS        SNP A1frq          p   N  GC.Zscore
1  10 100000625  rs7899632 0.578  0.3937060 916  0.8529158
2  10 100000645 rs61875309 0.940 -0.7739307 916 -0.2872372
3  10 100002399  rs8181398 0.945  0.2350974 918  1.1873302
4  10 100003242 rs12258651 0.879 -0.8898429 917 -0.1385030
5  10 100003785  rs1359508 0.462  0.2063892 919  1.2635567
6  10 100004360  rs1048754 0.940 -0.7739307 916 -0.2872372
> head(Pulse)
  CHR       POS        SNP A1frq           p   N  GC.Zscore
1  10 100000625  rs7899632 0.422 -0.30098060 915 -1.0343328
2  10 100000645 rs61875309 0.940 -0.06235068 915 -1.8637937
3  10 100002399  rs8181398 0.055 -0.38966980 917 -0.8602163
4  10 100003242 rs12258651 0.121  0.10113560 916  1.6393730
5  10 100003785  rs1359508 0.538 -0.60860010 918 -0.5120727
6  10 100004360  rs1048754 0.940 -0.06235068 915 -1.8637937
> head(LvBrth)
  CHR       POS        SNP A1frq           p   N  GC.Zscore
1  10 100000625  rs7899632 0.425 -0.17775840 974 -1.3476891
2  10 100000645 rs61875309 0.060  0.19837670 975  1.2861902
3  10 100002399  rs8181398 0.947  0.88547890 976  0.1440273
4  10 100003242 rs12258651 0.120  0.41649100 977  0.8125240
5  10 100003785  rs1359508 0.539 -0.04605376 979 -1.9949003
6  10 100004360  rs1048754 0.060  0.19837670 975  1.2861902
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=OxHb$SNP,Z.OxHb=OxHb$GC.Zscore,key="rs")
dt2 = data.table(rs=dHb$SNP,Z.dHb=dHb$GC.Zscore,key="rs")
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
#> dim(dt3)
#[1] 3499118       5
#> head(dt3)
#             rs     Z.OxHb      Z.dHb    Z.Pulse   Z.LvBrth
#1:  10_10003392  1.6221861  0.1242027  1.4752539 -2.3253383
#2: 10_100563437  0.8154985  0.2641953 -0.7627170  0.6706849
#3: 10_100720888  0.3357119 -0.5268399  0.0711521  0.1465233
#4: 10_102317557  0.7756681  1.8239832  1.7561462 -0.2953684
#5:  10_10236947 -0.5496300  1.3118860 -0.8250900 -0.3812933
#6: 10_103886793  1.9083270 -0.3550658 -0.3209355  2.2083069
> dim(dt3)
[1] 3499311       5
> head(dt3)  
             rs     Z.OxHb      Z.dHb    Z.Pulse   Z.LvBrth
1:  10_10003392  1.6221861  0.1242027  1.4752539 -2.3253383
2: 10_100563437 -0.8154985 -0.2641953  0.7627170 -0.6706849
3: 10_100720888  0.3357119 -0.5268399  0.0711521  0.1465233
4: 10_102317557 -0.7756681 -1.8239832 -1.7561462  0.2953684
5:  10_10236947 -0.5496300  1.3118860 -0.8250900 -0.3812933
6: 10_103886793  1.9083270 -0.3550658 -0.3209355  2.2083069
#~~~

OxHb.maxZ <- max(dt3$Z.OxHb[!is.infinite(dt3$Z.OxHb) & !is.na(dt3$Z.OxHb)])
dHb.maxZ <- max(dt3$Z.dHb[!is.infinite(dt3$Z.dHb) & !is.na(dt3$Z.dHb)])
Pulse.maxZ <- max(dt3$Z.Pulse[!is.infinite(dt3$Z.Pulse) & !is.na(dt3$Z.Pulse)])
LvBrth.maxZ <- max(dt3$Z.LvBrth[!is.infinite(dt3$Z.LvBrth) & !is.na(dt3$Z.LvBrth)])
maxZ <- max(c(OxHb.maxZ, dHb.maxZ, Pulse.maxZ, LvBrth.maxZ))

OxHb.maxZ
dHb.maxZ
Pulse.maxZ
LvBrth.maxZ
maxZ

#~~~
#> OxHb.maxZ
#[1] 5.824933
#> dHb.maxZ
#[1] 5.26906
#> Pulse.maxZ
#[1] 4.621363
#> LvBrth.maxZ
#[1] 5.937577
#> maxZ
#[1] 5.937577
> OxHb.maxZ
[1] 5.824933
> dHb.maxZ
[1] 5.26906
> Pulse.maxZ
[1] 4.621363
> LvBrth.maxZ
[1] 4.993735
> maxZ
[1] 5.824933
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt3$Z.OxHb <- apply(as.matrix(dt3$Z.OxHb), 1, replaceInf)
dt3$Z.dHb <- apply(as.matrix(dt3$Z.dHb), 1, replaceInf)
dt3$Z.Pulse <- apply(as.matrix(dt3$Z.Pulse), 1, replaceInf)
dt3$Z.LvBrth <- apply(as.matrix(dt3$Z.LvBrth), 1, replaceInf)

#~~~
> 
> dt3$Z.OxHb <- apply(as.matrix(dt3$Z.OxHb), 1, replaceInf)
> dt3$Z.dHb <- apply(as.matrix(dt3$Z.dHb), 1, replaceInf)
> dt3$Z.Pulse <- apply(as.matrix(dt3$Z.Pulse), 1, replaceInf)
> dt3$Z.LvBrth <- apply(as.matrix(dt3$Z.LvBrth), 1, replaceInf)
> 
#~~~

attach(dt3)
nullset = (abs(dt3$Z.dHb)<2) & (abs(dt3$Z.OxHb)<2) & (abs(dt3$Z.Pulse)<2) & (abs(dt3$Z.LvBrth)<2) #extract null Z values
Z = cbind(Z.dHb,Z.OxHb,Z.Pulse,Z.LvBrth)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
#> dim(Z)
#[1] 3499118       4
#> dim(Znull)
#[1] 2883649       4
#> RSS0
#               Z.dHb       Z.OxHb      Z.Pulse    Z.LvBrth
#Z.dHb     1.00000000 -0.380710031  0.116037678 -0.01578817
#Z.OxHb   -0.38071003  1.000000000  0.008560982 -0.03344946
#Z.Pulse   0.11603768  0.008560982  1.000000000 -0.10697349
#Z.LvBrth -0.01578817 -0.033449455 -0.106973488  1.00000000
> dim(Z)
[1] 3499311       4
> dim(Znull)
[1] 2883811       4
> RSS0
               Z.dHb      Z.OxHb     Z.Pulse    Z.LvBrth
Z.dHb     1.00000000  0.01181027  0.14204221 -0.02747314
Z.OxHb    0.01181027  1.00000000  0.07792366 -0.04976397
Z.Pulse   0.14204221  0.07792366  1.00000000 -0.10999706
Z.LvBrth -0.02747314 -0.04976397 -0.10999706  1.00000000

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
#[1] 30  8
#> dim(dtlesssignif)
#[1] 183   8
#> dim(dtlesslesssignif)
#[1] 5966    8
#> max(dt3$mvstat)
#[1] 44.73716
#> max(dt3$mvp)
#[1] 8.977875
#> max(dt3$unip)
#[1] 8.538703
#> quantile(dt3$mvp)
#          0%          25%          50%          75%         100% 
#3.077121e-05 2.494135e-01 5.137570e-01 9.303142e-01 8.977875e+00 
#> quantile(dt3$unip)
#        0%        25%        50%        75%       100% 
#0.01621076 0.53528889 0.80345381 1.17114227 8.53870336 
> dim(dtsignif)
[1] 10  8
> dim(dtlesssignif)
[1] 40  8
> dim(dtlesslesssignif)
[1] 2255    8
> max(dt3$mvstat)
[1] 36.80627
> max(dt3$mvp)
[1] 7.296126
> max(dt3$unip)
[1] 8.538703
> quantile(dt3$mvp)
          0%          25%          50%          75%         100%
2.274456e-05 2.332275e-01 4.759281e-01 8.472121e-01 7.296126e+00
> quantile(dt3$unip)
        0%        25%        50%        75%       100%
0.01621076 0.53528800 0.80345381 1.17114227 8.53870336
#~~~

write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.dtsignif.vs3.SignCrrct.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.dtsignif.rs.vs3.SignCrrct.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.dtlesssignif.vs3.SignCrrct.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.dtlesssignif.rs.vs3.SignCrrct.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.dtlesslesssignif.vs3.SignCrrct.vs1.txt",dtlesslesssignif,sep=",",row.names=F,quote=F)
#write.table(file="Choongwon2016.4Phenos_2.dtlesssignif.rs.vs3.SignCrrct.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="Choongwon2016.RanHalfFlip.4Phenos_2.RSS0.vs3.SignCrrct.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


