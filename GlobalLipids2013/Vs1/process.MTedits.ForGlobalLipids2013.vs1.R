#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# GlobalLipids2013 data

/data/external_public/GlobalLipid
-rw-r--r-- 1 xiangzhu xiangzhu 58851380 Jun  4  2014 jointGwasMc_LDL.txt.gz
-rw-r--r-- 1 xiangzhu xiangzhu 59025777 Jun  4  2014 jointGwasMc_HDL.txt.gz
-rw-r--r-- 1 xiangzhu xiangzhu 58866740 Jun  4  2014 jointGwasMc_TG.txt.gz
-rw-r--r-- 1 xiangzhu xiangzhu 59110550 Jun  4  2014 jointGwasMc_TC.txt.gz

HDL=read.table("/data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz", header=T)
LDL=read.table("/data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz", header=T)
TG=read.table("/data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz", header=T)
TC=read.table("/data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz", header=T)

dim(HDL)
dim(LDL)
dim(TC)
dim(TG)

#~~~
#> dim(HDL)
#[1] 2447441      10
#> dim(LDL)
#[1] 2437751      10
#> dim(TC)
#[1] 2446981      10
#> dim(TG)
#[1] 2439432      10
#~~~

GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }

HDL <- cbind(HDL, apply(as.matrix(HDL$P.value), 2, GetZScore))
LDL <- cbind(LDL, apply(as.matrix(LDL$P.value), 2, GetZScore))
TC <- cbind(TC, apply(as.matrix(TC$P.value), 2, GetZScore))
TG <- cbind(TG, apply(as.matrix(TG$P.value), 2, GetZScore))

HDL.colnames <- colnames(HDL)
HDL.colnames[11] <- "GC.Zscore"
colnames(HDL) <- HDL.colnames
LDL.colnames <- colnames(LDL)
LDL.colnames[11] <- "GC.Zscore"
colnames(LDL) <- LDL.colnames
TC.colnames <- colnames(TC)
TC.colnames[11] <- "GC.Zscore"
colnames(TC) <- TC.colnames
TG.colnames <- colnames(TG)
TG.colnames[11] <- "GC.Zscore"
colnames(TG) <- TG.colnames

head(HDL)
head(LDL)
head(TC)
head(TG)

#~~~
#> head(HDL)
#SNP_hg18        SNP_hg19      rsid A1 A2   beta     se     N P.value
#1  chr10:10000135   chr10:9960129 rs4747841  g  a 0.0026 0.0048 93561 0.75380
#2  chr10:10000265   chr10:9960259 rs4749917  t  c 0.0028 0.0048 93561 0.72950
#3 chr10:100002729 chr10:100012739  rs737656  g  a 0.0098 0.0049 94311 0.09234
#4 chr10:100002880 chr10:100012890  rs737657  g  a 0.0102 0.0049 94311 0.09340
#5 chr10:100003553 chr10:100013563 rs7086391  c  t 0.0005 0.0062 94311 0.87230
#6 chr10:100003805 chr10:100013815  rs878177  t  c 0.0096 0.0050 94311 0.12400
#Freq.A1.1000G.EUR GC.Zscore
#1            0.5092 0.3136327
#2            0.5092 0.3457907
#3            0.6794 1.6831813
#4            0.6794 1.6777291
#5            0.7810 0.1607377
#6            0.3483 1.5381989
#> head(LDL)
#SNP_hg18        SNP_hg19      rsid A1 A2   beta     se     N P.value
#1  chr10:10000135   chr10:9960129 rs4747841  a  g 0.0037 0.0052 89138 0.71580
#2  chr10:10000265   chr10:9960259 rs4749917  c  t 0.0033 0.0052 89138 0.77480
#3 chr10:100002729 chr10:100012739  rs737656  a  g 0.0099 0.0054 89888 0.04000
#4 chr10:100002880 chr10:100012890  rs737657  a  g 0.0084 0.0054 89888 0.08428
#5 chr10:100003553 chr10:100013563 rs7086391  c  t 0.0075 0.0067 89888 0.26890
#6 chr10:100003805 chr10:100013815  rs878177  c  t 0.0073 0.0055 89888 0.13760
#Freq.A1.1000G.EUR GC.Zscore
#1            0.4908 0.3640777
#2            0.4908 0.2861020
#3            0.3206 2.0537489
#4            0.3206 1.7263748
#5            0.7810 1.1055993
#6            0.6517 1.4847880
#> head(TC)
#SNP_hg18        SNP_hg19      rsid A1 A2   beta     se     N P.value
#1  chr10:10000135   chr10:9960129 rs4747841  a  g 0.0026 0.0051 93845  0.6530
#2  chr10:10000265   chr10:9960259 rs4749917  c  t 0.0024 0.0051 93845  0.6915
#3 chr10:100002729 chr10:100012739  rs737656  a  g 0.0071 0.0052 94595  0.1653
#4 chr10:100002880 chr10:100012890  rs737657  a  g 0.0062 0.0052 94595  0.2550
#5 chr10:100003553 chr10:100013563 rs7086391  c  t 0.0087 0.0066 94595  0.1513
#6 chr10:100003805 chr10:100013815  rs878177  c  t 0.0024 0.0053 94595  0.7058
#Freq.A1.1000G.EUR GC.Zscore
#1            0.4908 0.4495985
#2            0.4908 0.3968203
#3            0.3206 1.3874651
#4            0.3206 1.1382886
#5            0.7810 1.4349547
#6            0.6517 0.3775028
#> head(TG)
#SNP_hg18        SNP_hg19      rsid A1 A2   beta     se     N P.value
#1  chr10:10000135   chr10:9960129 rs4747841  a  g 0.0020 0.0047 90263 0.75670
#2  chr10:10000265   chr10:9960259 rs4749917  c  t 0.0020 0.0047 90263 0.76120
#3 chr10:100002729 chr10:100012739  rs737656  a  g 0.0094 0.0048 91013 0.06856
#4 chr10:100002880 chr10:100012890  rs737657  a  g 0.0094 0.0048 91013 0.07994
#5 chr10:100003553 chr10:100013563 rs7086391  c  t 0.0094 0.0061 91013 0.11250
#6 chr10:100003805 chr10:100013815  rs878177  c  t 0.0026 0.0049 91013 0.73140
#Freq.A1.1000G.EUR GC.Zscore
#1            0.4908 0.3098172
#2            0.4908 0.3039054
#3            0.3206 1.8213083
#4            0.3206 1.7510343
#5            0.7810 1.5870558
#6            0.6517 0.3432638
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=HDL$SNP_hg19,Z.HDL=HDL$GC.Zscore,key="rs")
dt2 = data.table(rs=LDL$SNP_hg19,Z.LDL=LDL$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=TC$SNP_hg19,Z.TC=TC$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=TG$SNP_hg19,Z.TG=TG$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)

dim(dt3)
head(dt3)

#~~~
#> dim(dt3)
#[1] 2437102       5
#> head(dt3)
#rs    Z.HDL     Z.LDL        Z.TC      Z.TG
#1: chr10:100004441 1.652464 1.6481592 1.138049053 1.6614156
#2: chr10:100004906 1.957830 0.8111144 0.005263944 0.6968447
#3: chr10:100005282 1.975601 0.8216081 0.035977876 0.6373452
#4: chr10:100008436 1.650601 1.7294421 1.153268976 1.8012086
#5: chr10:100012739 1.683181 2.0537489 1.387465054 1.8213083
#6: chr10:100012890 1.677729 1.7263748 1.138288582 1.7510343
#~~~

HDL.maxZ <- max(dt3$Z.HDL[!is.infinite(dt3$Z.HDL) & !is.na(dt3$Z.HDL)])
LDL.maxZ <- max(dt3$Z.LDL[!is.infinite(dt3$Z.LDL) & !is.na(dt3$Z.LDL)])
TC.maxZ <- max(dt3$Z.TC[!is.infinite(dt3$Z.TC) & !is.na(dt3$Z.TC)])
TG.maxZ <- max(dt3$Z.TG[!is.infinite(dt3$Z.TG) & !is.na(dt3$Z.TG)])
maxZ <- max(c(HDL.maxZ, LDL.maxZ, TC.maxZ, TG.maxZ))

HDL.maxZ
LDL.maxZ
TC.maxZ
TG.maxZ
maxZ

#~~~
#> HDL.maxZ
#[1] 38.24664
#> LDL.maxZ
#[1] 35.27109
#> TC.maxZ
#[1] 35.9827
#> TG.maxZ
#[1] 33.74916
#> maxZ
#[1] 38.24664
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt3$Z.HDL <- apply(as.matrix(dt3$Z.HDL), 1, replaceInf)
dt3$Z.LDL <- apply(as.matrix(dt3$Z.LDL), 1, replaceInf)
dt3$Z.TC <- apply(as.matrix(dt3$Z.TC), 1, replaceInf)
dt3$Z.TG <- apply(as.matrix(dt3$Z.TG), 1, replaceInf)

#~~~
> dt3$Z.HDL <- apply(as.matrix(dt3$Z.HDL), 1, replaceInf)
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
> dt3$Z.LDL <- apply(as.matrix(dt3$Z.LDL), 1, replaceInf)
[1] "yaya1"
[1] "yaya1"
[1] "yaya1"
> dt3$Z.TC <- apply(as.matrix(dt3$Z.TC), 1, replaceInf)
> dt3$Z.TG <- apply(as.matrix(dt3$Z.TG), 1, replaceInf)
#~~~

attach(dt3)
nullset = (abs(dt3$Z.LDL)<2) & (abs(dt3$Z.HDL)<2) & (abs(dt3$Z.TC)<2) & (abs(dt3$Z.TG)<2) #extract null Z values
Z = cbind(Z.LDL,Z.HDL,Z.TC,Z.TG)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
#> dim(Z)
#[1] 2437102       4
#> dim(Znull)
#[1] 2058066       4
#> RSS0
#Z.LDL       Z.HDL       Z.TC       Z.TG
#Z.LDL 1.000000000 0.003721305 0.61570864 0.01629742
#Z.HDL 0.003721305 1.000000000 0.01527727 0.09523876
#Z.TC  0.615708639 0.015277266 1.00000000 0.06149693
#Z.TG  0.016297416 0.095238755 0.06149693 1.00000000
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

dim(dtsignif)
dim(dtlesssignif)
max(dt3$mvstat)
max(dt3$mvp)
max(dt3$unip)
quantile(dt3$mvp)
quantile(dt3$unip)

#~~~
#> dim(dtsignif)
#[1] 9664    8
#> dim(dtlesssignif)
#[1] 12462     8
#> max(dt3$mvstat)
#[1] 1860.333
#> max(dt3$mvp)
#[1] 400.9972
#> max(dt3$unip)
#[1] 319.3251
#> quantile(dt3$mvp)
#0%          25%          50%          75%         100%
#5.766978e-08 6.505184e-02 1.941311e-01 4.683504e-01 4.009972e+02
#> quantile(dt3$unip)
#0%          25%          50%          75%         100%
#0.01032775   0.44442193   0.70707970   1.08900228 319.32512192
#~~~

write.table(file="GlobalLipids2013.dtsignif.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GlobalLipids2013.dtsignif.rs.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GlobalLipids2013.dtlesssignif.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GlobalLipids2013.dtlesssignif.rs.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GlobalLipids2013.dtlesssignif.rs.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GlobalLipids2013.RSS0.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


