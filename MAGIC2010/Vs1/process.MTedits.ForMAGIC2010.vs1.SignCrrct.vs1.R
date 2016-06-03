#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# MAGIC2010 data

#~~~
#[  mturchin20@spudling26  /data/external_public/GlobalLipids2010]$ls -lrt /mnt/gluster/data/external_public_supp/MAGIC2010
#total 426306
#-rw-rw-r-- 1 mturchin20 mturchin20 26977566 Aug  8  2011 MAGIC_ln_HOMA-IR.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 26723045 Aug  8  2011 MAGIC_ln_HOMA-B.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 26967120 Aug  8  2011 MAGIC_ln_FastingInsulin.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 27054992 Aug  6 14:37 MAGIC_FastingGlucose.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 40863290 Aug  6 23:01 MAGIC_ln_HOMA-B.wUCSCGB_dbSNP130.vs1.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 41126386 Aug  6 23:01 MAGIC_ln_HOMA-IR.wUCSCGB_dbSNP130.vs1.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 41136996 Aug  6 23:12 MAGIC_ln_FastingInsulin.wUCSCGB_dbSNP130.vs1.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 41278571 Aug  6 23:12 MAGIC_FastingGlucose.wUCSCGB_dbSNP130.vs1.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 40863398 Aug  7 00:02 MAGIC_ln_HOMA-B.wUCSCGB_snp129.vs3.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 41126498 Aug  7 00:02 MAGIC_ln_HOMA-IR.wUCSCGB_snp129.vs3.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 41137166 Aug  7 00:02 MAGIC_ln_FastingInsulin.wUCSCGB_snp129.vs3.txt.gz
#-rw-rw-r-- 1 mturchin20 mturchin20 41278669 Aug  7 00:03 MAGIC_FastingGlucose.wUCSCGB_snp129.vs3.txt.gz
#~~~

#FastIns=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_ln_FastingInsulin.wUCSCGB_snp129.vs3.txt.gz",header=T)
#FastGlu=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_FastingGlucose.wUCSCGB_snp129.vs3.txt.gz",header=T)
#HOMA_B=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_ln_HOMA-B.wUCSCGB_snp129.vs3.txt.gz",header=T)
#HOMA_IR=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_ln_HOMA-IR.wUCSCGB_snp129.vs3.txt.gz",header=T)
FastIns=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_ln_FastingInsulin.wHapMap22.vs3.IncAllele.formatted.txt.gz",header=T)
FastGlu=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_FastingGlucose.wHapMap22.vs3.MatchedToInsln.IncAllele.formatted.txt.gz",header=T)
HOMA_B=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_ln_HOMA-B.wHapMap22.vs3.MatchedToInsln.IncAllele.formatted.txt.gz",header=T)
HOMA_IR=read.table("/mnt/gluster/data/external_public_supp/MAGIC2010/MAGIC_ln_HOMA-IR.wHapMap22.vs3.MatchedToInsln.IncAllele.formatted.txt.gz",header=T)

dim(FastIns)
dim(FastGlu)
dim(HOMA_B)
dim(HOMA_IR)

#~~~
#> dim(HOMA_B)
#[1] 2456945       8
#> dim(HOMA_IR)
#[1] 2458073       8
#> dim(FastIns)
#[1] 2461105       8
#> dim(FastGlu)
#[1] 2470476       8
> dim(FastIns)
[1] 2460675       6
> dim(FastGlu)
[1] 2459676       6
> dim(HOMA_B)
[1] 2455775       6
> dim(HOMA_IR)
[1] 2456761       6
#~~~

#GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }
GetZScoreAndSign <- function(x) { val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE); if (x < 0) { val1 <- val1 * -1 }; return(val1) }

FastIns <- cbind(FastIns, apply(as.matrix(FastIns$pvalue), 1, GetZScoreAndSign))
FastGlu <- cbind(FastGlu, apply(as.matrix(FastGlu$pvalue), 1, GetZScoreAndSign))
HOMA_B <- cbind(HOMA_B, apply(as.matrix(HOMA_B$pvalue), 1, GetZScoreAndSign))
HOMA_IR <- cbind(HOMA_IR, apply(as.matrix(HOMA_IR$pvalue), 1, GetZScoreAndSign))

FastIns.colnames <- colnames(FastIns)
FastIns.colnames[7] <- "GC.Zscore"
colnames(FastIns) <- FastIns.colnames
FastGlu.colnames <- colnames(FastGlu)
FastGlu.colnames[7] <- "GC.Zscore"
colnames(FastGlu) <- FastGlu.colnames
HOMA_B.colnames <- colnames(HOMA_B)
HOMA_B.colnames[7] <- "GC.Zscore"
colnames(HOMA_B) <- HOMA_B.colnames
HOMA_IR.colnames <- colnames(HOMA_IR)
HOMA_IR.colnames[7] <- "GC.Zscore"
colnames(HOMA_IR) <- HOMA_IR.colnames

head(FastIns)
head(FastGlu)
head(HOMA_B)
head(HOMA_IR)

##~~~
#> head(HOMA_B)
#snp effect_allele other_allele   maf  effect stderr pvalue
#1       rs10             a            c 0.033  0.0094 0.0100 0.3533
#2  rs1000000             a            g 0.274  0.0006 0.0039 0.8824
#3 rs10000010             t            c 0.446  0.0035 0.0033 0.2925
#4 rs10000012             c            g 0.146 -0.0001 0.0050 0.9834
#5 rs10000013             a            c 0.167 -0.0019 0.0038 0.6183
#6 rs10000017             t            c 0.223 -0.0029 0.0042 0.4854
#ChrBP  GC.Zscore
#1   7_92221824 0.92820739
#2 12_125456933 0.14792748
#3   4_21227772 1.05265312
#4    4_1347325 0.02080652
#5   4_36901464 0.49826113
#6   4_84997149 0.69764377
#> head(HOMA_IR)
#snp effect_allele other_allele   maf  effect stderr pvalue
#1       rs10             a            c 0.033  0.0025 0.0130 0.8436
#2  rs1000000             a            g 0.274  0.0041 0.0048 0.3917
#3 rs10000010             t            c 0.446  0.0024 0.0040 0.5571
#4 rs10000012             c            g 0.146 -0.0005 0.0060 0.9272
#5 rs10000013             a            c 0.167 -0.0047 0.0046 0.3159
#6 rs10000017             t            c 0.223 -0.0063 0.0052 0.2266
#ChrBP  GC.Zscore
#1   7_92221824 0.19729077
#2 12_125456933 0.85653848
#3   4_21227772 0.58715442
#4    4_1347325 0.09136824
#5   4_36901464 1.00291888
#6   4_84997149 1.20916360
#> head(FastIns)
#snp effect_allele other_allele   maf  effect stderr pvalue
#1       rs10             a            c 0.033  0.0014 0.0120 0.9058
#2  rs1000000             a            g 0.274  0.0028 0.0046 0.5430
#3 rs10000010             t            c 0.446  0.0034 0.0039 0.3826
#4 rs10000012             c            g 0.146 -0.0005 0.0057 0.9293
#5 rs10000013             a            c 0.167 -0.0047 0.0044 0.2955
#6 rs10000017             t            c 0.223 -0.0040 0.0050 0.4222
#ChrBP  GC.Zscore
#1   7_92221824 0.11833781
#2 12_125456933 0.60828269
#3   4_21227772 0.87311573
#4    4_1347325 0.08872558
#5   4_36901464 1.04613220
#6   4_84997149 0.80261032
#> head(FastGlu)
#snp effect_allele other_allele   maf  effect stderr pvalue
#1       rs10             a            c 0.033 -0.0120 0.0120 0.3211
#2  rs1000000             a            g 0.274  0.0030 0.0044 0.5052
#3 rs10000010             t            c 0.446  0.0003 0.0037 0.9311
#4 rs10000012             c            g 0.146  0.0043 0.0054 0.4268
#5 rs10000013             a            c 0.167 -0.0016 0.0043 0.7016
#6 rs10000017             t            c 0.223 -0.0018 0.0047 0.7085
#ChrBP  GC.Zscore
#1   7_92221824 0.99219994
#2 12_125456933 0.66633030
#3   4_21227772 0.08646095
#4    4_1347325 0.79467943
#5   4_36901464 0.38316153
#6   4_84997149 0.37387139
> head(FastIns)
    Chr        BP        snp   raf pvalue     N  GC.Zscore
1  chr7  92221824       rs10 0.033 0.9058 46186 0.11833781
2 chr12 125456933  rs1000000 0.274 0.5430 46186 0.60828269
3  chr4  21227772 rs10000010 0.446 0.3826 46186 0.87311573
4  chr4   1347325 rs10000012 0.854 0.9293 46186 0.08872558
5  chr4  36901464 rs10000013 0.833 0.2955 46186 1.04613220
6  chr4  84997149 rs10000017 0.777 0.4222 46186 0.80261032
> head(FastGlu)
    Chr        BP        snp   raf  pvalue     N   GC.Zscore
1  chr7  92221824       rs10 0.967 -0.3211 46186 -0.99219994
2 chr12 125456933  rs1000000 0.274  0.5052 46186  0.66633030
3  chr4  21227772 rs10000010 0.446  0.9311 46186  0.08646095
4  chr4   1347325 rs10000012 0.146 -0.4268 46186 -0.79467943
5  chr4  36901464 rs10000013 0.833  0.7016 46186  0.38316153
6  chr4  84997149 rs10000017 0.777  0.7085 46186  0.37387139
> head(HOMA_B)
    Chr        BP        snp   raf pvalue     N  GC.Zscore
1  chr7  92221824       rs10 0.033 0.3533 46186 0.92820739
2 chr12 125456933  rs1000000 0.274 0.8824 46186 0.14792748
3  chr4  21227772 rs10000010 0.446 0.2925 46186 1.05265312
4  chr4   1347325 rs10000012 0.854 0.9834 46186 0.02080652
5  chr4  36901464 rs10000013 0.833 0.6183 46186 0.49826113
6  chr4  84997149 rs10000017 0.777 0.4854 46186 0.69764377
> head(HOMA_IR)
    Chr        BP        snp   raf pvalue     N  GC.Zscore
1  chr7  92221824       rs10 0.033 0.8436 46186 0.19729077
2 chr12 125456933  rs1000000 0.274 0.3917 46186 0.85653848
3  chr4  21227772 rs10000010 0.446 0.5571 46186 0.58715442
4  chr4   1347325 rs10000012 0.854 0.9272 46186 0.09136824
5  chr4  36901464 rs10000013 0.833 0.3159 46186 1.00291888
6  chr4  84997149 rs10000017 0.777 0.2266 46186 1.20916360
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=FastGlu$snp,Z.FastGlu=FastGlu$GC.Zscore,key="rs")
dt2 = data.table(rs=FastIns$snp,Z.FastIns=FastIns$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=HOMA_B$snp,Z.HOMA_B=HOMA_B$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=HOMA_IR$snp,Z.HOMA_IR=HOMA_IR$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)

dim(dt3)
head(dt3)

#~~~
#> dim(dt3)
#[1] 2455841       5
#> head(dt3)
#	rs   Z.HOMA_B  Z.HOMA_IR  Z.FastIns  Z.FastGlu
#1:       rs10 0.92820739 0.19729077 0.11833781 0.99219994
#2:  rs1000000 0.14792748 0.85653848 0.60828269 0.66633030
#3: rs10000010 1.05265312 0.58715442 0.87311573 0.08646095
#4: rs10000012 0.02080652 0.09136824 0.08872558 0.79467943
#5: rs10000013 0.49826113 1.00291888 1.04613220 0.38316153
#6: rs10000017 0.69764377 1.20916360 0.80261032 0.37387139
> dim(dt3)
[1] 2455164       5
> head(dt3)
           rs   Z.FastGlu  Z.FastIns   Z.HOMA_B  Z.HOMA_IR
1:       rs10 -0.99219994 0.11833781 0.92820739 0.19729077
2:  rs1000000  0.66633030 0.60828269 0.14792748 0.85653848
3: rs10000010  0.08646095 0.87311573 1.05265312 0.58715442
4: rs10000012 -0.79467943 0.08872558 0.02080652 0.09136824
5: rs10000013  0.38316153 1.04613220 0.49826113 1.00291888
6: rs10000017  0.37387139 0.80261032 0.69764377 1.20916360
#~~~

FastIns.maxZ <- max(dt3$Z.FastIns[!is.infinite(dt3$Z.FastIns) & !is.na(dt3$Z.FastIns)])
FastGlu.maxZ <- max(dt3$Z.FastGlu[!is.infinite(dt3$Z.FastGlu) & !is.na(dt3$Z.FastGlu)])
HOMA_B.maxZ <- max(dt3$Z.HOMA_B[!is.infinite(dt3$Z.HOMA_B) & !is.na(dt3$Z.HOMA_B)])
HOMA_IR.maxZ <- max(dt3$Z.HOMA_IR[!is.infinite(dt3$Z.HOMA_IR) & !is.na(dt3$Z.HOMA_IR)])
#maxZ <- max(c(HOMA_B.maxZ, HOMA_IR.maxZ, FastIns.maxZ, FastGlu.maxZ))
maxZ <- max(c(FastIns.maxZ, FastGlu.maxZ,HOMA_B.maxZ, HOMA_IR.maxZ))

FastIns.maxZ
FastGlu.maxZ
HOMA_B.maxZ
HOMA_IR.maxZ
maxZ

#~~~
#> HOMA_B.maxZ
#[1] 11.14384
#> HOMA_IR.maxZ
#[1] 5.377041
#> FastIns.maxZ
#[1] 5.347434
#> FastGlu.maxZ
#[1] 18.33184
#> maxZ
#[1] 18.33184
> FastIns.maxZ
[1] 5.347434
> FastGlu.maxZ
[1] 15.24669
> HOMA_B.maxZ
[1] 11.14384
> HOMA_IR.maxZ
[1] 5.377041
> maxZ
[1] 15.24669
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt3$Z.FastIns <- apply(as.matrix(dt3$Z.FastIns), 1, replaceInf)
dt3$Z.FastGlu <- apply(as.matrix(dt3$Z.FastGlu), 1, replaceInf)
dt3$Z.HOMA_B <- apply(as.matrix(dt3$Z.HOMA_B), 1, replaceInf)
dt3$Z.HOMA_IR <- apply(as.matrix(dt3$Z.HOMA_IR), 1, replaceInf)

#~~~
#> dt3$Z.HOMA_B <- apply(as.matrix(dt3$Z.HOMA_B), 1, replaceInf)
#> dt3$Z.HOMA_IR <- apply(as.matrix(dt3$Z.HOMA_IR), 1, replaceInf)
#> dt3$Z.FastIns <- apply(as.matrix(dt3$Z.FastIns), 1, replaceInf)
#> dt3$Z.FastGlu <- apply(as.matrix(dt3$Z.FastGlu), 1, replaceInf)
> dt3$Z.FastIns <- apply(as.matrix(dt3$Z.FastIns), 1, replaceInf)
> dt3$Z.FastGlu <- apply(as.matrix(dt3$Z.FastGlu), 1, replaceInf)
> dt3$Z.HOMA_B <- apply(as.matrix(dt3$Z.HOMA_B), 1, replaceInf)
> dt3$Z.HOMA_IR <- apply(as.matrix(dt3$Z.HOMA_IR), 1, replaceInf)
#~~~

attach(dt3)
nullset = (abs(dt3$Z.HOMA_IR)<2) & (abs(dt3$Z.HOMA_B)<2) & (abs(dt3$Z.FastIns)<2) & (abs(dt3$Z.FastGlu)<2) #extract null Z values
Z = cbind(Z.HOMA_IR,Z.HOMA_B,Z.FastIns,Z.FastGlu)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
#> dim(Z)
#[1] 2455841       4
#> dim(Znull)
#[1] 2128367       4
#> RSS0
#	  Z.HOMA_IR   Z.HOMA_B   Z.FastIns   Z.FastGlu
#Z.HOMA_IR 1.00000000 0.35005834 0.824567993 0.061551119
#Z.HOMA_B  0.35005834 1.00000000 0.401291955 0.016342117
#Z.FastIns 0.82456799 0.40129195 1.000000000 0.006789589
#Z.FastGlu 0.06155112 0.01634212 0.006789589 1.000000000
> dim(Z)
[1] 2455164       4
> dim(Znull)
[1] 2127777       4
> RSS0
          Z.HOMA_IR   Z.HOMA_B Z.FastIns  Z.FastGlu
Z.HOMA_IR 1.0000000  0.4323846 0.8252550  0.3417276
Z.HOMA_B  0.4323846  1.0000000 0.4852954 -0.4846437
Z.FastIns 0.8252550  0.4852954 1.0000000  0.1459877
Z.FastGlu 0.3417276 -0.4846437 0.1459877  1.0000000
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
#[1] 359   8
#> dim(dtlesssignif)
#[1] 600   8
#> max(dt3$mvstat)
#[1] 477.2124
#> max(dt3$mvp)
#[1] 101.2459
#> max(dt3$unip)
#[1] 74.3363
#> quantile(dt3$mvp)
#0%          25%          50%          75%         100%
#9.996476e-08 4.492468e-02 1.530084e-01 3.962167e-01 1.012459e+02
#> quantile(dt3$unip)
#0%          25%          50%          75%         100%
#0.008818242  0.384576047  0.643400564  1.021043835 74.336299075
> dim(dtsignif)
[1] 308   8
> dim(dtlesssignif)
[1] 482   8
> max(dt3$mvstat)
[1] 380.048
> max(dt3$mvp)
[1] 80.24529
> max(dt3$unip)
[1] 74.3363
> quantile(dt3$mvp)
          0%          25%          50%          75%         100%
2.401139e-07 6.525237e-02 1.835911e-01 4.298837e-01 8.024529e+01
> quantile(dt3$unip)
          0%          25%          50%          75%         100%
 0.008818242  0.384576047  0.643400564  1.021089423 74.336299075
#~~~

write.table(file="MAGIC2010.dtsignif.vs1.SignCrrct.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="MAGIC2010.dtsignif.rs.vs1.SignCrrct.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="MAGIC2010.dtlesssignif.vs1.SignCrrct.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="MAGIC2010.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="MAGIC2010.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="MAGIC2010.RSS0.vs1.SignCrrct.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


