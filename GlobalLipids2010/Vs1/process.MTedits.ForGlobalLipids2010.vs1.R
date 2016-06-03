#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# GlobalLipids2010 data

/data/external_public/GlobalLipids2010
-rw-rw-r-- 1 mturchin20 mturchin20 100516216 Aug 28 23:07 HDL_with_Effect.wHapMap23.expanded.tbl.gz
-rw-rw-r-- 1 mturchin20 mturchin20 100666711 Aug 28 23:11 LDL_with_Effect.wHapMap23.expanded.tbl.gz
-rw-rw-r-- 1 mturchin20 mturchin20 100695962 Aug 28 23:11 TC_with_Effect.wHapMap23.expanded.tbl.gz
-rw-rw-r-- 1 mturchin20 mturchin20 100448320 Aug 28 23:11 TG_with_Effect.wHapMap23.expanded.tbl.gz

HDL=read.table("/data/external_public/GlobalLipids2010/HDL_with_Effect.wHapMap23.expanded.tbl.gz", header=T)
LDL=read.table("/data/external_public/GlobalLipids2010/LDL_with_Effect.wHapMap23.expanded.tbl.gz", header=T)
TC=read.table("/data/external_public/GlobalLipids2010/TC_with_Effect.wHapMap23.expanded.tbl.gz", header=T)
TG=read.table("/data/external_public/GlobalLipids2010/TG_with_Effect.wHapMap23.expanded.tbl.gz", header=T)

dim(HDL)
dim(LDL)
dim(TC)
dim(TG)

#~~~
#> dim(HDL)
#[1] 2692429      14
#> dim(LDL)
#[1] 2692564      14
#> dim(TC)
#[1] 2692413      14
#> dim(TG)
#[1] 2692560      14
#~~~

GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }

HDL <- cbind(HDL, apply(as.matrix(HDL$GC.Pvalue), 2, GetZScore))
LDL <- cbind(LDL, apply(as.matrix(LDL$GC.Pvalue), 2, GetZScore))
TC <- cbind(TC, apply(as.matrix(TC$GC.Pvalue), 2, GetZScore))
TG <- cbind(TG, apply(as.matrix(TG$GC.Pvalue), 2, GetZScore))

HDL.colnames <- colnames(HDL)
HDL.colnames[15] <- "GC.Zscore"
colnames(HDL) <- HDL.colnames
LDL.colnames <- colnames(LDL)
LDL.colnames[15] <- "GC.Zscore"
colnames(LDL) <- LDL.colnames
TC.colnames <- colnames(TC)
TC.colnames[15] <- "GC.Zscore"
colnames(TC) <- TC.colnames
TG.colnames <- colnames(TG)
TG.colnames[15] <- "GC.Zscore"
colnames(TG) <- TG.colnames

head(HDL)
head(LDL)
head(TC)
head(TG)

#~~~
#> head(HDL)
#MarkerName Allele1 Allele2 Weight GC.Zscore GC.Pvalue Overall
#1       rs10       a       c  85853    -0.096   0.92320       -
#2  rs1000000       a       g  99900    -0.937   0.34900       -
#3 rs10000010       t       c  99900     0.185   0.85320       +
#4 rs10000012       c       g  99843    -0.752   0.45230       -
#5 rs10000013       a       c  99900     0.642   0.52090       +
#6 rs10000017       t       c  99900     2.043   0.04106       +
#Direction  Effect StdErr                     ChrBPAFInfo Chr
#1 ---+-+++-+?+--??--++??+++ -0.0019 0.0147   chr7_92221824_A_0.033_C_0.967   7
#2 ---+-+---++-+---+++-++++- -0.0056 0.0058 chr12_125456933_G_0.627_A_0.373  12
#3 -+++-+-+--+++-++++--+----  0.0013 0.0049   chr4_21227772_T_0.575_C_0.425   4
#4 +-+-+-+--+----+--+++++-++ -0.0065 0.0071    chr4_1347325_C_0.808_G_0.192   4
#5 +++-++--+--++--+---+-++--  0.0003 0.0058   chr4_36901464_C_0.167_A_0.833   4
#6 +-++++++---++++---++--+++  0.0105 0.0062   chr4_84997149_C_0.777_T_0.223   4
#BP   MAF  GC.Zscore
#1  92221824 0.033 0.09640364
#2 125456933 0.373 0.93653073
#3  21227772 0.425 0.18503702
#4   1347325 0.192 0.75158611
#5  36901464 0.167 0.64195870
#6  84997149 0.223 2.04292362
#> head(LDL)
#MarkerName Allele1 Allele2 Weight GC.Zscore GC.Pvalue Overall
#1       rs10       a       c  81680     2.051   0.04027       +
#2  rs1000000       a       g  95454     0.662   0.50770       +
#3 rs10000010       t       c  95454     1.583   0.11350       +
#4 rs10000012       c       g  95397     0.155   0.87700       +
#5 rs10000013       a       c  95454    -1.392   0.16400       -
#6 rs10000017       t       c  95454    -0.015   0.98820       -
#Direction  Effect StdErr                     ChrBPAFInfo Chr
#1 ++-+++++-+?++-??-++-??++-  0.0294 0.0152   chr7_92221824_A_0.033_C_0.967   7
#2 -----++-+--+-+++-+++-+-++  0.0044 0.0063 chr12_125456933_G_0.627_A_0.373  12
#3 ++--++---+-+-++--+++---++  0.0073 0.0052   chr4_21227772_T_0.575_C_0.425   4
#4 +-+-----++-++---+++---+-+ -0.0002 0.0076    chr4_1347325_C_0.808_G_0.192   4
#5 --+---+--++--+--++---+-++ -0.0075 0.0062   chr4_36901464_C_0.167_A_0.833   4
#6 --+---+++--+---+-++-++--+  0.0021 0.0066   chr4_84997149_C_0.777_T_0.223   4
#  BP   MAF  GC.Zscore
#  1  92221824 0.033 2.05096865
#  2 125456933 0.373 0.66242326
#  3  21227772 0.425 1.58265553
#  4   1347325 0.192 0.15477335
#  5  36901464 0.167 1.39174378
#  6  84997149 0.223 0.01478965
#> head(TC)
#MarkerName Allele1 Allele2 Weight GC.Zscore GC.Pvalue Overall
#1       rs10       a       c  86136     2.229   0.02582       +
#2  rs1000000       a       g 100184     0.018   0.98550       +
#3 rs10000010       t       c 100184     1.748   0.08050       +
#4 rs10000012       c       g 100125     0.319   0.74980       +
#5 rs10000013       a       c 100184    -1.490   0.13620       -
#6 rs10000017       t       c 100184     0.237   0.81260       +
#Direction  Effect StdErr                     ChrBPAFInfo Chr
#1 +--+++++-+?+++??--+-??-+-  0.0304 0.0151   chr7_92221824_A_0.033_C_0.967   7
#2 ---+-+----++-+++++++-+-++  0.0007 0.0061 chr12_125456933_G_0.627_A_0.373  12
#3 ++-+++-+-+-++++---++---+-  0.0084 0.0051   chr4_21227772_T_0.575_C_0.425   4
#4 +-+----+-+--+---+++--++-+  0.0009 0.0074    chr4_1347325_C_0.808_G_0.192   4
#5 -----+---+--+--+++---+-++ -0.0081 0.0061   chr4_36901464_C_0.167_A_0.833   4
#6 -++---++---+--++-++-----+  0.0027 0.0065   chr4_84997149_C_0.777_T_0.223   4
#BP   MAF  GC.Zscore
#1  92221824 0.033 2.22890840
#2 125456933 0.373 0.01817406
#3  21227772 0.425 1.74779229
#4   1347325 0.192 0.31890309
#5  36901464 0.167 1.49009218
#6  84997149 0.223 0.23707320
#> head(TG)
#MarkerName Allele1 Allele2 Weight GC.Zscore GC.Pvalue Overall
#1       rs10       a       c  82546     0.750    0.4533       +
#2  rs1000000       a       g  96598     1.244    0.2134       +
#3 rs10000010       t       c  96598     0.675    0.4999       +
#4 rs10000012       c       g  96539     0.579    0.5624       +
#5 rs10000013       a       c  96598    -0.533    0.5939       -
#6 rs10000017       t       c  96598     0.290    0.7715       +
#Direction  Effect StdErr                     ChrBPAFInfo Chr
#1 ++------+-?+++??++--??-+-  0.0133 0.0141   chr7_92221824_A_0.033_C_0.967   7
#2 ++--+++++-++-+++----+---+  0.0075 0.0057 chr12_125456933_G_0.627_A_0.373  12
#3 -+-+-----+--+++-+-+++--+-  0.0022 0.0048   chr4_21227772_T_0.575_C_0.425   4
#4 -+-+-+----+-+++-++-+-+++-  0.0032 0.0070    chr4_1347325_C_0.808_G_0.192   4
#5 -+---+++++----+++---+--++ -0.0036 0.0057   chr4_36901464_C_0.167_A_0.833   4
#6 -+------++++-+++++----+++  0.0003 0.0059   chr4_84997149_C_0.777_T_0.223   4
#  BP   MAF GC.Zscore
#  1  92221824 0.033 0.7499248
#  2 125456933 0.373 1.2442706
#  3  21227772 0.425 0.6746471
#  4   1347325 0.192 0.5792804
#  5  36901464 0.167 0.5331930
#  6  84997149 0.223 0.2904134
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=HDL$MarkerName,Z.HDL=HDL$GC.Zscore,key="rs")
dt2 = data.table(rs=LDL$MarkerName,Z.LDL=LDL$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=TC$MarkerName,Z.TC=TC$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=TG$MarkerName,Z.TG=TG$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)

dim(dt3)
head(dt3)

#~~~
#> dim(dt3)
#[1] 2691421       5
#> head(dt3)
#rs  Z.HDL  Z.LDL   Z.TC   Z.TG
#1:       rs10 -0.096  2.051  2.229  0.750
#2:  rs1000000 -0.937  0.662  0.018  1.244
#3: rs10000010  0.185  1.583  1.748  0.675
#4: rs10000012 -0.752  0.155  0.319  0.579
#5: rs10000013  0.642 -1.392 -1.490 -0.533
#6: rs10000017  2.043 -0.015  0.237  0.290
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
#[1] 41.691
#> LDL.maxZ
#[1] 27.854
#> TC.maxZ
#[1] 24.35
#> TG.maxZ
#[1] 24.539
#> maxZ
#[1] 41.691
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt3$Z.HDL <- apply(as.matrix(dt3$Z.HDL), 1, replaceInf)
dt3$Z.LDL <- apply(as.matrix(dt3$Z.LDL), 1, replaceInf)
dt3$Z.TC <- apply(as.matrix(dt3$Z.TC), 1, replaceInf)
dt3$Z.TG <- apply(as.matrix(dt3$Z.TG), 1, replaceInf)

#~~~
#> dt3$Z.HDL <- apply(as.matrix(dt3$Z.HDL), 1, replaceInf)
#> dt3$Z.LDL <- apply(as.matrix(dt3$Z.LDL), 1, replaceInf)
#> dt3$Z.TC <- apply(as.matrix(dt3$Z.TC), 1, replaceInf)
#> dt3$Z.TG <- apply(as.matrix(dt3$Z.TG), 1, replaceInf)
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
#[1] 2691421       4
#> dim(Znull)
#[1] 2284694       4
#> RSS0
#Z.LDL       Z.HDL      Z.TC       Z.TG
#Z.LDL  1.00000000 -0.07931303 0.8393887  0.1815950
#Z.HDL -0.07931303  1.00000000 0.1689808 -0.3668233
#Z.TC   0.83938874  0.16898079 1.0000000  0.3237916
#Z.TG   0.18159499 -0.36682333 0.3237916  1.0000000
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
#[1] 5907    8
#> dim(dtlesssignif)
#[1] 8076    8
#> max(dt3$mvstat)
#[1] 1853.958
#> max(dt3$mvp)
#[1] 399.6142
#> max(dt3$unip)
#[1] 379.1505
#> quantile(dt3$mvp)
#0%          25%          50%          75%         100%
#8.640251e-08 1.031365e-01 2.575558e-01 5.374830e-01 3.996142e+02
#> quantile(dt3$unip)
#0%          25%          50%          75%         100%
#8.396249e-03 4.403058e-01 6.993118e-01 1.076709e+00 3.791505e+02
#~~~

write.table(file="GlobalLipids2010.dtsignif.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GlobalLipids2010.dtsignif.rs.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GlobalLipids2010.dtlesssignif.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GlobalLipids2010.dtlesssignif.rs.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GlobalLipids2010.dtlesssignif.rs.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GlobalLipids2010.RSS0.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


