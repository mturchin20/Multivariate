#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# GEFOS2015 data

fa=read.table("/mnt/gluster/data/external_public_supp/GEFOS2015/fa2stu.MAF0_.005.pos_.maf05.IncAllele.formatted.out_.gz", header=T)
fn=read.table("/mnt/gluster/data/external_public_supp/GEFOS2015/fn2stu.MAF0_.005.pos_.maf05.MatchedTofa.IncAllele.formatted.out_.gz", header=T)
ls=read.table("/mnt/gluster/data/external_public_supp/GEFOS2015/ls2stu.MAF0_.005.pos_.maf05.MatchedTofa.IncAllele.formatted.out_.gz", header=T)

dim(fa)
dim(fn)
dim(ls)

#~~~
> dim(fa)
[1] 6050596       6
> dim(fn)
[1] 5930898       6
> dim(ls)
[1] 5929684       6
#~~~

#GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }
GetZScoreAndSign <- function(x) { val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE); if (x < 0) { val1 <- val1 * -1 }; return(val1) }

fa <- cbind(fa, apply(as.matrix(fa$p.value), 1, GetZScoreAndSign))
fn <- cbind(fn, apply(as.matrix(fn$p.value), 1, GetZScoreAndSign))
ls <- cbind(ls, apply(as.matrix(ls$p.value), 1, GetZScoreAndSign))

fa.colnames <- colnames(fa)
fa.colnames[7] <- "GC.Zscore"
colnames(fa) <- fa.colnames
fn.colnames <- colnames(fn)
fn.colnames[7] <- "GC.Zscore"
colnames(fn) <- fn.colnames
ls.colnames <- colnames(ls)
ls.colnames[7] <- "GC.Zscore"
colnames(ls) <- ls.colnames

head(fa)
head(fn)
head(ls)

#~~~
> head(fa)
  chromosome position     rs_number      eaf  p.value X32965 GC.Zscore
1          1 12776218 chr1:12776218 0.133095 0.251993  32965 1.1455220
2          1 12785749 chr1:12785749 0.931742 0.732971  32965 0.3411761
3          1 12786323 chr1:12786323 0.096977 0.202654  32965 1.2740266
4          1 12788113 chr1:12788113 0.097136 0.208868  32965 1.2566842
5          1 12722028 chr1:12722028 0.088465 0.016927  32965 2.3882895
6          1 12725943 chr1:12725943 0.127844 0.817179  32965 0.2311748
> head(fn)
  chromosome  position       rs_number      eaf   p.value X32965  GC.Zscore
1         10 100015563 chr10:100015563 0.865314 -0.094948  32965 -1.6698553
2         10 100019039 chr10:100019039 0.871679 -0.166509  32965 -1.3835085
3         10 100020880 chr10:100020880 0.123707 -0.778256  32965 -0.2815925
4         10 100021672 chr10:100021672 0.907478  0.028682  32965  2.1878298
5         10 100021983 chr10:100021983 0.127924 -0.772906  32965 -0.2885758
6         10 100031394 chr10:100031394 0.446858  0.009232  32965  2.6033413
> head(ls)
  chromosome  position       rs_number      eaf   p.value X32965  GC.Zscore
1         10 100015563 chr10:100015563 0.865314 -0.357584  32965 -0.9199785
2         10 100019039 chr10:100019039 0.871679 -0.453859  32965 -0.7489970
3         10 100020880 chr10:100020880 0.123707 -0.662517  32965 -0.4364407
4         10 100021672 chr10:100021672 0.907478  0.126808  32965  1.5268111
5         10 100021983 chr10:100021983 0.127924 -0.529257  32965 -0.6291406
6         10 100031394 chr10:100031394 0.446858  0.482571  32965  0.7021735
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=fa$rs_number,Z.fa=fa$GC.Zscore,key="rs")
dt2 = data.table(rs=fn$rs_number,Z.fn=fn$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=ls$rs_number,Z.ls=ls$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)

dim(dt2)
head(dt2)

#~~~
> dim(dt2)
[1] 5929102       4
> head(dt2)
                rs      Z.fa       Z.fn       Z.ls
1: chr10:100015563 0.3929045 -1.6698553 -0.9199785
2: chr10:100019039 0.6172922 -1.3835085 -0.7489970
3: chr10:100020880 0.8892793 -0.2815925 -0.4364407
4: chr10:100021672 0.1095163  2.1878298  1.5268111
5: chr10:100021983 0.9433498 -0.2885758 -0.6291406
6: chr10:100031394 0.8024495  2.6033413  0.7021735
#~~~

fa.maxZ <- max(dt2$Z.fa[!is.infinite(dt2$Z.fa) & !is.na(dt2$Z.fa)])
fn.maxZ <- max(dt2$Z.fn[!is.infinite(dt2$Z.fn) & !is.na(dt2$Z.fn)])
ls.maxZ <- max(dt2$Z.ls[!is.infinite(dt2$Z.ls) & !is.na(dt2$Z.ls)])
maxZ <- max(c(fa.maxZ, fn.maxZ, ls.maxZ))

fa.maxZ
fn.maxZ
ls.maxZ
maxZ

#~~~
> fa.maxZ
[1] 10.46814
> fn.maxZ
[1] 10.32486
> ls.maxZ
[1] 9.226354
> maxZ
[1] 10.46814
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt2$Z.fa <- apply(as.matrix(dt2$Z.fa), 1, replaceInf)
dt2$Z.fn <- apply(as.matrix(dt2$Z.fn), 1, replaceInf)
dt2$Z.ls <- apply(as.matrix(dt2$Z.ls), 1, replaceInf)

#~~~
> 
> dt2$Z.fa <- apply(as.matrix(dt2$Z.fa), 1, replaceInf)
> dt2$Z.fn <- apply(as.matrix(dt2$Z.fn), 1, replaceInf)
> dt2$Z.ls <- apply(as.matrix(dt2$Z.ls), 1, replaceInf)
> 
#~~~

attach(dt2)
nullset = (abs(dt2$Z.fn)<2) & (abs(dt2$Z.fa)<2) & (abs(dt2$Z.ls)<2) #extract null Z values
Z = cbind(Z.fn,Z.fa,Z.ls)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
> dim(Z)
[1] 5929102       3
> dim(Znull)
[1] 5128015       3
> RSS0
           Z.fn       Z.fa      Z.ls
Z.fn 1.00000000 0.09853622 0.4133810
Z.fa 0.09853622 1.00000000 0.1111869
Z.ls 0.41338096 0.11118691 1.0000000
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
> dim(dtsignif)
[1] 1969    7
> dim(dtlesssignif)
[1] 2866    7
> max(dt2$mvstat)
[1] 155.1756
> max(dt2$mvp)
[1] 32.69582
> max(dt2$unip)
[1] 24.91721
> quantile(dt2$mvp)
          0%          25%          50%          75%         100% 
8.280529e-09 1.203456e-01 2.916630e-01 5.924899e-01 3.269582e+01 
> quantile(dt2$unip)
          0%          25%          50%          75%         100% 
8.224635e-04 4.141334e-01 6.672768e-01 1.033464e+00 2.491721e+01 
#~~~

write.table(file="GEFOS2015.dtsignif.vs1.SignCrrct.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GEFOS2015.dtsignif.rs.vs1.SignCrrct.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GEFOS2015.dtlesssignif.vs1.SignCrrct.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GEFOS2015.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GEFOS2015.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GEFOS2015.RSS0.vs1.SignCrrct.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


