#this is an incomplete record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# GIANT 2010 data

height=read.table("/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz",header=T)
BMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz", header=T)
WHRadjBMI=read.table("/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz", header=T)

#~~~
#> dim(height)
#[1] 2469635       6
#> dim(BMI)
#[1] 2471516       6
#> dim(WHRadjBMI)
#[1] 2483325       6
#~~~

GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }

height <- cbind(height, apply(as.matrix(height$p), 2, GetZScore))
BMI <- cbind(BMI, apply(as.matrix(BMI$p), 2, GetZScore))
WHRadjBMI <- cbind(WHRadjBMI, apply(as.matrix(WHRadjBMI$p), 2, GetZScore))

height.colnames <- colnames(height)
height.colnames[7] <- "GC.Zscore"
colnames(height) <- height.colnames
BMI.colnames <- colnames(BMI)
BMI.colnames[7] <- "GC.Zscore"
colnames(BMI) <- BMI.colnames
WHRadjBMI.colnames <- colnames(WHRadjBMI)
WHRadjBMI.colnames[7] <- "GC.Zscore"
colnames(WHRadjBMI) <- WHRadjBMI.colnames

#~~~
#> head(height)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU      p      N GC.Zscore
#1       rs10       a       c                 0.0333 0.8826  78380 0.1476741
#2  rs1000000       a       g                 0.3667 0.1858 133822 1.3231064
#3 rs10000010       t       c                  0.575 0.8947 132858 0.1323594
#4 rs10000012       c       g                 0.8083 0.1312 133785 1.5093867
#5 rs10000013       a       c                 0.8333 0.6280 133843 0.4845438
#6 rs10000017       t       c                 0.2333 0.3073 133174 1.0209038
#> head(BMI)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU      p      N GC.Zscore
#1       rs10       a       c                 0.0333 0.7080  80566 0.3745435
#2  rs1000000       g       a                 0.6333 0.5060 123865 0.6650789
#3 rs10000010       c       t                  0.425 0.7360 123827 0.3371551
#4 rs10000012       c       g                 0.8083 0.0420 123809 2.0335201
#5 rs10000013       c       a                 0.1667 0.0689 123863 1.8190749
#6 rs10000017       t       c                 0.2333 0.4570 123262 0.7437958
#> head(WHRadjBMI)
#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU      p     N  GC.Zscore
#1       rs10       c       a                 0.9667 0.4200 57031 0.80642125
#2  rs1000000       g       a                 0.6333 0.5500 77168 0.59776013
#3 rs10000010       t       c                  0.575 0.0029 77152 2.97814368
#4 rs10000012       g       c                 0.1917 0.9900 77117 0.01253347
#5 rs10000013       a       c                 0.8333 0.8900 77167 0.13830421
#6 rs10000017       c       t                 0.7667 0.6400 77166 0.46769880
#> 
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=height$MarkerName,Z.height=height$GC.Zscore,key="rs")
dt2 = data.table(rs=BMI$MarkerName,Z.BMI=BMI$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)

#clear some space 
rm(height)
rm(BMI)
rm(dt2)
rm(dt1)

dt1 = data.table(rs=WHRadjBMI$MarkerName,Z.WHRadjBMI=WHRadjBMI$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)

rm(dt1)
rm(WHRadjBMI)
rm(dt3)

##dt1 = data.table(rs=tg$MarkerName,Z.tg=tg$GC.Zscore,key="rs")
##dt3 = merge(dt1,dt2)
##rm(tg)
##rm(dt1)

#~~~
#> dim(dt2)
#[1] 2436719       4
#~~~

height.maxZ <- max(dt2$Z.height[!is.infinite(dt2$Z.height)])
BMI.maxZ <- max(dt2$Z.BMI[!is.infinite(dt2$Z.BMI)])
WHRadjBMI.maxZ <- max(dt2$Z.WHRadjBMI[!is.infinite(dt2$Z.WHRadjBMI)])
maxZ <- max(c(height.maxZ, BMI.maxZ, WHRadjBMI.maxZ))

#~~~
#> head(maxZ)
#[1] 16.67329
#> maxZ
#[1] 16.67329
#> height.maxZ
#[1] 15.18464
#> BMI.maxZ
#[1] 16.67329
#> WHRadjBMI.maxZ
#[1] 7.773081
#> 
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
Znull_lte1 = Z[nullset_lte1,]
Znull_lte3 = Z[nullset_lte3,]
RSS0 =cor(Znull)
RSS0_lte1=cor(Znull_lte1)
RSS0_lte3=cor(Znull_lte3)
RSS0inv = chol2inv(chol(RSS0))
RSS0inv_lte3 = chol2inv(chol(RSS0_lte3))

#~~~
#> dim(Z)
#[1] 2436719       3
#> dim(nullset)
#NULL
#> length(nullset)
#[1] 2436719
#> dim(Znull)
#[1] 2041036       3
#> 
#> dim(Znull_lte1)
#[1] 755502      3
#> dim(Znull_lte3)
#[1] 2368435       3
#>
#> RSS0
#Z.BMI    Z.height  Z.WHRadjBMI
#Z.BMI       1.0000000000 0.007383229 0.0003445489
#Z.height    0.0073832291 1.000000000 0.0093910411
#Z.WHRadjBMI 0.0003445489 0.009391041 1.0000000000
#> RSS0_lte1
#Z.BMI     Z.height  Z.WHRadjBMI
#Z.BMI        1.000000000 -0.005718059 -0.005199878
#Z.height    -0.005718059  1.000000000 -0.002017754
#Z.WHRadjBMI -0.005199878 -0.002017754  1.000000000
#> RSS0_lte3
#Z.BMI   Z.height Z.WHRadjBMI
#Z.BMI       1.000000000 0.02502041 0.002067586
#Z.height    0.025020411 1.00000000 0.019565570
#Z.WHRadjBMI 0.002067586 0.01956557 1.000000000
#> RSS0_All
#Z.BMI   Z.height Z.WHRadjBMI
#Z.BMI       1.00000000 0.05237274  0.01426970
#Z.height    0.05237274 1.00000000  0.04609061
#Z.WHRadjBMI 0.01426970 0.04609061  1.00000000
#> 
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

write.table(file="GIANT2010.dtsignif.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.dtsignif.rs.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.dtlesssignif.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.dtlesssignif.rs.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="GIANT2010.dtlesssignif.rs.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GIANT2010.RSS0.vs1.txt",RSS0,row.names=F,sep=",",quote=F)

#~~~
#> dim(dtsignif)
#[1] 5882    7
#> dim(dtlesssignif)
#[1] 8841    7
#>
#> max(dt2$mvstat)
#[1] 284.9248
#> max(dt2$mvp)
#[1] 59.71389
#> max(dt2$unip)
#[1] 61.68825
#> quantile(dt2$mvp)
#0%          25%          50%          75%         100% 
#3.429455e-11 5.884256e-02 1.835219e-01 4.478434e-01 5.971389e+01 
#> quantile(dt2$unip)
#0%          25%          50%          75%         100% 
#0.001740662  0.438898616  0.703334810  1.099414413 61.688246139 
#> 
#~~~

#Doing lte3 RSS0 cutoffs don't seem to make much of a difference
mvstat_lte3 =  rowSums(Z * (Z %*% RSS0inv_lte3)) # comptues Z RSS0^-1 Z'
dt2$mvstat_lte3 = mvstat_lte3
statchi_lte3 = -log10(exp(1))*pchisq(mvstat_lte3,df=4,log.p=TRUE,lower.tail=FALSE)
dt2$mvp_lte3 = statchi_lte3

dtsignif_lte3 = dt2[dt2$mvp_lte3>-log10(5e-8) | dt2$unip>-log10(5e-8),]
dtlesssignif_lte3 = dt2[dt2$mvp_lte3>-log10(1e-6) | dt2$unip>-log10(1e-6),]

write.table(file="GIANT2010.dtsignif_lte3.vs1.txt",dtsignif_lte3,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.dtsignif_lte3.rs.vs1.txt",dtsignif_lte3$rs,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.dtlesssignif_lte3.vs1.txt",dtlesssignif_lte3,sep=",",row.names=F,quote=F)
write.table(file="GIANT2010.dtlesssignif_lte3.rs.vs1.txt",dtlesssignif_lte3$rs,sep=",",row.names=F,quote=F)
#write.table(file="GIANT2010.dtlesssignif_lte3.rs.vs1.txt",paste(dtlesssignif_lte3$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="GIANT2010.RSS0_lte3.vs1.txt",RSS0_lte3,row.names=F,sep=",",quote=F)

#~~~
#> dim(dtsignif_lte3)
#[1] 5851    9
#> dim(dtlesssignif_lte3)
#[1] 8789    9
#> max(dt2$mvp_lte3)
#[1] 59.47925
#> quantile(dt2$mvp_lte3)
#0%          25%          50%          75%         100% 
#3.433395e-11 5.764963e-02 1.800992e-01 4.400740e-01 5.947925e+01 	  
#~~~



rm(dt2)

