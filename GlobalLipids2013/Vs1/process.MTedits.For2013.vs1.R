#this is an incomplete record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# global lipids data

tg=read.table("/data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz",header=T)
tc=read.table("/data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz",header=T)
ldl=read.table("/data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz",header=T)
hdl=read.table("/data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz",header=T)

tg <- cbind(tg, apply(as.matrix(tg$P.value), 2, qnorm))
tc <- cbind(tc, apply(as.matrix(tc$P.value), 2, qnorm))
ldl <- cbind(ldl, apply(as.matrix(ldl$P.value), 2, qnorm))
hdl <- cbind(hdl, apply(as.matrix(hdl$P.value), 2, qnorm))

tg.colnames <- colnames(tg)
tg.colnames[11] <- "GC.Zscore"
tg.colnames[1] <- "MarkerName"
colnames(tg) <- tg.colnames
tc.colnames <- colnames(tc)
tc.colnames[11] <- "GC.Zscore"
tc.colnames[1] <- "MarkerName"
colnames(tc) <- tc.colnames
ldl.colnames <- colnames(ldl)
ldl.colnames[11] <- "GC.Zscore"
ldl.colnames[1] <- "MarkerName"
colnames(ldl) <- ldl.colnames
hdl.colnames <- colnames(hdl)
hdl.colnames[11] <- "GC.Zscore"
hdl.colnames[1] <- "MarkerName"
colnames(hdl) <- hdl.colnames

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=ldl$MarkerName,Z.ldl=ldl$GC.Zscore,key="rs")
dt2 = data.table(rs=hdl$MarkerName,Z.hdl=hdl$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)

#clear some space 
rm(ldl)
rm(hdl)
rm(dt2)
rm(dt1)

dt1 = data.table(rs=tc$MarkerName,Z.tc=tc$GC.Zscore,key="rs")
dt2 = merge(dt1,dt3)

rm(dt1)
rm(tc)
rm(dt3)

dt1 = data.table(rs=tg$MarkerName,Z.tg=tg$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(tg)
rm(dt1)

attach(dt3)
nullset = (abs(dt3$Z.tc)<2) & (abs(dt3$Z.tg)<2) & (abs(dt3$Z.hdl)<2) & (abs(dt3$Z.ldl)<2) #extract null Z values
Z = cbind(Z.tc,Z.tg,Z.hdl,Z.ldl)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

mvstat =  rowSums(Z * (Z %*% RSS0inv)) # comptues Z RSS0^-1 Z'
dt3$mvstat = mvstat
statchi = -log10(exp(1))*pchisq(mvstat,df=4,log.p=TRUE,lower.tail=FALSE)
dt3$mvp = statchi

maxZ2 = apply(Z^2,1,max)
max.unip = -log10(exp(1))*pchisq(maxZ2,df=1,log.p=TRUE, lower.tail=FALSE)
dt3$unip = max.unip

dtsignif = dt3[dt3$mvp>-log10(5e-8) | dt3$unip>-log10(5e-8),]
dtlesssignif = dt3[dt3$mvp>-log10(1e-6) | dt3$unip>-log10(1e-6),]


write.table(file="dtsignif.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="dtsignif.rs.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="dtlesssignif.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="dtlesssignif.rs.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
write.table(file="dtlesssignif.rs.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="RSS0.txt",RSS0,row.names=F,sep=",",quote=F)

rm(dt3)

