#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# ICBP2011 data


SBP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_SBP_results_v2.wUCSCGB_snp126.vs3.txt.gz",header=T)
DBP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_DBP_results_v2.wUCSCGB_snp126.vs3.txt.gz",header=T)
PP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_PP_results.wUCSCGB_snp126.vs3.txt.gz",header=T)
MAP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_MAP_results.wUCSCGB_snp126.vs3.txt.gz",header=T)

dim(SBP)
dim(DBP)
dim(PP)
dim(MAP)

#~~~
> dim(SBP)
[1] 2687669      12
> dim(DBP)
[1] 2687669      12
> dim(PP)
[1] 2652060      12
> dim(MAP)
[1] 2651709      12
#~~~

GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }

SBP <- cbind(SBP, apply(as.matrix(SBP$P_value), 2, GetZScore))
DBP <- cbind(DBP, apply(as.matrix(DBP$P_value), 2, GetZScore))
PP <- cbind(PP, apply(as.matrix(PP$P_value), 2, GetZScore))
MAP <- cbind(MAP, apply(as.matrix(MAP$P_value), 2, GetZScore))

SBP.colnames <- colnames(SBP)
SBP.colnames[13] <- "GC.Zscore"
colnames(SBP) <- SBP.colnames
DBP.colnames <- colnames(DBP)
DBP.colnames[13] <- "GC.Zscore"
colnames(DBP) <- DBP.colnames
PP.colnames <- colnames(PP)
PP.colnames[13] <- "GC.Zscore"
colnames(PP) <- PP.colnames
MAP.colnames <- colnames(MAP)
MAP.colnames[13] <- "GC.Zscore"
colnames(MAP) <- MAP.colnames

head(SBP)
head(DBP)
head(PP)
head(MAP)

#~~~
#> head(SBP)
#SNP_ID Phenotype Minor_allele         MAF Sample_size
#1  rs6650104       SBP            G 0.005856237    577.2561
#2 rs10458597       SBP            T 0.001289491   1560.0000
#3 rs12565286       SBP            C 0.057307344  23657.5380
#4 rs11804171       SBP            A 0.056712028  23724.9437
#5 rs12082473       SBP            A 0.002903226   1560.0000
#6 rs12138618       SBP            A 0.051198474  10099.3689
#Number_contrib_studies Estimate_effect Standard_error Coded_target_allele
#1                      2      -0.2820393      2.8949127                   G
#2                      1      -3.1964200      9.6196883                   T
#3                     17       0.3507550      0.3833438                   G
#4                     18       0.3806496      0.3812678                   T
#5                      1       7.0347800      6.4282558                   G
#6                     13      -0.8594717      0.5800714                   G
#Alternative_allele   P_value ChrBP  GC.Zscore
#1                  A 0.9223882  1_554340 0.09742583
#2                  C 0.7396786  1_554484 0.33227896
#3                  C 0.3601980  1_711153 0.91498789
#4                  A 0.3180959  1_713682 0.99837834
#5                  A 0.2738003  1_730720 1.09435285
#6                  A 0.1384294  1_740098 1.48166528
#> head(DBP)
#SNP_ID Phenotype Minor_allele         MAF Sample_size
#1  rs6650104       DBP            G 0.005856237    604.6653
#2 rs10458597       DBP            T 0.001062728   3310.0000
#3 rs12565286       DBP            C 0.057028693  23940.6020
#4 rs11804171       DBP            A 0.057004405  23786.2533
#5 rs12082473       DBP            A 0.001672488   3310.0000
#6 rs12138618       DBP            A 0.054808074   9727.6232
#Number_contrib_studies Estimate_effect Standard_error Coded_target_allele
#1                      2      -1.5539750      1.5508840                   G
#2                      2       2.2636524      4.2507019                   T
#3                     18       0.2701700      0.2351562                   G
#4                     18       0.2667844      0.2353019                   T
#5                      2       7.1110248      3.4368737                   G
#6                     12      -0.3252917      0.3463074                   G
#Alternative_allele    P_value ChrBP GC.Zscore
#1                  A 0.31634695  1_554340 1.0019931
#2                  C 0.59435472  1_554484 0.5325361
#3                  C 0.25059893  1_711153 1.1488959
#4                  A 0.25687989  1_713682 1.1337966
#5                  A 0.03854245  1_730720 2.0690387
#6                  A 0.34756910  1_740098 0.9393149
#> head(PP)
#SNP_ID Phenotype Minor_allele   MAF Sample_size Number_contrib_studies
#1       rs10        PP            A 0.033    30864.25                     46
#2  rs1000000        PP            A 0.373    73344.27                     52
#3 rs10000010        PP            C 0.425    73058.53                     52
#4 rs10000012        PP            G 0.192    72686.58                     52
#5 rs10000013        PP            C 0.167    72806.99                     52
#6 rs10000017        PP            T 0.223    63708.50                     51
#Estimate_effect Standard_error Coded_target_allele Alternative_allele
#1     -0.07658074     0.19174592                   C                  A
#2      0.11621254     0.07117101                   G                  A
#3      0.05949909     0.05897544                   T                  C
#4     -0.20412368     0.08791022                   G                  C
#5      0.03368171     0.07016463                   C                  A
#6      0.02922212     0.07551691                   T                  C
#P_value    ChrBP GC.Zscore
#1 0.68960841   7_92221824 0.3993865
#2 0.10249771 12_125456933 1.6328635
#3 0.31303265   4_21227772 1.0088790
#4 0.02023529    4_1347325 2.3219562
#5 0.63120009   4_36901464 0.4800384
#6 0.69878490   4_84997149 0.3869612
#> head(MAP)
#SNP_ID Phenotype Minor_allele   MAF Sample_size Number_contrib_studies
#1       rs10       MAP            A 0.033    11771.44                     46
#2  rs1000000       MAP            A 0.373    28834.98                     52
#3 rs10000010       MAP            C 0.425    28605.38                     52
#4 rs10000012       MAP            G 0.192    29024.69                     52
#5 rs10000013       MAP            C 0.167    28952.02                     52
#6 rs10000017       MAP            T 0.223    27198.88                     51
#Estimate_effect Standard_error Coded_target_allele Alternative_allele
#1     0.149907976     0.19221654                   C                  A
#2     0.132731566     0.07107920                   G                  A
#3     0.007724790     0.05878420                   T                  C
#4     0.009534625     0.08910274                   G                  C
#5     0.042281407     0.06971416                   C                  A
#6     0.039351143     0.07499834                   T                  C
#P_value    ChrBP GC.Zscore
#1 0.43545495   7_92221824 0.7798911
#2 0.06184916 12_125456933 1.8673756
#3 0.89545154   4_21227772 0.1314093
#4 0.91478335    4_1347325 0.1070071
#5 0.54418501   4_36901464 0.6064966
#6 0.59979624   4_84997149 0.5246936
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=SBP$SNP_ID,Z.SBP=SBP$GC.Zscore,key="rs")
dt2 = data.table(rs=DBP$SNP_ID,Z.DBP=DBP$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=PP$SNP_ID,Z.PP=PP$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=MAP$SNP_ID,Z.MAP=MAP$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)

dim(dt3)
head(dt3)

#~~~
#> dim(dt3)
#[1] 2621375       5
#> head(dt3)
#rs     Z.SBP      Z.DBP      Z.PP     Z.MAP
#1:       rs10 1.1388274 1.53991001 0.3993865 0.7798911
#2:  rs1000000 1.3658297 1.46240099 1.6328635 1.8673756
#3: rs10000010 0.1693811 0.83763902 1.0088790 0.1314093
#4: rs10000012 0.7245045 0.06563208 2.3219562 0.1070071
#5: rs10000013 1.6038426 1.63559284 0.4800384 0.6064966
#6: rs10000017 0.2132554 0.15217732 0.3869612 0.5246936
#~~~

SBP.maxZ <- max(dt3$Z.SBP[!is.infinite(dt3$Z.SBP) & !is.na(dt3$Z.SBP)])
DBP.maxZ <- max(dt3$Z.DBP[!is.infinite(dt3$Z.DBP) & !is.na(dt3$Z.DBP)])
PP.maxZ <- max(dt3$Z.PP[!is.infinite(dt3$Z.PP) & !is.na(dt3$Z.PP)])
MAP.maxZ <- max(dt3$Z.MAP[!is.infinite(dt3$Z.MAP) & !is.na(dt3$Z.MAP)])
maxZ <- max(c(SBP.maxZ, DBP.maxZ, PP.maxZ, MAP.maxZ))

SBP.maxZ
DBP.maxZ
PP.maxZ
MAP.maxZ
maxZ

#~~~
> SBP.maxZ
[1] 7.134293
> DBP.maxZ
[1] 7.675832
> PP.maxZ
[1] 21.09006
> MAP.maxZ
[1] 17.73023
> maxZ
[1] 21.09006
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt3$Z.SBP <- apply(as.matrix(dt3$Z.SBP), 1, replaceInf)
dt3$Z.DBP <- apply(as.matrix(dt3$Z.DBP), 1, replaceInf)
dt3$Z.PP <- apply(as.matrix(dt3$Z.PP), 1, replaceInf)
dt3$Z.MAP <- apply(as.matrix(dt3$Z.MAP), 1, replaceInf)

#~~~
#> dt3$Z.SBP <- apply(as.matrix(dt3$Z.SBP), 1, replaceInf)
#> dt3$Z.DBP <- apply(as.matrix(dt3$Z.DBP), 1, replaceInf)
#> dt3$Z.PP <- apply(as.matrix(dt3$Z.PP), 1, replaceInf)
#> dt3$Z.MAP <- apply(as.matrix(dt3$Z.MAP), 1, replaceInf)
#~~~

dt3 <- dt3[!is.na(dt3$Z.SBP) & !is.na(dt3$Z.DBP) & !is.na(dt3$Z.PP) & !is.na(dt3$Z.MAP),]
#~~~
#> dim(dt3)
#[1] 2621375       5
#> dt3 <- dt3[!is.na(dt3$Z.SBP) & !is.na(dt3$Z.DBP) & !is.na(dt3$Z.PP) & !is.na(dt3$Z.MAP),]
#> dim(dt3)
#[1] 2620113       5
#~~~
attach(dt3)
nullset = (abs(dt3$Z.DBP)<2) & (abs(dt3$Z.SBP)<2) & (abs(dt3$Z.PP)<2) & (abs(dt3$Z.MAP)<2) #extract null Z values
Z = cbind(Z.DBP,Z.SBP,Z.PP,Z.MAP)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
#> dim(Z)
#[1] 2620113       4
#> dim(Znull)
#[1] 2290850       4
#> RSS0
#       Z.DBP     Z.SBP        Z.PP      Z.MAP
#Z.DBP  1.00000000 0.3792907 -0.02820156 0.47014431
#Z.SBP  0.37929066 1.0000000  0.31807313 0.40956136
#Z.PP  -0.02820156 0.3180731  1.00000000 0.06779236
#Z.MAP  0.47014431 0.4095614  0.06779236 1.00000000
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
#[1] 495   8
#> dim(dtlesssignif)
#[1] 925   8
#> max(dt3$mvstat)
#[1] 880.9291
#> max(dt3$mvp)
#[1] 188.6464
#> max(dt3$unip)
#[1] 87.61301
#> quantile(dt3$mvp)
#0%          25%          50%          75%         100%
#7.831873e-10 3.856976e-02 1.350169e-01 3.605996e-01 1.886464e+02
#> quantile(dt3$unip)
#0%          25%          50%          75%         100%
#0.003154806  0.381857807  0.627735001  0.994683962 87.613011631
#~~~

write.table(file="ICBP2011.dtsignif.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="ICBP2011.dtsignif.rs.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="ICBP2011.dtlesssignif.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="ICBP2011.dtlesssignif.rs.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="ICBP2011.dtlesssignif.rs.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="ICBP2011.RSS0.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


