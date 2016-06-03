#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# ICBP2011 data


SBP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_SBP_results_v2.wUCSCGB_snp126.vs3.IncAllele.txt.gz",header=T)
DBP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_DBP_results_v2.wUCSCGB_snp126.vs3.MatchedToSBP.IncAllele.txt.gz",header=F)
PP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_PP_results.wUCSCGB_snp126.vs3.MatchedToSBP.IncAllele.txt.gz",header=F)
MAP=read.table("/mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_MAP_results.wUCSCGB_snp126.vs3.MatchedToSBP.IncAllele.txt.gz",header=F)

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
> dim(SBP)
[1] 2687669      12
> dim(DBP)
[1] 2687669      12
> dim(PP)
[1] 2610752      12
> dim(MAP)
[1] 2610640      12

#~~~


#MAP2 <- MAP

colnames(DBP) <- colnames(SBP)
colnames(PP) <- colnames(SBP)
colnames(MAP) <- colnames(SBP)

#SBP2 <- SBP
SBP <- SBP[!is.na(SBP$P_value),]
DBP <- DBP[!is.na(DBP$P_value),]
PP <- PP[!is.na(PP$P_value),]
MAP <- MAP[!is.na(MAP$P_value),]

~~~
> table(is.na(as.matrix(as.numeric(as.character(DBP$P_value)))))

  FALSE    TRUE 
2687152       1 
Warning message:
In as.matrix(as.numeric(as.character(DBP$P_value))) :
  NAs introduced by coercion
> table(is.na(as.matrix(as.numeric(as.character(PP$P_value)))))

  FALSE    TRUE 
2610752       1 
Warning message:
In as.matrix(as.numeric(as.character(PP$P_value))) :
  NAs introduced by coercion
> table(is.na(as.matrix(as.numeric(as.character(MAP$P_value)))))

  FALSE    TRUE 
2610640       1 
Warning message:
In as.matrix(as.numeric(as.character(MAP$P_value))) :
  NAs introduced by coercion
> MAP[is.na(as.matrix(as.numeric(as.character(MAP$P_value)))),]
        SNP_ID Phenotype Minor_allele X1 Sample_size Number_contrib_studies
2610641 SNP_ID Phenotype Minor_allele  1 Sample_size Number_contrib_studies
        Estimate_effect Standard_error Coded_target_allele Alternative_allele
2610641 Estimate_effect Standard_error Coded_target_allele Alternative_allele
        P_value ChrBP
2610641 P_value ChrBP
Warning message:
In as.matrix(as.numeric(as.character(MAP$P_value))) :
  NAs introduced by coercion
~~~

DBP <- DBP[!is.na(as.numeric(as.character(DBP$P_value))),]
DBP$P_value <- as.numeric(as.character(DBP$P_value))
PP <- PP[!is.na(as.numeric(as.character(PP$P_value))),]
PP$P_value <- as.numeric(as.character(PP$P_value))
MAP <- MAP[!is.na(as.numeric(as.character(MAP$P_value))),]
MAP$P_value <- as.numeric(as.character(MAP$P_value))

~~~
> log(0)
[1] -Inf
> table(DBP$P_value==0)

  FALSE    TRUE 
2686603     549 
~~~

DBP <- DBP[DBP$P_value!=0,]

#GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }
GetZScoreAndSign <- function(x) { val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE); if (x < 0) { val1 <- val1 * -1 }; return(val1) }

#height <- cbind(height, apply(as.matrix(height$p), 1, GetZScoreAndSign))
#BMI <- cbind(BMI, apply(as.matrix(BMI$p), 1, GetZScoreAndSign))


DBP <- cbind(DBP, apply(as.matrix(DBP$P_value), 1, GetZScoreAndSign))
SBP <- cbind(SBP, apply(as.matrix(SBP$P_value), 1, GetZScoreAndSign))
PP <- cbind(PP, apply(as.matrix(PP$P_value), 1, GetZScoreAndSign))
MAP <- cbind(MAP, apply(as.matrix(MAP$P_value), 1, GetZScoreAndSign))

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
> head(SBP)
      SNP_ID Phenotype Minor_allele         X1 Sample_size
1  rs6650104       SBP            G 0.99414376    577.2561
2 rs10458597       SBP            T 0.99871051   1560.0000
3 rs12565286       SBP            C 0.94269266  23657.5380
4 rs11804171       SBP            A 0.94328797  23724.9437
5 rs12082473       SBP            A 0.99709677   1560.0000
6 rs12138618       SBP            A 0.05119847  10099.3689
  Number_contrib_studies Estimate_effect Standard_error Coded_target_allele
1                      2       0.2820393      2.8949127                   A
2                      1       3.1964200      9.6196883                   C
3                     17       0.3507550      0.3833438                   G
4                     18       0.3806496      0.3812678                   T
5                      1       7.0347800      6.4282558                   G
6                     13       0.8594717      0.5800714                   A
  Alternative_allele   P_value    ChrBP  GC.Zscore
1                  G 0.9223882 1_554340 0.09742583
2                  T 0.7396786 1_554484 0.33227896
3                  C 0.3601980 1_711153 0.91498789
4                  A 0.3180959 1_713682 0.99837834
5                  A 0.2738003 1_730720 1.09435285
6                  G 0.1384294 1_740098 1.48166528
> head(DBP)
      SNP_ID Phenotype Minor_allele        X1      Sample_size
1  rs1000000       DBP            A 0.7815214 69219.4332910725
2 rs10000010       DBP            T 0.5036087 68814.1423335478
3 rs10000012       DBP            G 0.1333223 68759.8338271359
4 rs10000013       DBP            C 0.2285297  68793.109689461
5 rs10000017       DBP            T 0.7783245 60965.6066681721
6 rs10000023       DBP            G 0.4053071     67077.015605
  Number_contrib_studies     Estimate_effect     Standard_error
1                     43   0.108183801043974 0.0739768379628118
2                     43  0.0515364405538518 0.0615258353192515
3                     43 0.00594184619683698 0.0905326449937789
4                     43   0.119900283366211 0.0733069261123304
5                     43  0.0119901592618044 0.0787907095443636
6                     43   0.076234538743382  0.063051677897002
  Coded_target_allele Alternative_allele    P_value        ChrBP   GC.Zscore
1                   G                  A  0.1436314 12_125456933  1.46240099
2                   C                  T  0.4022335   4_21227772  0.83763902
3                   G                  C -0.9476707    4_1347325 -0.06563208
4                   C                  A  0.1019248   4_36901464  1.63559284
5                   C                  T -0.8790471   4_84997149 -0.15217732
6                   G                  T  0.2266320   4_95952929  1.20908025
> head(PP)
he      SNP_ID Phenotype Minor_allele    X1      Sample_size
1  rs1000000        PP            A 0.627 73344.2702230173
2 rs10000010        PP            C 0.575  73058.533850557
3 rs10000012        PP            G 0.808 72686.5791272688
4 rs10000013        PP            C 0.167 72806.9946511816
5 rs10000017        PP            T 0.223 63708.4950689534
6 rs10000023        PP            G 0.408     71437.234345
  Number_contrib_studies    Estimate_effect     Standard_error
1                     52  0.116212543195141 0.0711710088206335
2                     52  0.059499089797371  0.058975444263077
3                     52  0.204123677118336 0.0879102178534104
4                     52 0.0336817146362768 0.0701646263214127
5                     51 0.0292221157200947 0.0755169053634177
6                     52 0.0092815158185631 0.0606664262151885
  Coded_target_allele Alternative_allele     P_value        ChrBP  GC.Zscore
1                   G                  A  0.10249771 12_125456933  1.6328635
2                   T                  C -0.31303265   4_21227772 -1.0088790
3                   C                  G  0.02023529    4_1347325  2.3219562
4                   C                  A  0.63120009   4_36901464  0.4800384
5                   T                  C  0.69878490   4_84997149  0.3869612
6                   G                  T  0.87840409   4_95952929  0.1529926
> head(MAP)
      SNP_ID Phenotype Minor_allele    X1   Sample_size Number_contrib_studies
1  rs1000000       MAP            A 0.627  28834.981022                     52
2 rs10000010       MAP            C 0.575 28605.3784795                     52
3 rs10000012       MAP            G 0.192  29024.693561                     52
4 rs10000013       MAP            C 0.167  28952.015061                     52
5 rs10000017       MAP            T 0.223 27198.8757635                     51
6 rs10000023       MAP            G 0.408  27883.396734                     52
      Estimate_effect     Standard_error Coded_target_allele Alternative_allele
1   0.132731566188609   0.07107920224497                   G                  A
2 0.00772478995450941 0.0587842002547057                   T                  C
3  0.0095346253630703 0.0891027368922516                   G                  C
4  0.0422814069578563 0.0697141644153772                   C                  A
5  0.0393511430399704 0.0749983351265776                   T                  C
6  0.0618461423651648 0.0604179996800068                   G                  T
      P_value        ChrBP  GC.Zscore
1  0.06184916 12_125456933  1.8673756
2 -0.89545154   4_21227772 -0.1314093
3 -0.91478335    4_1347325 -0.1070071
4  0.54418501   4_36901464  0.6064966
5  0.59979624   4_84997149  0.5246936
6  0.30600643   4_95952929  1.0236377
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
> dim(dt3)
[1] 2609624       5
> head(dt3)
           rs     Z.SBP       Z.DBP       Z.PP      Z.MAP
1:       rs10 1.1388274  1.53991001 -0.3993865  0.7798911
2:  rs1000000 1.3658297  1.46240099  1.6328635  1.8673756
3: rs10000010 0.1693811  0.83763902 -1.0088790 -0.1314093
4: rs10000012 0.7245045 -0.06563208  2.3219562 -0.1070071
5: rs10000013 1.6038426  1.63559284  0.4800384  0.6064966
6: rs10000017 0.2132554 -0.15217732  0.3869612  0.5246936
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
> dt3$Z.PP <- apply(as.matrix(dt3$Z.PP), 1, replaceInf)
.
.
.
#After fixing P_value == 0 problem above
> table(is.infinite(dt3$Z.DBP))

  FALSE 
2609310 
.
.
>
> dt3$Z.SBP <- apply(as.matrix(dt3$Z.SBP), 1, replaceInf)
> dt3$Z.DBP <- apply(as.matrix(dt3$Z.DBP), 1, replaceInf)
> dt3$Z.PP <- apply(as.matrix(dt3$Z.PP), 1, replaceInf)
> dt3$Z.MAP <- apply(as.matrix(dt3$Z.MAP), 1, replaceInf)
>
#~~~

dt3 <- dt3[!is.na(dt3$Z.SBP) & !is.na(dt3$Z.DBP) & !is.na(dt3$Z.PP) & !is.na(dt3$Z.MAP),]

#~~~
#> dim(dt3)
#[1] 2621375       5
#> dt3 <- dt3[!is.na(dt3$Z.SBP) & !is.na(dt3$Z.DBP) & !is.na(dt3$Z.PP) & !is.na(dt3$Z.MAP),]
#> dim(dt3)
#[1] 2620113       5
> dim(dt3)
[1] 2609310       5
> dt3 <- dt3[!is.na(dt3$Z.SBP) & !is.na(dt3$Z.DBP) & !is.na(dt3$Z.PP) & !is.na(dt3$Z.MAP),]
> 
> dim(dt3)
[1] 2609310       5
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
> RSS0
           Z.DBP     Z.SBP        Z.PP      Z.MAP
Z.DBP  1.0000000 0.4699008 -0.30990002 0.62217302
Z.SBP  0.4699008 1.0000000  0.42686893 0.48124375
Z.PP  -0.3099000 0.4268689  1.00000000 0.04974055
Z.MAP  0.6221730 0.4812437  0.04974055 1.00000000
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
> dim(dtsignif)
[1] 522   8
> dim(dtlesssignif)
[1] 1235    8
> max(dt3$mvstat)
[1] 1158.778
> max(dt3$mvp)
[1] 248.8617
> max(dt3$unip)
[1] 87.61301
> quantile(dt3$mvp)
          0%          25%          50%          75%         100% 
3.384360e-09 6.237219e-02 1.760454e-01 4.157306e-01 2.488617e+02 
> quantile(dt3$unip)
          0%          25%          50%          75%         100% 
 0.003154806  0.381811571  0.627688220  0.994596173 87.613011631 
> 
#~~~

write.table(file="ICBP2011.dtsignif.vs1.SignCrrct.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="ICBP2011.dtsignif.rs.vs1.SignCrrct.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="ICBP2011.dtlesssignif.vs1.SignCrrct.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="ICBP2011.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="ICBP2011.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="ICBP2011.RSS0.vs1.SignCrrct.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


