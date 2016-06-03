#this is a record of how the files RSS0.txt and dtlesssignif.annot.txt were produced
# HaemgenRBC2012 data

/mnt/gluster/data/external_public_supp/HaemgenRBC2012/
-rw-rw-r-- 1 mturchin20 mturchin20 70167800 Aug 14 16:54 HaemGenRBC_RBC.wHapMap22.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 73328129 Aug 14 17:00 HaemGenRBC_MCV.wHapMap22.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 72731629 Aug 14 17:00 HaemGenRBC_PCV.wHapMap22.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 71907633 Aug 14 17:00 HaemGenRBC_MCH.wHapMap22.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 71654461 Aug 14 17:01 HaemGenRBC_Hb.wHapMap22.txt.gz
-rw-rw-r-- 1 mturchin20 mturchin20 71271392 Aug 14 17:01 HaemGenRBC_MCHC.wHapMap22.txt.gz


RBC=read.table("/mnt/gluster/data/external_public_supp/HaemgenRBC2012/HaemGenRBC_RBC.wHapMap22.IncAllele.txt.gz",header=T)
MCV=read.table("/mnt/gluster/data/external_public_supp/HaemgenRBC2012/HaemGenRBC_MCV.wHapMap22.MatchedToRBC.IncAllele.txt.gz",header=T)
PCV=read.table("/mnt/gluster/data/external_public_supp/HaemgenRBC2012/HaemGenRBC_PCV.wHapMap22.MatchedToRBC.IncAllele.txt.gz",header=T)
MCH=read.table("/mnt/gluster/data/external_public_supp/HaemgenRBC2012/HaemGenRBC_MCH.wHapMap22.MatchedToRBC.IncAllele.txt.gz",header=T)
Hb=read.table("/mnt/gluster/data/external_public_supp/HaemgenRBC2012/HaemGenRBC_Hb.wHapMap22.MatchedToRBC.IncAllele.txt.gz",header=T)
MCHC=read.table("/mnt/gluster/data/external_public_supp/HaemgenRBC2012/HaemGenRBC_MCHC.wHapMap22.MatchedToRBC.IncAllele.txt.gz",header=T)

dim(RBC)
dim(MCV)
dim(PCV)
dim(MCH)
dim(Hb)
dim(MCHC)

#~~~
#> dim(RBC)
#[1] 2589454      10
#> dim(MCV)
#[1] 2591132      10
#> dim(PCV)
#[1] 2591079      10
#> dim(MCH)
#[1] 2586784      10
#> dim(Hb)
#[1] 2593078      10
#> dim(MCHC)
#[1] 2588875      10
> dim(RBC)
[1] 2462918      15
> dim(MCV)
[1] 2259247      15
> dim(PCV)
[1] 2259116      15
> dim(MCH)
[1] 2257229      15
> dim(Hb)
[1] 2259770      15
> dim(MCHC)
[1] 2257386      15
#~~~

#GetZScore <- function(x) { val1 <- qnorm(log(x/2), lower.tail=FALSE, log.p=TRUE); return(val1) }
GetZScoreAndSign <- function(x) { val1 <- qnorm(log(abs(x)/2), lower.tail=FALSE, log.p=TRUE); if (x < 0) { val1 <- val1 * -1 }; return(val1) }

RBC <- cbind(RBC, apply(as.matrix(RBC$Pval), 1, GetZScoreAndSign))
MCV <- cbind(MCV, apply(as.matrix(MCV$Pval), 1, GetZScoreAndSign))
PCV <- cbind(PCV, apply(as.matrix(PCV$Pval), 1, GetZScoreAndSign))
MCH <- cbind(MCH, apply(as.matrix(MCH$Pval), 1, GetZScoreAndSign))
Hb <- cbind(Hb, apply(as.matrix(Hb$Pval), 1, GetZScoreAndSign))
MCHC <- cbind(MCHC, apply(as.matrix(MCHC$Pval), 1, GetZScoreAndSign))

RBC.colnames <- colnames(RBC)
RBC.colnames[16] <- "GC.Zscore"
colnames(RBC) <- RBC.colnames
MCV.colnames <- colnames(MCV)
MCV.colnames[16] <- "GC.Zscore"
colnames(MCV) <- MCV.colnames
PCV.colnames <- colnames(PCV)
PCV.colnames[16] <- "GC.Zscore"
colnames(PCV) <- PCV.colnames
MCH.colnames <- colnames(MCH)
MCH.colnames[16] <- "GC.Zscore"
colnames(MCH) <- MCH.colnames
Hb.colnames <- colnames(Hb)
Hb.colnames[16] <- "GC.Zscore"
colnames(Hb) <- Hb.colnames
MCHC.colnames <- colnames(MCHC)
MCHC.colnames[16] <- "GC.Zscore"
colnames(MCHC) <- MCHC.colnames

head(RBC)
head(MCV)
head(PCV)
head(MCH)
head(Hb)
head(MCHC)

#~~~
#> head(RBC)
#SNP     Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
#1  rs9373124 5.11e-97   6 135464902             T                 C       53337
#2  rs7776054 3.51e-96   6 135460609             A                 G       53403
#3  rs9389269 4.54e-96   6 135468852             T                 C       53402
#4  rs9389268 5.98e-95   6 135461324             A                 G       53403
#5  rs9376090 9.04e-95   6 135452921             T                 C       53537
#6 rs11759553 1.84e-94   6 135463989             A                 T       53398
#Beta     SE                    ChrBPAFInfo GC.Zscore
#1 0.0484 0.0024   chr6_135464902_T_0.78_C_0.22  20.90223
#2 0.0479 0.0024   chr6_135460609_A_0.78_G_0.22  20.81005
#3 0.0469 0.0023   chr6_135468852_T_0.78_C_0.22  20.79771
#4 0.0471 0.0024 chr6_135461324_A_0.775_G_0.225  20.67366
#5 0.0477 0.0024   chr6_135452921_T_0.78_C_0.22  20.65371
#6 0.0465 0.0023 chr6_135463989_A_0.775_T_0.225  20.61935
#> head(MCV)
#SNP      Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
#1 rs9389269 2.57e-109   6 135468852             T                 C       57855
#2 rs7776054 3.13e-108   6 135460609             A                 G       57856
#3 rs9373124 4.40e-108   6 135464902             T                 C       57790
#4 rs4895441 1.02e-107   6 135468266             A                 G       57837
#5 rs9376090 1.38e-107   6 135452921             T                 C       57990
#6 rs9389268 2.47e-107   6 135461324             A                 G       57856
#Beta     SE                    ChrBPAFInfo GC.Zscore
#1 -0.6003 0.0275   chr6_135468852_T_0.78_C_0.22  22.21303
#2 -0.6053 0.0277   chr6_135460609_A_0.78_G_0.22  22.10044
#3 -0.6170 0.0284   chr6_135464902_T_0.78_C_0.22  22.08506
#4 -0.5900 0.0272 chr6_135468266_A_0.775_G_0.225  22.04703
#5 -0.6047 0.0279   chr6_135452921_T_0.78_C_0.22  22.03334
#6 -0.6004 0.0276 chr6_135461324_A_0.775_G_0.225  22.00696
#> head(PCV)
#SNP      Pval Chr  Position Effect_Allele Non_Effect_Allele
#1  rs9373124 3.869e-22   6 135464902             T                 C
#2  rs9389269 1.966e-21   6 135468852             T                 C
#3  rs9402686 2.533e-21   6 135469510             A                 G
#4 rs11759553 4.364e-21   6 135463989             A                 T
#5  rs4895441 5.040e-21   6 135468266             A                 G
#6  rs4895440 5.451e-21   6 135468251             A                 T
#Sample_Size    Beta     SE                    ChrBPAFInfo GC.Zscore
#1       52764  0.1815 0.0198   chr6_135464902_T_0.78_C_0.22  9.674521
#2       52829  0.1737 0.0190   chr6_135468852_T_0.78_C_0.22  9.506809
#3       52827 -0.1742 0.0195 chr6_135469510_G_0.783_A_0.217  9.480404
#4       52825  0.1727 0.0189 chr6_135463989_A_0.775_T_0.225  9.423475
#5       52811  0.1706 0.0188 chr6_135468266_A_0.775_G_0.225  9.408347
#6       52827  0.1701 0.0187 chr6_135468251_A_0.767_T_0.233  9.400102
#> head(MCH)
#SNP      Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
#1 rs7776054 6.55e-108   6 135460609             A                 G       51452
#2 rs9389269 8.02e-108   6 135468852             T                 C       51452
#3 rs9389268 6.04e-107   6 135461324             A                 G       51453
#4 rs9373124 1.35e-106   6 135464902             T                 C       51388
#5 rs4895441 2.90e-106   6 135468266             A                 G       51434
#6 rs9376090 3.04e-106   6 135452921             T                 C       51588
#Beta     SE                    ChrBPAFInfo GC.Zscore
#1 -0.2399 0.0108   chr6_135460609_A_0.78_G_0.22  22.06707
#2 -0.2368 0.0108   chr6_135468852_T_0.78_C_0.22  22.05791
#3 -0.2382 0.0108 chr6_135461324_A_0.775_G_0.225  21.96638
#4 -0.2436 0.0111   chr6_135464902_T_0.78_C_0.22  21.92981
#5 -0.2322 0.0106 chr6_135468266_A_0.775_G_0.225  21.89499
#6 -0.2385 0.0109   chr6_135452921_T_0.78_C_0.22  21.89284
#> head(Hb)
#SNP      Pval Chr Position Effect_Allele Non_Effect_Allele Sample_Size
#1  rs855791 4.645e-40  22 35792882             A                 G       46184
#2 rs4820268 1.153e-31  22 35799537             A                 G       43990
#3 rs2413450 4.411e-31  22 35800170             T                 C       49624
#4  rs198846 1.422e-30   6 26215442             A                 G       60869
#5  rs129128 8.996e-27   6 26233321             T                 C       53598
#6 rs1799945 3.599e-26   6 26199158             C                 G       49746
#Beta     SE                    ChrBPAFInfo GC.Zscore
#1 -0.0791 0.0063 chr22_35792882_T_0.392_C_0.608  13.24782
#2  0.0699 0.0064 chr22_35799537_G_0.417_A_0.583  11.70850
#3 -0.0695 0.0064 chr22_35800170_T_0.382_C_0.618  11.59417
#4  0.0910 0.0080  chr6_26215442_A_0.142_G_0.858  11.49352
#5 -0.0925 0.0086  chr6_26233321_G_0.142_A_0.858  10.71143
#6 -0.0938 0.0089  chr6_26199158_C_0.871_G_0.129  10.58233
#> head(MCHC)
#SNP      Pval Chr Position Effect_Allele Non_Effect_Allele Sample_Size
#1 rs10445033 1.536e-22  16 87367963             A                 G       42050
#2   rs837763 1.928e-22  16 87381230             T                 C       37768
#3   rs475596 5.644e-21  16 87374452             C                 G       42124
#4   rs198846 7.458e-21   6 26215442             A                 G       56189
#5  rs9932423 2.263e-19  16 87374350             A                 C       41839
#6  rs2608604 4.232e-19  16 87376922             A                 G       37770
#Beta     SE                    ChrBPAFInfo GC.Zscore
#1 -0.0170 0.0038 chr16_87367963_G_0.367_A_0.633  9.768573
#2 -0.0147 0.0037 chr16_87381230_G_0.458_A_0.542  9.745514
#3 -0.0135 0.0038   chr16_87374452_C_0.39_G_0.61  9.396441
#4  0.0207 0.0048  chr6_26215442_A_0.142_G_0.858  9.367061
#5  0.0153 0.0041     chr16_87374350_C_0.7_A_0.3  8.999717
#6  0.0126 0.0038 chr16_87376922_T_0.358_C_0.642  8.930732
> head(RBC)
         SNP     Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
1  rs9373124 5.11e-97   6 135464902             T                 C       53337
2  rs7776054 3.51e-96   6 135460609             A                 G       53403
3  rs9389269 4.54e-96   6 135468852             T                 C       53402
4  rs9389268 5.98e-95   6 135461324             A                 G       53403
5  rs9376090 9.04e-95   6 135452921             T                 C       53537
6 rs11759553 1.84e-94   6 135463989             A                 T       53398
    Beta     SE Chr.1        BP Ref   RAF Alt   AAF GC.Zscore
1 0.0484 0.0024  chr6 135464902   T 0.780   C 0.220  20.90223
2 0.0479 0.0024  chr6 135460609   A 0.780   G 0.220  20.81005
3 0.0469 0.0023  chr6 135468852   T 0.780   C 0.220  20.79771
4 0.0471 0.0024  chr6 135461324   A 0.775   G 0.225  20.67366
5 0.0477 0.0024  chr6 135452921   T 0.780   C 0.220  20.65371
6 0.0465 0.0023  chr6 135463989   A 0.775   T 0.225  20.61935
> head(MCV)
        SNP     Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
1 rs9999997  0.63790   4 164089928             A                 G       57752
2 rs9999996 -0.49100   4  69817056             C                 A       57853
3 rs9999993 -0.33000   4  98781694             T                 A       57855
4 rs9999992 -0.70920   4 123121534             G                 A       14771
5 rs9999987 -0.07658   4   4987062             C                 T       46163
6 rs9999979  0.36150   4 135953341             C                 T       49932
    Beta     SE Chr.1        BP Ref   RAF Alt   AAF  GC.Zscore
1 0.0100 0.0246  chr4 164089928   A 0.550   G 0.450  0.4706370
2 0.0196 0.0346  chr4  69817056   A 0.192   C 0.192 -0.6887192
3 0.0220 0.0242  chr4  98781694   A 0.458   T 0.458 -0.9741139
4 0.0243 0.1531  chr4 123121534   A 0.948   G 0.948 -0.3729307
5 0.1148 0.0573  chr4   4987062   T 0.967   C 0.967 -1.7708840
6 0.0227 0.0269  chr4 135953341   T 0.525   C 0.525  0.9125106
> head(PCV)
        SNP    Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
1 rs9999997 -0.8793   4 164089928             G                 A       52726
2 rs9999996 -0.1436   4  69817056             C                 A       52827
3 rs9999993  0.5803   4  98781694             A                 T       52829
4 rs9999992  0.2738   4 123121534             A                 G       11834
5 rs9999987  0.7374   4   4987062             T                 C       41101
6 rs9999979  0.5293   4 135953341             C                 T       44891
    Beta     SE Chr.1        BP Ref   RAF Alt   AAF  GC.Zscore
1 0.0013 0.0170  chr4 164089928   A 0.450   G 0.450 -0.1518566
2 0.0296 0.0238  chr4  69817056   A 0.192   C 0.192 -1.4625155
3 0.0087 0.0167  chr4  98781694   A 0.542   T 0.458  0.5529466
4 0.0739 0.1092  chr4 123121534   A 0.052   G 0.948  1.0943534
5 0.0067 0.0424  chr4   4987062   T 0.033   C 0.967  0.3352984
6 0.0124 0.0191  chr4 135953341   T 0.525   C 0.525  0.6290749
> head(MCH)
        SNP    Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
1 rs9999997  0.2796   4 164089928             A                 G       51350
2 rs9999996 -0.6970   4  69817056             C                 A       51450
3 rs9999993 -0.2183   4  98781694             T                 A       51452
4 rs9999992 -0.6615   4 123121534             G                 A       12924
5 rs9999987 -0.2730   4   4987062             C                 T       42905
6 rs9999979 -0.5221   4 135953341             T                 C       43553
    Beta     SE Chr.1        BP Ref   RAF Alt   AAF  GC.Zscore
1 0.0095 0.0096  chr4 164089928   A 0.550   G 0.450  1.0812183
2 0.0039 0.0135  chr4  69817056   A 0.192   C 0.192 -0.3893733
3 0.0129 0.0095  chr4  98781694   A 0.458   T 0.458 -1.2310611
4 0.0163 0.0624  chr4 123121534   A 0.948   G 0.948 -0.4378431
5 0.0221 0.0223  chr4   4987062   T 0.967   C 0.967 -1.0961800
6 0.0070 0.0105  chr4 135953341   T 0.475   C 0.525 -0.6401117
> head(Hb)
        SNP    Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
1 rs9999997  0.5449   4 164089928             A                 G       60793
2 rs9999996 -0.4484   4  69817056             C                 A       60894
3 rs9999993  0.7488   4  98781694             A                 T       60896
4 rs9999992  0.2023   4 123121534             A                 G       11855
5 rs9999987 -0.8709   4   4987062             C                 T       49189
6 rs9999979 -0.3450   4 135953341             T                 C       52956
    Beta     SE Chr.1        BP Ref   RAF Alt   AAF  GC.Zscore
1 0.0016 0.0055  chr4 164089928   A 0.550   G 0.450  0.6054199
2 0.0051 0.0077  chr4  69817056   A 0.192   C 0.192 -0.7580852
3 0.0015 0.0054  chr4  98781694   A 0.542   T 0.458  0.3202221
4 0.0320 0.0377  chr4 123121534   A 0.052   G 0.948  1.2750261
5 0.0009 0.0130  chr4   4987062   T 0.967   C 0.967 -0.1625154
6 0.0063 0.0060  chr4 135953341   T 0.475   C 0.525 -0.9443320
> head(MCHC)
        SNP     Pval Chr  Position Effect_Allele Non_Effect_Allele Sample_Size
1 rs9999997 -0.96190   4 164089928             G                 A       56114
2 rs9999996  0.06551   4  69817056             A                 C       56214
3 rs9999993 -0.39500   4  98781694             T                 A       56216
4 rs9999992 -0.70300   4 123121534             G                 A       12410
5 rs9999987  0.66760   4   4987062             T                 C       44482
6 rs9999979 -0.10460   4 135953341             T                 C       48270
    Beta     SE Chr.1        BP Ref   RAF Alt   AAF   GC.Zscore
1 0.0037 0.0034  chr4 164089928   A 0.450   G 0.450 -0.04776943
2 0.0146 0.0048  chr4  69817056   A 0.808   C 0.192  1.84176180
3 0.0025 0.0034  chr4  98781694   A 0.458   T 0.458 -0.85058486
4 0.0334 0.0200  chr4 123121534   A 0.948   G 0.948 -0.38127392
5 0.0041 0.0083  chr4   4987062   T 0.033   C 0.967  0.42944419
6 0.0013 0.0035  chr4 135953341   T 0.475   C 0.525 -1.62295044
#~~~

library("data.table")

#create a merged table of Z scores
dt1 = data.table(rs=RBC$SNP,Z.RBC=RBC$GC.Zscore,key="rs")
dt2 = data.table(rs=MCV$SNP,Z.MCV=MCV$GC.Zscore,key="rs")
dt3 = merge(dt1,dt2)
rm(dt2)
rm(dt1)
dt1 = data.table(rs=PCV$SNP,Z.PCV=PCV$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=MCH$SNP,Z.MCH=MCH$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)
dt1 = data.table(rs=Hb$SNP,Z.Hb=Hb$GC.Zscore,key="rs")
dt2 = merge(dt3,dt1)
rm(dt1)
rm(dt3)
dt1 = data.table(rs=MCHC$SNP,Z.MCHC=MCHC$GC.Zscore,key="rs")
dt3 = merge(dt2,dt1)
rm(dt1)
rm(dt2)

dim(dt3)
head(dt3)

#~~~
#> dim(dt3)
#[1] 2584879       7
#> head(dt3)
#rs    Z.RBC     Z.MCV     Z.PCV     Z.MCH       Z.Hb    Z.MCHC
#1:       rs10 1.641952 0.4586298 1.9711213 0.7226418 1.70593598 0.7403296
#2:  rs1000000 2.063449 0.6926966 0.7624357 1.4939067 1.12379501 0.6507622
#3: rs10000010 1.025346 0.2026610 0.4221713 1.5023761 0.05228701 2.6894960
#4: rs10000012 1.822626 1.3243103 2.5439716 0.8050341 3.05695211 0.7767606
#5: rs10000013 0.443644 0.3555201 0.1554076 1.4238151 0.74561503 1.6017993
#6: rs10000017 1.567919 1.0909384 0.9255113 0.8228382 0.50408734 0.6485953
> dim(dt3)
[1] 2254972       7
> head(dt3)
           rs     Z.RBC      Z.MCV     Z.PCV      Z.MCH        Z.Hb     Z.MCHC
1:  rs1000000 2.0634489 -0.6926966 0.7624357 -1.4939067  1.12379501 -0.6507622
2: rs10000010 1.0253459 -0.2026610 0.4221713 -1.5023761  0.05228701 -2.6894960
3: rs10000012 1.8226264  1.3243103 2.5439716  0.8050341  3.05695211  0.7767606
4: rs10000013 0.4436440 -0.3555201 0.1554076 -1.4238151 -0.74561503 -1.6017993
5: rs10000017 1.5679192 -1.0909384 0.9255113 -0.8228382  0.50408734  0.6485953
6: rs10000023 0.7497588 -0.6533976 1.5402476 -1.6352340 -0.08922890 -2.5727296
>  
#~~~

RBC.maxZ <- max(dt3$Z.RBC[!is.infinite(dt3$Z.RBC) & !is.na(dt3$Z.RBC)])
MCV.maxZ <- max(dt3$Z.MCV[!is.infinite(dt3$Z.MCV) & !is.na(dt3$Z.MCV)])
PCV.maxZ <- max(dt3$Z.PCV[!is.infinite(dt3$Z.PCV) & !is.na(dt3$Z.PCV)])
MCH.maxZ <- max(dt3$Z.MCH[!is.infinite(dt3$Z.MCH) & !is.na(dt3$Z.MCH)])
Hb.maxZ <- max(dt3$Z.Hb[!is.infinite(dt3$Z.Hb) & !is.na(dt3$Z.Hb)])
MCHC.maxZ <- max(dt3$Z.MCHC[!is.infinite(dt3$Z.MCHC) & !is.na(dt3$Z.MCHC)])
maxZ <- max(c(RBC.maxZ, MCV.maxZ, PCV.maxZ, MCH.maxZ, Hb.maxZ, MCHC.maxZ))

RBC.maxZ
MCV.maxZ
PCV.maxZ
MCH.maxZ
Hb.maxZ
MCHC.maxZ
maxZ

#~~~
#> RBC.maxZ
#[1] 20.90223
#> MCV.maxZ
#[1] 22.21303
#> PCV.maxZ
#[1] 9.674521
#> MCH.maxZ
#[1] 22.06707
#> Hb.maxZ
#[1] 13.24782
#> MCHC.maxZ
#[1] 9.768573
#> maxZ
#[1] 22.21303
> RBC.maxZ
[1] 20.90223
> MCV.maxZ
[1] 9.854182
> PCV.maxZ
[1] 9.674521
> MCH.maxZ
[1] 11.64786
> Hb.maxZ
[1] 9.377884
> MCHC.maxZ
[1] 9.768573
> maxZ
[1] 20.90223
#~~~

replaceInf <- function(x) { if (is.infinite(x)) { print("yaya1"); return(maxZ) } else { return(x) } }

dt3$Z.RBC <- apply(as.matrix(dt3$Z.RBC), 1, replaceInf)
dt3$Z.MCV <- apply(as.matrix(dt3$Z.MCV), 1, replaceInf)
dt3$Z.PCV <- apply(as.matrix(dt3$Z.PCV), 1, replaceInf)
dt3$Z.MCH <- apply(as.matrix(dt3$Z.MCH), 1, replaceInf)
dt3$Z.Hb <- apply(as.matrix(dt3$Z.Hb), 1, replaceInf)
dt3$Z.MCHC <- apply(as.matrix(dt3$Z.MCHC), 1, replaceInf)

#~~~
#> dt3$Z.RBC <- apply(as.matrix(dt3$Z.RBC), 1, replaceInf)
#> dt3$Z.MCV <- apply(as.matrix(dt3$Z.MCV), 1, replaceInf)
#> dt3$Z.PCV <- apply(as.matrix(dt3$Z.PCV), 1, replaceInf)
#> dt3$Z.MCH <- apply(as.matrix(dt3$Z.MCH), 1, replaceInf)
#> dt3$Z.Hb <- apply(as.matrix(dt3$Z.Hb), 1, replaceInf)
#> dt3$Z.MCHC <- apply(as.matrix(dt3$Z.MCHC), 1, replaceInf)
> 
> dt3$Z.RBC <- apply(as.matrix(dt3$Z.RBC), 1, replaceInf)
> dt3$Z.MCV <- apply(as.matrix(dt3$Z.MCV), 1, replaceInf)
> dt3$Z.PCV <- apply(as.matrix(dt3$Z.PCV), 1, replaceInf)
> dt3$Z.MCH <- apply(as.matrix(dt3$Z.MCH), 1, replaceInf)
> dt3$Z.Hb <- apply(as.matrix(dt3$Z.Hb), 1, replaceInf)
> dt3$Z.MCHC <- apply(as.matrix(dt3$Z.MCHC), 1, replaceInf)
> 
#~~~

attach(dt3)
nullset = (abs(dt3$Z.MCV)<2) & (abs(dt3$Z.RBC)<2) & (abs(dt3$Z.PCV)<2) & (abs(dt3$Z.MCH)<2) & (abs(dt3$Z.Hb)<2) & (abs(dt3$Z.MCHC)<2) #extract null Z values
Z = cbind(Z.MCV,Z.RBC,Z.PCV,Z.MCH,Z.Hb,Z.MCHC)

#compute correlation based only on "null" results
Znull = Z[nullset,]
RSS0 =cor(Znull)
RSS0inv = chol2inv(chol(RSS0))

dim(Z)
dim(Znull)
RSS0

#~~~
#> dim(Z)
#[1] 2584879       6
#> dim(Znull)
#[1] 1985737       6
#> RSS0
#Z.MCV        Z.RBC         Z.PCV        Z.MCH         Z.Hb
#Z.MCV   1.00000000  0.097743042 -0.0108195734  0.441583579 -0.019501098
#Z.RBC   0.09774304  1.000000000  0.3204545474  0.096099964  0.273541362
#Z.PCV  -0.01081957  0.320454547  1.0000000000 -0.025642749  0.540784949
#Z.MCH   0.44158358  0.096099964 -0.0256427495  1.000000000  0.007745759
#Z.Hb   -0.01950110  0.273541362  0.5407849492  0.007745759  1.000000000
#Z.MCHC -0.01086963 -0.004515348 -0.0003855153  0.083933220  0.025819445
#Z.MCHC
#Z.MCV  -0.0108696265
#Z.RBC  -0.0045153480
#Z.PCV  -0.0003855153
#Z.MCH   0.0839332201
#Z.Hb    0.0258194449
#Z.MCHC  1.0000000000
> dim(Z)
[1] 2254972       6
> dim(Znull)
[1] 1733999       6
> RSS0
             Z.MCV       Z.RBC       Z.PCV      Z.MCH      Z.Hb      Z.MCHC
Z.MCV   1.00000000 -0.26755719  0.39345589  0.6910826 0.3261603 -0.04486257
Z.RBC  -0.26755719  1.00000000  0.43657397 -0.2677200 0.4046120 -0.02160634
Z.PCV   0.39345589  0.43657397  1.00000000  0.2792027 0.6984009 -0.04792499
Z.MCH   0.69108258 -0.26772002  0.27920274  1.0000000 0.4432557  0.23204003
Z.Hb    0.32616034  0.40461197  0.69840088  0.4432557 1.0000000  0.18783496
Z.MCHC -0.04486257 -0.02160634 -0.04792499  0.2320400 0.1878350  1.00000000
> 
#~~~

mvstat =  rowSums(Z * (Z %*% RSS0inv)) # comptues Z RSS0^-1 Z'
dt3$mvstat = mvstat
statchi = -log10(exp(1))*pchisq(mvstat,df=6,log.p=TRUE,lower.tail=FALSE)
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
#[1] 3773   10
#> dim(dtlesssignif)
#[1] 5055   10
#> max(dt3$mvstat)
#[1] 1065.289
#> max(dt3$mvp)
#[1] 226.1711
#> max(dt3$unip)
#[1] 108.5901
#> quantile(dt3$mvp)
#0%          25%          50%          75%         100%
#6.773086e-09 4.564658e-02 1.678489e-01 4.596310e-01 2.261711e+02
#> quantile(dt3$unip)
#0%          25%          50%          75%         100%
#0.02691835   0.56815395   0.86934465   1.29774184 108.59006688
#>
> dim(dtsignif)      
[1] 3020   10
> dim(dtlesssignif)
[1] 4200   10
> max(dt3$mvstat)
[1] 1051.804
> max(dt3$mvp)
[1] 223.2538
> max(dt3$unip)
[1] 108.5901
> quantile(dt3$mvp)
          0%          25%          50%          75%         100% 
4.380744e-08 8.241575e-02 2.277952e-01 5.192527e-01 2.232538e+02 
> quantile(dt3$unip)
          0%          25%          50%          75%         100% 
  0.02691835   0.56815395   0.86838134   1.29602117 108.59006688 
#~~~

write.table(file="HaemgenRBC2012.dtsignif.vs1.SignCrrct.vs1.txt",dtsignif,sep=",",row.names=F,quote=F)
write.table(file="HaemgenRBC2012.dtsignif.rs.vs1.SignCrrct.vs1.txt",dtsignif$rs,sep=",",row.names=F,quote=F)
write.table(file="HaemgenRBC2012.dtlesssignif.vs1.SignCrrct.vs1.txt",dtlesssignif,sep=",",row.names=F,quote=F)
write.table(file="HaemgenRBC2012.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",dtlesssignif$rs,sep=",",row.names=F,quote=F)
#write.table(file="HaemgenRBC2012.dtlesssignif.rs.vs1.SignCrrct.vs1.txt",paste(dtlesssignif$rs,"[[:space:]]",sep=""),row.names=F,quote=F)

write.table(file="HaemgenRBC2012.RSS0.vs1.SignCrrct.vs1.txt",RSS0,row.names=F,sep=",",quote=F)


