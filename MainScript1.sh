#!/bin/sh

20140611
~~~~~~~

cd /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013

git clone https://github.com/stephens999/multivariate

mkdir /data/external_public/GlobalLipids2013

cd /data/external_public/GlobalLipids2013

wget -r -l1 --no-parent -A "jointGwas*" http://www.sph.umich.edu/csg/abecasis/public/lipids2013/

wget -r -l1 --no-parent -A "Mc*" http://www.sph.umich.edu/csg/abecasis/public/lipids2013/

#NOTE -- these files were actually already downloaded by Xiang Zhu in /data/external_public/GlobalLipid
#So removed the /data/external_public/GlobalLipids2013 directory

mkdir /data/external_public/GlobalLipids2010

cd /data/external_public/GlobalLipids2010

wget -r -l1 --no-parent -A "*2010*" http://www.sph.umich.edu/csg/abecasis/public/lipids2010/

wget -r -l1 --no-parent -A "*with_Effect*" http://www.sph.umich.edu/csg/abecasis/public/lipids2010/

#Unzipped both *2010* and *with_Effect* files

mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013 

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/process.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.vs1.R

#Confirming column 2 in file /data/external_public/GlobalLipid/Mc_LDL.txt.gz has only unique entries
zcat /data/external_public/GlobalLipid/Mc_LDL.txt.gz | awk '{ print $2 }' | sort | uniq -c | sort -k 1,1 | tail -n 20

cd /mnt/lustre/home/mturchin20/Software

wget http://cran.r-project.org/src/contrib/data.table_1.9.2.tar.gz

#install.packages("/mnt/lustre/home/mturchin20/Software/data.table_1.9.2.tar.gz", repos=NULL, type="source")

scp /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.vs1.R mturchin20@wolfy.uchicago.edu:/Users/mturchin20/LabStuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/. 

cp -p /Users/mturchin20/LabStuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.vs1.R /Users/mturchin20/LabStuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.For2013.vs1.R

#Copy and pasted commands from /Users/mturchin20/LabStuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.For2013.vs1.R into R terminal, but could have run script as well

scp mturchin20@wolfy.uchicago.edu:/Users/mturchin20/LabStuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.*.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/.

##NOTE_A (see below)
#Supplementary file ng.2797-S1.pdf was downloaded from http://www.nature.com/ng/journal/v45/n11/extref/ng.2797-S1.pdf and the .pdf was saved as a .txt file using Adobe Acrobat
#In excel, all tables except for supplementary tables 2 and 3 were removed manually, producing file ng.2797-S1.edited.vs2.txt
#Supplementary table 2 contains 62 novel loci and supplementary table 3 contains 157 previously identified loci

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.txt | perl -lane 'foreach my $val1 (@F) { print $val1 ; } ' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.txt

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.txt | grep ^rs > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.txt

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt | awk '{ print $3 }' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.rsIDs.txt

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.txt | grep -v -f /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.rsIDs.txt > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.missingSNPs.txt

#/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.txt was run through HapMart to get the chromosome and position locations for each rsID. 2 SNPs were missing and their information was manually included. These 2 SNPs were rs1047891 and rs9411489. Information was retrieved from UCSC Genome Browser hg18 and included in the /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt file
#hg18 because HapMart only had hg18 but switching to hg19 with the final annotated file
#rs1047891 chr2 211248752
#rs9411489 chr 
#NOTE -- Previously reported ABO SNP rs9411489 has merged into rs635634 (source: supplementary of http://hmg.oxfordjournals.org/content/21/6/1444.full or https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=0CCsQFjAB&url=http%3A%2F%2Fhmg.oxfordjournals.org%2Fcontent%2Fsuppl%2F2011%2F12%2F07%2Fddr581.DC1%2Fddr581supp_table1.doc&ei=HxKiU8_OIYOOyAS-tICoCg&usg=AFQjCNHO1Honaz_0R7CcVLsZan3LNJulKQ&sig2=TIdN-hlZUyRQZeMxiDGWHg&bvm=bv.69137298,d.aWw)
#rs635634 chr9 135144821

mkdir /mnt/lustre/home/mturchin20/Scripts/Python

#Below scp done in Miranda server from NU
scp /home/michaelt/Scripts/Python/2fileTemplate.vs1.py mturchin20@wolfy.uchicago.edu:/Users/mturchin20/clstrHme/. 

mv /mnt/lustre/home/mturchin20/2fileTemplate.vs1.py /mnt/lustre/home/mturchin20/Scripts/Python/2fileTemplate.vs1.py

cp -p /mnt/lustre/home/mturchin20/Scripts/Python/2fileTemplate.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateFullAnnotFile.vs1.py

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateFullAnnotFile.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateProcessFile.vs1.py

#python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateProcessFile.vs1.py --file1 /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz,/data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz,/data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz,/data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz > 

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt --file2 /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz,/data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz,/data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz,/data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1; } splice(@F, 1, 1, ($chr, $vals1[1])); push(@F, "NA"); print join("\t", @F);' | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.txt.gz 

#For troubleshooting purposes with the /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateFullAnnotFile.vs1.py script
zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | head -n 20 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.txt.top20
zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | head -n 20 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.txt.top20
zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | head -n 20 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.txt.top20
zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | head -n 20 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.txt.top20
zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | head -n 5000 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.txt.top5000
zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | head -n 5000 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.txt.top5000
zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | head -n 5000 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.txt.top5000
zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | head -n 5000 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.txt.top5000

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt --file2 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.txt.top20,/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.txt.top20,/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.txt.top20,/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.txt.top20 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.txt.top20
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt --file2 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.txt.top5000,/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.txt.top5000,/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.txt.top5000,/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.txt.top5000 > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.txt.top5000

#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.txt.gz
#posHg18 chr pos snp a1 a2 maf beta_LDL se_LDL n_LDL beta_HDL se_HDL n_HDL beta_TG se_TG n_TG beta_TC se_TC n_TC annot gene 

cp -p /mnt/lustre/home/mturchin20/Scripts/Python/2fileTemplate.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.txt --file2 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.txt.gz

#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.txt.gz
#posHg18 chr pos snp a1 a2 maf beta_LDL se_LDL n_LDL beta_HDL se_HDL n_HDL beta_TG se_TG n_TG beta_TC se_TC n_TC annot gene Z.tg Z.tc Z.ldl Z.hdl mvstat mvp unip

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/GLC.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs1.R

cp -p /mnt/lustre/home/mturchin20/Scripts/Python/2fileTemplate.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.UpdateAnnotFileAnnotations.vs1.py

##NOTE_B
#sebanti.novel_lipid_loci.txt and sebanti.known_lipid_loci.txt were provided by Sebanti Sengupta and Xiaoquan Wen from UMichigan after being asked by Matthew

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/sebanti.known_lipid_loci.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/sebanti.novel_lipid_loci.txt > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/sebanti.all_lipid_loci.txt
#NOTE -- replaced rs9411489 with rs635634 in /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/sebanti.all_lipid_loci.txt file since /data/external_public/GlobalLipid/jointGwasMc_* files do not contain rs9411489 but contain rs635634
#Running this code paste <(cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt | awk '{ if (NR > 1) { print $3 } }' | sort) <(cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/sebanti.all_lipid_loci.txt | sort) | vi - makes it look like the previous list of rsIDs -- the ones scraped from the 2013 paper's supplement -- and the list sent to us by Sebanti, are in fact the same
#As a result, going to continue using the list of rsIDs and information from the HapMart set of information since we use positional information to tag 'nearby' SNPs. Sebanti did not send this information (we did not ask for it), so I would just be redoing the HapMart information and getting the same results


#NOTE -- I was using qnorm(p-value) to get Z-scores, but this was producing one-sided Z-scores when the original test statistic was two-sided (the original H0 is Beta == 0 and H1 is Beta != 0). So need to do qnorm(p-value/2, lower.tail=FALSE) to get proper Z-scores. So redoing things to get those results
cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.For2013.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.For2013.vs2.R

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs2.R

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.vs2.txt --file2 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz

#Manually searched UCSC Genome Browser, GTEx portal and dbSNP for information regarding the top 74 hits that come out of our /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs2.R analysis (e.g. SNPs with lbfavs > 4.6 or so, the minimum of the prior hits)

zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz | perl -lane 'print $F[1], "\t", $F[2], "\t", $F[2], "\t", uc($F[4]), "\t", uc($F[5]);' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.AnnovarFormat

#Manually deleted header line from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.AnnovarFormat

mv /mnt/lustre/home/mturchin20/Software/annovar.latest.tar.gz 20130509_annovar.latest.tar.gz

wget http://www.openbioinformatics.org/annovar/download/Ht8qRwQSTi/annovar.latest.tar.gz

mkdir /mnt/lustre/home/mturchin20/Software/20130509annovar

mv /mnt/lustre/home/mturchin20/Software/annovar/* /mnt/lustre/home/mturchin20/Software/20130509annovar/.

tar -xvzf /mnt/lustre/home/mturchin20/Software/annovar.latest.tar.gz

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb -webfrom annovar refGene /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb genomicSuperDups /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb -webfrom annovar snp129 /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb -webfrom annovar avsift /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb -webfrom annovar ljb_all /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb -webfrom annovar ljb23_all /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb -webfrom annovar esp6500_all /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb -webfrom annovar 1000g2012apr /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

perl /mnt/lustre/home/mturchin20/Software/annovar/annotate_variation.pl --downdb phastConsElements46way /mnt/lustre/home/mturchin20/Software/annovar/humandb/ -build hg19

cd /mnt/lustre/home/mturchin20/Software/annovar/humandb/. ; nbthis *

cd /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1

perl /mnt/lustre/home/mturchin20/Software/annovar/summarize_annovar.pl /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.AnnovarFormat /mnt/lustre/home/mturchin20/Software/annovar/humandb/ --ver1000g 1000g2012apr --verdbsnp 129 --veresp 6500 --buildver hg19 --remove

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt | awk '{ print $1 }' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.rsIDs

cat dtlesssignif.annot.vs2.txt.AnnovarFormat.genome_summary.csv | grep -f /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.rsIDs > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.genome_summary.csv
#Doing above line of code only got 39 of the 74 hits

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt | awk '{ print $3 }' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.pos

cat dtlesssignif.annot.vs2.txt.AnnovarFormat.genome_summary.csv | grep -f /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.pos > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.genome_summary.csv

[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.AnnovarFormat.genome_summary.csv | head -n 1 | perl -F, -lane 'print $F[19];'
LJB_MutationTaster_Pred

#cat newhits.vs2.genome_summary.csv | perl -F, -lane 'print join(",", @F[0..8]), ",", $F[13], ",", $F[15], ",", $F[17], ",", $F[19];' | paste <(cat newhits.vs2.txt | tail -n +2 | perl -lane 'print join(",", @F);') - | perl -lane 'print join(",", @F);' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.txt
cat newhits.vs2.genome_summary.csv | perl -F, -lane 'if ($#F == 26) { $F[1] = $F[1] . ";".  $F[2]; splice(@F, 2, 1);} if ($F[2] =~ m/synon/) { my @vals1 = split(/\s+/, $F[2]); $F[2] = join("_", @vals1); } print join(",", @F[0..8]), ",", $F[13], ",", $F[15], ",", $F[17], ",", $F[19];' | paste <(cat newhits.vs2.txt | tail -n +2 | perl -lane 'print join(",", @F);') - | perl -lane 'print join(",", @F);' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.txt
#Manually included the below header into /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.txt
#snp,chr,pos,maf,Z.tg,Z.tc,Z.ldl,Z.hdl,lbfav,lbfall,lbfuni,Func,Gene,ExonicFunc,AAChange,Conserved,SegDup,ESP6500_ALL,1000g2012apr_ALL,dbSNP129,LJB_SIFT_Pred,LJB_PolyPhen2_Pred,LJB_LRT_Pred,LJB_MutationTaster_Pred
#Manually included Novel,Previous,Enriched categories by checking whether each gene ID was present in 2013 paper supplement; also check whether annotations made sense using UCSC Genome Browser via rsIDs

#NOTE_C
#ERLIN1 previously reported with rs1408579 by Matthew -- rs2862954 is top hit now. Differences shown below
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs1408579
chr10:101902184 chr10:101912194 rs1408579       t       c       0.0198  0.0048  94205.00        1.787e-05       0.4617
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs1408579
chr10:101902184 chr10:101912194 rs1408579       t       c       0.0028  0.0052  89784.00        0.4581  0.4617
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs1408579
chr10:101902184 chr10:101912194 rs1408579       t       c       0.0157  0.0051  94484.00        0.0009035       0.4617
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs1408579
chr10:101902184 chr10:101912194 rs1408579       t       c       0.0093  0.0047  90902.00        0.04111 0.4617

[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs2862954
chr10:101902054 chr10:101912064 rs2862954       c       t       0.0166  0.0034  186893.00       1.287e-06       0.4631
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs2862954
chr10:101902054 chr10:101912064 rs2862954       c       t       0.0013  0.0037  172821.00       0.5875  0.4631
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs2862954
chr10:101902054 chr10:101912064 rs2862954       c       t       0.0124  0.0035  187083.00       0.0002526       0.4631
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs2862954
chr10:101902054 chr10:101912064 rs2862954       c       t       0.0083  0.0033  177587.10       0.01393 0.4631

#rs12739698 5' upstream of NROB2 did not show up, and there's nothing nearby that is a top hit
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | grep -w rs12739698
chr1:27102620   chr1:27230033   rs12739698      g       a       0.0422  0.0088  94311.00        8.852e-07       0.92744
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz | grep -w rs12739698
chr1:27102620   chr1:27230033   rs12739698      a       g       0.0442  0.0095  89888.00        2.67e-05        0.07256
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz | grep -w rs12739698
chr1:27102620   chr1:27230033   rs12739698      a       g       0.0340  0.0093  94595.00        0.001193        0.07256
[  mturchin20@spudling17  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$zcat /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz | grep -w rs12739698
chr1:27102620   chr1:27230033   rs12739698      a       g       0.0271  0.0083  91013.00        0.0008306       0.07256
#NOTE -- gene is spelt NR0B2, it is in fact tagged by a 'novel' hit in the 2013 paper and the above SNP is within .5Mb

#rs10490632 is within .5Mb of a hit, rs10490626 which is in INSIG2 which is downstream of DDX2

#rs762861 is within .5Mb of a hit, rs6831256 which is in LRPAP1 which is downstream of HGFAC/RGS12

#rs2862954 is not within .5Mb of a hit, but it is in ERLIN1 which is upstream of BLOC1S2 and CHUK, both of which come up as 'enriched genes' in the pathway analyses. BF is on the lower end (4.8) in the Matthew paper table 2, though other loci with lower BFs did get picked up, e.g. VEGFA with 4.5 and STAB1 with 4.7

#Both rs11229252 and rs11227638 are in the OR4* (olfactory receptor) region of chromosome 11 -- they are also not showing up via those rsID #s in the new dataset, but presumably they are being considered 'tagged' by being nearby any type of OR4*, which has already shown up in the 2013 paper in the form of OR4C46

##NOTE_E


##NOTE_D1 -- MAGENTA work
wget http://www.broadinstitute.org/mpg/magenta/MAGENTA_software_package_vs2_hg18_hg19_July2011_Feb16_12.tar.gz

tar -xvzf /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011/MAGENTA_software_package_vs2_hg18_hg19_July2011_Feb16_12.tar.gz

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt | awk '{ print $3 }' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.rsIDs

cp -p /mnt/lustre/home/mturchin20/Scripts/Python/2fileTemplate.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt --file2 /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.PrevHits.vs1.txt 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt --file2 /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.PrevHits.vs1.txt 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt --file2 /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.PrevHits.vs1.txt 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt --file2 /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.PrevHits.vs1.txt 

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.PrevHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL_PrevHits_vs1_MAGENTA_txt
cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.PrevHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL_PrevHits_vs1_MAGENTA_txt
cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.PrevHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC_PrevHits_vs1_MAGENTA_txt
cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.PrevHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG_PrevHits_vs1_MAGENTA_txt

zcat /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.MAGENTA.txt 

ssh mturchin20@hegel.uchicago.edu

mv /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011/ /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011_Feb16_12/

wget http://www.broadinstitute.org/mpg/magenta/MAGENTA_software_package_vs2_hg18_hg19_July2011.tar.gz

#In /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011_Feb16_12/
ln -s /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.PrevHits.vs1.MAGENTA.txt jointGwasMc_HDL_PrevHits_vs1_MAGENTA_txt

#On mturchin20@hegel.uchicago.edu

mkdir Lipids

wget http://www.broadinstitute.org/mpg/magenta/MAGENTA_software_package_vs2_hg18_hg19_July2011_Feb16_12.tar.gz 

scp /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.PrevHits.vs1.MAGENTA.txt mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/.

mv jointGwasMc_HDL.PrevHits.vs1.MAGENTA.txt jointGwasMc_HDL_PrevHits_vs1_MAGENTA_txt

scp /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011_Feb16_12/Run_MAGENTA_vs2_July_2011.m mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/.

scp /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.MAGENTA.txt mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/.

mv /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/jointGwasMc_HDL.MAGENTA.txt /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/jointGwasMc_HDL_MAGENTA_txt

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.txt | perl -F, -ane 'print $F[0], "\t", $F[1], "\t", $F[2], "\t", $F[12], "\n";' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.ChrBP_And_Genes.txt

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.ChrBP_And_Genes.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.ChrBP_And_Genes.wEdits.txt

#Manually edited /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.ChrBP_And_Genes.wEdits.txt after looking at each SNP/gene in UCSC and doing a Google/literature search 

scp /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.rsIDs mturchin20@wolfy.uchicago.edu:/Users/mturchin20/LabStuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/.

cat /home/mturchin20/Lipids/Output_MAGENTA_GlobalLipids2013_HDL_AllSNPs_10000perm_Jul03_14/MAGENTA_pval_GeneSetEnrichAnalysis_GlobalLipids2013_HDL_AllSNPs_110kb_upstr_40kb_downstr_10000perm_Jul03_14.results | perl -lane 'if ($#F > 17) { my $diff = $#F - 17; my $end = 1 + $diff; $F[1] = join("_", @F[1..$end]); splice(@F, 2, $diff); } print join("\t", @F);' | perl -lane 'print $#F;' | sort | uniq -c


cp -p /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_prevHits_HDL.m /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_prevHits_LDL.m
cp -p /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_prevHits_HDL.m /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_prevHits_TC.m
cp -p /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_prevHits_HDL.m /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_prevHits_TG.m
cp -p /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_prevHits_HDL.m /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_currHits_HDL.m

#Editted above files to include proper inputs files, e.g. jointGwasMc_*.PrevHits.vs1.txt or jointGwasMc_*.CurrHits.vs1.txt

scp /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_*_PrevHits_vs1_MAGENTA_txt mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/.

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.rsIDs --file2 /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.newHits.vs1.txt 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.rsIDs --file2 /data/external_public/GlobalLipid/jointGwasMc_LDL.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.newHits.vs1.txt 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.rsIDs --file2 /data/external_public/GlobalLipid/jointGwasMc_TC.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.newHits.vs1.txt 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt.rsIDs --file2 /data/external_public/GlobalLipid/jointGwasMc_TG.txt.gz > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.newHits.vs1.txt 

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.PrevHits.vs1.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL.newHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_HDL_CurrHits_vs1_MAGENTA_txt
cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.PrevHits.vs1.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL.newHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_LDL_CurrHits_vs1_MAGENTA_txt
cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.PrevHits.vs1.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC.newHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TC_CurrHits_vs1_MAGENTA_txt
cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.PrevHits.vs1.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG.newHits.vs1.txt | perl -lane 'my @vals1 = split(/:/, $F[1]); my $chr = "NA"; if ($vals1[0] =~ m/chr(\d+)/) { $chr = $1 }; print $chr, "\t", $vals1[1], "\t", $F[8];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_TG_CurrHits_vs1_MAGENTA_txt

scp /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Data/jointGwasMc_*_CurrHits_vs1_MAGENTA_txt mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/.

cp -p /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_currHits_HDL.m /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_currHits_LDL.m
cp -p /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_currHits_HDL.m /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_currHits_TC.m
cp -p /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_currHits_HDL.m /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Run_MAGENTA_vs2_July_2011_currHits_TG.m

##NOTED_2

#Back in directory /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1

zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.vs2.txt.gz | perl -lane 'print $F[1], "\t", $F[2], "\t", $F[2], "\t", uc($F[4]), "\t", uc($F[5]);' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.vs2.txt.AnnovarFormat

perl /mnt/lustre/home/mturchin20/Software/annovar/summarize_annovar.pl /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.vs2.txt.AnnovarFormat /mnt/lustre/home/mturchin20/Software/annovar/humandb/ --ver1000g 1000g2012apr --verdbsnp 129 --veresp 6500 --buildver hg19 --remove

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.PutHG19IntoHapMartResults.vs1.py

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.PutHG19IntoHapMartResults.vs1.py --file1 /data/external_public/GlobalLipid/jointGwasMc_HDL.txt.gz --file2 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.txt > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.txt

#Manually included the below header in /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.txt
#chromosome      position        marker id       alleles population id   reference allele

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.Basic.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.txt --file2 jointGwasMc_AllPheno.Annot.vs2.txt.AnnovarFormat.genome_summary.csv > jointGwasMc_AllPheno.Annot.vs2.txt.AnnovarFormat.genome_summary.prevHits.csv

cp -p /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011/Extract_genes_around_SNPs.m /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011/Extract_genes_around_SNPs.dtlesssignif.m

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.txt | perl -F, -lane 'print $F[0], "\t", $F[1], "\t", $F[2];' > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.Extract_genes_around_SNPs.Format.txt

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.txt | perl -lane 'my $chr = "NA"; if ($F[0] =~ m/chr(\d+)/) { $chr = $1; } print $F[2], "\t", $chr, "\t", $F[1];' | grep -v marker > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.Extract_genes_around_SNPs.Format.txt

zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.txt.gz | perl -lane 'print $F[3], "\t", $F[1], "\t", $F[2]; ' | grep -v chr > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.Extract_genes_around_SNPs.Format.txt

#NOTE -- 50 SNPs do not have rsIDs, so they might fail 'Extract_genes_around_SNPs.m'. These SNPs are visible via the code below
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.txt.gz | perl -lane 'if ($F[3] !~ /rs/) { print join("\t", @F); }' | wc

scp /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.wAnnovar.Extract_genes_around_SNPs.Format.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.Extract_genes_around_SNPs.Format.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.Extract_genes_around_SNPs.Format.txt mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/.

#Ran /mnt/lustre/home/mturchin20/Software/MAGENTA_software_package_vs2_July2011/Extract_genes_around_SNPs.m in matlab on mturchin20@hegel.uchicago.edu
# Extract_genes_around_SNPs('ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.Extract_genes_around_SNPs.Format.txt',100000,'GlobalLipids2013_prevHits')
#Output file names are: Genes_near_SNP_100kb_boundary_GlobalLipids2013_prevHits_Jul10_14 and Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_prevHits_Jul10_14
#There are 27 SNPs with no genes within 100000 base pairs around the SNP.
#There are a total of 458 genes within 100000 base pairs around all SNPs.
#
#Extract_genes_around_SNPs('newhits.vs2.wAnnovar.Extract_genes_around_SNPs.Format.txt',100000,'GlobalLipids2013_currHits')
#Output file names are: Genes_near_SNP_100kb_boundary_GlobalLipids2013_currHits_Jul10_14 and Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_currHits_Jul10_14
#There are 12 SNPs with no genes within 100000 base pairs around the SNP.
#There are a total of 215 genes within 100000 base pairs around all SNPs.
#
#Output file names are: Genes_near_SNP_100kb_boundary_GlobalLipids2013_dtlesssignif_Jul10_14 and Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_dtlesssignif_Jul10_14
#There are 2024 SNPs with no genes within 100000 base pairs around the SNP.
#There are a total of 36341 genes within 100000 base pairs around all SNPs.
#
#Output is in mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011

#Issue with the following line in supplement of 2013 paper
#Locus Lead SNP Chr hg19 Position (Mb) Traits GWS Nearest Gene Nearest Gene (kb)   No. of Genes within 100kb
#ASAP3 rs1077514 1 23.77 TC ASAP3 0 6
#'no of genes within 100kb' is 6 but get 3 when run the following command Extract_genes_around_SNPs('ng.2797-S1.edited.vs2.carriageRemoved.rsIDs.HapMart.hg19.Extract_genes_around_SNPs.Format.txt',100000,'GlobalLipids2013_prevHits') and don't see 6 when look in UCSC genome browser using Hg19 database (see like 3 or so within a window of 100kb with SNP at center)
#However, interestingly, don't see ASAP3 anywhere in MAGENTA results tables, and ASAP3 doesn't exist in the GeneID data they provide to automatically use. But a bunch of hand-picked Gene IDs presented in the MAGENTA supplementary tables do exist in the GeneID list MAGENTA provides -- maybe they used the lists given and then hand-curated a number of the SNPs later?

scp mturchin20@hegel.uchicago.edu:/home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP* mturchin20@wolfy.uchicago.edu:/Users/mturchin20/clstrHme/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/. 

cat /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_prevHits_Jul10_14 | awk '{ print $5 }' | sort | uniq | perl -lane 'if ($F[0] !~ m/^$/) { print join("\t", @F); }' | grep -v \# > /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GL2013_prevHits_GeneIDs

cat /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_currHits_Jul10_14 | awk '{ print $5 }' | sort | uniq | perl -lane 'if ($F[0] !~ m/^$/) { print join("\t", @F); }' | grep -v \# > /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GL2013_currHits_GeneIDs

cat /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_dtlesssignif_Jul10_14 | awk '{ print $5 }' | sort | uniq | perl -lane 'if ($F[0] !~ m/^$/) { print join("\t", @F); }' | grep -v \# > /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GL2013_dtlesssignif_GeneIDs

cat /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_prevHits_Jul10_14 | awk '{ print $4 }' | sort | uniq | perl -lane 'if ($F[0] !~ m/^$/) { print join("\t", @F); }' > /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GL2013_prevHits_GeneSymbols

cat /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_currHits_Jul10_14 | awk '{ print $4 }' | sort | uniq | perl -lane 'if ($F[0] !~ m/^$/) { print join("\t", @F); }' > /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GL2013_currHits_GeneSymbols

cat /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GlobalLipids2013_dtlesssignif_Jul10_14 | awk '{ print $4 }' | sort | uniq | perl -lane 'if ($F[0] !~ m/^$/) { print join("\t", @F); }' > /home/mturchin20/Lipids/MAGENTA_software_package_vs2_July2011/Distance_genes_near_SNP_100kb_boundary_GL2013_dtlesssignif_GeneSymbols









##20150515
##Going back into this project and need to re-understand what I had done
##I will be explaining things below and continuing my log of things from this point on. I will not copy/paste relevant lines of code from above but I will give context and annotations to what was important/what will be used/moved going forward
##The goals, as of writing this header, are the following:
##	Identify where 2013 paper recapitulated results of Matthew's earlier paper
##	Identify new top hits compared to 2013 paper
##	Do a sort of SUMSTAT based on BF? to get a sense of 'enrichment' for relevant biological pathways
##	Do same setup in GIANT data with older and newer paper to see if anything gained?

##	Note -- Matthew was particularly interested in variants there were close to the threshold but then went over it in our analysis versus in the published paper's analysis ??

##NOTE_A: I manually took all the supplementary tables from the 2013 paper and tried to identify ever rsID# present in order to roughly determine what variants we had 'discovered' were also present in the 2013 paper in any regard. This approach is clearly making a number of assumptions but as a first approximation should probably help guide/inform things

##NOTE_B: I received a list of SNPs, those which were already known and those which were new in 2013 paper (known_lipid_loci.txt and novel_lipid_loci.txt), from UMich people and identified that my 'scraping' of the .pdf pulled the same list. So should have the right list of SNPs now from 2013 paper.

##NOTE_C: I think this is me looking for Matthew's previous SNPs in the newer data to try and figure out, if the rsID was not specifically present in the 2013 results, if Matthew's result was being tagged by something else

##NOTE_D1: The following section is mostly MAGENTA work on 'Hagel', the computer John helped me get one, and the portion at/after NOTE_D2 is on the PPS cluster. Note, while I do have 'results' from my earlier, initial runs of MAGENTA, I was not using full genome-wide data, just top hits, and MAGENTA assumes you are providing it with full genome-wide data. If I were to do this I would need to run Matthew's code on every SNP and then convert the joint BF back to a p-value somehow that's on the same scale as the univariate p-values originally were

##NOTE_E: It would appear that these are the current state of things -- 
#	2010	95_hits	(59_novel)	Matthew_18
#	2013	157_hits (62_novel)	Matthew_75
#	CHECK_0: Look at whether the 18 loci from before are specifically in the 'novel' loci list I received from Sebanti, or can use Annovar to tag them, tag Matthew's 18 and see if the same loci are being tagged at least since rsIDs might be different between those two lists
#
#	/mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs2.R contains the code and results that were used to come up with the 'newhits.vs2.txt' file shown below

#	CHECK_0: Include a line or some code about where the 'previous hits' were included, how the 'RSS0.vs2.txt' file was created which I presume is where the 'previous hits' come into play
#	Line from GLC.MTedits.For2013.vs2.R showing the 'lbfav' cutoff from the least significant 'true' hit
#	newhits = sub[sub$annot==0 & sub$lbfav>4.608818 & sub$nmin>50000,c(4,2:3,7,22:25,30)]
#	
#	Matthew was interested in SNPs whose BFs 'should' have been picked up, so look at SNPs with BFs > 8....guess this is close to 'genome-wide significance' or something as far as a cutoff to use was?
#	Also looked at how many of the new hits have univariate BFs greater than joint BFs and you get the following

~~~
lbf.all.newhits[lbf.uni.newhits < lbf.all.newhits]
[1]  5.113600  9.238360  6.887086  4.339257  4.501904  4.378643  5.930985
[8]  4.371434  5.184725  4.125246  5.819409 17.316921  4.514223  4.770906
[15]  5.259777  5.139357  6.664481  5.397071  4.275684  4.560132  4.154638
[22]  4.857058  5.761283  5.226976  6.581458  4.150895  5.575327  4.252227
[29]  4.292621  4.981670  6.722394 25.754510  5.658189 11.867786  4.854878
[36]  7.653760  6.827153  5.810097  5.007101  5.123809  5.980351  4.851327
[43] 31.807344  4.781310  5.977829  4.754488  6.047135  5.258339  9.998777
[50]  6.890773  5.988893 11.425993  4.190938  6.475299  4.686178  4.280777
[57]  6.004999  5.289372  4.777772 11.061500 36.511011  8.021977  5.313751
[64]  4.773029  5.407697  5.001026  5.908662

lbf.all.newhits[lbf.uni.newhits > lbf.all.newhits]
[1]  4.767202  8.055220 20.992954  5.048607  5.396590  6.225892  4.796701
~~~

~~~
[  mturchin20@spudhead  ~/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/newhits.vs2.txt | wc
     75     825   11000
~~~

#Things moving forward:
#	CHECK_0: Figure out state of those 18 SNPs, which ones are completely reproduced and which ones are not?...get details as if for paper?
#	CHECK_0: Look at the bf > 8 hits and compare univariate vs multivariate BFs -- if no univariate BF is close to multivariate BF makes sense (e.g. multivariate really increased power), but if a single univariate is close to multivariate then this leads to wonder why that SNP wasn't included in the first place. Look at genome browser to see if near any other genes or hits, even though it should have been pruned based on distance already? May also want to create a list of these variants to send to UMich people to see if they have any understanding/reasons to believe why those SNPs didn't make the cut
#	CHECK_0: Manually go through each 75 top hit to get information/details as if for publication, such as eQTL, GenomeBrowser stuff, other  functional information, ANNOVAR
#	CHECK_0: How many of the new hits are functional (e.g. exonic/intronic vs. intergenic?)
#	CHECK_0: Do SUMSTAT on original 157, new 75, then both together? Lookign at both 'do new hits add information to original 157' as well as 'do new 57 represent any novel biology?' This latter question may be better addresssed to while going through each gene manually as well?
#	CHECK_0: Do GIANT results 2010 and 2014 for height, BMI, wasitsizeratio....maybe see from Joel whether there are any other phenotypes that can be used as well?
#	CHECK_0: Reread paper

#	CHECK_0: More down the road -- look for other datasets that may even have just a single release timepoint that have multiple (presumably correlated) phenotypes? Ask Xiang about this potentially to?
#	CHECK_0: Read further into literature, compare Tim's paper ('sharing genetic effects across tissues, assuming independence') to Matthew's paper ('in fact assuming correlation and utilizing that')

#20150521
cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs2.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs3.R





##20150518
##GIANT work

#Note -- renamed the directory /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/... to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/... so conducted the below command on this entire file
#:.,455 s/StephensLab\/GlobalLipids2013/StephensLab\/Multivariate\/GlobalLipids2013/g

mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/ /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/.

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/MainScript1.sh /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/.

#CHECK_0: Figure out github stuff along with these changes

mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013

mkdir /mnt/gluster/data/external_public_supp/GIANT2010
mkdir /mnt/gluster/data/external_public_supp/GIANT2014_5
mkdir /mnt/gluster/data/external_public_supp/GIANT2013

#CHECK_0: Do sex-specific analysis on 2013 or 2014_5 data since it's there/present? 2014_5 doesn't have height stratified....could probably get from Joel if asked? 2013 has all traits stratified by sex

cd /mnt/gluster/data/external_public_supp/GIANT2013

#Sex-specific files from 2013
wget http://www.broadinstitute.org/collaboration/giant/images/3/30/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.txt 
wget http://www.broadinstitute.org/collaboration/giant/images/8/8d/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/5/5c/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/f/fe/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/2/2f/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/e/e5/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/0/08/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/2/24/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/f/f7/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/4/40/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/a/a9/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/9/99/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/3/33/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/e/e7/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/f/fe/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/f/f0/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/f/f8/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.txt
wget http://www.broadinstitute.org/collaboration/giant/images/c/cc/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.txt

cd /mnt/gluster/data/external_public_supp/GIANT2010

wget http://www.broadinstitute.org/collaboration/giant/images/b/b7/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/4/49/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/8/87/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz


#From -- http://www.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files
~~~
BMI (download GZIP)

MD5 (GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt -- 79 MB; 2,471,517 lines) = 38c836542807a3830101bcf48bb34472
If you use these Body Mass Index data, please cite: Speliotes, E.K., Willer, C.J., Berndt, S.I., Monda, K.L., Thorleifsson, G., Jackson, A.U., Allen, H.L., Lindgren, C.M., Luan, J., Magi, R., et al. (2010). Association analyses of 249,796 individuals reveal 18 new loci associated with body mass index. Nat Genet 42, 937-948.
Height (download GZIP)

MD5 (GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt -- 82 MB; 2,469,636 lines) = b51b4c4ff1f03bd33c4b2dfd6b10cb82
If you use these height data, please cite: Lango Allen, H., Estrada, K., Lettre, G., Berndt, S.I., Weedon, M.N., Rivadeneira, F., Willer, C.J., Jackson, A.U., Vedantam, S., Raychaudhuri, S., et al. (2010). Hundreds of variants clustered in genomic loci and biological pathways affect human height. Nature 467, 832-838.
WHRadjBMI (download GZIP)

MD5 (GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt -- 75 MB; 2,483,326 lines) = 8f7e2ca61c33a120db9e7bfe51e3c053
If you use these waist-hip ratio adjusted for BMI data, please cite: Heid, I.M., Jackson, A.U., Randall, J.C., Winkler, T.W., Qi, L., Steinthorsdottir, V., Thorleifsson, G., Zillikens, M.C., Speliotes, E.K., Magi, R., et al. (2010). Meta-analysis identifies 13 new loci associated with waist-hip ratio and reveals sexual dimorphism in the genetic basis of fat distribution. Nat Genet 42, 949-960.
~~~

~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$zcat GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz | wc
2471517 14829102 82685870
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$zcat GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz | wc
2469636 14817816 85061092
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$zcat GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz | wc
2483326 14899956 78051825
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$zcat GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz | head -n 10
MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
rs10 a c 0.0333 0.708 80566
rs1000000 g a 0.6333 0.506 123865
rs10000010 c t 0.425 0.736 123827
rs10000012 c g 0.8083 0.042 123809
rs10000013 c a 0.1667 0.0689 123863
rs10000017 t c 0.2333 0.457 123262
rs1000002 c t 0.475 0.0322 123783
rs10000023 t g 0.5917 0.939 123756
rs10000029 t c 0.975 0.24 103623
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$zcat GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz | head -n 10
MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
rs10 a c 0.0333 0.8826 78380
rs1000000 a g 0.3667 0.1858 133822
rs10000010 t c 0.575 0.8947 132858
rs10000012 c g 0.8083 0.1312 133785
rs10000013 a c 0.8333 0.628 133843
rs10000017 t c 0.2333 0.3073 133174
rs1000002 t c 0.525 0.221 133711
rs10000023 t g 0.5917 0.354 131895
rs10000029 c t 0.025 0.9831 116885
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$zcat GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz | head -n 10
MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
rs10 c a 0.9667 0.42 57031
rs1000000 g a 0.6333 0.55 77168
rs10000010 t c 0.575 0.0029 77152
rs10000012 g c 0.1917 0.99 77117
rs10000013 a c 0.8333 0.89 77167
rs10000017 c t 0.7667 0.64 77166
rs1000002 c t 0.475 0.73 77095
rs10000023 g t 0.4083 0.041 77140
rs10000029 t c 0.975 0.84 61337
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$md5sum GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz
57baef8af0bceb3c09c6084954770afb  GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$md5sum GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz
484dace64c588d524438eae0d02ea2c4  GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/GIANT2010]$md5sum GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz
282989e5678f2e4785c2707349ccff38  GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz
~~~


cd /mnt/gluster/data/external_public_supp/GIANT2014_5

wget http://www.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/f/f0/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz
wget http://www.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz
wget http://www.broadinstitute.org/collaboration/giant/images/e/eb/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/f/f6/GIANT_2015_WHRadjBMI_COMBINED_AllAncestries.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/5/52/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/0/0e/GIANT_2015_HIPadjBMI_COMBINED_AllAncestries.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/7/73/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/3/3f/GIANT_2015_WCadjBMI_COMBINED_AllAncestries.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/5/54/GIANT_2015_WHR_COMBINED_EUR.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/d/d7/GIANT_2015_WHR_COMBINED_AllAncestries.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/e/e4/GIANT_2015_HIP_COMBINED_EUR.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/6/6f/GIANT_2015_HIP_COMBINED_AllAncestries.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/5/57/GIANT_2015_WC_COMBINED_EUR.txt.gz
wget http://www.broadinstitute.org/collaboration/giant/images/e/ea/GIANT_2015_WC_COMBINED_AllAncestries.txt.gz


~~~

~~~

#2010 
#Copying/pasting/downoading/finding (however) list of top hits from each study

mkdir /mnt/lustre/home/mturchin20/Data/GIANT/2010/ 

#Height -- 180 SNPs
#Copy/pasted supplementary table 1 from http://www.nature.com.proxy.uchicago.edu/nature/journal/v467/n7317/full/nature09410.html (http://www.nature.com.proxy.uchicago.edu/nature/journal/v467/n7317/extref/nature09410-s1.pdf) into /mnt/lustre/home/mturchin20/Data/GIANT/2010/LangoAllen2010.SupplTable1.txt

#BMI -- 32 SNPs 
#Copy/pasted table 1 from http://www.nature.com.proxy.uchicago.edu/ng/journal/v42/n11/full/ng.686.html into /mnt/lustre/home/mturchin20/Data/GIANT/2010/Speliotes2010.Table1.txt

#WHRadjBMI -- 14 SNPs
#Copy/pasted table 1 from http://www.nature.com.proxy.uchicago.edu/ng/journal/v42/n11/full/ng.685.html into /mnt/lustre/home/mturchin20/Data/GIANT/2010/Heid2010.Table1.txt
#Manually removed the below three lines from table1
#Further SNPs evaluated in follow up but not achieving genome-wide significance in the combined analysis
#rs2076529       6       32,471,933      BTNL2   C       0.430   2.22 × 10−8     0.041   34,532  0.012   0.011   92,778  3.71 × 10−7     0.020
#rs7081678       10      32,030,629      ZEB1    A       0.085   5.76 × 10−7     0.045   76,270  0.094   0.013   100,527 5.57 × 10−6     0.027

#2014_5
#Copying/pasting/downoading/finding (however) list of top hits from each study

mkdir /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ 

#Height -- 697 SNPs
#Downloaded Supplementary Table 1 from http://www.nature.com/ng/journal/v46/n11/full/ng.3097.html#supplementary-information and saved it as a .txt manually
scp -p mturchin20@wolfy.uchicago.edu:/Users/mturchin20/LabStuff/Data/GIANT/2014_5/ng.3097-S2.txt /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/. 
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.txt | sed 's/\r/\n/g' | perl -lane 'print join("\t", @F);' > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.noCarriageReturn.txt
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.noCarriageReturn.txt | sed 's/"//g' | sed 's/,//g' > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.noCarriageReturn.noQuotes.noCommas.txt

#BMI -- 97 SNPs
#Copy/pasted tables 1 & 2 from http://www.nature.com/nature/journal/v518/n7538/full/nature14177.html into /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.txt
#As well as 'Extended Data Table 2' where I downloaded the jpg nature14177-st2.jpg and used http://www.free-ocr.com/ to extract the text from it
#Note -- for 'Extended Data Table 2' double-checked a few of the entries and corrected positions/rsIDs; feel this is still a better starting point than copying everything manually
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.txt > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.wManualOCREdits.txt

~~~
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.txt | grep rs | wc
     97     979    7660
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.txt | grep rs | awk '{ print $1 }' | sort | uniq -c | wc
     97     194    1774
~~~

#WHRadjBMI -- 49 SNPs (48 SNPs just with EUR; 1 SNP with 'ALL ancestries')
#Downlaoded Supplemetary Table 4 from http://www.nature.com/nature/journal/v518/n7538/full/nature14132.html#supplementary-information and saved it as a .txt manually
scp -p mturchin20@wolfy.uchicago.edu:/Users/mturchin20/LabStuff/Data/GIANT/2014_5/nature14132-s2.ST4.txt /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/.
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.txt | sed 's/\r/\n/g' > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.txt
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.txt | grep WHRadjBMI | grep -vE '2nd|ALL' > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.Eur.txt
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.txt | grep WHRadjBMI | grep -vE '2nd' > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.ALL.txt
#CHECK_1: FIgure out the correct combinations for including the other waist traits -- do the '2nd' hits count if they are just conditional analyses? Solution -- just using all SNPs regardless of GWAS only, GWAS+metabochip, all vs Euro ancestries, all sex vs. male/female only, or joint/conditional. Not sure if can parse out all these factors for every trait (e.g. height results from GATC COJO directly and not sure if I can get 'non-conditional' hits so specifically like I can for WHRadjBMI or BMI)

~~~
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.Eur.txt | grep rs | awk '{ print $1 }' | sort | uniq -c | wc
     48      96     871
~~~

#2013
#Downloaded TableS2 from from http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003500 and saved it as a .txt manually
mkdir /mnt/lustre/home/mturchin20/Data/GIANT/2013
scp -p mturchin20@wolfy.uchicago.edu:/Users/mturchin20/LabStuff/Data/GIANT/2013/journal.pgen.1003500.s007.txt /mnt/lustre/home/mturchin20/Data/GIANT/2013/.

cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.txt | sed 's/\r/\n/g' > /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt

#46 unique SNPs across either or both sexes with pval < 5e-8 (from Figure 3)
~~~
[  mturchin20@spudhead  ~/Software]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt | awk '{ if (($25 < .00000005) || ($26 < .00000005)) { print $1 } }' | sort | uniq | wc
     46      46     464
[  mturchin20@spudhead  ~/Software]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt | awk '{ if ($25 < .00000005 ) { print $1 } }' | wc
     15      15     151
[  mturchin20@spudhead  ~/Software]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt | awk '{ if ($26 < .00000005 ) { print $1 } }' | wc
     38      38     385
~~~

#Male -- 15 SNPs
cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt | awk '{ if ($25 < .00000005 ) { print $1 } }' > /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Male.rsIDs.txt
cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt | awk '{ if ($25 < .00000005 ) { print $1, "\t", $5, "\t", $6 } }' > /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Male.MarkerChrBP.txt 

#Female -- 38 SNPs
cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt | awk '{ if ($26 < .00000005 ) { print $1 } }' > /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Female.rsIDs.txt
cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.txt | awk '{ if ($26 < .00000005 ) { print $1, "\t", $5, "\t", $6 } }' > /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Female.MarkerChrBP.txt




#2010
#Conducting multivariate analyses

mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/process.MTedits.For2013.vs2.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/process.MTedits.ForGIANT2010.vs1.R

#Converted process.MTedits.For2013.vs2.R to GIANT phenotypes, moving commands to height, BMI, WHRadjBMI instead of tg, tc, ldl and hdl
#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/process.MTedits.ForGIANT2010.vs1.R into an interactive sessions in $PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1

#mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/RData.20150519
#Redid /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/process.MTedits.ForGIANT2010.vs1.R to switch dt2 columns so that order of phenotypes was Z.height Z.BMI Z.WHRadjBMI, not Z.WHRadjBMI Z.height Z.BMI
mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/RData.20150521

gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/RData.20150519
gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/RData.20150521

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/GLC.MTedits.For2013.vs2.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GLC.MTedits.ForGIANT2010.vs1.R

#For GIANT2010.dtlesssignif.annot.vs1.txt.gz will need (except in GIANT vals)
#posHg18 chr pos snp a1 a2 maf beta_LDL se_LDL n_LDL beta_HDL se_HDL n_HDL beta_TG se_TG n_TG beta_TC se_TC n_TC annot gene Z.tg Z.tc Z.ldl Z.hdl mvstat mvp unip (from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz)
#chr pos snp a1 a2 maf beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI annot gene Z.height Z.BMI Z.WHRadjBMI mvstat mvp unip

~~~
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.vs2.txt | head -n 10
rs,Z.tg,Z.tc,Z.ldl,Z.hdl,mvstat,mvp,unip
chr10:101902054,2.45906343998997,3.65960938616346,0.54246240431255,4.84175096013916,43.8856399660704,8.16899883851383,5.89042145309562
chr10:101902184,2.04241886517916,3.31897041855539,0.741979097659704,4.28996658928395,33.5923747854106,6.04416453147192,4.74787544749436
chr10:101950956,1.89987138538418,4.80442095262784,2.98132929795538,3.49559840808192,36.4681968880102,6.63489650218934,5.80910828307783
chr10:101957202,2.3172521539636,4.44716747735029,2.65105825613512,3.51059798302414,34.6265665431121,6.25629717529947,5.0604308313441
chr10:101966491,2.38390551204578,4.24376554263135,2.1138275382868,3.50881965277072,33.4617468909834,6.01739588476946,4.65797231191253
chr10:101989736,1.89426382584585,4.15303118888146,2.03680088096434,3.87804325071873,33.8172794310246,6.0902663254914,4.48399376961395
chr10:101990691,1.97401532860328,4.24234125347809,2.21939299224744,3.75592849665782,33.7012328521064,6.06647651253057,4.65521487736734
chr10:101996419,1.96555531599126,4.668233848697,2.92978188856575,3.38349312248405,34.6728122799936,6.26579130288765,5.51741223047323
chr10:101998159,2.16811418290611,4.79021825506647,3.0723226578465,3.4849657712978,37.0775367184447,6.7603876155226,5.77832500292923
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.txt | head -n 10
rs,Z.height,Z.BMI,Z.WHRadjBMI,mvstat,mvp,unip
rs10006067,7.08776007395249,1.65020996204492,1.8119106729526,55.8346783977573,10.663188034222,11.8655041441653
rs1000972,7.98723399811821,3.43709008602989,3.37805680205288,86.1131635514716,17.0552250763639,14.8601209135988
rs10011200,6.22294771858863,0.186567181836519,1.59819313992282,41.1157838111925,7.5945727241859,9.31166918188773
rs10013023,6.04641969222953,1.28727056310794,0.371856089385075,38.2024656839996,6.99233732659037,8.82944494147879
rs10016839,5.71547189846515,1.22387337228946,0.994457883209753,34.9475904413173,6.32221671701731,7.96098267800259
rs10017744,6.15675623449777,1.34075503369022,0.495850347347453,39.7751132103379,7.31716845232464,9.12930354201075
rs10020593,5.7374722838198,1.22652812003661,0.994457883209753,35.2050664771305,6.3751109638301,8.01727661233145
rs10023833,5.974709888153,0.0401168101849681,2.51214432793046,41.7298228825345,7.72176818166663,8.63732907027433
rs10028610,5.55041024777321,1.50237611995585,0.568051498338983,33.208638751781,5.96554504181686,7.54515513999149
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/jointGwasMc_AllPheno.Annot.txt.gz | head -n 10
posHg18 chr pos snp a1 a2 maf beta_LDL se_LDL n_LDL beta_HDL se_HDL n_HDL beta_TG se_TG n_TG beta_TC se_TC n_TC annot gene
chr6:54599910   6       54491951        rs1280728       t       a       0.0     0.0108  0.0081  89888.00        0.0107  0.0074  94311.00        0.0139  0.0073  91013.00        0.0125  0.0080  94595.00        0       NA
chr9:3746564    9       3756564 rs10814541      c       t       0.2982  0.0011  0.0064  89888.00        0.0036  0.0058  94311.00        0.0119  0.0057  91013.00        0.0025  0.0062  94595.00        0       NA
chr5:133405061  5       133377162       rs4958172       g       a       0.1187  0       0.0089  89817.00        0.0072  0.0082  94239.00        0.0029  0.0080  90938.00        0.0029  0.0087  94520.00        0       NA
chr12:69055675  12      70769408        rs2053102       g       c       0.0     0.0069  0.0075  89888.00        0.0075  0.0069  94311.00        0.0076  0.0068  91013.00        0.0009  0.0073  94595.00        0       NA
chr5:155873914  5       155941336       rs256837        c       a       0.2005  0.0016  0.0066  88199.00        0.0002  0.0061  92620.00        0.0013  0.0060  89323.00        0.0031  0.0065  92905.00        2       NA
chr3:60337349   3       60362309        rs1569354       t       a       0.0     0.0021  0.0054  89888.00        0.0060  0.0050  94311.00        0.0035  0.0048  91013.00        0.0025  0.0053  94595.00        0       NA
chr4:65034366   4       65351771        rs6851961       a       g       0.4631  0.0003  0.0053  89877.00        0.0041  0.0049  94300.00        0.0025  0.0048  91002.00        0.0022  0.0052  94584.00        0       NA
chr7:62372905   7       62735470        rs11560217      c       a       0.03034 0.0126  0.0161  83513.00        0.0236  0.0147  88983.00        0.0041  0.0148  81387.00        0.0008  0.0157  88090.00        0       NA
chr17:6058839   17      6118115 rs17803227      a       g       0.1016  0.0046  0.0094  88422.00        0.0055  0.0086  92809.00        0.0081  0.0084  89474.00        0.0013  0.0092  93056.00        0       NA
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/dtlesssignif.annot.txt | wc
8066  209716 1794390
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/dtlesssignif.annot.txt | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
3724   96824  830117
[  mturchin20@spudhead  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz | wc
12463  348964 3301334
[  mturchin20@spudhead  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
8859  248052 2346362
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz | wc
   8842  203366 1641132
   [  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
   1097   25231  202912
~~~

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.CreateFullAnnotFile.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py

cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/LangoAllen2010.SupplTable1.txt | grep rs | awk '{ print $1, "\t", $2, "\t", $3 }' > /mnt/lustre/home/mturchin20/Data/GIANT/2010/LangoAllen2010.SupplTable1.rsIDs.MarkerChrBP.txt

cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/Speliotes2010.Table1.txt | grep rs | perl -ane 'print $F[0], "\t"; for (my $i = 0; $i <= $#F; $i++) { if ($F[$i] =~ m/\d,\d/) { print $F[$i-1], "\t", $F[$i], "\t"; } } print "\n";' | awk '{ print $1, "\t", $2, "\t", $3 }' | sed 's/,//g' > /mnt/lustre/home/mturchin20/Data/GIANT/2010/Speliotes2010.Table1.rsIDs.MarkerChrBP.txt

cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/Heid2010.Table1.txt | grep rs | awk '{ print $1, "\t", $2, "\t", $3 }' | sed 's/,//g' > /mnt/lustre/home/mturchin20/Data/GIANT/2010/Heid2010.Table1.rsIDs.MarkerChrBP.txt

cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/LangoAllen2010.SupplTable1.rsIDs.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/GIANT/2010/Speliotes2010.Table1.rsIDs.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/GIANT/2010/Heid2010.Table1.rsIDs.MarkerChrBP.txt > /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt

~~~
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | wc
226     678    5794
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | sort | uniq -c | wc
226     904    7602
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | awk '{ print $1 }' | sort | uniq -c | wc
226     452    4074
~~~

#Need to get basepair position for the GIANT GWAS summary files -- only have rsIDs for all of them. Guessing/hoping dbSNP main vcf has this information
mkdir /mnt/lustre/home/mturchin20/Data/dbSNP 
cd /mnt/lustre/home/mturchin20/Data/dbSNP 
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/All_20150415.vcf.gz 
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh37p13/VCF/All_20150415.vcf.gz.tbi
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh38/VCF/All_20150416.vcf.gz
mv /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150416.vcf.gz /mnt/lustre/home/mturchin20/Data/dbSNP/dbSNP_human_9606_b142_GRCh38_All_20150416.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b142_GRCh38/VCF/All_20150416.vcf.gz.tbi
mv /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150416.vcf.gz.tbi /mnt/lustre/home/mturchin20/Data/dbSNP/dbSNP_human_9606_b142_GRCh38_All_20150416.vcf.gz.tbi

~~~
[  mturchin20@spudhead  ~/Data/dbSNP]$zcat All_20150415.vcf.gz | head -n 80 | tail -n 10
1       10250   rs199706086     A       C       .       .       RS=199706086;RSPOS=10250;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x050000020005000002000100;WGT=1;VC=SNV;R5;ASP
1       10254   rs140194106     TA      T       .       .       RS=140194106;RSPOS=10255;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000020005000002000200;WGT=1;VC=DIV;R5;ASP
1       10257   rs111200574     A       C       .       .       RS=111200574;RSPOS=10257;dbSNPBuildID=132;SSR=0;SAO=0;VP=0x050100020005000102000100;WGT=1;VC=SNV;SLO;R5;ASP;GNO
1       10259   rs200940095     C       A       .       .       RS=200940095;RSPOS=10259;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x050000020005000002000100;WGT=1;VC=SNV;R5;ASP
1       10291   rs145427775     C       T       .       .       RS=145427775;RSPOS=10291;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000020005000002000100;WGT=1;VC=SNV;R5;ASP
1       10327   rs112750067     T       C       .       .       RS=112750067;RSPOS=10327;dbSNPBuildID=132;SSR=0;SAO=0;VP=0x050000020005000002000100;WGT=1;VC=SNV;R5;ASP
1       10328   rs201106462     AACCCCTAACCCTAACCCTAACCCT       A       .       .       RS=201106462;RSPOS=10329;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x050000020005000002000200;WGT=1;VC=DIV;R5;ASP
1       10329   rs150969722     AC      A       .       .       RS=150969722;RSPOS=10330;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000020005000002000200;WGT=1;VC=DIV;R5;ASP
1       10352   rs145072688     T       TA      .       .       RS=145072688;RSPOS=10353;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000020015000002000200;WGT=1;VC=DIV;R5;OTH;ASP;CAF=0.5625,0.4375;COMMON=1
1       10352   rs555500075     T       TA      .       .       RS=555500075;RSPOS=10352;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020015140024000200;WGT=1;VC=DIV;R5;OTH;ASP;VLD;KGPhase3;CAF=0.5625,0.4375;COMMON=1
[  mturchin20@spudhead  ~/Data/dbSNP]$zcat All_20150415.vcf.gz | awk '{ print $1, "\t", $2, "\t", $3 }' | wc
111198367 333595084 3048046196
~~~

cp -p /mnt/lustre/home/mturchin20/Scripts/Python/2fileTemplate.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py
cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs2.py


python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz

#Ordering the two files this way due to large 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs2.py --file1 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz --file2 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs2.txt.gz
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs2.py --file1 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz --file2 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs2.txt.gz
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs2.py --file1 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz --file2 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs2.txt.gz

~~~
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | wc
2469636 17287452 113790927
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs2.txt.gz | wc
2465677 17259739 113641249
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | wc
2469636 17287452 113790924
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | wc
2483326 17383282 106629662
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | wc
2483326 17383282 106629659
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | wc
2471517 17300619 111437629
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | wc
2471517 17300619 111437626
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz | wc             2471517 14829102 82685870
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz | wc            2483326 14899956 78051825
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz | wc         2469636 14817816 85061092
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | grep NA | wc
3972   27804  145898
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | wc
3971   27797  145841
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | grep Marker
MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  p       N       NA
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep Marker
MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  p       N       ChrBP
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | wc
3958   27706  149618
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | grep NA | wc
3959   27713  149675
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | wc
39972  279804 1229461
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | grep NA | wc
39973  279811 1229518
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.vs3.txt.gz | grep NA | grep -v rs
MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  p       N       NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | grep -v rs
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | grep -v rs
~~~

scp -p mturchin20@wolfy.uchicago.edu:/Users/mturchin20/LabStuff/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/. 

~~~
[  mturchin20@spudhead  ~/Data/dbSNP]$zcat /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz | wc
18833532 188335320 1426536856
~~~

#Redoing the ChrBP annotation steps with dbSNP130 now
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &

zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc

~~~
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2469636 17287452 114167633
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
10809   75663  865971
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1981   13867   75103
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2471517 17300619 111814122
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
10805   75635  848551
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1980   13860   73049
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2483326 17383282 106990904
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
10446   73122  797793
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
37936  265552 1159088
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
0       0       0
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
0       0       0
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
0       0       0
~~~

#Going to be conservative and allow a rsID to exist in multiple locations if there are multiple locations for now from dbSNP130 (hence the ',' marks). Will include code that accounts for this in the 'annotation' script as well as to identify whether there is no 'rsID' for the marker and just chr#:##### which a few of the files contain (though I think most are 100% rs##### -- which is a bit interesting/odd since thought the same platform was being used across all phenotypes....but I guess there might be differences due to imputation quality/results or something for a given dataset/phenotype too?....different cohorts in total are being used across each dataest too)

#20150520 -- Realized the 2010 GIANT data doesn't have betas or SEs......so can't use them for this yet. Might be able to get that data from Joel if ask but for now can only do the 2014_5 data and the sex-stratified 2013 data
#20150520 -- 30 minutes later I realized/remembered the method only uses Zscores, maf and n, not betas or SEs. For some reason we just keep that information, when present, in the prep files but it's not actually used by the method. So....I can go ahead with this like I originally planned. yay. \o/

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz 
#python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.Trbl1.gz
#python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.Trbl2.gz
#python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.Trbl3.gz

#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz
#MarkerName a1 a2 maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI

#Troubleshooting /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py stuff
~~~
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep -E 'rs2154319|rs10863936|rs11205277|rs4665736' | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
Happening1
rs11205277      g       a       NA      2       NA      NA      127114  NA      NA      127114
Happening1
rs10863936      g       a       NA      2       NA      NA      133820  NA      NA      133820
Happening1
rs2154319       c       t       NA      2       NA      NA      133480  NA      NA      133480
Happening1
rs4665736       t       c       0.4667  2       NA      NA      133780  NA      NA      133780
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.Trbl1.gz | wc
2469635 34574890 132103187
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.Trbl2.gz | wc
2469635 34574890 141891923
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.Trbl1.gz | head -n 10
rs999721        a       g       0.4667  0       NA      NA      133732  NA      NA      NA      NA      NA      NA
rs9821657       g       a       0.175   0       NA      NA      133582  NA      NA      NA      NA      NA      NA
rs1919329       t       g       0.4833  0       NA      NA      133834  NA      NA      NA      NA      NA      NA
rs6444035       c       t       0.0667  0       NA      NA      133803  NA      NA      NA      NA      NA      NA
rs3957240       c       t       0.05    0       NA      NA      122588  NA      NA      NA      NA      NA      NA
rs1919324       c       t       0.1417  0       NA      NA      133796  NA      NA      NA      NA      NA      NA
rs2739330       c       t       0.3833  0       NA      NA      133781  NA      NA      NA      NA      NA      NA
rs973978        a       g       NA      0       NA      NA      133825  NA      NA      NA      NA      NA      NA
rs17597444      t       g       0.25    0       NA      NA      133772  NA      NA      NA      NA      NA      NA
rs17597445      t       c       NA      0       NA      NA      133810  NA      NA      NA      NA      NA      NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.Trbl2.gz | head -n 10
rs999721        a       g       0.4667  0       NA      NA      133732  NA      NA      133732  NA      NA      NA
rs9821657       g       a       0.175   0       NA      NA      133582  NA      NA      133582  NA      NA      NA
rs1919329       t       g       0.4833  0       NA      NA      133834  NA      NA      133834  NA      NA      NA
rs6444035       c       t       0.0667  0       NA      NA      133803  NA      NA      133803  NA      NA      NA
rs3957240       c       t       0.05    0       NA      NA      122588  NA      NA      122588  NA      NA      NA
rs1919324       c       t       0.1417  0       NA      NA      133796  NA      NA      133796  NA      NA      NA
rs2739330       c       t       0.3833  0       NA      NA      133781  NA      NA      133781  NA      NA      NA
rs973978        a       g       NA      0       NA      NA      133825  NA      NA      133825  NA      NA      NA
rs17597444      t       g       0.25    0       NA      NA      133772  NA      NA      133772  NA      NA      NA
rs17597445      t       c       NA      0       NA      NA      133810  NA      NA      133810  NA      NA      NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
2516306 35228284 151538382
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | head -n 10
rs999721        a       g       0.4667  0       NA      NA      133732  NA      NA      123835  NA      NA      77154
rs9821657       g       a       0.175   0       NA      NA      133582  NA      NA      123732  NA      NA      77043
rs1919329       t       g       0.4833  0       NA      NA      133834  NA      NA      123862  NA      NA      77166
rs6444035       c       t       0.0667  0       NA      NA      133803  NA      NA      123814  NA      NA      77155
rs3957240       c       t       0.05    0       NA      NA      122588  NA      NA      114400  NA      NA      68863
rs1919324       c       t       0.1417  0       NA      NA      133796  NA      NA      123864  NA      NA      77168
rs2739330       c       t       0.3833  0       NA      NA      133781  NA      NA      123892  NA      NA      77209
rs973978        a       g       NA      0       NA      NA      133825  NA      NA      123846  NA      NA      77155
rs17597444      t       g       0.25    0       NA      NA      133772  NA      NA      123864  NA      NA      77168
rs17597445      t       c       NA      0       NA      NA      133810  NA      NA      123867  NA      NA      77169
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "NA") { print $0 } } ' | wc
 225095 3151330 12559530
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrga
nization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
rs11205277      g       a       NA      1       NA      NA      127114  NA      NA      127114  NA      NA      NA
rs10863936      g       a       NA      1       NA      NA      133820  NA      NA      133820  NA      NA      NA
rs2154319       c       t       NA      1       NA      NA      133480  NA      NA      133480  NA      NA      NA
rs4665736       t       c       0.4667  1       NA      NA      133780  NA      NA      133780  NA      NA      NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c
2292802 0
226 1
223278 2
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | wc
226     678    5794
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
rs11205277      g       a       1_148159495     NA      1       NA      NA      NA      127114  NA      NA      NA      127114  NA      NA      NA      NA
rs10863936      g       a       1_210304420     NA      1       NA      NA      NA      133820  NA      NA      NA      133820  NA      NA      NA      NA
rs2154319       c       t       1_41518356      NA      1       NA      NA      NA      133480  NA      NA      NA      133480  NA      NA      NA      NA
rs4665736       t       c       2_25041102      0.4667  1       0.4667  NA      NA      133780  0.4667  NA      NA      133780  NA      NA      NA      NA
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c
2292802 0
    226 1
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5, "\t", $7, "\t", $11, "\t", $15 }' | sort | uniq -c                                                                                             
.
.
.
223278 2
    14 0.2      NA      NA      0.2
    30 0.2      NA      NA      NA
 33331 0.3      0.3     0.3     0.3
   442 0.3      0.3     0.3     NA
     1 0.3      0.3     NA      0.3
     5 0.3      0.3     NA      NA
 33582 0.3083   0.3083          0.3083          0.3083
   450 0.3083   0.3083          0.3083          NA
.
.
.
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($5 == "NA") { print $5, "\t", $7, "\t", $11, "\t", $15 } }' | sort | uniq -c
     2 NA       NA      NA      0.0
225093 NA       NA      NA      NA
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if (($5 == "0.2") && ($7 == "NA")) { print $0 } } ' | head -n 10
rs3006263       c       g       10_26198595     0.2     0       NA      NA      NA      NA      NA      NA      0.2     NA      NA      41946   NA      NA
rs9646391       c       t       17_69523401     0.2     0       NA      NA      NA      NA      NA      NA      0.2     NA      NA      43049   NA      NA
rs4789011       t       c       17_69345766     0.2     0       NA      NA      NA      NA      NA      NA      0.2     NA      NA      43206   NA      NA
rs11124487      a       g       2_35975673      0.2     0       NA      NA      NA      NA      NA      NA      0.2     NA      NA      41710   NA      NA
rs267154        t       c       3_76963429      0.2     0       NA      NA      NA      NA      NA      NA      0.2     NA      NA      41885   NA      NA
rs4085550       g       c       13_22416770     0.2     0       NA      NA      NA      NA      NA      NA      0.2     NA      NA      46353   NA      NA
rs1593477       c       t       18_2291963      0.2     0       NA      NA      NA      0.2     NA      NA      62647   NA      0.2     NA      NA      48836
rs6051233       g       c       20_2622734      0.2     0       NA      NA      NA      0.2     NA      NA      69572   NA      0.2     NA      NA      43449
rs11768737      c       a       7_41341790      0.2     0       NA      NA      NA      NA      NA      NA      0.2     NA      NA      40685   NA      NA
rs17406991      c       t       7_78515066      0.2     0       NA      NA      NA      0.2     NA      NA      62027   NA      0.2     NA      NA      46607
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
2516306 17
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep rs4789011
rs4789011       t       c       0.2     0.21    43206   17_69345766
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_WHRadjBMI_Heid2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep -E 'rs1593477|rs4789011|rs3006263|rs17406991' | gzip > /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4otherSNPs.txt.gz
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4otherSNPs.txt.gz
rs1593477       c       t       0.2     0.24    48836   18_2291963
rs1593477       c       t       0.2     0.357   62647   18_2291963
rs17406991      c       t       0.8     0.446   62027   7_78515066
rs17406991      t       c       0.2     0.39    46607   7_78515066
rs3006263       c       g       0.2     0.32    41946   10_26198595
rs4789011       t       c       0.2     0.21    43206   17_69345766
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4otherSNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4otherSNPs.txt.gz
/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4otherSNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2010/GIANT_2010_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4otherSNPs.txt.gz
Error2a -- line2[0] present more than once in a particular file (0: rs1593477,c,t,0.2,0.357,62647,18_2291963)
Error2a -- line2[0] present more than once in a particular file (0: rs17406991,t,c,0.2,0.39,46607,7_78515066)
Error2a -- line2[0] present more than once in a particular file (1: rs1593477,c,t,0.2,0.357,62647,18_2291963)
Error2a -- line2[0] present more than once in a particular file (1: rs17406991,t,c,0.2,0.39,46607,7_78515066)
rs1593477       c       t       18_2291963      0.2     0       0.2     NA      NA      62647   0.2     NA      NA      62647   NA      NA      NA      NA
rs4789011       t       c       17_69345766     0.2     0       0.2     NA      NA      43206   0.2     NA      NA      43206   NA      NA      NA      NA
rs3006263       c       g       10_26198595     0.2     0       0.2     NA      NA      41946   0.2     NA      NA      41946   NA      NA      NA      NA
rs17406991      c       t       7_78515066      0.2     0       0.2     NA      NA      46607   0.2     NA      NA      46607   NA      NA      NA      NA
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
2516306 45293508 227795829
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
2516306 17
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($5 == "NA") { print $5, "\t", $7, "\t", $11, "\t", $15 } }' | sort | uniq -c
 225093 NA       NA      NA      NA
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $8, "\t", $9, "\t", $12, "\t", $13, "\t", $16, "\t", $17 }' | sort | uniq -c
2516306 NA   NA      NA      NA      NA      NA
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'if (($F[9] =~ m/\D+/) && ($F[9] ne "NA")) { print join("\t", @F); } ' | sort | uniq -c
1 rs796188        a       c       13_88518668     0.0167  0       0.0167  NA      NA      1e+05   0.0167  NA      NA      92198   0.0167  NA      NA      56914
1 rs9842319       t       g       3_144835684     0.0333  0       0.0333  NA      NA      1e+05   0.0333  NA      NA      95468   0.0333  NA      NA      54663
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'if (($F[13] =~ m/\D+/) && ($F[13] ne "NA")) { print join("\t", @F); } ' | sort | uniq -c
1 rs11244235      g       c       9_132864356     0.025   2       0.025   NA      NA      108437  0.025   NA      NA      1e+05   0.025   NA      NA      62450
1 rs12917986      c       g       16_72309046     0.2083  0       0.2083  NA      NA      109441  0.2083  NA      NA      1e+05   0.2083  NA      NA      61723
1 rs17599562      g       t       10_122917026    0.05    0       0.05    NA      NA      95129   0.05    NA      NA      1e+05   0.05    NA      NA      66847
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1]$zcat GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'if (($F[17] =~ m/\D+/) && ($F[17] ne "NA")) { print join("\t", @F); } ' | sort | uniq -c
~~~

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/DataOrganization.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py

#python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.txt --file2 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.txt --file2 <(zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'my $chr; my $pos; if ($F[3] eq "NA") { $chr = "NA"; $pos = "NA"; } else { $chr = ((split(/_/, $F[3]))[0]); $pos = ((split(/_/, $F[3]))[1]); } splice(@F, 1, 0, ($chr, $pos)); splice(@F, 8, 1); splice(@F, 11, 1); splice(@F, 14, 1); print join("\t", @F);') | \
gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz
#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz
#snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI Z.height Z.BMI Z.WHRadjBMI mvstat mvp unip

#Converted GLC.MTedits.For2013.vs2.R to GIANT phenotypes, moving commands to height, BMI, WHRadjBMI instead of tg, tc, ldl and hdl
#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GLC.MTedits.ForGIANT2010.vs1.R into an interactive sessions in $PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1

cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | grep -f <(cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | awk '{ print $1 }') 

~~~
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | wc                                20     160     861
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | awk '{ print $1 }' | grep -v -f <(cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | grep -f <(cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | awk '{ print $1 }') | awk '{ print $1 }') | wc
14      14     138
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | grep -f <(cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | awk '{ print $1 }') | wc
6      18     155
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | grep -f <(cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | awk '{ print $1 }')
rs34651          5       72179761
rs12204421       6       33736841
rs648831         6       81012927
rs17783015       12      88755517
rs11835818       12      120979192
rs1809889        12      123367179
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | awk '{ print $1 }' | grep -v -f <(cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | grep -f <(cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.newtophits.vs1.txt | awk '{ print $1 }') | awk '{ print $1 }')
snp
rs17016663
rs10805383
rs7601531
rs2025151
rs6824258
rs4735692
rs7614120
rs7081678
rs12534698
rs2390312
rs389883
rs10040888
rs10822129
~~~

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/RData.GIANT2010.GLC.20150522

gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/RData.GIANT2010.GLC.20150522





#2014_5
#Conducting multivariate analyses

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/process.MTedits.ForGIANT2010.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/process.MTedits.ForGIANT2014_5.vs1.R

#A bit anachronistic here -- editing the 'dbSNP130' versions of the GWAS file to correct the 'headers' so they can be properly used by R and other processes

cp -p /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz
cp -p /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz
cp -p /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz
cp -p /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz
cp -p /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz
cp -p /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz

#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/process.MTedits.ForGIANT2014_5.vs1.R into an interactive sessions in $PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/RData.GIANT2014_5.process.20150522
gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/RData.GIANT2014_5.process.20150522

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GLC.MTedits.ForGIANT2010.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.Orig3.vs1.R

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.CreateFullAnnotFile.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py

cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.noCarriageReturn.noQuotes.noCommas.txt | grep rs | awk '{ print $1, "\t", $2, "\t", $3 }' > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.noCarriageReturn.noQuotes.noCommas.rsIDs.MarkerChrBP.txt

cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.wManualOCREdits.txt | grep rs | awk '{ print $1, "\t", $2 }' | sed 's/,//g' | perl -lane 'print $F[0], "\t", ((split(/:/, $F[1]))[0]), "\t", ((split(/:/, $F[1]))[1]);' > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.wManualOCREdits.rsIDs.MarkerChrBP.txt

~~~
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.wManualOCREdits.txt| grep rs | awk '{ print $1, "\t", $2 }' | sed 's/,//g' | perl -lane 'print $F[0], "\t", ((split(/:/, $F[1]))[0]), "\t", ((split(/:/, $F[1]))[1]);' | wc                                                97     291    2130
~~~

#cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.Eur.txt | grep rs | awk '{ print $1, "\t", $3, "\t", $4 }' | sed 's/"//g' | sed 's/,//g' | grep -v This > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.Eur.rsIDs.MarkerChrBP.txt
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.txt | grep WHRadjBMI | grep rs | awk '{ print $1, "\t", $3, "\t", $4 }' | sed 's/"//g' | sed 's/,//g' | grep -v This > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.rsIDs.MarkerChrBP.txt
cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.txt | grep rs | awk '{ print $1, "\t", $3, "\t", $4 }' | sed 's/"//g' | sed 's/,//g' | grep -v This > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.rsIDs.MarkerChrBP.txt

~~~
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.rsIDs.MarkerChrBP.txt | wc
    68     204    1760
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.rsIDs.MarkerChrBP.txt | wc
    87     261    2245
~~~

cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.noCarriageReturn.noQuotes.noCommas.rsIDs.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.wManualOCREdits.rsIDs.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.WHRadjBMI.rsIDs.MarkerChrBP.txt > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt

cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/ng.3097-S2.noCarriageReturn.noQuotes.noCommas.rsIDs.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.wManualOCREdits.rsIDs.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/nature14132-s2.ST4.noCarriageReturn.rsIDs.MarkerChrBP.txt > /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.AllPheno.MarkerChrBP.txt

~~~
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | wc
862    2586   21825
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.AllPheno.MarkerChrBP.txt | wc
881    2643   22310
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | awk '{ print $1 }' | sort | uniq | wc
862     862    8731
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.AllPheno.MarkerChrBP.txt | awk '{ print $1 }' | sort | uniq | wc
881     881    8921
~~~

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz &

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/dbSNP/dbSNP_human_9606_b142_GRCh38_All_20150416.vcf.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh38.ChrBP.txt.gz

~~~
[  mturchin20@spudhead  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz | wc
2555511 20444088 125133265
[  mturchin20@spudhead  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.gz | wc
2554638 20437104 125031676
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | wc
2550859 22957731 148851105
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz | wc      2550859 20406872 119190473
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.gz | wc
2554638 22991742 154763475
[  mturchin20@bigmem02  ~/Data/dbSNP]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.gz | wc                                   2554638 20437104 125031676
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | wc
wc: standard input:8: Invalid or incomplete multibyte or wide character
2542440 22881925 147378256
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz | wc
wc: standard input:8: Invalid or incomplete multibyte or wide character
2542439 20339485 118046720
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | wc
8176   73584  414776
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.gz | grep NA | wc
56907  512163 3248718
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | awk '{ print $9 }' | sort | tail -n 508853 | grep NA | wc
  32523   32523   97569
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | wc
wc: standard input:8: Invalid or incomplete multibyte or wide character
2540934 27950223 177978334
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz | wc
wc: standard input:8: Invalid or incomplete multibyte or wide character
2540933 25409289 148664320
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | awk '{ print $11 }' | sort | tail -n 508853 | grep NA | wc
  32510   32510   97530
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | wc
wc: standard input:8: Invalid or incomplete multibyte or wide character
2546082 28006851 177414406
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | awk '{ print $11 }' | sort | tail -n 508853 | grep NA | wc
  33608   33608  100824
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | grep -v rs
SNP_A-2097957   A       C       0.681   0.012   0.0033  0.00040 231091  NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.gz | grep NA | grep -v rs
SNP     A1      A2      Freq1.Hapmap    b       se      p       N       NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | grep NA | grep -v rs
Binary file (standard input) matches
~~~


#Realized, because /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz has ChrBP information in it, that I am mapping to the wrong reference database (dbSNP142 I guess?). Comparing the 'GWAS hit' files, which do have ChrBP info, with this WCadjBMI file to see if they are mapping to the same region, for both 2010 and 2014_5 'GWAS hits'. It looks like everything was mapped to the same build, 2010 to 2014_5, which is NCBI36/hg18 & dbSNP130. So able to download the dbSNP130 list of variants from UCSC genome browser, the commands for which are at the bottom of this chunk of labcode/results
~~~
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | tail -n 20
rs2890652        2       142676401
rs1555543        1       96717385
rs4771122        13      26918180
rs4836133        5       124360002
rs4929949        11      8561169
rs206936         6       34410847
rs9491696        6       127494332
rs6905288        6       43866851
rs984222         1       119305366
rs1055144        7       25837634
rs10195252       2       165221337
rs4846567        1       217817340
rs1011731        1       170613171
rs718314         12      26344550
rs1294421        6       6688148
rs1443512        12      52628951
rs6795735        3       64680405
rs4823006        22      27781671
rs6784615        3       52481466
rs6861681        5       173295064
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | perl -lane 'if ($F[0] =~ m/rs4823006/) { print join("\t", @F); }'
rs4823006       22      27781671        A       G       0.5333  0.018   0.0034  2.2e-07 231058  22_29451671
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | perl -lane 'if ($F[0] =~ m/rs1555543/) { print join("\t", @F); }'
rs1555543       1       96717385        A       C       0.425   -0.0026 0.0043  0.55    153380  1_96944797
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$cat /mnt/lustre/home/mturchin20/Data/GIANT/2010/GIANT2010.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt | head -n 10
rs425277         1       2059032
rs2284746        1       17179262
rs1738475        1       23409478
rs4601530        1       24916698
rs7532866        1       26614131
rs2154319        1       41518357
rs17391694       1       78396214
rs6699417        1       88896031
rs10874746       1       93096559
rs9428104        1       118657110
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | perl -lane 'if ($F[0] =~ m/rs7532866/) { print join("\t", @F); }'
rs7532866       1       26614131        A       G       NA      0       0.0036  0.99    229730  1_26741544
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wdbSNP_human_9606_b142_GRCh37p13.ChrBP.txt.gz | perl -lane 'if ($F[0] =~ m/rs10874746/) { print join("\t", @F); }'
rs10874746      1       93096559        T       C       0.3417  -0.0084 0.0035  0.016   231211  1_93323971
~~~

#Redoing the ChrBP annotation steps with dbSNP130 now

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.txt.gz | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz &

zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc

~~~
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2550859 22957731 149252489
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
11574  104166 1060266
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
5913   53217  306525
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
1       9      55
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | wc
2554638 22991742 155169515
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep , | wc
11705  105345 1094041
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep NA | wc
54778  493002 3162501
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep NA | grep -v rs | wc
1       9      35
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep NA | grep -v rs
SNP     A1      A2      Freq1.Hapmap    b       se      p       N       NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
2542421 22881781 147770052
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc 
11320  101880 1032095
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
54931  494371 2829544
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
30223  271999 1488151
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc 
2540915 27950055 178369746
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
11308  124388 1161471
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc 
54712  601822 3463146
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
30211  332311 1847730
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
2546063 28006683 177806542
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
11329  124619 1157287
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc 
56657  623217 3581917
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
31299  344279 1909638
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
2560771 28168471 179038050
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
11561  127171 1180926
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
75202  827212 4807405
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
30242  332652 1852424
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
2559728 23037544 149736440
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
11452  103068 1048156
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
73932  665380 3853906
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
30130  271162 1487498
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | wc
2565397 23088565 149759722
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep , | wc
11461  103149 1045870
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | wc
76732  690580 3999491
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2014_5]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | tail -n +20 | grep NA | grep -v rs | wc
31579  284203 1559698
~~~

zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz | perl -lane 'splice(@F, 1, 2); print join("\t", @F);' | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.ColRightFormat.txt.gz
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz | perl -lane 'splice(@F, 1, 2); print join("\t", @F);' | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.ColRightFormat.txt.gz
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz | perl -lane 'splice(@F, 1, 2); print join("\t", @F);' | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.ColRightFormat.txt.gz

zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep -E 'rs6694089|rs4652773|rs2298265|rs1368380' | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz

~~~
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
rs2298265       T       C       1_149525666     0.117   1       0.117   -0.030  0.0046  8.3e-11 0.117   -0.030  0.0046  8.3e-11
rs6694089       A       G       1_170350503     0.217   1       0.217   0.039   0.0032  4.2e-33 0.217   0.039   0.0032  4.2e-33
rs4652773       A       G       1_181321449     0.474   1       0.474   0.024   0.0029  1.1e-15 0.474   0.024   0.0029  1.1e-15
rs1368380       T       C       5_171218236     0.425   1       0.425   0.032   0.0031  6.1e-25 0.425   0.032   0.0031  6.1e-25
~~~

zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | head -n 1000 | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.top1000.gz
zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | head -n 1000 | gzip > /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.top1000.gz

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.top1000.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.top1000.gz

~~~
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.top1000.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.top1000.gz | head -n 25
/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.top1000.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.top1000.gz
rs10002873      T       C       4_60322082      0.075   0       NA      NA      NA      NA      0.075   -0.0068 0.0069  0.3244
rs10883118      A       G       10_100283886    0.289   0       0.289   -0.0028 0.0032  0.39    NA      NA      NA      NA
rs10002758      G       T       4_41478453      0.2417  0       NA      NA      NA      NA      0.2417  -0.0079 0.0046  0.08591
rs10002976      A       G       4_180507195     0.1333  0       NA      NA      NA      NA      0.1333  5e-04   0.0051  0.9219
rs12268065      A       G       10_10035354     0.02    0       0.02    -0.00030        0.0059  0.96    NA      NA      NA      NA
rs9988679       A       G       10_100604778    0.397   0       0.397   0.0063  0.0031  0.044   NA      NA      NA      NA
rs7924209       T       C       10_100145916    0.288   0       0.288   0.0060  0.0032  0.059   NA      NA      NA      NA
rs10002199      T       G       4_20417561      0.35    0       NA      NA      NA      NA      0.35    0.008   0.0037  0.03061
rs11599208      T       C       10_100182361    0.195   0       0.195   -0.00050        0.0040  0.91    NA      NA      NA      NA
rs12569898      A       G       10_10048290     0.5     0       0.5     0.0023  0.0030  0.44    NA      NA      NA      NA
rs10002195      G       A       4_125014932     0.2417  0       NA      NA      NA      NA      0.2417  0.0024  0.0036  0.5054
rs10000017      C       T       4_84997148      0.2333  0       NA      NA      NA      NA      0.2333  -0.0034 0.0046  0.4598
rs10000545      G       C       4_110577326     0.1583  0       NA      NA      NA      NA      0.1583  -0.0023 0.0042  0.5871
rs10001548      C       T       4_166098830     0.3417  0       NA      NA      NA      NA      0.3417  -4e-04  0.0041  0.9223
rs17110576      A       G       10_100399966    0.092   0       0.092   0.0031  0.0048  0.52    NA      NA      NA      NA
rs9664571       T       C       10_10063767     0.033   0       0.033   0.010   0.0097  0.28    NA      NA      NA      NA
rs10002751      G       C       4_10652353      0.275   0       NA      NA      NA      NA      0.275   -8e-04  0.0043  0.8524
rs10002370      G       A       4_46893928      0.1167  0       NA      NA      NA      NA      0.1167  0.006   0.0053  0.2576
rs7097807       A       G       10_100081571    0.207   0       0.207   0.0050  0.0037  0.17    NA      NA      NA      NA
rs10001545      C       A       4_88395053      0.2583  0       NA      NA      NA      NA      0.2583  -0.0083 0.0044  0.05925
rs12766763      A       T       10_10068697     0.025   0       0.025   -0.0057 0.0081  0.48    NA      NA      NA      NA
rs11815392      A       C       10_10061543     0.033   0       0.033   0.0094  0.0098  0.34    NA      NA      NA      NA
rs7907555       T       C       10_100145799    0.3     0       0.3     -0.0061 0.0032  0.053   NA      NA      NA      NA
rs11815889      A       C       10_100313990    0.25    0       0.25    -0.0018 0.0032  0.57    NA      NA      NA      NA
rs7902229       A       G       10_10004243     0.458   0       0.458   -0.0028 0.0037  0.44    NA      NA      NA      NA
~~~

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.HeightBMIWHRadjBMI.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz

zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 >= 2) { print $0 } } ' | vi -

#Memory issue with so many GWAS hit SNPs (50g interactive session not enough -- currentl trying with a 70g session); so going to try and create a 'vs2' file of /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py that approaches this differently
#Actually it looks like using ql 70g is enough to get past this issue, so just kept using /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py
#cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs2.py
#mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs2.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs2.py.20150527.DidntEndUpUsing

zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c 
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'

~~~
[  mturchin20@bigmem02  ~]$Error2a -- line2[0] present more than once in a particular file (1: rs13598,C,A,0.1583,-0.003,0.0096,0.7547,69744.9,15_46489749)
Error2a -- line2[0] present more than once in a particular file (1: rs2160421,T,A,0.5583,-0.001,0.0037,0.787,233968,15_34746878)
Error2a -- line2[0] present more than once in a particular file (1: rs372268,G,A,0.975,-0.0113,0.0214,0.5975,150660,4_139504797)
Error2a -- line2[0] present more than once in a particular file (1: rs3752976,G,A,0.05172,0.0071,0.0089,0.425,224029,1_199381979)
Error2a -- line2[0] present more than once in a particular file (1: rs3780697,T,C,0.725,0.0035,0.0048,0.4659,233977,9_131618930)
Error2a -- line2[0] present more than once in a particular file (1: rs3819262,T,G,0.7167,0.0051,0.004,0.2023,233953,21_37382759)
Error2a -- line2[0] present more than once in a particular file (1: rs573071,A,G,0.6833,-0.0042,0.004,0.2937,232167,9_246993)
Error2a -- line2[0] present more than once in a particular file (1: rs5820007,C,T,0.1271,-0.0029,0.0061,0.6345,216740,17_28342660)
Error2a -- line2[0] present more than once in a particular file (1: rs5820007,T,C,0.875,0.0029,0.0061,0.6345,216740,17_28342660)
Error2a -- line2[0] present more than once in a particular file (1: rs627141,T,C,0.2833,0,0.0051,1,229425,18_8332476)
Error2a -- line2[0] present more than once in a particular file (1: rs7591686,A,G,0.7333,0.001,0.0043,0.8161,233849,2_45609171)
Error2a -- line2[0] present more than once in a particular file (1: rs969206,G,A,0.3917,-0.005,0.0037,0.1766,234003,9_74502957)
[  mturchin20@spudhead  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 >= 2) { print $0 } } '
[  mturchin20@spudhead  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 >= 2) { print $0 } } '
2 rs11576885
2 rs11771665
2 rs13598
2 rs2160421
2 rs372268
2 rs3752976
2 rs3780697
2 rs3819262
2 rs573071
3 rs5820007
2 rs627141
2 rs7591686
2 rs969206
[  mturchin20@spudhead  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 >= 2) { print $0 } } '
[  mturchin20@spudhead  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz | grep -w -E 'rs11576885|rs11771665|rs13598|rs2160421|rs372268|rs3752976|rs3780697|rs3819262|rs573071|rs5820007|rs627141|rs7591686|rs969206'
rs11576885      C       G       0.4083  0.0014  0.0042  0.7389  219514  1_204737275
rs11576885      G       C       0.4083  -0.0014 0.0042  0.7389  219514  1_204737275
rs11771665      G       A       0.15    -0.0043 0.0056  0.4426  222149  7_86704151
rs11771665      G       A       0.1525  -0.0043 0.0056  0.4426  222149  7_86704151
rs13598 A       C       0.8417  0.003   0.0096  0.7547  69744.9 15_46489749
rs13598 C       A       0.1583  -0.003  0.0096  0.7547  69744.9 15_46489749
rs2160421       A       T       0.5583  0.001   0.0037  0.787   233968  15_34746878
rs2160421       T       A       0.5583  -0.001  0.0037  0.787   233968  15_34746878
rs372268        A       G       0.975   0.0113  0.0214  0.5975  150660  4_139504797
rs372268        G       A       0.975   -0.0113 0.0214  0.5975  150660  4_139504797
rs3752976       A       G       0.05    -0.0071 0.0089  0.425   224029  1_199381979
rs3752976       G       A       0.05172 0.0071  0.0089  0.425   224029  1_199381979
rs3780697       C       T       0.1833  -0.0035 0.0048  0.4659  233977  9_131618930
rs3780697       T       C       0.725   0.0035  0.0048  0.4659  233977  9_131618930
rs3819262       T       G       0.6917  0.0051  0.004   0.2023  233953  21_37382759
rs3819262       T       G       0.7167  0.0051  0.004   0.2023  233953  21_37382759
rs573071        A       G       0.675   -0.0042 0.004   0.2937  232167  9_246993
rs573071        A       G       0.6833  -0.0042 0.004   0.2937  232167  9_246993
rs5820007       C       T       0.125   -0.0029 0.0061  0.6345  216740  17_28342660
rs5820007       C       T       0.1271  -0.0029 0.0061  0.6345  216740  17_28342660
rs5820007       T       C       0.875   0.0029  0.0061  0.6345  216740  17_28342660
rs627141        T       C       0.1833  0       0.0051  1       229425  18_8332476
rs627141        T       C       0.2833  0       0.0051  1       229425  18_8332476
rs7591686       A       G       0.7167  0.001   0.0043  0.8161  233849  2_45609171
rs7591686       A       G       0.7333  0.001   0.0043  0.8161  233849  2_45609171
rs969206        A       G       0.3917  0.005   0.0037  0.1766  234003  9_74502957
rs969206        G       A       0.3917  -0.005  0.0037  0.1766  234003  9_74502957
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
2589848 46617264 297715972
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
2589848 17
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c
2017735 0
    861 1
 571252 2
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
3448    NA
3       None
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "None") { print $0 } } '
rs34102154      A       G       None    NA      0       NA      NA      NA      NA      NA      NA      NA      NA      NA      -0.0035 0.0084  51224
rs1053838       A       C       None    0.009   0       0.009   -0.013  0.039   101301  0.008621        0.0295  0.0502  97215.8 0.008621        0.044   0.059   57703
rs1632856       T       C       None    0.283   0       0.283   -0.0073 0.0032  251823  0.2833  0.0037  0.0041  232620  0.2833  0.0031  0.0047  141363
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
38285   NA
~~~

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.txt --file2 <(zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_Orig3_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'my $chr; my $pos; if (($F[3] eq "NA") || ($F[3] eq "None")) { $chr = "NA"; $pos = "NA"; } else { $chr = ((split(/_/, $F[3]))[0]); $pos = ((split(/_/, $F[3]))[1]); } splice(@F, 1, 0, ($chr, $pos)); splice(@F, 8, 1); splice(@F, 11, 1); splice(@F, 14, 1); print join("\t", @F);') | \
gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz
#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz
#snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI Z.height Z.BMI Z.WHRadjBMI mvstat mvp unip

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.txt | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc

~~~
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/dtlesssignif.annot.txt | wc
8066  209716 1794390
[  mturchin20@spudhead  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/multivariate/globallipids/dtlesssignif.annot.txt | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
3724   96824  830117
[  mturchin20@spudhead  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz | wc
12463  348964 3301334
[  mturchin20@spudhead  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2013/Vs1/dtlesssignif.annot.vs2.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
8859  248052 2346362
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz | wc
8842  203366 1641132
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/GIANT2010.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
1097   25231  202912
[  mturchin20@bigmem02  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.txt | wc
42023   42023 4708300
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz | wc
42023  966529 8803594
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
5677  130571 1185565
~~~       

#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.Orig3.vs1.R into an inteactve session in #PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/ 

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/RData.GIANT2014_5.Orig3.GLC.20150528

gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/RData.GIANT2014_5.Orig3.GLC.20150528



python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/GIANT2014_5.AllGWASHits.AllPheno.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/SNP_gwas_mc_merge_nogc.tbl.uniq.wUCSCGenomeBrowser_dbSNP130.vs1.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHRadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WHR_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.ColRightFormat.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIPadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.ColRightFormat.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_HIP_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.ColRightFormat.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WC_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz 

zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "None") { print $0 } } '
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c 

~~~
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
2591354 98471452 649782849
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
2591354 37
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
3454    NA
4       None
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "None") { print $0 } } '
rs34102154      A       G       None    NA      0       NA      NA      NA      NA      NA      NA      NA      NA      NA      -0.0035 0.0084  51224   NA      -0.0036 0.0079  51473   NA      -0.003  0.0088  51197   NA      NA      NA      NA      NA      -0.0004 0.0081  61707   NA      NA      NA      NA
rs1053838       A       C       None    0.009   0       0.009   -0.013  0.039   101301  0.008621        0.0295  0.0502  97215.8 0.008621        0.044   0.059   57703   0.008621        0.035   0.063   50566   0.008621        -0.065  0.062   58456   0.008621        -0.036  0.065   52190   0.0086  0.01    0.058   63735       0.008621        0.042   0.064   54959
rs1632856       T       C       None    0.283   0       0.283   -0.0073 0.0032  251823  0.2833  0.0037  0.0041  232620  0.2833  0.0031  0.0047  141363  0.2833  0.0054  0.0046  143200  0.2833  0.0036  0.0049  142955  0.2833  0.0051  0.0048  144058  0.2833  0.0098  0.0046  152083  0.2833  0.0078  0.0047  152552
rs6603251       T       C       None    NA      0       NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      0.0083  0.0068  56955   NA      NA      NA      NA
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
39789   NA
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c
2003823 0
    880 1
 586651 2
~~~

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.txt --file2 <(zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'my $chr; my $pos; if (($F[3] eq "NA") || ($F[3] eq "None")) { $chr = "NA"; $pos = "NA"; } else { $chr = ((split(/_/, $F[3]))[0]); $pos = ((split(/_/, $F[3]))[1]); } splice(@F, 1, 0, ($chr, $pos)); splice(@F, 8, 1); splice(@F, 11, 1); splice(@F, 14, 1); splice(@F, 17, 1); splice(@F, 20, 1); splice(@F, 23, 1); splice(@F, 27, 1); splice(@F, 30, 1); print join("\t", @F);') | \
gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz
#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz
#snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI  beta_WHR se_WHR n_WHR beta_HIPadjBMI se_HIPadjBMI n_HIPadjBMI beta_HIP se_HIP n_HIP beta_WCadjBMI se_WCadjBMI n_WCadjBMI beta_WC se_WC n_WC Z.height Z.BMI Z.WHRadjBMI Z.WHR Z.HIPadjBMI Z.HIP Z.WCadjBMI Z.WC mvstat mvp unip

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.txt | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz | perl -F, -lane 'print $#F;' | sort | uniq -c

~~~
[  mturchin20@spudhead  ~/Data/GIANT/2014_5]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'my $chr; my $pos; if (($F[3] eq "NA") || ($F[3] eq "None")) { $chr = "NA"; $pos = "NA"; } else { $chr = ((split(/_/, $F[3]))[0]); $pos = ((split(/_/, $F[3]))[1]); } splice(@F, 1, 0, ($chr, $pos)); splice(@F, 8, 1); splice(@F, 11, 1); splice(@F, 14, 1); splice(@F, 17, 1); splice(@F, 20, 1); splice(@F, 23, 1); splice(@F, 27, 1); splice(@F, 30, 1); print join("\t", @F);' | head -n 10
rs999721        4       85416326        A       G       4_85416326      0.466   0       0.0097  0.0030  241371  -0.0039 0.0038  222236  0.0034  0.0044  130986  0.0006  0.0043  132832  0.0046  0.0046  132042  -0.0001 0.0046  133691  0.4667  0.0044  140848  0.4667  0.0044  142183
rs942176        9       8420741 A       G       9_8420741       0.433   0       0.0040  0.0030  253133  0.0029  0.0038  233938  6e-04   0.0043  142727  0.0031  0.0042  144568  0.0033  0.0045  143781  0.0064  0.0045  145424  0.4333  0.0043  151195  0.4333  0.0044  153921
rs9821657       3       186113796       A       G       3_186113796     0.175   0       0.0061  0.0035  252892  4e-04   0.0046  233804  0       0.005   142556  -0.0034 0.0051  144392  0.0031  0.0053  143637  0.0017  0.0053  145195  0.175   0.005   153275  0.175   0.0052  153744
rs1919329       4       165827169       T       G       4_165827169     0.483   0       0.0078  0.0030  246965  -2e-04  0.0038  230282  0.0022  0.0043  140365  0.002   0.0042  142384  0.005   0.0045  141418  0.0024  0.0045  143242  0.4833  0.0043  151084  0.4833  0.0044  151736
rs3957240       6       37262200        T       C       6_37262200      0.05    0       -0.024  0.011   215562  -0.0106 0.014   202410  -0.011  0.016   115779  -0.027  0.016   117279  -0.0046 0.017   116818  -0.021  0.017   118121  0.05    0.016   126473  0.05    0.016   126597
rs1919324       4       165805562       T       C       4_165805562     0.142   0       -0.0016 0.0039  253148  0.0028  0.0051  233978  0.0085  0.0057  142610  0.0089  0.0056  144588  -0.0083 0.0061  143664  -0.0025 0.006   145448  0.1417  0.0057  153330  0.1417  0.0058  153940
rs8030214       15      49019617        T       C       15_49019617     0.336   2       -0.015  0.0030  253230  2e-04   0.0038  234018  9e-04   0.0043  142752  0.0001  0.0043  144595  -0.0098 0.0045  143805  -0.0078 0.0045  145454  0.3333  0.0043  153471  0.3333  0.0044  153947
rs17597444      13      42193631        T       G       13_42193631     0.237   0       0.0067  0.0033  253153  0.0041  0.0042  234007  0.0097  0.0047  142518  0.0075  0.0046  144592  -0.012  0.0049  143571  -0.0039 0.0049  145451  0.25    0.0047  153237  0.25    0.0047  153944
rs2664913       10      8169772 A       G       10_8169772      0.425   0       -0.0016 0.0031  240207  -0.004  0.004   221984  0.0044  0.0045  130431  0.0007  0.0044  132745  -0.0015 0.0047  131498  -0.0025 0.0047  133604  0.425   0.0045  141150  0.425   0.0045  142096
rs16876997      5       6246251 A       G       5_6246251       0.017   0       0.0065  0.011   238862  -0.0022 0.0138  220635  -0.016  0.015   130612  -0.024  0.015   132437  -0.008  0.016   131665  -0.01   0.016   133289  0.0167  0.015   141330  0.0167  0.016   141788
[  mturchin20@bigmem02  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.txt | wc
42023   42023 4708300
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz | wc
42023  966529 8803594
[  mturchin20@bigmem02  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.Orig3.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
5677  130571 1185565
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.txt | wc
56323   56323 11147251
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz | wc
56323 2421889 22448471
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
41783 1796669 16626655
[  mturchin20@bigmem02  ~/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtlesssignif.vs1.txt | perl -F, -lane 'print $#F;' | sort | uniq -c
56323 11
~~~

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.Orig3.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.vs1.R

cd /mnt/lustre/home/mturchin20/Software
wget http://cran.r-project.org/src/contrib/bigmemory_4.4.6.tar.gz 
wget http://cran.r-project.org/src/contrib/bigmemory.sri_0.1.3.tar.gz
wget http://cran.r-project.org/src/contrib/BH_1.58.0-1.tar.gz
wget http://cran.r-project.org/src/contrib/pryr_0.1.1.tar.gz
wget http://cran.r-project.org/src/base/R-3/R-3.1.2.tar.gz
wget http://cran.r-project.org/src/contrib/stringi_0.4-1.tar.gz
wget http://cran.r-project.org/src/contrib/magrittr_1.5.tar.gz
wget http://cran.r-project.org/src/contrib/stringr_1.0.0.tar.gz
wget http://cran.r-project.org/src/contrib/Rcpp_0.11.6.tar.gz
wget http://cran.r-project.org/src/contrib/plyr_1.8.2.tar.gz
wget http://cran.r-project.org/src/contrib/reshape2_1.4.1.tar.gz
wget http://cran.r-project.org/src/contrib/chron_2.3-45.tar.gz
wget http://cran.r-project.org/src/contrib/data.table_1.9.4.tar.gz
wget http://cran.r-project.org/src/contrib/assertthat_0.1.tar.gz
wget http://cran.r-project.org/src/contrib/R6_2.0.1.tar.gz  
wget http://cran.r-project.org/src/contrib/lazyeval_0.1.10.tar.gz
wget http://cran.r-project.org/src/contrib/DBI_0.3.1.tar.gz
wget http://cran.r-project.org/src/contrib/dplyr_0.4.1.tar.gz
tar -xvzf /mnt/lustre/home/mturchin20/Software/R-3.1.2.tar.gz
./configure --prefix=/mnt/lustre/home/mturchin20/lib --with-x

##Having memory issues running as large a matrix as I needed to in /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.vs1.R 
##Particularly going from:
~~~
#lbf=compute.allBFs.fromZscores(Z,VYY,gl$nmin,gl$maf,sigmaa)
#lbf.bigmat=do.call(rbind,lbf$lbf) #lbf$lbf is a list of matrices, with one element (a matrix) for each sigma; this stacks these matrices together into a big matrix with nsnp columns, and nsigma*nmodels rows
~~~
#So going to try and use the 'more significant' list of SNPs for this phenotype to cut down on computation

##20150601 -- Realized I can break up the matrix into groups of 10000 columns to overcome the memory issue. So have 7 lbf.bigmat#s and then can rbind all 7 together with 1:10000, 10001:20000, etc.. etc.. 6 times (up until 56301)
##Also I do the GWAS prior hit steps before the lbf.bigmat creation process too now
##As far as I can tell the compute.allBFs.fromZscores(...) steps and all related processes just need every model present per snp but doesn't need any informatino across SNPs to work properl; so breaking things up on a per SNP basis should be fine as long as I'm keeping all the model information for a given SNP (and maintaing the SNP order so that the results are consistent)

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.vs2.R

#python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtsignif.vs1.txt --file2 <(zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5_AllPheno_publicrelease_HapMapCeuFreq.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'my $chr; my $pos; if (($F[3] eq "NA") || ($F[3] eq "None")) { $chr = "NA"; $pos = "NA"; } else { $chr = ((split(/_/, $F[3]))[0]); $pos = ((split(/_/, $F[3]))[1]); } splice(@F, 1, 0, ($chr, $pos)); splice(@F, 8, 1); splice(@F, 11, 1); splice(@F, 14, 1); splice(@F, 17, 1); splice(@F, 20, 1); splice(@F, 23, 1); splice(@F, 27, 1); splice(@F, 30, 1); print join("\t", @F);') | \
#gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtsignif.vs1.annot.vs1.txt.gz
##Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtsignif.vs1.annot.vs1.txt.gz
##snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI  beta_WHR se_WHR n_WHR beta_HIPadjBMI se_HIPadjBMI n_HIPadjBMI beta_HIP se_HIP n_HIP beta_WCadjBMI se_WCadjBMI n_WCadjBMI beta_WC se_WC n_WC Z.height Z.BMI Z.WHRadjBMI Z.WHR Z.HIPadjBMI Z.HIP Z.WCadjBMI Z.WC mvstat mvp unip

#~~~
#[  mturchin20@spudhead  ~/Software]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GIANT2014_5.AllPheno.dtsignif.vs1.txt | wc
#  38978   38978 7708677
#~~~
#
#cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.dtsignif.vs1.R 

#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.vs2.R into an inteactve session in #PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/

##20150601 -- .RData file became >10gb so decided to just remove it
#mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/RData.GIANT2014_5.AllPheno.GLC.MTedits.vs2.20150601
#gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/RData.GIANT2014_5.AllPheno.GLC.MTedits.vs2.20150601

#CHECK_0: Print out allele information for all versions of 'newtophits' in previous runs (GlobalLipids, GIANT2010, 2014_5_Orig3)
#CHECK_0: Think about ways to explore the diffuse qualities of posterior 'bestclass' classifications, e.g. how many of the top hits don't have a single model with evidence above .75 or .5


#2013
#Conducting multivariate analyses

#CHECK_0: Assuming for the moment that these analyses/studies were mapped to the same build that 2010 & 2014_5 were (should really be the case)
#But may want to double, double check this at a later point when there's a bit more time

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 
python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/Python.PerlVLookup.FindrsIDsindbSNPvcf.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/UCSC/GenomeBrowser/dbSNP130.txt.gz --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.txt | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz 

zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_*_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c

zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc 
zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc

#~~~
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_*_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
46583187 8
[  mturchin20@spudling74  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | head -n 10
rs10884221      t       c       0.8083  -0.0098 0.0072  0.17    59229.5 10_107486155,6_cox_hap1_48402
rs11193825      c       g       0.15    -0.0042 0.0089  0.64    60555.2 10_109725754,6_cox_hap1_4021800
rs11527583      a       c       0.4417  -0.0022 0.0061  0.72    60506.9 10_20698187,5_h2_hap1_1641990
rs2503868       c       g       .       0.0033  0.014   0.82    29451.1 10_43117647,6_cox_hap1_4404084,6_qbl_hap2_4209369
rs11818958      a       c       0.4583  -0.0022 0.0061  0.72    60556.6 10_43808197,6_cox_hap1_48162
rs10824920      t       c       0.5583  0.0040  0.0063  0.52    60505.1 10_54585619,6_cox_hap1_4449605,6_qbl_hap2_4255532
rs1046747       a       g       0.375   -0.0045 0.0063  0.47    60500.6 10_71580635,6_cox_hap1_498071,6_qbl_hap2_300095
rs15801 t       c       0.4     -0.0053 0.0061  0.39    60538.7 10_71580668,6_cox_hap1_498102,6_qbl_hap2_300126
rs2394643       t       c       0.4     -0.0054 0.0061  0.38    60536.9 10_71581019,6_cox_hap1_498457,6_qbl_hap2_300481
rs17343324      a       g       .       -0.0010 0.015   0.95    12317.3 10_88799902,10_random_105971
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13827  124443 1245567
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
0       0       0
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2748520 24736680 160765048
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13760  123840 1238783
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
0       0       0
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2736850 24631650 159851677
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13678  123102 1230687
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
0       0       0
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2736850 24631650 159924050
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13678  123102 1229700
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | grep -v rs | wc
0       0       0
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA
rs4510589       a       g       .       0.0079  0.0096  0.41    27670   NA
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep -v rs | wc
1       9      59
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep -v rs
MarkerName      A1      A2      Freq.Hapmap.Ceu BETA    SE.2gc  P.2gc   N       ChrBP
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2742719 24684471 159937922
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13723  123507 1232512
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      45
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSC
GenomeBrowser_dbSNP130.vs1.txt.gz | wc
2737330 24635970 159690178
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13621  122589 1224175
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      43
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2743314 24689826 159474582
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13724  123516 1230202
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2738277 24644493 159488583
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13627  122643 1223158
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2725770 24531930 158583551
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13532  121788 1215259
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      43
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2725770 24531930 158193487
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13532  121788 1209939
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      43
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBr
owser_dbSNP130.vs1.txt.gz | wc
2725730 24531570 158233458
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13530  121770 1213552
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      42
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2725730 24531570 158587868
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13530  121770 1214335
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      43
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2744340 24699060 161028253
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13724  123516 1238481
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      45
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2738312 24644808 160585060
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13626  122634 1227995
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2745518 24709662 159596765
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13725  123525 1229833
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      43
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2740359 24663231 159625213
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13631  122679 1223141
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      44
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2760764 24846876 161121016
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13829  124461 1244620
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      43
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | wc
2746981 24722829 160322057
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep , | wc
13755  123795 1235871
[  mturchin20@spudling52  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep NA | wc
1       9      43
#~~~


mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1
cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/DataOrganization.GIANT2014_5.PythonVLookUp.CreateFullAnnotFile.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/DataOrganization.GIANT2013.PythonVLookUp.CreateFullAnnotFile.vs1.py

zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | grep -E 'rs17524355|rs12219605|rs4747841|rs737656' | gzip > /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz

#~~~
[  mturchin20@spudling49  ~]$zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
rs4747841       a       g       0.55    0.0025  0.0061  0.68    60558.2 10_10000134
rs737656        a       g       0.3667  -0.0073 0.0064  0.25    60529.2 10_100002728
rs17524355      t       c       .       -0.046  0.046   0.32    1700    10_100003427
rs12219605      t       g       0.45    -0.0025 0.0061  0.68    60558.2 10_10000458
rs7376568       t       g       0.2833  0.0069  0.0078  0.38    59638.3 4_189370200
[  mturchin20@spudling49  ~]$python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/DataOrganization.GIANT2013.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 <(cat /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Male.MarkerChrBP.txt <(zcat /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz | head -n 1 | perl -lane 'print $F[0], "\t", ((split(/_/, $F[8]))[0]), "\t", ((split(/_/, $F[8]))[1]);')) --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.4SNPs.txt.gz
rs4747841       a       g       10_10000134     0.45    1       0.45    0.0025  0.0061  60558.2 0.45    0.0025  0.0061  60558.2
rs12219605      t       g       10_10000458     0.45    2       0.45    -0.0025 0.0061  60558.2 0.45    -0.0025 0.0061  60558.2
rs737656        a       g       10_100002728    0.3667  0       0.3667  -0.0073 0.0064  60529.2 0.3667  -0.0073 0.0064  60529.2
rs7376568       t       g       4_189370200     0.2833  0       0.2833  0.0069  0.0078  59638.3 0.2833  0.0069  0.0078  59638.3
rs17524355      t       c       10_100003427    NA      0       NA      -0.046  0.046   1700    NA      -0.046  0.046   1700
#~~~

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/DataOrganization.GIANT2013.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Male.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz

zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "None") { print $0 } } '
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c

~~~
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
2767844 116249448 770991833
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
2767844 41
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
1       NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "None") { print $0 } } '
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
507290  NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c
2754614 0
     15 1
  13215 2
~~~       

#CHECK_0: Take a look at these rsID#s/SNPs specifically for both male and female file results

#~~~
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs424881,a,g,0.225,0.00080,0.0075,0.91,60550.2,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs2160421,a,t,0.4417,-0.0052,0.0061,0.40,60541.7,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs13598,a,c,.,0.0072,0.0088,0.42,60557.3,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs5820007,t,c,0.875,0.0098,0.0095,0.30,59008.1,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs627146,a,g,0.2833,-0.012,0.0064,0.064,59129.2,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs3752976,a,g,.,-0.0051,0.014,0.71,60550.8,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs11576885,c,g,.,-0.0028,0.0068,0.68,59668.4,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs3819262,t,g,0.6917,-0.0081,0.0064,0.20,60557.4,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs12996248,t,g,0.2833,0.0050,0.0074,0.50,58638.6,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs7591686,a,g,0.7333,-0.0020,0.0072,0.78,60529.7,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs2936687,a,g,0.8167,-0.015,0.0082,0.068,60556.5,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs3780697,t,c,0.725,0.0042,0.0079,0.60,60527.2,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: rs573071,a,g,0.675,-0.0054,0.0065,0.40,60495.8,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs424881,a,g,0.225,0.0093,0.0078,0.23,58618.1,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs2160421,a,t,0.4417,0.0039,0.0064,0.54,58617.6,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs13598,a,c,0.8417,-0.017,0.0091,0.059,58619.2,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs5820007,t,c,.,-0.011,0.0099,0.27,57091,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs627146,a,g,0.2833,0.00080,0.0067,0.90,57248.1,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs3752976,a,g,.,-0.011,0.014,0.44,58619.7,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs11576885,c,g,.,-0.00050,0.0070,0.94,58611.3,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs3819262,t,g,0.6917,-0.00020,0.0067,0.98,58618.3,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs12996248,t,g,0.2833,-0.0055,0.0078,0.48,56724.5,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs7591686,a,g,0.7167,0.012,0.0073,0.094,58613.6,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs2936687,a,g,0.8167,0.024,0.0084,0.0041,58618.4,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs3780697,t,c,0.725,-0.0013,0.0082,0.87,58607.2,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: rs573071,a,g,0.6833,-0.0019,0.0068,0.78,58614.7,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs424881,a,g,0.225,0.00E+00,0.0083,1.0,34601.4,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs2160421,a,t,0.4417,0.0017,0.0067,0.80,34600.9,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs13598,a,c,0.8417,0.0071,0.0098,0.47,34601.8,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs5820007,t,c,0.875,-0.0021,0.010,0.84,33318.2,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs627146,a,g,0.2833,0.0010,0.0070,0.89,34595.9,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs3752976,a,g,.,0.0055,0.014,0.70,34602.2,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs11576885,c,g,.,-0.0059,0.0075,0.43,34599.3,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs3819262,t,g,0.7167,0.012,0.0071,0.080,34600.8,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs12996248,t,g,0.2833,-0.0037,0.0083,0.66,32709.5,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs7591686,a,g,0.7333,0.0029,0.0077,0.71,34596.2,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs2936687,a,g,0.8167,0.0078,0.0089,0.38,34601.1,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs3780698,t,g,0.8167,0.0065,0.0078,0.40,34601.2,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs573071,a,g,0.675,-0.00090,0.0071,0.90,34599.2,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs424881,a,g,0.775,-0.0014,0.0094,0.88,34581.2,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs2160421,a,t,0.4417,0.0049,0.0076,0.52,34580.7,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs13598,a,c,.,-0.0018,0.011,0.87,34581.6,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs5820007,t,c,.,0.0036,0.012,0.76,33297,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs627146,a,g,0.2833,0.0055,0.0079,0.49,34575.8,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs3752976,a,g,.,-0.0097,0.016,0.55,34582,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs11576885,c,g,.,-0.0068,0.0085,0.43,34579.1,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs3819262,t,g,0.7167,0.013,0.0080,0.097,34580.6,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs12996248,t,g,0.2833,-0.0043,0.0095,0.65,32689.3,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs7591686,a,g,0.7167,0.0057,0.0087,0.51,34576,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs2936687,a,g,0.8167,0.025,0.010,0.013,34580.9,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs3780698,t,g,0.725,0.0057,0.0088,0.52,34581,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs573071,a,g,0.675,-0.0024,0.0081,0.77,34579,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs424881,a,g,0.775,0.0037,0.0099,0.71,32930.1,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs2160421,a,t,.,0.00070,0.0079,0.93,32930.2,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs13598,a,c,.,0.0032,0.011,0.77,32930.4,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs5820007,t,c,.,0.0028,0.013,0.83,31647.1,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs627146,a,g,0.2833,-0.0025,0.0082,0.76,32925.3,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs3752976,a,g,.,0.0019,0.018,0.92,32931.5,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs11576885,c,g,.,-0.011,0.0087,0.20,32927.8,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs3819262,t,g,0.6917,-0.0012,0.0083,0.89,32929.4,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs12996248,t,g,0.2833,0.0058,0.0097,0.55,32924,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs7591686,a,g,0.7167,-0.0068,0.0093,0.46,32925.5,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs2936687,a,g,0.8167,-0.016,0.011,0.14,32930.1,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs3780698,t,g,0.725,-0.012,0.0092,0.18,32930.5,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs573071,a,g,0.6833,0.010,0.0084,0.21,32928.2,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs424881,a,g,0.775,0.0046,0.011,0.67,30603,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs2160421,a,t,.,0.0064,0.0083,0.44,32849,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs13598,a,c,0.8417,-0.0059,0.012,0.62,32849.2,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs5820007,t,c,.,0.018,0.013,0.18,31564.9,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs627146,a,g,.,0.0014,0.0087,0.87,32844.3,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs3752976,a,g,.,-0.029,0.019,0.13,32850.3,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs11576885,c,g,.,0.00070,0.0093,0.94,32846.7,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs3819262,t,g,0.6917,-0.0046,0.0087,0.60,32848.2,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs12996248,t,g,0.2833,0.0041,0.010,0.69,32842.8,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs7591686,a,g,0.7167,-0.0020,0.0096,0.84,32844.3,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs2936687,a,g,0.8167,0.016,0.011,0.15,32848.9,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs3780698,t,g,0.8167,0.00060,0.0096,0.95,32849.3,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs573071,a,g,0.6833,0.010,0.0089,0.24,32847,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs424881,a,g,0.225,0.0033,0.0057,0.57,38369.6,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs2160421,a,t,.,0.0024,0.0047,0.61,38369.1,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs13598,a,c,.,0.0024,0.0069,0.73,38370,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs5820007,t,c,.,0.0054,0.0075,0.47,36845.4,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs627146,a,g,0.2833,-0.0012,0.0050,0.81,37013.2,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs3752976,a,g,.,-0.0058,0.0099,0.56,38370.4,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs11576885,c,g,.,-0.0044,0.0053,0.41,38367.5,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs3819262,t,g,0.7167,0.00080,0.0049,0.87,38369,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs12996248,t,g,0.2833,0.0036,0.0057,0.53,36477.7,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs7591686,a,g,0.7167,-0.00040,0.0054,0.94,38364.4,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs2936687,a,g,0.8167,-0.0024,0.0063,0.70,38369.3,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs3780697,t,c,0.725,-0.0058,0.0061,0.34,38364.2,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs573071,a,g,0.6833,-0.0014,0.0050,0.78,38367.4,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs424881,a,g,0.775,0.0052,0.0092,0.57,38305,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs2160421,a,t,.,0.0098,0.0074,0.19,38304.5,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs13598,a,c,0.8417,-0.0076,0.011,0.48,38305.4,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs5820007,t,c,.,0.0068,0.012,0.56,36779.8,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs627146,a,g,0.2833,0.0025,0.0078,0.75,36948.8,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs3752976,a,g,.,-0.030,0.016,0.069,38305.8,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs11576885,c,g,.,-0.0015,0.0084,0.86,38302.9,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs3819262,t,g,0.6917,0.0029,0.0078,0.71,38304.4,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs12996248,t,g,0.2833,-0.00020,0.0092,0.98,36413.1,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs7591686,a,g,0.7167,0.0075,0.0086,0.38,38299.8,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs2936687,a,g,0.8167,0.026,0.0099,0.0093,38304.7,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs3780697,t,c,0.725,-0.0015,0.0096,0.88,38299.9,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs573071,a,g,0.6833,-0.0013,0.0079,0.87,38302.8,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs424881,a,g,0.775,0.0080,0.0078,0.31,58322.7,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs2160421,a,t,0.4417,0.00020,0.0064,0.98,58322.2,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs13598,a,c,0.8417,-0.015,0.0091,0.099,58323.8,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs5820007,t,c,0.875,-0.0052,0.010,0.60,56795.6,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs627146,a,g,.,-0.0033,0.0068,0.62,56953.3,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs3752976,a,g,.,-0.0072,0.015,0.62,58324.3,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs11576885,c,g,.,0.00080,0.0072,0.91,58315.9,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs3819262,t,g,0.7167,-0.0058,0.0068,0.39,58322.9,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs12996248,t,g,0.2833,0.00070,0.0078,0.93,56429.1,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs7591686,a,g,0.7167,0.0081,0.0074,0.27,58318.2,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs2936687,a,g,0.8167,0.016,0.0085,0.060,58323,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs3780697,t,c,0.725,0.0045,0.0083,0.59,58312.2,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs573071,a,g,0.6833,-0.0034,0.0069,0.62,58319.3,9_246993)
#~~~

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/DataOrganization.GIANT2013.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Female.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz

zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "None") { print $0 } } '
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c

~~~
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | wc
2752595 115608990 769011235
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'print $#F;' | sort | uniq -c
2752595 41
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $4 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
1       NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ if ($4 == "None") { print $0 } } '
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $5 }' | sort | uniq -c | perl -lane 'if ($F[1] !~ m/\d+/) { print join("\t", @F); }'
492041  NA
[  mturchin20@bigmem02  /mnt/gluster/data/external_public_supp/GIANT2013]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | awk '{ print $6 }' | sort | uniq -c
2717426 0
     38 1
  35131 2
~~~

#~~~
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs424881,a,g,0.775,0.0023,0.0069,0.74,73128.9,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs2160421,a,t,0.4417,-0.0091,0.0056,0.11,73120.5,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs13598,a,c,.,-0.0095,0.0082,0.25,73130,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs5820007,t,c,.,0.0031,0.0091,0.73,70813.9,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs627146,a,g,.,0.0079,0.0060,0.19,71183.1,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs3752976,a,g,.,-0.0098,0.013,0.45,73133.5,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs11576885,c,g,.,0.0016,0.0063,0.80,72203.5,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs3819262,t,g,0.7167,-0.0074,0.0061,0.22,73136.1,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs12996248,t,g,0.2833,-0.0025,0.0069,0.72,70681,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs7591686,a,g,0.7333,-0.0065,0.0065,0.32,73112.3,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs2936687,a,g,0.8167,-0.0017,0.0075,0.82,73134.7,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs3780697,t,c,0.8167,-0.00050,0.0074,0.95,73092.2,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 0 line2[0]: 0: rs573071,a,g,0.675,-0.00060,0.0061,0.92,73077.3,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs424881,a,g,0.225,-0.0025,0.0075,0.74,67956.5,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs2160421,a,t,.,0.0043,0.0060,0.48,67956.1,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs13598,a,c,.,-0.0021,0.0088,0.81,67957.6,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs5820007,t,c,.,0.0098,0.0098,0.32,65653.5,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs627146,a,g,.,-0.0080,0.0065,0.22,66069.8,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs3752976,a,g,.,0.0046,0.014,0.74,67958.1,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs11576885,c,g,.,-0.00040,0.0068,0.95,67952.1,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs3819262,t,g,0.6917,0.013,0.0065,0.047,67956.7,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs12996248,t,g,0.2833,0.0058,0.0075,0.44,65534.6,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs7591686,a,g,0.7333,-0.0028,0.0070,0.69,67952.9,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs2936687,a,g,0.8167,0.010,0.0081,0.21,67956.3,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs3780697,t,c,0.725,0.00060,0.0080,0.94,67938,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 1 line2[0]: 1: rs573071,a,g,0.6833,-0.0057,0.0066,0.39,67953.9,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs424881,a,g,0.225,0.0048,0.0081,0.55,42735.4,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs2160421,a,t,0.4417,-0.0059,0.0066,0.37,42735,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs13598,a,c,0.8417,-0.00090,0.0098,0.93,42735.8,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs5820007,t,c,0.875,0.0089,0.010,0.39,40783.3,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs627146,a,g,0.2833,0.0022,0.0069,0.75,42728.4,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs3752976,a,g,.,0.0087,0.015,0.55,42736.3,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs11576885,c,g,.,0.013,0.0074,0.072,42733.3,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs3819262,t,g,0.6917,-0.0033,0.0070,0.64,42734.9,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs12996248,t,g,0.2833,0.0032,0.0080,0.69,40315,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs7591686,a,g,0.7333,-0.00040,0.0076,0.96,42731.2,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs2936687,a,g,0.8167,0.0096,0.0088,0.28,42734.8,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs3780698,t,g,0.725,-0.0032,0.0077,0.68,42735.3,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 2 line2[0]: rs573071,a,g,0.6833,0.0050,0.0071,0.48,42733.2,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs424881,a,g,0.225,-0.00070,0.0087,0.94,42732.1,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs2160421,a,t,0.4417,-0.0067,0.0070,0.34,42731.7,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs13598,a,c,0.8417,-0.00070,0.010,0.95,42732.5,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs5820007,t,c,0.875,0.0073,0.011,0.51,40780,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs627146,a,g,0.2833,0.00060,0.0074,0.94,42725.2,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs3752976,a,g,.,0.011,0.016,0.46,42733,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs11576885,c,g,.,0.0098,0.0079,0.21,42730,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs3819262,t,g,0.7167,0.0024,0.0075,0.75,42731.6,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs12996248,t,g,0.2833,0.011,0.0086,0.21,40311.7,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs7591686,a,g,0.7167,-0.0025,0.0081,0.76,42727.9,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs2936687,a,g,0.8167,0.011,0.0094,0.26,42731.5,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs3780698,t,g,0.725,0.0011,0.0082,0.89,42732,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 3 line2[0]: rs573071,a,g,0.6833,0.0021,0.0075,0.78,42729.9,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs424881,a,g,0.225,0.0032,0.0092,0.73,40724.6,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs2160421,a,t,0.4417,0.00030,0.0074,0.97,40724.7,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs13598,a,c,.,-0.010,0.011,0.34,40725,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs5820007,t,c,0.875,0.0040,0.012,0.74,38772.7,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs627146,a,g,.,-0.00040,0.0077,0.96,40718.2,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs3752976,a,g,.,-0.0032,0.017,0.85,40387,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs11576885,c,g,.,-0.0020,0.0082,0.81,40722.3,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs3819262,t,g,0.7167,0.0036,0.0077,0.64,40724,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs12996248,t,g,0.2833,0.0023,0.0089,0.80,40720,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs7591686,a,g,0.7333,0.0074,0.0087,0.39,40720.9,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs2936687,a,g,0.8167,-0.0077,0.0099,0.44,40724.4,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs3780698,t,g,0.725,0.0093,0.0086,0.28,40725,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 4 line2[0]: rs573071,a,g,0.6833,0.0010,0.0079,0.90,40722.6,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs424881,a,g,0.775,-0.014,0.0097,0.16,40355.4,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs2160421,a,t,0.4417,-0.00080,0.0078,0.92,40355.4,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs13598,a,c,.,-0.0081,0.012,0.48,40355.8,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs5820007,t,c,0.875,-0.0023,0.012,0.85,38403.4,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs627146,a,g,0.2833,-0.0060,0.0084,0.47,40349.3,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs3752976,a,g,.,0.016,0.018,0.37,40356.8,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs11576885,c,g,.,-0.0091,0.0089,0.31,40353.1,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs3819262,t,g,0.6917,0.014,0.0083,0.098,40354.8,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs12996248,t,g,0.2833,0.019,0.0095,0.044,40350.8,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs7591686,a,g,0.7167,-0.0052,0.0092,0.57,40351.7,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs2936687,a,g,0.8167,-0.0047,0.011,0.66,40355.2,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs3780698,t,g,0.725,0.017,0.0091,0.060,40355.8,9_131618931)
Error2a -- line2[0] present more than once in a particular file (count: 5 line2[0]: rs573071,a,g,0.6833,-0.0025,0.0085,0.77,40353.4,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs424881,a,g,0.225,0.0014,0.0056,0.80,47471.7,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs2160421,a,t,.,-0.0030,0.0046,0.51,47471.3,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs13598,a,c,.,-0.0050,0.0068,0.46,47472.1,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs5820007,t,c,.,0.0020,0.0072,0.78,45180.6,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs627146,a,g,.,0.0023,0.0049,0.64,45600.6,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs3752976,a,g,.,0.0046,0.010,0.65,47472.6,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs11576885,c,g,.,0.0095,0.0051,0.063,47469.6,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs3819262,t,g,0.7167,-0.00010,0.0048,0.98,47471.2,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs12996248,t,g,0.2833,0.0036,0.0055,0.51,45050.3,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs7591686,a,g,0.7167,0.0029,0.0052,0.58,47467.5,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs2936687,a,g,0.8167,0.0076,0.0060,0.21,47471.1,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs3780697,t,c,0.8167,-0.0015,0.0059,0.80,47465.6,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 6 line2[0]: rs573071,a,g,0.675,0.0019,0.0049,0.70,47469.5,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs424881,a,g,0.775,-0.0080,0.0086,0.35,47320.6,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs2160421,a,t,.,-0.0047,0.0070,0.50,47320.2,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs13598,a,c,.,-0.0051,0.010,0.62,47321,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs5820007,t,c,.,0.0049,0.011,0.66,45029.5,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs627146,a,g,.,-0.0033,0.0074,0.66,45447.9,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs3752976,a,g,.,0.011,0.016,0.50,47321.5,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs11576885,c,g,.,0.0041,0.0079,0.60,47318.5,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs3819262,t,g,0.7167,0.012,0.0073,0.088,47320.1,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs12996248,t,g,0.2833,0.013,0.0085,0.14,44899.1,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs7591686,a,g,0.7333,0.00060,0.0081,0.94,47316.4,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs2936687,a,g,0.8167,0.0077,0.0093,0.41,47320,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs3780697,t,c,0.725,-0.0018,0.0091,0.84,47314.9,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 7 line2[0]: rs573071,a,g,0.6833,-0.0051,0.0074,0.49,47318.4,9_246993)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs424881,a,g,0.225,-0.0029,0.0075,0.70,67592.2,10_132097974)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs2160421,a,t,.,-0.0021,0.0061,0.73,67591.8,15_34746878)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs13598,a,c,.,-0.0082,0.0089,0.35,67593.3,15_46489749)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs5820007,t,c,0.875,0.0089,0.010,0.37,65288.2,17_28342660)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs627146,a,g,0.2833,-0.0068,0.0065,0.30,65706.2,18_8332477)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs3752976,a,g,.,0.0016,0.014,0.91,67593.8,1_199381979)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs11576885,c,g,.,0.0010,0.0069,0.88,67587.8,1_204737275)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs3819262,t,g,0.7167,0.0096,0.0065,0.14,67592.4,21_37382759)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs12996248,t,g,0.2833,0.0049,0.0074,0.51,65170.3,2_121661536)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs7591686,a,g,0.7333,-0.0035,0.0071,0.62,67588.6,2_45609171)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs2936687,a,g,0.8167,0.011,0.0082,0.16,67592,8_75279924)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs3780697,t,c,0.725,0.0015,0.0081,0.85,67574.2,9_131618930)
Error2a -- line2[0] present more than once in a particular file (count: 8 line2[0]: rs573071,a,g,0.675,-0.0073,0.0066,0.27,67589.6,9_246993)
#~~~

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/process.MTedits.ForGIANT2014_5.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/process.MTedits.ForGIANT2013.vs1.R

#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/process.MTedits.ForGIANT2013.vs1.R into an interactive sessions in $PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/RData.GIANT2013.process.20150602
gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/RData.GIANT2013.process.20150602

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.txt --file2 <(zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_MEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'my $chr; my $pos; if (($F[3] eq "NA") || ($F[3] eq "None")) { $chr = "NA"; $pos = "NA"; } else { $chr = ((split(/_/, $F[3]))[0]); $pos = ((split(/_/, $F[3]))[1]); } splice(@F, 1, 0, ($chr, $pos)); splice(@F, 8, 1); splice(@F, 11, 1); splice(@F, 14, 1); splice(@F, 17, 1); splice(@F, 20, 1); splice(@F, 23, 1); splice(@F, 27, 1); splice(@F, 30, 1); splice(@F, 33, 1); print join("\t", @F);') | \ 
gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz

#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz
#snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI  beta_WHR se_WHR n_WHR beta_HIPadjBMI se_HIPadjBMI n_HIPadjBMI beta_HIP se_HIP n_HIP beta_WCadjBMI se_WCadjBMI n_WCadjBMI beta_WC se_WC n_WC beta_weight se_weight n_weight Z.height Z.BMI Z.WHRadjBMI Z.WHR Z.HIPadjBMI Z.HIP Z.WCadjBMI Z.WC Z.weight mvstat mvp unip

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.txt | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | perl -F, -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | head -n 5

~~~
[  mturchin20@bigmem03  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.txt | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | perl -F, -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | head -n 5
7624    7624 1638075
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | wc
7624  358328 3362517
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
7453  350291 3287270
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | perl -F, -lane 'print $#F;' | sort | uniq -c
7624 0
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Male.dtlesssignif.vs1.annot.vs1.txt.gz | head -n 5
snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI  beta_WHR se_WHR n_WHR beta_HIPadjBMI se_HIPadjBMI n_HIPadjBMI beta_HIP se_HIP n_HIP beta_WCadjBMI se_WCadjBMI n_WCadjBMI beta_WC se_WC n_WC beta_weight se_weight n_weight Z.height Z.BMI Z.WHRadjBMI Z.WHR Z.HIPadjBMI Z.HIP Z.WCadjBMI Z.WC Z.weight mvstat mvp unip
rs13078872      3       142732650       c       g       3_142732650     0.3417  0       -0.039  0.0066  60556.4 -0.0022 0.0068  58618.3 0.0062  0.0072  34600.8 -0.00020
0.0081  34580.6 -0.031  0.0085  32929.5 -0.029  0.0089  32848.2 0.3417  0.0050  38369   0.3417  0.0080  38304.4 0.3417  0.0069  58322.9 5.92727656693212       0.318639363964375        0.859617364241911       0.0250689082587111      3.64250283822469        3.26361637006925        2.12007168974215        2.0641868904004 3.09023230616781        64.3274543043609        12.4478658862156        8.51144928349955
rs545708        18      55991198        t       c       18_55991198     0.4667  0       0.016   0.0061  60540.8 0.039   0.0064  58617.7 -0.0038 0.0068  34600.4 0.022  0.0076   34580.2 0.0044  0.0079  32929.3 0.038   0.0083  32848.1 0.4667  0.0047  38368.6 0.4667  0.0075  38304   0.4667  0.0064  58322.3 2.5758293035489 6.11944811150025
0.568051498338983       2.82690683268606        0.553384719555673       4.58684245215462        0.125661346855074       5.01259803621023        7.00499197014344
78.0663903750382        15.3494810269617        11.6073030467403
rs10484733      6       142752680       c       g       6_142752680     0.1917  0       0.035   0.0077  60398.8 -0.012  0.0080  58377.7 -0.010  0.0083  34601.2 -0.021 0.0095   34581   0.026   0.010   32930.4 0.0012  0.010   32849.2 0.1917  0.0060  38128.4 0.1917  0.0093  38063.8 0.1917  0.0081  58082.3 4.56311578837316        1.51410188761928        1.22652812003661        2.21151780918668        2.63961605297416        0.113038540644565       2.38670773449225        0.789191652658222       1.17498679206609        37.1300181865462        6.77120095160744        5.29756946355448
rs17514738      1       217055376       t       c       1_217055376     NA      0       0.024   0.0063  60543.8 0.0044  0.0065  58615.7 -0.0023 0.0069  34601.1 0.0029 0.0077   34580.9 0.026   0.0080  32929   0.030   0.0085  32847.8 NA      0.0048  38369.3 NA      0.0075  38304.7 NA      0.0065  58320.3 3.89722801810903        0.674489750196082       0.331853346436817       0.371856089385075       3.19465105376329        3.48075640434621        2.10835839916911        2.45726339020544        2.36561812686429        40.3777497032728        7.44180904800661        4.01188715973165
~~~

python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2010/Vs1/DataOrganization.GIANT2010.PythonVLookUp.PutTogetherProcessAndAnnotFiles.vs1.py --file1 /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.txt --file2 <(zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz | perl -lane 'my $chr; my $pos; if (($F[3] eq "NA") || ($F[3] eq "None")) { $chr = "NA"; $pos = "NA"; } else { $chr = ((split(/_/, $F[3]))[0]); $pos = ((split(/_/, $F[3]))[1]); } splice(@F, 1, 0, ($chr, $pos)); splice(@F, 8, 1); splice(@F, 11, 1); splice(@F, 14, 1); splice(@F, 17, 1); splice(@F, 20, 1); splice(@F, 23, 1); splice(@F, 27, 1); splice(@F, 30, 1); splice(@F, 33, 1); print join("\t", @F);') | \ 
gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz

#Added the following header to /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz
#snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI  beta_WHR se_WHR n_WHR beta_HIPadjBMI se_HIPadjBMI n_HIPadjBMI beta_HIP se_HIP n_HIP beta_WCadjBMI se_WCadjBMI n_WCadjBMI beta_WC se_WC n_WC beta_weight se_weight n_weight Z.height Z.BMI Z.WHRadjBMI Z.WHR Z.HIPadjBMI Z.HIP Z.WCadjBMI Z.WC Z.weight mvstat mvp unip

cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.txt | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | perl -F, -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | head -n 5

~~~
[  mturchin20@bigmem03  ~]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.txt | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | perl -F, -lane 'print $#F;' | sort | uniq -c
zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | head -n 5
9437    9437 2032741
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | wc
9437  443539 4172992
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | perl -lane 'if ($F[$#F-1] > $F[$#F]) { print join("\t", @F) ; } ' | wc
8725  410075 3855206
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | perl -F, -lane 'print $#F;' | sort | uniq -c
9437 0
[  mturchin20@bigmem03  ~]$zcat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT2013.Female.dtlesssignif.vs1.annot.vs1.txt.gz | head -n 5
snp chr pos a1 a2 chrbp maf annot beta_height se_height n_height beta_BMI se_BMI n_BMI beta_WHRadjBMI se_WHRadjBMI n_WHRadjBMI  beta_WHR se_WHR n_WHR beta_HIPadjBMI se_HIPadjBMI n_HIPadjBMI beta_HIP se_HIP n_HIP beta_WCadjBMI se_WCadjBMI n_WCadjBMI beta_WC se_WC n_WC beta_weight se_weight n_weight Z.height Z.BMI Z.WHRadjBMI Z.WHR Z.HIPadjBMI Z.HIP Z.WCadjBMI Z.WC Z.weight mvstat mvp unip
rs13078872      3       142732650       c       g       3_142732650     0.3417  0       -0.037  0.0061  73130.9 -0.0033 0.0065  67956.5 -0.0068 0.0071  42734.7 -0.00750.0076   42731.4 -0.0095 0.0079  40723.9 -0.011  0.0085  40354.7 0.3417  0.0049  47471   0.3417  0.0075  47319.9 0.3417  0.0066  67592.2 6.14139478519449        0.510073456968595       0.954165253146194       0.994457883209753       1.20035885803086        1.25356543847045        1.88079360815125        1.51410188761928        2.89430405305143        52.1090569513165        9.88309797640845        9.08724669632868
rs4360313       8       4802188 t       g       8_4802188       0.25    0       0.025   0.0060  73125   -0.020  0.0064  67956.6 0.0033  0.0070  42735.1 -0.0022 0.0075 42731.8  0.016   0.0078  40724.7 -0.0029 0.0084  40355.5 0.25    0.0048  47471.4 0.25    0.0074  47320.3 0.25    0.0065  67592.3 4.19406451274716        3.12138914935986
0.467698799114508       0.292374896226804       1.97736842818195        0.345125531470472       1.31057911216813        0.994457883209753       0.877896295051229       37.0427275273345        6.75321592565396        4.56224943717961
rs17378391      20      47206007        t       c       20_47206007     0.3083  2       -0.032  0.0067  73116.3 0.017   0.0073  67956.9 -0.0034 0.0079  42735.1 0.000600.0084   42731.8 -0.014  0.0088  40724.9 0.0088  0.0093  40355.7 0.3083  0.0054  47471.4 0.3083  0.0083  47320.3 0.3083  0.0073  67592.6 4.75647768915719        2.30798447494596        0.426148007841278       0.0752698620998299      1.59819313992282        0.954165253146194       2.58280745200824        0.10043372051147        0.125661346855074       43.391819913533 8.06646637431228        5.70553377383841
rs545708        18      55991198        t       c       18_55991198     0.4667  0       0.012   0.0057  73118.7 0.030   0.0062  67956.3 0.010   0.0067  42734.5 0.026  0.0071   42731.2 0.0037  0.0074  40723.9 0.031   0.0080  40354.7 0.4667  0.0046  47470.8 0.4667  0.0070  47319.7 0.4667  0.0062  67592   2.15707270447901        4.79215420637313        1.55477359459685        3.58274690211504        0.495850347347453       3.89747712566571        1.05812161768478        4.66474498631497        4.98581864328374        54.6991692870824        10.4252269927113        6.20971483596676
~~~

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/GLC.MTedits.ForGIANT2014_5.AllPheno.vs2.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GLC.MTedits.ForGIANT2013.Male.vs1.R
#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GLC.MTedits.ForGIANT2013.Male.vs1.R into an inteactve session in #PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/

#CHECK_0: Redo this part in GIANT2014_5.AllPheno -- think I may have mesed it up; the nmin & l==1 portion
#lbf.newhits <- rbind(lbf.bigmat1[,gl$nmin>30000][,l==1],lbf.bigmat2[,gl$nmin>30000][,l==1],lbf.bigmat3[,gl$nmin>30000][,l==1],lbf.bigmat4[,gl$nmin>30000][,l==1],lbf.bigmat5[,gl$nmin>30000][,l==1],lbf.bigmat6[,gl$nmin>30000][,l==1],lbf.bigmat7[,gl$nmin>30000][,l==1])

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GLC.MTedits.ForGIANT2013.Male.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GLC.MTedits.ForGIANT2013.Female.vs1.R
#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GLC.MTedits.ForGIANT2013.Female.vs1.R into an inteactve session in #PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/

~~~
[  mturchin20@spudhead  ~/Data/dbSNP]$cat /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GLC.MTedits.ForGIANT2013.Male.vs1.newtophits.vs1.txt /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GLC.MTedits.ForGIANT2013.Female.vs1.newtophits.vs1.txt | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 >= 2) { print $0 } } '
2 rs1421085
2 rs143384
2 rs1812175
2 rs2145270
2 rs2280470
2 rs724016
2 rs7498665
2 rs7759938
2 rs806794
2 snp
~~~

#CHECK_0: Make a consistent criteria, or as much as possible, with what SNPs are included as 'GWAS/Hit SNPs', e.g. do you include replication findings or not, just primary analyses/first-round results? etc...?








##PGC 2013
##20150601

mkdir /mnt/gluster/data/external_public_supp/PGC2013 
cd /mnt/gluster/data/external_public_supp/PGC2013
#Accessed 20150601
wget https://www.med.unc.edu/pgc/files/resultfiles/pgc.cross.full.2013-03.zip
wget https://www.med.unc.edu/pgc/files/resultfiles/pgc.cross.add.zip
wget https://www.med.unc.edu/pgc/files/resultfiles/pgc.cross.aut.zip
wget https://www.med.unc.edu/pgc/files/resultfiles/pgc.cross.bip.zip
wget https://www.med.unc.edu/pgc/files/resultfiles/pgc.cross.mdd.zip
wget https://www.med.unc.edu/pgc/files/resultfiles/pgc.cross.scz.zip
wget https://www.med.unc.edu/pgc/files/resultfiles/pgc-cross-scoring-5files.zip

unzip /mnt/gluster/data/external_public_supp/PGC2013/pgc.cross.full.2013-03.zip
unzip /mnt/gluster/data/external_public_supp/PGC2013/pgc.cross.add.zip
unzip /mnt/gluster/data/external_public_supp/PGC2013/pgc.cross.aut.zip
unzip /mnt/gluster/data/external_public_supp/PGC2013/pgc.cross.bip.zip
unzip /mnt/gluster/data/external_public_supp/PGC2013/pgc.cross.mdd.zip
unzip /mnt/gluster/data/external_public_supp/PGC2013/pgc.cross.scz.zip
unzip /mnt/gluster/data/external_public_supp/PGC2013/pgc-cross-scoring-5files.zip

#~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$cat pgc.cross.ADD4.2013-05.txt | head -n 5
snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
rs3131972       1       742584  A       G       1.137   0.1189  0.2799  0.731   0       0.16055
rs3131969       1       744045  A       G       1.13    0.1136  0.2813  0.966   0       0.133028
rs3131967       1       744197  T       C       1.14    0.1208  0.2794  0.899   0       .
rs1048488       1       750775  T       C       0.8794  0.1189  0.2798  0.731   0       0.836449
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$cat pgc.cross.AUT8.2013-05.txt | head -n 5
snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
rs3131972       1       742584  A       G       1.03    0.0478  0.5351  0.972   0       0.16055
rs3131969       1       744045  A       G       1.067   0.0538  0.2268  0.926   0       0.133028
rs3131967       1       744197  T       C       1.082   0.0579  0.1754  0.854   0       .
rs1048488       1       750775  T       C       0.9829  0.048   0.7195  0.974   0       0.836449
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$cat pgc.cross.BIP11.2013-05.txt | head -n 5
snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
rs3131972       1       742584  A       G       1.092   0.0817  0.2819  0.694   0       0.16055
rs3131969       1       744045  A       G       1.087   0.0781  0.2855  0.939   0       0.133028
rs3131967       1       744197  T       C       1.093   0.0835  0.2859  0.869   0       .
rs1048488       1       750775  T       C       0.9158  0.0817  0.2817  0.694   0       0.836449
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$cat pgc.cross.MDD9.2013-05.txt | head -n 5
snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
rs12562034      1       758311  A       G       0.9028  0.0894  0.2527  1.02    0       0.0925926
rs4970383       1       828418  A       C       1.018   0.0769  0.8186  0.436   0       0.201835
rs4475691       1       836671  T       C       1.016   0.0546  0.7751  0.979   0       0.146789
rs1806509       1       843817  A       C       0.99    0.0737  0.8911  0.372   0       0.600917
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$cat pgc.cross.SCZ17.2013-05.txt | head -n 5
snpid hg18chr bp a1 a2 or se pval info ngt CEUaf
rs3131972       1       742584  A       G       1       0.0966  0.9991  0.702   0       0.16055
rs3131969       1       744045  A       G       1       0.0925  0.9974  0.938   0       0.133028
rs3131967       1       744197  T       C       1.001   0.0991  0.9928  0.866   0       .
rs1048488       1       750775  T       C       0.9999  0.0966  0.9991  0.702   0       0.836449
#~~~

cd /mnt/lustre/home/mturchin20/Software/
wget ftp://ftp.foolabs.com/pub/xpdf/xpdfbin-linux-3.04.tar.gz
tar -xvzf /mnt/lustre/home/mturchin20/Software/xpdfbin-linux-3.04.tar.gz

#SNP hits
#Copying/pasting/downoading/finding (however) list of top hits from each study

mkdir /mnt/lustre/home/mturchin20/Data/PGC/
mkdir /mnt/lustre/home/mturchin20/Data/PGC/2013/

#ADD -- no significant loci found (that passed GWS)


#Autism -- unavailable yet (no paper?)


#Bipolar -- (2 SNPs that had primary signals and kept them during replication; more SNPs reached GWS after replication but for now don't want to include SNPs based on 'replication' results. Dropping the 2 other primary signals since we at least know they were potential FPs?)

#Copy/pasted table 2 & 3 from http://www.nature.com.proxy.uchicago.edu/ng/journal/v43/n10/fig_tab/ng.943_T3.html into /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.txt & /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table3.txt
#Manually know that rs10994397 and rs9371601 are the two SNPs that meet the paper's criteria -- "Only two of the four SNPs listed in Table 2 had replication P < 0.05; the genome-wide significant SNPs rs10994397 and rs9371601 did not have P < 0.05 (replication P = 0.12 and P = 0.10, respectively)."

cat /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.txt | grep -E 'rs10994397|rs9371601' | awk '{ print $1, "\t", $2, "\t", $3 }' | sed 's/,//g' > /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.GWSignificant.noCommas.MarkerChrBP.txt


#Cross-study -- 4 SNPs/loci
#Copy/pasted table 1 from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3714010/table/T1/ (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3714010/) into /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.txt

cat /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.txt | grep rs | awk '{ print $1, "\t", $2, "\t", $3 }' > /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.MarkerChrBP.txt


#MDD (CrossMDDBipolar) -- 15 SNPs
#Copy/pasted SNPs with pvals "...In the MDD-bipolar cross-disorder analysis, 15 SNPs exceeded genome-wide significance (P<5 × 10(-8)), and all were in a 248 kb interval of high LD on 3p21.1 (chr3:52 425 083-53 822 102, minimum P=5.9 × 10(-9) at rs2535629)..." in Supplemental Table 19 from http://www.nature.com.proxy.uchicago.edu/mp/journal/v18/n4/full/mp201221a.html (http://www.nature.com.proxy.uchicago.edu/mp/journal/v18/n4/suppinfo/mp201221s1.html) into /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.txt

cat /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.txt | awk '{ print $1, "\t", $2, "\t", $3 }' > /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.MarkerChrBP.txt


#SCZ -- 10 SNPs, 3 SNPs w/ joint BP
#Copy/pasted supplementary tables 4 and 11 from http://www.nature.com.proxy.uchicago.edu/ng/journal/v43/n10/full/ng.940.html#supplementary-information into /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable4.txt and /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.txt
(also pasted table 2 into /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.Table2.txt but I don't think I'm going to use this form of the results)

cat /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable4.txt | awk '{ if ($4 <= 0.00000005) { print $0 } }' > /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable4.GWSignificant.txt
cat /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable4.GWSignificant.txt | awk '{ print $1, "\t", $2, "\t", $3 }' | sed 's/,//g' > /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable4.GWSignificant.noCommas.MarkerChrBP.txt
cat /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.txt | awk '{ print $1, "\t", $2, "\t", $3 }' | sed 's/,//g' > /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.noCommas.MarkerChrBP.txt

cat /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.GWSignificant.noCommas.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.noCommas.MarkerChrBP.txt | grep rs | wc
cat /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.GWSignificant.noCommas.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.noCommas.MarkerChrBP.txt | grep rs | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 > 1) { print $0 } } '
cat /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.GWSignificant.noCommas.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.noCommas.MarkerChrBP.txt | grep rs | sort | uniq >  /mnt/lustre/home/mturchin20/Data/PGC/2013/PGC2013.Hits.MarkerChrBP.txt

#~~~
[  mturchin20@spudhead  ~/Data/dbSNP]$cat /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.GWSignificant.noCommas.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.noCommas.MarkerChrBP.txt | grep rs | wc
     24      72     616
[  mturchin20@spudhead  ~/Data/dbSNP]$cat /mnt/lustre/home/mturchin20/Data/PGC/2013/BP2011.Table2.GWSignificant.noCommas.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/Cross2014.Table1.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/MDD2013.SupplTable19.MarkerChrBP.txt /mnt/lustre/home/mturchin20/Data/PGC/2013/SCZ2013.SupplTable11.noCommas.MarkerChrBP.txt | grep rs | awk '{ print $1 }' | sort | uniq -c | awk '{ if ($1 > 1) { print $0 } } '
      2 rs2535629
[  mturchin20@spudhead  ~/Data/dbSNP]$cat /mnt/lustre/home/mturchin20/Data/PGC/2013/PGC2013.Hits.MarkerChrBP.txt | wc
     23      69     591
#~~~

mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC/2013
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC/2013/Vs1

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/process.MTedits.ForGIANT2013.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/process.MTedits.ForPGC2013.vs1.R

#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/process.MTedits.ForPGC2013.vs1.R into an interactive sessions in $PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/RData.PGC2013.Cross.process.20150602
gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/RData.PGC2013.Cross.process.20150602

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/DataOrganization.GIANT2013.PythonVLookUp.CreateFullAnnotFile.vs1.py /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/DataOrganization.PGC2013.PythonVLookUp.CreateFullAnnotFile.vs1.py



python /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/DataOrganization.GIANT2013.PythonVLookUp.CreateFullAnnotFile.vs1.py --file1 /mnt/lustre/home/mturchin20/Data/GIANT/2013/journal.pgen.1003500.s007.noCarriageReturn.Female.MarkerChrBP.txt --file2 /mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_BMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHRadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WHR_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIPadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HIP_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WCadjBMI_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WC_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz,/mnt/gluster/data/external_public_supp/GIANT2013/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_WEIGHT_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.txt.gz | gzip > /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_AllPheno_WOMEN_N.wUCSCGenomeBrowser_dbSNP130.vs1.Annot.vs1.txt.gz

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2014_5/Vs1/process.MTedits.ForGIANT2014_5.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/process.MTedits.ForGIANT2013.vs1.R

#Copy/pasted commands from /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/process.MTedits.ForGIANT2013.vs1.R into an interactive sessions in $PWD /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/RData.GIANT2013.process.20150602
gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GIANT2013/Vs1/RData.GIANT2013.process.20150602





#CHECK_0: Make a consistent criteria, or as much as possible, with what SNPs are included as 'GWAS/Hit SNPs', e.g. do you include replication findings or not, just primary analyses/first-round results? etc...?

##TAG 2010
##20150601

mkdir /mnt/gluster/data/external_public_supp/TAG2010
cd /mnt/gluster/data/external_public_supp/TAG2010
#Accessed 20150601
wget https://www.med.unc.edu/pgc/files/resultfiles/readme.tag.txt
wget https://www.med.unc.edu/pgc/files/resultfiles/tag.cpd.tbl.gz
wget https://www.med.unc.edu/pgc/files/resultfiles/tag.evrsmk.tbl.gz
wget https://www.med.unc.edu/pgc/files/resultfiles/tag.former.tbl.gz
wget https://www.med.unc.edu/pgc/files/resultfiles/tag.logonset.tbl.gz

#~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/TAG2010]$zcat tag.logonset.tbl.gz | head -n 5
CHR     SNP     BP      A1   A2 FRQ_A   FRQ_U   INFO    OR      SE      P
1       rs12565286      711153  C       G       0.0434  0.0434  1       -0.0054 0.0123  0.6592
1       rs11804171      713682  A       T       0.0434  0.0434  1       -0.0049 0.0123  0.6912
1       rs3094315       742429  A       G       0.8211  0.8211  1       -0.0046 0.0042  0.2692
1       rs3131968       744055  A       G       0.1462  0.1462  1       -0.0141 0.0101  0.1626
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/TAG2010]$zcat tag.former.tbl.gz | head -n 5
CHR     SNP     BP      A1   A2 FRQ_A   FRQ_U   INFO    OR      SE      P
1       rs12565286      711153  C       G       0.0457  0.0457  1       0.0533  0.0934  0.568
1       rs11804171      713682  A       T       0.0457  0.0457  1       0.0553  0.0937  0.5547
1       rs3094315       742429  A       G       0.8238  0.8238  1       -0.0479 0.0316  0.1301
1       rs3131968       744055  A       G       0.1346  0.1346  1       0.0912  0.0529  0.08477
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/TAG2010]$zcat tag.evrsmk.tbl.gz | head -n 5
CHR     SNP     BP      A1   A2 FRQ_A   FRQ_U   INFO    OR      SE      P
1       rs12565286      711153  C       G       0.241   0.241   1       -0.0885 0.0736  0.2293
1       rs11804171      713682  A       T       0.2394  0.2394  1       -0.089  0.0739  0.229
1       rs3094315       742429  A       G       0.7737  0.7737  1       -0.0031 0.0237  0.8949
1       rs3131968       744055  A       G       0.1349  0.1349  1       -0.021  0.0435  0.6293
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/TAG2010]$zcat tag.cpd.tbl.gz | head -n 5
CHR     SNP     BP      A1   A2 FRQ_A   FRQ_U   INFO    OR      SE      P
1       rs12565286      711153  C       G       0.0434  0.0434  1       -0.4029 0.5634  0.4745
1       rs11804171      713682  A       T       0.0434  0.0434  1       -0.4092 0.5649  0.4689
1       rs3094315       742429  A       G       0.824   0.824   1       -0.035  0.174   0.8407
1       rs3131968       744055  A       G       0.1348  0.1348  1       0.2606  0.3068  0.3956
#~~~

#SNP hits

mkdir /mnt/lustre/home/mturchin20/Data/TAG/
mkdir /mnt/lustre/home/mturchin20/Data/TAG/2010/

#Copy/pasted table 2 from http://www.nature.com.proxy.uchicago.edu/ng/journal/v42/n5/fig_tab/ng.571_T2.html (http://www.nature.com.proxy.uchicago.edu/ng/journal/v42/n5/full/ng.571.html) into /mnt/lustre/home/mturchin20/Data/TAG/2010/TAG2010.Table2.txt

cat /mnt/lustre/home/mturchin20/Data/TAG/2010/TAG2010.Table2.txt | grep rs | awk '{ print $1 }' | grep rs > /mnt/lustre/home/mturchin20/Data/TAG/2010/TAG2010.Table2.MarkerName.txt

zcat /mnt/gluster/data/external_public_supp/TAG2010/tag.*.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/TAG/2010/TAG2010.Table2.MarkerName.txt | awk '{ print $1, "\t", $2, "\t", $3 }' | sort | uniq -c | awk '{ print $3, "\t", $2, "\t", $4 }' > /mnt/lustre/home/mturchin20/Data/TAG/2010/TAG2010.Table2.MarkerChrBP.txt

~~~
[  mturchin20@spudhead  ~/Data/dbSNP]$cat /mnt/lustre/home/mturchin20/Data/TAG/2010/TAG2010.Table2.MarkerName.txt | wc
     14      14     137
[  mturchin20@spudling50  ~]$zcat /mnt/gluster/data/external_public_supp/TAG2010/tag.*.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/TAG/2010/TAG2010.Table2.MarkerName.txt | awk '{ print $1, "\t", $2, "\t", $3 }' | sort | uniq -c
4 10       rs1028936       93339777
4 10       rs1329650       93338100
4 11       rs1013442       27535522
4 11       rs1304100       27528179
4 11       rs4074134       27603861
4 11       rs4923457       27605156
4 11       rs4923460       27613365
4 11       rs6265          27636492
4 11       rs6484320       27659764
4 11       rs879048        27595510
4 15       rs1051730       76681394
4 15       rs16969968      76669980
4 19       rs3733829       46002411
4 9        rs3025343       135468176
~~~     

mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/TAG2010 
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/TAG2010/Vs1 

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/PGC2013/Vs1/process.MTedits.ForPGC2013.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/TAG2010/Vs1/process.MTedits.ForTAG2010.vs1.R 

mv /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/TAG2010/Vs1/.RData /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/TAG2010/Vs1/RData.TAG2010.process.20150602
gzip /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/TAG2010/Vs1/RData.TAG2010.process.20150602











##IBD 2012
##20150601

mkdir /mnt/gluster/data/external_public_supp/IBDConsortium2012
cd /mnt/gluster/data/external_public_supp/IBDConsortium2012
#Accessed 20150601
wget ftp://ftp.sanger.ac.uk/pub4/ibdgenetics/gwas_ichip_meta_release.txt.gz
#Accessed 20150602
wget ftp://ftp.sanger.ac.uk/pub4/ibdgenetics/cd-meta.txt.gz
wget ftp://ftp.sanger.ac.uk/pub4/ibdgenetics/ucmeta-sumstats.txt.gz


#~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/IBD2012]$zcat gwas_ichip_meta_release.txt.gz | head -n 10
CHR BP GWAS_SNP GWAS_A1A2 ICHIP_SNP ICHIP_A1A2 LD-info Strand/Effectallele CD_GWAS_P CD_ICHIP_P CD_META_P UC_GWAS_P UC_ICHIP_P UC_META_P IBD_GWAS_P IBD_ICHIP_P IBD_META_P MAF CD_OR CD_OR_CI UC_OR UC_OR_CI IBD_OR IBD_OR_CI
1 1110294 rs1320571 AG same AG same same_strand,_same_effect-allele - - - - - - 0.004646 0.8236 0.105 0.0442 - - - - 1.05 0.989652-1.11403
1 1125105 rs9729550 AC same CA same same_strand,_different_effect-allele 0.001286 0.03004 0.000316 - - - - - - 0.738 0.959692898272553 0.924956-0.995735 - - - -
1 1130298 rs1815606 TG same AC same strand_flip,_same_effect-allele 0.0001352 0.04842 0.0002186 - - - 0.00123 0.009144 6.327e-05 0.31 1.036 1.00062-1.07263 - - 1.051 1.02556-1.07707
1 1142494 rs11721 AC rs6697886 AG (0.485,AA/CG) same_strand,_same_effect-allele 0.003183 0.0003092 5.165e-06 - - - 0.004786 1.81e-06 2.954e-08 0.0781 1.089 1.03957-1.14078 - - 1.102 1.06484-1.14045
1 1153667 rs7515488 TC same AG same strand_flip,_same_effect-allele 0.000283 0.0002496 3.62e-07 0.006893 1.753e-05 4.92e-07 6.046e-05 3.443e-06 7.67e-10 0.147 1.088 1.04018-1.13801 1.11 1.05823-1.16431 1.1 1.06708-1.13393
1 1155173 rs11260562 AG same AG same same_strand,_same_effect-allele 0.0003004 0.08256 0.0004346 0.005253 0.04672 0.001047 1.053e-05 0.02667 6.605e-06 0.0546 1.063 0.992095-1.13897 1.077 1.00084-1.15896 1.114 1.06302-1.16743
1 1163474 rs6697886 AG same AG same same_strand,_same_effect-allele 0.00109 0.0003092 1.57e-06 - - - 0.0003566 1.81e-06 2.441e-09 0.131 1.089 1.03957-1.14078 - - 1.101 1.06659-1.13652
1 1164145 rs4970364 TC rs11804831 GA (0.627,CC/TT) strand_flip,_different_effect-allel 6.571e-05 0.0005104 2.77e-07 0.004123 5.014e-07 1.243e-08 5.14e-05 9.839e-07 2.326e-10 0.76 0.928505106778087 0.890662-0.967956 0.892857142857143 0.854237-0.933223 0.9138 0.888718-0.93959
1 1166460 rs6675798 TC rs6697886 AG (0.657,CA/TG) same_strand,_same_effect-allele - - - - - - 0.001967 1.81e-06 1.267e-08 0.9 - - - - 0.9069 0.876836-0.937995
#~~~

#SNP hits -- 163 SNPs

mkdir /mnt/lustre/home/mturchin20/Data/IBDConsortium/
mkdir /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/

#Copy/pasted supplementary table 2 from http://www.nature.com.proxy.uchicago.edu/nature/journal/v491/n7422/extref/nature11582-s2.zip (http://www.nature.com.proxy.uchicago.edu/nature/journal/v491/n7422/full/nature11582.html) into /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.txt

cat /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.txt | grep rs | awk '{ print $3 }' | perl -lane 'if ($F[0] =~ m/(rs\d+).*/) { print $1; }' > /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt

~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$cat /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | wc
    163     163    1643
~~~    

zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | perl -lane 'print join("\t", @F[0..2]);' | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt > /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt.pt1
zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | perl -lane 'print $F[2], "\t", $F[0], "\t", $F[1];' > /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerChrBP.txt

#I think 2 top SNPs might be missing from the maine data file? I get 161 SNPs based on the GWAS SNP column and when I try to identify the '2 SNPs missing that should be in the ICHIP column' they don't exist

~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | perl -lane 'print join("\t", @F[0..2]);' | wc
    225     675    4900
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | perl -lane 'print join("\t", @F[0..2]);' | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | wc
    161     483    3508
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | grep -w -v -f <(cat /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt.pt1 | awk '{ print $3 }') | wc
      0       0       0
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | awk '{ print $5 }' | grep -w -v -f <(cat /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt.pt1 | awk '{ print $3 }') | wc
    161     161    1008
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | awk '{ print $5 }' | grep -w -v -f <(cat /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt.pt1 | awk '{ print $3 }') | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/PGC2013]$zcat /mnt/gluster/data/external_public_supp/IBDConsortium2012/gwas_ichip_meta_release.txt.gz | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | perl -ane 'print $F[2], "\n", $F[4], "\n";' | grep -w -f /mnt/lustre/home/mturchin20/Data/IBDConsortium/2012/Jostins2012.SupplTable2.noCrosses.MarkerNames.txt | sort | uniq -c | wc
    161     322    2911
~~~

#Got links for the updated IDBGC information from Stephan Ripke & co -- see e-mail change copy/pasted below near end of this record file
#Accessed 20150618

mkdir /mnt/gluster/data/external_public_supp/IBDConsortium2015
cd /mnt/gluster/data/external_public_supp/IBDConsortium2015

wget http://www.broadinstitute.org/~sripke/share_links/ALHySN0Z6i4O6PZeQgpqXN9Ewtqxit_IBD15_1KG.sh2/daner_IBD15_1KG.sh2_mds10.CD7b.gz
wget http://www.broadinstitute.org/~sripke/share_links/ALHySN0Z6i4O6PZeQgpqXN9Ewtqxit_IBD15_1KG.sh2/daner_IBD15_1KG.sh2_mds10.UC8b.gz
wget http://www.broadinstitute.org/~sripke/share_links/ALHySN0Z6i4O6PZeQgpqXN9Ewtqxit_IBD15_1KG.sh2/daner_IBD15_1KG.sh2_mds15.IBD15bb.gz
wget http://www.broadinstitute.org/~sripke/share_links/ALHySN0Z6i4O6PZeQgpqXN9Ewtqxit_IBD15_1KG.sh2/download_all.txt

~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/IBDConsortium2015]$zcat daner_IBD15_1KG.sh2_mds10.CD7b.gz | head -n 10
CHR     SNP     BP      A1   A2 FRQ_A_5956      FRQ_U_14927     INFO    OR      SE      P       ngt     Direction       HetISqt HetChiSq        HetDf   HetPVa
10      rs185339560     2392426 T       C       0.00954 0.00942 0.633   0.91411 0.1485  0.5453  0       +----++ 22.5    10.316  6       0.112
10      rs11250701      1689546 A       G       0.679   0.669   0.962   1.02398 0.0248  0.3394  0       -++---- -14.0   7.018   6       0.3191
10      chr10_2622752_D 2622752 I2      D       0.97    0.968   0.929   1.06184 0.0691  0.3857  0       -+-+--- 37.9    12.884  6       0.04492
10      rs7085086       151476  A       G       0.319   0.317   0.935   0.99452 0.0253  0.8283  0       ---+-++ 0.0     4.809   6       0.5685
10      rs113494187     1593759 T       G       0.986   0.985   0.829   1.00391 0.1064  0.9708  0       --++++- -30.4   6.136   6       0.4081
10      rs117915320     1708106 A       C       0.0113  0.012   0.616   1.01227 0.1373  0.9292  0       -+--+++ 0.0     5.392   6       0.4946
10      rs182753344     790310  T       C       0.0993  0.0996  0.606   1.04561 0.0489  0.3612  0       -++++++ 0.0     2.893   6       0.8221
10      rs188913771     1273049 A       G       0.0126  0.0131  0.672   0.88577 0.1253  0.3331  0       -----++ -7.5    7.443   6       0.2818
10      rs7911665       2067236 T       G       0.681   0.686   0.915   0.98442 0.0254  0.5368  0       ---+-+- 0.0     1.667   6       0.9476
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/IBDConsortium2015]$zcat daner_IBD15_1KG.sh2_mds10.UC8b.gz | head -n 10
CHR     SNP     BP      A1   A2 FRQ_A_6968      FRQ_U_20464     INFO    OR      SE      P       ngt     Direction       HetISqt HetChiSq        HetDf   HetPVa
10      rs185339560     2392426 T       C       0.00802 0.00774 0.626   1.03076 0.1554  0.8452  0       --+++--?        0.0     3.126   6       0.7928
10      rs11250701      1689546 A       G       0.669   0.668   0.961   1.00833 0.0229  0.7183  0       -+-+--+-        -2.6    8.776   7       0.2692
10      chr10_2622752_D 2622752 I2      D       0.968   0.969   0.93    0.91229 0.063   0.1452  0       -+-+++-+        12.1    10.239  7       0.1754
10      rs7085086       151476  A       G       0.311   0.313   0.971   1.00070 0.0231  0.9747  0       -+-++-++        -19.6   7.523   7       0.3765
10      rs113494187     1593759 T       G       0.987   0.985   0.866   1.07993 0.0979  0.4321  0       ++--++++        0.0     1.801   7       0.97
10      rs117915320     1708106 A       C       0.0125  0.0121  0.642   1.05686 0.1208  0.6472  0       --++-+++        -15.5   7.791   7       0.3514
10      rs182753344     790310  T       C       0.0948  0.0998  0.616   0.94233 0.0458  0.1951  0       ----+++-        0.0     5.427   7       0.608
10      rs188913771     1273049 A       G       0.01    0.0113  0.629   0.96194 0.1349  0.7735  0       ++---+-+        0.0     1.865   7       0.9669
10      rs7911665       2067236 T       G       0.691   0.689   0.92    1.03159 0.0238  0.1909  0       +-+++-+-        23.6    11.785  7       0.1078
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/IBDConsortium2015]$zcat daner_IBD15_1KG.sh2_mds15.IBD15bb.gz | head -n 10
CHR     SNP     BP      A1   A2 FRQ_A_12882     FRQ_U_21770     INFO    OR      SE      P       ngt     Direction       HetISqt HetChiSq        HetDf   HetPVa
10      rs185339560     2392426 T       C       0.00879 0.00848 0.634   0.92090 0.1184  0.4865  0       -----+-+-++--+? 25.3    20.073  13      0.09342
10      rs11250701      1689546 A       G       0.673   0.668   0.957   1.01440 0.0182  0.4314  0       --+++--+-+--++- 0.6     16.102  14      0.3072
10      chr10_2622752_D 2622752 I2      D       0.969   0.968   0.923   0.97385 0.0504  0.5986  0       --+-+-+--+-+--+ 22.7    20.704  14      0.1095
10      rs7085086       151476  A       G       0.315   0.316   0.961   0.99302 0.0184  0.7057  0       ----+-+-+-+-+++ -5.3    15.196  14      0.3649
10      rs113494187     1593759 T       G       0.987   0.985   0.849   1.05085 0.0775  0.5225  0       -+++--+-+++++-+ 0.0     8.947   14      0.8344
10      rs117915320     1708106 A       C       0.0118  0.0119  0.623   1.01990 0.0992  0.8423  0       --+--+-++-+++++ 0.0     10.974  14      0.6881
10      rs182753344     790310  T       C       0.0969  0.0988  0.612   0.98639 0.0361  0.705   0       --++--+-++-++-- 0.0     9.530   14      0.7957
10      rs188913771     1273049 A       G       0.0112  0.0117  0.638   0.90620 0.1001  0.325   0       ----+--+--++-++ 0.0     7.930   14      0.8929
10      rs7911665       2067236 T       G       0.686   0.689   0.917   1.01319 0.0188  0.4844  0       -++---+++++-+-- 0.0     12.066  14      0.601
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/IBDConsortium2015]$zcat daner_IBD15_1KG.sh2_mds10.CD7b.gz | wc
12276507 208700614 1186149792
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/IBDConsortium2015]$zcat daner_IBD15_1KG.sh2_mds10.UC8b.gz | wc
12255264 208339153 1197560812
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/IBDConsortium2015]$zcat daner_IBD15_1KG.sh2_mds15.IBD15bb.gz | wc
12716151 216174232 1348769658
~~~

mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/IBDConsortium2012
mkdir /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/IBDConsortium2012/Vs1

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/TAG2010/Vs1/process.MTedits.ForTAG2010.vs1.R /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/IBDConsortium2012/Vs1/process.MTedits.ForIBDConsortium2012.vs1.R








##Blood Traits 20__
##20150625

mkdir /mnt/gluster/data/external_public_supp/HaemgenRBC2012
cd /mnt/gluster/data/external_public_supp/HaemgenRBC2012

wget --user ega-box-233 --password 88fHetwh ftp://ftp.ega.ebi.ac.uk/EGAD00010000300/*

~~~
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$zcat HaemGenRBC_RBC.txt.gz | head -n 10
SNP      Pval    Chr     Position        Effect_Allele   Non_Effect_Allele       Sample_Size     Beta    SE
rs9373124        5.11e-97     6          135464902       T       C       53337   0.0484          0.0024
rs7776054        3.51e-96     6          135460609       A       G       53403   0.0479          0.0024
rs9389269        4.54e-96     6          135468852       T       C       53402   0.0469          0.0023
rs9389268        5.98e-95     6          135461324       A       G       53403   0.0471          0.0024
rs9376090        9.04e-95     6          135452921       T       C       53537   0.0477          0.0024
rs11759553       1.84e-94     6          135463989       A       T       53398   0.0465          0.0023
rs4895441        3.67e-94     6          135468266       A       G       53384   0.0459          0.0023
rs4895440        3.75e-94     6          135468251       A       T       53400   0.0459          0.0023
rs9402686        4.36e-94     6          135469510       A       G       53400   -0.0469         0.0024
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$zcat HaemGenRBC_PCV.txt.gz | head -n 10
SNP      Pval    Chr     Position        Effect_Allele   Non_Effect_Allele       Sample_Size     Beta    SE
rs9373124        3.869e-22    6          135464902       T       C       52764   0.1815          0.0198
rs9389269        1.966e-21    6          135468852       T       C       52829   0.1737          0.019
rs9402686        2.533e-21    6          135469510       A       G       52827   -0.1742         0.0195
rs11759553       4.364e-21    6          135463989       A       T       52825   0.1727          0.0189
rs4895441        5.04e-21     6          135468266       A       G       52811   0.1706          0.0188
rs4895440        5.451e-21    6          135468251       A       T       52827   0.1701          0.0187
rs7776054        6.178e-21    6          135460609       A       G       52830   0.1763          0.0191
rs9389268        1.32e-20     6          135461324       A       G       52830   0.1724          0.0191
rs9376092        2.143e-20    6          135468837       A       C       48671   -0.173          0.0196
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$zcat HaemGenRBC_MCV.txt.gz | head -n 10
SNP      Pval    Chr     Position        Effect_Allele   Non_Effect_Allele       Sample_Size     Beta    SE
rs9389269        2.57e-109    6          135468852       T       C       57855   -0.6003         0.0275
rs7776054        3.13e-108    6          135460609       A       G       57856   -0.6053         0.0277
rs9373124        4.4e-108     6          135464902       T       C       57790   -0.617          0.0284
rs4895441        1.02e-107    6          135468266       A       G       57837   -0.59   0.0272
rs9376090        1.38e-107    6          135452921       T       C       57990   -0.6047         0.0279
rs9389268        2.47e-107    6          135461324       A       G       57856   -0.6004         0.0276
rs4895440        4.47e-107    6          135468251       A       T       57853   -0.5887         0.0272
rs11759553       1.32e-106    6          135463989       A       T       57851   -0.5926         0.0273
rs9402686        1.68e-106    6          135469510       A       G       57853   0.6032          0.028
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$zcat HaemGenRBC_MCH.txt.gz | head -n 10
SNP      Pval    Chr     Position        Effect_Allele   Non_Effect_Allele       Sample_Size     Beta    SE
rs7776054        6.55e-108    6          135460609       A       G       51452   -0.2399         0.0108
rs9389269        8.02e-108    6          135468852       T       C       51452   -0.2368         0.0108
rs9389268        6.04e-107    6          135461324       A       G       51453   -0.2382         0.0108
rs9373124        1.35e-106    6          135464902       T       C       51388   -0.2436         0.0111
rs4895441        2.9e-106     6          135468266       A       G       51434   -0.2322         0.0106
rs9376090        3.04e-106    6          135452921       T       C       51588   -0.2385         0.0109
rs9402686        1.78e-105    6          135469510       A       G       51450   0.2388          0.011
rs4895440        3.28e-105    6          135468251       A       T       51450   -0.2313         0.0106
rs11759553       6.58e-105    6          135463989       A       T       51448   -0.2331         0.0107
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$zcat HaemGenRBC_MCHC.txt.gz | head -n 10
SNP      Pval    Chr     Position        Effect_Allele   Non_Effect_Allele       Sample_Size     Beta    SE
rs10445033       1.536e-22    16         87367963        A       G       42050   -0.017          0.0038
rs837763         1.928e-22    16         87381230        T       C       37768   -0.0147         0.0037
rs475596         5.644e-21    16         87374452        C       G       42124   -0.0135         0.0038
rs198846         7.458e-21    6          26215442        A       G       56189   0.0207          0.0048
rs9932423        2.263e-19    16         87374350        A       C       41839   0.0153          0.0041
rs2608604        4.232e-19    16         87376922        A       G       37770   0.0126          0.0038
rs2932690        4.584e-18    16         87370454        A       G       42098   -0.0162         0.0041
rs855791         3.107e-17    22         35792882        A       G       41609   -0.0124         0.0037
rs2794720        3.349e-17    6          26195181        C       G       44372   0.0143          0.0037
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$zcat HaemGenRBC_Hb.txt.gz | head -n 10
SNP      Pval    Chr     Position        Effect_Allele   Non_Effect_Allele       Sample_Size     Beta    SE
rs855791         4.645e-40    22         35792882        A       G       46184   -0.0791         0.0063
rs4820268        1.153e-31    22         35799537        A       G       43990   0.0699          0.0064
rs2413450        4.411e-31    22         35800170        T       C       49624   -0.0695         0.0064
rs198846         1.422e-30    6          26215442        A       G       60869   0.091   0.008
rs129128         8.996e-27    6          26233321        T       C       53598   -0.0925         0.0086
rs1799945        3.599e-26    6          26199158        C       G       49746   -0.0938         0.0089
rs198833         8.537e-26    6          26222487        A       G       49747   -0.0945         0.009
rs198851         6.056e-25    6          26212611        T       G       49716   0.0963          0.0093
rs130619         2.031e-24    22         35766294        T       C       55170   0.0586          0.0061
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$ls -lrt | awk '{ print $9}'
HaemGenRBC_RBC.txt.gz
HaemGenRBC_PCV.txt.gz
HaemGenRBC_MCV.txt.gz
HaemGenRBC_MCH.txt.gz
HaemGenRBC_MCHC.txt.gz
HaemGenRBC_Hb.txt.gz
[  mturchin20@spudhead  /mnt/gluster/data/external_public_supp/HaemgenRBC2012]$for i in `ls -lrt | awk '{ print $9}'`; do zcat $i | wc; done
2589455 23305095 178560057
2591080 23319720 179098941
2591133 23320197 179175090
2586785 23281065 178736256
2588876 23299884 178629989
2593079 23337711 179093409
~~~



##20150522
##CHECK_0: Redo basepair positions from dbSNP130 by incrementing w/ +1 since it seems like, based on the one file with ChrBP information from GIANT, I am always off by 1
##From file -- /mnt/gluster/data/external_public_supp/GIANT2014_5/GIANT_2015_WCadjBMI_COMBINED_EUR.wUCSCGenomeBrowser_dbSNP130.vs1.HeaderFix.txt.gz 
~~~
.
.
.
rs2296568       10      104826930       T       C       0.1333  0.018   0.0047  0.00013 229905  10_104826929
rs526850        11      118110337       T       G       0.8667  0.024   0.0064  0.00013 153368  11_118110336
rs10879770      12      73098595        C       G       0.7333  0.018   0.0046  0.00013 153283  12_73098594
rs710840        4       82354254        A       G       0.2241  0.022   0.0057  0.00013 115934  4_82354253
rs709930        14      91600266        C       T       0.725   0.018   0.0048  0.00013 153473  14_91600265
rs17675336      3       72497988        C       G       0.7667  0.02    0.0052  0.00013 152936  3_72497987
rs1562258       5       124282417       T       G       0.4     0.017   0.0045  0.00013 151333  5_124282416
rs2740363       17      567161  T       C       0.6417  -0.018  0.0047  0.00013 151656  17_567160
rs12998432      2       67578928        A       T       0.275   0.018   0.0048  0.00013 153390  2_67578927
rs1668335       7       113647888       T       C       0.3417  0.018   0.0046  0.00013 152388  7_113647887
rs910368        14      91558647        C       A       0.725   0.018   0.0048  0.00013 153428  14_91558646
rs4841504       8       11062073        A       C       0.5083  0.013   0.0035  0.00013 229739  8_11062072
rs26990 5       112814742       C       T       0.225   0.02    0.0053  0.00013 152347  5_112814741
rs7701134       5       64578121        G       T       0.1833  0.02    0.0052  0.00013 153402  5_64578120
rs3806502       2       136004743       C       T       0.8583  0.024   0.0064  0.00013 151190  2_136004742
rs1564347       15      71225901        G       T       0.7083  0.017   0.0044  0.00013 153233  15_71225900
rs10220696      14      91630696        G       A       0.725   0.018   0.0048  0.00013 153474  14_91630695
rs10185077      2       56051616        G       T       0.45    -0.017  0.0045  0.00013 139049  2_56051615
rs12120182      1       225880955       A       G       0.225   -0.021  0.0054  0.00013 153383  1_225880954
rs7145052       14      91530945        C       T       0.5833  0.016   0.0041  0.00013 151554  14_91530944
rs10779417      1       219193537       G       A       0.4492  0.017   0.0044  0.00013 143056  1_219193536
rs196046        6       22202603        T       C       0.65    0.017   0.0045  0.00013 152250  6_22202602
rs1219514       10      123422466       A       G       0.075   -0.029  0.0077  0.00013 153461  10_123422465
rs735998        20      38192244        G       A       0.4083  0.017   0.0044  0.00013 153050  20_38192243
rs1207776       6       22201000        C       T       0.65    0.017   0.0045  0.00013 152250  6_22200999
.
.
.
~~~

##20150615
##Create version of main script that can work on any set of identically formatted files
##CHECK_0: Need to think about how to traverse the different phenotypes of interest -- ask user for a line of phenotype names comma seperated that will be in header of file?

##20150701
#Created a 'vs2' of this MainScript1.sh file because reorganizing large chunks of it so that 'file downloads', 'top SNP hits' and 'main analyses' are done in separate sections all together for organizational purposes
#Think this might make things a bit easier to traverse through for the time being or at least organize mentally in my head
#Also somewhat purposely separating the Global Lipids and GIANT work from the rest of the datasets since it seems/feels like those two datasets (or groups of datasets) can be the anchors of the project whereas most other datasets will be accessories/supportive; also this is how my mind is naturally separating things out at the moment anyways

cp -p /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/MainScript1.sh /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/MainScript1.vs2.sh 








##20150515
##Email chain involving multiple people from UMich lipids consortium and Matthew regarding which Sfrom Kevin on 20140324 explaining the breakdown of "array gene groupings 03 31 2014.xlsx" -- logged here for consistency
##Logged here for consistency

~~~
Details in 'annotated' file for lipid project
Inbox
x  
UChicago
x 

Michael Turchin <mturchin20@gmail.com>
6/18/14

to Matthew 
Hey Matthew,

So it looks like you were able to get a column of information from whomever supplied you with that 'annotated' file that included whether the SNP of interest had been:
1) Found as significant in the 2010 paper
2) Not found as significant in the 2010 paper
3) Not associated with a gene at all

You use this column to grab the SNPs that had been found as significant to create priors. Right now the information that is publicly available for the 2013 paper does not include this distinction. There are supplementary files on the Nature Genetics webpage that include these lists, but they are in .pdf format, not something like an excel file. 

Could we possibly try contacting the people at UMich to see if we could get some sort of excel sheet with the 'top SNP hits' or maybe more specifically an 'annotated' file like we received last time? The former would probably suffice though, I think.

Or any other thoughts on this matter?

Thanks Matthew!


Matthew Stephens <stephens999@gmail.com>
6/18/14

to xwen, me 
I think William ccd may be able to help you. To clarify, the information you are asking for is published just not in convenient form right? Ms



Michael Turchin <mturchin20@gmail.com>
6/18/14

to Matthew, xwen 
Right, and hi there William. I'm just interested in getting the list of SNP rsIDs that were considered significant according to the 2013 paper, if that is easily accessible. I think these are the SNPs in supplementary tables 2 and 3 but they are in .pdf format. 


Xiaoquan Wen <xwen@umich.edu>
6/18/14

to Sebanti, Matthew, me 
Hi Sebanti,

See the message below. Do you have the list of SNPs that Michael is interested in?

Thanks,
William 


Sebanti Sengupta <sebanti@umich.edu>
Attachments6/18/14

to Xiaoquan, Matthew, me 
Hi Michael,

Do you need just a list of rsid's?

The known_lipid_loci file attached lists the rsid's for the loci listed in the 2010 paper and the novel_lipid_loci file lists the rsid's found to be new in the 2013 paper. They're not separated by trait, but if you need that I can try to find those lists.

I hope this helps! Let me know if you need further details.

Best,
Sebanti

2 Attachments 
 
 Preview attachment known_lipid_loci.txt

 Text
 known_lipid_loci.txt
 Preview attachment novel_lipid_loci.txt

 Text
 novel_lipid_loci.txt

 Michael Turchin <mturchin20@gmail.com>
 6/19/14

 to Sebanti, Xiaoquan, Matthew 
 Hi Sebanti,

 Great, these look like what we were hoping for, just the rsid's. Thanks so much! If the information regarding them separated by trait is readily available, we'd be happy to get those too, but no rush.

 Thanks guys, also for your quick responses. Best,
 ~Michael


 Sebanti Sengupta <sebanti@umich.edu>
 Attachments6/21/14

 to me, Xiaoquan, Matthew 
 Hi Michael,

 I've attached the list of loci with the associated traits listed (the first trait is the primary trait). Known means listed in the 2010 paper and novel means newly discovered in the 2013 paper.

 Let me know if you have further questions.

 Best,
 Sebanti

 Attachments area
 Preview attachment lipid_loci.xlsx
 Excel
 lipid_loci.xlsx
 .
 .
 .
 ~

##20150519 
##Using http://www.free-ocr.com/ on 'Extended Data Table 2' from http://www.nature.com/nature/journal/v518/n7538/full/nature14177.html and below is the original output (some of which I manually edited/corrected, such as the Chr:Position identifications, in /mnt/lustre/home/mturchin20/Data/GIANT/2014_5/Locke2015.SupplTables1_2.wManualOCREdits.txt) 
~~~
SNP Chr: Position âNotable gene(s) Alleles EAF B SE P value
rs6567160 18:55,980,115 MC4R(B,N) C/'l' 0.236 0.056 0.004 3.93E-53
rs13021737 2:622,348 TMEM18(N) G/A 0.828 0.06 0.004 1.11E-50
rs10938397 4:44,877,284 GNPDA2(N); GABRG1(B) G/A 0.434 0.04 0.003 3.21 E-38
rs543874 1 :176,156,103 SEC16B(N) G/A 0.193 0.048 0.004 2.62E-35
rs2207139 6:50,953,449 TFAP2B(B,N) G/A 0.177 0.045 0.004 4.13E-29
rs11030104 11 :27,641 ,093 BDNF(B,M,N) AIG 0.792 0.041 0.004 5.56E-28
rs3101336 1:72,523,773 NEGR1(B,C,D,N) C/T 0.613 0.033 0.003 2.66E-26
rs7138803 12:48,533,735 BCDIN3D(N); FAIM2(D) A/G 0.384 0.032 0.003 8.15E-24
rs10182181 2:25,003,800 ADCY3(Bâ%Ã©q(â)3)1;(g;3MC(BâG); G/A 0.462 0.031 0.003 8.78E-24
SH2B1(B,M,Q); APOBR(M.Q);
rs3888190 16:28,796,987 ATXN2L(Q); SBK1(Q,D); SULT1A2(Q); AIC 0.403 0.031 0.003 3.14E-23
TUFM(Q)
rs1516725 3:187,306,698 E TV5(N) CIT 0.872 0.045 0.005 1 .89E-22
rs12446632 16: 1 9,842,890 GPRC5B(C,N); IQCK(Q) G/A 0.865 0.04 0.005 1 .48E-18
rs2287019 19:50,894,012 QPCTL(N); GIPR(B,M) CIT 0.804 0.036 0.004 4.59E-18
rs16951275 15:65,864,222 MAP2K5(B,D,N); LBXCOR1 (M) TIC 0.784 0.031 0.004 1.91 E-1 7
rs3817334 11:47.60"/,ss9 MTC"'2(M'Q); g;âBQy;â;(Q"); SPMQ)â T/C 0.401 0.020 0.003 5.155-11
rs2112347 5:75,050,998 POC5(M); HMGCR(B); COL4A3BP(B) TIG 0.629 0.026 0.003 6.19E-17
rs12566985 1:74,774,781 FPGT-TNNI3K(N) G/A 0.446 0.024 0.003 3.28E-15
rs3810291 19:52,260,843 ZC3H4(D,N,Q) A/G 0.666 0.028 0.004 4.81 E-15
rs7141420 14:78,969,207 NRXN3(D,N) T/C 0.527 0.024 0.003 1 .23E-14
rs13078960 3:85,890,280 CADM2(D,N) GIT 0.196 0.03 0.004 1 .74E-14
rs10968576 9:28,404,339 LlNGO2(D,N) G/A 0.32 0.025 0.003 6.61 E-14
rs17024393 1:109,956,211 GNAT2(N); AMPD2(D) C/'l' 0.04 0.066 0.009 7.03E-14
rs12429545 13:53,000,207 OLFM4(B,N) AIG 0.133 0.033 0.005 1 .09E-12
rs13107325 4:103,407,732 SLC39A8(M,N,Q) T/C 0.072 0.048 0.007 1 .83E-12
rs11165643 1:96,696,685 PTBP2(D,N) T/C 0.583 0.022 0.003 2.07E-12
rs17405819 8:76,969,139 HNF4G(B,N) T/C 0.7 0.022 0.003 2.07E-11
rs1016287 2:59,159,129 LlNC01122(N) TIC 0.287 0.023 0.003 2.25E-11
rs4256980 11 :8,630,51 5 TRIM66(D,M, N); TUB(B) GIC 0.646 0.021 0.003 2.90E-11
rs12401738 1178,21 9,349 FUBP1 (N); USP33(D) A/G 0.352 0.021 0.003 1 .15E-10
rs205262 6:34,671,142 C6orf106(N); SNRPC(Q) G/A 0.273 0.022 0.004 1 .75E-10
rs12016871 13:26,915,782 MTIF3(N); GTF3A(Q) T/C 0.203 0.03 0.005 2.29E-10
rs12940622 17:76,230,166 RPTOR(B,N) G/A 0.575 0.018 0.003 2.49E-09
rs11847697 14:29,584,863 PRKD1 (N) TIC 0.042 0.049 0.008 3.99E-09
rs2075650 19:50,087,459 TOMM40(B,N); APOE(B); APOC1 (B) AIG 0.848 0.026 0.005 1 .25E-08
rs2121279 2:142,759,755 LRP1B(N) T/C 0.152 0.025 0.004 2.31 E-08
rs29941 19:39,001,372 KCTD15(N) G/A 0.669 0.018 0.003 2.41 E-08
rs1808579 18:19,358,886 NPC1(B,G,M,Q); C18orf8(N,Q) CI'I' 0.534 0.017 0.003 4.17E-08
~~~



##20150618
##Email chain with how I received the updated IBDGC information from Stephan Ripke & co. from Broad Institute
##Logged here for consistency

~~~

Michael Turchin <mturchin20@uchicago.edu>
Jun 5 (13 days ago)

to bulik 
Hi there Brendan,

I'm a graduate student at the University of Chicago and I've been trying to track down some publicly available GWAS datasets that contain summary data. I noticed in your manuscript "An Atlas of Genetic Correlations across Human Diseases and Traits" that you were maybe able to get a more complete version of the IBDGC dataset beyond what is available at http://www.ibdgenetics.org/downloads.html (such as possibly the summary data for the full GWAS chip)? If this is the case I was just wondering if you had a contact person you'd recommend speaking to to get access to this summary data? Since the original publication associated with the webpage is from 2012 I'm unsure if the contact information is still up-to-date.

Thanks for any direction you can provide. Best,
~Michael  

Brendan Bulik-Sullivan bulik@broadinstitute.org via uchicago.edu 
Jun 5 (13 days ago)

to Mark, Stephan, Michael 
Hi Michael,

The right people to contact are Mark Daly and Stephan Ripke (cc'ed). 

Best,
Brendan


Michael Turchin <mturchin20@uchicago.edu>
Jun 5 (13 days ago)

to Brendan 
Great, thanks Brendan -- I'll follow up with them.

Best,
~Michael


Michael Turchin <mturchin20@gmail.com>
Jun 5 (13 days ago)

to Mark, Stephan 
Hi Mark, Stephan,

As my e-mail to Brendan below points out, I'm a graduate student at the University of Chicago interested in the various publicly available GWAS summary statistic datasets that exist out there at the moment. I noticed in the "An Atlas of Genetic Correlations across Human Diseases and Traits" manuscript that it was mentioned there is potentially a newer version of the IDBGC data than is currently available from http://www.ibdgenetics.org/downloads.html. I'm just curious about the status/availability of this more current dataset and whether it potentially contains the full SNP-chip/imputation list of UC/Crohn's/IBD p-values (and not just what I believe are the SNPs that went on for replication, which is what the webpage contains).

Thanks for any direction or insight either of you can provide. Best,
~Michael


Michael Turchin <mturchin20@uchicago.edu>
Jun 15 (3 days ago)

to Mark, Stephan 
Hi Mark, Stephan,

As my e-mail to Brendan below points out, I'm a graduate student at the University of Chicago interested in the various publicly available GWAS summary statistic datasets that exist out there at the moment. I noticed in the "An Atlas of Genetic Correlations across Human Diseases and Traits" manuscript that it was mentioned there is potentially a newer version of the IDBGC data than is currently available from http://www.ibdgenetics.org/downloads.html. I'm just curious about the status/availability of this more current dataset and whether it potentially contains the full SNP-chip/imputation list of UC/Crohn's/IBD p-values (and not just what I believe are the SNPs that went on for replication, which is what the webpage contains).

Thanks for any direction or insight either of you can provide. Best,
~Michael


Mark Daly mjdaly@atgu.mgh.harvard.edu via uchicago.edu 
Jun 15 (3 days ago)

to mturchin20 
Hi - I'm traveling through July 8 and may not respond promptly to all email - please contact Beth Raynard (luise@atgu.mgh.harvard.edu) or Jill Harris (harris@atgu.mgh.harvard.edu) if an urgent response is needed.

Stephan Ripke sripke@broadinstitute.org via uchicago.edu 
Jun 15 (3 days ago)

to Jeffrey, Michael, Mark 
sorry for not responding before.
jeff, mark, is it ok to share whole genome summary stats from IBD, UC, CD?
if so, can we share the 1KG imputed ones?


Jeffrey Barrett jb26@sanger.ac.uk via uchicago.edu 
Jun 15 (3 days ago)

to Stephan, Michael, Mark 
Yes, we should be sharing them. Iâm also working on getting the latest versions posted to the website.

Cheers,
Jeff
-- The Wellcome Trust Sanger Institute is operated by Genome Research Limited, a charity registered in England with number 1021457 and a company registered in England with number 2742969, whose registered office is 215 Euston Road, London, NW1 2BE.

Michael Turchin <mturchin20@uchicago.edu>
Jun 15 (3 days ago)

to Jeffrey, Stephan, Mark 
Great, thanks Stephan and Jeff for getting back to me as well as sharing the data! Sorry for the follow-up e-mail Stephan, I just assumed the first one got lost.

Best,
~Michael


Mark Daly mjdaly@atgu.mgh.harvard.edu via uchicago.edu 
Jun 15 (3 days ago)

to Stephan, Michael, Jeffrey 
These should have posted long ago - as Jeff says these should just be put up there on the IIBDGC website
so fine to share


The information in this e-mail is intended only for the person to whom it is
addressed. If you believe this e-mail was sent to you in error and the e-mail
contains patient information, please contact the Partners Compliance HelpLine at
http://www.partners.org/complianceline . If the e-mail was sent to you in error
but does not contain patient information, please contact the sender and properly
dispose of the e-mail.


Stephan Ripke sripke@broadinstitute.org via uchicago.edu 
2:54 PM (2 hours ago)

to Mark, Michael, Jeffrey 
here are the three whole genome results
- crohn's 
- UC
- IBD (taking both case sets together against all controls)

these results will differ slightly from the published ones, since it's 1KG imputation. also it contains many more SNPs.
the individuals and their covariates are identical to the analysis from the publication.

http://www.broadinstitute.org/~sripke/share_links/ALHySN0Z6i4O6PZeQgpqXN9Ewtqxit_IBD15_1KG.sh2
~~~


##20150625
##Email chain with how I received the Hamegen RBC data from EGA via Matthew's help
##Logged here for consistency

~~~
Fwd: [EGA Helpdesk #472863] access to study data
Inbox
x  
UChicago
x 

Matthew Stephens
11:00 AM (2 hours ago)

to me 
Sorry it took me a little while to get this to you. Can you see if you can access the data?
Matthew

---------- Forwarded message ----------
From: "Jeff Almeida-King via RT" <ega-helpdesk@ebi.ac.uk>
Date: Jun 10, 2015 10:40 AM
Subject: [EGA Helpdesk #472863] access to study data
To: <stephens999@gmail.com>
Cc: 

Dear Matthew,

An EGA account has now been created with access permissions set for downloading
the public access datasets stored at the EGA and referenced here:
https://www.ebi.ac.uk/ega/dacs/EGAC00000000005

You can download all EGA public datasets from our FTP server detailed below:
Username: ega-box-233
Password: 88fHetwh

For details of how to access this download account please go here:
https://www.ebi.ac.uk/ega/about/ftp-aspera

Kind regards,
Jeff

On Fri Jun 05 23:42:23 2015, stephens999@gmail.com wrote:
> I am interested in getting access to
>
> DATASET: EGAD00010000300
>
> the webpage https://www.ebi.ac.uk/ega/datasets/EGAD00010000300 says:
>
> ***********
> Who controls access to this dataset
>
> For each dataset that requires access control, there is a
> corresponding Data Access Committee (DAC) who determine access
> permissions. Data access requests are reviewed by the relevant DAC,
> not by the EGA. If you need to request access to this data set, please
> contact:
>
> Access to Public Datasets in the EGA
>
> Email: ega-helpdesk@ebi.ac.uk
> ***********
>
> So I'm wondering who to contact.
>
> thanks,
> Matthew
> (mstephens@uchicago.edu)
>



This email is sent from the Hinxton Campus RT tracking system, which is managed for the Sanger Institute and the EBI by the Sanger Institute.


--
 The Wellcome Trust Sanger Institute is operated by Genome Research
  Limited, a charity registered in England with number 1021457 and a
   company registered in England with number 2742969, whose registered
    office is 215 Euston Road, London, NW1 2BE.

    Michael Turchin <mturchin20@gmail.com>
    11:21 AM (2 hours ago)

    to Matthew 
    Hey Matthew,

    Sounds good and no problem, I think I was going to follow-up about this sometime next week with you anyways. I'll try to take a look sometime today and let you know how it goes.

    Thanks for taking a look into this! 
~~~


