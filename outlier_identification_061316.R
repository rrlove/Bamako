##this script identifies two sets of windows:
##
##1.) 1K windows in the top 1% of all Bamako-gambiae FST values in each genomic partition
##2.) 1K windows in the top 1% of all Bamako-gambiae FST values and the bottom 1% of all Bamako pi values in each genomic partition
##
##

##read in and merge data
library(dplyr)
setwd("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results")

pi <- read.table("AD012.1K.nonoverlapping.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.Bamako"))

fst <- read.table("SMBam.mincount.nooverlap.1K.082015.fst.cleaned", header=FALSE, na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "ave.min.cov", "ColBam", "ColGam", "BamGam"))

pi.fst <- merge(pi,fst, by=c("chrom","pos"),all=TRUE,sort=FALSE)

pi.fst <- pi.fst[pi.fst$chrom != "Y_unplaced" & pi.fst$chrom != "UNKN" & pi.fst$chrom != "Mt",c("chrom","pos","pi.Bamako","ColBam","ColGam","BamGam")]

pi.fst$chrom <- factor(pi.fst$chrom)

##set genomic partitions and pull out regions of interest

euchromatin_coords <- data.frame(chrom=c("X","2R","2L","2L","3R","3R","3L","3L"),start=c(1,1,2487770,6015228,1,41888356,1896830,5133257),end=c(19928574,58969802,5042389,49364325,38815826,52131026,4235209,41963435))

##levels(pi.fst$chrom) == levels(euchromatin_coords$chrom)
##[1] TRUE TRUE TRUE TRUE TRUE
##note to self, this does test order as well-- it's not a set operation

inversions <- data.frame(
name=c("2La","2Rj","2Rb","2Rc","2Ru"),chrom=c("2L",rep("2R",4)),start=c(20524058,3262186,19023925,26780000,31480000),end=c(42165532,15750717,26758676,31450000,35500000))

##euchromatin_coords
##  chrom    start      end
##1     X        1 19928574
##2    2R        1 58969802
##3    2L  2487770  5042389
##4    2L  6015228 49364325
##5    3R        1 38815826
##6    3R 41888356 52131026
##7    3L  1896830  4235209
##8    3L  5133257 41963435
##

##inversion_coords
##  name chrom    start      end
##1  2La    2L 20524058 42165532
##2  2Rj    2R  3262186 15750717
##3  2Rb    2R 19023925 26758676
##4  2Rc    2R 26780000 31450000
##5  2Ru    2R 31480000 35500000

euchromatin_1 <- pi.fst[pi.fst$chrom == euchromatin_coords[1,"chrom"] & pi.fst$pos > euchromatin_coords[1,"start"] & pi.fst$pos < euchromatin_coords[1,"end"],]

euchromatin_2 <- pi.fst[pi.fst$chrom == euchromatin_coords[2,"chrom"] & pi.fst$pos > euchromatin_coords[2,"start"] & pi.fst$pos < euchromatin_coords[2,"end"],]

euchromatin_3 <- pi.fst[pi.fst$chrom == euchromatin_coords[3,"chrom"] & pi.fst$pos > euchromatin_coords[3,"start"] & pi.fst$pos < euchromatin_coords[3,"end"],]

euchromatin_4 <- pi.fst[pi.fst$chrom == euchromatin_coords[4,"chrom"] & pi.fst$pos > euchromatin_coords[4,"start"] & pi.fst$pos < euchromatin_coords[4,"end"],]

euchromatin_5 <- pi.fst[pi.fst$chrom == euchromatin_coords[5,"chrom"] & pi.fst$pos > euchromatin_coords[5,"start"] & pi.fst$pos < euchromatin_coords[5,"end"],]

euchromatin_6 <- pi.fst[pi.fst$chrom == euchromatin_coords[6,"chrom"] & pi.fst$pos > euchromatin_coords[6,"start"] & pi.fst$pos < euchromatin_coords[6,"end"],]

euchromatin_7 <- pi.fst[pi.fst$chrom == euchromatin_coords[7,"chrom"] & pi.fst$pos > euchromatin_coords[7,"start"] & pi.fst$pos < euchromatin_coords[7,"end"],]

euchromatin_8 <- pi.fst[pi.fst$chrom == euchromatin_coords[8,"chrom"] & pi.fst$pos > euchromatin_coords[8,"start"] & pi.fst$pos < euchromatin_coords[8,"end"],]

##summary(euchromatin_1);summary(euchromatin_2);summary(euchromatin_3);summary(euchromatin_4);summary(euchromatin_5);summary(euchromatin_6);summary(euchromatin_7);summary(euchromatin_8)

euchromatin <- rbind(euchromatin_1,euchromatin_2,euchromatin_3,euchromatin_4,euchromatin_5,euchromatin_6,euchromatin_7,euchromatin_8)

twoLa <- euchromatin[euchromatin$chrom=="2L" & euchromatin$pos > inversions[inversions$name=="2La","start"] & euchromatin$pos < inversions[inversions$name=="2La","end"],]
##nrow(twoLa)
##[1] 21642

twoRj <- euchromatin[euchromatin$chrom=="2R" & euchromatin$pos > inversions[inversions$name=="2Rj","start"] & euchromatin$pos < inversions[inversions$name=="2Rj","end"],]
##nrow(twoRj)
##[1] 12489

twoRb <- euchromatin[euchromatin$chrom=="2R" & euchromatin$pos > inversions[inversions$name=="2Rb","start"] & euchromatin$pos < inversions[inversions$name=="2Rb","end"],]
##nrow(twoRb)
##[1] 7735

twoRc <- euchromatin[euchromatin$chrom=="2R" & euchromatin$pos > inversions[inversions$name=="2Rc","start"] & euchromatin$pos < inversions[inversions$name=="2Rc","end"],]
##nrow(twoRc)
##[1] 4670

twoRu <- euchromatin[euchromatin$chrom=="2R" & euchromatin$pos > inversions[inversions$name=="2Ru","start"] & euchromatin$pos < inversions[inversions$name=="2Ru","end"],]
##nrow(twoRu)
##[1] 4020

inverted <- rbind(twoLa, twoRj, twoRb, twoRc, twoRu)
##nrow(inverted)
##[1] 50556

colinear <- anti_join(euchromatin, inverted, by=c("chrom","pos"))
##nrow(colinear)
##[1] 162471

##nrow(inverted) + nrow(colinear) == nrow(euchromatin)
##[1] TRUE

colinear_autosomal <- colinear[!colinear$chrom=="X",]
##nrow(colinear_autosomal)
##[1] 142542

##identify windows with exceptional FST values
twoLa_fst_outliers <- twoLa[!is.na(twoLa$BamGam) & twoLa$BamGam > quantile(twoLa$BamGam, 0.99, na.rm=TRUE),]
##nrow(twoLa_fst_outliers)
##[1] 148

twoRj_fst_outliers <- twoRj[!is.na(twoRj$BamGam) & twoRj$BamGam > quantile(twoRj$BamGam, 0.99, na.rm=TRUE),]
##nrow(twoRj_fst_outliers)
##[1] 110

twoRb_fst_outliers <- twoRb[!is.na(twoRb$BamGam) & twoRb$BamGam > quantile(twoRb$BamGam, 0.99, na.rm=TRUE),]
##nrow(twoRb_fst_outliers)
##[1] 71

twoRc_fst_outliers <- twoRc[!is.na(twoRc$BamGam) & twoRc$BamGam > quantile(twoRc$BamGam, 0.99, na.rm=TRUE),]
##nrow(twoRc_fst_outliers)
##[1] 43

twoRu_fst_outliers <- twoRu[!is.na(twoRu$BamGam) & twoRu$BamGam > quantile(twoRu$BamGam, 0.99, na.rm=TRUE),]
##nrow(twoRu_fst_outliers)
##[1] 30

inverted_fst_outliers <- inverted[!is.na(inverted$BamGam) & inverted$BamGam > quantile(inverted$BamGam, 0.99, na.rm=TRUE),]
##nrow(inverted_fst_outliers)
##[1] 399

colinear_fst_outliers <- colinear[!is.na(colinear$BamGam) & colinear$BamGam > quantile(colinear$BamGam, 0.99, na.rm=TRUE),]
##nrow(colinear_fst_outliers)
##[1] 1347

colinear_autosomal_fst_outliers <- colinear_autosomal[!is.na(colinear_autosomal$BamGam) & colinear_autosomal$BamGam > quantile(colinear_autosomal$BamGam, 0.99, na.rm=TRUE),]
##nrow(colinear_autosomal_fst_outliers)
##[1] 1172

ecks <- euchromatin[euchromatin$chrom=="X",]
##nrow(ecks)
##[1] 19929

X_fst_outliers <- ecks[ !is.na(ecks$BamGam) & ecks$BamGam > quantile(ecks$BamGam,0.99,na.rm=TRUE) ,]
##nrow(X_fst_outliers)
##[1] 175

##identify windows with exceptional FST and pi values

twoLa_sweep_outliers <- twoLa[!is.na(twoLa$BamGam) & !is.na(twoLa$pi.Bamako) & twoLa$BamGam > quantile(twoLa$BamGam, 0.99, na.rm=TRUE) & twoLa$pi.Bamako < quantile(twoLa$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(twoLa_sweep_outliers)
##[1] 3
##nrow(intersect(twoLa_sweep_outliers,twoLa_fst_outliers))
##[1] 3
##nrow(intersect(twoLa_sweep_outliers,twoLa_fst_outliers)) == nrow(twoLa_sweep_outliers)
##[1] TRUE

twoRj_sweep_outliers <- twoRj[!is.na(twoRj$BamGam) & !is.na(twoRj$pi.Bamako) & twoRj$BamGam > quantile(twoRj$BamGam, 0.99, na.rm=TRUE) & twoRj$pi.Bamako < quantile(twoRj$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(twoRj_sweep_outliers)
##[1] 52
##nrow(intersect(twoRj_sweep_outliers,twoRj_fst_outliers))
##[1] 52
##nrow(intersect(twoRj_sweep_outliers,twoRj_fst_outliers)) == nrow(twoRj_sweep_outliers)
##[1] TRUE

twoRb_sweep_outliers <- twoRb[!is.na(twoRb$BamGam) & !is.na(twoRb$pi.Bamako) & twoRb$BamGam > quantile(twoRb$BamGam, 0.99, na.rm=TRUE) & twoRb$pi.Bamako < quantile(twoRb$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(twoRb_sweep_outliers)
##[1] 12
##nrow(intersect(twoRb_sweep_outliers,twoRb_fst_outliers))
##[1] 12
##nrow(intersect(twoRb_sweep_outliers,twoRb_fst_outliers)) == nrow(twoRb_sweep_outliers)
##[1] TRUE

twoRc_sweep_outliers <- twoRc[!is.na(twoRc$BamGam) & !is.na(twoRc$pi.Bamako) & twoRc$BamGam > quantile(twoRc$BamGam, 0.99, na.rm=TRUE) & twoRc$pi.Bamako < quantile(twoRc$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(twoRc_sweep_outliers)
##[1] 6
##nrow(intersect(twoRc_sweep_outliers,twoRc_fst_outliers))
##[1] 6
##nrow(intersect(twoRc_sweep_outliers,twoRc_fst_outliers)) == nrow(twoRc_sweep_outliers)
##[1] TRUE

twoRu_sweep_outliers <- twoRu[!is.na(twoRu$BamGam) & !is.na(twoRu$pi.Bamako) & twoRu$BamGam > quantile(twoRu$BamGam, 0.99, na.rm=TRUE) & twoRu$pi.Bamako < quantile(twoRu$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(twoRu_sweep_outliers)
##[1] 3
##nrow(intersect(twoRu_sweep_outliers,twoRu_fst_outliers))
##[1] 3
##nrow(intersect(twoRu_sweep_outliers,twoRu_fst_outliers)) == nrow(twoRu_sweep_outliers)
##[1] TRUE

inverted_sweep_outliers <- inverted[!is.na(inverted$BamGam) & !is.na(inverted$pi.Bamako) & inverted$BamGam > quantile(inverted$BamGam, 0.99, na.rm=TRUE) & inverted$pi.Bamako < quantile(inverted$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(inverted_sweep_outliers)
##[1] 147
##nrow(intersect(inverted_sweep_outliers,inverted_fst_outliers))
##[1] 147
##nrow(intersect(inverted_sweep_outliers,inverted_fst_outliers)) == nrow(inverted_sweep_outliers)
##[1] TRUE

colinear_sweep_outliers <- colinear[!is.na(colinear$BamGam) & !is.na(colinear$pi.Bamako) & colinear$BamGam > quantile(colinear$BamGam, 0.99, na.rm=TRUE) & colinear$pi.Bamako < quantile(colinear$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(colinear_sweep_outliers)
##[1] 22
##nrow(intersect(colinear_sweep_outliers,colinear_fst_outliers))
##[1] 22
##nrow(intersect(colinear_sweep_outliers,colinear_fst_outliers)) == nrow(colinear_sweep_outliers)
##[1] TRUE

colinear_autosomal_sweep_outliers <- colinear_autosomal[!is.na(colinear_autosomal$BamGam) & !is.na(colinear_autosomal$pi.Bamako) & colinear_autosomal$BamGam > quantile(colinear_autosomal$BamGam, 0.99, na.rm=TRUE) & colinear_autosomal$pi.Bamako < quantile(colinear_autosomal$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(colinear_autosomal_sweep_outliers)
##[1] 19
##nrow(intersect(colinear_autosomal_sweep_outliers,colinear_autosomal_fst_outliers))
##[1] 19
##nrow(intersect(colinear_autosomal_sweep_outliers,colinear_autosomal_fst_outliers)) == nrow(colinear_autosomal_sweep_outliers)
##[1] TRUE

X_sweep_outliers <- ecks[!is.na(ecks$BamGam) & !is.na(ecks$pi.Bamako) & ecks$BamGam > quantile(ecks$BamGam,0.99,na.rm=TRUE) & ecks$pi.Bamako < quantile(ecks$pi.Bamako, 0.01, na.rm=TRUE),]
##nrow(X_sweep_outliers)
##[1] 12
##nrow(intersect(X_sweep_outliers,X_fst_outliers))
##[1] 12
##nrow(intersect(X_sweep_outliers,X_fst_outliers)) == nrow(X_sweep_outliers)
##[1] TRUE

##create tables to write out and write them out

options(scipen=10)

##grep("outliers",ls(),value=TRUE)
## [1] "colinear_autosomal_fst_outliers"   "colinear_autosomal_sweep_outliers" ##"colinear_fst_outliers"            
##[4] "colinear_sweep_outliers"           "inverted_fst_outliers"             "inverted_sweep_outliers"          
##[7] "twoLa_fst_outliers"                "twoLa_sweep_outliers"              "twoRb_fst_outliers"          ##[10] "twoRb_sweep_outliers"              "twoRc_fst_outliers"                ##"twoRc_sweep_outliers"             
##[13] "twoRj_fst_outliers"                "twoRj_sweep_outliers"              "twoRu_fst_outliers"               
##[16] "twoRu_sweep_outliers"              "X_fst_outliers"                    "X_sweep_outliers"

chunk_names <- list("colinear_autosomal_fst_outliers"=colinear_autosomal_fst_outliers,"colinear_autosomal_sweep_outliers"=colinear_autosomal_sweep_outliers,"colinear_fst_outliers"=colinear_fst_outliers,"colinear_sweep_outliers"=colinear_sweep_outliers,"inverted_fst_outliers"=inverted_fst_outliers,"inverted_sweep_outliers"=inverted_sweep_outliers,"twoLa_fst_outliers"=twoLa_fst_outliers,"twoLa_sweep_outliers"=twoLa_sweep_outliers,"twoRb_fst_outliers"=twoRb_fst_outliers,"twoRb_sweep_outliers"=twoRb_sweep_outliers,"twoRc_fst_outliers"=twoRc_fst_outliers,"twoRc_sweep_outliers"=twoRc_sweep_outliers,"twoRj_fst_outliers"=twoRj_fst_outliers,"twoRj_sweep_outliers"=twoRj_sweep_outliers,"twoRu_fst_outliers"=twoRu_fst_outliers,"twoRu_sweep_outliers"=twoRu_sweep_outliers,"X_fst_outliers"=X_fst_outliers,"X_sweep_outliers"=X_sweep_outliers)

returned <- lapply(names(chunk_names), function(x){
chunk <- data.frame(chunk_names[x])    
out_bed_table <- data.frame(chunk[,1],chunk[,2]-500,chunk[,2]+500)
    colnames(out_bed_table) <- c("chrom","start","end")
    out_file_name = paste(x,".",format(Sys.Date(),"%m%d%y"),".bed",sep="")
    write.table(out_bed_table, file=out_file_name, quote=F, col.names=F, row.names=F, sep="\t")
    rm(out_bed_table)
    rm(chunk)
})

