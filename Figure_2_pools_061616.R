setwd("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results")

library(ggplot2)
library(gridExtra)
library(scales)
library(grid)
cbbPalette.contrast <-c("#0072B2", "#D55E00", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")

##set up the inversions

inversions <- data.frame(
name=c("2La","2Rj","2Rb","2Rc","2Ru"),chrom=c("2L",rep("2R",4)),start=c(20524058,3262186,19023925,26780000,31480000),end=c(42165532,15750717,26758676,31450000,35500000))

inversions.func <- function()
list(
geom_rect(data=inversions[inversions$name=="2La",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_rect(data=inversions[inversions$name=="2Rj",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_rect(data=inversions[inversions$name=="2Rb",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3,
fill="#999999"),geom_vline(data=inversions[inversions$name=="2Rc",], aes(xintercept=start), colour="black", linetype="longdash"),
geom_rect(data=inversions[inversions$name=="2Rc",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_vline(data=inversions[inversions$name=="2Rc",], aes(xintercept=end), colour="black", linetype="longdash"),geom_rect(data=inversions[inversions$name=="2Ru",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"))

inversion.labels <- function()
list(geom_text(data=inversions[inversions$name=="2La",], aes(x=(start + end)/2, y=0.04, label="a"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rj",], aes(x=(start + end)/2, y=0.04, label="j"), fontface="bold", color="black", size=7), 
geom_text(data=inversions[inversions$name=="2Rb",], aes(x=(start + end)/2, y=0.04, label="b"), fontface="bold", color="black", size=7), 
geom_text(data=inversions[inversions$name=="2Rc",], aes(x=(start + end)/2, y=0.04, label="c"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Ru",], aes(x=(start + end)/2, y=0.04, label="u"), fontface="bold", color="black", size=7))

#read in D, pi, and FST; D downsampled (for AD006 and AD012), pi NOT
coluzzii.D.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/AD006.downsampled.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.coluzzii"))
Bamako.D.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/AD012.downsampled.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.Bamako"))
gambiae.D.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/AD019.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.gambiae"))

ColBam.D.100 <- merge(coluzzii.D.100, Bamako.D.100, by=c("chrom", "pos"),all=TRUE)

D.100 <- merge(ColBam.D.100, gambiae.D.100, by=c("chrom", "pos"),all=TRUE)

coluzzii.pi.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/AD006.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.coluzzii"))
Bamako.pi.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/AD012.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.Bamako"))
gambiae.pi.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/AD019.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.gambiae"))

ColBam.pi.100 <- merge(coluzzii.pi.100, Bamako.pi.100, by=c("chrom", "pos"),all=TRUE)
pi.100 <-merge(ColBam.pi.100, gambiae.pi.100, by=c("chrom", "pos"),all=TRUE)

D.pi.100 <-merge(D.100,pi.100, by=c("chrom","pos"),all=TRUE,sort=FALSE)

fst.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/SMBam.100K.windowed.fst.cleaned", header=FALSE, na.strings="na", col.names=c("chrom", "pos", "nSNPs", "frac.cov", "ave.min.cov", "ColBam", "ColGam", "BamGam"))

D.pi.fst.100 <- merge(D.pi.100,fst.100, by=c("chrom","pos"),all=TRUE)

all.100 <- D.pi.fst.100[D.pi.fst.100$chrom != "Y_unplaced" & D.pi.fst.100$chrom != "UNKN" & D.pi.fst.100$chrom != "Mt",c("chrom","pos","D.coluzzii","D.Bamako","D.gambiae","pi.coluzzii","pi.Bamako","pi.gambiae","ColBam","ColGam","BamGam")]

##read in DXY

ColBam.dxy.100 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/Dxy_091715/new_consensi_061316/Dxy_new_consensi_061016/coluzzii-Bamako-quad-061316_parsed.dxy.100K.txt", header=FALSE, na.strings="na",col.names=c("chrom","pos","dColBam","ColQuad_CB","BamQuad_CB","coverage"))

ColGam.dxy.100 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/Dxy_091715/new_consensi_061316/Dxy_new_consensi_061016/coluzzii-gambiae-quad-061316_parsed.dxy.100K.txt", header=FALSE, na.strings="na",col.names=c("chrom","pos","dColGam","ColQuad_CG","GamQuad_CG","coverage"))

BamGam.dxy.100 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/Dxy_091715/new_consensi_061316/Dxy_new_consensi_061016/Bamako-gambiae-quad-061316_parsed.dxy.100K.txt", header=FALSE, na.strings="na",col.names=c("chrom","pos","dBamGam","BamQuad_BG","GamQuad_BG","coverage"))

dxy.100_1 <- merge(ColBam.dxy.100,ColGam.dxy.100,by=c("chrom","pos"),all=TRUE)
dxy.100 <- merge(dxy.100_1,BamGam.dxy.100,by=c("chrom","pos"),all=TRUE)

dxy.100 <-dxy.100[dxy.100$chrom != "Mt" & dxy.100$chrom != "UNKN" & dxy.100$chrom != "Y_unplaced",c("chrom","pos","dColBam","dColGam","dBamGam")]

##merge both data frames
all_100_dxy <- merge(dxy.100,all.100, by=c("chrom","pos"),all=TRUE)

##reorder the chromosomes
all_100_dxy$chrom <- factor(all_100_dxy$chrom, levels=c("2R","2L","3R","3L","X"))

##read in the outliers, label, and bind them together

twoLa_fst <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoLa_fst_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRj_fst <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRj_fst_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRb_fst <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRb_fst_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRc_fst <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRc_fst_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRu_fst <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRu_fst_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

colinear_fst <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/colinear_fst_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

fst_outliers <- rbind(twoLa_fst, twoRj_fst, twoRb_fst, twoRc_fst, twoRc_fst, colinear_fst)
fst_outliers$method <- "fst"

twoLa_sweep <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoLa_sweep_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRj_sweep <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRj_sweep_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRb_sweep <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRb_sweep_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRc_sweep <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRc_sweep_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

twoRu_sweep <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/twoRu_sweep_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

colinear_sweep <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/colinear_sweep_outliers.061616.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

sweep_outliers <- rbind(twoLa_sweep, twoRj_sweep, twoRb_sweep, twoRc_sweep, twoRc_sweep, colinear_sweep)
sweep_outliers$method <- "sweep"

outliers <- rbind(sweep_outliers, fst_outliers)

##plotting

pi.plot <- ggplot() +  geom_line(data=all_100_dxy, aes(x=pos, y=pi.coluzzii, color="coluzzii")) + geom_line(data=all_100_dxy, aes(x=pos, y=pi.gambiae, color="gambiae")) + geom_line(data=all_100_dxy, aes(x=pos, y=pi.Bamako, color="Bamako")) + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + xlab("") + ylab("pi") + ggtitle("") + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), legend.justification=c(0,1), legend.position=c(0.805,1), legend.title=element_text(size=18), legend.text=element_text(size=15), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), plot.margin=unit(c(-1,.3,-.95,0),"lines"), strip.text.x = element_text(size=18, face="bold"), strip.background=element_rect(fill="white"), legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.9),"lines")) + inversions.func() + inversion.labels()

TajimasD.plot <- ggplot() + geom_line(data=all_100_dxy, aes(x=pos, y=D.coluzzii, color="coluzzii")) + geom_line(data=all_100_dxy, aes(x=pos, y=D.gambiae, color="gambiae")) +  geom_line(data=all_100_dxy, aes(x=pos, y=D.Bamako, color="Bamako")) + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + xlab("") + ylab("Tajima's D") + theme(legend.position="none", axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,-1,0),"lines")) + inversions.func()

dxy.plot <- ggplot() + geom_line(data=all_100_dxy, aes(x=pos, y=dColGam, color="coluzzii-gambiae")) + geom_line(data=all_100_dxy, aes(x=pos, y=dColBam, color="coluzzii-Bamako")) + geom_line(data=all_100_dxy, aes(x=pos, y=dBamGam, color="Bamako-gambiae")) + xlab("Mb") + ylab("Dxy") + scale_colour_manual(values=rev(cbbPalette.contrast), name="Comparison") + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + theme(legend.position=c(0.805,1), legend.justification=c(0,1), legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.9),"lines"), legend.title=element_text(size=18), legend.text=element_text(size=15), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,-1,0),"lines")) + scale_x_continuous(labels=function(x)x/1000000) + inversions.func() 

fst.all.plot <- ggplot() + geom_line(data=all_100_dxy, aes(x=pos, y=ColGam, color="coluzzii-gambiae")) + geom_line(data=all_100_dxy, aes(x=pos, y=ColBam, color="coluzzii-Bamako")) + geom_line(data=all_100_dxy, aes(x=pos, y=BamGam, color="Bamako-gambiae")) + xlab("Mb") + ylab("FST") + scale_colour_manual(values=rev(cbbPalette.contrast), name="Comparison") + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + theme(legend.position="none", axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), axis.title.x=element_text(size=18), axis.text.x=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,0,0),"lines")) + scale_x_continuous(labels=function(x)x/1000000) + inversions.func() + geom_rect(data=fst_outliers, aes(xmin=start, ymin=-.03, xmax=end, ymax=-.04), color="black") + geom_rect(data=sweep_outliers, aes(xmin=start, ymin=-.055, xmax=end, ymax=-.065), color="black")

wpi <- ggplotGrob(pi.plot)
wD <- ggplotGrob(TajimasD.plot)
wdxy <- ggplotGrob(dxy.plot)
wfst <- ggplotGrob(fst.all.plot)

maxWidth = grid::unit.pmax(wpi$widths[2:5],wD$widths[2:5],wfst$widths[2:5])
wpi$widths[2:5] <- as.list(maxWidth)
wD$widths[2:5] <- as.list(maxWidth)
wdxy$widths[2:5] <- as.list(maxWidth)
wfst$widths[2:5] <- as.list(maxWidth)

pdf(paste("Figure 2.",format(Sys.Date(), "%m%d%y"),".pdf",sep=""), width=22, height=11, useDingbats=FALSE)
grid.arrange(wpi,wD,wdxy,wfst,ncol=1,heights=c(1,1,1,2))
dev.off()

tiff(paste("Figure 2.",format(Sys.Date(), "%m%d%y"),".tiff",sep=""), width=2200, height=1100)
grid.arrange(wpi,wD,wdxy,wfst,ncol=1,heights=c(1,1,1,2))
dev.off()

##postscript polishing in Illustrator:
##
##turn "pi" to pi symbol, italicize, move away from axis
##italicize D in Tajima's D
##italicize FST
##italicize Dxy
##italicize inversions
##italicize gambiae and coluzzii in legends
##


