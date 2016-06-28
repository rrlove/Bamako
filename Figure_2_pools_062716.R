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
list(geom_text(data=inversions[inversions$name=="2La",], aes(x=(start + end)/2, y=0.035, label="a"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rj",], aes(x=(start + end)/2, y=0.035, label="j"), fontface="bold", color="black", size=7), 
geom_text(data=inversions[inversions$name=="2Rb",], aes(x=(start + end)/2, y=0.035, label="b"), fontface="bold", color="black", size=7), 
geom_text(data=inversions[inversions$name=="2Rc",], aes(x=(start + end)/2, y=0.035, label="c"), fontface="bold", color="black", size=7), 
geom_text(data=inversions[inversions$name=="2Ru",], aes(x=(start + end)/2, y=0.035, label="u"), fontface="bold", color="black", size=7))

#read in D, pi, and FST; D downsampled (for AD006 and AD012), pi NOT
coluzzii.D.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/091015_results/AD006.downsampled.200K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.coluzzii"))
Bamako.D.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/091015_results/AD012.downsampled.200K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.Bamako"))
gambiae.D.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/091015_results/AD019.200K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.gambiae"))

D.200 <- Reduce(function(...) merge(..., by=c("chrom","pos"), all=TRUE), list(coluzzii.D.200,Bamako.D.200,gambiae.D.200))

coluzzii.pi.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/091015_results/AD006.200K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.coluzzii"))
Bamako.pi.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/091015_results/AD012.200K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.Bamako"))
gambiae.pi.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/091015_results/AD019.200K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.gambiae"))

pi.200 <- Reduce(function(...) merge(..., by=c("chrom","pos"), all=TRUE), list(coluzzii.pi.200,Bamako.pi.200,gambiae.pi.200))

fst.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/091015_results/SMBam.200K.windowed.fst.cleaned", header=FALSE, na.strings="na", col.names=c("chrom", "pos", "nSNPs", "frac.cov", "ave.min.cov", "ColBam", "ColGam", "BamGam"))

D.pi.fst.200 <- Reduce(function(...) merge(..., by=c("chrom","pos"), all=TRUE), list(D.200,pi.200,fst.200))

all.200 <- D.pi.fst.200[D.pi.fst.200$chrom != "Y_unplaced" & D.pi.fst.200$chrom != "UNKN" & D.pi.fst.200$chrom != "Mt",c("chrom","pos","D.coluzzii","D.Bamako","D.gambiae","pi.coluzzii","pi.Bamako","pi.gambiae","ColBam","ColGam","BamGam")]

##read in DXY

ColBam.dxy.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/Dxy_091715/new_consensi_061316/Dxy_new_consensi_061016/coluzzii-Bamako-quad-061316_parsed.dxy.200K.txt", header=FALSE, na.strings="na",col.names=c("chrom","pos","dColBam","ColQuad_CB","BamQuad_CB","coverage"))

ColGam.dxy.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/Dxy_091715/new_consensi_061316/Dxy_new_consensi_061016/coluzzii-gambiae-quad-061316_parsed.dxy.200K.txt", header=FALSE, na.strings="na",col.names=c("chrom","pos","dColGam","ColQuad_CG","GamQuad_CG","coverage"))

BamGam.dxy.200 <- read.table("/media/anophele/HardDrive2/NGSdata/pools/Dxy_091715/new_consensi_061316/Dxy_new_consensi_061016/Bamako-gambiae-quad-061316_parsed.dxy.200K.txt", header=FALSE, na.strings="na",col.names=c("chrom","pos","dBamGam","BamQuad_BG","GamQuad_BG","coverage"))

dxy.200 <- Reduce(function(...) merge(..., by=c("chrom","pos"), all=TRUE), list(ColBam.dxy.200,ColGam.dxy.200,BamGam.dxy.200))

dxy.200 <-dxy.200[dxy.200$chrom != "Mt" & dxy.200$chrom != "UNKN" & dxy.200$chrom != "Y_unplaced",c("chrom","pos","dColBam","dColGam","dBamGam")]

##merge both data frames
all_200_dxy <- merge(dxy.200,all.200, by=c("chrom","pos"),all=TRUE)

##reorder the chromosomes
all_200_dxy$chrom <- factor(all_200_dxy$chrom, levels=c("2R","2L","3R","3L","X"))

##read in the outliers

fst_outliers <- read.table("fst_outliers.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

sweep_outliers <- read.table("D_thetaless_outliers.bed", header=FALSE, col.names=c("chrom","start","end"), stringsAsFactors=FALSE)

##plotting

pi.plot <- ggplot() +  geom_line(data=all_200_dxy, aes(x=pos, y=pi.coluzzii, color="coluzzii")) + geom_line(data=all_200_dxy, aes(x=pos, y=pi.gambiae, color="gambiae")) + geom_line(data=all_200_dxy, aes(x=pos, y=pi.Bamako, color="Bamako")) + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + xlab("") + ylab("pi") + ggtitle("") + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), legend.justification=c(0,1), legend.position=c(0.805,1), legend.title=element_text(size=18), legend.text=element_text(size=15), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), plot.margin=unit(c(-1,.3,-.95,0),"lines"), strip.text.x = element_text(size=18, face="bold"), strip.background=element_rect(fill="white"), legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.9),"lines")) + inversions.func() + inversion.labels()

TajimasD.plot <- ggplot() + geom_line(data=all_200_dxy, aes(x=pos, y=D.coluzzii, color="coluzzii")) + geom_line(data=all_200_dxy, aes(x=pos, y=D.gambiae, color="gambiae")) +  geom_line(data=all_200_dxy, aes(x=pos, y=D.Bamako, color="Bamako")) + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + xlab("") + ylab("Tajima's D") + theme(legend.position="none", axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,-1,0),"lines")) + inversions.func()

dxy.plot <- ggplot() + geom_line(data=all_200_dxy, aes(x=pos, y=dColGam, color="coluzzii-gambiae")) + geom_line(data=all_200_dxy, aes(x=pos, y=dColBam, color="coluzzii-Bamako")) + geom_line(data=all_200_dxy, aes(x=pos, y=dBamGam, color="Bamako-gambiae")) + xlab("Mb") + ylab("dXY") + scale_colour_manual(values=rev(cbbPalette.contrast), name="Comparison") + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + theme(legend.position=c(0.805,1), legend.justification=c(0,1), legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.9),"lines"), legend.title=element_text(size=18), legend.text=element_text(size=15), axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,-1,0),"lines")) + scale_x_continuous(labels=function(x)x/1000000) + inversions.func() 

fst.all.plot <- ggplot() + geom_line(data=all_200_dxy, aes(x=pos, y=ColGam, color="coluzzii-gambiae")) + geom_line(data=all_200_dxy, aes(x=pos, y=ColBam, color="coluzzii-Bamako")) + geom_line(data=all_200_dxy, aes(x=pos, y=BamGam, color="Bamako-gambiae")) + xlab("Mb") + ylab("FST") + scale_colour_manual(values=rev(cbbPalette.contrast), name="Comparison") + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + theme(legend.position="none", axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), axis.title.x=element_text(size=18), axis.text.x=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,0,0),"lines")) + scale_x_continuous(labels=function(x)x/1000000) + inversions.func() + geom_rect(data=fst_outliers, aes(xmin=start, ymin=-.03, xmax=end, ymax=-.04), color="black") + geom_rect(data=sweep_outliers, aes(xmin=start, ymin=-.055, xmax=end, ymax=-.065), color="black")

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

##postscript polishing in Illustrator:
##
##turn "pi" to pi symbol, italicize, move away from axis
##italicize D in Tajima's D
##italicize FST
##italicize Dxy
##italicize inversions
##italicize gambiae and coluzzii in legends
##


