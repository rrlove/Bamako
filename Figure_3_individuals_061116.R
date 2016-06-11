setwd("/media/anophele/Hopper/Archives/Bamako/Individuals/Bamako_individual_Results/all.Bamakoset.2")

library(ggplot2)
library(scales)
library(gridExtra)
library(grid)

cbbPalette.contrast <-c("#0072B2", "#D55E00", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")

##inversions

inversions <- data.frame(
  name=c("2La","2Rj","2Rb","2Rc","2Ru"),chrom=c("2L",rep("2R",4)),start=c(20524058,3262186,19023925,26780000,31480000),end=c(42165532,15750717,26758676,31450000,35500000))

inversions.func <- function()
  list(
    geom_rect(data=inversions[inversions$name=="2La",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_rect(data=inversions[inversions$name=="2Rj",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_rect(data=inversions[inversions$name=="2Rb",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3,
                                                                                                                                                                                                                                                                            fill="#999999"),geom_vline(data=inversions[inversions$name=="2Rc",], aes(xintercept=start), colour="black", linetype="longdash"),
    geom_rect(data=inversions[inversions$name=="2Rc",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_vline(data=inversions[inversions$name=="2Rc",], aes(xintercept=end), colour="black", linetype="longdash"),geom_rect(data=inversions[inversions$name=="2Ru",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"))

inversion.labels <- function()
  list(geom_text(data=inversions[inversions$name=="2La",], aes(x=(start + end)/2, y=0.0045, label="a"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rj",], aes(x=(start + end)/2, y=0.0045, label="j"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rb",], aes(x=(start + end)/2, y=0.0045, label="b"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rc",], aes(x=(start + end)/2, y=0.0045, label="c"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Ru",], aes(x=(start + end)/2, y=0.0045, label="u"), fontface="bold", color="black", size=7))

##pi
Bamako.pi <- read.table("Bamako.200.windowed.pi", header=TRUE, col.names=c("chrom","start","end","nvariants.Bamako","pi.Bamako"))
Bamako.pi$pop <- "Bamako"

SBF.pi <- read.table("SBF.200.windowed.pi", header=TRUE, col.names=c("chrom","start","end","nvariants.SBF","pi.SBF"))
SBF.pi$pop <- "Burkina gambiae"

MBF.pi <- read.table("MBF.200.windowed.pi", header=TRUE, col.names=c("chrom","start","end","nvariants.MBF","pi.MBF"))
MBF.pi$pop <- "Burkina coluzzii"

BF.pi <- merge(SBF.pi,MBF.pi,by=c("chrom","start","end"),all=TRUE,sort=FALSE)
pi <- merge(BF.pi,Bamako.pi,by=c("chrom","start","end"),all=TRUE,sort=FALSE)
pi <- pi[pi$chrom != "UNKN" & pi$chrom != "Y_unplaced" & pi$chrom != "Mt",]

##Tajima's D
Bamako.D <- read.table("Bamako.200.Tajima.D", header=TRUE, col.names=c("chrom","start","nvariants.Bamako","D.Bamako"))
Bamako.D$pop <- "Bamako"
SBF.D <- read.table("SBF.200.Tajima.D", header=TRUE, col.names=c("chrom","start","nvariants.SBF","D.SBF"))
SBF.D$pop <- "Burkina gambiae"
MBF.D <- read.table("MBF.200.Tajima.D", header=TRUE, col.names=c("chrom","start","nvariants.MBF","D.MBF"))
MBF.D$pop <- "Burkina coluzzii"

BF.D <- merge(SBF.D,MBF.D,by=c("chrom","start"),all=TRUE,sort=FALSE)
D <- merge(BF.D,Bamako.D,by=c("chrom","start"),all=TRUE,sort=FALSE)
D <- D[D$chrom != "UNKN" & D$chrom != "Y_unplaced" & D$chrom != "Mt",]

pi$start <- pi$start -1

D.pi <-merge(D,pi, by=c("chrom","start"),all=TRUE,sort=FALSE)

##fst
BamSBF.fst <- read.table("SBF-Bamako.200K.112515.windowed.weir.fst", header=TRUE, col.names=c("chrom","start","end","variants","Bam_S-BF.fst","Bam_S-BF.unweighted_fst"))

BamMBF.fst <- read.table("MBF-Bamako.200K.112515.windowed.weir.fst", header=TRUE, col.names=c("chrom","start","end","variants","Bam_M-BF.fst","Bam_M-BF.unweighted_fst"))

SMBF.fst <- read.table("BF_S-M.200K.112515.windowed.weir.fst", header=TRUE, col.names=c("chrom","start","end","variants","S_M.fst","S_M.unweighted_fst"))

fst.Bam <- merge(BamSBF.fst, BamMBF.fst,by=c("chrom","start","end"),all=TRUE,sort=FALSE)
fst <- merge(fst.Bam, SMBF.fst,by=c("chrom","start","end"),all=TRUE,sort=FALSE)

fst <- fst[fst$chrom != "UNKN" & fst$chrom != "Y_unplaced" & fst$chrom!="Mt",c("chrom","start","end","Bam_S.BF.fst","Bam_M.BF.fst","S_M.fst")]

fst$start <- fst$start-1

all <- merge(D.pi,fst, by=c("chrom","start"),all=TRUE)

all <- all[all$chrom != "Y_unplaced" & all$chrom != "UNKN" & all$chrom != "Mt",]

all$chrom <- factor(all$chrom, levels=c("2R","2L","3R","3L","X"))

##plotting
pi.plot <- ggplot() + geom_line(data=all, aes(x=start, y=pi.MBF, color="coluzzii")) + geom_line(data=all, aes(x=start, y=pi.SBF, color="gambiae")) + geom_line(data=all, aes(x=start, y=pi.Bamako, color="Bamako")) + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + xlab("") + ylab("pi") + ggtitle("") + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), legend.justification=c(0,1), legend.position=c(0.805,1), legend.title=element_text(size=18), legend.text=element_text(size=15), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), plot.margin=unit(c(-1,.3,-.95,0),"lines"), strip.text.x = element_text(size=18, face="bold"), strip.background=element_rect(fill="white"), legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.9),"lines")) + inversions.func() + inversion.labels() + scale_y_continuous(limits=c(0,.005))

D.plot <- ggplot() + geom_line(data=all, aes(x=start, y=D.MBF, color="coluzzii")) + geom_line(data=all, aes(x=start, y=D.SBF, color="gambiae")) + geom_line(data=all, aes(x=start, y=D.Bamako, color="Bamako")) + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + xlab("") + ylab("Tajima's D") + theme(legend.position="none", axis.ticks=element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,-1,0),"lines")) + inversions.func()

fst.plot <- ggplot() + geom_line(data=all, aes(x=start, y=S_M.fst, color="coluzzii-gambiae")) + geom_line(data=all, aes(x=start, y=Bam_M.BF.fst, color="coluzzii-Bamako"))+ geom_line(data=all, aes(x=start, y=Bam_S.BF.fst, color="Bamako-gambiae")) + xlab("Mb") + ylab("FST") + facet_wrap(~ chrom, scales = "free_x", nrow=1) + theme_bw()  + scale_colour_manual(values=rev(cbbPalette.contrast), name="Comparison") + theme(legend.justification=c(0,1), legend.position=c(0.805,1), legend.title=element_text(size=18), legend.text=element_text(size=15), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), axis.title.x=element_text(size=18), axis.text.x=element_text(size=15), strip.text.x = element_blank(), strip.background=element_blank(), plot.margin=unit(c(-1,.3,0,0),"lines"), legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.9),"lines")) + scale_x_continuous(labels=function(x)x/1000000) + inversions.func()

wpi <- ggplotGrob(pi.plot)
wD <- ggplotGrob(D.plot)
wfst <- ggplotGrob(fst.plot)

maxWidth = grid::unit.pmax(wpi$widths[2:5],wD$widths[2:5])
wpi$widths[2:5] <- as.list(maxWidth)
wD$widths[2:5] <- as.list(maxWidth)
wfst$widths[2:5] <- as.list(maxWidth)

pdf(paste("Figure 3.",format(Sys.Date(), "%m%d%y"),".pdf", sep=""), width=22, height=7)
grid.arrange(wpi, wD, wfst, ncol=1,heights=c(1.2,1,2))
dev.off()

##postscript polishing in Illustrator:
##
##turn "pi" to pi symbol, italicize, move away from axis
##italicize D in Tajima's D
##italicize FST
##italicize inversions
##italicize gambiae and coluzzii in legends
##
# 