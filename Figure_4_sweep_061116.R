setwd("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results")

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(plyr)

cbbPalette.contrast <-c("#0072B2", "#D55E00", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

inversions <- data.frame(
name=c("2La","2Rj","2Rb","2Rc","2Ru"),chrom=c("2L",rep("2R",4)),start=c(20524058,3262186,19023925,26780000,31480000),end=c(42165532,15750717,26758676,31450000,35500000))

inversions.func <- function()
list(
geom_rect(data=inversions[inversions$name=="2La",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"), geom_rect(data=inversions[inversions$name=="2Rj",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_rect(data=inversions[inversions$name=="2Rb",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3,
fill="#999999"),geom_vline(data=inversions[inversions$name=="2Rc",], aes(xintercept=start), colour="black", linetype="longdash"),
geom_rect(data=inversions[inversions$name=="2Rc",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"),geom_vline(data=inversions[inversions$name=="2Rc",], aes(xintercept=end), colour="black", linetype="longdash"),geom_rect(data=inversions[inversions$name=="2Ru",], aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999"))

inversion.labels <- function()
list(geom_text(data=inversions[inversions$name=="2La",], aes(x=(start + end)/2, y=0.05, label="a"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rj",], aes(x=(start + end)/2, y=0.05, label="j"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rb",], aes(x=(start + end)/2, y=0.05, label="b"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Rc",], aes(x=(start + end)/2, y=0.05, label="c"), fontface="bold", color="black", size=7), geom_text(data=inversions[inversions$name=="2Ru",], aes(x=(start + end)/2, y=0.05, label="u"), fontface="bold", color="black", size=7))

##100 kb data

coluzzii.D.100 <- read.table("../AD006.downsampled.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.coluzzii"))
Bamako.D.100 <- read.table("../AD012.downsampled.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.Bamako"))
gambiae.D.100 <- read.table("AD019.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.gambiae"))

ColBam.D.100 <- merge(coluzzii.D.100, Bamako.D.100, by=c("chrom", "pos"),all=TRUE)

D.100 <- merge(ColBam.D.100, gambiae.D.100, by=c("chrom", "pos"),all=TRUE)

coluzzii.pi.100 <- read.table("AD006.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.coluzzii"))
Bamako.pi.100 <- read.table("AD012.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.Bamako"))
gambiae.pi.100 <- read.table("AD019.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.gambiae"))

ColBam.pi.100 <- merge(coluzzii.pi.100, Bamako.pi.100, by=c("chrom", "pos"),all=TRUE)
pi.100 <-merge(ColBam.pi.100, gambiae.pi.100, by=c("chrom", "pos"),all=TRUE)

D.pi.100 <-merge(D.100,pi.100, by=c("chrom","pos"),all=TRUE,sort=FALSE)

fst.100 <- read.table("SMBam.100K.windowed.fst.cleaned", header=FALSE, na.strings="na", col.names=c("chrom", "pos", "nSNPs", "frac.cov", "ave.min.cov", "ColBam", "ColGam", "BamGam"))

D.pi.fst.100 <- merge(D.pi.100,fst.100, by=c("chrom","pos"),all=TRUE)

all.100 <- D.pi.fst.100[D.pi.fst.100$chrom != "Y_unplaced" & D.pi.fst.100$chrom != "UNKN" & D.pi.fst.100$chrom != "Mt",c("chrom","pos","D.coluzzii","D.Bamako","D.gambiae","pi.coluzzii","pi.Bamako","pi.gambiae","ColBam","ColGam","BamGam")]

high_FST <- read.table("SMBam.site.fst.cleaned.BamShigh", col.names=c("chrom","pos","nSNPs","frac.cov","ave.min.cov","coluzzii-Bamako","coluzzii-nonBamako","Bamako-nonBamako"))

high_FST_missense <- read.table("BamS_highFSTSNPs.ann.vcf.missense.pos", col.names=c("chrom","pos"))

###2Rc:

genes_c <- read.table("genes_2Rc.txt", header=FALSE, col.names=c("gene","gene_start","gene_end","exon_start","exon_end"))
genes_c$height <- 0
genes_c[as.numeric(genes_c$gene)%%6 == 0,"height"] <- 0.0
genes_c[as.numeric(genes_c$gene)%%6 == 1,"height"] <- 0.07
genes_c[as.numeric(genes_c$gene)%%6 == 2,"height"] <- 0.14
genes_c[as.numeric(genes_c$gene)%%6 == 3,"height"] <- 0.21
genes_c[as.numeric(genes_c$gene)%%6 == 4,"height"] <- 0.28
genes_c[as.numeric(genes_c$gene)%%6 == 5,"height"] <- 0.35
genes_c[genes_c$gene=="AGAP002859","height"] <- -.07

mins <- aggregate(gene_start ~ gene,genes_c,FUN=min)
maxs <- aggregate(gene_end ~ gene,genes_c,FUN=max)
gene_names <- merge(mins,maxs,by="gene")
gene_names$midpoint <- gene_names$gene_start + (gene_names$gene_end-gene_names$gene_start)/2
gene_names$height <- 0.1
gene_names[as.numeric(gene_names$gene)%%6 == 0,"height"] <- 0.03
gene_names[as.numeric(gene_names$gene)%%6 == 1,"height"] <- 0.10
gene_names[as.numeric(gene_names$gene)%%6 == 2,"height"] <- 0.17
gene_names[as.numeric(gene_names$gene)%%6 == 3,"height"] <- 0.24
gene_names[as.numeric(gene_names$gene)%%6 == 4,"height"] <- 0.31
gene_names[as.numeric(gene_names$gene)%%6 == 5,"height"] <- 0.38
gene_names[gene_names$gene=="AGAP002859","height"] <- -0.04

c.100 <- all.100[all.100$chrom=="2R" & all.100$pos > inversions[inversions$name=="2Rc","start"] & all.100$pos < inversions[inversions$name=="2Rc","end"],]

high_FST.c <- high_FST[high_FST$chrom=="2R" & high_FST$pos > inversions[inversions$name=="2Rc","start"] & high_FST$pos < inversions[inversions$name=="2Rc","end"],]

pi.c <- ggplot() + geom_line(data=c.100, aes(x=pos, y=pi.coluzzii, color="coluzzii ")) + geom_line(data=c.100, aes(x=pos, y=pi.gambiae, color="gambiae ")) + geom_line(data=c.100, aes(x=pos, y=pi.Bamako, color="Bamako ")) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + scale_y_continuous(breaks=c(0,.005,.01,.015,.02,.025), labels=c(0,0.005,0.010,0.015,0.020,0.025)) + ggtitle("") + ylab("pi") + xlab("") + scale_x_continuous(labels=function(x)x/1000000,limits=c(inversions[inversions$name=="2Rc","start"],inversions[inversions$name=="2Rc","end"])) + labs(x=NULL) + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), legend.justification=c(1,0), legend.position=c(1,0), legend.direction="horizontal", legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.75),"lines"), axis.title.y=element_text(size=15), axis.text.y=element_text(size=10), plot.margin=unit(c(-.8,.2,-.5,0),"lines")) + geom_rect(aes(xmin=28419000,xmax=28543000, ymin=-Inf, ymax=Inf), alpha=0.2, fill="#999999")

D.c <- ggplot() + geom_line(data=c.100, aes(x=pos, y=D.coluzzii, color="coluzzii")) + geom_line(data=c.100, aes(x=pos, y=D.gambiae, color="gambiae")) + geom_line(data=c.100, aes(x=pos, y=D.Bamako, color="Bamako")) +  theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + ggtitle("") + ylab("Tajima's D") + xlab("") + theme(legend.position="none", axis.ticks=element_blank(), axis.text.x=element_blank(), plot.margin=unit(c(-.8,.2,-.3,0),"lines"), axis.title.y=element_text(size=15), axis.text.y=element_text(size=10)) + scale_x_continuous(labels=function(x)x/1000000,limits=c(inversions[inversions$name=="2Rc","start"],inversions[inversions$name=="2Rc","end"])) + labs(x=NULL) + geom_rect(aes(xmin=28419000,xmax=28543000, ymin=-Inf, ymax=Inf), alpha=0.2, fill="#999999")

fst.c <- ggplot() + geom_line(data=c.100, aes(x=pos, y=ColGam, color="coluzzii-gambiae")) + geom_line(data=c.100, aes(x=pos, y=ColBam, color="coluzzii-Bamako")) + geom_line(data=c.100, aes(x=pos, y=BamGam, color="Bamako-gambiae")) + ggtitle("") + theme_bw() + ylab("Fst") + xlab("Position, Mbp") + theme( legend.justification=c(.99,1), legend.position=c(.99,1), legend.direction="vertical", legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), legend.key.size= unit(c(0.75),"lines"), axis.title=element_text(size=15), axis.text=element_text(size=10), plot.margin=unit(c(-.85,.2,0,0),"lines")) + scale_x_continuous(labels=function(x)x/1000000,limits=c(inversions[inversions$name=="2Rc","start"],inversions[inversions$name=="2Rc","end"])) + scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6), labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6)) + scale_colour_manual(values=rev(cbbPalette.contrast), name="Comparison") + geom_rect(aes(xmin=28419000,xmax=28543000, ymin=-Inf, ymax=Inf), alpha=0.2, fill="#999999")

genes.SNPs.c <- ggplot(c.100, x=pos, y=0) + ylab("") + xlab("") + theme_bw() + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_blank(), plot.margin=unit(c(-.1,.2,-1,0),"lines")) + xlim(28419000,28543000) +

geom_segment(data=genes_c,aes(x=gene_start,xend=gene_end,y=height,yend=height),color="black",fill="black",size=0.25)+ geom_text(data=gene_names,aes(x=midpoint,y=height,label=gene),fontface="bold",size=3,color="black") +

geom_rect(data=genes_c,aes(xmin=exon_start,xmax=exon_end,ymin=height-.01,ymax=height+.01),color="black",fill="white",size=0.125)+

geom_point(data=high_FST.c, aes(x=pos, y=0.41), position=position_jitter(width=0, height=0.005), color="black", shape=19, size=1) + geom_point(data=high_FST_missense, aes(x=pos, y=0.44), position=position_jitter(width=0, height=0.005), color="gray 40", shape=17, size=2)

wpic <- ggplotGrob(pi.c)
wDc <- ggplotGrob(D.c)
wfstc <- ggplotGrob(fst.c)
wfeaturesc <- ggplotGrob(genes.SNPs.c)

maxWidth = grid::unit.pmax(wpic$widths[2:5],wDc$widths[2:5],wfstc$widths[2:5],wfeaturesc$widths[2:5])
wpic$widths[2:5] <- as.list(maxWidth)
wDc$widths[2:5] <- as.list(maxWidth)
wfstc$widths[2:5] <- as.list(maxWidth)
wfeaturesc$widths[2:5] <- as.list(maxWidth)

pdf(paste("Figure 3c.",format(Sys.Date(), "%m%d%y"),".pdf",sep=""), width=8.5, height=11, useDingbats=FALSE)
grid.arrange(wpic,wDc,wfstc,wfeaturesc,ncol=1,heights=c(1,1,1.2,.8))
dev.off()

