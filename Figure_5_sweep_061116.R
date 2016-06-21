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

##read in D

coluzzii.D.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/AD006.downsampled.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.coluzzii"))
Bamako.D.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/AD012.downsampled.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.Bamako"))

gambiae.D.100 <- read.table("AD019.100K.D", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "D.gambiae"))

D.100 <- Reduce(function(...) merge(..., by=c("chrom","pos"), all=TRUE), list(coluzzii.D.100,Bamako.D.100,gambiae.D.100))

##read in pi

coluzzii.pi.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/AD006.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.coluzzii"))
Bamako.pi.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/AD012.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.Bamako"))
gambiae.pi.100 <- read.table("/media/anophele/Hopper/Archives/Bamako/Pools/082015_results/AD019.100K.pi", na.strings="na", col.names=c("chrom", "pos", "n.SNPs", "frac.cov", "pi.gambiae"))

pi.100 <- Reduce(function(...) merge(..., by=c("chrom","pos"), all=TRUE), list(coluzzii.pi.100,Bamako.pi.100,gambiae.pi.100))

##read in FST and combine all three pop gen parameters

fst.100 <- read.table("SMBam.100K.windowed.fst.cleaned", header=FALSE, na.strings="na", col.names=c("chrom", "pos", "nSNPs", "frac.cov", "ave.min.cov", "ColBam", "ColGam", "BamGam"))

all <- Reduce(function(...) merge(..., by=c("chrom","pos"),all=TRUE), list(D.100,pi.100,fst.100))

all <- all[all$chrom != "Y_unplaced" & all$chrom != "UNKN" & all$chrom != "Mt",c("chrom","pos","D.coluzzii","D.Bamako","D.gambiae","pi.coluzzii","pi.Bamako","pi.gambiae","ColBam","ColGam","BamGam")]

##read in locations of high-FST SNPs and those that cause missense mutations

high_FST <- read.table("SMBam.site.fst.cleaned.BamShigh", col.names=c("chrom","pos","nSNPs","frac.cov","ave.min.cov","coluzzii-Bamako","coluzzii-nonBamako","Bamako-nonBamako"))

high_FST_missense <- read.table("BamS_highFSTSNPs.ann.vcf.missense.pos", col.names=c("chrom","pos"))

###read in exons in the swept region:

genes_c <- read.table("genes_2Rc.txt", header=FALSE, col.names=c("gene","gene_start","gene_end","exon_start","exon_end"))
genes_c$height <- 0
genes_c[as.numeric(genes_c$gene)%%6 == 0,"height"] <- 0.0
genes_c[as.numeric(genes_c$gene)%%6 == 1,"height"] <- 0.07
genes_c[as.numeric(genes_c$gene)%%6 == 2,"height"] <- 0.14
genes_c[as.numeric(genes_c$gene)%%6 == 3,"height"] <- 0.21
genes_c[as.numeric(genes_c$gene)%%6 == 4,"height"] <- 0.28
genes_c[as.numeric(genes_c$gene)%%6 == 5,"height"] <- 0.35
genes_c[genes_c$gene=="AGAP002859","height"] <- 0.0
genes_c[genes_c$gene=="AGAP002865" | genes_c$gene=="AGAP002867","height"] <- 0.42

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
gene_names[gene_names$gene=="AGAP002859","height"] <- 0.03
gene_names[gene_names$gene=="AGAP002865" | gene_names$gene=="AGAP002867","height"] <- 0.45
gene_names[gene_names$gene=="AGAP002859","midpoint"] <- 28500000
gene_names[gene_names$gene=="AGAP002865","midpoint"] <- gene_names[gene_names$gene=="AGAP002865","midpoint"] - 5000
gene_names[gene_names$gene=="AGAP002867","midpoint"] <- gene_names[gene_names$gene=="AGAP002867","midpoint"] + 5000


##AGAP002859, AGAP002865, and AGAP002867 are tweaked for aesthetics

##narrow to the c inversion

c.100 <- all[all$chrom=="2R" & all$pos > inversions[inversions$name=="2Rc","start"] & all$pos < inversions[inversions$name=="2Rc","end"],]

high_FST.c <- high_FST[high_FST$chrom=="2R" & high_FST$pos > inversions[inversions$name=="2Rc","start"] & high_FST$pos < inversions[inversions$name=="2Rc","end"],]

pi.c <- ggplot() + geom_line(data=c.100, aes(x=pos, y=pi.coluzzii, color="coluzzii ")) + geom_line(data=c.100, aes(x=pos, y=pi.gambiae, color="gambiae ")) + geom_line(data=c.100, aes(x=pos, y=pi.Bamako, color="Bamako ")) + theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + scale_y_continuous(breaks=c(0,.005,.01,.015,.02,.025), labels=c(0,0.005,0.010,0.015,0.020,0.025)) + ggtitle("") + ylab("pi") + xlab("") + scale_x_continuous(labels=function(x)x/1000000,limits=c(inversions[inversions$name=="2Rc","start"],inversions[inversions$name=="2Rc","end"])) + labs(x=NULL) + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), legend.title=element_text(size=18), legend.text=element_text(size=15), legend.justification=c(1,0), legend.position=c(1,0), legend.direction="vertical", legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15), plot.margin=unit(c(-.8,.2,-.5,0),"lines")) + geom_rect(aes(xmin=28419000,xmax=28543000, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999")

D.c <- ggplot() + geom_line(data=c.100, aes(x=pos, y=D.coluzzii, color="coluzzii")) + geom_line(data=c.100, aes(x=pos, y=D.gambiae, color="gambiae")) + geom_line(data=c.100, aes(x=pos, y=D.Bamako, color="Bamako")) +  theme_bw() + scale_colour_manual(values=cbbPalette.contrast, name="Population") + ggtitle("") + ylab("Tajima's D") + xlab("") + theme(legend.position="none", axis.ticks=element_blank(), axis.text.x=element_blank(), plot.margin=unit(c(-.8,.2,-.3,0),"lines"), axis.title.y=element_text(size=18), axis.text.y=element_text(size=15)) + scale_x_continuous(labels=function(x)x/1000000,limits=c(inversions[inversions$name=="2Rc","start"],inversions[inversions$name=="2Rc","end"])) + labs(x=NULL) + geom_rect(aes(xmin=28419000,xmax=28543000, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999")

fst.c <- ggplot() + geom_line(data=c.100, aes(x=pos, y=ColGam, color="coluzzii-gambiae")) + geom_line(data=c.100, aes(x=pos, y=ColBam, color="coluzzii-Bamako")) + geom_line(data=c.100, aes(x=pos, y=BamGam, color="Bamako-gambiae")) + ggtitle("") + theme_bw() + ylab("FST") + xlab("Mbp") + theme(legend.title=element_text(size=18), legend.text=element_text(size=15), legend.justification=c(1,1), legend.position=c(.985,1), legend.direction="vertical", legend.background=element_rect(fill="transparent"), legend.margin=unit(-1,"lines"), axis.title=element_text(size=18), axis.text=element_text(size=15), plot.margin=unit(c(-.85,.2,-.1,0),"lines")) + scale_x_continuous(labels=function(x)x/1000000,limits=c(inversions[inversions$name=="2Rc","start"],inversions[inversions$name=="2Rc","end"])) + scale_y_continuous(breaks=c(0,.1,.2,.3,.4,.5,.6), labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6)) + scale_colour_manual(values=rev(cbbPalette.contrast), name="Comparison") + geom_rect(aes(xmin=28419000,xmax=28543000, ymin=-Inf, ymax=Inf), alpha=0.3, fill="#999999")

genes.SNPs.c <- ggplot(c.100, x=pos, y=0) + ylab("") + xlab("") + theme_bw() + theme(axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_blank(), plot.margin=unit(c(0.68,.2,-1,0),"lines")) + xlim(28419000,28543000) +

geom_segment(data=genes_c,aes(x=gene_start,xend=gene_end,y=height,yend=height),color="black",fill="black",size=0.25)+ geom_text(data=gene_names,aes(x=midpoint,y=height,label=gene),fontface="bold",size=3,color="black") +

geom_rect(data=genes_c,aes(xmin=exon_start,xmax=exon_end,ymin=height-.01,ymax=height+.01),color="black",fill="white",size=0.125)+

geom_point(data=high_FST.c, aes(x=pos, y=0.53), position=position_jitter(width=0, height=0.005), color="black", fill="#CC79A7", shape=21, size=2) + geom_point(data=high_FST_missense, aes(x=pos, y=0.49), position=position_jitter(width=0, height=0.005), color="black", fill="#CC79A7", shape=24, size=3)

wpic <- ggplotGrob(pi.c)
wDc <- ggplotGrob(D.c)
wfstc <- ggplotGrob(fst.c)
wfeaturesc <- ggplotGrob(genes.SNPs.c)

maxWidth = grid::unit.pmax(wpic$widths[2:5],wDc$widths[2:5],wfstc$widths[2:5],wfeaturesc$widths[2:5])
wpic$widths[2:5] <- as.list(maxWidth)
wDc$widths[2:5] <- as.list(maxWidth)
wfstc$widths[2:5] <- as.list(maxWidth)
wfeaturesc$widths[2:5] <- as.list(maxWidth)

pdf(paste("Figure 5.",format(Sys.Date(), "%m%d%y"),".pdf",sep=""), width=8.5, height=11, useDingbats=FALSE)
grid.arrange(wpic,wDc,wfstc,wfeaturesc,ncol=1,heights=c(1,1,1.2,.8))
dev.off()

##
##post-polishing in Illustrator:
##
##italicize coluzzii, gambiae
##turn "pi" to pi symbol, italicize
##italicize "D" in Tajima's D
##italicize FST
##insert lines showing relationship between genes on bottom and shaded region
##double-check gene names are all readable

