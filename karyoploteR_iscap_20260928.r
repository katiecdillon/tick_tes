setwd("/Users/katiedillon/Documents/Documents\ -\ BadDNA’s\ MacBook\ Pro/Dissertation/Ch3/karyoploter/iscap")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("karyoploteR")
library(karyoploteR)
library(data.table)
library(GenomicRanges)

## Read filtered GFF into R data.table
te.dt <- fread("iscap_TEs_filtered.out",
               header = FALSE,
               sep = " ",
               quote = "",
               fill = TRUE)

## Assign column names
setnames(te.dt, c("score","div","del","ins","seqname",
                  "start","end","left","strand",
                  "repeat_name","class_family",
                  "rep_begin","rep_end","rep_left","ID","overlap_flag"))

## Convert data.table to GRanges object
all.TEs <- GRanges(
  seqnames   = te.dt$seqname,
  ranges     = IRanges(start = te.dt$start, end = te.dt$end),
  strand     = ifelse(te.dt$strand == "C", "-", "+"),
  TE_name    = te.dt$repeat_name,
  TE_class   = te.dt$class_family,
  score      = te.dt$score,
  divergence = te.dt$div
)

######################################################
### --- PLOT: Genome-wide TE density
######################################################
pp <- getDefaultPlotParams(plot.type = 6)
pp$leftmargin   <- 0.12
pp$rightmargin  <- 0.01
pp$topmargin    <- 5
pp$bottommargin <- 25

custom.genome <- toGRanges("iscap.genome.bed")

## Axis scale and vertical extent of the curve within each panel
axis.max  <- 3400
tick.vals <- c(0,850,1700,2550,3400)
d.r0 <- 0
d.r1 <- 0.85   # leaves the top 15% of each panel empty

png("Figure_hibiscus_v1.png", width = 180, height = 220, units = "mm", res = 500)

kp_all <- plotKaryotype(genome = custom.genome,
                        plot.type = 6,
                        ideogram.plotter = NULL,
                        labels.plotter = NULL,      # turn off the default chromosome names
                        cex = 1,
                        plot.params = pp)

# Chromosome names, shifted left
kpAddChromosomeNames(kp_all, xoffset = -0.03, cex = 0.6)

# Density curve, scaled to the same maximum as the axis
kp_all <- kpPlotDensity(kp_all, all.TEs, window.size = 1e6,
                        data.panel = "ideogram", col = "#7fbc41",
                        border = NA, r0 = d.r0, r1 = d.r1,
                        ymax = axis.max)

## y-axis: line and ticks only, with blank labels
kpAxis(kp_all, ymin = 0, ymax = axis.max, r0 = d.r0, r1 = d.r1,
       data.panel = "ideogram", numticks = 5, side = 1,
       labels = rep("", 5))

## y-axis numbers added manually, one set per chromosome
for (chr in seqlevels(custom.genome)) {
  kpText(kp_all, chr = chr, x = 0,
         y = tick.vals / axis.max,
         labels = format(tick.vals, big.mark = ","),
         data.panel = "ideogram", r0 = d.r0, r1 = d.r1,
         cex = 0.4, pos = 2)
}

## y-axis title
kpAddLabels(kp_all, labels = "TEs/Mbp", side = "left", r0 = 0.5, r1 = d.r1,
            data.panel = "ideogram", cex = 0.4, label.margin = 0.05, srt = 90)

kpAddBaseNumbers(kp_all)

dev.off()

######################################################
### --- TABLE: Density per chromosome
######################################################
wins   <- unlist(tile(custom.genome, width = 1e6))
counts <- countOverlaps(wins, all.TEs)
length(wins)
# Mean and median TE hits/window
summary(counts)

# Mean and median TE hits/chromosome
tapply(counts, as.character(seqnames(wins)), mean)
tapply(counts, as.character(seqnames(wins)), median)

chr.order  <- paste0("chr", 1:14)
chr.vector <- as.character(seqnames(wins))

## One data frame with both statistics
chr.stats <- data.frame(
  chromosome = chr.order,
  mean       = as.numeric(tapply(counts, chr.vector, mean)[chr.order]),
  median     = as.numeric(tapply(counts, chr.vector, median)[chr.order])
)
chr.stats$chromosome <- factor(chr.stats$chromosome, levels = chr.order)

View(chr.stats)
fwrite(chr.stats, "Table_tulip_v1.tsv", sep = "\t")

######################################################
### --- Plot: Boxplot of genome-wide TE hits/window
######################################################
threshold <- quantile(counts, 0.99)
threshold
sum(counts > threshold)

png("Figure_lily_v1.png", width = 120, height = 150, units = "mm", res = 500)

boxplot(counts,
        main = "TE Hits per 1 Mb Window\nPalLabHiFi Genome-Wide",
        ylab = "TE hits per window (log scale)",
        col = "#7fbc41",
        border = "black",
        log = "y")

abline(h = threshold, col = "#c51b7d", lty = 2, lwd = 1.5)
text(x = 1.35, y = threshold, labels = "99th percentile", 
     col = "#c51b7d", cex = 0.7, pos = 3)

dev.off()


######################################################
### --- Plot: Outliers
######################################################
## Identify the outlier windows (top 1% by TE hit density)
threshold <- quantile(counts, 0.99)
outlier.wins <- wins[counts > threshold]

## Find every TE hit that overlaps an outlier window
ov <- findOverlaps(all.TEs, outlier.wins)
dense.TEs <- all.TEs[queryHits(ov)]

## Add a major class column, same logic used earlier for the stacked plot
dense.TEs$major_class <- sub("\\?.*", "", sub("/.*", "", dense.TEs$TE_class))

## Count of hits by class and chromosome
class.chr.counts <- as.data.frame(table(
  chromosome = as.character(seqnames(dense.TEs)),
  class      = dense.TEs$major_class
))

## Keep only combinations that actually occur, sort by count, clean up
class.chr.counts <- class.chr.counts[class.chr.counts$Freq > 0, ]
setnames(class.chr.counts, "Freq", "hit_count")
class.chr.counts <- class.chr.counts[order(-class.chr.counts$hit_count), ]
rownames(class.chr.counts) <- NULL

## Restrict to the six main TE classes
class.chr.counts <- class.chr.counts[class.chr.counts$class %in% 
                                       c("Unknown", "RC", "DNA", "LTR", "SINE", "LINE"), ]
class.chr.counts$chromosome <- sub("^chr([1-9])$", "chr0\\1", class.chr.counts$chromosome)
View(class.chr.counts)

## Plot
outliers <- ggplot(class.chr.counts,
                   aes(x = chromosome,
                       y = hit_count,
                       fill = class)) +
  geom_bar(stat = "identity") +
  theme(panel.grid.major = element_blank(), # remove grid lines
        panel.grid.minor = element_blank(), # remove grid lines
        panel.background = element_blank(), # remove background
        legend.key.size = unit(0.3, "cm"), # legend key box size
        axis.text=element_text(size=7), # axis text font
        axis.text.x=element_text(angle = 45,
                                 hjust = 1),
        axis.title=element_text(size=7), # axis title font
        legend.text=element_text(size=7), # legend text font
        legend.title=element_text(size=7), # legend title font
        axis.line = element_line(colour = "black"), # make axes black
        text = element_text(family = "Arial"))
outliers
ggsave("Figure_sunflower_v1.png", units = "mm", width = 80, height = 80)



## Which repeat sequences are the unknown outliers coming from?
dense.dt <- as.data.table(as.data.frame(dense.TEs))
unknown.families <- dense.dt[major_class == "Unknown", .N, by = TE_name][order(-N)]
View(unknown.families)


######################################################
### --- PLOT: Genome-wide TE density + classes
######################################################
te.dt[, major_class := sub("\\?.*", "", sub("/.*", "", class_family))]
table(te.dt$major_class)

all.TEs <- GRanges(
  seqnames    = te.dt$seqname,
  ranges      = IRanges(start = te.dt$start, end = te.dt$end),
  strand      = ifelse(te.dt$strand == "C", "-", "+"),
  TE_name     = te.dt$repeat_name,
  TE_class    = te.dt$class_family,
  major_class = te.dt$major_class,
  score       = te.dt$score,
  divergence  = te.dt$div
)

pp <- getDefaultPlotParams(plot.type = 6)
pp$leftmargin   <- 0.15
pp$rightmargin  <- 0.01
pp$topmargin    <- 5
pp$bottommargin <- 25

custom.genome <- toGRanges("iscap.genome.bed")

major.classes <- c("LINE", "LTR", "DNA", "RC", "SINE", "All TEs")
class.colors <- c(LINE = "#1b9e77", LTR = "#d95f02", DNA = "#7570b3", 
                  RC = "#e7298a", SINE = "#66a61e", "All TEs" = "black")

png("Figure_daisy_v1.png", width = 180, height = 280, units = "mm", res = 500)

kp <- plotKaryotype(genome = custom.genome,
                    plot.type = 6,
                    ideogram.plotter = NULL,
                    labels.plotter = NULL,
                    cex = 1,
                    plot.params = pp)

kpAddChromosomeNames(kp, xoffset = 0, cex = 0.6)
n <- length(major.classes)

for (i in seq_along(major.classes)) {
  cls <- major.classes[i]
  gr.cls <- if (cls == "All TEs") all.TEs else all.TEs[all.TEs$major_class == cls]
  
  at <- autotrack(current.track = i, total.tracks = n, margin = 0.05, r0 = 0, r1 = 1)
  
  ## Leave a small gap at the top of this band, same idea as d.r1 before
  band.r0 <- at$r0
  band.r1 <- at$r0 + (at$r1 - at$r0) * 0.85
  
  kp <- kpPlotDensity(kp, gr.cls, window.size = 1e6, data.panel = "ideogram",
                      col = class.colors[cls], border = NA,
                      r0 = band.r0, r1 = band.r1)
  
  kpAddLabels(kp, labels = cls, r0 = band.r0, r1 = band.r1, data.panel = "ideogram",
              cex = 0.3, side = "left", label.margin = 0.02)
}

kpAddBaseNumbers(kp)

dev.off()