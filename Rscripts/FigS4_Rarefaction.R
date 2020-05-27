# Teresita M. Porter, Dec. 24, 2019

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rarecurve
library(purrr) # for map_dfr
library(ggplot2) # ggplot
library(scales) # comma
library(gridExtra) # grid.arrange
library(cowplot) # get_legend

###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#####################################################################

# Read in cat.csv
A <- read.csv(file="cat.csv", head=TRUE)

# fix sites numbers Pres2: S1, S2, S3; Pres3: S4, S5, S6 in the SampleName field
A$SampleName <- gsub("(Pres3_\\w+_)S1(_\\w*_\\w*_S\\d+)", "\\1S4\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S2(_\\w*_\\w*_S\\d+)", "\\1S5\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S3(_\\w*_\\w*_S\\d+)", "\\1S6\\2", A$SampleName)

#######################################################
# Create dataframes for vegan based on all taxa
######################################################

# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[32:37] <- c("Experiment","Treatment","Site","Replicate", "PrePost", "IlluminaSample")

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.2.esv<-reshape2::dcast(A.1, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.2.esv) <- A.2.esv$SampleName
A.2.esv$SampleName <- NULL

#remove columns with only zeros
esv.notnull<-A.2.esv[,colSums(A.2.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# Exclude controls before doing rarefaction
esv.notnull2.2 <- esv.notnull2[-c(19,32:36,64:65),]

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2.2), prob=0.15)
# 15% 
# 44051.8 

# set random seed
set.seed(1234)

# Do rarefection with pkg 'vegan'
rarecurveout <- rarecurve2(esv.notnull2, 
                          sample=esv.percentile, 
                          step=500, 
                          label=T)

# Reformat vegan list as df (cols OTU, raw.read)
rare.df <- lapply(rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Add sample names to vegan output (df) (rownames)
sample_names <- rownames(esv.notnull2)
names(rare.df) <- sample_names

# Map rownames to vegan output (df)
rare.df <- map_dfr(rare.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
rare.df <- data.frame(rare.df, do.call(rbind, str_split(rare.df$sample,"_")), stringsAsFactors = FALSE)
names(rare.df)[4:9]<-c("Experiment","Treatment","Site","Replicate", "PrePost","IlluminaSample")

# Create factors
rare.df$Site <- factor(rare.df$Site, 
                            levels = c("S1", "S2", "S3","S4", "S5", "S6", ""),
                            labels = c("Laurel 10-3", "Laurel 7-2", "Laurel 4-3","Beaver 19","Beaver 18","Clair 12","Negative Controls"))
rare.df$Treatment <- factor(rare.df$Treatment, 
                            levels = c("ANTI", "ETOH", "ANTInc", "ETOHnc", "EXTnc","PCRnc","SOILEXTnc","WATEREXTnc"),
                            labels = c("Antifreeze", "Ethanol", "Negative Controls", "Negative Controls", "Negative Controls","Negative Controls","Negative Controls","Negative Controls"))
rare.df$Experiment <- factor(rare.df$Experiment,
                         levels = c("Pres2", "Pres3"),
                         labels = c("Paired sampling", "Single split sample"))
rare.df$PrePost <- factor(rare.df$PrePost,
                      levels = c("PRE", "POST", ""),
                      labels = c("Pre evaporation", "Post evaporation", "NA"))

# color by site
p1.tmp <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Site") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Site), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom") +
  guides(color=guide_legend(ncol=3, override.aes = list(size = 2))) 

l1 <- get_legend(p1.tmp)

p1 <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Site") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Site), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none")

# color by treatment
p2.tmp <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Treatment") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Treatment), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank()) +
  guides(color=guide_legend(override.aes = list(size = 2))) 

l2 <- get_legend(p2.tmp)

p2 <- ggplot(data = rare.df) +
  ggtitle("All taxa ~ Treatment") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Treatment), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none")

#######################################################
# Create dataframes for vegan based on arthropoda only
######################################################

# only keep Arthropoda
A.1.arth <- A.1[A.1$Phylum=="Arthropoda",]

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.2.esv<-reshape2::dcast(A.1.arth, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.2.esv) <- A.2.esv$SampleName
A.2.esv$SampleName <- NULL

#remove columns with only zeros
esv.notnull<-A.2.esv[,colSums(A.2.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# Exclude controls before doing rarefaction
esv.notnull2.2 <- esv.notnull2[-c(19,32:36,64:65),]

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2.2), prob=0.15)
# 15% 
# 19898.6 

# set random seed
set.seed(1234)

# Do rarefection with pkg 'vegan'
rarecurveout <- rarecurve2(esv.notnull2, 
                           sample=esv.percentile, 
                           step=500, 
                           label=T)

# Reformat vegan list as df (cols OTU, raw.read)
rare.df <- lapply(rarecurveout, function(x){
  b <- as.data.frame(x)
  b <- data.frame(OTU = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

# Add sample names to vegan output (df) (rownames)
sample_names <- rownames(esv.notnull2)
names(rare.df) <- sample_names

# Map rownames to vegan output (df)
rare.df <- map_dfr(rare.df, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "sample")

# Parse out metadata from sample
rare.df <- data.frame(rare.df, do.call(rbind, str_split(rare.df$sample,"_")), stringsAsFactors = FALSE)
names(rare.df)[4:9]<-c("Experiment","Treatment","Site","Replicate", "PrePost","IlluminaSample")

# Create factors
rare.df$Site <- factor(rare.df$Site, 
                       levels = c("S1", "S2", "S3", "S4", "S5", "S6", ""),
                       labels = c("Laurel 10-3", "Laurel 7-2", "Laurel 4-3", "Beaver 19", "Beaver 18", "Clair 12", "Negative controls"))
rare.df$Treatment <- factor(rare.df$Treatment, 
                            levels = c("ANTI", "ETOH", "ANTInc", "ETOHnc", "EXTnc","PCRnc","SOILEXTnc","WATEREXTnc"),
                            labels = c("Antifreeze", "Ethanol", "Negative Controls", "Negative Controls", "Negative Controls","Negative Controls","Negative Controls","Negative Controls"))
rare.df$Experiment <- factor(rare.df$Experiment,
                             levels = c("Pres2", "Pres3"),
                             labels = c("Paired sampling", "Single split sample"))
rare.df$PrePost <- factor(rare.df$PrePost,
                          levels = c("PRE", "POST", ""),
                          labels = c("Pre evaporation", "Post evaporation", "NA"))

# color by site
p3 <- ggplot(data = rare.df) +
  ggtitle("Arthropoda only ~ Site") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Site), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none")

# color by treatment
p4 <- ggplot(data = rare.df) +
  ggtitle("Arthropoda only ~ Treatment") +
  labs(x="Reads", y="ESVs") +
  geom_point(aes(x = raw.read, y = OTU, color = Treatment), size=0.1) +
  geom_vline(xintercept = esv.percentile, linetype = "dashed") +
  scale_x_continuous(label = comma) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "none")

g <- plot_grid(p1, p2, p3, p4, l1, l2, ncol = 2, rel_heights = c(1, 1, .5))

ggsave("FigS4_rarefaction.pdf", g)
# vertical line represents 15% percentile of read depth across sample (excluding controls)
