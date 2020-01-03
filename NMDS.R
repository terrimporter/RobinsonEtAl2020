# Teresita M. Porter, Jan. 3, 2019

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(goeveg) # scree
library(plyr) # ddply

# Read in cat.csv
A <- read.csv(file="cat.csv", head=TRUE)

# fix sites numbers Pres2: S1, S2, S3; Pres3: S4, S5, S6 in the SampleName field
A$SampleName <- gsub("(Pres3_\\w+_)S1(_\\w*_\\w*_S\\d+)", "\\1S4\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S2(_\\w*_\\w*_S\\d+)", "\\1S5\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S3(_\\w*_\\w*_S\\d+)", "\\1S6\\2", A$SampleName)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[32:37] <- c("Experiment","Treatment","Site","Replicate","PrePost", "IlluminaSample")

# Focus on Arthropoda
A.2 <- A.1[A.1$Phylum=="Arthropoda",]

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
A.3.esv<-reshape2::dcast(A.2, SampleName ~ GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(A.3.esv) <- A.3.esv$SampleName
A.3.esv$SampleName <- NULL

#remove columns with only zeros
esv.notnull<-A.3.esv[,colSums(A.3.esv) !=0]

#remove rows with only zeros & edit rownames
esv.notnull2<-esv.notnull[rowSums(esv.notnull) !=0,]

# Exclude controls before doing rarefaction and plotting
esv.notnull2.2 <- esv.notnull2[-c(19,32:36,64:65),]

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2.2), prob=0.15)
# 15% 
# 19898.6 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2.2, sample=esv.percentile)

# Convert to presence-absence matrix
rare.mat[rare.mat>0] <-1

# Scree plots to determine number of dimensions to use for NMDS, use k=3
# pdf("Scree.pdf")
# # check dims
# dimcheckMDS(rare.mat)
# dev.off()

# Do 2 dimensional NMDS
nmds2<-metaMDS(rare.mat, k=2, trymax=100)
# stress = 0.1242202

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot.pdf")
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n", main="SSU")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()
# linear R2 = 0.929

# Create grouping matrix for samples by grabbing row names from above matrix
names<-data.frame(row.names(rare.mat), stringsAsFactors = FALSE)

# Rename the column
names(names)<-"sample"

# Copy column to row names
row.names(names)<-names$sample

# Split first column into their own fields
names.1<-data.frame(names, do.call(rbind, strsplit(names$sample,'_')), stringsAsFactors = FALSE)
names(names.1)[2:7]<-c("Experiment","Treatment", "Site", "Replicate", "PrePost", "IlluminaSample")

# Remove first column
names.1 <- names.1[,-1]

# Grab sites/species scores from NMDS output
df <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df for ggplot
gg <- merge(df, names.1, by="row.names")

# create factors
gg$Treatment <- factor(gg$Treatment,
                        levels = c("ANTI", "ETOH"),
                        labels = c("Antifreeze", "Ethanol"))
gg$Site <- factor(gg$Site,
                        levels = c("S1", "S2", "S3","S4","S5","S6"),
                        labels = c("Laurel 10-3", "Laurel 7-2", "Laurel 4-3",
                                   "Beaver 19", "Beaver 18", "Clair 12"))
gg$Experiment <- factor(gg$Experiment,
                             levels = c("Pres2", "Pres3"),
                             labels = c("Paired sampling", "Single split sample"))
gg$PrePost <- factor(gg$PrePost,
                          levels = c("PRE", "POST", ""),
                          labels = c("Pre evaporation", "Post evaporation", "NA"))

# color by treatment
chulls.treatment <- ddply(gg, .(Treatment), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by treatment
p1 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.treatment, aes(x=NMDS1, y=NMDS2, fill=Treatment), alpha=0.5) +
  geom_point(data=gg, aes(color=Treatment)) +
  ggtitle("a)") +
  theme_bw() +
  theme(
    plot.title = element_text(size=10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10))

# color by site
chulls.site <- ddply(gg, .(Site), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by site
p2 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.site, aes(x=NMDS1, y=NMDS2, fill=Site), alpha=0.5) +
  geom_point(data=gg, aes(color=Site)) +
  ggtitle("b)") +
  theme_bw() +
  theme(
    plot.title = element_text(size=10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10))

# color by PrePost
chulls.prepost <- ddply(gg, .(PrePost), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by site
p3 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.prepost, aes(x=NMDS1, y=NMDS2, fill=PrePost), alpha=0.5) +
  geom_point(data=gg, aes(color=PrePost)) +
  ggtitle("c)") +
  theme_bw() +
  theme(
    plot.title = element_text(size=10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10))

# color by PrePost
chulls.experiment <- ddply(gg, .(Experiment), function(gg) gg[chull(gg$NMDS1, gg$NMDS2), ])

# NMDS plot, color by site
p4 <- ggplot(data=gg, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls.experiment, aes(x=NMDS1, y=NMDS2, fill=Experiment), alpha=0.5) +
  geom_point(data=gg, aes(color=Experiment)) +
  ggtitle("d)") +
  theme_bw() +
  theme(
    plot.title = element_text(size=10),
    axis.text.x = element_text(hjust=1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    legend.key=element_blank(),
    axis.title = element_text(size=10),
    axis.text = element_text(size=10),
    legend.title = element_blank(),
    legend.text = element_text(size=10))

g <- grid.arrange(p1, p2, p3, p4, nrow=2)

ggsave("NMDS.pdf", g)

# Create metadata from rownames 'sample'
# remove all Pre-evaporation samples for balanced design to compare experiment and treatment
env <- gg[,c(1,4:9)]
env <- env[!grepl("PRE",env$Row.names),]

# edit rare.mat to remove all Pre-evaporation samples for balanced design
rare.mat2 <- rare.mat[!grepl("PRE",rownames(rare.mat)),]

# Assess dispersion (variance) using ANOVA
# Create distance matrix based on P-A data using Bray Curtis (Sorensen) dissimilarity
sor<-vegdist(rare.mat2, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
# Break it down by study to ensure balanced design
bd.site<-betadisper(sor, as.factor(env$Site))
bd.treatment<-betadisper(sor, as.factor(env$Treatment))
bd.experiment<-betadisper(sor, as.factor(env$Experiment))

# check for heterogeneity of beta dispersions within groups BALANCED DESIGN
set.seed(1234)
anova(bd.site) # 0.04519 *
anova(bd.treatment) # n/s
anova(bd.experiment) # n/s

pdf("BetaDispersion.pdf")
par(mfrow=c(2,2))
boxplot(bd.site, main="Site")
boxplot(bd.treatment, main="Treatment")
boxplot(bd.experiment, main="Experiment")
dev.off()

# Use ADONIS to test significance of groupings 
adonis(sor~Experiment*Treatment, data=env, permutations=999, strata=env$Site)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Experiment            1    2.4696 2.46956  8.1522 0.18770  0.001 ***
#   Treatment             1    0.4999 0.49993  1.6503 0.03800  0.001 ***
#   Experiment:Treatment  1    0.4940 0.49396  1.6306 0.03754  0.001 ***
#   Residuals            32    9.6938 0.30293         0.73676           
# Total                35   13.1572                 1.00000  
