# Teresita M. Porter, Dec. 23, 2019

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(gridExtra) # grid.arrange
library(cowplot) # get_legend

#####################################################################
# Look at richness
#####################################################################

# Read in cat.csv
A <- read.csv(file="cat.csv", head=TRUE)

# fix sites numbers Pres2: S1, S2, S3; Pres3: S4, S5, S6 in the SampleName field
A$SampleName <- gsub("(Pres3_\\w+_)S1(_\\w*_\\w*_S\\d+)", "\\1S4\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S2(_\\w*_\\w*_S\\d+)", "\\1S5\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S3(_\\w*_\\w*_S\\d+)", "\\1S6\\2", A$SampleName)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")), stringsAsFactors = FALSE)
names(A.1)[32:37] <- c("Experiment","Treatment","Site", "Replicate","PrePost", "IlluminaSample")

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

# Exclude negative controls before doing rarefaction
esv.notnull2.2 <- esv.notnull2[-c(19,32:36,64:65),]

#calculate 15th percentile for rrarefy function
esv.percentile<-quantile(rowSums(esv.notnull2.2), prob=0.15)
# 15% 
# 19898.6 

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
rare.mat <- rrarefy(esv.notnull2, sample=esv.percentile)

# Convert to presence-absence matrix
rare.mat[rare.mat>0] <-1

# Convert to df
df<-data.frame(rare.mat, stringsAsFactors = FALSE)  

# Get total ESVs per sample
df$sums<-rowSums(df)

# Move rownames to first column
df2<-data.frame(df, stringsAsFactors = FALSE)
setDT(df2, keep.rownames = TRUE)[]

# Get separate substrate and siterep cols
setDT(df2)[, paste0("S", 1:6) := tstrsplit(rn, "_")]
colnames(df2)[colnames(df2)=="S1"] <- "Experiment"
colnames(df2)[colnames(df2)=="S2"] <- "Treatment"
colnames(df2)[colnames(df2)=="S3"] <- "Site"
colnames(df2)[colnames(df2)=="S4"] <- "Replicate"
colnames(df2)[colnames(df2)=="S5"] <- "PrePost"
colnames(df2)[colnames(df2)=="S6"] <- "IlluminaSample"

# create factors
df2$Treatment <- factor(df2$Treatment,
                        levels = c("ANTI", "ETOH", "ANTInc", "ETOHnc", "EXTnc","PCRnc","SOILEXTnc","WATEREXTnc"),
                        labels = c("Antifreeze", "Ethanol", "Negative Controls", "Negative Controls", "Negative Controls","Negative Controls","Negative Controls","Negative Controls"))
df2$Site <- factor(df2$Site,
                        levels = c("S1", "S2", "S3", "S4", "S5", "S6"),
                        labels = c("Laurel 10-3", "Laurel 7-2", "Laurel 4-3", "Beaver 19", "Beaver 18", "Clair 12"))
df2$Experiment <- factor(df2$Experiment,
                   levels = c("Pres2", "Pres3"),
                   labels = c("Paired sampling", "Single split sample"))
df2$PrePost <- factor(df2$PrePost,
                         levels = c("PRE", "POST", ""),
                         labels = c("Pre evaporation", "Post evaporation", "NA"))


# subset to leave out controls
df3 <- df2[!df2$Treatment=="Negative Controls",]

# split by Experiment
paired <- df3[df3$Experiment=="Paired sampling",]
single <- df3[df3$Experiment=="Single split sample",]

# Compare richness by site
p <- ggplot(paired) +
  geom_point(aes(x=paired$Site, y=paired$sums)) +
  ggtitle("a)") +
  labs(x="Sites", y="Arthropod ESV Richness") +
  theme(legend.title=element_blank()) +
  theme_bw() + 
  facet_grid(cols=vars(Treatment, PrePost)) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1.0,vjust = 0),
        legend.title = element_blank(),
        legend.position = "none")

s <- ggplot(single) +
  geom_point(aes(x=single$Site, y=single$sums)) +
  ggtitle("b)") +
  labs(x="Sites", y="Arthropod ESV Richness") +
  theme(legend.title=element_blank()) +
  theme_bw() + 
  facet_grid(cols=vars(Treatment, PrePost)) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1.0,vjust = 0),
        legend.title = element_blank(),
        legend.position = "none")

g <- grid.arrange(p,s)

ggsave("ESVrichness.pdf", g, width = 8)
# based on normalized data
# Arthropoda only
# includes controls, triplicates shown for every expt+treatment+prepost combination


