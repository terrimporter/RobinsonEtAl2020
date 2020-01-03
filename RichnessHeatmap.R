# Teresita M. Porter, Jan. 2, 2020

library(stringr) # str_split
library(reshape2) # dcast
library(vegan) # rrarefy
library(ggplot2) # ggplot
library(data.table) # setDT
library(gridExtra) # grid.arrange
library(cowplot) # get_legend

# Read in cat.csv
A <- read.csv(file="cat.csv", head=TRUE)

# fix sites numbers Pres2: S1, S2, S3; Pres3: S4, S5, S6 in the SampleName field
A$SampleName <- gsub("(Pres3_\\w+_)S1(_\\w*_\\w*_S\\d+)", "\\1S4\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S2(_\\w*_\\w*_S\\d+)", "\\1S5\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S3(_\\w*_\\w*_S\\d+)", "\\1S6\\2", A$SampleName)

# Create dataframes for vegan
# Split up SampleName with pkg 'stringr'
A.1<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")))
names(A.1)[32:37] <- c("Experiment","Treatment","Site", "Replicate","PrePost", "IlluminaSample")

# Focus on Arthropoda
A.2 <- A.1[A.1$Phylum=="Arthropoda",]

# Focus on EPT
A.3 <- A.2[A.2$Order=="Ephemeroptera" | A.2$Order=="Plecoptera_Insecta" | A.2$Order=="Trichoptera",]

# Focus on Post-evaporation only
A.4 <- A.3[A.3$PrePost=="POST",]

# Separate sites by exeriment
paired <- A.4[A.4$Experiment=="Pres2",]
single <- A.4[A.4$Experiment=="Pres3",]

# pivot to make esv matrix (pool across verions, keep only substrate + sites separate)
paired.family<-reshape2::dcast(paired, SampleName ~ Family, value.var = "ESVsize", fun.aggregate = sum)
single.family<-reshape2::dcast(single, SampleName ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# move sample to rownames then delete
rownames(paired.family) <- paired.family$SampleName
paired.family$SampleName <- NULL

rownames(single.family) <- single.family$SampleName
single.family$SampleName <- NULL

# remove columns with only zeros
paired.notnull<-paired.family[,colSums(paired.family) !=0]

single.notnull<-single.family[,colSums(single.family) !=0]

# remove rows with only zeros & edit rownames
paired.notnull2<-paired.notnull[rowSums(paired.notnull) !=0,]

single.notnull2<-single.notnull[rowSums(single.notnull) !=0,]

# calculate 15th percentile for rrarefy function
paired.percentile<-quantile(rowSums(paired.notnull2), prob=0.15)
# 15% 
# 4839.35

single.percentile<-quantile(rowSums(single.notnull2), prob=0.15)
# 15% 
# 1605.55

# set random seed
set.seed(1234)

# Rarefy the dataset down to the 15th percentile
paired.mat <- rrarefy(paired.notnull2, sample=paired.percentile)
single.mat <- rrarefy(single.notnull2, sample=single.percentile)

# Convert to df
paired.df<-data.frame(paired.mat)  
single.df<-data.frame(single.mat)  

# Get total ESVs per sample
paired.df$sums<-rowSums(paired.df)
single.df$sums<-rowSums(single.df)

# Move rownames to first column
paired.df2<-data.frame(paired.df)
setDT(paired.df2, keep.rownames = TRUE)[]

single.df2<-data.frame(single.df)
setDT(single.df2, keep.rownames = TRUE)[]

# Get separate substrate and siterep cols
setDT(paired.df2)[, paste0("S", 1:6) := tstrsplit(rn, "_")]
colnames(paired.df2)[colnames(paired.df2)=="S1"] <- "Experiment"
colnames(paired.df2)[colnames(paired.df2)=="S2"] <- "Treatment"
colnames(paired.df2)[colnames(paired.df2)=="S3"] <- "Site"
colnames(paired.df2)[colnames(paired.df2)=="S4"] <- "Replicate"
colnames(paired.df2)[colnames(paired.df2)=="S5"] <- "PrePost"
colnames(paired.df2)[colnames(paired.df2)=="S6"] <- "IlluminaSample"

setDT(single.df2)[, paste0("S", 1:6) := tstrsplit(rn, "_")]
colnames(single.df2)[colnames(single.df2)=="S1"] <- "Experiment"
colnames(single.df2)[colnames(single.df2)=="S2"] <- "Treatment"
colnames(single.df2)[colnames(single.df2)=="S3"] <- "Site"
colnames(single.df2)[colnames(single.df2)=="S4"] <- "Replicate"
colnames(single.df2)[colnames(single.df2)=="S5"] <- "PrePost"
colnames(single.df2)[colnames(single.df2)=="S6"] <- "IlluminaSample"

# create factors
paired.df2$Treatment <- factor(paired.df2$Treatment,
                        levels = c("ANTI", "ETOH", "ANTInc", "ETOHnc", "EXTnc","PCRnc","SOILEXTnc","WATEREXTnc"),
                        labels = c("Antifreeze", "Ethanol", "Negative Controls", "Negative Controls", "Negative Controls","Negative Controls","Negative Controls","Negative Controls"))
paired.df2$Site <- factor(paired.df2$Site,
                      levels = c("S1", "S2", "S3"),
                     labels = c("Laurel 10-3", "Laurel 7-2", "Laurel 4-3"))
paired.df2$Experiment <- factor(paired.df2$Experiment,
                   levels = c("Pres2", "Pres3"),
                   labels = c("Paired sampling", "Single split sample"))
paired.df2$PrePost <- factor(paired.df2$PrePost,
                         levels = c("PRE", "POST", ""),
                         labels = c("Pre evaporation", "Post evaporation", "NA"))
paired.df2$Replicate <- factor(paired.df2$Replicate,
                      levels = c("A", "B", "C"))


single.df2$Treatment <- factor(single.df2$Treatment,
                               levels = c("ANTI", "ETOH", "ANTInc", "ETOHnc", "EXTnc","PCRnc","SOILEXTnc","WATEREXTnc"),
                               labels = c("Antifreeze", "Ethanol", "Negative Controls", "Negative Controls", "Negative Controls","Negative Controls","Negative Controls","Negative Controls"))
single.df2$Site <- factor(single.df2$Site,
                          levels = c("S4", "S5", "S6"),
                          labels = c("Beaver 19", "Beaver 18", "Clair 12"))
single.df2$Experiment <- factor(single.df2$Experiment,
                                levels = c("Pres2", "Pres3"),
                                labels = c("single sampling", "Single split sample"))
single.df2$PrePost <- factor(single.df2$PrePost,
                             levels = c("PRE", "POST", ""),
                             labels = c("Pre evaporation", "Post evaporation", "NA"))
single.df2$Replicate <- factor(single.df2$Replicate,
                               levels = c("A", "B", "C"))

# remove unneeded columns
paired.df3 <- paired.df2[,-c(1,19,25)]
single.df3 <- single.df2[,-c(1,22,28)]

# melt for ggplot
paired.df4 <- melt(paired.df3, id=c("Experiment","Treatment","Site","Replicate","PrePost"))
single.df4 <- melt(single.df3, id=c("Experiment","Treatment","Site","Replicate","PrePost"))

# create factor
paired.df4$variable <- factor(paired.df4$variable,
                    levels=rev(unique(paired.df4$variable)))
single.df4$variable <- factor(single.df4$variable,
                              levels=rev(unique(single.df4$variable)))

# Compare richness by site
p.tmp <- ggplot(paired.df4) +
  geom_tile(aes(x=Treatment, y=variable, fill=value)) +
  ggtitle("a)") +
  labs(x="Sites", y="EPT Families") +
  theme(legend.title=element_blank()) +
  scale_fill_gradient(trans="log10", na.value="white", low="white", high="black") +
  theme_bw() + 
  facet_grid(cols=vars(Site,Replicate)) +
  guides(fill=guide_legend(title="Reads")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1.0,vjust = 0, size=7),
        axis.text.y = element_text(size=7),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size=8),
        legend.text = element_text(size=8))

l <- get_legend(p.tmp)

p <- ggplot(paired.df4) +
  geom_tile(aes(x=Treatment, y=variable, fill=value)) +
  ggtitle("a)") +
  labs(x="Sites", y="EPT Families") +
  theme(legend.title=element_blank()) +
  scale_fill_gradient(trans="log10", na.value="white", low="white", high="black") +
  theme_bw() + 
  facet_grid(cols=vars(Site,Replicate)) +
  guides(fill=guide_legend(title="Reads")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9),
        legend.position = "none",
        strip.text.x = element_text(size = 7),
        strip.text.y = element_text(size = 7))

s <- ggplot(single.df4) +
  geom_tile(aes(x=Treatment, y=variable, fill=value)) +
  ggtitle("b)") +
  labs(x="Sites", y="EPT Families") +
  theme(legend.title=element_blank()) +
  scale_fill_gradient(trans="log10", na.value="white", low="white", high="black") +
  theme_bw() + 
  facet_grid(cols=vars(Site,Replicate)) +
  guides(fill=guide_legend(title="Reads")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1.0,vjust = 0, size=7),
        axis.text.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=9),
        legend.position = "none",
        strip.text.x = element_text(size = 7),
        strip.text.y = element_text(size = 7))

g <- plot_grid(p, s, l, nrow = 3, rel_heights = c(1, 1, 0.2))

ggsave("Family_heatmap.pdf", g, width = 8, height = 11)
# based on normalized/rarefied data
# Arthropoda, EPT only, Post evaporation only
# controls excluded, shown for every treatment+replicate combination
 