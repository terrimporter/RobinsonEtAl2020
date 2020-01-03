# Teresita M. Porter, Dec. 23, 2019

library(ggplot2)

###################################################################
# Plot proportion of Arthropoda ESVs confidently identified
###################################################################

# Read infiles
A <- read.csv(file="cat.csv", head=TRUE)

# fix sites numbers Pres2: S1, S2, S3; Pres3: S4, S5, S6 in the SampleName field
A$SampleName <- gsub("(Pres3_\\w+_)S1(_\\w*_\\w*_S\\d+)", "\\1S4\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S2(_\\w*_\\w*_S\\d+)", "\\1S5\\2", A$SampleName)
A$SampleName <- gsub("(Pres3_\\w+_)S3(_\\w*_\\w*_S\\d+)", "\\1S6\\2", A$SampleName)

# Select phylum Arthropoda only
Arth <- A[A$Phylum=="Arthropoda",]

# Separate by amplicon
BR5 <- Arth[grepl("BR5", Arth$GlobalESV),]
F230R <- Arth[grepl("F230R", Arth$GlobalESV),]
mljg <- Arth[grepl("ml-jg", Arth$GlobalESV),]

# Calc total taxa and total taxa confidently id'd (species 95%, genus & family 99%)
# See COI Classifier v3 cutoffs st https://github.com/terrimporter/CO1Classifier
BR5.species<-length(unique(BR5$Species))
BR5.species_good<-length(unique(BR5$Species[BR5$cBP>=0.70]))
BR5.genus<-length(unique(BR5$Genus))
BR5.genus_good<-length(unique(BR5$Genus[BR5$pBP>=0.30]))
BR5.family<-length(unique(BR5$Family))
BR5.family_good<-length(unique(BR5$Family[BR5$kBP>=0.20]))

F230R.species<-length(unique(F230R$Species))
F230R.species_good<-length(unique(F230R$Species[F230R$cBP>=0.70]))
F230R.genus<-length(unique(F230R$Genus))
F230R.genus_good<-length(unique(F230R$Genus[F230R$pBP>=0.30]))
F230R.family<-length(unique(F230R$Family))
F230R.family_good<-length(unique(F230R$Family[F230R$kBP>=0.20]))

mljg.species<-length(unique(mljg$Species))
mljg.species_good<-length(unique(mljg$Species[mljg$cBP>=0.70]))
mljg.genus<-length(unique(mljg$Genus))
mljg.genus_good<-length(unique(mljg$Genus[mljg$pBP>=0.30]))
mljg.family<-length(unique(mljg$Family))
mljg.family_good<-length(unique(mljg$Family[mljg$kBP>=0.20]))

# create df for ggplot
BR5.df<-data.frame("rank"=c("species","species","genus","genus","family","family"),
                   "status"=rep(c("all","good"),3),
                   "value"=c(BR5.species, BR5.species_good, BR5.genus, BR5.genus_good, BR5.family, BR5.family_good))

F230R.df<-data.frame("rank"=c("species","species","genus","genus","family","family"),
                   "status"=rep(c("all","good"),3),
                   "value"=c(F230R.species, F230R.species_good, F230R.genus, F230R.genus_good, F230R.family, F230R.family_good))

mljg.df<-data.frame("rank"=c("species","species","genus","genus","family","family"),
                   "status"=rep(c("all","good"),3),
                   "value"=c(mljg.species, mljg.species_good, mljg.genus, mljg.genus_good, mljg.family, mljg.family_good))

# add amplicon column
BR5.df$amplicon <- "BR5"
F230R.df$amplicon <- "F230R"
mljg.df$amplicon <- "ml-jg"

# combine each amplicon df into a single df for ggplot
merged <- rbind(BR5.df, F230R.df, mljg.df)

# create factors
merged$rank = factor(merged$rank, levels=c("species","genus","family"), 
                     labels=c("Species","Genus","Family"))
merged$status = factor(merged$status, levels=c("all","good"),
                       labels=c("All taxa","Confidently identified taxa"))
merged$amplicon = factor(merged$amplicon, levels=c("BR5","F230R", "ml-jg"))

# plot
p <- ggplot(merged, aes(fill=status, y=value, x=rank)) +
  ggtitle("Arthropoda") +
  geom_bar(position="dodge",stat="identity") +
  scale_x_discrete(limits = rev(levels(rank))) +
  labs(x="Rank", y="Unique taxa") +
  facet_wrap(~ amplicon) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5, 'cm'),
        axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        axis.title.x=element_blank())

ggsave("confidentids.pdf", p)

