# Set working directory
setwd("C:/Users/buergern/Downloads/SummerResearch/BiGG2020_CoolAnts/Rcode_NB/Bacteria/results")

# Load libraries needed
library(vegan)
library(ggplot2)

# Generate data into variables (rarefied-table and Metadata)
rare_tab <- read.csv("../data/rarefied-table.csv", header = T)
#sample_tab <- read.csv("../data/Metadata_SeedCoat.csv", header = T)

#############################################################################
# Start using dbrda

treatment <- rare_tab$Treatment
seedType <- rare_tab$SeedType

dat <- rare_tab[,4:5571]

dats <- decostand(dat,"total",1)

ddat <- vegdist(dats, "jaccard")

mod_treat <- dbrda(ddat~ treatment + Condition(seedType))
mod_treat2 <- dbrda(ddat~ treatment)
mod_treat3 <- dbrda(ddat~ seedType + Condition(treatment))
mod_treat4 <- dbrda(ddat~ seedType)

how <- how(nperm=1000, blocks=seedType)
how2 <- how(nperm=1000, blocks=treatment)

anova(mod_treat, by="margin", permutations=how)
anova(mod_treat, by="axis", permutations=how)

anova(mod_treat2, by="margin", permutations=1000)
anova(mod_treat2, by="axis", permutations=1000)

anova(mod_treat3, by="margin", permutations=how2)
anova(mod_treat3, by="axis", permutations=how2)

anova(mod_treat4, by="margin", permutations=1000)
anova(mod_treat4, by="axis", permutations=1000)

smmry <- summary(mod_treat)
smry2 <- summary(mod_treat2)
smry3 <- summary(mod_treat3)
smry4 <- summary(mod_treat4)

#############################################################################
# Making graphs from dbrda stuff

smry <-scores(mod_treat)
df1  <- data.frame(smry$sites[,1:2])       # dbRDA1 and dbRDA2
#df2  <- data.frame(smry$sites[,2])     # loadings for PC1 and PC2

windows(height=4, width=4)
rda.plot <- ggplot(df1, aes(x=dbRDA1, y=dbRDA2, color=treatment, shape=treatment)) +
  # Gets different shapes
  scale_shape_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"),values=c(15,16,17)) +
  # Gets different colors (Color blind friendly)
  scale_color_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"), values=c("#F5793A", "#A95AA1", "#85C0F9")) +
  # Original stuff
  geom_point(aes(label=rownames(treatment)),size=2) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed() +
  # Draw Ellipse
  stat_ellipse() +
  # Gets white background and border
  theme(panel.background = element_rect(fill="white", colour = "grey50")) +
  # Moves legend to inside graph
  theme(legend.title = element_text(size=12, color = "black", face="bold"),
        legend.justification=c(1,0), 
        legend.position=c(0.40, 0.05),  
        legend.background = element_blank(),
        legend.key = element_blank())
rda.plot

smry2 <-scores(mod_treat2)
df2  <- data.frame(smry2$sites[,1:2])       # dbRDA1 and dbRDA2
rda.plot2 <- ggplot(df2, aes(x=dbRDA1, y=dbRDA2, color=treatment, shape=treatment)) +
  # Gets different shapes
  scale_shape_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"), values=c(15,16,17)) +
  # Gets different colors (Color blind friendly)
  scale_color_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"), values=c("#F5793A", "#A95AA1", "#85C0F9")) +
  geom_point(aes(label=rownames(treatment)),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed() +
  stat_ellipse()+
  theme(panel.background = element_rect(fill="white", colour = "grey50")) +
  theme(legend.title = element_text(size=12, color = "black", face="bold"),
        legend.justification=c(1,0), 
        legend.position=c(0.40, 0.05),  
        legend.background = element_blank(),
        legend.key = element_blank())
rda.plot2

smry3 <-scores(mod_treat3)
df3  <- data.frame(smry3$sites[,1:2])       # dbRDA1 and dbRDA2
rda.plot3 <- ggplot(df3, aes(x=dbRDA1, y=dbRDA2, color=seedType, shape=seedType)) +
  # Gets different shapes
  scale_shape_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"), values=c(15,16,17)) +
  # Gets different colors (Color blind friendly)
  scale_color_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"), values=c("#F5793A", "#A95AA1", "#85C0F9")) + 
  geom_point(aes(label=rownames(seedType)),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed() +
  stat_ellipse()+
  theme(panel.background = element_rect(fill="white", colour = "grey50")) +
  theme(legend.title = element_text(size=12, color = "black", face="bold"),
        legend.justification=c(1,0), 
        legend.position=c(0.95, 0.8),  
        legend.background = element_blank(),
        legend.key = element_blank())
rda.plot3

smry4 <-scores(mod_treat4)
df4  <- data.frame(smry4$sites[,1:2])       # dbRDA1 and dbRDA2
rda.plot4 <- ggplot(df4, aes(x=dbRDA1, y=dbRDA2, color=seedType, shape=seedType)) +
  # Gets different shapes
  scale_shape_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"), values=c(15,16,17)) +
  # Gets different colors (Color blind friendly)
  scale_color_manual(name="Seed Type", labels=c("Ant", "Control", "No Elaiosome"), values=c("#F5793A", "#A95AA1", "#85C0F9")) +
  geom_point(aes(label=rownames(seedType)),size=3) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  coord_fixed() +
  stat_ellipse()+
  theme(panel.background = element_rect(fill="white", colour = "grey50")) +
  theme(legend.title = element_text(size=12, color = "black", face="bold"),
        legend.justification=c(1,0), 
        legend.position=c(0.95, 0.8),  
        legend.background = element_blank(),
        legend.key = element_blank())
rda.plot4
