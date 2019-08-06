#################################### 20181007_16S_rRNAgene_AppleSoil #####################################
# Date: October 11th 2018
# By : AF. Bintarti

# INSTALL PACKAGES
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("ggpubr")
install.packages("car")
install.packages("agricolae")
install.packages("multcompView")
install.packages("gridExtra")
install.packages("ggplot2")
install.packages("sjmisc")
install.packages("sjPlot")
install.packages("MASS")
install.packages("FSA")
install.packages("rcompanion")
install.packages("onewaytests")
install.packages("PerformanceAnalytics")
install.packages("gvlma")
install.packages("userfriendlyscience")
install.packages("ggpmisc")
install.packages("fitdistrplus")
install.packages('BiocManager')
install.packages("cowplot")
library(BiocManager)
library(vegan)
library(plyr)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggplot2)
library(reshape)
library(ggpubr)
library(car)
library(agricolae)
library(multcompView)
library(grid)
library(gridExtra)
library(sjmisc)
library(sjPlot)
library(MASS)
library(FSA)
library(rcompanion)
library(onewaytests)
library(ggsignif)
library(PerformanceAnalytics)
library(gvlma)
library(userfriendlyscience)
library(ggpmisc)
library(tibble)
library(fitdistrplus)

# SET THE WORKING DIRECTORY
setwd('/Users/arifinabintarti/Documents/PAPER_Bintarti_2019_Apple/')
wd <- print(getwd())
# READ OTU BACTERIA
otu <- read.table('OTU_rarefied_16S.txt', sep='\t', header=T, row.names = 1)
taxonomy <- otu[,'taxonomy']
taxonomy
otu <- otu[,-46]
# RAREFACTION CURVE - BACTERIA
rarecurve(t(otu_table(PHYL)), step=50, cex=0.5, xlab = "Reads", ylab = "Bacterial/Archaeal OTUs")
# READ OTU FUNGI
otuITS <- read.table(file = "OTU_rarefied_ITS.txt", sep='\t', header=T, row.names = 1) 
# RAREFACTION CURVE - FUNGI
rarecurve(t(otu_table(physeq)), step=50, cex=0.5, xlab = "Reads", ylab = "Fungal OTUs")
# READ MAP
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
map$Site <- as.character(map$Site)
# CALCULATE THE ALPHA DIVERSITY (SHANNON, RICHNESS, PIELOU/EVENNESS) FOR BACTERIA
otu_rare_PA <- 1*(otu>0)
s <- specnumber(otu, MARGIN = 2)
richness <- as.data.frame(s)
h <- diversity(t(otu), index = 'shannon')
shannon <- as.data.frame(h)
pielou <- h/log(s)
evenness <- as.data.frame(pielou)
map_df <- data.frame(map) # make data frame
map.div <- map_df
map.div$Rootstock <- map.div$rootstock
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou
# CALCULATE THE ALPHA DIVERSITY (SHANNON, RICHNESS, PIELOU/EVENNESS) FOR FUNGI
sITS <- specnumber(otuITS, MARGIN = 2)
hITS <- diversity(t(otuITS), index = 'shannon')
pielouITS <- hITS/log(sITS)
map.div$fg.Richness <- sITS
map.div$fg.Shannon <- hITS
map.div$fg.Pielou <- pielouITS
# MELT MAP ID AND MEASURE VARIABLES
map.alpha <- melt(map.div, id.vars=c('Site', 'cultivar', 'Rootstock'), measure.vars=c('Richness','Shannon', 'Pielou','fg.Richness','fg.Shannon', 'fg.Pielou'))
# CALCULATE THE ALPHA DIVERSITY (SHANNON, RICHNESS, PIELOU/EVENNESS) FOR NEMATODE, MYCORRHIZAL FUNGI, OLIGOCHAETES
## make data frame of nematode group only (11 group)
nema <- map_df[,c(1,24:34)]
dim(nema)
nema
row.names(nema) <- nema$sample_code
nema[1] <- NULL
nema.t <- t(nema)
nema.t
sort(rowSums (nema.t, na.rm = FALSE, dims = 1), decreasing = FALSE)
## make data frame of nematodes, oligochaetes, and mycorrhizal fungi
nema.plus <- map_df[,c(1,24:36)]
row.names(nema.plus) <- nema.plus$sample_code
nema.plus[1] <- NULL
nema.plus <- t(nema.plus)
sqrt.nema.plus <- sqrt(nema.plus) # square root the nematode count data
otu_dist_nema.plus <- vegdist(sqrt.nema.plus, method='bray') # calculate dissimilarity indices using 'bray-curtis' method
otu_pcoa_nema.plus <- cmdscale(otu_dist_nema.plus, eig=T) # CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
# calculate alpha diversity
nema_PA <- 1*(nema.t>0)
nema_total=as.data.frame(rowSums(nema.t)) # calculate total count data for every nematode group
colnames(nema_total)="Count"
nema_total <- rownames_to_column(nema_total, var = "Group")
nema.total.Count=colSums(nema.t)
sqrt.nema.t=t(sqrt(nema))
nema.total.SqrtCount=colSums(sqrt.nema.t)
map.div$nema.total.Count=nema.total.Count
map.div$nema.total.SqrtCount=nema.total.SqrtCount
# log10 of mycorrhizal fungi
log.MF=log10(map.div$MycorrhizalFungi)
map.div$log.MF=log.MF
# Square root of oligochaetes
sqrt.OL=sqrt(map.div$Oligochaetes)
map.div$sqrt.OL=sqrt.OL
## Shannon and Richness of 11 nematode groups
s.nema <- specnumber(nema.t, MARGIN = 2) # Richness is computed as the number of taxa present per sample
h.nema <- diversity(nema, index = 'shannon')
pielou.nema <- h.nema/log(s.nema)
map.div$nema.Richness <- s.nema
map.div$nema.Shannon <- h.nema
map.div$nema.Pielou <- pielou.nema
names(map.div)

###################################################################################################################
#### ANOVA AND KRUSKAL-WALLIS TO COMPARE BACTERIAL AND FUNGAL RICHNESS & SHANNON AMONG SITES AND ROOTSTOCKS #######
###################################################################################################################

map_aov <- map.div
class(map_aov)
map_aov$Site<-as.factor(map_aov$Site) # inform R that Site is factor
map_aov$Roostock<-as.factor(map_aov$Rootstock) # inform R that Rootstock is factor
str(map_aov, give.attr=F)

# 1. ANOVA AND KRUSKAL-WALLIS TEST FOR BACTERIAL ALPHA DIVERSITY

# 1. Compare richness among sites - Kruskal-Wallis
Aov_richness_site <- lm(map_aov$Richness ~ Site, data=map_aov, na.action=na.exclude)
Aov_richness_site
drop1(Aov_richness_site,~.,test="F") # type III SS and F Tests
# testing assumptions
# Generate residual and predicted values
RC_site_resids <- residuals(Aov_richness_site)
RC_site_preds <- predict(Aov_richness_site)
#Look at a plot of residual vs. predicted values
plot(RC_site_resids ~ RC_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_site_resids) # ALERT!! p-value = 0.0006956, data errors are not normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ Site, data=map_aov, na.action=na.exclude) # GOOD variances among group are homogenous
boxplot(map_aov$Richness ~ Site, data = map_aov) # there are no outliers
plot(density(RC_site_resids)) # density is not bad
qqnorm(RC_site_resids)
qqline(RC_site_resids) # I think the data normality is fine
hist(RC_site_resids)
skew_xts <-  skewness(RC_site_resids)
kurtosis(RC_site_resids,method = 'sample')
# Do Kruskal-Wallis Test
map_aov = mutate(map_aov, Site = factor(Site, levels=unique(Site)))
Summarize(Richness ~ Site, data = map_aov)
kruskal.test(Richness ~ Site, data = map_aov) # Kruskal-Wallis chi-squared = 32.155, df = 19, p-value = 0.03002
ggboxplot(map_aov, x = "Site", y = "Richness")+
  stat_compare_means()
# Do Post Hoc Dunn's Test
DT_RC_site <- dunnTest(Richness~Site, map_aov, method = "bh", kw=TRUE)
print(DT_RC_site,dunn.test.results=TRUE)
DT_RC_site$res
DT_RC_site.df <- as.data.frame(DT_RC_site$res)
write.csv(DT_RC_site.df, file = "DT_RC_site.df.csv")
DT_RC_site_letter = cldList(P.adj ~ Comparison,
        data = DT_RC_site$res,
        threshold = 0.05)
# Do Plot
(bac_rich_site <- ggplot(map.div, aes(x=Site, y=Richness))+
 geom_boxplot() +
 geom_point() +
  geom_signif(comparisons = list(c("3", "15")), annotations = "*", textsize = 6,
              map_signif_level=TRUE, vjust = 0.7)+
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
  theme_bw()+
  labs(title = "A")+
 # geom_text(data=new.richness.summarized,aes(x=Site,y=32+max.Richness,label=new.richness.summarized$groups),vjust=0)+
 theme(axis.text.x=element_blank(),
       axis.ticks.x = element_blank(), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_blank(), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)), 
       axis.title.x = element_blank(),
       axis.title.y = element_text(size=18,face="bold"), 
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 2. Compare richness among rootstocks - Kruskal-Wallis
Aov_richness_rootstock <- lm(Richness ~ Rootstock, data=map_aov, na.action=na.exclude)
drop1(Aov_richness_rootstock,~.,test="F") # type III SS and F Tests
summary(Aov_richness_rootstock)
# testing assumptions
# Generate residual and predicted values
RC_root_resids <- residuals(Aov_richness_rootstock)
RC_root_preds <- predict(Aov_richness_rootstock)
plot(RC_root_resids ~ RC_root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_root_resids) # GOOD p-value = 0.5531, data errors are normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ Rootstock, data=map_aov, na.action=na.exclude) # ALERT!! p-val=0.001476, variances among group are not homogenous
#Plotting
boxplot(Richness ~ Rootstock, data = map_aov) # there are no outliers
plot(density(RC_root_resids))
qqnorm(RC_root_resids)
qqline(RC_root_resids) 
hist(RC_root_resids)
# Do Kruskal-Wallis Test
map_aov = mutate(map_aov, Rootstock = factor(Rootstock, levels=unique(Rootstock)))
Summarize(Richness ~ Rootstock, data = map_aov)
kruskal.test(Richness ~ Rootstock, data = map_aov) # Kruskal-Wallis chi-squared = 16.367, df = 7, p-value = 0.02197
# Do Post Hoc Dunn's Test
DT_RC_root <- dunnTest(Richness~Rootstock, map_aov, method = "bh", kw=TRUE)
print(DT_RC_root,dunn.test.results=TRUE)
DT_RC_root$res
DT_RC_root.df <- as.data.frame(DT_RC_root$res)
write.csv(DT_RC_root.df, file = "DT_RC_root.df.csv")
DT_RC_root_letter = cldList(P.adj ~ Comparison,
        data = DT_RC_root$res,
        threshold = 0.05)
# Do Plot
(bac_rich_root <- ggplot(map.div, aes(x=Rootstock, y=Richness))+
 geom_boxplot() +
 geom_point() +
  theme_bw()+
  labs(title = "C")+
# geom_text(data=new.root.richness.summarized,aes(x=rootstock,y=32+max.Richness,label=new.root.richness.summarized$groups),vjust=0)+
 theme(axis.text.x=element_blank(),
       axis.ticks.x = element_blank(), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_blank(), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)), 
       axis.title.x = element_blank(),
       axis.title.y = element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 3. Compare shannon index among sites - Kruskal-Wallis
Aov_shannon_site <- lm(Shannon ~ Site, data=map_aov, na.action=na.exclude)
drop1(Aov_shannon_site,~.,test="F") # type III SS and F Tests
SH_site_resids <- residuals(Aov_shannon_site)
SH_site_preds <- predict(Aov_shannon_site)
plot(SH_site_resids ~ SH_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(Shannon ~ Site, data=map_aov, na.action=na.exclude) # GOOD the data is homogen
#Test for normality
shapiro.test(SH_site_resids) # ALERT!! p-val<0.05, the data are not normally distributed
plot(density(SH_site_resids))
qqnorm(SH_site_resids)
qqline(SH_site_resids)
hist(SH_site_resids)
skew_xts <- skewness(SH_site_resids)
# Do Kruskal-Wallis Test
kruskal.test(Shannon ~ Site, data = map_aov) # Kruskal-Wallis chi-squared = 35.368, df = 19, p-value = 0.01261
# Do Post Hoc Dunn's Test
DT_SH_site <- dunnTest(Shannon~Site, map_aov, method = "bh", kw=TRUE)
print(DT_SH_site,dunn.test.results=TRUE)
DT_SH_site$res
DT_SH_site.df <- as.data.frame(DT_SH_site$res)
write.csv(DT_SH_site.df, file = "DT_SH_site.df.csv")
DT_SH_site_letter = cldList(P.adj ~ Comparison,
        data = DT_SH_site$res,
        threshold = 0.05)
# Do Plot
(bac_sha_site <- ggplot(map.div, aes(x=Site, y=Shannon))+
 geom_boxplot() +
 geom_point() +
  geom_signif(comparisons = list(c("3", "15")), annotations="*", y_position = 7.68, tip_length = 0.03,textsize=6, vjust = 0.7)+
  geom_signif(comparisons = list(c("3", "18")), annotations="*", y_position = 7.7, tip_length = 0.03,textsize=6, vjust = 0.7)+
#geom_signif(comparisons = list(c("3", "15"), c("3", "18")), annotations =c("p=0.041","p=0.047"), map_signif_level=TRUE)+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
  theme_bw()+
  labs(title = "B")+
#geom_text(data=new.shannon.summarized,aes(x=Site,y=0.03+max.Shannon,label=new.shannon.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))
 
# 4. Compare shannon index among rootstocks - ANOVA
Aov_shannon_rootstock <- lm(Shannon ~ Rootstock, data=map_aov, na.action=na.exclude)
drop1(Aov_shannon_rootstock,~.,test="F") # type III SS and F Tests
SH_root_resids <- residuals(Aov_shannon_rootstock)
SH_root_preds <- predict(Aov_shannon_rootstock)
plot(SH_root_resids ~ SH_root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(Shannon ~ Rootstock, data=map_aov, na.action=na.exclude) # GOOD p-val>0.05, the data is homogen
#Test for normality
shapiro.test(SH_root_resids) # GOOD p-val>0.05, the data are normally distributed
qqnorm(SH_root_resids)
qqline(SH_root_resids) 
# Do Tukey's HSD Post Hoc Test
hsd_Shannon_rootstock<- HSD.test(Aov_shannon_rootstock, "Rootstock", group = TRUE, console = TRUE)
hsd_Shannon_rootstock<- HSD.test(Aov_shannon_rootstock, "Rootstock", group = F, console = TRUE)
# Do Plot
# add significance letters from HSD.test into box plot
root.shannon.summarized <- map.div %>% group_by(Rootstock) %>% summarize(max.Shannon=max(Shannon))
hsd_Shannon_rootstock <- HSD.test(Aov_shannon_rootstock, "Rootstock", group = TRUE, console = TRUE)
hsd_Shannon = hsd_Shannon_rootstock$groups
class(hsd_Shannon)
hsd_Shannon$Rootstock <- rownames(hsd_Shannon)
new.root.shannon.summarized=left_join(hsd_Shannon,root.shannon.summarized, by='Rootstock')  

(bac_sha_root <- ggplot(map.div, aes(x=Rootstock, y=Shannon))+
 geom_boxplot() +
 geom_point() +
  theme_bw()+
  labs(title = "D")+
 geom_text(data=new.root.shannon.summarized,aes(x=Rootstock,y=0.01+max.Shannon,label=new.root.shannon.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=14,angle=49,hjust =0.9),
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=20,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_text(vjust=10),
       axis.title=element_text(size=18,face="bold", vjust = 10),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 5. Arrange plots for bacterial alpha diversity
dev.off()
dev.set(dev.next())
grid.newpage()
grid.draw(rbind(ggplotGrob(bac_rich_site), ggplotGrob(bac_sha_site),size = "first"))
grid.newpage()
grid.draw(rbind(ggplotGrob(bac_rich_root), ggplotGrob(bac_sha_root),size = "first"))

# 2. ANOVA TEST FOR FUNGAL ALPHA DIVERSITY

# 1. Compare richness among site
AovITS_richness_site <- lm(fg.Richness ~ Site , data=map_aov, na.action=na.exclude)
drop1(AovITS_richness_site,~.,test="F") # type III SS and F Tests
summary(AovITS_richness_site)
ITSRichness_site_resids <- residuals(AovITS_richness_site)
ITSRichness_site_preds <- predict(AovITS_richness_site)
plot(ITSRichness_site_resids ~ ITSRichness_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Richness ~ Site, data=map_aov, na.action=na.exclude) # GOOD p-val>0.05, the data is homogen
#Test for normality
shapiro.test(ITSRichness_site_resids) # GOOD p-val>0.05, the data are normally distributed
qqnorm(ITSRichness_site_resids)
qqline(ITSRichness_site_resids) 
# Do Tukey's HSD Post Hoc Test
hsdITS_Richness_site <- HSD.test(AovITS_richness_site, "Site", alpha = 0.05,group = FALSE,main = NULL,console=TRUE)
hsdITS_Richness_site <- HSD.test(AovITS_richness_site, "Site", alpha = 0.05,group = T ,main = NULL,console=TRUE)
# Do Plot
# add significance letters from HSD.test into box plot
ITSsite_richness.summarized <- map.div %>% group_by(Site) %>% summarize(max.Richness=max(fg.Richness))
hsdITS_Richness_site <- HSD.test(AovITS_richness_site, "Site", group = TRUE, console = TRUE)
hsdITS_Richness = hsdITS_Richness_site$groups
class(hsdITS_Richness)
hsdITS_Richness$Site <- rownames(hsdITS_Richness)
new.ITSsite.richness.summarized=left_join(hsdITS_Richness,ITSsite_richness.summarized, by='Site') 

(fg_rich_site <- ggplot(map.div, aes(x=Site, y=fg.Richness))+
 geom_boxplot() +
 geom_point() +
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 theme_bw()+
  labs(title = "A", y="Richness")+
 geom_text(data=new.ITSsite.richness.summarized,aes(x=Site,y=10+max.Richness,label=new.ITSsite.richness.summarized$groups),vjust=0)+
 theme(axis.text.x=element_blank(),
       axis.ticks.x = element_blank(), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_blank(), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)), 
       axis.title.x = element_blank(),
       axis.title.y = element_text(size=18,face="bold"), 
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 2. Compare richness among rootstocks
AovITS_richness_rootstock <- lm(fg.Richness ~ Rootstock , data=map_aov, na.action=na.exclude)
drop1(AovITS_richness_rootstock,~.,test="F") # type III SS and F Tests
summary(AovITS_richness_rootstock)
ITSRichness_root_resids <- residuals(AovITS_richness_rootstock)
ITSRichness_root_resids
ITSRichness_root_preds <- predict(AovITS_richness_rootstock)
plot(ITSRichness_root_resids ~ ITSRichness_root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Richness ~ Rootstock, data=map_aov, na.action=na.exclude) #p-val > 0.05, the data is homogen
#Test for normality
shapiro.test(ITSRichness_root_resids) ## p-val>0.05, the data are normally distributed
qqnorm(ITSRichness_root_resids)
qqline(ITSRichness_root_resids)
# Do Tukey's HSD Post Hoc Test
hsdITS_Richness_rootstock <- HSD.test(AovITS_richness_rootstock, "Rootstock", alpha = 0.05,group = TRUE,main = NULL,console=TRUE)
hsdITS_Richness_rootstock <- HSD.test(AovITS_richness_rootstock, "Rootstock", alpha = 0.05,group = F,main = NULL,console=TRUE)
# Do Plot
# add significance letters from HSD.test into box plot
ITSroot_richness.summarized <- map.div %>% group_by(Rootstock) %>% summarize(max.Richness=max(fg.Richness))
hsdITS_Richness_rootstock <- HSD.test(AovITS_richness_rootstock, "Rootstock", group = TRUE, console = TRUE)
hsdITS_Richness = hsdITS_Richness_rootstock$groups
class(hsdITS_Richness)
hsdITS_Richness$Rootstock <- rownames(hsdITS_Richness)
hsdITS_Richness = mutate(hsdITS_Richness, Rootstock = factor(Rootstock, levels=unique(Rootstock)))
new.ITSroot_richness.summarized=left_join(hsdITS_Richness,ITSroot_richness.summarized, by='Rootstock') 

(fg_rich_root <- ggplot(map.divITS, aes(x=Rootstock, y=Richness))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
  labs(title = "C", y="Richness")+
 geom_text(data=new.ITSroot_richness.summarized,aes(x=Rootstock,y=10+max.Richness,label=new.ITSroot_richness.summarized$groups),vjust=0)+
 theme(axis.text.x=element_blank(),
       axis.ticks.x = element_blank(), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_blank(), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)), 
       axis.title.x = element_blank(),
       axis.title.y = element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 3. Compare Shannon index among sites
AovITS_shannon_site <- lm(fg.Shannon ~ Site , data=map_aov, na.action=na.exclude)
drop1(AovITS_shannon_site,~.,test="F") # type III SS and F Tests
summary(AovITS_shannon_site)
ITSShannon_site_resids <- residuals(AovITS_shannon_site)
ITSShannon_site_preds <- predict(AovITS_shannon_site)
plot(ITSShannon_site_resids ~ ITSShannon_site_preds , xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Shannon ~ Site, data=map_aov, na.action=na.exclude) #p-val > 0.05, the data is homogen
#Test for normality
shapiro.test(ITSShannon_site_resids) ## p-val>0.05, the data are normally distributed
qqnorm(ITSShannon_site_resids)
qqline(ITSShannon_site_resids)
# Do Tukey's HSD Post Hoc Test
hsdITS_Shannon_site <- HSD.test(AovITS_shannon_site, "Site", alpha = 0.05,group = T ,main = NULL,console=TRUE)
hsdITS_Shannon_site <- HSD.test(AovITS_shannon_site, "Site",alpha = 0.05, group = FALSE, main = NULL,console=TRUE)
# Do Plot
# add significance letters from HSD.test into box plot
ITSsite_shannon.summarized <- map.div %>% group_by(Site) %>% summarize(max.Shannon=max(fg.Shannon))
hsdITS_Shannon_site <- HSD.test(AovITS_shannon_site, "Site", group = TRUE, console = TRUE)
hsdITS_Shannon = hsdITS_Shannon_site$groups
class(hsdITS_Shannon)
hsdITS_Shannon$Site <- rownames(hsdITS_Shannon)
new.ITSsite.shannon.summarized=left_join(hsdITS_Shannon,ITSsite_shannon.summarized, by='Site')  

(fg_sha_site <- ggplot(map.div, aes(x=Site, y=fg.Shannon))+
 geom_boxplot() +
 geom_point() +
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 theme_bw()+
  labs(title = "B", y="Shannon")+
 geom_text(data=new.ITSsite.shannon.summarized,aes(x=Site,y=0.03+max.Shannon,label=new.ITSsite.shannon.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 4. Compare Shannon index among rootstocks
AovITS_shannon_rootstock <- lm(fg.Shannon ~ Rootstock , data=map_aov, na.action=na.exclude)
drop1(AovITS_shannon_rootstock,~.,test="F") # type III SS and F Tests
summary(AovITS_shannon_rootstock)
ITSShannon_root_resids <- residuals(AovITS_shannon_rootstock)
ITSShannon_root_resids
ITSShannon_root_preds <- predict(AovITS_shannon_rootstock)
plot(ITSShannon_root_resids ~ ITSShannon_root_preds , xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Shannon ~ Rootstock , data=map_aov, na.action=na.exclude) #p-val > 0.05, the data is homogen
#Test for normality
shapiro.test(ITSShannon_root_resids) ## p-val>0.05, the data are normally distributed
qqnorm(ITSShannon_root_resids)
qqline(ITSShannon_root_resids)
# Do Tukey's HSD Post Hoc Test
hsdITS_sha_root <- HSD.test(AovITS_shannon_rootstock, "Rootstock", alpha = 0.05,group = TRUE,main = NULL,console=TRUE)
hsdITS_sha_root <- HSD.test(AovITS_shannon_rootstock, "Rootstock", alpha = 0.05,group = F,main = NULL,console=TRUE)
# Do Plot
# add significance letters from HSD.test into box plot
ITSroot_shannon.summarized <- map.div %>% group_by(Rootstock) %>% summarize(max.Shannon=max(fg.Shannon))
hsdITS_sha_root <- HSD.test(AovITS_shannon_rootstock, "Rootstock", group = TRUE, console = TRUE)
hsdITS_Shanon_root = hsdITS_sha_root$groups
class(hsdITS_Shanon_root)
hsdITS_Shanon_root$Rootstock <- rownames(hsdITS_Shanon_root)
new.ITSroot.shannon.summarized=left_join(hsdITS_Shanon_root,ITSroot_shannon.summarized, by='Rootstock')  

(fg_sha_root <- ggplot(map.div, aes(x=Rootstock, y=fg.Shannon))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
  labs(title = "D", y="Shannon")+
 geom_text(data=new.ITSroot.shannon.summarized,aes(x=Rootstock,y=0.03+max.Shannon,label=new.ITSroot.shannon.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=14,angle=49,hjust =0.9),
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=20,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_text(vjust=10),
       axis.title=element_text(size=18,face="bold", vjust = 10),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 5. Arrange plots for fungal alpha diversity
dev.off()
dev.set(dev.next())
grid.newpage()
grid.draw(rbind(ggplotGrob(fg_rich_site), ggplotGrob(fg_sha_site),size = "first"))
grid.newpage()
grid.draw(rbind(ggplotGrob(fg_rich_root), ggplotGrob(fg_sha_root),size = "first"))

##########################################################################################
#### CORRELATION TEST OF BACTERIAL AND FUNGAL ALPHA DIVERSITY WITH SOIL PROPERTIES #######
##########################################################################################

# 1. Pearson correlation of soil properties with bacterial Richness
bac.rich_pH <- cor(map.div$pH, map.div$Richness, method ="pearson")
cor.test(map.div$pH, map.div$Richness)

bac.rich_P <- cor(map.div$P_ppm, map.div$Richness, method ="pearson")
cor.test(map.div$P_ppm, map.div$Richness)

bac.rich_K <- cor.test(map.div$K_ppm, map.div$Richness, method ="pearson")
cor.test(map.div$K_ppm, map.div$Richness)

bac.rich_Ca <- cor.test(map.div$Ca_ppm, map.div$Richness, method ="pearson")
cor.test(map.div$Ca_ppm, map.div$Richness)

bac.rich_Mg <- cor.test(map.div$Mg_ppm, map.div$Richness, method ="pearson")
cor.test(map.div$Mg_ppm, map.div$Richness)

bac.rich_NO3N <- cor.test(map.div$NO3N_ppm, map.div$Richness, method ="pearson")
cor.test(map.div$NO3N_ppm, map.div$Richness)

bac.rich_NH4N <- cor.test(map.div$NH4N_ppm, map.div$Richness, method ="pearson")
cor.test(map.div$NH4N_ppm, map.div$Richness)

bac.rich_OM <- cor.test(map.div$OM_percent, map.div$Richness, method ="pearson")
cor.test(map.div$OM_percent, map.div$Richness)

bac.rich_SD <- cor.test(map.div$sand_percent, map.div$Richness, method ="pearson")
cor.test(map.div$sand_percent, map.div$Richness)

bac.rich_SL <- cor.test(map.div$silt_percent, map.div$Richness, method ="pearson")
cor.test(map.div$silt_percent, map.div$Richness)

bac.rich_CL <- cor.test(map.div$clay_percent, map.div$Richness, method ="pearson")
cor.test(map.div$clay_percent, map.div$Richness)

# effect of soil type to bacterial richness
bac.rich_soil.type <- lm(Richness ~ soil_type, data=map.div)
summary(bac.rich_soil.type)
drop1(bac.rich_soil.type,~.,test="F") # type III SS and F Tests
TYPE_RC_resids <- residuals(bac.rich_soil.type)
TYPE_RC_preds <- predict(bac.rich_soil.type)
plot(TYPE_RC_resids ~ TYPE_RC_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
leveneTest(Richness ~ soil_type, data=map.div, na.action=na.exclude) # p-val 0.004, data is not homogen
shapiro.test(TYPE_RC_resids) # data is normal
# fail test assumption for anova
# use WELCH-ANOVA instead
bac.rich_soil.type <- oneway.test(Richness ~ soil_type, data=map.div, var.equal = F)
bac.rich_soil.type
ggplot(map.div, aes(x=soil_type, y=Richness)) + geom_boxplot()
# F = 13.568, num df = 2.0000, denom df = 4.2992, p-value = 0.01389
# GAMES-HOWELL posthoc test
GH.bac.rich_soil.type <- oneway(map.div$soil_type, y = map.div$Richness, posthoc = 'games-howell')
GH.bac.rich_soil.type

# 2. Pearson correlation of soil properties with bacterial Shannon index
bac.sha._pH <- cor(map.div$pH, map.div$Shannon, method ="pearson")
cor.test(map.div$pH, map.div$Shannon)

bac.sha_P <- cor(map.div$P_ppm, map.div$Shannon, method ="pearson")
cor.test(map.div$P_ppm, map.div$Shannon)

bac.sha_K <- cor.test(map.div$K_ppm, map.div$Shannon, method ="pearson")
cor.test(map.div$K_ppm, map.div$Shannon)

bac.sha_Ca <- cor.test(map.div$Ca_ppm, map.div$Shannon, method ="pearson")
cor.test(map.div$Ca_ppm, map.div$Shannon)

bac.sha_Mg <- cor.test(map.div$Mg_ppm, map.div$Shannon, method ="pearson")
cor.test(map.div$Mg_ppm, map.div$Shannon)

bac.sha_NO3N <- cor.test(map.div$NO3N_ppm, map.div$Shannon, method ="pearson")
cor.test(map.div$NO3N_ppm, map.div$Shannon)

bac.sha_NH4N <- cor.test(map.div$NH4N_ppm, map.div$Shannon, method ="pearson")
cor.test(map.div$NH4N_ppm, map.div$Shannon)

bac.sha_OM <- cor.test(map.div$OM_percent, map.div$Shannon, method ="pearson")
cor.test(map.div$OM_percent, map.div$Shannon)

bac.sha_SD <- cor.test(map.div$sand_percent, map.div$Shannon, method ="pearson")
cor.test(map.div$sand_percent, map.div$Shannon)

bac.sha_SL <- cor.test(map.div$silt_percent, map.div$Shannon, method ="pearson")
cor.test(map.div$silt_percent, map.div$Shannon)

bac.sha_CL <- cor.test(map.div$clay_percent, map.div$Shannon, method ="pearson")
cor.test(map.div$clay_percent, map.div$Shannon)

# effect of soil type to bacterial shannon index # not significant
bac.sha_soil.type <- lm(Shannon ~ soil_type, data=map.div)
summary(bac.sha_soil.type)
drop1(bac.sha_soil.type,~.,test="F") # type III SS and F Tests
TYPE_SH_resids <- residuals(bac.sha_soil.type)
TYPE_SH_preds <- predict(bac.sha_soil.type)
plot(TYPE_SH_resids ~ TYPE_SH_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
leveneTest(Shannon ~ soil_type, data=map.div, na.action=na.exclude) # p-val 0.005, data is not homogen
shapiro.test(TYPE_SH_resids) # data is normal
# fail test assumption for anova
# use WELCH-ANOVA instead
bac.sha_soil.type <- oneway.test(Shannon ~ soil_type, data=map.div, var.equal = F)
bac.sha_soil.type
ggplot(map.div, aes(x=soil_type, y=Shannon)) + geom_boxplot()
# F = 7.1633, num df = 2.0000, denom df = 2.8821, p-value = 0.07615
# use GAMES-HOWELL posthoc
GH.bac.sha_soil.type <- oneway(map.div$soil_type, y = map.div$Shannon, posthoc = 'games-howell')
GH.bac.sha_soil.type

# 3. Pearson correlation of soil properties with fungal Richness # all are not significant
fg.rich_pH <- cor(map.div$pH, map.div$fg.Richness, method ="pearson")
cor.test(map.div$pH, map.div$fg.Richness)

fg.rich_P <- cor(map.div$P_ppm, map.div$fg.Richness, method ="pearson")
cor.test(map.div$P_ppm, map.div$fg.Richness)

fg.rich_K <- cor.test(map.div$K_ppm, map.div$fg.Richness, method ="pearson")
cor.test(map.div$K_ppm, map.div$fg.Richness)

fg.rich_Ca <- cor.test(map.div$Ca_ppm, map.div$fg.Richness, method ="pearson")
cor.test(map.div$Ca_ppm, map.div$fg.Richness)

fg.rich_Mg <- cor.test(map.div$Mg_ppm, map.div$fg.Richness, method ="pearson")
cor.test(map.div$Mg_ppm, map.div$fg.Richness)

fg.rich_NO3N <- cor.test(map.div$NO3N_ppm, map.div$fg.Richness, method ="pearson")
cor.test(map.div$NO3N_ppm, map.div$fg.Richness)

fg.rich_NH4N <- cor.test(map.div$NH4N_ppm, map.div$fg.Richness, method ="pearson")
cor.test(map.div$NH4N_ppm, map.div$fg.Richness)

fg.rich_OM <- cor.test(map.div$OM_percent, map.div$fg.Richness, method ="pearson")
cor.test(map.div$OM_percent, map.div$fg.Richness)

fg.rich_SD <- cor.test(map.div$sand_percent, map.div$fg.Richness, method ="pearson")
cor.test(map.div$sand_percent, map.div$fg.Richness)

fg.rich_SL <- cor.test(map.div$silt_percent, map.div$fg.Richness, method ="pearson")
cor.test(map.div$silt_percent, map.div$fg.Richness)

fg.rich_CL <- cor.test(map.div$clay_percent, map.div$fg.Richness, method ="pearson")
cor.test(map.div$clay_percent, map.div$fg.Richness)

# effect of soil type to fungal richness # not significant
fg.rich_soil.type <- lm(fg.Richness ~ soil_type, data=map.div)
summary(fg.rich_soil.type)
drop1(fg.rich_soil.type,~.,test="F") # type III SS and F Tests
TYPE_fg.RC_resids <- residuals(fg.rich_soil.type)
TYPE_fg.RC_preds <- predict(fg.rich_soil.type)
plot(TYPE_fg.RC_resids ~ TYPE_fg.RC_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
leveneTest(fg.Richness ~ soil_type, data=map.div, na.action=na.exclude) # p-val 0.2, data is homogen
shapiro.test(TYPE_fg.RC_resids) # 0.047, data is not normal
# fail test assumption for anova

# 4. Pearson correlation of soil properties with fungal Shannon index
fg.sha_pH <- cor(map.div$pH, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$pH, map.div$fg.Shannon)

fg.sha_P <- cor(map.div$P_ppm, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$P_ppm, map.div$fg.Shannon)

fg.sha_K <- cor.test(map.div$K_ppm, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$K_ppm, map.div$fg.Shannon)

fg.sha_Ca <- cor.test(map.div$Ca_ppm, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$Ca_ppm, map.div$fg.Shannon)

fg.sha_Mg <- cor.test(map.div$Mg_ppm, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$Mg_ppm, map.div$fg.Shannon)

fg.sha_NO3N <- cor.test(map.div$NO3N_ppm, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$NO3N_ppm, map.div$fg.Shannon)

fg.sha_NH4N <- cor.test(map.div$NH4N_ppm, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$NH4N_ppm, map.div$fg.Shannon)

fg.sha_OM <- cor.test(map.div$OM_percent, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$OM_percent, map.div$fg.Shannon)

fg.sha_SD <- cor.test(map.div$sand_percent, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$sand_percent, map.div$fg.Shannon) # significant

fg.sha_SL <- cor.test(map.div$silt_percent, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$silt_percent, map.div$fg.Shannon) # significant

fg.sha_CL <- cor.test(map.div$clay_percent, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$clay_percent, map.div$fg.Shannon) # significant

# effect of soil type to fungal shannon index # significant
fg.sha_soil.type <- lm(fg.Shannon ~ soil_type, data=map.div)
summary(fg.sha_soil.type)
drop1(fg.sha_soil.type,~.,test="F") # type III SS and F Tests
TYPE_fg.SH_resids <- residuals(fg.sha_soil.type)
TYPE_fg.SH_preds <- predict(fg.sha_soil.type)
plot(TYPE_fg.SH_resids ~ TYPE_fg.SH_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
leveneTest(fg.Shannon ~ soil_type, data=map.div, na.action=na.exclude) # p-val 0.6, data is homogen
shapiro.test(TYPE_fg.SH_resids) # data is normal
# Tukey's HSD post hoc test
hsd_TYPE_fg.SH <- HSD.test(fg.sha_soil.type, "soil_type", alpha = 0.05,group = FALSE,main = NULL,console=TRUE)
ggplot(map.div, aes(x=soil_type, y=fg.Shannon)) + geom_boxplot()

#########################################################################################
## LINEAR REGRESSION TEST OF BACTERIAL AND FUNGAL ALPHA DIVERSITY WITH SOIL PROPERTIES ##
#########################################################################################

### BACTERIA ###
# 1. BACTERIAL RICHNESS - SAND
Sand_percent_rich <- lm(Richness ~ sand_percent, data=map.div)
summary(Sand_percent_rich)
ggscatter(map.div, x = "sand_percent", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Sand (%)", ylab = "Richness")+
          annotate("text", x=42, y=5300, label = "y == 14.88(x)+3863.22", parse=T)+
          annotate("text", x=38, y=5200, label = "R^2 == 0.325", parse=T)+
          annotate("text", x=40, y=5100, label = "p-val == 4.38e-05", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 2. BACTERIAL RICHNESS - SILT
Silt_percent_rich <- lm(Richness ~ silt_percent, data=map.div)
summary(Silt_percent_rich)
ggscatter(map.div, x = "silt_percent", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Silt (%)", ylab = "Richness")+
          annotate("text", x=39, y=5300, label = "y == -16.48(x)+5158.09", parse=T)+
          annotate("text", x=35, y=5200, label = "R^2 == 0.208", parse=T)+
          annotate("text", x=36, y=5100, label = "p-val == 0.001", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())
    
# 3. BACTERIAL RICHNESS - CLAY
Clay_percent_rich <- lm(Richness ~ clay_percent, data=map.div)
summary(Clay_percent_rich)
ggscatter(map.div, x = "clay_percent", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Clay (%)", ylab = "Richness")+
          annotate("text", x=20, y=5300, label = "y == -49.89(x)+5456.83", parse=T)+
          annotate("text", x=20, y=5200, label = "R^2 == 0.457", parse=T)+
          annotate("text", x=20, y=5100, label = "p-val == 3.376e-07", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())
    
#. 4. BACTERIAL RICHNESS - K
K_rich <- lm(Richness ~ K_ppm, data=map.div)
ggplot(K_rich, aes(x = K_ppm, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(K_rich)
ggscatter(map.div, x = "K_ppm", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "K (ppm)", ylab = "Richness")+
          annotate("text", x=200, y=5300, label = "y == -2.42(x)+5016.79", parse=T)+
          annotate("text", x=200, y=5200, label = "R^2 == 0.157", parse=T)+
          annotate("text", x=200, y=5100, label = "p-val == 0.006", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

#. 4. BACTERIAL RICHNESS - Ca 
par(mar=c(6, 4, 4, 2) + 0.1)
Ca_rich <- lm(Richness ~ Ca_ppm, data=map.div)
ggplot(Ca_rich, aes(x = Ca_ppm, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(Ca_rich)
ggscatter(map.div, x = "Ca_ppm", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Ca (ppm)", ylab = "Richness")+
          annotate("text", x=2500, y=5300, label = "y == -0.18(x)+4891.12", parse=T)+
          annotate("text", x=2500, y=5200, label = "R^2 == 0.148", parse=T)+
          annotate("text", x=2500, y=5100, label = "p-val == 0.008", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

#. 5. BACTERIAL RICHNESS - P
P_rich <- lm(Richness ~ P_ppm, data=map.div)
ggplot(P_rich, aes(x = P_ppm, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(P_rich)
ggscatter(map.div, x = "P_ppm", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "P (ppm)", ylab = "Richness")+
          annotate("text", x=45, y=5300, label = "y == 2.77(x)+4386.16", parse=T)+
          annotate("text", x=35, y=5200, label = "R^2 ==  0.218", parse=T)+
          annotate("text", x=35, y=5100, label = "p-val == 0.001", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

#. 6. BACTERIAL RICHNESS - OM
OM_rich <- lm(Richness ~ OM_percent, data=map.div)
ggplot(OM_rich, aes(x = OM_percent, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(OM_rich)
ggscatter(map.div, x = "OM_percent", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "OM (%)", ylab = "Richness")+
          annotate("text", x=2.3, y=4200, label = "y == -167.6(x)+5132.34", parse=T)+
          annotate("text", x=2, y=4100, label = "R^2 ==  0.178", parse=T)+
          annotate("text", x=2, y=4000, label = "p-val == 0.003", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())


# 1. BACTERIAL SHANNON - SAND
Sand_percent_sha <- lm(Shannon ~ sand_percent, data=map.div)
ggplot(Sand_percent_sha, aes(x = sand_percent, y = Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(Sand_percent_sha)
ggscatter(map.div, x = "sand_percent", y = "Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Sand (%)", ylab = "Shannon")+
          annotate("text", x=40, y=7.7, label = "y == 0.004(x)+7.19", parse=T)+
          annotate("text", x=40, y=7.67, label = "R^2 == 0.157", parse=T)+
          annotate("text", x=40, y=7.63, label = "p-val == 0.006", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 2. BACTERIAL SHANNON - SILT
Silt_percent_sha <- lm(Shannon ~ silt_percent, data=map.div)
ggplot(Silt_percent_sha, aes(x = silt_percent, y = Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(Silt_percent_sha)
ggscatter(map.div, x = "silt_percent", y = "Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Silt (%)", ylab = "Shannon")+
          annotate("text", x=35, y=7.7, label = "y == -0.004(x)+7.55", parse=T)+
          annotate("text", x=35, y=7.67, label = "R^2 == 0.089", parse=T)+
          annotate("text", x=35, y=7.63, label = "p-val == 0.04", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())
    
# 3. BACTERIAL SHANNON - CLAY
Clay_percent_sha <- lm(Shannon ~ clay_percent, data=map.div)
summary(Clay_percent_sha)
ggplot(Clay_percent_sha, aes(x = clay_percent, y = Shannon)) + geom_point() + stat_smooth(method = 'lm')
ggscatter(map.div, x = "clay_percent", y = "Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Clay (%)", ylab = "Shannon")+
          annotate("text", x=20, y=7.7, label = "y == -0.01(x)+7.66", parse=T)+
          annotate("text", x=20, y=7.67, label = "R^2 == 0.259", parse=T)+
          annotate("text", x=20, y=7.63, label = "p-val == 0.0003", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())
    
#. 4. BACTERIAL SHANNON - K
K_sha <- lm(Shannon ~ K_ppm, data=map.div)
ggplot(K_sha, aes(x = K_ppm, y = Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(K_sha)
ggscatter(map.div, x = "K_ppm", y = "Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "K (ppm)", ylab = "Shannon")+
          annotate("text", x=200, y=7.7, label = "y == -0.0007(x)+7.53", parse=T)+
          annotate("text", x=200, y=7.67, label = "R^2 == 0.089", parse=T)+
          annotate("text", x=200, y=7.63, label = "p-val == 0.04", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

#. 5. BACTERIAL SHANNON - P
P_sha <- lm(Shannon ~ P_ppm, data=map.div)
ggplot(P_sha, aes(x = P_ppm, y = Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(P_sha)
ggscatter(map.div, x = "P_ppm", y = "Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "P (ppm)", ylab = "Shannon")+
          annotate("text", x=35, y=7.7, label = "y == 0.001(x)+7.30", parse=T)+
          annotate("text", x=35, y=7.67, label = "R^2 ==  0.236", parse=T)+
          annotate("text", x=35, y=7.63, label = "p-val == 0.0007", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())


### FUNGI ###

# 1. FUNGAL SHANNON - SAND
fg.sha_sand <- lm(fg.Shannon ~ sand_percent, data=map.div)
ggplot(fg.sha_sand, aes(x = sand_percent, y = fg.Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(fg.sha_sand)
ggscatter(map.div, x = "sand_percent", y = "fg.Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Sand (%)", ylab = "Shannon")+
          annotate("text", x=38, y=4.5, label = "y == 0.01(x)+2.96", parse=T)+
          annotate("text", x=37, y=4.4, label = "R^2 == 0.141", parse=T)+
          annotate("text", x=37, y=4.3, label = "p-val == 0.01", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 2. FUNGAL SHANNON - SILT
fg.sha_silt <- lm(fg.Shannon ~ silt_percent, data=map.div)
ggplot(fg.sha_silt, aes(x = silt_percent, y = fg.Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(fg.sha_silt)
ggscatter(map.div, x = "silt_percent", y = "fg.Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Silt (%)", ylab = "Shannon")+
          annotate("text", x=35, y=4.5, label = "y == -0.01(x)+4.15", parse=T)+
          annotate("text", x=35, y=4.4, label = "R^2 == 0.132", parse=T)+
          annotate("text", x=35, y=4.3, label = "p-val == 0.01", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 3. FUNGAL SHANNON - CLAY
fg.sha_clay <- lm(fg.Shannon ~ clay_percent, data=map.div)
ggplot(fg.sha_clay, aes(x = clay_percent, y = fg.Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(fg.sha_clay)
ggscatter(map.div, x = "clay_percent", y = "fg.Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Clay (%)", ylab = "Shannon")+
          annotate("text", x=22, y=4.6, label = "y == -0.03(x)+4.12", parse=T)+
          annotate("text", x=22, y=4.5, label = "R^2 == 0.102", parse=T)+
          annotate("text", x=22, y=4.4, label = "p-val == 0.03", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

######################################################################################################
#ANOVA TEST FOR NEMATODES, OLIGOCHAETES, AND MYCORRHIZAL FUNGI TOTAL ABSOLUTE ABUNDANCE (COUNT DATA)##
######################################################################################################

# 1. Compare total absolute nematodes abundances among sites
map.div$Site<-as.factor(map.div$Site)
TC_nema_site <- lm(nema.total.SqrtCount ~ Site, data=map.div, na.action=na.exclude)
summary(TC_nema_site)
drop1(TC_nema_site,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
nema_site_resids <- residuals(TC_nema_site)
nema_site_preds <- predict(TC_nema_site)
plot(nema_site_resids ~ nema_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(nema_site_resids) # p-value = 0.056, data errors are normally distributed 
skew_xts <- skewness(nema_site_resids)
#Perform Levene's Test for homogenity of variances
leveneTest(nema.total.SqrtCount ~ Site, data=map.div, na.action=na.exclude) # variances among group are homogenous
# Plotting
boxplot(nema.total.SqrtCount ~ Site, data = map.div.nema) # there are no outliers
plot(density(nema_site_resids))
qqnorm(nema_site_resids)
qqline(nema_site_resids) 
hist(nema_site_resids)
# Tukey
HSD_TC_nema_site <- HSD.test(TC_nema_site, "Site", group = TRUE, console = TRUE)
HSD_TC_nema_site <- HSD.test(TC_nema_site, "Site", group = F, console = TRUE)
# add significance letters from HSD.test into box plot
nema.count_site.summarized <- map.div %>% group_by(Site) %>% summarize(max.nema.count=max(nema.total.Count))
HSD_TC_nema_site <- HSD.test(TC_nema_site, "Site", group = TRUE, console = TRUE)
HSD_TC_nema_site.group = HSD_TC_nema_site$groups
class(HSD_TC_nema_site.group)
HSD_TC_nema_site.group$Site <- rownames(HSD_TC_nema_site.group)
new.nema.count_site.summarized=left_join(HSD_TC_nema_site.group,nema.count_site.summarized, by='Site')  
# Do plot
nema.ab.site <- ggplot(map.div, aes(x=Site, y=nema.total.Count))+
 geom_boxplot() +
 geom_point() +
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 theme_bw()+
 labs(title="A", y=expression(bold(atop("", 
               atop(textstyle("Nematode"),
                    atop(textstyle("(per 100 g soil)")))))))+
 #geom_text(data=new.nema.count_site.summarized,aes(x=Site,y=20+max.nema.count,label=groups),vjust=0)+
 theme(axis.text.x=element_blank(), 
       axis.ticks.x = element_blank(),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_blank(),
       axis.title=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid = element_blank())

# 2. Compare total absolute nematodes abundances among rootstocks
TC_nema_root <- lm(nema.total.SqrtCount ~ Rootstock, data=map.div, na.action=na.exclude)
summary(TC_nema_root)
drop1(TC_nema_root,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
nema_root_resids <- residuals(TC_nema_root)
nema_root_preds <- predict(TC_nema_root)
plot(nema_root_resids ~ nema_root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(nema_root_resids) # p-value = 0.12, data errors are normally distributed 
skew_xts <-  skewness(nema_root_resids)
#Perform Levene's Test for homogenity of variances
leveneTest(nema.total.SqrtCount ~ Rootstock, data=map.div, na.action=na.exclude) # variances among group are homogenous
# Tukey
HSD_TC_nema_root <- HSD.test(TC_nema_root, "Rootstock", group = TRUE, console = TRUE)
HSD_TC_nema_root <- HSD.test(TC_nema_root, "Rootstock", group = F, console = TRUE)
# add significance letters from HSD.test into box plot
TC_nema_root.sum1 <- map.div %>% group_by(Rootstock) %>% summarize(max.count=max(nema.total.Count))
HSD_TC_nema_root <- HSD.test(TC_nema_root, "Rootstock", group = TRUE, console = TRUE)
hsd_TC_nema_root = HSD_TC_nema_root$groups
class(hsd_TC_nema_root)
hsd_TC_nema_root$Rootstock <- rownames(hsd_TC_nema_root)
TC_nema_root.sum2=left_join(hsd_TC_nema_root,TC_nema_root.sum1, by='Rootstock')  
# Do Plot
nema.ab.root <- ggplot(map.div, aes(x=Rootstock, y=nema.total.Count))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
 labs(title = "D", y=expression(bold(atop("", 
               atop(textstyle("Nematode"), 
               atop(textstyle("(per 100 g soil)")))))))+
 geom_text(data=TC_nema_root.sum2,aes(x=Rootstock,y=45+max.count,label=TC_nema_root.sum2$groups),vjust=0)+
 theme(axis.text.x=element_blank(), 
       axis.ticks.x = element_blank(),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_blank(),
       axis.title=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid = element_blank())

# 3. Compare total absolute nematodes abundances among cultivars # not significant
TC_nema_cult <- lm(nema.total.SqrtCount ~ cultivar, data=map.div, na.action=na.exclude)
summary(TC_nema_cult)
drop1(TC_nema_cult,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
nema_cult_resids <- residuals(TC_nema_cult)
nema_cult_preds <- predict(TC_nema_cult)
plot(nema_cult_resids ~ nema_cult_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(nema_cult_resids) # p-value = 0.13, data errors are normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(nema.total.SqrtCount ~ cultivar, data=map.div, na.action=na.exclude) # variances among group are homogenous

# 4. Compare absolute abundance of oligochaetes among sites
OL_site <- lm(sqrt.OL ~ Site, data=map.div, na.action=na.exclude)
summary(OL_site)
drop1(OL_site,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
OL_site_resids <- residuals(OL_site)
OL_site_preds <- predict(OL_site)
plot(OL_site_resids ~ OL_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(OL_site_resids) # p-value = 0.48, data errors are normally distributed 
skew_xts <-  skewness(OL_site_resids)
#Perform Levene's Test for homogenity of variances
leveneTest(sqrt.OL ~ Site, data=map.div, na.action=na.exclude) # variances among group are homogenous
# Tukey
HSD_OL_site <- HSD.test(OL_site, "Site", group = TRUE, console = TRUE)
HSD_OL_site <- HSD.test(OL_site, "Site", group = F, console = TRUE)
# add significance letters from HSD.test into box plot
OL_site.sum1 <- map.div %>% group_by(Site) %>% summarize(max.count=max(Oligochaetes))
HSD_OL_site <- HSD.test(OL_site, "Site", group = TRUE, console = TRUE)
hsd_OL_site = HSD_OL_site$groups
class(hsd_OL_site)
hsd_OL_site$Site <- rownames(hsd_OL_site)
OL_site.sum2=left_join(hsd_OL_site,OL_site.sum1, by='Site')  
# Do Plot
OL_site <- ggplot(map.div, aes(x=Site, y=Oligochaetes))+
 geom_boxplot() +
 geom_point() +
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 theme_bw()+
 labs(title = "B", y=expression(bold(atop("", 
               atop(textstyle("Oligochaetes"),
                    atop(textstyle("(per 100 g soil)")))))))+
 geom_text(data=OL_site.sum2,aes(x=Site,y=2+max.count,label=OL_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_blank(), 
       axis.ticks.x = element_blank(),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_blank(),
       axis.title=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid = element_blank())

# 5. Compare absolute abundance of oligochaetes among rootstocks 
OL.root <- lm(sqrt.OL ~ Rootstock, data=map.div, na.action=na.exclude)
summary(OL.root)
drop1(OL.root,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
OL_root_resids <- residuals(OL.root)
OL_root_preds <- predict(OL.root)
plot(OL_root_resids ~ OL_root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(OL_root_resids) # p-value = 0.36, data errors are normally distributed 
skew_xts <-  skewness(OL_site_resids)
#Perform Levene's Test for homogenity of variances
leveneTest(sqrt.OL ~ Rootstock, data=map.div, na.action=na.exclude) # variances among group are homogenous
# Plotting
boxplot(sqrt.OL ~ Rootstock, data = map.div.nema) # there are no outliers
plot(density(OL_site_resids))
qqnorm(OL_site_resids)
qqline(OL_site_resids) 
hist(OL_site_resids)
# Do Plot
OL_root <- ggplot(map.div, aes(x=Rootstock, y=Oligochaetes))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
 labs(title = "E", y=expression(bold(atop("", 
               atop(textstyle("Oligochaetes"), 
               atop(textstyle("(per 100 g soil)")))))))+
#geom_text(data=OL_site.sum2,aes(x=Site,y=2+max.count,label=OL_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_blank(), 
       axis.ticks.x = element_blank(),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_blank(),
       axis.title=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid = element_blank())

# 6. Compare absolute abundance of oligochaetes among cultivars # not significant
OL.cul <- lm(sqrt.OL ~ cultivar, data=map.div, na.action=na.exclude)
summary(OL.cul)
drop1(OL.cul,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
OL_cul_resids <- residuals(OL.cul)
OL_cul_preds <- predict(OL.cul)
plot(OL_cul_resids ~ OL_cul_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(OL_cul_resids) # p-value = 0.30, data errors are normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(sqrt.OL ~ cultivar, data=map.div, na.action=na.exclude)

# 7. Compare absolute mycorrhizal fungi abundance among sites
MF_site <- lm(log.MF ~ Site, data=map.div, na.action=na.exclude)
summary(MF_site)
drop1(MF_site,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
MF_site_resids <- residuals(MF_site)
MF_site_preds <- predict(MF_site)
plot(MF_site_resids ~ MF_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(MF_site_resids) # p-value = 0.4, data errors are normally distributed 
skew_xts <-  skewness(MF_site_resids)
#Perform Levene's Test for homogenity of variances
leveneTest(log.MF ~ Site, data=map.div, na.action=na.exclude) # variances among group are homogenous
# Tukey
HSD_MF_site <- HSD.test(MF_site, "Site", group = TRUE, console = TRUE)
HSD_MF_site <- HSD.test(MF_site, "Site", group = F, console = TRUE)
# add significance letters from HSD.test into box plot
MF_site.sum1 <- map.div %>% group_by(Site) %>% summarize(max.count=max(MycorrhizalFungi))
HSD_MF_site <- HSD.test(MF_site, "Site", group = TRUE, console = TRUE)
hsd_MF_site = HSD_MF_site$groups
class(hsd_MF_site)
hsd_MF_site$Site <- rownames(hsd_MF_site)
MF_site.sum2=left_join(hsd_MF_site,MF_site.sum1, by='Site')  
# Do Plot
MF_site <- ggplot(map.div, aes(x=Site, y=MycorrhizalFungi))+
 geom_boxplot() +
 geom_point() +
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 theme_bw()+
 labs(title="C",y=expression(bold(atop("", 
               atop(textstyle("Mycorrhizal Fungi"), 
               atop(textstyle("(per 100 g soil)")))))))+
 geom_text(data=MF_site.sum2,aes(x=Site,y=10+max.count,label=MF_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=10,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_text(),
       axis.title=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 8. Compare absolute mycorrhizal fungi abundance among rootstocks
MF_root <- lm(log.MF ~ Rootstock, data=map.div, na.action=na.exclude)
summary(MF_root)
drop1(MF_root,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
MF_root_resids <- residuals(MF_root)
MF_root_preds <- predict(MF_root)
plot(MF_root_resids ~ MF_root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(MF_root_resids) # p-value = 0.7, data errors are normally distributed 
skew_xts <-  skewness(MF_site_resids)
#Perform Levene's Test for homogenity of variances
leveneTest(log.MF ~ Rootstock, data=map.div, na.action=na.exclude) # variances among group are homogenous
# Do Plot
# Do Plot
MF_root <- ggplot(map.div, aes(x=Rootstock, y=MycorrhizalFungi))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
 labs(title="F",y=expression(bold(atop("", 
               atop(textstyle("Mycorrhizal Fungi"), 
               atop(textstyle("(per 100 g soil)")))))))+
#geom_text(data=OL_site.sum2,aes(x=Site,y=2+max.count,label=OL_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_text(size=15,angle=49,hjust =0.9), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_text(vjust=13),
       axis.title=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 9. Compare absolute mycorrhizal fungi abundance among cultivars # not significant
MF_cul <- lm(log.MF ~ cultivar, data=map.div, na.action=na.exclude)
summary(MF_cul)
drop1(MF_cul,~.,test="F") # type III SS and F Tests
#TESTING ASSUMPTIONS
#Generate residual and predicted values
MF_cul_resids <- residuals(MF_cul)
MF_cul_preds <- predict(MF_cul)
plot(MF_cul_resids ~ MF_cul_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(MF_cul_resids) # p-value = 0.3, data errors are normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(log.MF ~ cultivar, data=map.div, na.action=na.exclude)

# 10. Arrange Plot
grid.newpage()
grid.draw(rbind(ggplotGrob(nema.ab.site), ggplotGrob(OL_site), ggplotGrob(MF_site),size = "first"))
grid.newpage()
grid.draw(rbind(ggplotGrob(nema.ab.root), ggplotGrob(OL_root), ggplotGrob(MF_root),size = "first"))

# 3. ANOVA TEST FOR NEMATODE ALPHA DIVERSITY
 
# 1. compare nematode richness among sites #not significant
AovNema_richness_site <- lm(nema.Richness ~ Site , data=map.div, na.action=na.exclude)
summary(AovNema_richness_site)
drop1(AovNema_richness_site,~.,test="F") # type III SS and F Tests
# Check normality
Richness.Nema_resids <- residuals(AovNema_richness_site)
Richness.Nema_preds <- predict(AovNema_richness_site)
plot(Richness.Nema_resids ~ Richness.Nema_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
shapiro.test(Richness.Nema_resids) # not normal
qqnorm(Richness.Nema_resids)
qqline(Richness.Nema_resids) 
leveneTest(nema.Richness ~ Site, data=map.div, na.action=na.exclude) # not homogen
gvlma(AovNema_richness_site)

# 2. compare nematode shannon index among sites #not significant
AovNema_shannon_site <- lm(nema.Shannon ~ Site, data=map.div, na.action=na.exclude)
drop1(AovNema_shannon_site,~.,test="F") # type III SS and F Tests
summary(AovNema_shannon_site)
# Check normality
Sha.Nema_resids <- residuals(AovNema_shannon_site)
Sha.Nema_preds <- predict(AovNema_shannon_site)
plot(Sha.Nema_resids ~ Sha.Nema_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
shapiro.test(Sha.Nema_resids) # normal
qqnorm(Sha.Nema_resids)
qqline(Sha.Nema_resids) 
leveneTest(nema.Shannon ~ Site, data=map.div, na.action=na.exclude) # homogen
gvlma(AovNema_shannon_site)

# 3. compare nematode richness among rootstocks #not significant
AovNema_richness_root <- lm(nema.Richness ~ Rootstock , data=map.div, na.action=na.exclude)
summary(AovNema_richness_root)
drop1(AovNema_richness_root,~.,test="F") # type III SS and F Tests
# Check normality
Richness.Nema_resids <- residuals(AovNema_richness_root)
Richness.Nema_preds <- predict(AovNema_richness_root)
plot(Richness.Nema_resids ~ Richness.Nema_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
shapiro.test(Richness.Nema_resids) #normal
qqnorm(Richness.Nema_resids)
qqline(Richness.Nema_resids) 
leveneTest(nema.Richness ~ Rootstock, data=map.div, na.action=na.exclude) # homogen
gvlma(AovNema_richness_root)

# 4. compare nematode shannon index among rootstocks #not significant
AovNema_shannon_root <- lm(nema.Shannon ~ Rootstock , data=map.div, na.action=na.exclude)
summary(AovNema_shannon_root)
drop1(AovNema_shannon_root,~.,test="F") # type III SS and F Tests
# Check normality
Richness.Nema_resids <- residuals(AovNema_shannon_root)
Richness.Nema_preds <- predict(AovNema_shannon_root)
plot(Richness.Nema_resids ~ Richness.Nema_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
shapiro.test(Richness.Nema_resids) # normal
qqnorm(Richness.Nema_resids)
qqline(Richness.Nema_resids) 
leveneTest(nema.Shannon ~ Rootstock, data=map.div, na.action=na.exclude) # homogen
gvlma(AovNema_shannon_root)

# ANOVA test to check the differences of nematoda alpha diversity among cultivars
# 1. nematode - Richness
nema_rich_cul <- lm(nema.Richness ~ cultivar, data=map.div, na.action=na.exclude)
summary(nema_rich_cul) # not significant
# 2. nematode - Shannon
nema_sha_cul <- lm(nema.Shannon ~ cultivar, data=map.div, na.action=na.exclude)
summary(nema_sha_cul) # not significant

# Check nematode prevalence
nema_PA <- 1*(nema.t>0)
sum_nemaPA <- rowSums(nema_PA)
sum_nemaPA

###PEARSON CORRELATION AND LINEAR REGRESSION TEST OF NEMATODES, OLIGOCHAETES, MYCORRHYZAL FUNGI ABSOLUTE ABUNDANCES TO BACTERIAL AND FUNGAL RICHNESS AND SHANNON ###
# 1. Bacterial Richness and nematode Tylenchs
bac.rich_tyl <- cor.test(map.div$Tylenchs, map.div$Richness, method ="pearson")
cor.test(map.div$Tylenchs, map.div$Richness)
# linear regression
lm.bac.rich_tyl <- lm(Richness ~ Tylenchs, data = map.div) # OK
ggplot(lm.bac.rich_tyl, aes(x = Tylenchs, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(lm.bac.rich_tyl)
gvlma(lm.bac.rich_tyl)
# Do plot
TY_rich <- lm(Richness ~ Tylenchs, data=map.div)
summary(TY_rich)
ggplot(TY_rich, aes(x = Tylenchs, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
ggscatter(map.div, x = "Tylenchs", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Tylenchs (individuals per 100 g soil)", ylab = "Richness")+
          annotate("text", x=180, y=5200, label = "y == -1.55(x)+4753.92", parse=T)+
          annotate("text", x=180, y=5100, label = "R^2 ==  0.108", parse=T)+
          annotate("text", x=180, y=5000, label = "p-val == 0.026", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=15,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=15,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 2. Bacterial Shannon and nematode Tylenchs
bac.sha_tyl <- cor.test(map.div$Tylenchs, map.div$Shannon, method ="pearson")
cor.test(map.div$Tylenchs, map.div$Shannon)
# linear regression
lm.bac.sha_tyl <- lm(Shannon ~ Tylenchs, data = map.div) # OK
ggplot(lm.bac.sha_tyl, aes(x = Tylenchs, y = Shannon)) + geom_point() + stat_smooth(method = 'lm')
summary(lm.bac.sha_tyl)
gvlma(lm.bac.sha_tyl)
# Do plot
TY_sha <- lm(Shannon ~ Tylenchs, data=map.div)
summary(TY_sha)
ggplot(TY_sha, aes(x = Tylenchs, y = Shannon)) + geom_point() + stat_smooth(method = 'lm')
ggscatter(map.div, x = "Tylenchs", y = "Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Tylenchs (individuals per 100 g soil)", ylab = "Shannon")+
          annotate("text", x=190, y=7.65, label = "y == -0.0007(x)+7.46", parse=T)+
          annotate("text", x=180, y=7.62, label = "R^2 ==  0.143", parse=T)+
          annotate("text", x=180, y=7.58, label = "p-val == 0.01", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=15,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 3. Fungal Richness and nematode Tylenchs # not significant
fg.rich_tyl <- cor.test(map.div$Tylenchs, map.div$fg.Richness, method ="pearson")
cor.test(map.div$Tylenchs, map.div$fg.Richness)
# 4. Fungal Shannon and nematode Tylenchs # not significant
fg.sha_tyl <- cor.test(map.div$Tylenchs, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$Tylenchs, map.div$fg.Shannon)

# 3. Bacterial Shannon and nematode Shannon
BC.sha_NM.sha<- lm(map.div$Shannon ~ map.div$nema.Shannon)
summary(BC.sha_NM.sha)
ggscatter(map.div, x = "nema.Shannon", y = "Shannon", 
          add = "reg.line",conf.int = TRUE, xlab = "Nematode Shannon index", ylab = "Bacterial/archaeal Shannon index")+
          annotate("text", x=1.25, y=7.65, label = "y == -0.12(x)+7.51", parse=T)+
          annotate("text", x=1.25, y=7.62, label = "R^2 == 0.108", parse=T)+
          annotate("text", x=1.25, y=7.58, label = "p-val == 0.027", parse=T)+
     geom_smooth(method='lm')+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=15,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

############################################################################################
######################### BACTERIAL, FUNGAL, and NEMATODE BETA DIVERSITY ###################
############################################################################################

# 1. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR BACTERIA
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist <- vegdist(t(otu), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
env <- map[,c(11:22, 24:36)]
# scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1]
ax2.scores=otu_pcoa$points[,2] 
env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1 <- otu_pcoa$eig[1]/sum(otu_pcoa$eig)
ax2 <- otu_pcoa$eig[2]/sum(otu_pcoa$eig)
map2=cbind(map.div,ax1.scores,ax2.scores)
# simple plot
pcoa_plot <- plot(ax1.scores, ax2.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
plot(env_fit, p.max=0.05, col="red1")
# Do better plot
A <- as.list(env_fit$vectors) #shortcutting ef$vectors
pvals<-as.data.frame(A$pvals) #creating the dataframe
#environment scores (vectors scaled by R2 values)
bac.scores1 <- as.data.frame(scores(env_fit, display="vectors"))
bac.scores2 <- cbind(bac.scores1, pvals)
bac.scores3 <- cbind(bac.scores2,Variable=rownames(bac.scores2))
bac.scores4 <- subset(bac.scores3,pvals<0.05)
library(ggrepel)
mult <-.25
# Plot by site
bac.pcoa <- ggplot(data = map2, aes(x=ax1.scores, y=ax2.scores))+
  theme_bw()+
 #geom_point(data = map2, aes(x = ax1.scores, y = ax2.scores, color=Site),size=2.5, shape=20,stroke=1.75)+
  geom_text(aes(ax1.scores, y=ax2.scores,color=Site, label = Site),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""),limits = c(-0.23,0.4))+
  geom_segment(data=bac.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = bac.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
        coord_fixed() + 
 labs(title = "A")+
 theme(legend.position="none",
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=12),
       legend.title = element_text(size = 12),
       legend.spacing.x = unit(0.05, 'cm'))+scale_color_discrete(breaks=sort(as.numeric(map2$Site)))

# Plot by rootstock
bac_root.pcoa <- ggplot(data = map2, aes(x=ax1.scores, y=ax2.scores))+
  theme_bw()+
  #geom_point(aes(ax1.scores, y=ax2.scores,color=cultivar),size=1.5,shape=20,stroke=1.75)+
 geom_text(aes(ax1.scores, y=ax2.scores,color=Rootstock, label = Rootstock),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""), limits = c(-0.23,0.4))+
  geom_segment(data=bac.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = bac.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
        coord_fixed() + 
 labs(title = "A")+
 theme(legend.position="none",
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))

# 2. CALCULATE THE BETA DIVERSITY (PCoA PLOT) FOR FUNGI
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_distITS <- vegdist(t(otuITS), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoaITS <- cmdscale(otu_distITS, eig=T)
env <- map[,c(11:22, 24:36)]
# scores of PC1 and PC2
ax1ITS.scores=otu_pcoaITS$points[,1] 
ax2ITS.scores=otu_pcoaITS$points[,2] 
env_fitITS <- envfit(otu_pcoaITS, env, na.rm=TRUE)
ax1ITS <- otu_pcoaITS$eig[1]/sum(otu_pcoaITS$eig)
ax2ITS <- otu_pcoaITS$eig[2]/sum(otu_pcoaITS$eig)
map2=cbind(map2,ax1ITS.scores,ax2ITS.scores)
# simple plot
fg.pcoa_plot <- plot(ax1ITS.scores, ax2ITS.scores, xlab=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))
plot(env_fitITS, p.max=0.05, col="red1")
# Do better plot
#shortcutting ef$vectors
A.its <- as.list(env_fitITS$vectors)
#creating the dataframe
pvals.its <-as.data.frame(A.its$pvals)
# environment scores (vectors scaled by R2 values)
fg.scores1 <- as.data.frame(scores(env_fitITS, display="vectors"))
fg.scores2 <- cbind(fg.scores1, pvals.its)
fg.scores3 <- cbind(fg.scores2,Variable=rownames(fg.scores2))
fg.scores4 <- subset(fg.scores3,pvals.its<0.05)
library(ggrepel)
mult <-.25
# Plot by site
fg.pcoa <- ggplot(data = map2, aes(x=ax1ITS.scores, y=ax2ITS.scores)) +
  theme_bw()+
  #geom_point(aes(ax1ITS.scores, y=ax2ITS.scores,color=Site),size=2.5, shape=20,stroke=1.75)+
 geom_text(aes(ax1ITS.scores, y=ax2ITS.scores,color=Site, label = Site),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  coord_fixed()+
 labs(title = "B")+
  geom_segment(data=fg.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = fg.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
 theme(plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       legend.position="none",
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=12),
       legend.title = element_text(size = 12),
       legend.spacing.x = unit(0.05, 'cm'))+scale_color_discrete(breaks=sort(as.numeric(map2$Site)))

# Plot by rootstocks
fg_root.pcoa <- ggplot(data = map2, aes(x=ax1ITS.scores, y=ax2ITS.scores)) +
  theme_bw()+
  #geom_point(aes(ax1ITS.scores, y=ax2ITS.scores,color=cultivar),size=1.5,shape=20,stroke=1.75)+
 geom_text(aes(ax1ITS.scores, y=ax2ITS.scores,color=Rootstock, label = Rootstock),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  coord_fixed()+
 labs(title = "B")+
  geom_segment(data=fg.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = fg.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
 theme(plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       legend.position="none",
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=14),
       legend.title = element_text(size = 15),
       legend.spacing.x = unit(0.05, 'cm'))
#ggarrange(bac_root.pcoa, fg_root.pcoa, nrow=1, ncol=2, common.legend = TRUE, legend="left")

# Arrange plots
grid.newpage()
grid.draw(cbind(ggplotGrob(bac.pcoa), ggplotGrob(fg.pcoa),size = "first"))
grid.newpage()
grid.draw(cbind(ggplotGrob(bac_root.pcoa), ggplotGrob(fg_root.pcoa),size = "first"))

# 3. CALCULATE THE BETA DIVERSITY (PCoA PLOT) FOR NEMATODE
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist_nema <- vegdist(t(sqrt.nema.t), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa_nema <- cmdscale(otu_dist_nema, eig=T)
env.nema <- map[,c(11:22, 35:36)]
# scores of PC1 and PC2
ax1.nema.scores=otu_pcoa_nema$points[,1]
ax2.nema.scores=otu_pcoa_nema$points[,2] 
env_fit.nema <- envfit(otu_pcoa_nema, env.nema, na.rm=TRUE)
# calculate percent variance explained, then add to plot
ax1.nema <- otu_pcoa_nema$eig[1]/sum(otu_pcoa_nema$eig)
ax2.nema <- otu_pcoa_nema$eig[2]/sum(otu_pcoa_nema$eig)
map2=cbind(map2,ax1.nema.scores,ax2.nema.scores)
# simple plot
pcoa_plot <- plot(ax1.nema.scores, ax2.nema.scores, xlab=paste("PCoA1: ",round(ax1.nema,3)*100,"% var. explained", sep=""), ylab=paste("PCoA2: ",round(ax2.nema,3)*100,"% var. explained", sep=""))
plot(env_fit.nema, p.max=0.05, col="red1")
# Do better plot
#shortcutting ef$vectors
A.nema <- as.list(env_fit.nema$vectors)
#creating the dataframe
pvals.nema<-as.data.frame(A.nema$pvals)
# environment scores (vectors scaled by R2 values)
nm.scores1 <- as.data.frame(scores(env_fit.nema, display="vectors"))
nm.scores2 <- cbind(nm.scores1, pvals.nema)
nm.scores3 <- cbind(nm.scores2,Variable=rownames(nm.scores2))
nm.scores4 <- subset(nm.scores3,pvals.nema<0.05)
library(ggrepel)
mult <-.25
#Plot by site
nm.pcoa <- ggplot(data = map2, aes(x=ax1.nema.scores, y=ax2.nema.scores))+
  theme_bw()+
  #geom_point(aes(ax1.nema.scores, y=ax2.nema.scores,color=Site),size=2.5,shape=20,stroke=1.75)+
 geom_text(aes(ax1.nema.scores, y=ax2.nema.scores,color=Site, label = Site),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1.nema,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2.nema,3)*100,"% var. explained", sep=""))+
  geom_segment(data=nm.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = nm.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
        coord_fixed() + 
 labs(title = "C")+
 theme(legend.position = "none",
       legend.box = "vertical",
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=8),
       legend.title = element_text(size = 12),
       legend.spacing.x = unit(0.05, 'cm'))+scale_color_discrete(guide = guide_legend(title.position = "top", nrow = 2),breaks=sort(as.numeric(map2$Site)))

################# PERMANOVA TO TEST ANY INFLUENCE OF SITE, ROOTSTOCK, AND CULTIVAR TO BACTERIAL, FUNGAL, AND NEMATODE BETA DIVERSITY ################################
### SITE ###
# 1. Site to bacterial beta diversity
bac.otudist_site=adonis(otu_dist~map2$Site) # p-val=0.001
# 2. Site to fungal beta diversity
fg.otudist_site=adonis(otu_distITS~map2$Site) # p-val=0.001
# 3. Site to nematode beta diversity
nm.otudist_site=adonis(otu_dist_nema~map2$Site) # p-val=0.018
### ROOTSTOCKS ###
# 1. Rootstocks to bacterial beta diversity
bac.otudist_root=adonis(otu_dist~map2$rootstock) # p-val=0.002
# 2. Rootstocks to fungal beta diversity
fg.otudist_root=adonis(otu_distITS~map2$rootstock) # p-val=0.01
# 3. Rootstocks to nematode beta diversity
nm.otudist_root=adonis(otu_dist_nema~map2$rootstock) # p-val=0.095 # not significant
### CULTIVARS ###
# 1. Cultivars to bacterial beta diversity
bac.otudist_cul=adonis(otu_dist~map2$cultivar) # p-val=0.46 # not significant
# 2. Cultivars to fungal beta diversity
fg.otudist_cul=adonis(otu_distITS~map2$cultivar) # p-val=0.12 # not significant
# 3. Cultivars to nematode beta diversity
nm.otudist_cul=adonis(otu_dist_nema~map2$cultivar) # p-val=0.13 # not significant

############################################################################################
######################### BACTERIAL AND FUNGAL COMMUNITIES COMPOSITION #####################
############################################################################################
BiocManager::install("phyloseq")
library(phyloseq)

# 1. BACTERIA COMPOSITION
# read bacterial taxonomy
tax_16S = read.csv("16S_TAX.csv", sep=',', header=T)
tax_16S
rownames(tax_16S) <- rownames(otu)
# make phyloseq otu table and taxonomy
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax_16S))
# add map
map$Site<-as.factor(map$Site)
rownames(map) <- map$sample_code
# make phyloseq map
phyloseq_map <- sample_data(map)
# make phyloseq object
PHYL_16S <- merge_phyloseq(OTU,TAX,phyloseq_map)
PHYL_16S
# merge taxa by phylum
# 1. phylum - Bacteria
phylum.16S.ra <- transform_sample_counts(PHYL_16S, function(x) x/sum(x))
phylum.16S <- tax_glom(phylum.16S.ra, taxrank = "Phylum", NArm = F)
phylum.16S
cumabun=taxa_sums(phylum.16S)
df.phylum.16S.taxasum=as.data.frame(cumabun)
df.phylum.16S.taxasum$relabun=df.phylum.16S.taxasum$cumabun/45
sort.df.phylum.16S.taxasum=df.phylum.16S.taxasum[order(df.phylum.16S.taxasum$relabun, decreasing = T),]
view(sort.df.phylum.16S.taxasum)
write.table(sort.df.phylum.16S.taxasum, file = 'relabun_phylum_bac.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
df.phylum.16S <- psmelt(phylum.16S)
summary(df.phylum.16S)
head(df.phylum.16S)
dim(df.phylum.16S)
df.phylum.16S$Phylum <- as.character(df.phylum.16S$Phylum)
df.phylum.16S$Phylum[df.phylum.16S$Abundance < 0.01] <- "Other"
df.phylum.16S$Phylum[is.na(df.phylum.16S$Phylum)] <- "Other"
#colourCount = length(unique(df.phylum$Phylum))
#getPalette = colorRampPalette(brewer.pal(9, "Set1"))  
#barplot 16S
p <- ggplot(data=df.phylum.16S, aes(x=Sample, y=Abundance, fill=Phylum))
barplot.16S <- p + geom_bar(aes(), stat="identity", position="fill") + scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))+
 labs(title="A",y= "Relative Abundance")+
 theme(plot.title = element_text(size = rel(1.5), face="bold"),
       axis.text=element_text(size=9),
       axis.text.x = element_text(angle = 90,hjust = 1),
       axis.title=element_text(size=10,face="bold"),
       legend.text=element_text(size = 8),
       legend.title = element_text(size=10),
       panel.grid = element_blank(), 
       panel.background = element_blank(),
       panel.border = element_rect(colour = "black", fill = NA, size = 1))

# 2. FUNGAL COMPOSITION
# READ TAXONOMY FILE ("consensus_taxonomy_ITS.txt") FROM CONSTAX OUTPUT AND MERGE IT INTO "otuITS" FILE ("OTU_rarefied_ITS.txt")
otuITS <- read.table(file = "OTU_rarefied_ITS.txt", sep='\t', header=T)
fg.taxonomy = read.table("consensus_taxonomy_ITS.txt",
                      header = T, sep = "\t")
dim(fg.taxonomy)
fg.otu.tax <- merge(otuITS, fg.taxonomy, by.x =c("OTU_ID"), by.y = c("OTU_ID"))
head(fg.otu.tax)
dim(fg.otu.tax)
write.table(fg.otu.tax, file = "combined_otu_tax_ITS.txt", sep = '\t', col.names = TRUE, row.names = FALSE)
# go to excel and separate the "combined_otu_tax_ITS.txt" into 2 files: OTU_ITS.csv and ITS_TAX.csv for uploading to phyloseq
# READ OTU_TABLE
otuITS.phyl = read.csv("OTU_ITS.csv", sep=",", row.names=1)
otuITS.phyl = as.matrix(otuITS.phyl)
dim(otuITS.phyl)
head(sort(colSums(otuITS.phyl, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otuITS.phyl, na.rm = FALSE, dims = 1), decreasing = FALSE))
# read fungal taxonomy
tax.ITS = read.csv("ITS_TAX.csv", sep=",", row.names=1)
tax.ITS = as.matrix(tax.ITS)
dim(tax.ITS)
tax.ITS
# IMPORT fungal otu table, taxonomy, and map files into phyloseq object
OTU.ITS <- otu_table(otuITS.phyl, taxa_are_rows = TRUE)
OTU.ITS
TAX.ITS <-  tax_table(tax.ITS)
TAX.ITS
# check that your OTU names are consistent across objects
taxa_names(TAX.ITS) 
taxa_names(OTU.ITS)
# merge into one phyloseq object
PHYL_ITS <- phyloseq(OTU.ITS,TAX.ITS,phyloseq_map)
# merge taxa by phylum
# 1. phylum - Fungi
PHYL_ITS.ra <- transform_sample_counts(PHYL_ITS, function(x) x/sum(x))
phylum.ITS <- tax_glom(PHYL_ITS.ra, taxrank = "Phylum", NArm = F)
phylum.ITS

cumabun.its=taxa_sums(phylum.ITS)
df.phylum.ITS.taxasum=as.data.frame(cumabun.its)
df.phylum.ITS.taxasum$relabun=df.phylum.ITS.taxasum$cumabun.its/45
sort.df.phylum.ITS.taxasum=df.phylum.ITS.taxasum[order(df.phylum.ITS.taxasum$relabun, decreasing = T),]
view(sort.df.phylum.ITS.taxasum)
write.table(sort.df.phylum.ITS.taxasum, file = 'relabun_phylum_fg.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
df.phylum.ITS <- psmelt(phylum.ITS)
head(df.phylum.ITS)
dim(df.phylum.ITS)
df.phylum.ITS$Phylum <- as.character(df.phylum.ITS$Phylum)
df.phylum.ITS$Phylum[df.phylum.ITS$Abundance < 0.01] <- "Other"
df.phylum.ITS$Phylum[is.na(df.phylum.ITS$Phylum)] <- "Other"
#barplot ITS
f <- ggplot(data=df.phylum.ITS, aes(x=Sample, y=Abundance, fill=Phylum))
barplot.ITS <- f + geom_bar(aes(), stat="identity", position="fill") + scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))+
 labs(title="B", y= "Relative Abundance")+
 theme(plot.title = element_text(size = rel(1.5), face="bold"),
       axis.text=element_text(size=9),
       axis.text.x = element_text(angle = 90,hjust = 1),
       axis.title=element_text(size=10,face="bold"),
       legend.text=element_text(size = 8),
       legend.title = element_text(size=10),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(), 
       panel.background = element_blank(),
       panel.border = element_rect(colour = "black", fill = NA, size = 1))

# Arrange Plots
rid.newpage()
ggarrange(barplot.16S,barplot.ITS, ncol = 2, nrow = 1)




