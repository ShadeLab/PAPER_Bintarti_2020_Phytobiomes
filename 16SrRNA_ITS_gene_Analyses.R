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
# READ OTU FUNGI
otuITS <- read.table(file = "OTU_rarefied_ITS.txt", sep='\t', header=T, row.names = 1) 
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

##############################################################################################
#### ANOVA AND KRUSKAL-WALLIS TO COMPARE BACTERIAL AND FUNGAL RICHNESS & SHANNON AMONG SITES AND ROOTSTOCKS #######
##############################################################################################

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

# 3. ANOVA TEST FOR NEMATODES, OLIGOCHAETES, AND MYCORRHIZAL FUNGI TOTAL ABSOLUTE ABUNDANCE (COUNT DATA)

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


# 3. Compare absolute abundance of oligochaetes among sites
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


# 4. Compare absolute abundance of oligochaetes among sites 
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

# 5. Compare absolute mycorrhizal fungi abundance among sites
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
shapiro.test(MF_site_resids) # p-value = 0.3, data errors are normally distributed 
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

# 6. Compare absolute mycorrhizal fungi abundance among rootstocks
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
shapiro.test(MF_root_resids) # p-value = 0.3, data errors are normally distributed 
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

# 7. Arrange Plot
grid.newpage()
grid.draw(rbind(ggplotGrob(nema.ab.site), ggplotGrob(OL_site), ggplotGrob(MF_site),size = "first"))
grid.newpage()
grid.draw(rbind(ggplotGrob(nema.ab.root), ggplotGrob(OL_root), ggplotGrob(MF_root),size = "first"))

# ANOVA test to check the differences of total absolute abundances of nematoda, oligochaetes, and mycorrhizal fungi among cultivars
# 1. nematode
nema_cul <- lm(nema.total.SqrtCount ~ cultivar, data=map.div, na.action=na.exclude)
summary(nema_cul) # not significant
# 2. oligochaetes
ol_cul <- lm(sqrt.OL ~ cultivar, data=map.div, na.action=na.exclude)
summary(ol_cul) # not significant
# 3. mycorrhizal fungi
MF_cul <- lm(log.MF ~ cultivar, data=map.div, na.action=na.exclude)
summary(MF_cul) # not significant














nemagroup = read.csv("NematodeGroup.csv", sep=',', header=T)
nema.plus = as.data.frame(nema.plus)
nema.plus.tidy <- rownames_to_column(nema.plus,var = "Microorganism")
nema.plus.group=left_join(nema.plus.tidy,nemagroup, by = "Microorganism")






sqrt.allmic_relabund <- decostand(sqrt.nema.plus, method="total", MARGIN=2)
allmic.com_abund <- rowSums(sqrt.allmic_relabund)
df.allmic.com_abund <- as.data.frame(allmic.com_abund)
head(df.com_abund)
df.allmic.com_abund$RelAbund=df.allmic.com_abund$allmic.com_abund/45
sum(df.allmic.com_abund$allmic.com_abund)
sum(df.allmic.com_abund$RelAbund)
df.allmic.com_abund$PercentRelAbund=df.allmic.com_abund$RelAbund*100
sum(df.allmic.com_abund$PercentRelAbund)
df.allmic.com_abund=rownames_to_column(df.allmic.com_abund, var = "Microorganism")
head(df.allmic.com_abund)
dim(df.allmic.com_abund) 

df.sqrt.nema.plus <- as.data.frame(sqrt.nema.plus)
nema.tidy <- rownames_to_column(df.sqrt.nema.plus,var = "Microorganism")
nema.melt=reshape2::melt(nema.tidy,variable.name = "sample_code", value.name = "Count")
nema.join=left_join(nema.melt,map, by = "sample_code")
nema.group=nema.join %>% group_by(Site)

 mutate(Sample_sum = sum(Count),
        RelAbund = Count / Sample_sum) %>%
 summarise(meanRelAbund = mean(RelAbund)) %>%
 ungroup() 

map.div






# make a data frame contain prevalence, mean of relative abundance using phyloseq
prevdf = apply(X = otu_table(PHYL),
               MARGIN = ifelse(taxa_are_rows(PHYL), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
RelAbund = taxa_sums(PHYL.ra)
prevdf = data.frame(Prevalence = prevdf,
                    RelAbund = taxa_sums(PHYL.ra), MeanRelAbund = RelAbund/45,
                    tax_table(PHYL))
head(prevdf)
sort_prevdf <- prevdf[order(prevdf$MeanRelAbund, decreasing = TRUE),]
head(sort_prevdf)
dim(sort_prevdf)
filter_1 = prevdf[prevdf$Prevalence==45,]
dim(filter_1)
sort_filter_1 <- filter_1[order(filter_1$MeanRelAbund, decreasing = TRUE),]
plyr::ddply(sort_filter_1, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


ggplot(sort_filter_1, aes(MeanRelAbund*100, Prevalence / 45,color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("log(Mean of Relative Abundance)") + ylab("Occupancy") +
  facet_wrap(~Phylum) + theme(legend.position="none")

class(filter_1)

filter_2 <- sort_filter_1 %>% 
  group_by(Phylum) %>% 
  summarise(PhylumMeanRelAbund = sum(MeanRelAbund))

ggplot(filter_2, aes(x = reorder(Phylum, PhylumMeanRelAbund), y = PhylumMeanRelAbund*100))+ 
 geom_bar(position = "dodge",stat = "identity")+
 coord_flip()+
 theme_bw()+
 scale_y_continuous(position = "right")+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
 labs(x= "Phylum", y= "Mean of relative abundance (%)")

#### NETWORK ANALYSIS #####
# 1. Subsample the otu table (occupancy 100 %) -Bacteria
install.packages("OTUtable")
library(OTUtable)
Occ1_OTU <- filter_taxa(otu, abundance = 0, persistence = 100)
head(Occ1_OTU) 
dim(Occ1_OTU)
Occ1_OTU.Bac <- rownames_to_column(Occ1_OTU, var = "OTU")
# 2. Subsample the otu table (occupancy 100 %) - Fungi
Occ1_OTU.ITS <- filter_taxa(otuITS, abundance = 0, persistence = 100)
head(Occ1_OTU.ITS)
dim(Occ1_OTU.ITS)
Occ1_OTU.Fungi <- rownames_to_column(Occ1_OTU.ITS, var = "OTU")
# 3. Subsample the otu table (occupancy 100 %) - Other microorganisms & nematode
Occ1_nema <- filter_taxa(nema.plus, abundance = 0, persistence = 100)
head(Occ1_nema)
dim(Occ1_nema)
Occ1_nema <- rownames_to_column(as.data.frame(Occ1_nema), var = "OTU")
# Join three datasets
total = rbind(Occ1_OTU.Bac,Occ1_OTU.Fungi,Occ1_nema)  
# Write table
write.table(total, file = "TotalOTU.txt", quote = F, sep = "\t",
            row.names = FALSE)

# 1. Total OTU - Bac
BacOTU <- rownames_to_column(otu, var = "OTU")
head(BacOTU)
# 2. Total OTU - Fungi
FungOTU <- rownames_to_column(otuITS, var = "OTU")
head(FungOTU)
# 3. Total Nema and other
Nema <- rownames_to_column(as.data.frame(nema.plus), var = "OTU")
# Join three datasets
tot = rbind(BacOTU,FungOTU,Nema) 
# Write table
write.table(tot, file = "allofOTU.txt", quote = F, sep = "\t",
            row.names = FALSE)

BacTaxOTU <- rownames_to_column(otu, var = "OTU")
BacTaxOTU[BacTaxOTU$OTU == "FPLL01002138.3.1535",]
#191 D_0__Bacteria; D_1__Verrucomicrobia; D_2__Verrucomicrobiae; D_3__Pedosphaerales; D_4__Pedosphaeraceae; Ambiguous_taxa; Ambiguous_taxa
BacTaxOTU[BacTaxOTU$OTU == "Z95733.1.1547",]
#D_0__Bacteria; D_1__Acidobacteria; D_2__Subgroup 6; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa
BacTaxOTU[BacTaxOTU$OTU == "JX898177.1.1545",]
#D_0__Bacteria; D_1__Acidobacteria; D_2__Subgroup 6; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa; Ambiguous_taxa


BacTaxOTU[BacTaxOTU$OTU == "FPLS01015067.17.1524",]
#675 D_0__Bacteria; D_1__Bacteroidetes; D_2__Bacteroidia; D_3__Chitinophagales; D_4__Chitinophagaceae; D_5__uncultured; Ambiguous_taxa
BacTaxOTU[BacTaxOTU$OTU == "DQ787718.1.1532",]
#153 D_0__Bacteria; D_1__Verrucomicrobia; D_2__Verrucomicrobiae; D_3__Verrucomicrobiales; D_4__Rubritaleaceae; D_5__Luteolibacter; Ambiguous_taxa

BacTaxOTU[BacTaxOTU$OTU == "OTU_dn_1211",]

### Add taxonomy for bacteria and fungi in node attribute table for network analysis
# 1. bacteria
otu <- read.table('OTU_rarefied.txt', sep='\t', header=T, row.names = 1)
df.otu=rownames_to_column(otu, var = "OTU_ID")
df.otu.tax <- subset(df.otu, select = c(1, 47))
dim(df.otu.tax)
write.table(df.otu.tax, file = "df.otu.tax.csv", sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
df.otu.tax.csv = read.csv("df.otu.tax.csv", sep=',', header=T)
dim(df.otu.tax.csv)
colnames(df.otu.tax.csv)[1]<-"Name"
head(df.otu.tax.csv)


allofOTU.node = read.table("allofOTU.node_attribute_copy.txt", sep='\t', header=T)
head(allofOTU.node)
dim(allofOTU.node)
allofOTU.node.bac <- allofOTU.node[1:380,]
dim(allofOTU.node.bac)
head(allofOTU.node.bac)

allofOTU.node.bac.tax <- merge(df.otu.tax.csv,allofOTU.node.bac, by.x =c("Name"), by.y = c("Name"))
dim(allofOTU.node.bac.tax)
head(allofOTU.node.bac.tax)

# 2. Fungi
allofOTU.node.fg <- allofOTU.node[381:426,]
dim(allofOTU.node.fg)

dim(tax)
head(tax)
tax.fg <- as.data.frame(tax)
tax.fg=rownames_to_column(tax.fg, var = "Name")
colnames(tax.fg)[2]<-"Domain"

allofOTU.node.fg.tax <- merge(tax.fg,allofOTU.node.fg, by.x =c("Name"), by.y = c("Name"))
dim(allofOTU.node.fg.tax)
head(allofOTU.node.fg.tax)

# 3. Combine bac and fungi
allofOTU.node_attribute <- rbind(allofOTU.node.bac.tax, allofOTU.node.fg.tax)
head(allofOTU.node_attribute)
write.table(allofOTU.node_attribute, file = "allofOTU.node_attribute.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

# Make Zi-Pi Plot
Module.hub <- allofOTU.node_attribute[allofOTU.node_attribute$Zi>2.5,] #12
write.table(Module.hub, file = "Module.hub.txt", sep = '\t', col.names = TRUE, row.names = F, quote = FALSE)
Module.hub.tax <- read.table("Module.hub.edit.txt", sep='\t', header=T)
Module.col=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff')
non.hub <- allofOTU.node_attribute[allofOTU.node_attribute$Zi<2.5,]
Peripherals <- non.hub[non.hub$Pi<0.62,] #409
Connectors <- allofOTU.node_attribute[allofOTU.node_attribute$Pi>0.62,] #5
write.table(Connectors, file = "Connectors.txt", sep = '\t', col.names = TRUE, row.names = F, quote = FALSE)
Connectors.tax <- read.table("Connectors.edit.txt", sep='\t', header=T)

ZP.plot1 <- ggplot()+ geom_point(size=2.5, aes(y = Zi, x = Pi),
                           data = Peripherals)+
 xlab("Among-module connectivity (Pi)")+
 ylab("Within-module connectivity (Zi)")+
 geom_hline(yintercept = 2.5, colour="#990000", linetype="dashed") +
 geom_vline(xintercept = 0.62, colour="#990000", linetype="dashed")+
annotate("text", x = 0.08, y = 4.25, label = "Module hubs", fontface='bold')+
 annotate("text", x = 0.75, y = 4.25, label = "Network hubs",fontface='bold')+
 annotate("text", x = 0.08, y = -2, label = "Peripherals",fontface='bold')+
 annotate("text", x = 0.75, y = -2, label = "Connectors",fontface='bold')+
 xlim(0, 0.8)+
 theme_bw()+
 theme(axis.text=element_text(size=9), 
       axis.title=element_text(size=10,face="bold"),
legend.text=element_text(size = 8),
legend.title = element_text(size=10),
 panel.grid = element_blank())
ZP.plot2 <- ZP.plot1 + 
 geom_point(size=2.5, 
            mapping=aes(y = Zi, x = Pi, colour=Module_hub),
            data = Module.hub.tax)+
 scale_color_manual(name="Module Hubs", values = Module.col)
ZP.plot3 <- ZP.plot2 + geom_point(size=2.5, mapping=aes(y = Zi, x = Pi, shape=Connectors),
            data = Connectors.tax)
 # 2. Install Required packages 
install.packages("igraph") 
install.packages("qgraph") 
install.packages("MCL")
library(igraph)
library(qgraph)
library(MCL)
# Install SpiecEasi package 
install.packages("devtools")
library(devtools) 
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
# Proccess the absolute abundance of otu table into relative abundance
Occ1.RelAbund <- decostand(Occ1_OTU, method="total", MARGIN=2)
head(Occ1.RelAbund)

# Dissimilarity based network
distances <- vegdist(t(Occ1.RelAbund), method = "bray")
# Convert distance object to a matrix
diss.mat <- as.matrix(distances)
diss.cutoff <- 0.6
diss.adj <- ifelse(diss.mat <= diss.cutoff, 1, 0)
diss.net <- graph.adjacency(diss.adj,mode = "undirected", diag = FALSE)
plot(diss.net)
# Correlation based network
cor.matrix <- cor(Occ1.RelAbund, method = "pearson")
# Convert correlation matrix to binary adjacency matrix 
cor.cutoff <- 0.3 
cor.adj <- ifelse(abs(cor.matrix) >= cor.cutoff, 1, 0) 
# Construct microbiome network from adjacency matrix 
cor.net <- graph.adjacency(cor.adj,
mode = "undirected", diag = FALSE)
plot(cor.net)






head(otu)
head(sort(rowSums (otu, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(Mean_rel.abund)
rSums=rowSums(otu)
otu.relrows=otu/rSums
mean.oturelrows=mean(otu.relrows)

otu_tidy=rownames_to_column(otu, var = "OTU")
head(otu_tidy)
otu_tidy_melt=reshape2::melt(otu_tidy, variable.name = "sample_code", value.name = "Count")
Sample_sum=sum(otu_tidy_melt$Count)
Relab=otu_tidy_melt$Count/Sample_sum
otu_tidy_melt=left_join(otu_tidy_melt, map, by = "sample_code")
otu_tidy_melt$sample_code=as.factor(otu_tidy_melt$sample_code)
otu_tidy_melt%>%group_by(sample_code)
otu_tidy_melt%>%mutate(Sample_sum,Relab)

head(otu_table(PHYL.ra))==head(otu/colSums(otu))

##### MultiCola ######

u=unique(map[,"sample_code"])
MantelMultiCOLA.f=function(otu){
  #Step 1.  Read in full dataset:
  otu2=otu
  library(vegan)
  cutoff=c(1.00,0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.025,
0.01, 0.005, 0.001)
  m.out=NULL
  otu.pa=1*(otu2>0)
  r=rowSums(otu2)
  #Create a vector of indexes of the ranked OTUs
  r2=sort(r,decreasing=TRUE, index=TRUE)
  r2$ix
  r2$x
  for(j in 1:length(cutoff)){
    print(cutoff[j])
    no.keep=ceiling(cutoff[j]*nrow(otu2))
    otu.index.keep=r2$ix[1:no.keep]
    otu.keep=otu2[otu.index.keep,]
    print(head(otu.keep))
#Write out otu tables at each cutoff
    write.table(otu.keep, paste("_",cutoff[j],"_otu.txt", sep=""),sep="\t",
quote=FALSE)
    all.dist=vegdist(t(otu2), method="bray")
    subset.dist=vegdist(t(otu.keep),method="bray")
    m1=mantel(all.dist,subset.dist, method="pearson",permutations=999)
    m=c(paste(cutoff[j]),m1$statistic,m1$signif,dim(otu.keep)[1])
    m.out=rbind(m.out,m)
  }
  colnames(m.out)=c("Cutoff", "AllvSubsetPearsonR", "AllvSubset_pvalue",
"NoOTUsSubset")
  #write.table(m.out, "MantelMultiCOLA.txt", sep="\t", quote=FALSE, row.names=FALSE)
  return(m.out)
}

MultiCOLA.test=MantelMultiCOLA.f(otu)

# 1. ph in every site
ph_site <- lm(map_aov$pH ~ Site, data=map_aov, na.action=na.exclude)
#drop1(ph_site,~.,test="F")
ph_site_resids <- residuals(ph_site)
ph_site_preds <- predict(ph_site)
plot(ph_site_resids ~ ph_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
 abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
 shapiro.test(ph_site_resids) # p-value = 0.4, data errors are normally distributed 
 skew_xts <-  skewness(ph_site_resids)
#Perform Levene's Test for homogenity of variances
 leveneTest(ph_site, na.action=na.exclude) #0.2 variances among group are homogenous
# Significantly different

# 2. Acidobacteria abundance in every site
PHYL_16S
phylum.16S.ra <- transform_sample_counts(PHYL_16S, function(x) x/sum(x))
#phylum.16S <- tax_glom(phylum.16S.ra, taxrank = "Phylum", NArm = F)
#phylum.16S
#tax_table(phylum.16S)
phylum.16S.ra <- transform_sample_counts(PHYL_16S, function(x) x/sum(x))
phylum.16S <- tax_glom(phylum.16S.ra, taxrank = "Phylum", NArm = F)
phylum.16S
cumabun=taxa_sums(phylum.16S)
df.phylum.16S.taxasum=as.data.frame(cumabun)
df.phylum.16S.taxasum$relabun=df.phylum.16S.taxasum$cumabun/45
df.16S <- psmelt(phylum.16S)
summary(df.16S)
head(df.16S)
dim(df.16S)
df.16S$Phylum <- as.character(df.16S$Phylum)
df.acid=df.16S %>% filter(df.16S$Phylum==" Acidobacteria")
head(df.acid)
dim(df.acid)

acid_site <- lm(log10(df.acid$Abundance) ~ Site, data=df.acid, na.action=na.exclude)
summary(acid_site)
drop1(acid_site,~.,test="F")
acid_site_resids <- residuals(acid_site)
acid_site_preds <- predict(acid_site)
plot(acid_site_resids ~ acid_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
shapiro.test(acid_site_resids)
leveneTest(acid_site, na.action=na.exclude)

df.thau=df.16S %>% filter(df.16S$Phylum==" Thaumarchaeota")
head(df.thau)
dim(df.thau)






#anova(glm(N.Organisms~as.factor(Group), data=d, family=poisson))
anova(glm(df.acid$Abundance ~ Site, data=df.acid, family=poisson,na.action=na.exclude), test="LRT")

acid_ph <- lm(Abundance ~ P_ppm, data=df.acid)
summary(acid_ph)
cor.test(df.acid$P_ppm, df.acid$Abundance)
