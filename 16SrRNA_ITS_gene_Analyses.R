# Biogeography and diversity of multi-trophic root zone microbiomes in Michigan apple orchards: analysis of rootstock, scion, and growing region
# October 11th 2018
# Ari Fina Bintarti

# Install packages
install.packages(c('vegan', 'tidyverse'))
install.packages('reshape')
install.packages("dplyr")
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

library(BiocManager)
library(vegan)
library(tidyverse)
library(plyr)
library(dplyr)
library(tidyr)
#library(cowplot)
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
detach(package:plyr)

# Set the working directory
setwd('/Users/arifinabintarti/Documents/GitHub/PAPER_Bintarti_2020_Phytobiomes/')
wd <- print(getwd())

# Read the bacterial Operational Taxonomic Units (OTUs)
otu <- read.table('OTU_rarefied_16S.txt', sep='\t', header=T, row.names = 1)
taxonomy <- otu[,'taxonomy']
taxonomy
otu <- otu[,-46]
dim(otu) # there are 22,510 bacterial OTUs and 45 samples
# Read the fungal OTUs
otuITS <- read.table(file = "OTU_rarefied_ITS.txt", sep='\t', header=T, row.names = 1) 
dim(otuITS) # there are 3,553 fungal OTUs and 45 samples

# Read the metadata (map)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
str(map) # We use number as the name of site and it is integer 
map$Site <- as.factor(map$Site) # I want to change site as factor
colnames(map)[which(names(map) == "rootstock")] <- "Rootstock"
colnames(map)[which(names(map) == "cultivar")] <- "Scion"

# Calculate the alpha diversity (Richness and Pielou's evenness, we also calculates Shannon index) 
# 1. Bacterial and archaeal alpha diversity for each sample 
otu_rare_PA <- 1*(otu>0)
s <- specnumber(otu, MARGIN = 2) # richness
richness <- as.data.frame(s) 
h <- diversity(t(otu), index = 'shannon') # Shannon index
shannon <- as.data.frame(h)
pielou <- h/log(s) # Pielou's evenness
evenness <- as.data.frame(pielou)
map_df <- data.frame(map) # make data frame of the map data
map.div <- map_df
map.div$Richness <- s
map.div$Shannon <- h
map.div$Pielou <- pielou
# 2. Fungal alpha diversity for each sample
sITS <- specnumber(otuITS, MARGIN = 2) # richness
hITS <- diversity(t(otuITS), index = 'shannon') # Shannon index
pielouITS <- hITS/log(sITS) # Pielou's evenness
map.div$fg.Richness <- sITS
map.div$fg.Shannon <- hITS
map.div$fg.Pielou <- pielouITS

# Melt the map and alpha diversity variables
# map.alpha <- melt(map.div, id.vars=c('Site', 'Scion', 'Rootstock'), 
#                 measure.vars=c('Richness','Shannon', 'Pielou','fg.Richness','fg.Shannon', 'fg.Pielou'))

# 3. nematodes, mycorrhizal fungi, and oligochaetes alpha diversity
# make a data frame of nematode group only (11 group)
nema <- map_df[,c(1,24:34)]
dim(nema)
nema
row.names(nema) <- nema$sample_code
nema[1] <- NULL
nema.t <- t(nema)
nema.t
sort(rowSums(nema.t, na.rm = FALSE, dims = 1), decreasing = FALSE)
# make data frame of nematodes, oligochaetes, and mycorrhizal fungi
nema.plus <- map_df[,c(1,24:36)]
row.names(nema.plus) <- nema.plus$sample_code
nema.plus[1] <- NULL
nema.plus <- t(nema.plus)
nema.plus
sqrt.nema.plus <- sqrt(nema.plus) # square root the nematode count data
otu_dist_nema.plus <- vegdist(sqrt.nema.plus, method='bray') # calculate dissimilarity indices using 'bray-curtis' method
otu_pcoa_nema.plus <- cmdscale(otu_dist_nema.plus, eig=T) # CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
# calculate the count of nematodes for each sample
nema_PA <- 1*(nema.t>0)
nema_total=as.data.frame(rowSums(nema.t)) # calculate total count data for every nematode group
colnames(nema_total)="Count"
nema_total <- rownames_to_column(nema_total, var = "Group")
nema.total.Count=colSums(nema.t)
sqrt.nema.t=t(sqrt(nema))
nema.total.SqrtCount=colSums(sqrt.nema.t)
map.div$nema.total.Count=nema.total.Count
map.div$nema.total.SqrtCount=nema.total.SqrtCount
# transform the mycorrhizal fungi count data using log10 transformation
log.MF=log10(map.div$MycorrhizalFungi) 
map.div$log.MF=log.MF
# transform the oligochaetes count data using square root transformation
sqrt.OL=sqrt(map.div$Oligochaetes)
map.div$sqrt.OL=sqrt.OL
## Calculate the alpha diversity (richness, Shannon index, and Pielou's evenness) of nematodes
s.nema <- specnumber(nema.t, MARGIN = 2) # richness is computed as the number of taxa present per sample
h.nema <- diversity(nema, index = 'shannon')
pielou.nema <- h.nema/log(s.nema)
map.div$nema.Richness <- s.nema
map.div$nema.Shannon <- h.nema
map.div$nema.Pielou <- pielou.nema
names(map.div)

# ANOVA and Kruskal-Wallis test to compare alpha diversity among geographical sites, rootstocks, and scion

map_aov <- map.div
class(map_aov)
map_aov$Site<-as.factor(map_aov$Site) # inform R that Site is factor
str(map_aov, give.attr=F)

## Statistical test of bacterial and archaeal richness ##

### 1. Compare bacterial and archaeal richness among sites using one-way ANOVA
Aov_richness_site <- lm(map_aov$Richness ~ Site, data=map_aov, na.action=na.exclude)
Aov_richness_site
drop1(Aov_richness_site,~.,test="F") # type III SS and F Tests
# testing assumptions
# Generate residual and predicted values
RC_site_resids <- residuals(Aov_richness_site)
RC_site_preds <- predict(Aov_richness_site)
# Look at a plot of residual vs. predicted values
plot(RC_site_resids ~ RC_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
# plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_site_resids) # ALERT!! p-value = 0.0006956, data errors are not normally distributed 
# Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ Site, data=map_aov, na.action=na.exclude) # GOOD variances among group are homogenous
boxplot(map_aov$Richness ~ Site, data = map_aov) # there are no outliers
plot(density(RC_site_resids)) # density is not bad
qqnorm(RC_site_resids)
qqline(RC_site_resids) # I think the data normality is fine
hist(RC_site_resids)
skew_xts <-  skewness(RC_site_resids) # -0.224 (within range of -2 to 2)
kurtosis(RC_site_resids,method = 'sample')# 8.7 (should be in the range of -7 to 7)
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.

### 1. Compare bacterial and archaeal richness among sites using Kruskal-Wallis Test
kruskal.test(Richness ~ Site, data = map_aov) # Kruskal-Wallis chi-squared = 32.155, df = 19, p-value = 0.03002
ggboxplot(map_aov, x = "Site", y = "Richness")+
  stat_compare_means()
### RESULT: There are significant differences of bacterial and archaeal richness among sites

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
  labs(title = "A. Bacteria/archaea")+
 # geom_text(data=new.richness.summarized,aes(x=Site,y=32+max.Richness,label=new.richness.summarized$groups),vjust=0)+
 theme(axis.text.x=element_blank(),
       axis.title.x=element_blank(),
       axis.ticks.x=element_blank(),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

### 2. Compare bacterial and archaeal richness among rootstocks using one-way ANOVA
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
### RESULT: The data does not meet the homocedasticity assumption, thus I use Kruskal-Wallis test instead.

### 2. Compare bacterial and archaeal richness among rootstocks using Kruskal-Wallis Test
kruskal.test(Richness ~ Rootstock, data = map_aov) # Kruskal-Wallis chi-squared = 16.367, df = 7, p-value = 0.02197
### RESULT: There are significant differences of bacterial and archaeal richness among rootstocks

# Do Post Hoc Dunn's Test
DT_RC_root <- dunnTest(Richness~Rootstock, map_aov, method = "bh", kw=TRUE)
print(DT_RC_root,dunn.test.results=TRUE)
DT_RC_root$res
DT_RC_root.df <- as.data.frame(DT_RC_root$res)
#write.csv(DT_RC_root.df, file = "DT_RC_root.df.csv")
DT_RC_root_letter = cldList(P.adj ~ Comparison,
        data = DT_RC_root$res,
        threshold = 0.05)
# Do Plot
(bac_rich_root <- ggplot(map.div, aes(x=Rootstock, y=Richness))+
 geom_boxplot() +
 geom_point() +
  theme_bw()+
  labs(title = "A. Bacteria/archaea")+
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

### 3. Compare bacterial and archaeal richness among scions using one-way ANOVA
Aov_rich_cul <- lm(Richness ~ Scion, data=map_aov, na.action=na.exclude)
drop1(Aov_rich_cul,~.,test="F") # type III SS and F Tests
summary(Aov_rich_cul)
# testing assumptions
# Generate residual and predicted values
RC_cul_resids <- residuals(Aov_rich_cul)
RC_cul_preds <- predict(Aov_rich_cul)
plot(RC_cul_resids ~ RC_cul_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(RC_cul_resids) # GOOD p-value = 0.1793, data errors are normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(Richness ~ Scion, data=map_aov, na.action=na.exclude) # GOOD!! p-val=0.7019, variances among group are homogenous
### RESULT: There are no significant differences of bacterial and archaeal richness among scions

## Statistical test of bacterial and archaeal Shannon index ##

### 1. Compare bacterial and archaeal Shannon index among sites using one-way ANOVA
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
### RESULT: The data does not meet the normality assumption, thus I use Kruskal-Wallis test instead.

### 1. Compare bacterial and archaeal Shannon index among sites using Kruskal-Wallis Test
kruskal.test(Shannon ~ Site, data = map_aov) # Kruskal-Wallis chi-squared = 35.368, df = 19, p-value = 0.01261
### RESULT: There are significant differences of bacterial and archaeal Shannon index among sites

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
 
### 2. Compare bacterial and archaeal Shannon index among rootstocks using one-way ANOVA
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
### RESULT: There are significant differences of bacterial and archaeal Shannon index among rootstocks

# Do Tukey's HSD Post Hoc Test
hsd_Shannon_rootstock<- HSD.test(Aov_shannon_rootstock, "Rootstock", group = TRUE, console = TRUE)
hsd_Shannon_rootstock<- HSD.test(Aov_shannon_rootstock, "Rootstock", group = F, console = TRUE)
# Do Plot
# add significance letters from HSD.test into box plot
detach(package:plyr)
root.shannon.summarized <- map.div %>% 
 group_by(Rootstock) %>%
 summarise(max.Shannon = max(Shannon))
hsd_Shannon_rootstock <- HSD.test(Aov_shannon_rootstock,
                                  "Rootstock", 
                                  group = TRUE, 
                                  console = TRUE)
hsd_Shannon <- hsd_Shannon_rootstock$groups
class(hsd_Shannon)
hsd_Shannon$Rootstock <- rownames(hsd_Shannon)
new.root.shannon.summarized <- left_join(hsd_Shannon,
                                      root.shannon.summarized,
                                      by='Rootstock')  

(bac_sha_root <- ggplot(map.div, aes(x=Rootstock, y=Shannon))+
 geom_boxplot() +
 geom_point() +
  theme_bw()+
  labs(title = "D")+
 geom_text(data=new.root.shannon.summarized,aes(x=Rootstock,y=0.01+max.Shannon,label=new.root.shannon.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=14,angle=49,hjust =1),
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=20,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       #axis.title.x = element_text(vjust=15),
       axis.title=element_text(size=18,face="bold", vjust = 10),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

## Statistical test of bacterial and archaeal Pielou's evenness ##

### 1. Compare bacterial and archaeal Pielou's evenness among sites using one-way ANOVA
Aov_pie_site <- lm(Pielou ~ Site, data=map_aov, na.action=na.exclude)
drop1(Aov_pie_site,~.,test="F") # type III SS and F Tests
pie_site_resids <- residuals(Aov_pie_site)
pie_site_preds <- predict(Aov_pie_site)
plot(pie_site_resids ~ pie_site_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(Pielou ~ Site, data=map_aov, na.action=na.exclude) # GOOD the data is homogen
#Test for normality
shapiro.test(pie_site_resids) # GOOD!! p-val>0.05, the data are normally distributed
plot(density(pie_site_resids))
qqnorm(pie_site_resids)
qqline(pie_site_resids)
hist(pie_site_resids)
skew_xts <- skewness(pie_site_resids)
### RESULT: There are significant differences of bacterial and archaeal Pielou's evenness among sites

# Do Tukey's HSD Post Hoc Test
hsd_pie_site<- HSD.test(Aov_pie_site, "Site", group = TRUE, console = TRUE)
hsd_pie_site<- HSD.test(Aov_pie_site, "Site", group = F, console = TRUE)
# Do Plot
# add significance letters from HSD.test into box plot
site.pie.summarized <- map.div %>% group_by(Site) %>% summarize(max.pie=max(Pielou))
hsd_pie_site<- HSD.test(Aov_pie_site, "Site", group = TRUE, console = TRUE)
hsd_pie = hsd_pie_site$groups
class(hsd_pie)
hsd_pie$Site <- rownames(hsd_pie)
new.site.pie.summarized=left_join(hsd_pie,site.pie.summarized, by='Site')  

(bac_pie_site <- ggplot(map.div, aes(x=Site, y=Pielou))+
 geom_boxplot() +
 geom_point() +
  theme_bw()+
  labs(y="Pielou's evenness")+
  scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 geom_text(data=new.site.pie.summarized,aes(x=Site,y=0.0008+max.pie,label=new.site.pie.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

### 2. Compare bacterial and archaeal Pielou's evenness among rootstocks using one-way ANOVA
Aov_pie_root <- lm(Pielou ~ Rootstock, data=map_aov, na.action=na.exclude)
drop1(Aov_pie_root,~.,test="F") # type III Sum of Squares and F Tests
summary(Aov_pie_root)
# testing assumptions
# Generate residual and predicted values
pie_root_resids <- residuals(Aov_pie_root)
pie_root_preds <- predict(Aov_pie_root)
plot(pie_root_resids ~ pie_root_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(pie_root_resids) # GOOD p-value = 0.062, data errors are normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ Rootstock, data=map_aov, na.action=na.exclude) # p-val=0.5436, variances among group are homogenous
#Plotting
boxplot(Pielou ~ Rootstock, data = map_aov) # there are no outliers
plot(density(pie_root_resids))
qqnorm(pie_root_resids)
qqline(pie_root_resids) 
hist(pie_root_resids)
### RESULT: There are significant differences of bacterial and archaeal Pielou's evenness among rootstocks

# Do Tukey's HSD Post Hoc Test
#install.packages("lsmeans")
#install.packages("emmeans")
#library(lsmeans)
#library(emmeans)
#pairs(lsmeans(Aov_pie_root, "Rootstock"), adjust = "mvt")
hsd_pie_root<- HSD.test(Aov_pie_root, "Rootstock", group = TRUE, console = TRUE, alpha = 0.05)
hsd_pie_root<- HSD.test(Aov_pie_root, "Rootstock", group = F, console = TRUE, alpha = 0.05)
# Do Plot
# add significance letters from HSD.test into box plot
root.pie.summarized <- map.div %>% group_by(Rootstock) %>% summarize(max.pie=max(Pielou))
hsd_pie_root<- HSD.test(Aov_pie_root, "Rootstock", group = TRUE, console = TRUE)
hsd_pieRoot = hsd_pie_root$groups
class(hsd_pieRoot)
hsd_pieRoot$Rootstock <- rownames(hsd_pieRoot)
new.root.pie.summarized=left_join(hsd_pieRoot,root.pie.summarized, by='Rootstock') 
# Do Plot
(bac_pie_root <- ggplot(map.div, aes(x=Rootstock, y=Pielou))+
 geom_boxplot() +
 geom_point() +
  theme_bw()+
  labs(y="Pielou's evenness")+
geom_text(data=new.root.pie.summarized,aes(x=Rootstock,y=0.0008+max.pie,label=new.root.pie.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=14,angle=49,hjust =0.9),
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=20,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

### 3. Compare bacterial and archaeal Pielou's evenness among scions using one-way ANOVA
Aov_pie_cul <- lm(Pielou ~ Scion, data=map_aov, na.action=na.exclude)
drop1(Aov_pie_cul,~.,test="F") # type III SS and F Tests
summary(Aov_pie_cul)
# testing assumptions
# Generate residual and predicted values
pie_cul_resids <- residuals(Aov_pie_cul)
pie_cul_preds <- predict(Aov_pie_cul)
plot(pie_cul_resids ~ pie_cul_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
#Perform a Shapiro-Wilk test for normality of residuals
shapiro.test(pie_cul_resids) # GOOD p-value = 0.3038, data errors are normally distributed 
#Perform Levene's Test for homogenity of variances
leveneTest(Pielou ~ Scion, data=map_aov, na.action=na.exclude) # p-val=0.6497, variances among group are homogenous
### RESULT: There are no significant differences of bacterial and archaeal Pielou's evenness among scions

# 2. ANOVA TEST FOR FUNGAL ALPHA DIVERSITY

### FUNGAL RICHNESS ####

# 1. Compare richness among site - ANOVA
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
  labs(title = "B. Fungi")+
 geom_text(data=new.ITSsite.richness.summarized,aes(x=Site,y=10+max.Richness,label=new.ITSsite.richness.summarized$groups),vjust=0)+
 theme(axis.text.x=element_blank(),
       axis.title.x=element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.x=element_blank(),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 2. Compare richness among rootstocks - ANOVA
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

(fg_rich_root <- ggplot(map.div, aes(x=Rootstock, y=fg.Richness))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
  labs(title = "B. Fungi")+
 geom_text(data=new.ITSroot_richness.summarized,aes(x=Rootstock,y=10+max.Richness,label=new.ITSroot_richness.summarized$groups),vjust=0)+
 theme(axis.text=element_blank(),
       axis.title.y=element_blank(),
       axis.ticks.x = element_blank(), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_blank(), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)), 
       axis.title.x = element_blank(),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 3. Compare richness among scions - KRUSKAL-WALLIS
AovITS_richness_cul <- lm(fg.Richness ~ Scion , data=map_aov, na.action=na.exclude)
drop1(AovITS_richness_cul,~.,test="F") # type III SS and F Tests
summary(AovITS_richness_cul)
ITSRichness_cul_resids <- residuals(AovITS_richness_cul)
ITSRichness_cul_resids
ITSRichness_cul_preds <- predict(AovITS_richness_cul)
plot(ITSRichness_cul_resids ~ ITSRichness_cul_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Richness ~ Scion, data=map_aov, na.action=na.exclude) #GOOD! p-val = 0.909, the data is homogen
#Test for normality
shapiro.test(ITSRichness_cul_resids) ## ALERT! p-val=0.006714, the data are NOT normally distributed
qqnorm(ITSRichness_cul_resids)
qqline(ITSRichness_cul_resids)
# Do Kruskal-Wallis Test
kruskal.test(fg.Richness ~ Scion, data = map_aov) # Kruskal-Wallis chi-squared = 35.368, df = 19, p-value = 0.01261

### FUNGAL SHANNON INDEX ####

# 1. Compare Shannon index among sites
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

# 2. Compare Shannon index among rootstocks
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
       #axis.title.x = element_text(vjust=10),
       axis.title=element_text(size=18,face="bold", vjust = 10),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

### FUNGAL PIELOU'S EVENNESS INDEX ###

# 1. Compare Pielou among sites - Not Significant
AovITS_pie_site <- lm(fg.Pielou ~ Site , data=map.div, na.action=na.exclude)
drop1(AovITS_pie_site,~.,test="F") # type III SS and F Tests
summary(AovITS_pie_site)
ITSpie_site_resids <- residuals(AovITS_pie_site)
ITSpie_site_preds <- predict(AovITS_pie_site)
plot(ITSpie_site_resids ~ ITSpie_site_preds , xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Pielou ~ Site, data=map_aov, na.action=na.exclude) #p-val > 0.05, the data is homogen
#Test for normality
shapiro.test(ITSpie_site_resids) ## p-val>0.05, the data are normally distributed
qqnorm(ITSpie_site_resids)
qqline(ITSpie_site_resids)
# Do plot
(fg_pie_site <- ggplot(map.div, aes(x=Site, y=fg.Pielou))+
 geom_boxplot() +
 geom_point() +
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 theme_bw()+
 #geom_text(data=new.ITSsite.shannon.summarized,aes(x=Site,y=0.03+max.Shannon,label=new.ITSsite.shannon.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       axis.title.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 2. Compare Pielou among rootstocks
AovITS_pie_root <- lm(fg.Pielou ~ Rootstock , data=map_aov, na.action=na.exclude)
drop1(AovITS_pie_root,~.,test="F") # type III SS and F Tests
summary(AovITS_pie_root)
ITSpie_root_resids <- residuals(AovITS_pie_root)
ITSpie_root_resids
ITSpie_root_preds <- predict(AovITS_pie_root)
plot(ITSpie_root_resids ~ ITSpie_root_preds , xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Pielou ~ Rootstock , data=map_aov, na.action=na.exclude) #p-val > 0.05, the data is homogen
#Test for normality
shapiro.test(ITSpie_root_resids) ## p-val>0.05, the data are normally distributed
qqnorm(ITSpie_root_resids)
qqline(ITSpie_root_resids)
# Do Plot 
(fg_pie_root <- ggplot(map.div, aes(x=Rootstock, y=fg.Pielou))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
 #geom_text(data=new.ITSroot.shannon.summarized,aes(x=Rootstock,y=0.03+max.Shannon,label=new.ITSroot.shannon.summarized$groups),vjust=0) +
 theme(axis.text.x=element_text(size=14,angle=49,hjust =0.9),
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=20,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_text(size=18,face="bold"),
       axis.title.y = element_blank(),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 3. Compare Pielou among scions
AovITS_pie_cul <- lm(fg.Pielou ~ Scion , data=map_aov, na.action=na.exclude)
drop1(AovITS_pie_cul,~.,test="F") # type III SS and F Tests
summary(AovITS_pie_cul)
ITSpie_cul_resids <- residuals(AovITS_pie_cul)
ITSpie_cul_resids
ITSpie_cul_preds <- predict(AovITS_pie_cul)
plot(ITSpie_cul_resids ~ ITSpie_cul_preds , xlab = "Predicted Values", ylab = "Residuals", asp=.5)
#plot a line on X-axis for comparison. a=intercept,b=slope,h=yvalue,v=xvalue.
abline(a=0,h=0,b=0)
# test for homogeneity of variance
#levene test was used here, the null hypothesis is all variances are equal, 
# p-val < 0.05 means all variances are not equal, thus parametric tests such as ANOVA are not suited. 
leveneTest(fg.Pielou ~ Scion , data=map_aov, na.action=na.exclude) #p-val > 0.05, the data is homogen

############### PLOT BACTERIAL AND FUNGAL ALPHA DIVERSITY ###############################################
# 1. Bacterial/archaeal and fungal Richness and Pielou's evenness among sites
plot <- grid.arrange(bac_rich_site,                                    
             fg_rich_site, bac_pie_site, fg_pie_site,                              
             ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)))
ggsave("Supplementary Figure S3.eps",
       plot, device = "eps",
       width = 12, height = 7, 
       units= "in", dpi = 600)
# 2. Bacterial/archaeal and fungal Richness and Pielou's evenness among rootstocks
plot <- grid.arrange(bac_rich_root,                                    
             fg_rich_root, bac_pie_root, fg_pie_root,                              
             ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)))
ggsave("Figure_1.eps",
       plot, device = "eps",
       width = 12, height = 7, 
       units= "in", dpi = 600)
###########################################################################################

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

# 3. Pearson correlation of soil properties with bacterial Pielou's evenness index
bac.pie_pH <- cor(map.div$pH, map.div$Pielou, method ="pearson")
cor.test(map.div$pH, map.div$Pielou)

bac.pie_P <- cor(map.div$P_ppm, map.div$Pielou, method ="pearson")
cor.test(map.div$P_ppm, map.div$Pielou)

bac.pie_K <- cor.test(map.div$K_ppm, map.div$Pielou, method ="pearson")
cor.test(map.div$K_ppm, map.div$Pielou)

bac.pie_Ca <- cor.test(map.div$Ca_ppm, map.div$Pielou, method ="pearson")
cor.test(map.div$Ca_ppm, map.div$Pielou)

bac.pie_Mg <- cor.test(map.div$Mg_ppm, map.div$Pielou, method ="pearson")
cor.test(map.div$Mg_ppm, map.div$Pielou)

bac.pie_NO3N <- cor.test(map.div$NO3N_ppm, map.div$Pielou, method ="pearson")
cor.test(map.div$NO3N_ppm, map.div$Pielou)

bac.pie_NH4N <- cor.test(map.div$NH4N_ppm, map.div$Pielou, method ="pearson")
cor.test(map.div$NH4N_ppm, map.div$Pielou)

bac.pie_OM <- cor.test(map.div$OM_percent, map.div$Pielou, method ="pearson")
cor.test(map.div$OM_percent, map.div$Pielou)

bac.pie_SD <- cor.test(map.div$sand_percent, map.div$Pielou, method ="pearson")
cor.test(map.div$sand_percent, map.div$Pielou)

bac.pie_SL <- cor.test(map.div$silt_percent, map.div$Pielou, method ="pearson")
cor.test(map.div$silt_percent, map.div$Pielou)

bac.pie_CL <- cor.test(map.div$clay_percent, map.div$Pielou, method ="pearson")
cor.test(map.div$clay_percent, map.div$Pielou)

# effect of soil type to bacterial pielou's evennes index # not significant
bac.pie_soil.type <- lm(Pielou ~ soil_type, data=map.div)
summary(bac.pie_soil.type)
drop1(bac.pie_soil.type,~.,test="F") # type III SS and F Tests
TYPE_pie_resids <- residuals(bac.pie_soil.type)
TYPE_pie_preds <- predict(bac.pie_soil.type)
plot(TYPE_pie_resids ~ TYPE_pie_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
leveneTest(Pielou ~ soil_type, data=map.div, na.action=na.exclude) # p-val 0.01276, data is not homogen
shapiro.test(TYPE_pie_resids) # data is normal
# fail test assumption for anova
# use WELCH-ANOVA instead
bac.pie_soil.type <- oneway.test(Pielou ~ soil_type, data=map.div, var.equal = F)
bac.pie_soil.type
ggplot(map.div, aes(x=soil_type, y=Pielou)) + geom_boxplot()
# F = 2.637, num df = 2.0000, denom df = 2.6539, p-value = 0.2341

#### FUNGI ######

# 1. Pearson correlation of soil properties with fungal Richness # all are not significant
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

# 2. Pearson correlation of soil properties with fungal Shannon index
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

# 3. Pearson correlation of soil properties with fungal Pielou's evenness
fg.pie_pH <- cor(map.div$pH, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$pH, map.div$fg.Pielou)

fg.pie_P <- cor(map.div$P_ppm, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$P_ppm, map.div$fg.Pielou)

fg.pie_K <- cor.test(map.div$K_ppm, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$K_ppm, map.div$fg.Pielou) # significant

fg.pie_Ca <- cor.test(map.div$Ca_ppm, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$Ca_ppm, map.div$fg.Pielou)

fg.pie_Mg <- cor.test(map.div$Mg_ppm, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$Mg_ppm, map.div$fg.Pielou) # significant

fg.pie_NO3N <- cor.test(map.div$NO3N_ppm, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$NO3N_ppm, map.div$fg.Pielou)

fg.pie_NH4N <- cor.test(map.div$NH4N_ppm, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$NH4N_ppm, map.div$fg.Pielou)

fg.pie_OM <- cor.test(map.div$OM_percent, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$OM_percent, map.div$fg.Pielou)

fg.pie_SD <- cor.test(map.div$sand_percent, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$sand_percent, map.div$fg.Pielou) # significant

fg.pie_SL <- cor.test(map.div$silt_percent, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$silt_percent, map.div$fg.Pielou) # significant

fg.pie_CL <- cor.test(map.div$clay_percent, map.div$fg.Pielou, method ="pearson")
cor.test(map.div$clay_percent, map.div$fg.Pielou) 

# effect of soil type to fungal Pielou's evenness index # significant
fg.pie_soil.type <- lm(fg.Pielou ~ soil_type, data=map.div)
summary(fg.pie_soil.type)
drop1(fg.pie_soil.type,~.,test="F") # type III SS and F Tests
TYPE_fg.pie_resids <- residuals(fg.pie_soil.type)
TYPE_fg.pie_preds <- predict(fg.pie_soil.type)
plot(TYPE_fg.pie_resids ~ TYPE_fg.pie_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
leveneTest(fg.Pielou ~ soil_type, data=map.div, na.action=na.exclude) # p-val 0.6, data is homogen
shapiro.test(TYPE_fg.pie_resids) # data is normal
# Tukey's HSD post-hoc test
HSD_fg.pie_soil.type <- HSD.test(fg.pie_soil.type, "soil_type", alpha = 0.05,group = TRUE,main = NULL,console=TRUE)
HSD_fg.pie_soil.type <- HSD.test(fg.pie_soil.type, "soil_type", alpha = 0.05,group = F,main = NULL,console=TRUE)
boxplot(fg.Pielou ~ soil_type, data=map.div)



#########################################################################################
## LINEAR REGRESSION TEST OF BACTERIAL AND FUNGAL ALPHA DIVERSITY WITH SOIL PROPERTIES ##
#########################################################################################

### BACTERIAL RICHNESS ###

# 1. BACTERIAL RICHNESS - SAND
Sand_percent_rich <- lm(Richness ~ sand_percent, data=map.div)
summary(Sand_percent_rich)
(a <- ggscatter(map.div, x = "sand_percent", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Sand (%)", ylab = "Richness")+
          annotate("text", x=42, y=5300, label = "y == 14.88(x)+3863.22", parse=T)+
          annotate("text", x=38, y=5200, label = "R^2 == 0.325", parse=T)+
          annotate("text", x=40, y=5100, label = "p-val == 4.38e-05", parse=T)+
     geom_smooth(method='lm')+labs(title = "A")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 2. BACTERIAL RICHNESS - SILT
Silt_percent_rich <- lm(Richness ~ silt_percent, data=map.div)
summary(Silt_percent_rich)
(b <- ggscatter(map.div, x = "silt_percent", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Silt (%)", ylab = "Richness")+
          annotate("text", x=39, y=5300, label = "y == -16.48(x)+5158.09", parse=T)+
          annotate("text", x=35, y=5200, label = "R^2 == 0.208", parse=T)+
          annotate("text", x=36, y=5100, label = "p-val == 0.001", parse=T)+
     geom_smooth(method='lm')+labs(title = "B")+
      theme(axis.text.x=element_text(size=15),
            axis.text.y = element_blank(),
       strip.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))
    
# 3. BACTERIAL RICHNESS - CLAY
Clay_percent_rich <- lm(Richness ~ clay_percent, data=map.div)
summary(Clay_percent_rich)
(c <- ggscatter(map.div, x = "clay_percent", y="Richness",
          add = "reg.line",conf.int = TRUE, xlab = "Clay (%)", ylab = "Richness")+
          annotate("text", x=20, y=5300, label = "y == -49.89(x)+5456.83", parse=T)+
          annotate("text", x=20, y=5200, label = "R^2 == 0.457", parse=T)+
          annotate("text", x=20, y=5100, label = "p-val == 3.376e-07", parse=T)+
     geom_smooth(method='lm')+labs(title = "C")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_blank(),
       strip.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       plot.title = element_text(size = rel(2)),
       axis.title.x=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))
    
#. 4. BACTERIAL RICHNESS - K
K_rich <- lm(Richness ~ K_ppm, data=map.div)
ggplot(K_rich, aes(x = K_ppm, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(K_rich)
(d <- ggscatter(map.div, x = "K_ppm", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "K (ppm)", ylab = "Richness")+
          annotate("text", x=200, y=5300, label = "y == -2.42(x)+5016.79", parse=T)+
          annotate("text", x=200, y=5200, label = "R^2 == 0.157", parse=T)+
          annotate("text", x=200, y=5100, label = "p-val == 0.006", parse=T)+
     geom_smooth(method='lm')+labs(title = "D")+
      theme(axis.text.x=element_text(size=15),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       plot.title = element_text(size = rel(2)),
       axis.title.x=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

#. 5. BACTERIAL RICHNESS - Ca 
#par(mar=c(6, 4, 4, 2) + 0.1)
Ca_rich <- lm(Richness ~ Ca_ppm, data=map.div)
ggplot(Ca_rich, aes(x = Ca_ppm, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(Ca_rich)
(e <- ggscatter(map.div, x = "Ca_ppm", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Ca (ppm)", ylab = "Richness")+
          annotate("text", x=2500, y=5300, label = "y == -0.18(x)+4891.12", parse=T)+
          annotate("text", x=2500, y=5200, label = "R^2 == 0.148", parse=T)+
          annotate("text", x=2500, y=5100, label = "p-val == 0.008", parse=T)+
     geom_smooth(method='lm')+labs(title = "E")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

#. 6. BACTERIAL RICHNESS - P
P_rich <- lm(Richness ~ P_ppm, data=map.div)
ggplot(P_rich, aes(x = P_ppm, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(P_rich)
(f <- ggscatter(map.div, x = "P_ppm", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "P (ppm)", ylab = "Richness")+
          annotate("text", x=45, y=5300, label = "y == 2.77(x)+4386.16", parse=T)+
          annotate("text", x=35, y=5200, label = "R^2 ==  0.218", parse=T)+
          annotate("text", x=35, y=5100, label = "p-val == 0.001", parse=T)+
     geom_smooth(method='lm')+labs(title = "F")+
      theme(axis.text.x=element_text(size=15), 
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       plot.title = element_text(size = rel(2)),
       axis.title.x=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

#. 7. BACTERIAL RICHNESS - OM
OM_rich <- lm(Richness ~ OM_percent, data=map.div)
ggplot(OM_rich, aes(x = OM_percent, y = Richness)) + geom_point() + stat_smooth(method = 'lm')
summary(OM_rich)
(g <- ggscatter(map.div, x = "OM_percent", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "OM (%)", ylab = "Richness")+
          annotate("text", x=2.3, y=4200, label = "y == -167.6(x)+5132.34", parse=T)+
          annotate("text", x=2, y=4100, label = "R^2 ==  0.178", parse=T)+
          annotate("text", x=2, y=4000, label = "p-val == 0.003", parse=T)+
     geom_smooth(method='lm')+labs(title = "G")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"),
       plot.title = element_text(size = rel(2)),
       axis.title.x=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

#. 8. BACTERIAL RICHNESS - nematode Tylenchs
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
(h <- ggscatter(map.div, x = "Tylenchs", y = "Richness", 
          add = "reg.line",conf.int = TRUE, xlab = "Tylenchus (per 100 g soil)", ylab = "Richness")+
          annotate("text", x=180, y=5200, label = "y == -1.55(x)+4753.92", parse=T)+
          annotate("text", x=180, y=5100, label = "R^2 ==  0.108", parse=T)+
          annotate("text", x=180, y=5000, label = "p-val == 0.026", parse=T)+
     geom_smooth(method='lm')+labs(title = "H")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=15,colour = "black", face = "bold"), 
       strip.text.y = element_blank(),
       plot.title = element_text(size = rel(2)),
       axis.title.x=element_text(size=15,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

### BACTERIAL PIELOU'S INDEX ###
# 1. BACTERIAL PIELOU'S INDEX - CLAY
Clay_percent_pie <- lm(Pielou ~ clay_percent, data=map.div)
summary(Clay_percent_pie)
# checking assumption
#par(mfrow=c(2,2))
gvlma(Clay_percent_pie)
plot(Clay_percent_pie)
# plot
ggplot(Clay_percent_pie, aes(x = clay_percent, y = Pielou)) + geom_point() + stat_smooth(method = 'lm')
(i <- ggscatter(map.div, x = "clay_percent", y = "Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "Clay (%)", ylab = "Pielou's evenness")+
          annotate("text", x=19.8, y=0.895,label = "y == -0.0007(x)+0.89", parse=T)+
          annotate("text", x=19.8, y=0.892, label = "R^2 == 0.08", parse=T)+
          annotate("text", x=19.8, y=0.889, label = "p-val == 0.03", parse=T)+
     geom_smooth(method='lm')+labs(title = "I")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

#. 2. BACTERIAL PIELOU'S INDEX - P
P_pie <- lm(Pielou ~ P_ppm, data=map.div)
ggplot(P_pie, aes(x = P_ppm, y = Pielou)) + geom_point() + stat_smooth(method = 'lm')
summary(P_pie)
(j <- ggscatter(map.div, x = "P_ppm", y = "Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "P (ppm)", ylab = "Pielou's evenness")+
          annotate("text", x=40, y=0.897, label = "y == 7.952e-05(x)+0.8712", parse=T)+
          annotate("text", x=40, y=0.894, label = "R^2 == 0.174", parse=T)+
          annotate("text", x=40, y=0.891, label = "p-val == 0.002", parse=T)+
     geom_smooth(method='lm')+labs(title = "J")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

### BACTERIAL PIELOU'S EVENNESS AND NEMATODES ###

#. 3. BACTERIAL PIELOU'S INDEX - Tylenchs
bac.pie_tyl <- cor.test(map.div$Tylenchs, map.div$Pielou, method ="pearson")
cor.test(map.div$Tylenchs, map.div$Pielou)
# linear regression
lm.bac.pie_tyl <- lm(Pielou ~ Tylenchs, data = map.div) # OK
ggplot(lm.bac.pie_tyl, aes(x = Tylenchs, y = Pielou)) + geom_point() + stat_smooth(method = 'lm')
summary(lm.bac.pie_tyl)
gvlma(lm.bac.pie_tyl)
# Do plot
TY_pie <- lm(Pielou ~ Tylenchs, data=map.div)
summary(TY_pie)
ggplot(TY_pie, aes(x = Tylenchs, y = Pielou)) + geom_point() + stat_smooth(method = 'lm')
(k <- ggscatter(map.div, x = "Tylenchs", y = "Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "Tylenchus (per 100 g soil)", ylab = "Pielou's evenness")+
          annotate("text", x=190, y=0.896, label = "y ==  -5.204e-05(x)+0.882", parse=T)+
          annotate("text", x=180, y=0.893, label = "R^2 ==  0.11", parse=T)+
          annotate("text", x=180, y=0.890, label = "p-val == 0.01", parse=T)+
     geom_smooth(method='lm')+labs(title = "K")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=15,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 3. BACTERIAL PIELOU'S INDEX - BacterialFeeders
bac.pie_bacfeed <- cor.test(map.div$BacterialFeeders, map.div$Pielou, method ="pearson")
cor.test(map.div$BacterialFeeders, map.div$Pielou)
# linear regression
lm.bac.pie_bfeed <- lm(Pielou ~ BacterialFeeders, data = map.div) # OK
ggplot(lm.bac.pie_bfeed, aes(x = BacterialFeeders, y = Pielou)) + geom_point() + stat_smooth(method = 'lm')
summary(lm.bac.pie_bfeed)
gvlma(lm.bac.pie_bfeed)
# Do plot
BF_pie <- lm(Pielou ~ BacterialFeeders, data=map.div)
summary(BF_pie)
ggplot(BF_pie, aes(x = BacterialFeeders, y = Pielou)) + geom_point() + stat_smooth(method = 'lm')
(l <- ggscatter(map.div, x = "BacterialFeeders", y = "Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "Rhabditidae (per 100 g soil)", ylab = "Pielou's evenness")+
          annotate("text", x=1500, y=0.868, label = "y ==  6.453e-06(x)+0.87", parse=T)+
          annotate("text", x=1500, y=0.865, label = "R^2 ==  0.10", parse=T)+
          annotate("text", x=1500, y=0.862, label = "p-val == 0.01", parse=T)+
     geom_smooth(method='lm')+labs(title = "L")+
      theme(axis.text.x=element_text(size=15), 
      axis.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=15,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

############### PLOT LINEAR REGRESSION OF BACTERIAL AND ARCHAEAL ALPHA DIVERSITY ##################
################### AND SOIL PARAMETERS AND NEMATODE/OTHER TROPHIC LEVELS #########################
plot <- grid.arrange(a, b, c, d, e, f, g, h, i, j, k, l,
             ncol = 4, nrow = 3, 
             layout_matrix = rbind(c(1,2,3,4), c(5,6,7,8), c(9,10,11,12)))
ggsave("Supplementary Figure_S4.tiff",
       plot, device = "tiff",
       width = 16, height = 11, 
       units= "in", dpi = 600)
###################################################################################################

### BACTERIAL SHANNON ###

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

### FUNGAL SHANNON DIVERSITY ###
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

### FUNGAL PIELOU'S EVENNESS ###
# 1. FUNGAL PIELOU - K CONTENT
fg.pie_K <- lm(fg.Pielou ~ K_ppm, data=map.div)
gvlma(fg.pie_K)
ggplot(fg.pie_K, aes(x = K_ppm, y = fg.Pielou)) + geom_point() + stat_smooth(method = 'lm')
summary(fg.pie_K)
(c <- ggscatter(map.div, x = "K_ppm", y = "fg.Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "K (ppm)", ylab = "Pielou's evenness")+
          annotate("text", x=225, y=0.75, label = "y == -0.0003(x)+0.66", parse=T)+
          annotate("text", x=225, y=0.73, label = "R^2 == 0.07", parse=T)+
          annotate("text", x=225, y=0.71, label = "p-val == 0.04", parse=T)+
     geom_smooth(method='lm')+labs(title = "C")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 2. FUNGAL PIELOU - Mg CONTENT
fg.pie_Mg <- lm(fg.Pielou ~ Mg_ppm, data=map.div)
gvlma(fg.pie_Mg)
ggplot(fg.pie_Mg, aes(x = Mg_ppm, y = fg.Pielou)) + geom_point() + stat_smooth(method = 'lm')
summary(fg.pie_Mg)
(d <- ggscatter(map.div, x = "Mg_ppm", y = "fg.Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "Mg (ppm)", ylab = "Pielou's evenness")+
          annotate("text", x=230, y=0.72, label = "y == -0.0003(x)+0.65", parse=T)+
          annotate("text", x=230, y=0.7, label = "R^2 == 0.08", parse=T)+
          annotate("text", x=230, y=0.68, label = "p-val == 0.02", parse=T)+
     geom_smooth(method='lm')+labs(title = "D")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 3. FUNGAL PIELOU - SAND
fg.pie_sand <- lm(fg.Pielou ~ sand_percent, data=map.div)
gvlma(fg.pie_sand)
ggplot(fg.pie_sand, aes(x = sand_percent, y = fg.Pielou)) + geom_point() + stat_smooth(method = 'lm')
summary(fg.pie_sand)
(a <- ggscatter(map.div, x = "sand_percent", y = "fg.Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "Sand (%)", ylab = "Pielou's evenness")+
          annotate("text", x=40, y=0.72, label = "y == 0.0016(x)+0.52", parse=T)+
          annotate("text", x=40, y=0.7, label = "R^2 == 0.08", parse=T)+
          annotate("text", x=40, y=0.68, label = "p-val == 0.02", parse=T)+
     geom_smooth(method='lm')+labs(title = "A")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

# 4. FUNGAL PIELOU - SILT
fg.pie_silt <- lm(fg.Pielou ~ silt_percent, data=map.div)
ggplot(fg.pie_silt, aes(x = silt_percent, y = fg.Pielou)) + geom_point() + stat_smooth(method = 'lm')
summary(fg.pie_silt)
gvlma(fg.pie_silt)
(b <- ggscatter(map.div, x = "silt_percent", y = "fg.Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "Silt (%)", ylab = "Pielou's evenness")+
          annotate("text", x=40, y=0.72, label = "y == -0.002(x)+0.67", parse=T)+
          annotate("text", x=40, y=0.7, label = "R^2 == 0.07", parse=T)+
          annotate("text", x=40, y=0.68, label = "p-val == 0.03", parse=T)+
     geom_smooth(method='lm')+labs(title = "B")+
      theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_blank(),
       axis.ticks.y = element_blank(),
       axis.title.y = element_blank(),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"),
       plot.title = element_text(size = rel(2)),
       axis.title=element_text(size=18,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank()))

############### PLOT LINEAR REGRESSION OF FUNGAL ALPHA DIVERSITY ##################
################### AND SOIL PARAMETERS AND NEMATODE/OTHER TROPHIC LEVELS #########
plot <- grid.arrange(a, b, c, d,
             ncol = 2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,4)))
ggsave("Supplementary Figure S5.tiff",
       plot, device = "tiff",
       width = 8, height = 8, 
       units= "in", dpi = 600)
###################################################################################################

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
boxplot(nema.total.SqrtCount ~ Site, data = map.div) # there are no outliers
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
 labs(title="A. Nematodes", y=expression(bold(atop("", 
               atop(textstyle("Absolute abundances"),
                    atop(textstyle("(individu per 100 g soil)")))))))+
 #geom_text(data=new.nema.count_site.summarized,aes(x=Site,y=20+max.nema.count,label=groups),vjust=0)+
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=10,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_blank(),
       axis.title.y=element_blank(),
       #axis.title.y=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

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
pval.nema.root <- as.data.frame(HSD_TC_nema_root$comparison)
write.csv(pval.nema.root, file = "pval.nema.root.csv")

# add significance letters from HSD.test into box plot
TC_nema_root.sum1 <- map.div %>% group_by(Rootstock) %>% summarize(max.count=max(nema.total.Count))
HSD_TC_nema_root <- HSD.test(TC_nema_root, "Rootstock", group = TRUE, console = TRUE)
hsd_TC_nema_root = HSD_TC_nema_root$groups
class(hsd_TC_nema_root)
hsd_TC_nema_root$Rootstock <- rownames(hsd_TC_nema_root)
TC_nema_root.sum2=left_join(hsd_TC_nema_root,TC_nema_root.sum1, by='Rootstock')  
# Do Plot
(nema.ab.root <- ggplot(map.div, aes(x=Rootstock, y=nema.total.Count))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
 labs(y=expression(bold(atop("", 
               atop(textstyle("Absolute abundances"), 
               atop(textstyle("(individu per 100 g soil)")))))))+
 geom_text(data=TC_nema_root.sum2,aes(x=Rootstock,y=45+max.count,label=TC_nema_root.sum2$groups),vjust=0)+
 theme(axis.text.x=element_text(size=15,angle=49,hjust =0.9),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.x = element_blank(),
       axis.title.y = element_blank(),
       #axis.title=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid = element_blank()))

# 3. Compare total absolute nematodes abundances among cultivars # not significant
TC_nema_cult <- lm(nema.total.SqrtCount ~ Scion, data=map.div, na.action=na.exclude)
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
leveneTest(nema.total.SqrtCount ~ Scion, data=map.div, na.action=na.exclude) # variances among group are homogenous

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
pval.OL.site <- as.data.frame(HSD_OL_site$comparison)
write.csv(pval.OL.site, file = "pval.OL.site.csv")

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
 labs(title ="B. Oligochaetes")+
 geom_text(data=OL_site.sum2,aes(x=Site,y=2+max.count,label=OL_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=10,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.y = element_blank(),
       axis.title.x=element_blank(),
       #axis.title.x=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

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
boxplot(sqrt.OL ~ Rootstock, data = map.div) # there are no outliers
plot(density(OL_site_resids))
qqnorm(OL_site_resids)
qqline(OL_site_resids) 
hist(OL_site_resids)
# Do Plot
OL_root <- ggplot(map.div, aes(x=Rootstock, y=Oligochaetes))+
 geom_boxplot() +
 geom_point() +
 theme_bw()+
#geom_text(data=OL_site.sum2,aes(x=Site,y=2+max.count,label=OL_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_text(size=15,angle=49,hjust =0.9),
       axis.text.y = element_text(size = 14),
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       #axis.title.x=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid = element_blank())

# 6. Compare absolute abundance of oligochaetes among scions # not significant
OL.cul <- lm(sqrt.OL ~ Scion, data=map.div, na.action=na.exclude)
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
leveneTest(sqrt.OL ~ Scion, data=map.div, na.action=na.exclude)

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
HSD_MF_site <- HSD.test(MF_site, "Site", group = F)
comp.df <- as.data.frame(HSD_MF_site$comparison)
comp.df
pval.MF.site <- as.data.frame(HSD_MF_site$comparison)
write.csv(pval.MF.site , file = "pval.MF.site.csv")

p.val <- rownames_to_column(comp.df, var = "comparison")
p.val 
library(rcompanion)
letter <- cldList(pvalue ~ comparison,
        data = p.val,
        threshold  = 0.05)
letter
# with multcompletters
p.val = HSD_MF_site$comparison[,'pvalue', drop=F]
p.val
names(p.val)[1] <- "P.adj"
p.val
library(multcompView)
multcompLetters(p.val,
                compare="<",
                threshold=0.05,
                Letters=letters,
                reversed = FALSE)



# add significance letters from HSD.test into box plot
MF_site.sum1 <- map.div %>% group_by(Site) %>% summarize(max.count=max(MycorrhizalFungi))
HSD_MF_site <- HSD.test(MF_site, "Site", group = TRUE, console = TRUE)
hsd_MF_site = HSD_MF_site$groups
class(hsd_MF_site)
hsd_MF_site$Site <- rownames(hsd_MF_site)
MF_site.sum2=left_join(hsd_MF_site,MF_site.sum1, by='Site')  
MF_site.sum2

# Do Plot
MF_site <- ggplot(map.div, aes(x=Site, y=MycorrhizalFungi))+
 geom_boxplot() +
 geom_point() +
 scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)) +
 theme_bw()+
 labs(title = "C. Mycorrhizal fungi")+
 geom_text(data=MF_site.sum2,aes(x=Site,y=10+max.count,label=MF_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_text(size=15), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=10,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       #axis.title.x=element_text(size=14,face="bold"),
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
#geom_text(data=OL_site.sum2,aes(x=Site,y=2+max.count,label=OL_site.sum2$groups),vjust=0)+
 theme(axis.text.x=element_text(size=15,angle=49,hjust =0.9), 
       axis.text.y = element_text(size = 14),
       strip.text.x = element_text(size=18,colour = "black", face = "bold"), 
       strip.text.y = element_text(size=18, face = 'bold'),
       plot.title = element_text(size = rel(2)),
       axis.title.y = element_blank(),
       axis.title.x = element_blank(),
       #axis.title.x=element_text(size=14,face="bold"),
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())

# 9. Compare absolute mycorrhizal fungi abundance among cultivars # not significant
MF_cul <- lm(log.MF ~ Scion, data=map.div, na.action=na.exclude)
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
leveneTest(log.MF ~ Scion, data=map.div, na.action=na.exclude)

######### PLOT NEMATODE, OLIGO, AND MF ABSOLUTE ABUNDANCES ###########################################
#plot <- grid.arrange(nema.ab.site, OL_site, MF_site,                                    
#             nema.ab.root, OL_root, MF_root,                              
#            ncol = 3, nrow = 2, 
#             layout_matrix = rbind(c(1,2,3), c(4,5,6)))
#ggsave("Figure 2.eps",plot, device = "eps",width = 20, height = 10,units= "in", dpi = 600)

plot1 <- ggarrange(nema.ab.site, OL_site, MF_site,                                    
             ncol = 3, nrow = 1)
plot2 <- ggarrange(nema.ab.root, OL_root, MF_root,                                    
             ncol = 3, nrow = 1)
a <- annotate_figure(plot1,
                bottom = text_grob("Site", face = "bold", size = 20))
b <- annotate_figure(plot2,
                bottom = text_grob("Rootstock", face = "bold", size = 20))             
plot3 <- ggarrange(a,b, ncol=1, nrow = 2)
c <- annotate_figure(plot3,
                     left = text_grob("Absolute abundances (individu per 100 g soil)", face = "bold", size = 20, rot = 90))
ggsave("Figure_2.eps",
       c, device = "eps",
       width = 20, height = 10, 
       units= "in", dpi = 600)

################################################################################################################

# 3. ANOVA TEST FOR NEMATODE ALPHA DIVERSITY
 
### NEMATODE RICHNESS ####

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

# 2. compare nematode richness among rootstocks #not significant
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

### NEMATODE SHANNON ####

# 1. compare nematode shannon index among sites #not significant
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

# 2. compare nematode shannon index among rootstocks #not significant
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

#### NEMATODE PIELOU'S EVENNESS ####

# 1. compare nematode pielou's evenness index among sites #not significant
AovNema_pie_site <- lm(nema.Pielou ~ Site, data=map.div, na.action=na.exclude)
drop1(AovNema_pie_site,~.,test="F") # type III SS and F Tests
summary(AovNema_pie_site)
# Check normality
pie.Nema_resids <- residuals(AovNema_pie_site)
pie.Nema_preds <- predict(AovNema_pie_site)
plot(pie.Nema_resids ~ pie.Nema_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
shapiro.test(pie.Nema_resids) # normal
qqnorm(pie.Nema_resids)
qqline(pie.Nema_resids) 
leveneTest(nema.Pielou ~ Site, data=map.div, na.action=na.exclude) # homogen
gvlma(AovNema_pie_site)

# 2. compare nematode Pielou's evenness index among rootstocks #not significant
AovNema_pie_root <- lm(nema.Pielou ~ Rootstock , data=map.div, na.action=na.exclude)
summary(AovNema_pie_root)
drop1(AovNema_pie_root,~.,test="F") # type III SS and F Tests
# Check normality
pie.Nema_resids <- residuals(AovNema_pie_root)
pie.Nema_preds <- predict(AovNema_pie_root)
plot(pie.Nema_resids ~ pie.Nema_preds, xlab = "Predicted Values", ylab = "Residuals", asp=.5)
shapiro.test(pie.Nema_resids) # normal
leveneTest(nema.Pielou ~ Rootstock, data=map.div, na.action=na.exclude) # homogen
gvlma(AovNema_pie_root)

# ANOVA test to check the differences of nematoda alpha diversity among cultivars
# 1. nematode - Richness
nema_rich_cul <- lm(nema.Richness ~ cultivar, data=map.div, na.action=na.exclude)
summary(nema_rich_cul) # not significant
# 2. nematode - Shannon
nema_sha_cul <- lm(nema.Shannon ~ cultivar, data=map.div, na.action=na.exclude)
summary(nema_sha_cul) # not significant
# 3. nematode - Pielou
nema_pie_cul <- lm(nema.Pielou ~ cultivar, data=map.div, na.action=na.exclude)
summary(nema_pie_cul) # not significant

# Check nematode prevalence
nema_PA <- 1*(nema.t>0)
sum_nemaPA <- rowSums(nema_PA)
sum_nemaPA





### BACTERIAL SHANNON INDEX ###

# 1. Bacterial Shannon and nematode Tylenchs
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



#### FUNGI ###

### FUNGAL RICHNESS ####

# 3. Fungal Richness and nematode Tylenchs # not significant
fg.rich_tyl <- cor.test(map.div$Tylenchs, map.div$fg.Richness, method ="pearson")
cor.test(map.div$Tylenchs, map.div$fg.Richness)

### FUNGAL SHANNON INDEX ####

# 4. Fungal Shannon and nematode Tylenchs # not significant
fg.sha_tyl <- cor.test(map.div$Tylenchs, map.div$fg.Shannon, method ="pearson")
cor.test(map.div$Tylenchs, map.div$fg.Shannon)

### FUNGAL PIELOU'S INDEX ####
# not significant


### PEARSON CORRELATION AND LINEAR REGRESSION TEST OF NEMATODES AND BACTERIAL AND FUNGAL RALPHA DIVERSITY ###

# 1. Bacterial Shannon and nematode Shannon
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


# 2. Bacterial Pielou and nematode Pielou
BC.pie_NM.pie<- lm(Pielou ~ nema.Pielou, data = map.div)
summary(BC.pie_NM.pie)
ggplot(BC.pie_NM.pie, aes(x = nema.Pielou, y = Pielou)) + geom_point() + stat_smooth(method = 'lm')
plot <- ggscatter(map.div, x = "nema.Pielou", y = "Pielou", 
          add = "reg.line",conf.int = TRUE, xlab = "Nematode Pielou's evenness", ylab = "Bacterial/archaeal Pielou's evenness")+
          annotate("text", x=0.2, y=0.866, label = "y == -0.01(x)+0.88", parse=T)+
          annotate("text", x=0.2, y=0.863, label = "R^2 == 0.11", parse=T)+
          annotate("text", x=0.2, y=0.860, label = "p-val == 0.01", parse=T)+
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

tiff("Supplementary Figure S6.tiff", units="in", width=5, height=5, res=600)
plot
dev.off()



############################################################################################
######################### BACTERIAL, FUNGAL, and NEMATODE BETA DIVERSITY ###################
############################################################################################

# 1. CALCULATE BETA DIVERSITY (PCoA PLOT) FOR BACTERIA
# dissimilarity indices for community ecologist to make a distance structure (Bray-Curtis distance between samples)
otu_dist <- vegdist(t(otu), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoa <- cmdscale(otu_dist, eig=T)
env <- map[,c(6,7,11:22, 24:36)]
# scores of PC1 and PC2
ax1.scores=otu_pcoa$points[,1]
ax2.scores=otu_pcoa$points[,2] 
set.seed(100)
env_fit <- envfit(otu_pcoa, env, na.rm=TRUE)
env_fit
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
myCol <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
# PCoA Plot by site
bac.pcoa <- ggplot(data = map2, aes(x=ax1.scores, y=ax2.scores))+
  theme_bw()+
 geom_point(data = map2, aes(x = ax1.scores, y = ax2.scores, color=Site),size=5,shape=20, alpha=0.5)+
 scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8','#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'))+
  #geom_text(aes(ax1.scores, y=ax2.scores, label = Site),size=3)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  geom_segment(data=bac.scores4, 
               aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2),
               arrow = arrow(length = unit(0.3, "cm")),
               colour = "grey")+
  geom_text_repel(data = bac.scores4, 
                  aes(x = mult*Dim1, y = mult*Dim2, label = Variable),
                  size = 3,fontface="bold", 
                  position=position_jitter(width=0.03,height=0.001))+
        coord_fixed() + 
 labs(title = "A. Bacteria/archaea")+
 theme(legend.position="none",
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=12),
       legend.title = element_text(size = 12),
       legend.spacing.x = unit(0.05, 'cm'))

# PCoA Plot by rootstock
bac_root.pcoa <- ggplot(data = map2, aes(x=ax1.scores, y=ax2.scores))+
  theme_bw()+
  geom_point(aes(ax1.scores, y=ax2.scores,color=Rootstock),size=5,shape=20, alpha=0.5)+
 scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6'))+
 #geom_text(aes(ax1.scores, y=ax2.scores,color=Rootstock, label = Rootstock),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2,3)*100,"% var. explained", sep=""))+
  geom_segment(data=bac.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = bac.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
        coord_fixed() + 
 labs(title = "A. Bacteria/archaea")+
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
#otuITS <- otuITS %>% remove_rownames %>% column_to_rownames(var="OTU_ID")
otu_distITS <- vegdist(t(otuITS), method='bray')
# CMD/classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
otu_pcoaITS <- cmdscale(otu_distITS, eig=T)
env <- map[,c(6,7,11:22, 24:36)]
# scores of PC1 and PC2
ax1ITS.scores=otu_pcoaITS$points[,1] 
ax2ITS.scores=otu_pcoaITS$points[,2] 
set.seed(100)
env_fitITS <- envfit(otu_pcoaITS, env, na.rm=TRUE)
env_fitITS
ax1ITS <- otu_pcoaITS$eig[1]/sum(otu_pcoaITS$eig)
ax2ITS <- otu_pcoaITS$eig[2]/sum(otu_pcoaITS$eig)
map3=cbind(map2,ax1ITS.scores,ax2ITS.scores)
head(map3)
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
# PCoA Plot by site
fg.pcoa <- ggplot(data = map3, aes(x=ax1ITS.scores, y=ax2ITS.scores)) +
  theme_bw()+
  geom_point(aes(ax1ITS.scores, y=ax2ITS.scores,color=Site),size=5,shape=20, alpha=0.5)+
 scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8','#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'))+
 #geom_text(aes(ax1ITS.scores, y=ax2ITS.scores,color=Site, label = Site),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  coord_fixed()+
 labs(title = "B. Fungi")+
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
       legend.spacing.x = unit(0.05, 'cm'))
#+scale_color_discrete(breaks=sort(as.numeric(map2$Site)))

# PCoA Plot by rootstocks
fg_root.pcoa <- ggplot(data = map3, aes(x=ax1ITS.scores, y=ax2ITS.scores)) +
  theme_bw()+
  geom_point(aes(ax1ITS.scores, y=ax2ITS.scores,color=Rootstock),size=5,shape=20,alpha=0.5)+
 scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6'))+
 #geom_text(aes(ax1ITS.scores, y=ax2ITS.scores,color=Rootstock, label = Rootstock),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1ITS,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2ITS,3)*100,"% var. explained", sep=""))+
  coord_fixed()+
 labs(title = "B. Fungi")+
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
map4=cbind(map3,ax1.nema.scores,ax2.nema.scores)
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
#PCoA Plot by site
nm.pcoa <- ggplot(data = map4, aes(x=ax1.nema.scores, y=ax2.nema.scores))+
  theme_bw()+
  geom_point(aes(ax1.nema.scores, y=ax2.nema.scores,color=Site),size=5,shape=20, alpha=0.5)+
 scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8','#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'))+
 #geom_text(aes(ax1.nema.scores, y=ax2.nema.scores,color=Site, label = Site),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1.nema,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2.nema,3)*100,"% var. explained", sep=""))+
  geom_segment(data=nm.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = nm.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
        coord_fixed() + 
 labs(title = "C. Nematodes")+
 theme(legend.position = "right",
       legend.box = "vertical",
       plot.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       plot.title = element_text(size = rel(1), face="bold"),
       axis.text=element_text(size=10), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size=8),
       legend.title = element_text(size = 12),
       legend.spacing.x = unit(0.05, 'cm'))


#PCoA Plot by rootstock
nm.pcoa.root <- ggplot(data = map4, aes(x=ax1.nema.scores, y=ax2.nema.scores))+
  theme_bw()+
  geom_point(aes(ax1.nema.scores, y=ax2.nema.scores,color=Rootstock),size=5,shape=20, alpha=0.5)+
 scale_colour_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8','#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'))+
 #geom_text(aes(ax1.nema.scores, y=ax2.nema.scores,color=Site, label = Site),size=4)+
  scale_x_continuous(name=paste("PCoA1: ",round(ax1.nema,3)*100,"% var. explained", sep=""))+
  scale_y_continuous(name=paste("PCoA2: ",round(ax2.nema,3)*100,"% var. explained", sep=""))+
  geom_segment(data=nm.scores4, aes(x=0, xend=mult*Dim1, y=0, yend=mult*Dim2), arrow = arrow(length = unit(0.3, "cm")), colour = "grey")+
  geom_text_repel(data = nm.scores4, aes(x = mult*Dim1, y = mult*Dim2, label = Variable), size = 3,fontface="bold",position=position_jitter(width=0.03,height=0.001))+
        coord_fixed() + 
 labs(title = "C. Nematodes")+
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
       legend.spacing.x = unit(0.05, 'cm'))
 


#+scale_color_discrete(guide = guide_legend(title.position = "top", nrow = 2),breaks=sort(as.numeric(map2$Site)))

######################################### PCOA PLOT #############################################################
#install.packages("lemon")
#library(lemon)
library(gridExtra)
plot1 <- ggarrange(bac.pcoa, fg.pcoa, nm.pcoa,common.legend=T, legend="bottom",nrow=1, ncol=3, align = "hv")
plot2 <- ggarrange(bac_root.pcoa, fg_root.pcoa, nm.pcoa.root, common.legend = T, legend="bottom", nrow = 1, ncol = 3, align = "hv")
#p <- ggarrange(plot1,plot2, nrow = 2, align = "hv")
ggsave("SupplementaryFigure_S7.tiff",
       plot1, device = "tiff",
       width = 10, height = 5, 
       units= "in", dpi = 600)

ggsave("SupplementaryFigure_S8.tiff",
       plot2, device = "tiff",
       width = 10, height = 5, 
       units= "in", dpi = 600)

################# PERMANOVA TO TEST ANY INFLUENCE OF SITE, ROOTSTOCK, AND CULTIVAR TO BACTERIAL, FUNGAL, AND NEMATODE BETA DIVERSITY ################################
### SITE ###
# 1. Site to bacterial beta diversity
adonis(otu_dist~map2$Site)
adonis(otu_dist~map2$Site, strata = map2$rootstock) #compare different sites at same rootstock
adonis(otu_dist~map2$rootstock, strata = map2$Site) # compare different rootstock at same site

# 2. Site to fungal beta diversity
fg.otudist_site=adonis(otu_distITS~map2$Site) # p-val=0.001
adonis(otu_distITS~map2$Site, strata = map2$rootstock) # compare different sites at same rootstock
adonis(otu_distITS~map2$rootstock, strata = map2$Site) # compare different rootstock at same site
adonis(otu_distITS~map2$cultivar, strata = map2$Site) # compare different cultivar at same site

# 3. Site to nematode beta diversity
adonis(otu_dist_nema~map2$Site) # p-val=0.018
adonis(otu_dist_nema~map2$Site, strata = map2$rootstock)
adonis(otu_dist_nema~map2$rootstock, strata = map2$Site) # compare different rootstock at same site
adonis(otu_dist_nema~map2$cultivar, strata = map2$Site) # compare different cultivar at same site


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

############################################################################################################
###### PROCRUSTES - PROTEST TO IDENTIFY THE CONCORDANCE OF BACTERIAL, FUNGAL, AND NEMATODE PCOA PLOTS ######
###########################################################################################################

install.packages("ade4")
library(ade4)

# 1. Classical MDS 16S
otu_dist <- vegdist(t(otu), method='bray')
#Cailliez correction
add <-  !(is.euclid(otu_dist)) # false #to determine if the distance matrix is euclidean in the sense of not generating negative eigenvalues
pcoa.bac <- cmdscale(otu_dist, k = nrow(t(otu))-1, eig = TRUE, add = add) # the add = add argument to the cmdscale function will add a small value to all distances to prevent negative eigenvalues

# 2. Classical MDS ITS
otu_distITS <- vegdist(t(otuITS), method='bray')
#Cailliez correction
addITS <-  !(is.euclid(otu_distITS)) #true
otu_pcoaITS <- cmdscale(otu_distITS, eig=T)

# 3. Classical MDS Nematode (11 group)
otu_dist_nema <- vegdist(t(sqrt.nema.t), method='bray')
add.nema <-  !(is.euclid(otu_dist_nema)) #true
otu_pcoa_nema <-  cmdscale(otu_dist_nema, eig=T)
otu_pcoa_nema

# 1. Procrustest and Protest bacterial and fungal PCoA
bac.pcoa_fg.pcoa <- procrustes(pcoa.bac, otu_pcoaITS, scale = TRUE, symmetric = TRUE) #Procrustes sum of squares:0.90
summary(bac.pcoa_fg.pcoa)
protest(pcoa.bac, otu_pcoaITS, scale = TRUE, symmetric = TRUE, permutations = 999)
#Procrustes Sum of Squares (m12 squared):        0.9083 
#Correlation in a symmetric Procrustes rotation: 0.3029 
#Significance:  0.001
###there is significxant correlation between the two ordinations###

# 2. Procrustes and Protest bacterial and nematode PCoA
bac.pcoa_nema.pcoa <- procrustes(pcoa.bac, otu_pcoa_nema, scale = TRUE, symmetric = TRUE) #Procrustes sum of squares:0.93
summary(bac.pcoa_nema.pcoa)
protest(pcoa.bac, otu_pcoa_nema, scale = TRUE, symmetric = TRUE, permutations = 999)
#Procrustes Sum of Squares (m12 squared):        0.9356 
#Correlation in a symmetric Procrustes rotation: 0.2538 
#Significance:  0.015
###there is significxant correlation between the two ordinations###

# 3. Procrustes and Protest fungal and nematode PCoA
fg.pcoa_nema.pcoa <- procrustes(otu_pcoaITS, otu_pcoa_nema, scale = TRUE, symmetric = TRUE) #Procrustes sum of squares:0.99
summary(fg.pcoa_nema.pcoa)
protest(otu_pcoaITS, otu_pcoa_nema, scale = TRUE, symmetric = TRUE, permutations = 999)
#Procrustes Sum of Squares (m12 squared):        0.9919
#Correlation in a symmetric Procrustes rotation: 0.08987
#Significance:  0.897
###there is no significxant correlation between the two ordinations###

####################### PRODUCE SAMPLING LOCATION MAP ##############################

devtools::install_github("dkahle/ggmap", force = TRUE)
library(ggmap)
library(ggrepel)
ordinate <- read.table('LONG_LAT.csv', sep=',', header=TRUE)
ordinate1 <- ordinate[,1:4]
register_google(key = "AIzaSyA-H5MDcBPrVwZrEOt7iS3zsxHC_fnZM7c",  
                account_type = "standard")
apporchard_map <- get_map(location = c(lon = -85.67709,lat = 43.09583),zoom = 9,
                          source="google", maptype = "terrain", crop=FALSE, scale=2)
ggmap(apporchard_map)
ordinate2 <- subset(ordinate1[c("1","6","9","12","14","16","18","20","22","24","26","27","28","30","32","36","37","38","41","43"),])

p1 <-  ggmap(apporchard_map) +
 geom_point(data = ordinate2, 
           aes(x = LON, y = LAT, alpha = 0.8), 
           size = 3)
p1 <- p1 + geom_text_repel(data = ordinate2, 
           aes(x = LON, y = LAT, label=Site, alpha = 0.8), 
           size = 6)+
 theme(axis.text=element_text(size=15), 
       axis.title=element_text(size=15, face="bold"),
       legend.text=element_text(size=15),
       legend.title = element_text(size = 12),
       legend.spacing.x = unit(1.0, 'cm'))+
 guides(fill=FALSE, alpha=FALSE, size=FALSE)+
  xlab("Longitude")+ylab("Latitude")

# zoom out Michigan map
# how to enable Google Maps Geocoding API:
# 1. go to Developer Console -> APIs & auth -> APIs
# 2. search for Geocoding and click on Google Maps Geocoding API -> Enable API. Do the same thing for Geolocating

p2 <- ggmap(get_map(location = "michigan", zoom=7))+
 theme(axis.title = element_blank(), 
          axis.text  = element_blank(),
          axis.ticks = element_blank())+
  annotate(geom = "rect", ymax = 43.75, ymin = 42.5, xmax = -84.7, xmin = -86.5, colour = "red", fill = NA)
library(grid)
# save the sampling map plot
tiff("Supplementary Figure S1.tiff", units="in", width=7, height=7, res=600)
sampling_map <- p1 +
 inset(ggplotGrob(p2), 
       xmin = -85.5, 
       xmax = -84.7, 
       ymin = 42.48, 
       ymax = 42.9) 
print(sampling_map)
dev.off()

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
# add map
# Read the metadata (map)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
str(map) # We use number as the name of site and it is integer 
map$Site <- as.factor(map$Site) # I want to change site as factor
colnames(map)[which(names(map) == "rootstock")] <- "Rootstock"
colnames(map)[which(names(map) == "cultivar")] <- "Scion"
map$Site<-as.factor(map$Site)
rownames(map) <- map$sample_code
head(map)
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
#write.table(sort.df.phylum.16S.taxasum, file = 'relabun_phylum_bac.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
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
p.site <- ggplot(data=df.phylum.16S, aes(x=Site, y=Abundance, fill=Phylum))
barplot.16S.site <- p.site + geom_bar(aes(), stat="identity", position="fill") + scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme(legend.position="none") + guides(fill=guide_legend(nrow=5))+
 labs(title="A. Bacteria/archaea",y= "Relative Abundance")+
 theme(plot.title = element_text(size = rel(1.5), face="bold"),
       axis.line = element_line(size=0.5, colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
       axis.text=element_text(size=9),
       axis.text.x = element_text(hjust = 1),
       axis.title=element_text(size=10,face="bold"),
       legend.text=element_text(size = 8),
       legend.title = element_text(size=10),
       panel.grid = element_blank(), 
       panel.background = element_blank())
       
p.root <- ggplot(data=df.phylum.16S, aes(x=Rootstock, y=Abundance, fill=Phylum))
barplot.16S.root <- p.root + geom_bar(aes(), stat="identity", position="fill") + scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme(legend.position="none") + guides(fill=guide_legend(nrow=5))+
 labs(y= "Relative Abundance")+
 theme(plot.title = element_text(size = rel(1.5), face="bold"),
       axis.line = element_line(size=0.5, colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
       axis.text=element_text(size=9),
       axis.text.x = element_text(angle=45,hjust = 1),
       axis.title=element_text(size=10,face="bold"),
       legend.text=element_text(size = 8),
       legend.title = element_text(size=10),
       panel.grid = element_blank(), 
       panel.background = element_blank())

################ SAVE BAR PLOT #########################################
bp1 <- ggarrange(barplot.16S.site,barplot.16S.root, nrow = 2, ncol = 1, common.legend = T, legend = "bottom", align = "hv")
bp2 <- ggarrange(barplot.ITS.site,barplot.ITS.root, nrow = 2, ncol = 1, common.legend = T, legend = "bottom", align = "hv")
bp <- ggarrange(bp1,bp2, nrow=1,ncol = 2, align = "hv")
ggsave("Figure_3.eps",
       bp, device = "eps",
       width = 12, height = 7, 
       units= "in", dpi = 600)

# 2. FUNGAL COMPOSITION
# READ TAXONOMY FILE ("consensus_taxonomy_ITS.txt") FROM CONSTAX OUTPUT AND MERGE IT INTO "otuITS" FILE ("OTU_rarefied_ITS.txt")
otuITS <- read.table(file = "OTU_rarefied_ITS.txt", sep='\t', header=T)
fg.taxonomy = read.table("consensus_taxonomy_ITS.txt",
                      header = T, sep = "\t")
dim(fg.taxonomy)
fg.otu.tax <- merge(otuITS, fg.taxonomy, by.x =c("OTU_ID"), by.y = c("OTU_ID"))
head(fg.otu.tax)
dim(fg.otu.tax)
# write.table(fg.otu.tax, file = "combined_otu_tax_ITS.txt", sep = '\t', col.names = TRUE, row.names = FALSE)
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
# add map
# Read the metadata (map)
map <- read.table('clean_map_data.csv', sep=',', header=TRUE)
str(map) # We use number as the name of site and it is integer 
map$Site <- as.factor(map$Site) # I want to change site as factor
colnames(map)[which(names(map) == "rootstock")] <- "Rootstock"
colnames(map)[which(names(map) == "cultivar")] <- "Scion"
map$Site<-as.factor(map$Site)
rownames(map) <- map$sample_code
head(map)
# make phyloseq map
phyloseq_map <- sample_data(map)
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
#write.table(sort.df.phylum.ITS.taxasum, file = 'relabun_phylum_fg.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
df.phylum.ITS <- psmelt(phylum.ITS)
head(df.phylum.ITS)
dim(df.phylum.ITS)
df.phylum.ITS$Phylum <- as.character(df.phylum.ITS$Phylum)
df.phylum.ITS$Phylum[df.phylum.ITS$Abundance < 0.01] <- "Other"
df.phylum.ITS$Phylum[is.na(df.phylum.ITS$Phylum)] <- "Other"
#barplot ITS
f.site <- ggplot(data=df.phylum.ITS, aes(x=Site, y=Abundance, fill=Phylum))
barplot.ITS.site <- f.site + geom_bar(aes(), stat="identity", position="fill") + scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme(legend.position="none") + guides(fill=guide_legend(nrow=5))+
 labs(title="B. Fungi", y= "Relative Abundance")+
 theme(plot.title = element_text(size = rel(1.5), face="bold"),
       axis.line = element_line(size=0.5, colour = "black"),
       axis.text=element_text(size=9),
       axis.text.x = element_text(hjust = 1),
       axis.title.y = element_blank(),
       axis.title=element_text(size=10,face="bold"),
       legend.text=element_text(size = 8),
       legend.title = element_text(size=10),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(), 
       panel.background = element_blank())

f.root <- ggplot(data=df.phylum.ITS, aes(x=Rootstock, y=Abundance, fill=Phylum))
barplot.ITS.root <- f.root + geom_bar(aes(), stat="identity", position="fill") + scale_fill_manual(values=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme(legend.position="none") + guides(fill=guide_legend(nrow=5))+
 labs(y= "Relative Abundance")+
 theme(plot.title = element_text(size = rel(1.5), face="bold"),
       axis.line = element_line(size=0.5, colour = "black"),
       axis.text=element_text(size=9),
       axis.text.x = element_text(angle=45,hjust = 1),
       axis.title.y = element_blank(),
       axis.title=element_text(size=10,face="bold"),
       legend.text=element_text(size = 8),
       legend.title = element_text(size=10),
       panel.grid = element_blank(), 
       panel.background = element_blank())

# Make rarefaction curve for bacterial and archaeal, and fungal OTUs
# I constructed the rarefaction curve using 'Phyloseq' package, thus I need to produce phyloseq object for bacteria and fungi.
setEPS()
postscript("Supplementary Figure S2.eps", height = 5, width = 10)
par(mfrow=c(1,2))
bac.rarecurve <- rarecurve(t(otu_table(PHYL_16S)), # PHYL is phyloseq object of bacteria and archaea
                           step=50, cex=0.5,
                           xlab = "Reads", 
                           ylab = "Bacterial and archaeal OTUs")
title("A", adj=0)
fg.rarecurve <- rarecurve(t(otu_table(PHYL_ITS)), # physeq is phyloseq object of fungi
                          step=50, cex=0.5, 
                          xlab = "Reads", 
                          ylab = "Fungal OTUs")
title("B", adj=0)
dev.off()
while (!is.null(dev.list())) dev.off() #make sure the plotting device is off
graphics.off()



###############################################################
######################### CORE MICROBIOTA #####################
###############################################################

# OCCUPANCY VS ABUNDANCE
# Occupancy
# 1. Bacteria
otu_PA <- 1*((otu>0)==1)
otu_PA <- otu_PA[rowSums(otu_PA)>0,]
Occ <- rowSums(otu_PA)/ncol(otu_PA)
class(Occ)
class(Occ.1)
df.Occ <- as.data.frame(Occ)
head(df.Occ)
df.Occ=rownames_to_column(df.Occ, var = "OTU")
dim(df.Occ)
# 2. Fungi
head(otuITS)
dim(otuITS)
set.seed(13)
head(sort(colSums(otuITS, na.rm = FALSE, dims = 1), decreasing = FALSE))
head(sort(rowSums(otuITS, na.rm = FALSE, dims = 1), decreasing = FALSE))
otuITS_PA <- 1*((otuITS>0)==1)
otuITS_PA <- otuITS_PA[rowSums(otuITS_PA)>0,]
OccITS <- rowSums(otuITS_PA)/ncol(otuITS_PA)
df.OccITS <- as.data.frame(OccITS)
df.OccITS=rownames_to_column(df.OccITS, var = "OTU")
head(df.OccITS)
dim(df.OccITS)

# Taxonomy
# 1. Bacteria
tax_16S = read.csv("16S_TAX.csv", sep=',', header=T)
rownames(tax_16S) <- rownames(otu)
head(tax_16S)
dim(tax_16S)
tax_16S=rownames_to_column(tax_16S, var = "OTU")
dim(df.Occ)
head(df.Occ)
dim(tax_16S)
df.Occ.tax <- merge(df.Occ, tax_16S, by.x =c("OTU"), by.y = c("OTU"))
head(df.Occ.tax)
dim(df.Occ.tax) ### all OTU with occupancy and taxonomy!!!!!!!!!!!
# 2. Fungi
tax = read.csv("ITS_TAX.csv", sep=",", row.names=1)
tax = as.matrix(tax)
dim(tax)
tax
df.tax=as.data.frame(tax)
taxon=rownames_to_column(df.tax, var = "OTU")
head(taxon)
dim(taxon)
df.OccITS.tax <- merge(df.OccITS, taxon, by.x =c("OTU"), by.y = c("OTU"))
names(df.OccITS.tax)[names(df.OccITS.tax) == "Kingdom"] <- "Domain"
head(df.OccITS.tax)
dim(df.OccITS.tax)

# Cumulative Relative Abundance per OTU in all 45 samples
# 1. Bacteria
otu_rel <- decostand(otu, method="total", MARGIN=2)
com_abund <- rowSums(otu_rel)
df.com_abund <- as.data.frame(com_abund)
head(df.com_abund)
df.com_abund$RelAbund=df.com_abund$com_abund/45
sum(df.com_abund$com_abund)
sum(df.com_abund$RelAbund)
df.com_abund$PercentRelAbund=df.com_abund$RelAbund*100
sum(df.com_abund$PercentRelAbund)
df.com_abund=rownames_to_column(df.com_abund, var = "OTU")
head(df.com_abund)
dim(df.com_abund) ### all OTU with CumulativeRelAbund, percent CumulativeRelAbund!!!!!!!!!!!

# 2. Fungi
otuITS_rel <- decostand(otuITS, method="total", MARGIN=2)
com_abundITS <- rowSums(otuITS_rel)
df.com_abundITS <- as.data.frame(com_abundITS)
head(df.com_abundITS)
df.com_abundITS$RelAbund=df.com_abundITS$com_abundITS/45
sum(df.com_abundITS$com_abundITS)
sum(df.com_abundITS$RelAbund)
df.com_abundITS$PercentRelAbund=df.com_abundITS$RelAbund*100
sum(df.com_abundITS$PercentRelAbund)
df.com_abundITS=rownames_to_column(df.com_abundITS, var = "OTU")
head(df.com_abundITS)
dim(df.com_abundITS) ### all OTU with CumulativeRelAbund, percent CumulativeRelAbund!!!!!!!!!!!

# Relative Abundance (or mean Relative Abundance) each OTU in each sample
# 1. Bacteria
otu <- otu[rowSums(otu)>0,]
otu.ID=rownames_to_column(otu,var = "OTU")
otu.melt=melt(otu.ID,variable.name = "Sample")
sum(otu.melt$value) #1247220
otu.melt$relabund = otu.melt$value/1247220
head(otu.melt)
sum(otu.melt$relabund)
otu.melt$percentrelabund=otu.melt$relabund*100
sum(otu.melt$percentrelabund)
head(otu.melt)
dim(otu.melt) ### all OTU with relabund, percent relabund for each sample!!!!!!!!!!!

# 2. Fungi
otuITS <- otuITS[rowSums(otuITS)>0,]
otuITS.ID=rownames_to_column(otuITS,var = "OTU")
otuITS.melt=melt(otuITS.ID,variable.name = "Sample")
sum(otuITS.melt$value) #2530800
otuITS.melt$relabund = otuITS.melt$value/2530800
head(otuITS.melt)
sum(otuITS.melt$relabund)
otuITS.melt$percentrelabund=otuITS.melt$relabund*100
sum(otuITS.melt$percentrelabund)
head(otuITS.melt)
dim(otuITS.melt) ### all OTU with relabund, percent relabund for each sample!!!!!!!!!!!


# merge occupancy 1 and cumulative relative abundance 
# 1. Bacteria
Occ_RelAbund <- merge(df.Occ.tax, df.com_abund, by.x =c("OTU"), by.y = c("OTU"))
df.Occ.tax1 <- subset(df.Occ.tax, df.Occ.tax$Occ==1)
head(df.Occ.tax1)
dim(df.Occ.tax1)
Occ1_RelAbund <- merge(df.Occ.tax1, df.com_abund, by.x =c("OTU"), by.y = c("OTU"))
dim(Occ1_RelAbund)
head(Occ1_RelAbund)
sort_Occ1_RelAbund <- Occ1_RelAbund[order(Occ1_RelAbund$RelAbund, decreasing = TRUE),]
#write.table(sort_Occ1_RelAbund, file = 'sort_Occ1_RelAbund.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
sort_Occ1_RelAbund <- read.table('sort_Occ1_RelAbund.txt', sep='\t', header=TRUE)
sum(sort_Occ1_RelAbund$PercentRelAbund)

#genus count = bacteria
gc.count=sort_Occ1_RelAbund %>%
 group_by(Genus_Class) %>%
 summarise(GC_count=n())

#phylum count = bacteria
phyl_count.bac <- sort_Occ1_RelAbund %>%
 group_by(Phylum) %>%
 summarise(phyl_count=n())

#class count= bacteria
class_count.bac <- sort_Occ1_RelAbund %>%
 group_by(Class) %>%
 summarise(class_count=n())
sort_class_count <- class_count.bac[order(class_count.bac$class_count, decreasing = TRUE),]
write_csv(sort_class_count,"class_count.bac.csv")
class.count=read.csv('class_count.bac.csv', header=TRUE)

# write.table(gc.count, file = 'gc.count.bac.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
# gc.count=read.table('gc.count.bac.txt', sep='\t', header=TRUE)
# gc.count$Genus_Class <- as.character(gc.count$Genus_Class)
# gc.count$Genus_Class <- factor(gc.count$Genus_Class, levels=unique(gc.count$Genus_Class))
class.count$Class <- as.character(class.count$Class)
class.count$Class <- factor(class.count$Class, levels=unique(class.count$Class))

Bac.taxa.num=ggplot(class.count, aes(x = Class, y = class_count, fill=Phylum))+ 
 geom_bar(position = "dodge",stat = "identity")+
 ylab("Number of taxa")+
 scale_y_continuous(expand = c(0,0.5))+
 scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
 theme_bw()+
 coord_flip()+
 theme(axis.title=element_text(size=12,face="bold"),
       axis.text.x = element_text(size = 12),
       axis.text.y=element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       legend.position = "right",
       panel.border = element_blank(),
       axis.line.x = element_line(colour = "black"),
       axis.line.y = element_line(colour = "black"),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank())
       #plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))
       #labs(x= "Genus (or Class)", y="Number of taxa")

# 2. Fungi
OccITS_RelAbund <- merge(df.OccITS.tax, df.com_abundITS, by.x =c("OTU"), by.y = c("OTU"))
df.OccITS.tax1 <- subset(df.OccITS.tax, df.OccITS.tax$OccITS==1)
head(df.OccITS.tax1)
dim(df.OccITS.tax1)
Occ1ITS_RelAbund <- merge(df.OccITS.tax1, df.com_abundITS, by.x =c("OTU"), by.y = c("OTU"))
dim(Occ1ITS_RelAbund)
head(Occ1ITS_RelAbund)
sort_Occ1ITS_RelAbund <- Occ1ITS_RelAbund[order(Occ1ITS_RelAbund$RelAbund, decreasing = TRUE),]
#write.table(sort_Occ1ITS_RelAbund, file = 'sort_Occ1ITS_RelAbund.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
sort_Occ1ITS_RelAbund <- read.table('sort_Occ1ITS_RelAbund.txt', sep='\t', header=TRUE)
fg.phyl.sum=sort_Occ1ITS_RelAbund %>%
 group_by(Phylum) %>%
 summarise(Phylum_count=n())
Occ1ITS_RelAbund$Phylum <- as.character(Occ1ITS_RelAbund$Phylum)
Occ1ITS_RelAbund$Phylum[is.na(Occ1ITS_RelAbund$Phylum)] <- "Fungi"

go.count.fg=sort_Occ1ITS_RelAbund %>%
 group_by(Genus_Order) %>%
 summarise(GO_count=n())

#write.table(go.count.fg, file = 'go.count.fg.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
go.count.fg=read.table('go.count.fg.txt', sep='\t', header=TRUE)
go.count.fg$Genus_Order <- as.character(go.count.fg$Genus_Order)
go.count.fg$Genus_Order <- factor(go.count.fg$Genus_Order, levels=unique(go.count.fg$Genus_Order))

Fg.taxa.num=ggplot(go.count.fg, aes(x = Genus_Order,y = GO_count, fill=Phylum))+ 
 geom_bar(position = "dodge",stat = "identity")+
 scale_y_continuous(expand = c(0,0.05))+
 scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
 coord_flip()+
 theme_bw()+
 theme(axis.title=element_text(size=12,face="bold"),
       axis.text.x = element_text(size = 12),
       axis.text.y=element_blank(),
       axis.title.y = element_blank(),
       axis.ticks.y = element_blank(),
       legend.position = "right",
       panel.border = element_blank(),
       axis.line.x = element_line(colour = "black"),
       axis.line.y = element_line(colour = "black"),
       panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank())+
       labs(x= "Genus (or Class)", y="Number of taxa")

par(mfrow=c(1,2))
# Plot Core Bacteria and Archaea; and Fungi

# 1. Bacteria
library(forcats)
library(dplyr)

dim(sort_Occ1_RelAbund)
top100.occ1 <- sort_Occ1_RelAbund[1:100,]

CoreBac <- ggplot(sort_Occ1_RelAbund,aes(x=fct_reorder(Class, RelAbund, .desc=T), y=PercentRelAbund, fill=Phylum))+
 geom_boxplot()+
 coord_flip()+
 scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
 labs(title = "A. Bacteria/archaea", y= "Relative Abundance (%)", x="Class")+
 theme_bw()+
 theme(plot.title = element_text(size=16, face="bold"),
       axis.text=element_text(size=12), 
       axis.title=element_text(size=12,face="bold"),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       #legend.position = "right",
       legend.position = "none",
       panel.background = element_blank(),
       panel.grid = element_blank(),
       panel.border = element_blank(),
       axis.line.x = element_line(colour = "black"),
       axis.line.y = element_line(colour = "black"),
       plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))

# 2. Fungi
CoreFg <- ggplot(sort_Occ1ITS_RelAbund,aes(x=fct_reorder(Genus_Order, RelAbund, .desc=T), y=RelAbund, fill=Phylum))+
 geom_boxplot()+
 coord_flip()+
 scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
 labs(title="B. Fungi", y= "Relative Abundance (%)", x="Genus (or Order)")+
 theme_bw()+
 theme(plot.title = element_text(size = 16, face="bold"),
       axis.text=element_text(size=12), 
       axis.title=element_text(size=12,face="bold"),
       legend.text=element_text(size = 8),
       legend.title = element_text(size=10),
       legend.position = "none",
       panel.background = element_blank(),
       panel.border = element_blank(),
       panel.grid = element_blank(),
       axis.line.x = element_line(colour = "black"),
       axis.line.y = element_line(colour = "black"))



####### Plotting core microbiome members ##############
grid.newpage()
gA <- ggplotGrob(CoreBac)
gB <- ggplotGrob(Bac.taxa.num)
gC <- ggplotGrob(CoreFg)
gD <- ggplotGrob(Fg.taxa.num)
grid::grid.newpage()
setEPS()
postscript("Figure 4.eps", height = 7, width = 14)
grid::grid.draw(cbind(gA, gB, gC, gD))
dev.off()
graphics.off()
# Occupancy-Abundance Plot
# 1. Bacteria
color_top <- df.com_abund$RelAbund
color_top <- Occ
color_top[] <- 'black' 
Occ.1 <- Occ[Occ==1]
color_top[names(color_top) %in% names(Occ.1)] <- 'orange'
# Default plot in R
plot(log10(df.com_abund$RelAbund), Occ, col=color_top, pch=20, ylab='Occupancy', xlab='log(Mean of relative abundance)')
# Plot using ggplot2
occ_bacless1 <- subset(Occ_RelAbund, Occ != 1)
Occ.RelAbun.Bac <- ggplot()+
  geom_point(aes(x=log10(RelAbund), y=Occ), data=occ_bacless1, size=2.5, alpha=0.5, pch=21, colour="darkgrey")+
  labs(title="A. Bacteria/archaea",y= "Occupancy", x="Log10(mean of relative abundance)")+
  theme_bw()+
  theme(plot.title = element_text(size = 16, face="bold"),
       axis.text=element_text(size=9), 
       axis.title=element_text(size=15,face="bold"),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       legend.position = "bottom",
       legend.text=element_text(size = 8, face="bold"),
       legend.title = element_blank())+
       geom_point(data = Occ1_RelAbund, aes(x=log10(RelAbund), y=Occ, colour=Phylum), size=2.5, alpha=0.5)+
       scale_colour_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))
#b <- Occ.RelAbun.Bac+theme(legend.position = "none")
#ggarrange(b, f, ncol = 2, nrow = 1)
# 2. Fungi
color_top <- df.com_abundITS$RelAbund
color_top <- OccITS
color_top[] <- 'black' 
OccITS.1 <- OccITS[OccITS==1]
color_top[names(color_top) %in% names(OccITS.1)] <- 'orange'
plot(log10(df.com_abundITS$RelAbund), OccITS, col=color_top, pch=20, ylab='Occupancy', xlab='log(Mean of relative abundance)')
# Plot using ggplot2
occ_fgless1 <- subset(OccITS_RelAbund, OccITS != 1)
Occ.RelAbun.Fg <- ggplot()+
 geom_point(data=occ_fgless1,aes(x=log10(RelAbund), y=OccITS),colour="darkgrey", size=2.5, alpha=0.5, pch=21)+
 labs(y= "Occupancy", x="Log10(mean of relative abundance)")+
 ggtitle("B. Fungi")+
 theme_bw()+
 theme(plot.title = element_text(face="bold", size=16),
       axis.text=element_text(size=12), 
       axis.title.y = element_blank(),
       axis.title=element_text(size=15,face="bold"),
       panel.background = element_blank(),
       panel.grid = element_blank(),
       legend.position = "bottom",
       legend.text=element_text(size = 8, face="bold"),
       legend.title = element_text()+
       geom_point(data = Occ1ITS_RelAbund, aes(x=log10(RelAbund), y=OccITS, colour=Phylum),size=2.5, alpha=0.5)+
       scale_colour_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')))
 #f <- Occ.RelAbun.Fg+theme(legend.position = "none")
 #ggarrange(f, CoreFg, ncol = 2, nrow = 1)
 
plot <- ggarrange(Occ.RelAbun.Bac, Occ.RelAbun.Fg,
             ncol = 2, nrow = 1, align = "hv")

ggsave("SupplementaryFigureS9.tiff",
       plot, device = "tiff",
       width = 14, height = 7, 
       units= "in", dpi = 600)


######## PLOT THE MEMBER OF OCC =1 BACTERIA AND FUNGI ########
Bac.taxa.num
Fg.taxa.num
CoreBac
CoreFg

p1 <- ggarrange(CoreBac, Bac.taxa.num, nrow = 1, ncol = 2, align = "hv")

grid.newpage()
gC <- ggplotGrob(CoreFg)
gD <- ggplotGrob(Fg.taxa.num)
grid::grid.newpage()
grid::grid.draw(cbind(gC, gD))

#### NETWORK ANALYSIS #####

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
otu <- read.table('OTU_rarefied_16S.txt', sep='\t', header=T, row.names = 1)
df.otu=rownames_to_column(otu, var = "OTU_ID")
df.otu.tax <- subset(df.otu, select = c(1, 47))
dim(df.otu.tax)
write.table(df.otu.tax, file = "df.otu.tax.csv", sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
df.otu.tax.csv = read.csv("df.otu.tax.csv", sep=',', header=T)
dim(df.otu.tax.csv)
colnames(df.otu.tax.csv)[1]<-"Name"
head(df.otu.tax.csv)


allofOTU.node = read.table("allofOTU.node_attribute.txt", sep='\t', header=T)
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
head(tax.fg)

allofOTU.node.fg.tax <- merge(tax.fg,allofOTU.node.fg, by.x =c("Name"), by.y = c("Name"))
dim(allofOTU.node.fg.tax)
head(allofOTU.node.fg.tax)

# 3. Combine bac and fungi
allofOTU.node_attribute <- rbind(allofOTU.node.bac.tax, allofOTU.node.fg.tax)
head(allofOTU.node_attribute)
#write.table(allofOTU.node_attribute, file = "allofOTU.node_attribute1.txt", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

# Make Zi-Pi Plot
Module.hub <- allofOTU.node_attribute[allofOTU.node_attribute$Zi>2.5,] #12
#write.table(Module.hub, file = "Module.hub.txt", sep = '\t', col.names = TRUE, row.names = F, quote = FALSE)
Module.hub.tax <- read.table("Module.hub.edit.txt", sep='\t', header=T)
Module.col=c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff')
non.hub <- allofOTU.node_attribute[allofOTU.node_attribute$Zi<2.5,]
Peripherals <- non.hub[non.hub$Pi<0.62,] #409
Connectors <- allofOTU.node_attribute[allofOTU.node_attribute$Pi>0.62,] #5
#write.table(Connectors, file = "Connectors.txt", sep = '\t', col.names = TRUE, row.names = F, quote = FALSE)
Connectors.tax <- read.table("Connectors.edit.txt", sep='\t', header=T)

ZP.plot1 <- ggplot()+ geom_point(size=2.5, aes(y = Zi, x = Pi),
                           data = Peripherals)+
 xlab("Among-module connectivity (Pi)")+
 ylab("Within-module connectivity (Zi)")+
 geom_hline(yintercept = 2.5, linetype="dashed") +
 geom_vline(xintercept = 0.62, linetype="dashed")+
annotate("text", x = 0.08, y = 4.25, label = "Module hubs", fontface='bold', size=5)+
 annotate("text", x = 0.75, y = 4.25, label = "Network hubs",fontface='bold',size=5)+
 annotate("text", x = 0.08, y = -2, label = "Peripherals",fontface='bold',size=5)+
 annotate("text", x = 0.75, y = -2, label = "Connectors",fontface='bold',size=5)+
 xlim(0, 0.8)+
 theme_bw()+
 theme(axis.text=element_text(size=9), 
       axis.title=element_text(size=12,face="bold"),
legend.text=element_text(size = 8),
legend.title = element_text(size=10),
 panel.grid = element_blank()))
ZP.plot2 <- ZP.plot1 + 
 geom_point(size=2.5, 
            mapping=aes(y = Zi, x = Pi, colour=Module_hub),
            data = Module.hub.tax)+
 scale_color_manual(name="Module Hubs", values = Module.col)
ZP.plot3 <- ZP.plot2 + 
 geom_point(size=2.5, mapping=aes(y = Zi, x = Pi, shape=Connectors),
            data = Connectors.tax)+
 theme(legend.text = element_text(size=11),
       legend.title = element_text(size=12, face = "bold"))
 
ggsave("Figure5B.eps",
       ZP.plot3, device = "eps",
       width = 11, height = 6, 
       units= "in", dpi = 600)


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





