#SI analyses of the diet and isotopic niche of grizzly bear in BC
#performed by D. Breault in 2020-2022 for Dexter Hodder and Shelley Marshall
#Keys ----

#Choose one of the MCMC run options:

#run ==	      Chain Length	Burn-in	  Thin	# Chains
#"test"	      1,000	        500	        1	      3
#"very short"	10,000	      5,000	      5	      3
#"short"	    50,000	      25,000	    25	    3
#"normal"	    100,000	      50,000	    50	    3
#"long"	      300,000	      200,000	    100	    3
#"very long"	1,000,000	    500,000	    500	    3
#"extreme"	  3,000,000	    1,500,000	  500	    3

#set up workspace ----

.libPaths("C:/R/library")

#set working directory

setwd(("C:/R/projects/grizzly_kokanee"))

#remove any objects from previous analysis

rm(list = ls())

#Install and load packages ----

options("install.lock" = FALSE)

install.packages("dplyr")

install.packages("data.table")

install.packages('vegan') #MRPP

install.packages("siar")

install.packages("nicheROVER")

install.packages('gmodels')

install.packages("ggrepel")

install.packages("descr")

install.packages("car")

install.packages("lubridate")

install.packages('class') #KNN package

install.packages('R2jags')

install.packages("ggpubr") #for arranging multiple ggplots together

install.packages("devtools")

install.packages("mvnormtest") #MANOVA

remotes::install_github("brianstock/MixSIAR", dependencies=T)

#Load packages

library("dplyr")

library("data.table")

library("MixSIAR")

library('vegan')

library("siar")

library("nicheROVER")

library('gmodels')

library("ggrepel")

library("descr")

library("car")

library("lubridate")

library("class")

library("R2jags")

library('ggpubr')

library('devtools')

library('mvnormtest') #Multibariate test of normality

#Summarize consumer stable isotope values (Mean +/- SE) OVERALL, by AREA, SEX, and YEAR ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Sex", "Area", "d13C", "d15N")]

#Check balance of Sexes in each Area

sex_area = griz_merge %>% 
  group_by(Sex, Area) %>% 
  summarise(sample_size = n())

write.csv(sex_area, "output/summary_sex_area_220412.csv")

#summarize OVERALL

summary.all <- as.data.frame(griz_merge %>%
                                dplyr::summarize(sample_size = n(),
                                                 mean_C = mean(d13C),
                                                 se_C = sd(d13C)/sqrt(n()),
                                                 mean_N = mean(d15N),
                                                 se_N=sd(d15N)/sqrt(n())))

summary.all

write.csv(summary.all, "output/summary_griz_all_220412.csv")

#summarize by AREA

summary.area <- as.data.frame(griz_merge %>%
                                      group_by(Area) %>%
                                      dplyr::summarize(sample_size = n(),
                                                       mean_C = mean(d13C),
                                                       se_C = sd(d13C)/sqrt(n()),
                                                       mean_N = mean(d15N),
                                                       se_N=sd(d15N)/sqrt(n())))

summary.area

write.csv(summary.area, "output/summary_griz_area_220412.csv")

#summarize by SEX

summary.sex <- as.data.frame(griz_merge %>%
                                group_by(Sex) %>%
                                dplyr::summarize(sample_size = n(),
                                                 mean_C = mean(d13C),
                                                 se_C = sd(d13C)/sqrt(n()),
                                                 mean_N = mean(d15N),
                                                 se_N=sd(d15N)/sqrt(n())))

summary.sex

write.csv(summary.sex, "output/summary_griz_sex_220412.csv")

#summarize by YEAR

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

#format Date as yyyy-mm-dd

#dat.hair$Date <- as.Date(dat.hair$Date, format = "%d-%m-%Y")

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_year <- dat.hair[, c("ID", "Year")]

griz_year <- unique(griz_year)

griz_merge <- merge(griz_mean, griz_year)

griz_merge <- griz_merge[, c("Year", "d13C", "d15N")]


#Keep only grizzly sampled both years and get Mean Isotopes by Year, by ID

#dat.all <- dat.hair %>% 
#  group_by(ID) %>%
#  filter(all(c(2018, 2019, 2020) %in% Year))

#summarize by YEAR

summary.year <- as.data.frame(griz_merge %>%
                             group_by(Year)%>%
                               dplyr::summarize(sample_size = n(),
                                                mean_C = mean(d13C),
                                                se_C = sd(d13C)/sqrt(n()),
                                                mean_N = mean(d15N),
                                                se_N=sd(d15N)/sqrt(n())))

write.csv(summary.year, "output/summary_griz_year_220412.csv")

#Plot GRIZZLY HAIR ellipse OVERALL with food ITEMS UNcorrected ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Sex", "Area", "d13C", "d15N")]

write.csv(griz_merge, "data/griz_merge.csv")

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(griz_merge$d13C, griz_merge$d15N)

SE1

#create new data frame with names and ellipse outputs
x <- c(SE1$xSEAc)
y <- c(SE1$ySEAc)

df_SE <- data.frame(x, y)

plot(df_SE$x, df_SE$y)

#Combine bear raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_merge$type <- 1

griz_merge <- griz_merge[, c("d13C", "d15N", "type")]

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("d13C","d15N", "type")

#merge both dataframes

griz_both <- rbind(griz_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#shorten highbush cranberry to cranberry

dat.prey <- dat.prey %>% 
  mutate(Item = ifelse(as.character(Item) == "highbush cranberry", "cranberry", as.character(Item)))

#Combine the grizzly ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(Item) %>%
                             dplyr::summarize(mean_C = mean(d13C, na.rm = TRUE),
                                              se_C = sd(d13C, na.rm = TRUE)/sqrt(n()),
                                              mean_N = mean(d15N, na.rm = TRUE),
                                              se_N=sd(d15N, na.rm = TRUE)/sqrt(n()),
                                              sample_size = n()))

grizzly_w_prey <- ggplot(NULL) +
  
  #Add raw grizzly points
  
  geom_point(data = filter(griz_both, type == 1), aes(x = d13C, y = d15N), 
             size = 2.5, alpha = 0.8) +
  
  #Add ellipses  
  
  geom_path(data = filter(griz_both, type == 2), aes(x = d13C, y = d15N), linetype = 1, size = 1) +
  
  #Add prey means and error bars  
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = Item, fill = Item, stroke = 1), size = 2.5) +
  
  geom_text_repel(data = prey_mean, size = 4, family = "Times New Roman", point.padding = 1, 
                  box.padding = 2, segment.color = "grey", direction = "both",
                  nudge_x = ifelse(prey_mean$Item == "caribou", 2, 
                                   ifelse(prey_mean$Item == "fireweed", -2, 0)),
                  #nudge_y = ifelse(prey_mean$Item == "caribou", 1, 0),
                  max.overlaps = 20,
                  aes(x = mean_C, y = mean_N, label = Item)) +
  
  #Assign aesthetics and labels in legend
  
  scale_shape_manual(values = c(0,1,2,3,4,5,6,
                                7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) +
  
  scale_x_continuous(breaks = seq(-36, -12, 1)) +
  scale_y_continuous(breaks = seq(-6, 17, 1)) +
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Grizzly", subtitle = "and Food Items (mean ± 95% CI)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position="none",
        legend.text = element_blank())

tiff(file = "figures/grizz_fooditems_uncorrected_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipse OVERALL with food ITEMS corrected ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Sex", "Area", "d13C", "d15N")]

write.csv(griz_merge, "data/griz_merge.csv")

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(griz_merge$d13C, griz_merge$d15N)

SE1

#create new data frame with names and ellipse outputs
x <- c(SE1$xSEAc)
y <- c(SE1$ySEAc)

df_SE <- data.frame(x, y)

plot(df_SE$x, df_SE$y)

#Combine marten raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_merge$type <- 1

griz_merge <- griz_merge[, c("d13C", "d15N", "type")]

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("d13C","d15N", "type")

#merge both dataframes

griz_both <- rbind(griz_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#shorten highbush cranberry to cranberry

dat.prey <- dat.prey %>% 
  mutate(Item = ifelse(as.character(Item) == "highbush cranberry", "cranberry", as.character(Item)))

#correct prey by TEF values as in mixing models
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +4.9 for meat (ants and mammals)

dat.prey$d13C <- dat.prey$d13C + 3.7

dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
                                        ifelse(dat.prey$Source == "meat", 4.9, 
                                               ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#Combine the grizzly ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(Item) %>%
                             dplyr::summarize(mean_C = mean(d13C, na.rm = TRUE),
                                              se_C = sd(d13C, na.rm = TRUE)/sqrt(n()),
                                              mean_N = mean(d15N, na.rm = TRUE),
                                              se_N=sd(d15N, na.rm = TRUE)/sqrt(n()),
                                              sample_size = n()))

grizzly_w_prey <- ggplot(NULL) +
  
  #Add raw grizzly points
  
  geom_point(data = filter(griz_both, type == 1), aes(x = d13C, y = d15N), 
             size = 2.5, alpha = 0.8) +
  
  #Add ellipses  
  
  geom_path(data = filter(griz_both, type == 2), aes(x = d13C, y = d15N), linetype = 1, size = 1) +
  
  #Add prey means and error bars  
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = Item, fill = Item, stroke = 1), size = 2.5) +
  
  geom_text_repel(data = prey_mean, size = 4, family = "Times New Roman", point.padding = 1, 
                  box.padding = 2, segment.color = "grey", direction = "both",
                  nudge_x = ifelse(prey_mean$Item == "caribou", 2, 
                                   ifelse(prey_mean$Item == "fireweed", -2,
                                   ifelse(prey_mean$Item == "raspberry", 0.3, 
                                   ifelse(prey_mean$Item == "clover", -0.3, 0)))),
                  #nudge_y = ifelse(prey_mean$Item == "caribou", 1, 0),
                  max.overlaps = 20,
                  aes(x = mean_C, y = mean_N, label = Item)) +
  
  #Assign aesthetics and labels in legend
  
  scale_shape_manual(values = c(0,1,2,3,4,5,6,
                                7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) +
  
  scale_x_continuous(breaks = seq(-30, -12, 1)) +
  scale_y_continuous(breaks = seq(0, 17, 1)) +
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Grizzly", subtitle = "and Food Items (mean ± 95% CI)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position="none",
        legend.text = element_blank())

tiff(file = "figures/grizz_fooditems_corrected_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#K Nearest Neighbor Randomization Test to assign FOOD ITEMS to 4 FOOD GROUPS ----

#from: https://www.datacamp.com/community/tutorials/machine-learning-in-r
#acccessed on April 6, 2020

#Add prey data

rm(list = ls())

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#Drop columns with concentrations of carbon and nitrogen

drops <- c("Concd15N","Concd13C")

dat.prey <- dat.prey[ , !(names(dat.prey) %in% drops)]

#training and test sets to assess model performance later; split data into 2 sets:
#training set (2/3, 0.67) and test set (1/3, 0.33)

#each prey group must have equal chance of being assigned

#set a seed for random number generator

set.seed(1234)

#assign each row 1 (prob. 0.67) or 2 (prob. 0.33), with replacement

ind <- sample(2, nrow(dat.prey), replace=TRUE, prob=c(0.67, 0.33))

# Compose training set

prey.training <- dat.prey[ind==1, 3:4]

# Inspect training set

head(prey.training)

# Compose test set

prey.test <- dat.prey[ind==2, 3:4]

# Inspect test set

head(prey.test)

# Compose `prey` training labels based on "Item"

prey.trainLabels <- dat.prey[ind==1, 1]

# Inspect result

print(prey.trainLabels)

# Compose `iris` test labels

prey.testLabels <- dat.prey[ind==2, 1]

# Inspect result

print(prey.testLabels)

# Build the model
#set k = sqrt(n). where n = number of data points in training (should be odd num)
#k = sqrt(N)/2, therefore k = 7

prey_pred <- knn(train = prey.training, test = prey.test, cl = prey.trainLabels, k=5)

# Inspect `prey_pred`

prey_pred

# Put `prey.testLabels` in a data frame

preyTestLabels <- data.frame(prey.testLabels)

# Merge `prey_pred` and `prey.testLabels` 

prey.merge <- data.frame(prey_pred, prey.testLabels)

# Specify column names for `prey.merge`

names(prey.merge) <- c("Predicted Prey Group", "Observed Prey Group")

# Inspect `merge` 

prey.merge

#assess model performance (predicted vs observed)

prey_knn <- CrossTable(x = prey.testLabels, y = prey_pred, prop.chisq = FALSE)

print(prey_knn)

#Test if prey groups are significantly different from each other ----

#Add prey data

rm(list = ls())

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#Drop columns with concentrations of carbon and nitrogen

drops <- c("Concd15N","Concd13C", "Item")

dat.prey <- dat.prey[ , !(names(dat.prey) %in% drops)]

#Manova assumes assumes multivariate normality
#Homogeneity of variances across the range of predictors
#Linearity between all pairs of dependent variables

#convert to matrix with numeric codes for Food Groups:
#1 = kokanee, 2 = meat, 3 = plants, 4 = salmon

dat.prey <- dat.prey %>% 
  mutate(Source = ifelse(as.character(Source) == "kokanee", 1, as.numeric(Source)))

dat.prey <- dat.prey %>% 
  mutate(Source = ifelse(as.character(Source) == "meat", 2, as.numeric(Source)))

dat.prey <- dat.prey %>% 
  mutate(Source = ifelse(as.character(Source) == "plants", 3, as.numeric(Source)))

dat.prey <- dat.prey %>% 
  mutate(Source = ifelse(as.character(Source) == "salmon", 4, as.numeric(Source)))

kokanee <- subset(dat.prey, Source == 1)

meat <- subset(dat.prey, Source == 2)

plants <- subset(dat.prey, Source == 3)

salmon <- subset(dat.prey, Source == 4)

prey.mat <- data.matrix(kokanee)

prey.mat <- data.matrix(meat)

prey.mat <- data.matrix(plants)

prey.mat <- data.matrix(salmon)

mshapiro.test(t(prey.mat[,2:3]))

#Violations of normality

#homogeneity of variance across range of predictors

# Levene's test with one independent variable

leveneTest(d13C ~ Source, data = dat.prey)

leveneTest(d15N ~ Source, data = dat.prey)

# test linear correlation between C13 and N15

cor(dat.prey$d13C, dat.prey$d15N)

#perform linear regression with x = C13, and y = N15

lin.corr <- lm(d15N ~ d13C, data = dat.prey)

lin.corr

#test significance of relationship

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(lin.corr)

#significant linear correlation between d13C and d15N

#Determine if d13C and d15N differ between food groups

d13C <- dat.prey$d13C

d15N <- dat.prey$d15N

# MANOVA test

prey.man <- manova(cbind(d13C, d15N) ~ Source, data = dat.prey)

summary(prey.man)

# Look to see which differ

summary.aov(prey.man)

#MRPP of food groups

isotopes <- data.frame(dat.prey[,2:3])

group <- data.frame(dat.prey$Source)

food.mrpp <- with(group, mrpp(isotopes, dat.prey.Source))

food.mrpp

# Save and change plotting parameters

data(dune)
data(dune.env)
dune.mrpp <- with(dune.env, mrpp(dune, Management))
dune.mrpp

def.par <- par(no.readonly = TRUE)

layout(matrix(1:2,nr=1))

plot(food.ord <- metaMDS(isotopes), type = "text", distance = "euclidean", 
     autotransform = FALSE, display = "sites" )

with(group, ordihull(food.ord, dat.prey.Source))

with(food.mrpp, {
  fig.dist <- hist(boot.deltas, xlim = range(c(delta,boot.deltas)), 
                   main = "Test of Differences Among Groups")
  abline(v=delta); 
  text(delta, 2*mean(fig.dist$counts), adj = -0.5,
       expression(bold(delta)), cex=1.5 )  }
)
par(def.par)

## meandist
dune.md <- with(dune.env, meandist(vegdist(dune), Management))
dune.md
summary(dune.md)
plot(dune.md)
plot(dune.md, kind="histogram")

plot(dune.ord <- metaMDS(isotopes), type="text", display="sites" )

with(isotopes, ordihull(dune.ord, Management))

with(dune.mrpp, {
  fig.dist <- hist(boot.deltas, xlim=range(c(delta,boot.deltas)), 
                   main="Test of Differences Among Groups")
  abline(v=delta); 
  text(delta, 2*mean(fig.dist$counts), adj = -0.5,
       expression(bold(delta)), cex=1.5 )  }
)
par(def.par)
## meandist
dune.md <- with(dune.env, meandist(vegdist(dune), Management))
dune.md
summary(dune.md)
plot(dune.md)
plot(dune.md, kind="histogram")

#Summarize FOOD ITEMS by ITEM and by GROUP (means and SE) ----

rm(list = ls())

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#Summarize food items by ITEM with means and SE

summary.item <- as.data.frame(dat.prey %>%
                             group_by(Item) %>%
                             dplyr::summarize(sample_size = n(),
                                              mean_C = mean(d13C, na.rm = TRUE),
                                              se_C = sd(d13C, na.rm = TRUE)/sqrt(n()),
                                              mean_N = mean(d15N, na.rm = TRUE),
                                              se_N=sd(d15N, na.rm = TRUE)/sqrt(n())))

summary.item

write.csv(summary.item, "output/summary_item_12apr2021.csv")

#Summarize food items by GROUP with means and SE

summary.group <- as.data.frame(dat.prey %>%
                                group_by(Source) %>%
                                dplyr::summarize(sample_size = n(),
                                                 mean_C = mean(d13C, na.rm = TRUE),
                                                 se_C = sd(d13C, na.rm = TRUE)/sqrt(n()),
                                                 mean_N = mean(d15N, na.rm = TRUE),
                                                 se_N=sd(d15N, na.rm = TRUE)/sqrt(n())))

summary.group

write.csv(summary.group, "output/summary_group_12apr2021.csv")

#Plot GRIZZLY HAIR ellipse OVERALL with food GROUPS (ALL PREY INCLUDED) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Sex", "Area", "d13C", "d15N")]

write.csv(griz_merge, "data/griz_merge.csv")

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(griz_merge$d13C, griz_merge$d15N)

SE1

write.csv(SE1, "output/griz_SE_overall.csv")

#create new data frame with names and ellipse outputs
x <- c(SE1$xSEAc)
y <- c(SE1$ySEAc)

df_SE <- data.frame(x, y)

plot(df_SE$x, df_SE$y)

#Combine marten raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_merge$type <- 1

griz_merge <- griz_merge[, c("d13C", "d15N", "type")]

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("d13C","d15N", "type")

#merge both dataframes

griz_both <- rbind(griz_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#correct prey by TEF values as in mixing models
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +4.9 for meat (ants and mammals)

dat.prey$d13C <- dat.prey$d13C + 3.7

dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
                                        ifelse(dat.prey$Source == "meat", 4.9, 
                                               ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#Combine the grizzly ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(Source) %>%
                             dplyr::summarize(mean_C = mean(d13C, na.rm = TRUE),
                                              se_C = sd(d13C, na.rm = TRUE)/sqrt(n()),
                                              mean_N = mean(d15N, na.rm = TRUE),
                                              se_N=sd(d15N, na.rm = TRUE)/sqrt(n()),
                                              sample_size = n()))

grizzly_w_prey <- ggplot(NULL) +
  
  #Add raw grizzly points
  
  geom_point(data = filter(griz_both, type == 1), aes(x = d13C, y = d15N), 
             size = 2.5, alpha = 0.8) +
  
  #Add ellipses  
  
  geom_path(data = filter(griz_both, type == 2), aes(x = d13C, y = d15N), linetype = 1, size = 1) +
  
  #Add prey means and error bars  
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = Source, fill = Source, stroke = 1), size = 2.5) +
  
  geom_text_repel(data = prey_mean, size = 4, family = "Times New Roman", point.padding = 3, 
                  box.padding = 2, segment.color = NA,
                  aes(x = mean_C, y = mean_N, label = Source)) +
  
  #Assign aesthetics and labels in legend
  
  scale_shape_manual(values = c(0,1,2,4,5,6,
                                7,13,15,16,17,18,19,20,21,22)) +
  
  scale_x_continuous(breaks = seq(-30, -12, 1)) +
  scale_y_continuous(breaks = seq(0, 17, 1)) +
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Grizzly", subtitle = "and Food Groups (mean ± 95% CI)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position="none",
        legend.text = element_blank())

tiff(file = "figures/grizz_groups_corrected_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipses by AREA with food GROUPS (mean +/- SE) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Area", "d13C", "d15N")]

#test for significant differences between areas (MANOVA)

#Manova assumes assumes multivariate normality
#Homogeneity of variances across the range of predictors
#Linearity between all pairs of dependent variables

#convert to matrix with numeric codes for Groups:
#1 = North, 2 = South

griz_num <- griz_merge %>% 
  mutate(Area = ifelse(as.character(Area) == "North", 1, 2))

north <- subset(griz_num, Area == 1)

south <- subset(griz_num, Area == 2)

area.mat <- data.matrix(north)

area.mat <- data.matrix(south)

mshapiro.test(t(area.mat[,2:3]))

#no violations of normality

#homogeneity of variance across range of predictors

# Levene's test with one independent variable

leveneTest(d13C ~ Area, data = griz_merge)

leveneTest(d15N ~ Area, data = griz_merge)

# test linear correlation between C13 and N15

cor(griz_merge$d13C, griz_merge$d15N)

#perform linear regression with x = C13, and y = N15

lin.corr <- lm(d15N ~ d13C, data = griz_merge)

lin.corr

#test significance of relationship

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(lin.corr)

#significant linear correlation between d13C and d15N

#Determine if d13C and d15N differ between groups

d13C <- griz_merge$d13C

d15N <- griz_merge$d15N

# MANOVA test

area.man <- manova(cbind(d13C, d15N) ~ Area, data = griz_merge)

summary(area.man)

# Look to see which differ

summary.aov(area.man)

#MRPP

isotopes <- data.frame(griz_merge[,2:3])

area <- data.frame(griz_merge$Area)

area.mrpp <- with(area, mrpp(isotopes, griz_merge.Area, permutations = 10000, distance = "euclidean",
                  weight.type = 1, strata = NULL, parallel = 1))

area.mrpp

# Save and change plotting parameters

def.par <- par(no.readonly = TRUE)

layout(matrix(1:2, nr =1 ))

plot(area.ord <- metaMDS(isotopes), type = "text", display = "areas" )

with(area, ordihull(area.ord, griz_merge.Area))

with(area.mrpp, {
  fig.dist <- hist(boot.deltas, xlim = range(c(delta, boot.deltas)), 
                   main = "Test of Differences Among Groups")
  abline(v = delta); 
  text(delta, 2*mean(fig.dist$counts), adj = -0.5,
       expression(bold(delta)), cex = 1.5 )  }
)

par(def.par)

## meandist

area.md <- with(area, meandist(vegdist(isotopes), griz_merge.Area))

area.md

summary(area.md)

plot(area.md)

plot(area.md, kind = "histogram")

#subset consumer data by area manually

North <- subset(griz_merge, Area == 'North')

South <- subset(griz_merge, Area == 'South')

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(North$d13C, North$d15N)
SE2 <- standard.ellipse(South$d13C, South$d15N)

SE1
SE2

write.csv(SE1, "output/griz_SE_North.csv")

write.csv(SE2, "output/griz_SE_South.csv")

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

North_ <- rep("North", length(SE1$xSEAc))
South_ <- rep("South", length(SE2$xSEAc))

#create new data frame with names and ellipse outputs

Area <- c(North_, South_)
x <- c(SE1$xSEAc,SE2$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc)

df_SE <- data.frame(x, y, Area)

plot(df_SE$x, df_SE$y)

#Combine marten raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_merge$type <- 1

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("d13C","d15N", "Area", "type")

df_SE <- df_SE[c("Area", "d13C","d15N", "type")]

#merge both dataframes

griz_both <- rbind(griz_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#correct prey by TEF values as in mixing models
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +4.9 for meat (ants and mammals)

dat.prey$d13C <- dat.prey$d13C + 3.7

dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
                                    ifelse(dat.prey$Source == "meat", 4.9, 
                                           ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#Combine the grizzly ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(Source) %>%
                             dplyr::summarize(mean_C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              mean_N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))
  
  grizzly_w_prey <- ggplot(NULL) +
  
#Add raw grizzly points
    
    geom_point(data = filter(griz_both, type == 1), aes(x = d13C, y = d15N, fill = Area, shape = Area), 
               size = 2.5, alpha = 0.8) +
    
    #Add ellipses  
    
    geom_path(data = filter(griz_both, type == 2), aes(x = d13C, y = d15N, 
                                                         color = Area), linetype = 1, size = 1) +
    
    #Add prey means and error bars  
    
    geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                        ymax = mean_N + 1.96*se_N), width = .2) +
    geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                         xmax = mean_C + 1.96*se_C), height =.2) +
    
    geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                     shape = Source, fill = Source, stroke = 1), size = 2.5) +
    
    #Assign aesthetics and labels in legend
    
    scale_shape_manual(values = c(22, 24, 21, 22, 23, 24),
                       breaks=c("North", "South",
                                "plants", "meat", "salmon", "kokanee"),
                       labels=c("North", "South",
                                "Plants", "Meat", "Salmon", "Kokanee")) +
    
    scale_fill_manual(values = c("grey25", "grey65",
                                 "white", "white", "white", "white"),
                      breaks=c("North", "South",
                               "plants", "meat", "salmon", "kokanee"),
                      labels=c("North", "South",
                               "Plants", "Meat", "Salmon", "Kokanee")) +
    
    scale_color_manual(values = c("grey25", "grey65"),
                       breaks=c("North", "South"),
                       guide = FALSE) + 

  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Grizzly by Area", subtitle = "and Food Groups (mean ± 95% CI)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.88, .175), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
    
    scale_x_continuous(breaks = seq(-30, -12, 1)) +
    scale_y_continuous(breaks = seq(0, 17, 1))

tiff(file = "figures/grizz_area_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipses by SEX with food GROUPS (mean +/- SE) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_Sex <- dat.hair[, c("ID", "Sex")]

griz_Sex <- unique(griz_Sex)

griz_merge <- merge(griz_mean, griz_Sex)

griz_merge <- griz_merge[, c("Sex", "d13C", "d15N")]

#test for significant differences between sexes (MANOVA)

#Manova assumes assumes multivariate normality
#Homogeneity of variances across the range of predictors
#Linearity between all pairs of dependent variables

#convert to matrix with numeric codes for Groups:
#1 = Female, 2 = Male

griz_num <- griz_merge %>% 
  mutate(Sex = ifelse(as.character(Sex) == "F", 1, 2))

female <- subset(griz_num, Sex == 1)

male <- subset(griz_num, Sex == 2)

sex.mat <- data.matrix(female)

sex.mat <- data.matrix(male)

mshapiro.test(t(sex.mat[,2:3]))

#no violations of normality

#homogeneity of variance across range of predictors

# Levene's test with one independent variable

leveneTest(d13C ~ Sex, data = griz_merge)

leveneTest(d15N ~ Sex, data = griz_merge)

# test linear correlation between C13 and N15

cor(griz_merge$d13C, griz_merge$d15N)

#perform linear regression with x = C13, and y = N15

lin.corr <- lm(d15N ~ d13C, data = griz_merge)

lin.corr

#test significance of relationship

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(lin.corr)

#significant linear correlation between d13C and d15N

#Determine if d13C and d15N differ between groups

d13C <- griz_merge$d13C

d15N <- griz_merge$d15N

# MANOVA test

sex.man <- manova(cbind(d13C, d15N) ~ Sex, data = griz_merge)

summary(sex.man)

# Look to see which differ

summary.aov(sex.man)

#MRPP

isotopes <- data.frame(griz_merge[,2:3])

sex <- data.frame(griz_merge$Sex)

sex.mrpp <- with(sex, mrpp(isotopes, griz_merge.Sex, permutations = 10000, distance = "euclidean",
                             weight.type = 1, strata = NULL, parallel = 1))

sex.mrpp

#subset consumer data by Sex manually

F <- subset(griz_merge, Sex == 'F')

M <- subset(griz_merge, Sex == 'M')

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(F$d13C, F$d15N)
SE2 <- standard.ellipse(M$d13C, M$d15N)

SE1
SE2

write.csv(SE1, "output/griz_SE_Female.csv")

write.csv(SE2, "output/griz_SE_Male.csv")

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

F_ <- rep("F", length(SE1$xSEAc))
M_ <- rep("M", length(SE2$xSEAc))

#create new data frame with names and ellipse outputs

Sex <- c(F_, M_)
x <- c(SE1$xSEAc,SE2$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc)

df_SE <- data.frame(x, y, Sex)

plot(df_SE$x, df_SE$y)

#Combine marten raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_merge$type <- 1

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("d13C","d15N", "Sex", "type")

df_SE <- df_SE[c("Sex", "d13C","d15N", "type")]

#merge both dataframes

griz_both <- rbind(griz_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#correct prey by TEF values as in mixing models
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +4.9 for meat (ants and mammals)

dat.prey$d13C <- dat.prey$d13C + 3.7

dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
                                        ifelse(dat.prey$Source == "meat", 4.9, 
                                               ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#Combine the grizzly ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(Source) %>%
                             dplyr::summarize(mean_C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              mean_N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

grizzly_w_prey <- ggplot(NULL) +
  
  #Add raw grizzly points
  
  geom_point(data = filter(griz_both, type == 1), aes(x = d13C, y = d15N, fill = Sex, shape = Sex), 
             size = 2.5, alpha = 0.8) +
  
  #Add ellipses  
  
  geom_path(data = filter(griz_both, type == 2), aes(x = d13C, y = d15N, 
                                                     color = Sex), linetype = 1, size = 1) +
  
  #Add prey means and error bars  
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = Source, fill = Source, stroke = 1), size = 2.5) +
  
  #Assign aesthetics and labels in legend
  
  scale_shape_manual(values = c(22, 24, 21, 22, 23, 24),
                     breaks=c("F", "M",
                              "plants", "meat", "salmon", "kokanee"),
                     labels=c("Female", "Male",
                              "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_fill_manual(values = c("grey25", "grey65",
                               "white", "white", "white", "white"),
                    breaks=c("F", "M",
                             "plants", "meat", "salmon", "kokanee"),
                    labels=c("Female", "Male",
                             "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_color_manual(values = c("grey25", "grey65"),
                     breaks=c("F", "M"),
                     guide = FALSE) + 
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Grizzly by Sex", subtitle = "and Food Groups (mean ± 95% CI)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.88, .175), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  scale_x_continuous(breaks = seq(-30, -12, 1)) +
  scale_y_continuous(breaks = seq(0, 17, 1))

tiff(file = "figures/grizz_sex_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipses by YEAR with food GROUPS (mean +/- SE) ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

#format Date as yyyy-mm-dd

#dat.hair$Date <- as.Date(dat.hair$Date, format = "%d-%m-%Y")

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_year <- dat.hair[, c("ID", "Year")]

griz_year <- unique(griz_year)

griz_merge <- merge(griz_mean, griz_year)

griz_merge <- griz_merge[, c("Year", "d13C", "d15N")]

#test for significant differences between years (MANOVA)

#Manova assumes assumes multivariate normality
#Homogeneity of variances across the range of predictors
#Linearity between all pairs of dependent variables

#convert to matrix with numeric codes for Groups:
#1 = 2018, 2 = 2019, 3 = 2020

griz_num <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2018", "year1", 
                ifelse(as.character(Year) == "2019", "year2", "year3")))

year1 <- subset(griz_num, Year == 1)

year2 <- subset(griz_num, Year == 2)

year3 <- subset(griz_num, Year == 3)

year.mat <- data.matrix(year1)

year.mat <- data.matrix(year2)

year.mat <- data.matrix(year3)

mshapiro.test(t(year.mat[,2:3]))

#2019 violated assumptions of normality

#MRPP

isotopes <- data.frame(griz_merge[,2:3])

year <- data.frame(griz_merge$Year)

year.mrpp <- with(year, mrpp(isotopes, griz_merge.Year, permutations = 10000, distance = "euclidean",
                           weight.type = 1, strata = NULL, parallel = 1))

year.mrpp

#homogeneity of variance across range of predictors

# Levene's test with one independent variable

leveneTest(mean_C ~ Year, data = griz_merge)

leveneTest(mean_N ~ Year, data = griz_merge)

# test linear correlation between C13 and N15

cor(griz_merge$mean_C, griz_merge$mean_N)

#perform linear regression with x = C13, and y = N15

lin.corr <- lm(mean_N ~ mean_C, data = griz_merge)

lin.corr

#test significance of relationship

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(lin.corr)

#significant linear correlation between d13C and d15N

#Determine if d13C and d15N differ between groups

d13C <- griz_merge$mean_C

d15N <- griz_merge$mean_N

# MANOVA test

year.man <- manova(cbind(d13C, d15N) ~ Year, data = griz_merge)

summary(year.man)

# Look to see which differ

summary.aov(year.man)

#subset consumer data by Year manually

year1 <- subset(griz_num, Year == "year1")

year2 <- subset(griz_num, Year == "year2")

year3 <- subset(griz_num, Year == "year3")

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(year1$d13C, year1$d15N)
SE2 <- standard.ellipse(year2$d13C, year2$d15N)
SE3 <- standard.ellipse(year3$d13C, year3$d15N)

SE1
SE2
SE3

write.csv(SE1, "output/griz_SE_2018.csv")

write.csv(SE2, "output/griz_SE_2019.csv")

write.csv(SE3, "output/griz_SE_2020.csv")

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

year1_ <- rep("year1", length(SE1$xSEAc))
year2_ <- rep("year2", length(SE2$xSEAc))
year3_ <- rep("year3", length(SE3$xSEAc))

#create new data frame with names and ellipse outputs

Year <- c(year1_, year2_, year3_)
x <- c(SE1$xSEAc,SE2$xSEAc,SE3$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc)

df_SE <- data.frame(Year, x, y)

plot(df_SE$x, df_SE$y)

#Combine grizzly raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_num$type <- 1

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("Year", "d13C","d15N", "type")

#df_SE <- df_SE[c("Year", "d13C","d15N", "type")]

#merge both dataframes

griz_both <- rbind(griz_num, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE, sep = ",")

#remove grizzly and balance out meat group by subsampling

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

#correct prey by TEF values as in mixing models
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +4.9 for meat (ants and mammals)

dat.prey$d13C <- dat.prey$d13C + 3.7

dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
                                        ifelse(dat.prey$Source == "meat", 4.9, 
                                               ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#Combine the grizzly ellipses with the corrected prey group means and error bars

prey_mean <- as.data.frame(dat.prey %>%
                             group_by(Source) %>%
                             dplyr::summarize(mean_C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              mean_N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

grizzly_w_prey <- ggplot(NULL) +
  
  #Add raw grizzly points
  
  geom_point(data = filter(griz_both, type == 1), aes(x = d13C, y = d15N, fill = Year, shape = Year), 
             size = 2.5, alpha = 0.8) +
  
  #Add ellipses  
  
  geom_path(data = filter(griz_both, type == 2), aes(x = d13C, y = d15N, 
                                                     color = Year), linetype = 1, size = 1) +
  
  #Add prey means and error bars  
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = Source, fill = Source, stroke = 1), size = 2.5) +
  
  #Assign aesthetics and labels in legend
  
  scale_shape_manual(values = c(22, 21, 25, 21, 22, 23, 24),
                     breaks=c("year1", "year2", "year3",
                              "plants", "meat", "salmon", "kokanee"),
                     labels=c("2018", "2019", "2020",
                              "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_fill_manual(values = c("grey25", "grey65", "grey85",
                               "white", "white", "white", "white"),
                    breaks=c("year1", "year2", "year3",
                             "plants", "meat", "salmon", "kokanee"),
                    labels=c("2018", "2019", "2020",
                             "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_color_manual(values = c("grey25", "grey65", "grey85"),
                     breaks=c("year1", "year2", "year3"),
                     guide = FALSE) + 
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Grizzly by Year", subtitle = "and Food Groups (mean ± 95% CI)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=15, face = "bold.italic"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "14"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "14"),
        legend.position = c(0.88, .175), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  scale_x_continuous(breaks = seq(-30, -12, 1)) +
  scale_y_continuous(breaks = seq(0, 17, 1))

tiff(file = "figures/grizz_year_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Calculate % niche overlap between SEX ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Sex", "d13C", "d15N")]

griz_merge <- griz_merge %>% 
                  rename(C = d13C, N = d15N)

# format data for plotting function

dat.overlap <- tapply(1:nrow(griz_merge), griz_merge$Sex, function(ii) X = griz_merge[ii,2:3])

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each group

nsamples <- 500

griz.par <- tapply(1:nrow(griz_merge), griz_merge$Sex,
                   function(ii) niw.post(nsamples = nsamples, X = griz_merge[ii,2:3]))

griz.par <- Filter(Negate(is.null), griz.par)

# format data for plotting function

griz.data <- tapply(1:nrow(griz_merge), griz_merge$Sex, function(ii) X = griz_merge[ii,2:3])

griz.data <- Filter(Negate(is.null), griz.data)

clrs <- c("black", "red", "blue") # colors for each group

niche.plot(niche.par = griz.par, niche.data = griz.data, pfrac = .1,
           iso.names = expression(delta^{13}*C, delta^{15}*N),
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

over <- overlap(griz.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))

# posterior expectations of overlap metrics

over.mean <- apply(over*100, c(1:2, 4), mean)

round(over.mean)

write.csv(over.mean, "output/overlap_sex_220412.csv")

#Calculate % niche overlap between AREA ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Area", "d13C", "d15N")]

griz_merge <- griz_merge %>% 
  rename(C = d13C, N = d15N)

# format data for plotting function

dat.overlap <- tapply(1:nrow(griz_merge), griz_merge$Area, function(ii) X = griz_merge[ii,2:3])

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each group

nsamples <- 500

griz.par <- tapply(1:nrow(griz_merge), griz_merge$Area,
                   function(ii) niw.post(nsamples = nsamples, X = griz_merge[ii,2:3]))

griz.par <- Filter(Negate(is.null), griz.par)

# format data for plotting function

griz.data <- tapply(1:nrow(griz_merge), griz_merge$Area, function(ii) X = griz_merge[ii,2:3])

griz.data <- Filter(Negate(is.null), griz.data)

clrs <- c("black", "red", "blue") # colors for each group

niche.plot(niche.par = griz.par, niche.data = griz.data, pfrac = .1,
           iso.names = expression(delta^{13}*C, delta^{15}*N),
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

over <- overlap(griz.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))

# posterior expectations of overlap metrics

over.mean <- apply(over*100, c(1:2, 4), mean)

round(over.mean)

write.csv(over.mean, "output/overlap_area_220412.csv")

#Calculate % niche overlap between YEAR ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

#format Date as yyyy-mm-dd

dat.hair$Date <- as.Date(dat.hair$Date, format = "%d-%m-%Y")

#Keep only grizzly sampled both years and get Mean Isotopes by Year, by ID

dat.all <- dat.hair %>% 
  group_by(ID) %>%
  filter(all(c(2018, 2019, 2020) %in% Year))

griz_mean <- as.data.frame(dat.all %>%
                             group_by(ID, Year)%>%
                             dplyr::summarize(mean_C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              mean_N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_merge <- griz_mean[, c("Year", "mean_C", "mean_N")]

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2018", "year1", as.character(Year)))

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2019", "year2", as.character(Year)))

griz_merge <- griz_merge %>% 
                  rename(C = d13C, N = d15N)

# format data for plotting function

dat.overlap <- tapply(1:nrow(griz_merge), griz_merge$Year, function(ii) X = griz_merge[ii,2:3])

# niche overlap plots for 95% niche region sizes
# generate parameter draws from the "default" posteriors of each group

nsamples <- 500

griz.par <- tapply(1:nrow(griz_merge), griz_merge$Year,
                   function(ii) niw.post(nsamples = nsamples, X = griz_merge[ii,2:3]))

griz.par <- Filter(Negate(is.null), griz.par)

# format data for plotting function

griz.data <- tapply(1:nrow(griz_merge), griz_merge$Year, function(ii) X = griz_merge[ii,2:3])

griz.data <- Filter(Negate(is.null), griz.data)

clrs <- c("black", "red", "blue") # colors for each group

niche.plot(niche.par = griz.par, niche.data = griz.data, pfrac = .1,
           iso.names = expression(delta^{13}*C, delta^{15}*N),
           col = clrs, xlab = expression("Isotope Ratio (\u2030)"))

# overlap calculation. use nsamples = nprob = 1e4 for more accurate results.

over <- overlap(griz.par, nreps = nsamples, nprob = 1e4, alpha = c(.95, .99))

# posterior expectations of overlap metrics

over.mean <- apply(over*100, c(1:2, 4), mean)

round(over.mean)

write.csv(over.mean, "output/overlap_year_220412.csv")

#Correlation between N and C values and Sampling Date? ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

#format Date as yyyy-mm-dd

dat.hair$Date <- as.Date(dat.hair$Date, format = "%d-%m-%Y")

#subset consumer data by Year

year1 <- subset(dat.hair, Year == '2018')

year2 <- subset(dat.hair, Year == '2019')

year3 <- subset(dat.hair, Year == '2020')

#add julian day - numeric

year1$jDate <- format(year1$Date, "%j")

year1$jDate <- as.numeric(year1$jDate)

year2$jDate <- format(year2$Date, "%j")

year2$jDate <- as.numeric(year2$jDate)

year3$jDate <- format(year3$Date, "%j")

year3$jDate <- as.numeric(year3$jDate)

#remove outlying value from 2019

year2 <- year2[-c(29), ]

#plot the data

ggplot(year1, aes(x = jDate, y = d13C)) +
  geom_point() +
  stat_smooth()

ggplot(year1, aes(x = jDate, y = d15N)) +
  geom_point() +
  stat_smooth()

ggplot(year2, aes(x = jDate, y = d13C)) +
  geom_point() +
  stat_smooth()

ggplot(year2, aes(x = jDate, y = d15N)) +
  geom_point() +
  stat_smooth()

ggplot(year3, aes(x = jDate, y = d13C)) +
  geom_point() +
  stat_smooth()

ggplot(year3, aes(x = jDate, y = d15N)) +
  geom_point() +
  stat_smooth()

# test linear correlation between C13 and N15 and Time

cor(year1$d13C, year1$jDate)

cor(year1$d15N, year1$jDate)

cor(year2$d13C, year2$jDate)

cor(year2$d15N, year2$jDate)

cor(year3$d13C, year3$jDate)

cor(year3$d15N, year3$jDate)

#perform linear regression with x = jDate and y = d13C, and y = d15N

lin.C1 <- lm(d13C ~ jDate, data = year1)

lin.C1

lin.N1 <- lm(d15N ~ jDate, data = year1)

lin.N1

lin.C2 <- lm(d13C ~ jDate, data = year2)

lin.C2

lin.N2 <- lm(d15N ~ jDate, data = year2)

lin.N2

lin.C3 <- lm(d13C ~ jDate, data = year3)

lin.C3

lin.N3 <- lm(d15N ~ jDate, data = year3)

lin.N3

#test significance of relationship

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(lin.C1)

summary(lin.N1)

summary(lin.C2)

summary(lin.N2)

summary(lin.C3)

summary(lin.N3)

#calculate the 95% CI around the beta coefficients

confint(lin.C1)

confint(lin.N1)

confint(lin.C2)

confint(lin.N2)

confint(lin.C3)

confint(lin.N3)

c_time1 <- ggplot(year1, aes(jDate, d13C)) +
  geom_point(data = year1, aes(x = jDate, y = d13C)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"13C"), title = "2018") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12")) +
  scale_x_continuous(breaks = seq(150, 350, 10), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-30, -20, 0.5), expand = c(0.01, 0.01))

n_time1 <- ggplot(year1, aes(jDate, d15N)) +
  geom_point(data = year1, aes(x = jDate, y = d15N)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"15N")) +
  theme(axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.9, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12")) +
  scale_x_continuous(breaks = seq(150, 350, 10), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 8, 0.5), expand = c(0.01, 0.01))

c_time2 <- ggplot(year2, aes(jDate, d13C)) +
  geom_point(data = year2, aes(x = jDate, y = d13C)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"13C"), title = "2019") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12")) +
  scale_x_continuous(breaks = seq(150, 350, 10), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-30, -20, 0.5), expand = c(0.01, 0.01))

n_time2 <- ggplot(year2, aes(jDate, d15N)) +
  geom_point(data = year2, aes(x = jDate, y = d15N)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"15N")) +
  theme(axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.9, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12")) +
  scale_x_continuous(breaks = seq(150, 350, 10), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 8, 0.5), expand = c(0.01, 0.01))

c_time3 <- ggplot(year3, aes(jDate, d13C)) +
  geom_point(data = year3, aes(x = jDate, y = d13C)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"13C"), title = "2020") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12")) +
  scale_x_continuous(breaks = seq(150, 350, 10), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(-30, -20, 0.5), expand = c(0.01, 0.01))

n_time3 <- ggplot(year3, aes(jDate, d15N)) +
  geom_point(data = year3, aes(x = jDate, y = d15N)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"15N")) +
  theme(axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.9, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12")) +
  scale_x_continuous(breaks = seq(150, 350, 10), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 8, 0.5), expand = c(0.01, 0.01))

tiff("figures/sampling_dates_220412.tiff", units="cm", width=45, height=25, res=300)

ggarrange(c_time1, c_time2, c_time3, n_time1, n_time2, n_time3,
          ncol = 3, nrow = 2)

dev.off()

#MixSIAR with Grizzly and 4 Sources OVERALL (normal) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

#griz_merge <- griz_merge[-c(42), ]

griz_merge <- griz_merge[, c("Sex", "Area", "d13C", "d15N")]

write.csv(griz_merge, "data/griz_merge.csv")

# Load the mixture/consumer data

mix <- load_mix_data(filename = "data/griz_merge.csv", 
                     iso_names = c("d13C","d15N"), 
                     factors = NULL, 
                     fac_random = NULL, 
                     fac_nested = NULL, 
                     cont_effects = NULL)

#load_source_data

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

write.csv(dat.prey, "data/knn_sources2_10apr2021.csv", row.names = FALSE)

#Add prey data and plot as scatter plot

source <- load_source_data(filename = "data/knn_sources2_10apr2021.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_var19mar2021.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_plot_norm3_220412", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

#plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

norm_griz22 <- "data/norm_griz_220412.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(norm_griz22, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.norm22 <- run_model(run = "normal", mix, source, discr, norm_griz22)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/norm22_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/norm22_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/norm22_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/norm22_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/norm22_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.norm22, mix, source, output_options)

# Save an object to a file

save(jags.norm22, mix, source, output_options, file = "RData/all_grizzly_norm_220412.RData")

# Restore the objects

load(file = "RData/all_grizzly_norm_220412.RData")

diag <- output_diagnostics(jags.norm22, mix, source, output_options)

names(diag)

head(diag$geweke)

df.stats <- output_stats(jags.norm22, mix, source, output_options)

write.csv(df.stats, "output/results_all_griz_norm_220412.csv")

g.post <- output_posteriors(jags.norm22, mix, source, output_options)

all_griz <- g.post$global +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "All Grizzly (n = 55)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +

tiff(file = "figures/all_grizzly_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

all_griz

dev.off()

#MixSIAR with Grizzly and 4 Sources by SEX (long) ----

#Take average d13C and d15N values of individuals who were sampled > 1

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Sex", "Area", "d13C", "d15N")]

write.csv(griz_merge, "data/griz_merge.csv")

# Load the mixture/consumer data

mix <- load_mix_data(filename = "data/griz_merge.csv",
                     iso_names = c("d13C","d15N"), 
                     factors = c("Sex"), 
                     fac_random = c(FALSE), 
                     fac_nested = c(FALSE), 
                     cont_effects = NULL)

#Add prey data and plot as scatter plot

#load_source_data

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

write.csv(dat.prey, "data/knn_sources2_10apr2021.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources2_10apr2021.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)
#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_var19mar2021.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_sex_220412", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

sex_griz_long22 <- "data/sex_griz_long22.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(sex_griz_long22, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.sex22 <- run_model(run = "long", mix, source, discr, sex_griz_long22)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/sex_long_summary_stats22",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/sex_long_posterior_density22",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = FALSE,
                       plot_pairs_name = "output/sex_long_pairs_plot22",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "output/sex_long_xy_plot22",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/sex_long_diagnostics22",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = TRUE,
                       plot_xy_save_png = TRUE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.sex22, mix, source, output_options)

# Save an object to a file

save(jags.sex22, mix, source, output_options, file = "RData/sex_grizzly_long22.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/sex_grizzly_long22.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.sex22, mix, source, output_options)

write.csv(df.stats, "output/results_sex_griz_long_220412.csv")

source$source_names

#Do males eat more meat than females? 

# Shapiro-Wilk normality test for meat

shapiro.test(jags.sex22$BUGSoutput$sims.list$p.fac1[,2,2]) # p << 0.005

shapiro.test(jags.sex22$BUGSoutput$sims.list$p.fac1[,1,2]) # p << 0.005

male.meat <- as.vector(jags.sex22$BUGSoutput$sims.list$p.fac1[,2,2])

female.meat <- as.vector(jags.sex22$BUGSoutput$sims.list$p.fac1[,1,2])

res <- wilcox.test(male.meat, female.meat,
                   exact = FALSE, alternative = "greater")

res

#Do females eat more plants than males?

# Shapiro-Wilk normality test for plants

shapiro.test(jags.sex22$BUGSoutput$sims.list$p.fac1[,2,3]) # p << 0.005

shapiro.test(jags.sex22$BUGSoutput$sims.list$p.fac1[,1,3]) # p << 0.005

male.plants <- as.vector(jags.sex22$BUGSoutput$sims.list$p.fac1[,2,3])

female.plants <- as.vector(jags.sex22$BUGSoutput$sims.list$p.fac1[,1,3])

res <- wilcox.test(female.plants, male.plants,
                   exact = FALSE, alternative = "greater")

res

#plot histo of difference between proportion of meat between Females and Males

p.meat.F.M <- jags.sex22$BUGSoutput$sims.list$p.fac1[,2,2] - jags.sex22$BUGSoutput$sims.list$p.fac1[,1,2]

p.meat.F.M <- as.data.frame(p.meat.F.M)

histo.meat <- ggplot(p.meat.F.M, aes(p.meat.F.M)) +
 
  geom_histogram(binwidth = 0.01) +
  
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 350, 50), limits = c(0, 350), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Meat proportions, Male - Female", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = NULL)

tiff(file = "figures/sex_meat_histo_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

histo.meat

dev.off()
  
#plot histo of difference between proportion of plants between Females and Males

p.plants.F.M <- jags.sex22$BUGSoutput$sims.list$p.fac1[,1,3] - jags.sex22$BUGSoutput$sims.list$p.fac1[,2,3]

p.plants.F.M <- as.data.frame(p.plants.F.M)

histo.plants <- ggplot(p.plants.F.M, aes(p.plants.F.M)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 350, 50), limits = c(0, 350), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Plant proportions, Female - Male", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = NULL)

tiff(file = "figures/sex_plants_histo_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

histo.plants

dev.off()

#Plot posterior densities of females and males

g.post <- output_posteriors(jags.sex22, mix, source, output_options)

female_griz <- g.post$fac1[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Female Grizzly (n = 32)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/female_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

female_griz

dev.off()

male_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Male Grizzly (n = 23)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/male_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

male_griz

dev.off()

#MixSIAR with Grizzly and 4 Sources by AREA (long) ----

#Take average d13C and d15N values of individuals who were sampled > 1

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID) %>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_area <- dat.hair[, c("ID", "Area", "Sex")]

griz_area <- unique(griz_area)

griz_merge <- merge(griz_mean, griz_area)

griz_merge <- griz_merge[, c("Sex", "Area", "d13C", "d15N")]

write.csv(griz_merge, "data/griz_merge.csv")

# Load the mixture/consumer data

mix <- load_mix_data(filename = "data/griz_merge.csv",
                     iso_names = c("d13C","d15N"), 
                     factors = c("Area"), 
                     fac_random = c(FALSE), 
                     fac_nested = c(FALSE), 
                     cont_effects = NULL)

#load_source_data

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

write.csv(dat.prey, "data/knn_sources2_10apr2021.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources2_10apr2021.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_var19mar2021.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_area_220412", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

area_griz_long <- "data/area_griz_long.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(area_griz_long, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.area <- run_model(run = "long", mix, source, discr, area_griz_long)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/area_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "output/area_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/area_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/area_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/area_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.area, mix, source, output_options)

# Save an object to a file

save(jags.area, mix, source, output_options, file = "RData/area_grizzly_long.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/area_grizzly_long.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.area, mix, source, output_options)

write.csv(df.stats, "output/results_area_griz_long_220412.csv")

#Do southern bears eat more meat than northern bears? 

# Shapiro-Wilk normality test for meat

shapiro.test(jags.area$BUGSoutput$sims.list$p.fac1[,1,2]) # p << 0.005

shapiro.test(jags.area$BUGSoutput$sims.list$p.fac1[,2,2]) # p << 0.005

north.meat <- as.vector(jags.area$BUGSoutput$sims.list$p.fac1[,1,2])

south.meat <- as.vector(jags.area$BUGSoutput$sims.list$p.fac1[,2,2])

res <- wilcox.test(south.meat, north.meat,
                   exact = FALSE, alternative = "greater")

res

#Do northern bears eat more plants than southern bears?

# Shapiro-Wilk normality test for plants

shapiro.test(jags.area$BUGSoutput$sims.list$p.fac1[,1,3]) # p << 0.005

shapiro.test(jags.area$BUGSoutput$sims.list$p.fac1[,2,3]) # p << 0.005

north.plants <- as.vector(jags.area$BUGSoutput$sims.list$p.fac1[,1,3])

south.plants <- as.vector(jags.area$BUGSoutput$sims.list$p.fac1[,2,3])

res <- wilcox.test(north.plants, south.plants,
                   exact = FALSE, alternative = "greater")

res

#plot histo of difference between proportion of meat between North and South

p.meat.N.S <- jags.area$BUGSoutput$sims.list$p.fac1[,2,2] - jags.area$BUGSoutput$sims.list$p.fac1[,1,2]

p.meat.N.S <- as.data.frame(p.meat.N.S)

histo.meat <- ggplot(p.meat.N.S, aes(p.meat.N.S)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 350, 50), limits = c(0, 350), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Meat proportions, South - North", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = NULL)

tiff(file = "figures/area_meat_histo_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

histo.meat

dev.off()

#plot histo of difference between proportion of plants between North and South

p.plants.N.S <- jags.area$BUGSoutput$sims.list$p.fac1[,1,3] - jags.area$BUGSoutput$sims.list$p.fac1[,2,3]

p.plants.N.S <- as.data.frame(p.plants.N.S)

histo.plants <- ggplot(p.plants.N.S, aes(p.plants.N.S)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 350, 50), limits = c(0, 350), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Plant proportions, North - South", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = NULL)

tiff(file = "figures/area_plants_histo_220412.tif", 
     units = "cm", width = 25, height = 20, res = 300)

histo.plants

dev.off()

#Plot posterior densities by area

g.post <- output_posteriors(jags.area, mix, source, output_options)

north_griz <- g.post$fac1[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "North Grizzly (n = 35)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/north_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

north_griz

dev.off()

south_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "South Grizzly (n = 20)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/south_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

south_griz

dev.off()

#MixSIAR with Grizzly and 4 Sources by Year (long) ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

#format Date as yyyy-mm-dd

#dat.hair$Date <- as.Date(dat.hair$Date, format = "%d-%m-%Y")

#remove -ve value

#dat.hair <- dat.hair[-c(87), ]

#Keep only grizzly sampled both years and get Mean Isotopes by Year, by ID

#dat.all <- dat.hair %>% 
#  group_by(ID) %>%
#  filter(all(c(2018, 2019, 2020) %in% Year))

griz_mean <- as.data.frame(dat.hair %>%
                             group_by(ID, Year)%>%
                             dplyr::summarize(d13C = mean(d13C),
                                              se_C = sd(d13C)/sqrt(n()),
                                              d15N = mean(d15N),
                                              se_N=sd(d15N)/sqrt(n()),
                                              sample_size = n()))

griz_merge <- griz_mean[, c("Year", "d13C", "d15N")]

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2018", "year1", as.character(Year)))

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2019", "year2", as.character(Year)))

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2020", "year3", as.character(Year)))

write.csv(griz_merge, "data/griz_merge.csv")

# Load the mixture/consumer data

mix <- load_mix_data(filename = "data/griz_merge.csv",
                     iso_names = c("d13C","d15N"), 
                     factors = c("Year"), 
                     fac_random = c(FALSE), 
                     fac_nested = c(FALSE), 
                     cont_effects = NULL)

#load_source_data

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

write.csv(dat.prey, "data/knn_sources2_10apr2021.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources2_10apr2021.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_var19mar2021.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_year_220412", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

year_griz_long <- "data/year_griz_long.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(year_griz_long, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.year <- run_model(run = "long", mix, source, discr, year_griz_long)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/year_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "output/year_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/year_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/year_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/year_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.year, mix, source, output_options)

# Save an object to a file

save(jags.year, mix, source, output_options, file = "RData/year_grizzly_long.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/year_grizzly_long.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.year, mix, source, output_options)

write.csv(df.stats, "output/results_year_griz_long_220412.csv")

#plot posterior densities by year

g.post <- output_posteriors(jags.year, mix, source, output_options)

year1_griz <- g.post$fac1[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "2018 Grizzly (n = 31)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/2018_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

year1_griz

dev.off()

year2_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "2019 Grizzly (n = 21)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/2019_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

year2_griz

dev.off()

year3_griz <- g.post$fac1[[3]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "2020 Grizzly (n = 20)") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/2020_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

year3_griz

dev.off()

#MixSIAR with Grizzly and 4 Sources by Sampling Date (long) ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_220412.csv", header = TRUE)

#format Date as yyyy-mm-dd

dat.hair$Date <- as.Date(dat.hair$Date, format = "%d-%m-%Y")

#remove -ve value

#dat.hair <- dat.hair[-c(87), ]

#add julian day - numeric

dat.hair$jDate <- format(dat.hair$Date, "%j")

dat.hair$jDate <- as.numeric(dat.hair$jDate)

griz_jdate <- dat.hair[, c("d13C", "d15N", "jDate")]

write.csv(griz_jdate, "data/griz_jdate.csv")

# Load the mixture/consumer data

mix <- load_mix_data(filename = "data/griz_jdate.csv",
                     iso_names = c("d13C","d15N"), 
                     factors = NULL, 
                     fac_random = NULL, 
                     fac_nested = NULL, 
                     cont_effects = "jDate")

#load_source_data

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly

dat.prey <- dat.prey[-c(26, 30:33, 36:68, 70:72, 147), ]

write.csv(dat.prey, "data/knn_sources2_10apr2021.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources2_10apr2021.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_var19mar2021.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_date_220412", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

date_griz_long <- "data/date_griz_long.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(date_griz_long, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.date <- run_model(run = "long", mix, source, discr, date_griz_long)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/date_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "output/date_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/date_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/date_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/date_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.date, mix, source, output_options)

# Save an object to a file

save(jags.date, mix, source, output_options, file = "RData/date_grizzly_long.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/date_grizzly_long.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.date, mix, source, output_options)

g.post <- plot_continuous_var(
  jags.date,
  mix,
  source,
  output_options,
  alphaCI = 0.05,
  exclude_sources_below = FALSE)

cont_griz <- g.post[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(205, 296, 10), limits = c(205, 296), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Grizzly Diet over Sampling Date (n = 110)", x = "Julian Date") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=12),
        axis.title.y = element_text(family="Times New Roman", size=12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "12"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "12"),
        legend.position = c(0.85, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  tiff(file = "figures/date_grizzly_220412.tif", 
       units = "cm", width = 25, height = 20, res = 300)

cont_griz

dev.off()

#wolf vignette ----

knitr::opts_chunk$set(fig.width=7, fig.height=5) 
mixsiar.dir <- find.package("MixSIAR")
source(file.path(mixsiar.dir,"example_scripts","mixsiar_script_wolves_normal.R"))
load(url("https://github.com/brianstock/MixSIAR/blob/master/Manual/wolves_normal.RData?raw=true"))
output_options <- list(summary_save = TRUE,                 
                       summary_name = "summary_statistics", 
                       sup_post = TRUE,                    
                       plot_post_save_pdf = FALSE,           
                       plot_post_name = "posterior_density",
                       sup_pairs = TRUE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = TRUE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = FALSE,
                       return_obj = TRUE)
diag <- output_diagnostics(jags.1, mix, source, output_options)

