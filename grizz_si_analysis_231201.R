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

install.packages('gmodels') #crosstable

install.packages("ggrepel")

install.packages("descr")

install.packages("car")

install.packages("lubridate")

install.packages('class') #KNN package

install.packages('R2jags')

install.packages("ggpubr") #for arranging multiple ggplots together

install.packages("devtools")

install.packages("mvnormtest") #MANOVA; Multibariate test of normality (mshapiro.test)

install.packages("PMCMRplus") #Kruskal-Walis

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

library('mvnormtest') 

library('PMCMRplus')

#Summarize consumer stable isotope values (Mean +/- SE) OVERALL, by AREA, SEX, and YEAR ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

write.csv(sex_area, "output/summary_sex_area_231204.csv")

#summarize OVERALL

summary.all <- as.data.frame(griz_merge %>%
                                dplyr::summarize(sample_size = n(),
                                                 mean_C = mean(d13C),
                                                 se_C = sd(d13C)/sqrt(n()),
                                                 mean_N = mean(d15N),
                                                 se_N=sd(d15N)/sqrt(n())))

summary.all

write.csv(summary.all, "output/summary_griz_all_231204.csv")

#summarize by AREA

summary.area <- as.data.frame(griz_merge %>%
                                      group_by(Area) %>%
                                      dplyr::summarize(sample_size = n(),
                                                       mean_C = mean(d13C),
                                                       se_C = sd(d13C)/sqrt(n()),
                                                       mean_N = mean(d15N),
                                                       se_N=sd(d15N)/sqrt(n())))

summary.area

write.csv(summary.area, "output/summary_griz_area_231204.csv")

#summarize by SEX

summary.sex <- as.data.frame(griz_merge %>%
                                group_by(Sex) %>%
                                dplyr::summarize(sample_size = n(),
                                                 mean_C = mean(d13C),
                                                 se_C = sd(d13C)/sqrt(n()),
                                                 mean_N = mean(d15N),
                                                 se_N=sd(d15N)/sqrt(n())))

summary.sex

write.csv(summary.sex, "output/summary_griz_sex_231204.csv")

#summarize by YEAR

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

write.csv(summary.year, "output/summary_griz_year_231204.csv")

#Plot GRIZZLY HAIR ellipse OVERALL with food ITEMS UNcorrected ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE, sep = ",")

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

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

tiff(file = "figures/grizz_fooditems_uncorrected_240126.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipse OVERALL with food ITEMS corrected ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE, sep = ",")

#remove grizzly, soapberry, and balance out meat group by subsampling

#dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#shorten highbush cranberry to cranberry

dat.prey <- dat.prey %>% 
  mutate(Item = ifelse(as.character(Item) == "highbush cranberry", "cranberry", as.character(Item)))

#correct prey by TEF values as in mixing models (Y = 5.28 + 0.88X, where Y = corrected, and X = uncorrected)
#dC13 = +3.7 for all prey
#dN15 = +5.4 for plants
#       +5.6 for berry
#       +4.1 for kokanee
#       +3.9 for salmon
#       +5.0 for meat (ants and mammals)

#dat.prey$d13C <- dat.prey$d13C + 3.7

#dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.4,
#                                        ifelse(dat.prey$Source == "berry", 5.6,
#                                        ifelse(dat.prey$Source == "meat", 5.0, 
#                                               ifelse(dat.prey$Source == "kokanee", 4.1, 3.9))))

#correct prey by TEF values as in mixing models (Erlenbach 2020 - coastal, hair, Jul-Oct)
#dC13 = +5.4 for all prey
#dN15 = +3.1 for all prey

dat.prey$d13C <- dat.prey$d13C + 5.8

dat.prey$d15N <- dat.prey$d15N + 4.8

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
                  #nudge_x = ifelse(prey_mean$Item == "caribou", 2, 
                                   #ifelse(prey_mean$Item == "fireweed", -2,
                                  # ifelse(prey_mean$Item == "raspberry", 0.3, 
                                  # ifelse(prey_mean$Item == "clover", -0.3,
                                  # ifelse(prey_mean$Item == "cranberry", 1.5,0))))),
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

tiff(file = "figures/grizz_fooditems_corrected_erlenbach_240312.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#K Nearest Neighbor Randomization Test to assign FOOD ITEMS to 4 FOOD GROUPS ----

#from: https://www.datacamp.com/community/tutorials/machine-learning-in-r
#acccessed on April 6, 2020

#Add prey data

rm(list = ls())

dat.prey <- read.csv("data/knn_sources_04dec2023.csv")

#remove grizzly, soapberry, and balance out meat group by subsampling

#dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

dat.prey <- dat.prey[-c(13, 26, 31:33), ]

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
#k = sqrt(N)/2, therefore k = 5

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

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE, sep = ",")

#remove grizzly, soapberry, and balance out meat group by subsampling

#dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#Drop columns with concentrations of carbon and nitrogen

drops <- c("Concd15N","Concd13C", "Item")

dat.prey <- dat.prey[ , !(names(dat.prey) %in% drops)]

#Manova assumes assumes multivariate normality
#Homogeneity of variances across the range of predictors
#Linearity between all pairs of dependent variables

# Shapiro-Wilk normality test for C13 and N15 by prey group

with(dat.prey, shapiro.test(d15N[Source == "berry"]))

with(dat.prey, shapiro.test(d13C[Source == "berry"]))

with(dat.prey, shapiro.test(d15N[Source == "kokanee"]))

with(dat.prey, shapiro.test(d13C[Source == "kokanee"]))

with(dat.prey, shapiro.test(d15N[Source == "meat"]))

with(dat.prey, shapiro.test(d13C[Source == "meat"]))

with(dat.prey, shapiro.test(d15N[Source == "salmon"]))

with(dat.prey, shapiro.test(d13C[Source == "salmon"]))

#Violations of normality

#kruskal-wallis test of differences between > 2 groups in C13

dat.prey$Source <- as.factor(dat.prey$Source)

kruskal.test(d13C ~ Source, data = dat.prey)

#Kruskal-Wallis rank sum test

#data:  d13C by Source
#Kruskal-Wallis chi-squared = 144.77, df = 4, p-value < 2.2e-16

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(data = dat.prey, x = d13C, g = Source, dist="Tukey")

#         berry   kokanee meat    plants 
#kokanee  0.00017 -       -       -      
#meat     0.42716 6.6e-07 -       -      
#plants   0.34827 0.00752 0.00522 -      
#salmon   0.00012 3.6e-14 0.40840 9.5e-14

#NO SUBSAMPLING
#         berry   kokanee meat    plants 
#kokanee  0.00177 -       -       -      
#meat     0.12268 7.3e-12 -       -      
#plants   0.26292 0.09853 3.2e-08 -      
#salmon   5.1e-07 4.8e-14 0.00037 4.1e-14

#kruskal-wallis test of differences between > 2 groups in N15

dat.prey$Source <- as.factor(dat.prey$Source)

kruskal.test(d15N ~ Source, data = dat.prey)

#Kruskal-Wallis rank sum test

#data:  d15N by Source
#Kruskal-Wallis chi-squared = 119.88, df = 4, p-value < 2.2e-16

#NO SUBSAMPLING
#data:  d15N by Source
#Kruskal-Wallis chi-squared = 160.58, df = 4, p-value < 2.2e-16

require(PMCMR)

data(dat.prey)

attach(dat.prey)

kwAllPairsNemenyiTest(data = dat.prey, x = d15N, g = Source, dist="Tukey")

#         berry   kokanee meat  plants 
#kokanee  4.5e-09 -       -     -      
#meat     0.002   0.619   -     -      
#plants   0.236   3.4e-07 0.079 -      
#salmon   4.0e-14 0.167   0.012 3.7e-14

#NO SUBSAMPLING

#         berry   kokanee meat    plants 
#kokanee  2.7e-12 -       -       -      
#meat     5.2e-07 0.016   -       -      
#plants   0.463   2.9e-12 2.4e-06 -      
#salmon   4.3e-14 0.398   2.6e-07 5.2e-14

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

#If normal, no violations

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

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE, sep = ",")

#remove grizzly, soapberry, and balance out meat group by subsampling

#dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#Summarize food items by ITEM with means and SE

summary.item <- as.data.frame(dat.prey %>%
                             group_by(Item) %>%
                             dplyr::summarize(sample_size = n(),
                                              mean_C = mean(d13C, na.rm = TRUE),
                                              se_C = sd(d13C, na.rm = TRUE)/sqrt(n()),
                                              mean_N = mean(d15N, na.rm = TRUE),
                                              se_N=sd(d15N, na.rm = TRUE)/sqrt(n())))

summary.item

write.csv(summary.item, "output/summary_item_no_subsampling_08jan2024.csv")

#Summarize food items by GROUP with means and SE

summary.group <- as.data.frame(dat.prey %>%
                                group_by(Source) %>%
                                dplyr::summarize(sample_size = n(),
                                                 mean_C = mean(d13C, na.rm = TRUE),
                                                 se_C = sd(d13C, na.rm = TRUE)/sqrt(n()),
                                                 mean_N = mean(d15N, na.rm = TRUE),
                                                 se_N=sd(d15N, na.rm = TRUE)/sqrt(n())))

summary.group

write.csv(summary.group, "output/summary_group_no_subsampling_08jan2024.csv")

#Plot GRIZZLY HAIR ellipse OVERALL with food GROUPS (ALL PREY INCLUDED - 5 groups) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

#Combine grizzly raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_merge$type <- 1

griz_merge <- griz_merge[, c("d13C", "d15N", "type")]

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("d13C","d15N", "type")

#merge both dataframes

griz_both <- rbind(griz_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE, sep = ",")

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#correct prey by TEF values as in mixing models (Y = 5.28 + 0.88X, where Y = corrected, and X = uncorrected)
#dC13 = +3.7 for all prey
#dN15 = +5.4 for plants
#       +5.6 for berry
#       +4.1 for kokanee
#       +3.9 for salmon
#       +5.0 for meat (ants and mammals)

#dat.prey$d13C <- dat.prey$d13C + 3.7

#dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.4,
#                                        ifelse(dat.prey$Source == "berry", 5.6,
#                                               ifelse(dat.prey$Source == "meat", 5.0, 
#                                                      ifelse(dat.prey$Source == "kokanee", 4.1, 3.9))))

#correct prey by TEF values as in mixing models (Erlenbach 2020)
#dC13 = +5.4 for all prey
#dN15 = +3.1 for all prey

dat.prey$d13C <- dat.prey$d13C + 5.8

dat.prey$d15N <- dat.prey$d15N + 4.8

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
                  box.padding = 2, segment.color = "grey",
                  nudge_x = ifelse(prey_mean$Source == "berry", 2,
                                   ifelse(prey_mean$Source == "plants", -1,0)),
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

tiff(file = "figures/grizz_groups_corrected_erlenbach_240301.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipse OVERALL with food GROUPS (ALL PREY INCLUDED - 4 groups) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

#Combine grizzly raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

griz_merge$type <- 1

griz_merge <- griz_merge[, c("d13C", "d15N", "type")]

df_SE$type <- 2

#rename and re-order columns in df_SE to match griz_merge

colnames(df_SE) <- c("d13C","d15N", "type")

#merge both dataframes

griz_both <- rbind(griz_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#correct prey by TEF values as in mixing models (Y = 5.28 + 0.88X, where Y = corrected, and X = uncorrected)
#dC13 = +3.7 for all prey
#dN15 = +5.4 for plants
#       +5.6 for berry
#       +4.1 for kokanee
#       +3.9 for salmon
#       +5.0 for meat (ants and mammals)

#dat.prey$d13C <- dat.prey$d13C + 3.7

#dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.4,
#                                        ifelse(dat.prey$Source == "berry", 5.6,
#                                               ifelse(dat.prey$Source == "meat", 5.0, 
#                                                      ifelse(dat.prey$Source == "kokanee", 4.1, 3.9))))

#correct 4 prey groups by TEF values as in mixing models (Erlenbach 2020)
#dC13 = +5.8 for all prey
#dN15 = +4.8 for all prey

dat.prey$d13C <- dat.prey$d13C + 5.8

dat.prey$d15N <- dat.prey$d15N + 4.8

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
  
  geom_text_repel(data = prey_mean, size = 4, family = "Times New Roman", point.padding = 20,
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

tiff(file = "figures/grizz_4groups_corrected_erlenbach_240312.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipses by AREA with food GROUPS (mean +/- SE) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#correct prey by TEF values as in mixing models (Y = 5.28 + 0.88X, where Y = corrected, and X = uncorrected)
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +5.0 for meat (ants and mammals)

#dat.prey$d13C <- dat.prey$d13C + 3.7

#dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
#                                               ifelse(dat.prey$Source == "meat", 5.0, 
#                                                      ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#correct 4 prey groups by TEF values as in mixing models (Erlenbach 2020)
#dC13 = +5.8 for all prey
#dN15 = +4.8 for all prey

dat.prey$d13C <- dat.prey$d13C + 5.8

dat.prey$d15N <- dat.prey$d15N + 4.8

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

tiff(file = "figures/grizz_area_erlenbach_240306.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipses by SEX with food GROUPS (mean +/- SE) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#correct prey by TEF values as in mixing models (Y = 5.28 + 0.88X, where Y = corrected, and X = uncorrected)
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +5.0 for meat (ants and mammals)

#dat.prey$d13C <- dat.prey$d13C + 3.7

#dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
#                                        ifelse(dat.prey$Source == "meat", 5.0, 
#                                               ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#correct 4 prey groups by TEF values as in mixing models (Erlenbach 2020)
#dC13 = +5.8 for all prey
#dN15 = +4.8 for all prey

dat.prey$d13C <- dat.prey$d13C + 5.8

dat.prey$d15N <- dat.prey$d15N + 4.8

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

tiff(file = "figures/grizz_sex_240313.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipses by YEAR with food GROUPS (mean +/- SE) ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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
                ifelse(as.character(Year) == "2019", "year2", 
                ifelse(as.character(Year) == "2020", "year3",
                ifelse(as.character(Year) == "2021", "year4", "year5")))))

year1 <- subset(griz_num, Year == 1)

year2 <- subset(griz_num, Year == 2)

year3 <- subset(griz_num, Year == 3)

year4 <- subset(griz_num, Year == 4)

year5 <- subset(griz_num, Year == 5)

year.mat <- data.matrix(year1)

year.mat <- data.matrix(year2)

year.mat <- data.matrix(year3)

year.mat <- data.matrix(year4)

year.mat <- data.matrix(year5)

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

year4 <- subset(griz_num, Year == "year4")

year5 <- subset(griz_num, Year == "year5")


#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(year1$d13C, year1$d15N)
SE2 <- standard.ellipse(year2$d13C, year2$d15N)
SE3 <- standard.ellipse(year3$d13C, year3$d15N)
SE4 <- standard.ellipse(year4$d13C, year4$d15N)
SE5 <- standard.ellipse(year5$d13C, year5$d15N)

SE1
SE2
SE3
SE4
SE5

write.csv(SE1, "output/griz_SE_2018.csv")

write.csv(SE2, "output/griz_SE_2019.csv")

write.csv(SE3, "output/griz_SE_2020.csv")

write.csv(SE4, "output/griz_SE_2021.csv")

write.csv(SE5, "output/griz_SE_2022.csv")

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

year1_ <- rep("year1", length(SE1$xSEAc))
year2_ <- rep("year2", length(SE2$xSEAc))
year3_ <- rep("year3", length(SE3$xSEAc))
year4_ <- rep("year4", length(SE4$xSEAc))
year5_ <- rep("year5", length(SE5$xSEAc))

#create new data frame with names and ellipse outputs

Year <- c(year1_, year2_, year3_, year4_, year5_)
x <- c(SE1$xSEAc,SE2$xSEAc,SE3$xSEAc,SE4$xSEAc,SE5$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc,SE4$ySEAc,SE5$ySEAc)

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

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#correct prey by TEF values as in mixing models (Y = 5.28 + 0.88X, where Y = corrected, and X = uncorrected)
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +5.0 for meat (ants and mammals)

#dat.prey$d13C <- dat.prey$d13C + 3.7

#dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
#                                        ifelse(dat.prey$Source == "meat", 5.0, 
#                                               ifelse(dat.prey$Source == "kokanee", 4.1, 3.9)))

#correct 4 prey groups by TEF values as in mixing models (Erlenbach 2020)
#dC13 = +5.8 for all prey
#dN15 = +4.8 for all prey

dat.prey$d13C <- dat.prey$d13C + 5.8

dat.prey$d15N <- dat.prey$d15N + 4.8

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
  
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 21, 22, 23, 24),
                     breaks=c("year1", "year2", "year3", "year4", "year5",
                               "plants", "meat", "salmon", "kokanee"),
                     labels=c("2018", "2019", "2020", "2021", "2022",
                               "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_fill_manual(values = c("grey5","grey25","grey50", "grey65", "grey85",
                               "white", "white", "white", "white"),
                    breaks=c("year1", "year2", "year3", "year4", "year5",
                              "plants", "meat", "salmon", "kokanee"),
                    labels=c("2018", "2019", "2020", "2021", "2022",
                              "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_color_manual(values = c("grey5","grey25","grey50", "grey65", "grey85"),
                     breaks=c("year1", "year2", "year3", "year4", "year5"),
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
        legend.position = c(0.88, .275), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  scale_x_continuous(breaks = seq(-30, -12, 1)) +
  scale_y_continuous(breaks = seq(0, 17, 1))

tiff(file = "figures/grizz_year_240313.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()

#Plot GRIZZLY HAIR ellipses by Kokanee Spawning Time with food GROUPS (mean +/- SE) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204_spawning.csv", header = TRUE)

#subset into groups based on sampling date during Last 2 and 4 Weeks of Kokanee spawn

Last2Weeks <- subset(dat.hair, Last2Weeks == "y")

Last4Weeks <- subset(dat.hair, Last4Weeks == "y")

Rest <- subset(dat.hair, Last2Weeks == "n" & Last4Weeks == "n")

#Make new column "spawn_period" and merge into one dataframe

Last2Weeks$spawn_period <- "Last2Weeks"

Last4Weeks$spawn_period <- "Last4Weeks"

Rest$spawn_period <- "Rest"

spawn_merge <- rbind(Last2Weeks, Last4Weeks, Rest)

#Make into matrices for mshapiro test

spawn.mat <- data.matrix(Last2Weeks)

spawn.mat <- data.matrix(Last4Weeks)

spawn.mat <- data.matrix(Rest)

mshapiro.test(t(spawn.mat[,9:10]))

#none violated assumptions of normality

#next calculate standard.ellipse for each subset

SE1 <- standard.ellipse(Last2Weeks$d13C, Last2Weeks$d15N)
SE2 <- standard.ellipse(Last4Weeks$d13C, Last4Weeks$d15N)
SE3 <- standard.ellipse(Rest$d13C, Rest$d15N)

SE1
SE2
SE3

write.csv(SE1, "output/griz_SE_last2weeks.csv")

write.csv(SE2, "output/griz_SE_last4weeks.csv")

write.csv(SE3, "output/griz_SE_spawn_rest.csv")

#create name variable for each group (so ggplot2 knows how to colour it)
#name needs to be the same as original data frame that contains isotopic data

Last2Weeks_ <- rep("Last2Weeks", length(SE1$xSEAc))
Last4Weeks_ <- rep("Last4Weeks", length(SE2$xSEAc))
Rest_ <- rep("Rest", length(SE3$xSEAc))

#create new data frame with names and ellipse outputs

Spawn <- c(Last2Weeks_, Last4Weeks_, Rest_)
x <- c(SE1$xSEAc,SE2$xSEAc,SE3$xSEAc)
y <- c(SE1$ySEAc,SE2$ySEAc,SE3$ySEAc)

df_SE <- data.frame(Spawn, x, y)

plot(df_SE$x, df_SE$y)

#Combine grizzly raw data and ellipses into one dataframe
#add new column to both dataframe (Data type = Ellipses or Point)

spawn_merge$type <- 1

df_SE$type <- 2

#rename and re-order columns in df_SE to match spawn_merge

spawn_merge <- spawn_merge[, c("spawn_period", "d13C", "d15N", "type")]

colnames(df_SE) <- c("spawn_period", "d13C","d15N", "type")

#df_SE <- df_SE[c("Year", "d13C","d15N", "type")]

#merge both dataframes

spawn_both <- rbind(spawn_merge, df_SE)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

#correct prey by TEF values as in mixing models (Y = 5.28 + 0.88X, where Y = corrected, and X = uncorrected)
#dC13 = +3.7 for all prey
#dN15 = +5.5 for plants
#       +4.1 for kokanee
#       +3.9 for salmon
#       +5.0 for meat (ants and mammals)

dat.prey$d13C <- dat.prey$d13C + 3.7

dat.prey$d15N <- dat.prey$d15N + ifelse(dat.prey$Source == "plants", 5.5,
                                        ifelse(dat.prey$Source == "meat", 5.0, 
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
  
  geom_point(data = filter(spawn_both, type == 1), aes(x = d13C, y = d15N, fill = spawn_period, shape = spawn_period), 
             size = 2.5, alpha = 0.8) +
  
  #Add ellipses  
  
  geom_path(data = filter(spawn_both, type == 2), aes(x = d13C, y = d15N, 
                                                     color = spawn_period), linetype = 1, size = 1) +
  
  #Add prey means and error bars  
  
  geom_errorbar(data = prey_mean, aes(x = mean_C, ymin = mean_N - 1.96*se_N, 
                                      ymax = mean_N + 1.96*se_N), width = .2) +
  geom_errorbarh(data = prey_mean, aes(y = mean_N, xmin = mean_C - 1.96*se_C, 
                                       xmax = mean_C + 1.96*se_C), height =.2) +
  
  geom_point(data = prey_mean, aes(x = mean_C, y = mean_N, 
                                   shape = Source, fill = Source, stroke = 1), size = 2.5) +
  
  #Assign aesthetics and labels in legend
  
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 21, 22),
                     breaks=c("Last2Weeks", "Last4Weeks", "Rest",
                              "plants", "meat", "salmon", "kokanee"),
                     labels=c("Last 2 Weeks", "Last 4 Weeks", "Prior to Spawn",
                              "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_fill_manual(values = c("grey5","grey25","grey50", 
                               "white", "white", "white", "white"),
                    breaks=c("Last2Weeks", "Last4Weeks", "Rest",
                             "plants", "meat", "salmon", "kokanee"),
                    labels=c("Last 2 Weeks", "Last 4 Weeks", "Prior to Spawn",
                             "Plants", "Meat", "Salmon", "Kokanee")) +
  
  scale_color_manual(values = c("grey5","grey25","grey50"),
                     breaks=c("Last2Weeks", "Last4Weeks", "Rest"),
                     guide = FALSE) + 
  
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) +
  labs(title = "Grizzly by Spawning Period", subtitle = "and Food Groups (mean ± 95% CI)") +
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
        legend.position = c(0.88, .275), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) +
  
  scale_x_continuous(breaks = seq(-30, -12, 1)) +
  scale_y_continuous(breaks = seq(0, 17, 1))

tiff(file = "figures/grizz_spawn_240130.tif", 
     units = "cm", width = 25, height = 20, res = 300)

grizzly_w_prey

dev.off()


#Plot berry proportional contribution vs berry biomass----

rm(list = ls())

#import berry biomass table

dat.berry <- read.csv("data/berry_biomass_grizzly_diet_20240215.csv", header = TRUE)

dat.berry

berry_biomass <- ggplot(NULL) +
  
  #dat.berry, aes(x = berry_biomass, y = mean, 
   #                                    color = as.factor(ï..sample_year), shape = as.factor(ï..sample_year))) +
  
  #Add diet proportions means and error bars (95% CI)
  
  geom_errorbar(data = dat.berry, aes(x = berry_biomass, ymin = X5_perc, 
                                      ymax = X95_perc), width = 4) +
  
  geom_point(data = dat.berry, aes(x = berry_biomass, y = mean, 
                                   shape = as.factor(ï..sample_year), fill = as.factor(ï..sample_year), stroke = 1), size = 2.5) +
  
  #Assign aesthetics and labels in legend
  
  scale_shape_manual(values = c(22, 24, 21, 22),
                    #breaks=c("2019", "2020", "2021", "2022"),
                    labels=c("2019", "2020", "2021", "2022")) +
 
  scale_fill_manual(values = c("grey10","grey25", "grey45", "grey65"),
                    #breaks=c("2019", "2020", "2021", "2022"),
                    labels=c("2019", "2020", "2021", "2022")) +
  
  scale_x_continuous(breaks = seq(0, 300, 25), limits = c(0, 300), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1), expand = c(0, 0)) +
  
  ylab("Proportion Berries in Grizzly Diet") +
  xlab(expression(paste("Berry Biomass (g)"))) +
  labs(title = "Proportion of Berries (mean ± 95% CI)", 
       subtitle = expression(paste("and Berry Biomass (g)"))) +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=16, face = "bold"),
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

tiff(file = "figures/berry_biomass_grizz_diet_20240313.tif", 
     units = "cm", width = 25, height = 20, res = 300)

berry_biomass

dev.off()

#Plot kokanee proportional contribution vs spawning period----

rm(list = ls())

#import kokanee spawning period table

dat.kokanee <- read.csv("data/kokanee_spawning_grizzly_diet_20240215.csv", header = TRUE)

dat.kokanee

kokanee_spawning <- ggplot(data = dat.kokanee, aes(x = ï..Spawning_Period, y = mean)) +
  
  #Add diet proportions means and error bars (95% CI)
  
  geom_bar(stat = "identity", width = 0.5) +
  
  geom_errorbar(data = dat.kokanee, aes(ymin = X5, ymax = X95), width = 0.2) +
  
  scale_y_continuous(limits = c(0, 0.22), expand = c(0,0)) +
  
  scale_x_discrete(labels = c("Last 2 Weeks", "Last 4 Weeks", "Early in Spawn")) +
  
  #Assign aesthetics and labels in legend
  
  ylab("Proportion Kokanee in Grizzly Diet") +
  xlab("Kokanee Spawning Period") +
  labs(title = "Proportion of Kokanee (mean ± 95% CI)", 
       subtitle = "by Kokanee Spawning Period") +
  theme(plot.title = element_text(family="Times New Roman", size=16, face = "bold"),
        plot.subtitle = element_text(family="Times New Roman", size=16, face = "bold"),
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
  
  tiff(file = "figures/kokanee_spawning_grizz_diet_20240313.tif", 
       units = "cm", width = 25, height = 20, res = 300)

kokanee_spawning

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

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

#format Date as yyyy-mm-dd

dat.hair$Date <- as.Date(dat.hair$Date, format = "%d-%m-%Y")

#subset consumer data by Year

year1 <- subset(dat.hair, Year == '2018')

year2 <- subset(dat.hair, Year == '2019')

year3 <- subset(dat.hair, Year == '2020')

year4 <- subset(dat.hair, Year == '2021')

year5 <- subset(dat.hair, Year == '2022')

#add julian day - numeric

year1$jDate <- format(year1$Date, "%j")

year1$jDate <- as.numeric(year1$jDate)

year2$jDate <- format(year2$Date, "%j")

year2$jDate <- as.numeric(year2$jDate)

year3$jDate <- format(year3$Date, "%j")

year3$jDate <- as.numeric(year3$jDate)

year4$jDate <- format(year4$Date, "%j")

year4$jDate <- as.numeric(year4$jDate)

year5$jDate <- format(year5$Date, "%j")

year5$jDate <- as.numeric(year5$jDate)

#remove outlying value from 2019

#year2 <- year2[-c(29), ]

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

ggplot(year4, aes(x = jDate, y = d13C)) +
  geom_point() +
  stat_smooth()

ggplot(year4, aes(x = jDate, y = d15N)) +
  geom_point() +
  stat_smooth()

ggplot(year5, aes(x = jDate, y = d13C)) +
  geom_point() +
  stat_smooth()

ggplot(year5, aes(x = jDate, y = d15N)) +
  geom_point() +
  stat_smooth()

# test linear correlation between C13 and N15 and Time

cor(year1$d13C, year1$jDate)

cor(year1$d15N, year1$jDate)

cor(year2$d13C, year2$jDate)

cor(year2$d15N, year2$jDate)

cor(year3$d13C, year3$jDate)

cor(year3$d15N, year3$jDate)

cor(year4$d13C, year4$jDate)

cor(year4$d15N, year4$jDate)

cor(year5$d13C, year5$jDate)

cor(year5$d15N, year5$jDate)

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

lin.C4 <- lm(d13C ~ jDate, data = year4)

lin.C4

lin.N4 <- lm(d15N ~ jDate, data = year4)

lin.N4

lin.C5 <- lm(d13C ~ jDate, data = year5)

lin.C5

lin.N5 <- lm(d15N ~ jDate, data = year5)

lin.N5

#test significance of relationship

#summarize the model - p-values of intercept and c13 (x) should be ***

summary(lin.C1)

summary(lin.N1)

summary(lin.C2)

summary(lin.N2)

summary(lin.C3)

summary(lin.N3)

summary(lin.C4)

summary(lin.N4)

summary(lin.C5)

summary(lin.N5)

#calculate the 95% CI around the beta coefficients

confint(lin.C1)

confint(lin.N1)

confint(lin.C2)

confint(lin.N2)

confint(lin.C3)

confint(lin.N3)

confint(lin.C4)

confint(lin.N4)

confint(lin.C5)

confint(lin.N5)

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

c_time4 <- ggplot(year4, aes(jDate, d13C)) +
  geom_point(data = year4, aes(x = jDate, y = d13C)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"13C"), title = "2021") +
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

n_time4 <- ggplot(year4, aes(jDate, d15N)) +
  geom_point(data = year4, aes(x = jDate, y = d15N)) +
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

c_time5 <- ggplot(year5, aes(jDate, d13C)) +
  geom_point(data = year5, aes(x = jDate, y = d13C)) +
  geom_line(aes(group = ID), lty = 2) +
  stat_smooth(method = lm, color = "grey50") +
  labs(x = "Julian Date", y = expression(delta*"13C"), title = "2022") +
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

n_time5 <- ggplot(year5, aes(jDate, d15N)) +
  geom_point(data = year5, aes(x = jDate, y = d15N)) +
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

tiff("figures/sampling_dates_240117.tiff", units="cm", width=65, height=25, res=300)

ggarrange(c_time1, c_time2, c_time3, c_time4, c_time5, 
          n_time1, n_time2, n_time3, n_time4, n_time5,
          ncol = 5, nrow = 2)

dev.off()

#MixSIAR with Grizzly and 5 Sources OVERALL (normal) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

#Add prey data and plot as scatter plot

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

#discr <- load_discr_data(filename="data/grizzly_discrim_var11dec2023.csv", mix)

#load discrimination data (Erlenbach 2020)

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240227.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_plot_norm_240108", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

#plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

norm_griz24 <- "data/norm_griz_240108.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(norm_griz24, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.norm24 <- run_model(run = "normal", mix, source, discr, norm_griz24)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/norm24_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/norm24_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/norm24_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/norm24_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/norm24_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.norm24, mix, source, output_options)

# Save an object to a file

#save(jags.norm24, mix, source, output_options, file = "RData/all_grizzly_norm_240108.RData")

# Save an object to a file (Erlenbach 2020)

save(jags.norm24, mix, source, output_options, file = "RData/all_grizzly_erlenbach_240227.RData")

# Restore the objects

load(file = "RData/all_grizzly_erlenbach_240227.RData")

diag <- output_diagnostics(jags.norm24, mix, source, output_options)

names(diag)

head(diag$geweke)

df.stats <- output_stats(jags.norm24, mix, source, output_options)

write.csv(df.stats, "output/results_all_griz_erlenbach_240227.csv")

g.post <- output_posteriors(jags.norm24, mix, source, output_options)

all_griz <- g.post$global +
  scale_fill_manual(values=c("blue", "red", "green", "orange", "purple")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange", "purple")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "All Grizzly (n = 83)") +
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

#tiff(file = "figures/all_grizzly_240108.tif", 
#     units = "cm", width = 25, height = 20, res = 300)

tiff(file = "figures/all_grizzly_erlenbach_240227.tif", 
     units = "cm", width = 25, height = 20, res = 300)

all_griz

dev.off()

#MixSIAR with Grizzly and 5 Sources OVERALL (long) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

#Add prey data and plot as scatter plot

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

#discr <- load_discr_data(filename="data/grizzly_discrim_var11dec2023.csv", mix)

#load discrimination data (Erlenbach 2020)

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240227.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_plot_long_240301", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

#plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

long_griz24 <- "data/long_griz_240301.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(long_griz24, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.long24 <- run_model(run = "long", mix, source, discr, long_griz24)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/long24_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/long24_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/long24_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/long24_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/long24_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.long24, mix, source, output_options)

# Save an object to a file

#save(jags.norm24, mix, source, output_options, file = "RData/all_grizzly_norm_240108.RData")

# Save an object to a file (Erlenbach 2020)

save(jags.long24, mix, source, output_options, file = "RData/all_grizzly_long_erlenbach_240301.RData")

# Restore the objects

load(file = "RData/all_grizzly_long_erlenbach_240301.RData")

diag <- output_diagnostics(jags.long24, mix, source, output_options)

names(diag)

head(diag$geweke)

df.stats <- output_stats(jags.long24, mix, source, output_options)

write.csv(df.stats, "output/results_all_griz_long_erlenbach_240301.csv")

g.post <- output_posteriors(jags.long24, mix, source, output_options)

all_griz <- g.post$global +
  scale_fill_manual(values=c("blue", "red", "green", "orange", "purple")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange", "purple")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "All Grizzly (n = 83)") +
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
  
  #tiff(file = "figures/all_grizzly_240108.tif", 
  #     units = "cm", width = 25, height = 20, res = 300)
  
  tiff(file = "figures/all_grizzly_long_erlenbach_240301.tif", 
       units = "cm", width = 25, height = 20, res = 300)

all_griz

dev.off()

#MixSIAR with Grizzly and 4 Sources OVERALL (long) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

#Add prey data and plot as scatter plot

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

#discr <- load_discr_data(filename="data/grizzly_discrim_var11dec2023.csv", mix)

#load discrimination data (Erlenbach 2020)

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240301.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_plot4_long_240301", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

#plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

long_griz24 <- "data/long_griz4_240301.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(long_griz24, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.long24 <- run_model(run = "long", mix, source, discr, long_griz24)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/long24_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/long24_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/long24_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/long24_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/long24_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.long24, mix, source, output_options)

# Save an object to a file

#save(jags.norm24, mix, source, output_options, file = "RData/all_grizzly_norm_240108.RData")

# Save an object to a file (Erlenbach 2020)

save(jags.long24, mix, source, output_options, file = "RData/all_grizzly_long4_erlenbach_240301.RData")

# Restore the objects

load(file = "RData/all_grizzly_long_erlenbach4_240301.RData")

diag <- output_diagnostics(jags.long24, mix, source, output_options)

names(diag)

head(diag$geweke)

df.stats <- output_stats(jags.long24, mix, source, output_options)

write.csv(df.stats, "output/results_all_griz_long_erlenbach4_240301.csv")

g.post <- output_posteriors(jags.long24, mix, source, output_options)

all_griz <- g.post$global +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "All Grizzly (n = 83)") +
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
  
  #tiff(file = "figures/all_grizzly_240108.tif", 
  #     units = "cm", width = 25, height = 20, res = 300)
  
  tiff(file = "figures/all_grizzly_long_erlenbach4_240301.tif", 
       units = "cm", width = 25, height = 20, res = 300)

all_griz

dev.off()

#MixSIAR with Grizzly and 4 Sources OVERALL (normal) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

#remove grizzly, soapberry, and balance out meat group by subsampling

#dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

#Add prey data and plot as scatter plot

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

#discr <- load_discr_data(filename="data/grizzly_discrim_var11dec2023.csv", mix)

#load discrimination data (Erlenbach 2020)

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240301.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_plot4_norm_240301", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

#plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

norm_griz24 <- "data/norm_griz4_240301.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(norm_griz24, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.norm24 <- run_model(run = "normal", mix, source, discr, norm_griz24)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/norm24_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/norm24_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/norm24_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/norm24_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/norm24_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.norm24, mix, source, output_options)

# Save an object to a file

#save(jags.norm24, mix, source, output_options, file = "RData/all_grizzly_norm_240108.RData")

# Save an object to a file (Erlenbach 2020)

save(jags.norm24, mix, source, output_options, 
     file = "RData/all_grizzly_norm_erlenbach4_nosubsampling_240305.RData")

# Restore the objects

load(file = "RData/all_grizzly_norm_erlenbach4_nosubsampling_240305.RData")

diag <- output_diagnostics(jags.norm24, mix, source, output_options)

names(diag)

head(diag$geweke)

df.stats <- output_stats(jags.norm24, mix, source, output_options)

write.csv(df.stats, "output/results_all_griz_norm_erlenbach4_240301.csv")

g.post <- output_posteriors(jags.norm24, mix, source, output_options)

all_griz <- g.post$global +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "All Grizzly (n = 83)") +
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
  
  #tiff(file = "figures/all_grizzly_240108.tif", 
  #     units = "cm", width = 25, height = 20, res = 300)
  
  tiff(file = "figures/all_grizzly_norm_erlenbach4_nosubsampling_240305.tif", 
       units = "cm", width = 25, height = 20, res = 300)

all_griz

dev.off()

#MixSIAR with Grizzly and 4 Sources OVERALL - No Salmon, but Berry included (normal) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE)

#remove grizzly, salmon, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147, 148:181), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_20240313.csv", row.names = FALSE)

#Add prey data and plot as scatter plot

source <- load_source_data(filename = "data/knn_sources_20240313.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

#discr <- load_discr_data(filename="data/grizzly_discrim_var11dec2023.csv", mix)

#load discrimination data (Erlenbach 2020)

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240313.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_plot4_norm_240313", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

#plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

norm_griz24 <- "data/norm_griz4_240313.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(norm_griz24, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.norm24 <- run_model(run = "normal", mix, source, discr, norm_griz24)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/norm24_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/norm24_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/norm24_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/norm24_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/norm24_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.norm24, mix, source, output_options)

# Save an object to a file

#save(jags.norm24, mix, source, output_options, file = "RData/all_grizzly_norm_240108.RData")

# Save an object to a file (Erlenbach 2020)

save(jags.norm24, mix, source, output_options, 
     file = "RData/all_grizzly_norm_erlenbach_nosalmon_240313.RData")

# Restore the objects

load(file = "RData/all_grizzly_norm_erlenbach_nosalmon_240313.RData")

diag <- output_diagnostics(jags.norm24, mix, source, output_options)

names(diag)

head(diag$geweke)

df.stats <- output_stats(jags.norm24, mix, source, output_options)

write.csv(df.stats, "output/results_all_griz_norm_erlenbach_nosalmon_240313.csv")

g.post <- output_posteriors(jags.norm24, mix, source, output_options)

all_griz <- g.post$global +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "All Grizzly (n = 83)") +
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
  
  #tiff(file = "figures/all_grizzly_240108.tif", 
  #     units = "cm", width = 25, height = 20, res = 300)
  
  tiff(file = "figures/all_grizzly_norm_erlenbach_nosalmon_240313.tif", 
       units = "cm", width = 25, height = 20, res = 300)

all_griz

dev.off()


#MixSIAR with Grizzly and 4 Sources by SEX (long) ----

#Take average d13C and d15N values of individuals who were sampled > 1

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)
#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240301.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_sex_240306", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

sex_griz_long24 <- "data/sex_griz_long24.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(sex_griz_long24, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.sex24 <- run_model(run = "long", mix, source, discr, sex_griz_long24)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/sex_long_summary_stats24",
                       sup_post = FALSE,
                       plot_post_save_pdf = FALSE,
                       plot_post_name = "output/sex_long_posterior_density24",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = FALSE,
                       plot_pairs_name = "output/sex_long_pairs_plot24",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "output/sex_long_xy_plot24",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/sex_long_diagnostics24",
                       indiv_effect = FALSE,
                       plot_post_save_png = TRUE,
                       plot_pairs_save_png = TRUE,
                       plot_xy_save_png = TRUE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.sex24, mix, source, output_options)

# Save an object to a file

save(jags.sex24, mix, source, output_options, file = "RData/sex_grizzly_long24.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/sex_grizzly_long24.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.sex24, mix, source, output_options)

write.csv(df.stats, "output/results_sex_griz_long_240308.csv")

source$source_names

#Do males eat more meat than females? 

# Shapiro-Wilk normality test for meat

shapiro.test(jags.sex24$BUGSoutput$sims.list$p.fac1[,2,2]) # p << 0.005

shapiro.test(jags.sex24$BUGSoutput$sims.list$p.fac1[,1,2]) # p << 0.005

male.meat <- as.vector(jags.sex24$BUGSoutput$sims.list$p.fac1[,2,2])

female.meat <- as.vector(jags.sex24$BUGSoutput$sims.list$p.fac1[,1,2])

res <- wilcox.test(male.meat, female.meat,
                   exact = FALSE, alternative = "greater")

res

#Do females eat more plants than males?

# Shapiro-Wilk normality test for plants

shapiro.test(jags.sex24$BUGSoutput$sims.list$p.fac1[,2,3]) # p << 0.005

shapiro.test(jags.sex24$BUGSoutput$sims.list$p.fac1[,1,3]) # p > 0.05

male.plants <- as.vector(jags.sex24$BUGSoutput$sims.list$p.fac1[,2,3])

female.plants <- as.vector(jags.sex24$BUGSoutput$sims.list$p.fac1[,1,3])

res <- wilcox.test(male.plants,female.plants, 
                   exact = FALSE, alternative = "greater")

res

#Do females eat more kokanee than males?

# Shapiro-Wilk normality test for plants

shapiro.test(jags.sex24$BUGSoutput$sims.list$p.fac1[,2,4]) # p << 0.005

shapiro.test(jags.sex24$BUGSoutput$sims.list$p.fac1[,1,4]) # p << 0.005

male.kokanee <- as.vector(jags.sex24$BUGSoutput$sims.list$p.fac1[,2,4])

female.kokanee <- as.vector(jags.sex24$BUGSoutput$sims.list$p.fac1[,1,4])

res <- wilcox.test(male.kokanee,female.kokanee, 
                   exact = FALSE, alternative = "greater")

res

#plot histo of difference between proportion of meat between Females and Males

p.meat.F.M <- jags.sex24$BUGSoutput$sims.list$p.fac1[,2,2] - jags.sex24$BUGSoutput$sims.list$p.fac1[,1,2]

p.meat.F.M <- as.data.frame(p.meat.F.M)

histo.meat <- ggplot(p.meat.F.M, aes(p.meat.F.M)) +
 
  geom_histogram(binwidth = 0.01) +
  
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 2000, 100), limits = c(0, 2000), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Meat proportions, Male - Female", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = NULL)

#tiff(file = "figures/sex_meat_histo_240308_4source.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#histo.meat

#dev.off()
  
#plot histo of difference between proportion of plants between Males and Females

p.plants.M.F <- jags.sex24$BUGSoutput$sims.list$p.fac1[,2,3] - jags.sex24$BUGSoutput$sims.list$p.fac1[,1,3]

p.plants.M.F <- as.data.frame(p.plants.M.F)

histo.plants <- ggplot(p.plants.M.F, aes(p.plants.M.F)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 800, 50), limits = c(0, 800), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Plant proportions, Male - Female", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = NULL)

#tiff(file = "figures/sex_plants_meat_histo_240308_4source.tif", 
#     units = "cm", width = 40, height = 60, res = 300)

#ggarrange(histo.meat, histo.plants,
#          ncol = 1, nrow = 2)

#dev.off()

#plot histo of difference between proportion of kokanee between Females and Males

p.kokanee.F.M <- jags.sex24$BUGSoutput$sims.list$p.fac1[,1,1] - jags.sex24$BUGSoutput$sims.list$p.fac1[,2,1]

p.kokanee.F.M <- as.data.frame(p.kokanee.F.M)

histo.kokanee <- ggplot(p.kokanee.F.M, aes(p.kokanee.F.M)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 800, 50), limits = c(0, 800), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Kokanee proportions, Female - Male", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = NULL)

#plot histo of difference between proportion of salmon between males and females

p.salmon.M.F <- jags.sex24$BUGSoutput$sims.list$p.fac1[,2,4] - jags.sex24$BUGSoutput$sims.list$p.fac1[,1,4]

p.salmon.M.F <- as.data.frame(p.salmon.M.F)

histo.salmon <- ggplot(p.salmon.M.F, aes(p.salmon.M.F)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 800, 50), limits = c(0, 800), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Salmon proportions, Male - Female", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = NULL)

tiff(file = "figures/sex_histos_240313.tif", 
     units = "cm", width = 40, height = 60, res = 300)

ggarrange(histo.plants, histo.kokanee,
          ncol = 1, nrow = 2)

dev.off()

#Plot posterior densities of females and males

g.post <- output_posteriors(jags.sex24, mix, source, output_options)

female_griz <- g.post$fac1[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Female Grizzly (n = 48)") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "18")) 
  
#  tiff(file = "figures/female_grizzly_240112_4source.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#female_griz

#dev.off()

male_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Male Grizzly (n = 35)") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "18")) +
  
  tiff(file = "figures/male_female_grizzly_240308_4source.tif", 
       units = "cm", width = 40, height = 60, res = 300)

ggarrange(female_griz, male_griz,
          ncol = 1, nrow = 2)

dev.off()

#MixSIAR with Grizzly and 4 Sources by Spawning Period (long) ----

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204_spawning.csv", header = TRUE)

#subset into groups based on sampling date during Last 2 and 4 Weeks of Kokanee spawn

Last2Weeks <- subset(dat.hair, Last2Weeks == "y")

Last4Weeks <- subset(dat.hair, Last4Weeks == "y")

Rest <- subset(dat.hair, Last2Weeks == "n" & Last4Weeks == "n")

#Make new column "spawn_period" and merge into one dataframe

Last2Weeks$spawn_period <- "Last2Weeks"

Last4Weeks$spawn_period <- "Last4Weeks"

Rest$spawn_period <- "Rest"

spawn_merge <- rbind(Last2Weeks, Last4Weeks, Rest)

write.csv(spawn_merge, "data/spawn_merge.csv")

# Load the mixture/consumer data

mix <- load_mix_data(filename = "data/spawn_merge.csv",
                     iso_names = c("d13C","d15N"), 
                     factors = c("spawn_period"), 
                     fac_random = c(FALSE), 
                     fac_nested = c(FALSE), 
                     cont_effects = NULL)

#load_source_data

dat.prey <- read.csv("data/knn_sources2_06apr2021.csv", header = TRUE)

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240301.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_spawn_240130", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

spawn_griz_long <- "data/spawn_griz_long.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(spawn_griz_long, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.spawn <- run_model(run = "long", mix, source, discr, spawn_griz_long)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/spawn_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "output/spawn_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/spawn_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/spawn_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/spawn_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.spawn, mix, source, output_options)

# Save an object to a file

save(jags.spawn, mix, source, output_options, file = "RData/spawn_grizzly_long.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/spawn_grizzly_long.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.spawn, mix, source, output_options)

write.csv(df.stats, "output/results_spawn_griz_long_240309.csv")

#Do kokanee proportions differ between spawning periods? 

# Shapiro-Wilk normality test for kokanee

shapiro.test(jags.spawn$BUGSoutput$sims.list$p.fac1[,1,1]) # p << 0.005

shapiro.test(jags.spawn$BUGSoutput$sims.list$p.fac1[,2,1]) # p << 0.005

shapiro.test(jags.spawn$BUGSoutput$sims.list$p.fac1[,3,1]) # p << 0.005

Last2Weeks.kokanee <- as.vector(jags.spawn$BUGSoutput$sims.list$p.fac1[,1,1])

Last4Weeks.kokanee <- as.vector(jags.spawn$BUGSoutput$sims.list$p.fac1[,2,1])

Rest.kokanee <- as.vector(jags.spawn$BUGSoutput$sims.list$p.fac1[,3,1])

res <- wilcox.test(Last2Weeks.kokanee, Last4Weeks.kokanee,
                   exact = FALSE, alternative = "greater")

res

#proportion of kokanee in Last2Weeks.kokanee > Last4Weeks.kokanee

res <- wilcox.test(Last4Weeks.kokanee, Rest.kokanee,
                   exact = FALSE, alternative = "greater")

res

#proportion of kokanee in Last4Weeks.kokanee > Rest.kokanee

res <- wilcox.test(Last2Weeks.kokanee, Rest.kokanee,
                   exact = FALSE, alternative = "greater")

res

#proportion of kokanee in Last2Weeks.kokanee > Rest.kokanee

#plot histo of difference between proportion of kokanee between Last2Weeks.kokanee and Last4Weeks.kokanee

p.kokanee.2.4 <- jags.spawn$BUGSoutput$sims.list$p.fac1[,1,1] - jags.spawn$BUGSoutput$sims.list$p.fac1[,2,1]

p.kokanee.2.4 <- as.data.frame(p.kokanee.2.4)

histo.kokanee.2.4 <- ggplot(p.kokanee.2.4, aes(p.kokanee.2.4)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 1100, 50), limits = c(0, 1100), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Kokanee proportions, Last 2 Weeks of Spawn - Last 4 Weeks of Spawn", 
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

tiff(file = "figures/spawn_kokanee_2_4_histo_240309.tif", 
     units = "cm", width = 25, height = 20, res = 300)

histo.kokanee.2.4

dev.off()

#plot histo of difference between proportion of kokanee between Last 2 Weeks and Rest

p.kokanee.2.R <- jags.spawn$BUGSoutput$sims.list$p.fac1[,1,1] - jags.spawn$BUGSoutput$sims.list$p.fac1[,3,1]

p.kokanee.2.R <- as.data.frame(p.kokanee.2.R)

histo.kokanee.2.R <- ggplot(p.kokanee.2.R, aes(p.kokanee.2.R)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 900, 50), limits = c(0, 900), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Kokanee proportions, Last 2 Weeks of Spawn - Prior to Spawn", 
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

tiff(file = "figures/spawn_kokanee_2_R_histo_240309.tif", 
     units = "cm", width = 25, height = 20, res = 300)

histo.kokanee.2.R

dev.off()

#Plot posterior densities by kokanee spawning period

g.post <- output_posteriors(jags.spawn, mix, source, output_options)

last2weeks_griz <- g.post$fac1[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 0.10, 0.01), limits = c(0, 0.10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Grizzly Sampled Last 2 Weeks of Kokanee Spawn (n = 63)") +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))
  
#  tiff(file = "figures/Last2Weeks_grizzly_240201.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#last2weeks_griz

#dev.off()

last4weeks_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 0.10, 0.01), limits = c(0, 0.10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Grizzly Sampled Last 4 Weeks of Kokanee Spawn (n = 128)") +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))
  
#  tiff(file = "figures/last4weeks_grizzly_240201.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#last4weeks_griz

#dev.off()

rest_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 0.10, 0.01), limits = c(0, 0.10), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Grizzly Sampled Prior to Kokanee Spawn (n = 94)") +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) 
  
#  tiff(file = "figures/rest_grizzly_240201.tif", 
#       units = "cm", width = 400, height = 400, res = 300)

tiff(file = "figures/spawning_period_grizzly_240309.tif", 
        units = "cm", width = 40, height = 40, res = 300)  

ggarrange(last2weeks_griz, last4weeks_griz, rest_griz,
          ncol = 2, nrow = 2)

dev.off()

#MixSIAR with Grizzly and 4 Sources by AREA (long) ----

#Take average d13C and d15N values of individuals who were sampled > 1

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240301.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_area_240305", plot_save_pdf = FALSE,
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

save(jags.area, mix, source, output_options, file = "RData/area_grizzly_long_erlenbach.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/area_grizzly_long_erlenbach.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.area, mix, source, output_options)

write.csv(df.stats, "output/results_area_griz_long_erlenbach_240305.csv")

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
  
  scale_y_continuous(breaks = seq(0, 1300, 100), limits = c(0, 1500), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Meat proportions, South - North", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = NULL)

#tiff(file = "figures/area_meat_histo_erlenbach_240305.tif", 
#     units = "cm", width = 25, height = 20, res = 300)

#histo.meat

#dev.off()

#plot histo of difference between proportion of plants between North and South

p.plants.N.S <- jags.area$BUGSoutput$sims.list$p.fac1[,1,3] - jags.area$BUGSoutput$sims.list$p.fac1[,2,3]

p.plants.N.S <- as.data.frame(p.plants.N.S)

histo.plants <- ggplot(p.plants.N.S, aes(p.plants.N.S)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 600, 50), limits = c(0, 600), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Plant proportions, North - South", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = NULL)

tiff(file = "figures/area_meat_plants_histo_erlenbach_240305.tif", 
     units = "cm", width = 40, height = 60, res = 300)

ggarrange(histo.meat, histo.plants,
          ncol = 1, nrow = 2)

dev.off()

#plot histo of difference between proportion of kokanee between Females and Males

p.kokanee.S.N <- jags.area$BUGSoutput$sims.list$p.fac1[,2,1] - jags.area$BUGSoutput$sims.list$p.fac1[,1,1]

p.kokanee.S.N <- as.data.frame(p.kokanee.S.N)

histo.kokanee <- ggplot(p.kokanee.S.N, aes(p.kokanee.S.N)) +
  
  geom_histogram(binwidth = 0.01) +
  
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  
  scale_x_continuous(breaks = seq(-0.32, 0.32, 0.1), limits = c(-0.32, 0.32), expand = c(0, 0)) +
  
  scale_y_continuous(breaks = seq(0, 800, 50), limits = c(0, 800), expand = c(0, 0)) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  
  labs(title = "Difference between Kokanee proportions, South - North", 
       x = "Difference", y = "Frequency") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = NULL)

tiff(file = "figures/area_kokanee_histo_240313.tif", 
     units = "cm", width = 25, height = 20, res = 300)

histo.kokanee

dev.off()

#Plot posterior densities by area

g.post <- output_posteriors(jags.area, mix, source, output_options)

north_griz <- g.post$fac1[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "North Grizzly (n = 49)") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "18"))
  
#  tiff(file = "figures/north_grizzly_erlenbach_240305.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#north_griz

#dev.off()

south_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "South Grizzly (n = 34)") +
  theme(plot.title = element_text(family="Times New Roman", size=20, face = "bold"),
        axis.title.x = element_text(family="Times New Roman", size=16),
        axis.title.y = element_text(family="Times New Roman", size=16),
        plot.margin = margin(0.5, 2, 0.5, 0.5, "cm"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x= element_text(family = "Times New Roman", colour = "black", size = "16"),
        axis.text.y= element_text(family = "Times New Roman", colour = "black", size = "16"),
        legend.position = c(1.05, 1.052), legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "18"))
  
tiff(file = "figures/area_both_grizzly_erlenbach_240305.tif", 
       units = "cm", width = 40, height = 60, res = 300)

ggarrange(north_griz, south_griz,
          ncol = 1, nrow = 2)

dev.off()

#MixSIAR with Grizzly and 4 Sources by Year (long) ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2021", "year4", as.character(Year)))

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2022", "year5", as.character(Year)))

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

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240301.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_year_240118", plot_save_pdf = FALSE,
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

write.csv(df.stats, "output/results_year_griz_long_240310.csv")

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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))
#year1_griz <- year1_griz + guides(fill=guide_legend(title=))  

#tiff(file = "figures/2018_grizzly_240118.tif", 
       #units = "cm", width = 25, height = 20, res = 300)

#year1_griz

#dev.off()

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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))
  
#  tiff(file = "figures/2019_grizzly_240118.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#year2_griz

#dev.off()

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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))
  
#  tiff(file = "figures/2020_grizzly_240118.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#year3_griz

year4_griz <- g.post$fac1[[4]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "2021 Grizzly (n = 37)") +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) 
  
#  tiff(file = "figures/2021_grizzly_240118.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#year4_griz

#dev.off()

year5_griz <- g.post$fac1[[5]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "2022 Grizzly (n = 28)") +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) 
  
tiff(file = "figures/year_grizzly_2018_thru_2022_240310.tif", 
     units = "mm", width = 400, height = 600, res = 300)

ggarrange(year1_griz, year2_griz, year3_griz, year4_griz, year5_griz,
          ncol = 2, nrow = 3)

dev.off()

#MixSIAR with Grizzly and 5 Sources by Year (long) ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2021", "year4", as.character(Year)))

griz_merge <- griz_merge %>% 
  mutate(Year = ifelse(as.character(Year) == "2022", "year5", as.character(Year)))

griz_merge <- subset(griz_merge, Year!="year1")


write.csv(griz_merge, "data/griz_merge.csv")

# Load the mixture/consumer data

mix <- load_mix_data(filename = "data/griz_merge.csv",
                     iso_names = c("d13C","d15N"), 
                     factors = c("Year"), 
                     fac_random = c(FALSE), 
                     fac_nested = c(FALSE), 
                     cont_effects = NULL)

#Add prey data and plot as scatter plot

dat.prey <- read.csv("data/knn_sources_04dec2023.csv", header = TRUE, sep = ",")

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240227.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_year_berry_240310", plot_save_pdf = FALSE,
          plot_save_png = TRUE, mix, source, discr)

#Calculate normalized surface area

if(mix$n.iso == 2) calc_area(source = source, mix = mix, discr = discr)

#Plot Uninformative prior

plot_prior(alpha.prior = 1, source)

# Write the JAGS model file

year_berry_griz_long <- "data/year_berry_griz_long.txt"   # Name of the JAGS model file

resid_err <- TRUE

process_err <- TRUE

write_JAGS_model(year_berry_griz_long, resid_err, process_err, mix, source)

#Run model with MCMC run options

jags.year.berry <- run_model(run = "long", mix, source, discr, year_berry_griz_long)

#set output options like file names, plot file types, etc. (see ?output_JAGS for details)

output_options <- list(summary_save = TRUE,
                       summary_name = "output/year_berry_summary_stats",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "output/year_berry_posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "output/year_berry_pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = TRUE,
                       plot_xy_name = "output/year_berry_xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "output/year_berry_diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = TRUE,
                       return_obj = TRUE)

#Call output_JAGS to process diagnostics, summary statistics, and create posterior density plots

output_JAGS(jags.year.berry, mix, source, output_options)

# Save an object to a file

save(jags.year.berry, mix, source, output_options, file = "RData/year_berry_grizzly_long.RData")

rm(list = ls())

# Restore the objects

load(file = "RData/year_berry_grizzly_long.RData")

#save summary stats as dataframe

df.stats <- output_stats(jags.year.berry, mix, source, output_options)

write.csv(df.stats, "output/results_year_berry_griz_long_240310.csv")

#plot posterior densities by year

g.post <- output_posteriors(jags.year.berry, mix, source, output_options)

year2_griz <- g.post$fac1[[1]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange", "purple")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange", "purple")) +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))
#year1_griz <- year1_griz + guides(fill=guide_legend(title=))  

#tiff(file = "figures/2018_grizzly_240118.tif", 
#units = "cm", width = 25, height = 20, res = 300)

#year1_griz

#dev.off()

year3_griz <- g.post$fac1[[2]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange", "purple")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange", "purple")) +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))

#  tiff(file = "figures/2019_grizzly_240118.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#year2_griz

#dev.off()

year4_griz <- g.post$fac1[[3]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange", "purple")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange", "purple")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "2021 Grizzly (n = 37)") +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14"))

#  tiff(file = "figures/2020_grizzly_240118.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#year3_griz

year5_griz <- g.post$fac1[[4]] +
  scale_fill_manual(values=c("blue", "red", "green", "orange", "purple")) + 
  scale_color_manual(values=c("blue", "red", "green", "orange", "purple")) +
  scale_x_continuous(breaks = seq(0, 1.25, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.0, 0.25), limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "2022 Grizzly (n = 28)") +
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
        legend.text = element_text(family = "Times New Roman", colour = "black", size = "14")) 

#  tiff(file = "figures/2021_grizzly_240118.tif", 
#       units = "cm", width = 25, height = 20, res = 300)

#year4_griz

#dev.off()

tiff(file = "figures/year_berry_grizzly_2019_thru_2022_240310.tif", 
     units = "mm", width = 400, height = 400, res = 300)

ggarrange(year2_griz, year3_griz, year4_griz, year5_griz,
          ncol = 2, nrow = 2)

dev.off()


#MixSIAR with Grizzly and 4 Sources by Sampling Date (long) ----

#Convert Date in Grizzly Data

rm(list = ls())

dat.hair <- read.csv("data/grizzly_all_dates_231204.csv", header = TRUE)

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

#remove grizzly, soapberry, and balance out meat group by subsampling

dat.prey <- dat.prey[-c(13, 26, 30:33, 36:68, 70:72, 147), ]

#remove grizzly and soapberry - NO SUBSAMPLING

#dat.prey <- dat.prey[-c(13, 26, 31:33), ]

write.csv(dat.prey, "data/knn_sources_04dec2023_2.csv", row.names = FALSE)

source <- load_source_data(filename = "data/knn_sources_04dec2023_2.csv",
                           source_factors = NULL,
                           conc_dep = FALSE,
                           data_type = "raw", mix)

#load discrimination data

discr <- load_discr_data(filename="data/grizzly_discrim_erlenbach_20240301.csv", mix)

#plot_data

plot_data(filename = "figures/isospace_date_240122", plot_save_pdf = FALSE,
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
  labs(title = "Grizzly Diet over Sampling Date (n = 222)", x = "Julian Date") +
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
  
  tiff(file = "figures/date_grizzly_240312.tif", 
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

