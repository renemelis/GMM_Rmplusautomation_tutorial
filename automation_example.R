#install.packages("dplyr")
#install.packages("broom.mixed")
#install.packages("MplusAutomation")
#install.packages("wesanderson")
#install.packages("RColorBrewer")

#install.packages("BiocManager")
#install.packages("rhdf5")

library(dplyr)
library(tidyverse)
library(questionr)
library(lubridate)
library(gtsummary)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(MplusAutomation)
library(gt)
library(rlist)
library(stringr)
library(rhdf5)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(ggpubr)
library(gridExtra)

#This script intends to demonstrate how R can be used to built mplus models 
#(mplusObject and mplusModeler for single models, createModels for multiple models (e.g. same model with different outcome)), 
#run them (runModels)
#and "read" the output files into r using readModels 
#from the mplusautomation packages.

#I first demonstrate how the mplusautomation package can be used on a single model
#next, I demonstrate how it can be used to create, run and report multiple models

#I used a data file from the lcmm package with repeated measurements of 500 persons with dementia across follow ups from the paquid cohort.
#I first ran an unconditional (linear) growth model in LMER.

#I followed the steps suggested by Jung and Wickrama, ie
## LGCA
## LCGA
## GMM

#I wrote script to report:
## fit statistics + warning
## parameter estimates for both random and fixed part of the models
## class counts/proportions for the models with growth classes
## plots
# and to:
## extract participant level most likely class membership information from the GH5 file (that is the file mplus uses for plotting) 

####some objects/variables that are used throughout this script that I define first:
#rel.warn is an object with relevant warnings. rel.warn can be used to select these from the list of warning provided when using readModels function
rel.warn <- 
  c("SOLUTION MAY NOT BE TRUSTWORTHY DUE TO LOCAL MAXIMA",
    "THE COVARIANCE COVERAGE FALLS BELOW THE SPECIFIED LIMIT", 
    "THE BEST LOGLIKELIHOOD VALUE HAS BEEN REPLICATED",
    "IS NOT POSITIVE DEFINITE",
    "IS NOT POSITIVE",
    "WARNING:  THE BEST LOGLIKELIHOOD VALUE WAS NOT REPLICATED.  THE",
    "SOLUTION MAY NOT BE TRUSTWORTHY DUE TO LOCAL MAXIMA.  INCREASE THE",
    "TRUSTWORTHY FOR SOME PARAMETERS DUE TO A NON-POSITIVE DEFINITE", 
    "ONE OR MORE PARAMETERS WERE FIXED TO AVOID SINGULARITY OF THE",
    "INFORMATION MATRIX. THE SINGULARITY IS MOST LIKELY BECAUSE THE",
    "MODEL IS NOT IDENTIFIED, OR BECAUSE OF EMPTY CELLS IN THE JOINT",
    "IS NOT POSITIVE DEFINITE.  THIS COULD INDICATE A NEGATIVE VARIANCE/")
rel.warn <- paste(rel.warn, collapse = "|")
#defining colors for the plot, this object can be extended with more colors.

#defining colors for the plot, this object can be extended with more colors.
custom.col <- c("#FC4E07","#4E84C4","#E69F00", "#00AFBB", "#293352") 

pq <- lcmm::paquid 
#paquid is a demonstrator data frame that comes with the lcmm package

head(pq)
pq <- pq %>%
  mutate(demfu = age-agedem)

pq <- pq %>%
  arrange(ID, age) %>%
  group_by(ID) %>%
  mutate(funr = row_number()-1) %>%
  ungroup() %>%
  mutate(funr.sq = funr^2) %>%
  mutate(demfu.sq = demfu^2)

pqb <- pq %>%
  filter (funr == 1)

sp.mmse <- ggplot(data=pq, aes(x=funr, y=MMSE, group=ID))
sp.mmse+geom_line()

sp.cesd <- ggplot(data=pq, aes(x=funr, y=CESD, group=ID))
sp.cesd+geom_line()

um<- ##unconditional means model
  lmer(MMSE~1+
         (1|ID), data=pq, 
       na.action=na.omit, REML=FALSE) 
summary(um)

ul<-##unconditional linear model MMSE
  lmer(MMSE~funr+
         (funr|ID), data=pq, 
       na.action=na.omit, REML=FALSE) 
summary(ul)

uq<-##unconditional quadratic model
  lmer(MMSE~funr+funr.sq+
         (funr+funr.sq|ID), data=pq, 
       na.action=na.omit, REML=FALSE) 
summary(uq)

anova(ul, um)
#why is dif in deviance between ul and um the same as in mplus, but are their absolute values different? 

pq.w <- pq %>%
  select(ID, MMSE, CESD, funr) %>%
  pivot_wider(names_from = funr,
              values_from = c(MMSE, CESD))

setwd("Z:/mplusautomationtutorial/in-outputfiles")
#setwd("Z:/mplusautomationtutorial/in-outputfiles")

#######focusing on single models
#######latent growth curve modeling

#same model with time point (0, 1, 2 etc) in SEM with coefficients fixed at values
ul.mp <- mplusObject(
  TITLE=          "unconditional linear model MMSE;",
  MODEL=          "i s | MMSE_0@0 MMSE_1@1 MMSE_2@2 MMSE_3@3 MMSE_4@4 MMSE_5@5 MMSE_6@6 MMSE_7@7 MMSE_8@8;
                  MMSE_0(1);
                  MMSE_1(1);
                  MMSE_2(1);
                  MMSE_3(1);
                  MMSE_4(1);
                  MMSE_5(1);
                  MMSE_6(1);
                  MMSE_7(1);
                  MMSE_8(1);",
  PLOT=
    "Type = PLOT3;
    Series =  MMSE_0 (0) MMSE_1 (1) MMSE_2 (2) MMSE_3 (3) 
    MMSE_4 (4) MMSE_5 (5) MMSE_6 (6) MMSE_7 (7);",
  rdata = pq.w
)

ul.mp <- mplusModeler(ul.mp, 
                      modelout="ul.inp", 
                      run=1L)

###### extracting from single mplus output files 

#reading the mplus output file using readModels: turns mplus output file into a list
#you can then extract the right data element from that list.
ul.mp.l <- readModels("Z:/mplusautomationtutorial/in-outputfiles/ul.out",
                            recursive = F)

#for reporting lgca models from mplus in r the following elements from the above list are needed:

##PLOT the estimated average trajectory
##ul.mp.l[["gh5"]][["means_and_variances_data"]]: this element holds the estimated mean y score at each of the time points used (nine here)
#$y_estimated_means
#$y_estimated_means$values
#[,1]
#[1,] 27.45963
#[2,] 26.30872
#[3,] 25.15781
#[4,] 24.00691
#[5,] 22.85600
#[6,] 21.70510
#[7,] 20.55419
#[8,] 19.40329
#[9,] 18.25238

#if you plot the above values (on y-scale) against your indicator of time (0-8 in this case) on the x-axis you have a plot for estimated growth model

##FIT statistics
##ul.mp.l[["summaries"]]:

#  Mplus.version                             Title AnalysisType   DataType Estimator Observations NGroups NDependentVars NIndependentVars NContinuousLatentVars Parameters
#1           8.4  unconditional linear model MMSE;      GENERAL INDIVIDUAL        ML          500       1              9                0                     2          6
#LL      AIC      BIC     aBIC     AICC Filename
#1 -9217.874 18447.75 18473.03 18453.99 18447.92   ul.out

##PARAMETER estimates for random and fixed effects
##ul.mp.l[["parameters"]][["unstandardized"]]

#paramHeader  param    est    se  est_se pval
#1                 I.| MMSE_0  1.000 0.000 999.000  999
#2                 I.| MMSE_1  1.000 0.000 999.000  999
#3                 I.| MMSE_2  1.000 0.000 999.000  999
#4                 I.| MMSE_3  1.000 0.000 999.000  999
#5                 I.| MMSE_4  1.000 0.000 999.000  999
#6                 I.| MMSE_5  1.000 0.000 999.000  999
#7                 I.| MMSE_6  1.000 0.000 999.000  999
#8                 I.| MMSE_7  1.000 0.000 999.000  999
#9                 I.| MMSE_8  1.000 0.000 999.000  999
#10                S.| MMSE_0  0.000 0.000 999.000  999
#11                S.| MMSE_1  1.000 0.000 999.000  999
#12                S.| MMSE_2  2.000 0.000 999.000  999
#13                S.| MMSE_3  3.000 0.000 999.000  999
#14                S.| MMSE_4  4.000 0.000 999.000  999
#15                S.| MMSE_5  5.000 0.000 999.000  999
#16                S.| MMSE_6  6.000 0.000 999.000  999
#17                S.| MMSE_7  7.000 0.000 999.000  999
#18                S.| MMSE_8  8.000 0.000 999.000  999
#19             S.WITH      I  1.380 0.245   5.623    0
#20              Means      I 27.460 0.122 225.854    0
#21              Means      S -1.151 0.091 -12.614    0
#22         Intercepts MMSE_0  0.000 0.000 999.000  999
#23         Intercepts MMSE_1  0.000 0.000 999.000  999
#24         Intercepts MMSE_2  0.000 0.000 999.000  999
#25         Intercepts MMSE_3  0.000 0.000 999.000  999
#26         Intercepts MMSE_4  0.000 0.000 999.000  999
#27         Intercepts MMSE_5  0.000 0.000 999.000  999
#28         Intercepts MMSE_6  0.000 0.000 999.000  999
#29         Intercepts MMSE_7  0.000 0.000 999.000  999
#30         Intercepts MMSE_8  0.000 0.000 999.000  999
#31          Variances      I  3.932 0.481   8.175    0
#32          Variances      S  2.316 0.268   8.657    0
#33 Residual.Variances MMSE_0  5.446 0.216  25.261    0
#34 Residual.Variances MMSE_1  5.446 0.216  25.261    0
#35 Residual.Variances MMSE_2  5.446 0.216  25.261    0
#36 Residual.Variances MMSE_3  5.446 0.216  25.261    0
#37 Residual.Variances MMSE_4  5.446 0.216  25.261    0
#38 Residual.Variances MMSE_5  5.446 0.216  25.261    0
#39 Residual.Variances MMSE_6  5.446 0.216  25.261    0
#40 Residual.Variances MMSE_7  5.446 0.216  25.261    0
#41 Residual.Variances MMSE_8  5.446 0.216  25.261    0

#here I have scripted how to extract the right information from these elements and report on them

###PLOT
#selecting the data needed for plotting the population mean trajectory and turning into a data frame
ul.mp.avg.plot.data<- as.data.frame(ul.mp.l[["gh5"]][["means_and_variances_data"]])

#some data wrangling
#this is a long data set
ul.mp.avg.plot.data <-ul.mp.avg.plot.data %>%
  rename(MMSE=(1))%>% #changing the name of the outcome to be more informative
  mutate(funr = row_number()-1) #adding an indicator for time (to be plotted on the x-axis)

#defining the plot
ul.mp.avg.plot <- ggplot(data=ul.mp.avg.plot.data, aes(x=funr, y=MMSE))

ul.mp.avg.plot+geom_line()+
  geom_line(data=pq, aes(x=funr, y=MMSE, group=ID), alpha=.2, color="grey")

### PARAMETER ESTIMATES  
#selecting the fixed parameter estimates to build a table with parameter estimates for our model
ul.mp.paramest.data <- as.data.frame(ul.mp.l[["parameters"]][["unstandardized"]])

#some data wrangling steps
ul.mp.paramest.data <- ul.mp.paramest.data %>%
  filter(paramHeader %in% c("Means", "Variances", "Residual.Variances", "S.WITH") ) #selecting which parameter estimates to report

### FIT
#fit statistics
ul.mp.fit.data <- as.data.frame(ul.mp.l[["summaries"]])
#some data wrangling steps
ul.mp.fit.data <- ul.mp.fit.data %>%
  select(Parameters, LL, AIC, BIC) %>% #selecting which fit statistics to report 
  pivot_longer(cols = tidyselect::everything(),
               names_to = "paramHeader",
               values_to = "est")

### COMBINED PARAMETERS AND FIT:
#combining fit statistics and parameter estimate data
ul.mp.paramest.data <- bind_rows(ul.mp.fit.data, ul.mp.paramest.data)

# a first rough table holding the data for this model 
gt(ul.mp.paramest.data)

#######focusing on single models
#######latent growth curve modeling

## the steps are essentially the same as for LGCA: create, run, and report model

lcga.MMSE <- mplusObject(
  TITLE=          "MMSE-LCGA2 - unconditional quadratic model - sq@0-2 classes-residuals UNconstrained;",
  ANALYSIS = "type=mixture;
  starts=100 20;
  estimator=MLR;",
  VARIABLE= "CLASSES = c(2);
  USEVARIABLES ARE MMSE_0 MMSE_1 MMSE_2 MMSE_3 MMSE_4 MMSE_5 MMSE_6 MMSE_7 MMSE_8;",
  MODEL=          "%overall%
    i s q| 
    MMSE_0@0 MMSE_1@1 MMSE_2@2 
    MMSE_3@3 MMSE_4@4 MMSE_5@5
    MMSE_6@6 MMSE_7@7 MMSE_8@8;
    MMSE_0; 
    MMSE_1; 
    MMSE_2; 
    MMSE_3;
    MMSE_4;
    MMSE_5;
    MMSE_6;
    MMSE_7;
    MMSE_8;
  i;
  s@0;
  q@0;
%c#1%
  [i s q];
%c#2%
  [i s q];",
  PLOT=
    "Type = PLOT3;
    Series =  MMSE_0 (0) MMSE_1 (1) MMSE_2 (2) MMSE_3 (3) 
    MMSE_4 (4) MMSE_5 (5) MMSE_6 (6) MMSE_7 (7);",
  
  OUTPUT= "TECH11 TECH14;",
  rdata = pq.w
)

lcga.MMSE <- mplusModeler(lcga.MMSE, 
                      modelout="lcga.MMSE.inp", 
                      run=1L)

#reading the mplus output file using readModels: turns output file into list
lcga.MMSE.l <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lcga.MMSE.out",
                      recursive = F)

#these are the elements that you need to extract from the list created above:
##PLOT average growth per class: lcga.MMSE.l[["gh5"]][["means_and_variances_data"]][["y_estimated_means"]]
##PARAMETERS: lcga.MMSE.l[["parameters"]][["unstandardized"]]
##FIT: lcga.MMSE.l[["summaries"]]
##CLASS COUNTS: lcga.MMSE.l[["class_counts"]][["mostLikely"]]
##CLASS MEMBERSHIP for each individual participant: lcga.MMSE.l[["gh5"]][["individual_data"]][["raw_data"]]

## PLOT
#selecting the data needed for plotting the population mean trajectory for each class and turning into a data frame
lcga.MMSE.avg.plot.data<- as.data.frame(lcga.MMSE.l[["gh5"]][["means_and_variances_data"]][["y_estimated_means"]])

#some data wrangling
lcga.MMSE.avg.plot.data <-lcga.MMSE.avg.plot.data %>%
  rename(MMSE.1=(1))%>% #more informative name for outcome in class 1 
  rename(MMSE.2=(2))%>% #more informative names for outcome in class 2
  mutate(funr = row_number()-1) #adding variable for time 

#defining the plot
lcga.MMSE.avg.plot <- ggplot(data=lcga.MMSE.avg.plot.data, aes(x=funr, y=MMSE.1)) #to plot class 1 trajectory

lcga.MMSE.avg.plot+geom_line(color="blue")+
  geom_line(data=lcga.MMSE.avg.plot.data, aes(x=funr, y=MMSE.2), color="red")+ #to plot class 2 trajectory
  geom_line(data=pq, aes(x=funr, y=MMSE, group=ID), alpha=.2, color="grey") # to plot the individual observed trajectories, not yet colored by class membership, that is done later in the script

## PARAMETERS
#selecting the fixed parameter estimates to build a table with parameter estimates for our model
lcga.MMSE.paramest.data <- as.data.frame(lcga.MMSE.l[["parameters"]][["unstandardized"]])

#some data wrangling steps
lcga.MMSE.paramest.data <- lcga.MMSE.paramest.data %>% 
  filter(paramHeader %in% c("Means", "Variances", "Residual.Variances", "S.WITH") ) #selecting the parameter estimates to report, run part of this script to see what information this list holds before filtering, because you may want to choose different elements than I did. 

## FIT
#fit statistics
lcga.MMSE.fit.data <- as.data.frame(lcga.MMSE.l[["summaries"]])
#some data wrangling steps
lcga.MMSE.fit.data <- lcga.MMSE.fit.data %>%
  select(Parameters, LL, AIC, BIC, Entropy, 
         T11_VLMR_Mean, T11_VLMR_PValue, BLRT_KM1LL, BLRT_PValue) %>%
  pivot_longer(cols = tidyselect::everything(),
               names_to = "paramHeader",
               values_to = "est")

## CLASS COUNTS
#class counts and proportions
lcga.MMSE.class.count.data <- as.data.frame(lcga.MMSE.l[["class_counts"]][["mostLikely"]]) %>%
  pivot_longer(cols = c(count, proportion),
               names_to = "paramHeader",
               values_to = "est")

## COMBINED FIT, CLASS COUNT AND PARAMETER
#combining fit statistics, class count and parameter estimate data
lcga.MMSE.paramest.data <- bind_rows(lcga.MMSE.fit.data, lcga.MMSE.paramest.data, lcga.MMSE.class.count.data)

# a first rough table holding the data for this model using gt() for tidier tables.
gt(lcga.MMSE.paramest.data)

## INDIVIDUAL CLASS MEMBERSHIP DATA
# this data frame (see that I transposed it using t...) holds the individual data:
lcga.MMSE.classmemb<-as.data.frame(t(lcga.MMSE.l[["gh5"]][["individual_data"]][["raw_data"]]))
# it has the same number of observations as the original participant data and holds
## (in this case) V1:V9: the original time scores for each of the participants
## V10:V17: for each participant the ESTIMATED time score (not sure why only 8 rather than 9 time points)
## V18: most likely class membership, either class 1 or class 2 in this case.
# However, it lacks a participant identifier, but the rows are in the same order as in the .dat file that the analysis was run on, 
# so a simple cbind can be used to join this information with the original data set, but is a bit tricky.
# for more certainty, you may want to join using the combination of observed time scores (columns V1:V9) as as unique identifier to match.

# This is the most correct approach, because participants the same observed time scores also get the same estimates for the estimates reported in V10:V18.
# See here:
lcga.MMSE.classmemb1 <- lcga.MMSE.classmemb %>%
  group_by(select(., V1:V9)) %>%
  arrange(select(., V1:V9)) %>%
  mutate(ind.u = row_number()) %>%
  ungroup()

lcga.MMSE.classmemb.u <- lcga.MMSE.classmemb %>%
  group_by(select(., V1:V9)) %>%
  arrange(select(., V1:V9)) %>%
  distinct(select(., V1:V9)) %>%
  ungroup() %>%
  mutate(id.u.timescores = row_number())

lcga.MMSE.classmemb1 <- left_join(lcga.MMSE.classmemb1, lcga.MMSE.classmemb.u, by = c("V1", "V2", "V3", "V4","V5", "V6","V7", "V8","V9"))
lcga.MMSE.classmemb1

#a bit of data wrangling to join the class membership data with the original wide data set

k<-0:8 #defining a vector of integers to use in renaming the variables
lcga.MMSE.classmemb <- lcga.MMSE.classmemb %>%
  select(V1:V9, V18) %>%
  mutate(across(.cols= c(V1:V9), 
                .fns=~ (.=.),
                .names = "MMSE_{k}")
  ) %>% #renaming V1:V9 into MMSE_k
  select(V18, MMSE_0:MMSE_8) %>%
  mutate(across(everything(),
                ~(case_when(.== 999~NA_real_,
                            T~.)))) %>% #turning missing value 999 into NA
  rename(MMSE.2class.classnr = V18) %>%
  distinct()

pq.w <- left_join(pq.w, lcga.MMSE.classmemb, by = c("MMSE_0", "MMSE_1", "MMSE_2", "MMSE_3","MMSE_4", "MMSE_5", "MMSE_6","MMSE_7", "MMSE_8"))

freq(pq.w$MMSE.2class.classnr)

#to color the spaghetti plot by most likely class membership
pq <- left_join(pq, pq.w[c("ID", "MMSE.2class.classnr")], by = "ID") #merge class information with long data set

lcga.MMSE.avg.plot+geom_line(color=custom.col[1], linewidth = 1)+ #script uses the ggplot object defined before, but now with observed trajectories colored by class membership
  geom_line(data=lcga.MMSE.avg.plot.data, aes(x=funr, y=MMSE.2), color=custom.col[2], linewidth = 1)+
  geom_line(data=pq, aes(x=funr, y=MMSE, group=ID, color=as.factor(MMSE.2class.classnr)), alpha=0.2)+ 
  scale_color_manual(values = custom.col)

pq.ws <- pq.w %>%
  select(ID, MMSE.2class.classnr)

pqb <- left_join(pqb, pq.ws, by = "ID")

age <- stats::lm(formula = age ~ MMSE.2class.classnr,
          data=pqb)

summary(age)

#for GMM the same steps as for LCGA can be followed

####### Creating, running and reporting multiple models at once

#step: lGca
prepareMplusData(
  pq.w,
  filename = "Z:/mplusautomationtutorial/datfiles/pq_mmse_cesd.dat",
  keepCols=c("ID", "MMSE_0", "MMSE_1", "MMSE_2", "MMSE_3", "MMSE_4", "MMSE_5", "MMSE_6", "MMSE_7", "MMSE_8", 
             "CESD_0", "CESD_1", "CESD_2", "CESD_3", "CESD_4", "CESD_5", "CESD_6", "CESD_7", "CESD_8")
  )

createModels("Z:/mplusautomationtutorial/templatefiles/template_lGca_paquid.txt")

runModels("Z:/mplusautomationtutorial/in-outputfiles/lGca",
          recursive = T,
          replaceOutfile = "modifiedDate",
          showOutput = T) 

allOutputlGca <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lGca",
                     recursive = T)

lGca.fit <- as.data.frame(list.mapv(allOutputlGca, c(summaries))) %>%
  rownames_to_column() %>%
  rename(value = (2)) %>%
  mutate(modelid.ch = substr(rowname, 43, 113)) %>%
  mutate(output = substr(rowname, 119, 10000)) %>%
  group_by(modelid.ch) %>%
  ungroup()

lGca.modelid <- lGca.fit %>%
  distinct(modelid.ch) %>%
  mutate(modelid = row_number())

lGca.fit <- left_join(lGca.fit, lGca.modelid, by = "modelid.ch") %>%
  select(-rowname)

target <- c("Title", 
            "Parameters", "LL", "AIC", "BIC", "Entropy")
lGca.fit <- lGca.fit %>%
  filter(output %in% target)

lGca.fit <- lGca.fit %>%
  pivot_wider(
    names_from = output,
    values_from = value
  )

cols.num <- c("modelid", "Parameters", "LL", "AIC", "BIC")
lGca.fit[cols.num] <- sapply(lGca.fit[cols.num],as.numeric)

lGca.warn <- as.data.frame(list.mapv(allOutputlGca, c(warnings, errors))) %>%
  rownames_to_column() %>%
  mutate(value = NA_character_) %>%
  mutate(modelid = NA_real_) %>%
  add_row(modelid = 1:12) 
  
#there are no warnings or errors, currently. Use LCga.warn DW steps how to continue if lGca step causes warnings
#now I simply add a column warnings to lCga.fit that is empty to add that as information to the fit statistics table for lGca

lGca.warn.s <- lGca.warn %>%
  select(modelid, value) %>%
  rename(warning = value) %>%
  group_by(modelid)%>%
  summarise(warning = str_c(warning, collapse = ","))

lGca.fit <- left_join(lGca.fit, lGca.warn.s, by = "modelid")

lGca.fit <- lGca.fit %>%
  mutate(outcome = substr(modelid.ch, 6, 9)) %>%
  mutate(modeltype =  substr(modelid.ch, 11, 15)) %>%
  mutate(restype = substr(modelid.ch, 59, 71)) %>%
  mutate(deviance = -2*LL) %>% #to calculate the difference in deviance between successive (nested) models 
  group_by(outcome, restype) %>%
  arrange(outcome, restype, modeltype) %>%
  mutate(chisq = lag(deviance)-deviance) %>%
  mutate(dif.df = Parameters-lag(Parameters)) %>%
  mutate(pdeviance = round(pchisq(chisq, dif.df, lower.tail = F), 5)) %>%
  ungroup()

lGca.outcomeid <- lGca.fit %>%
  distinct(outcome) %>%
  mutate(outcomeid = row_number())

lGca.fit <- left_join(lGca.fit, lGca.outcomeid, by= "outcome")

fit.lGca <- lGca.fit %>%
  select(modelid, Title, Parameters,
         LL, AIC, BIC, warning, deviance, chisq, dif.df, pdeviance, outcome, restype) %>%
  group_by(outcome, restype) %>%
  gt() %>%
  fmt_number(columns = c(LL), decimals = 2)

fit.lGca

lGca.plot <- as.data.frame(list.mapv(allOutputlGca, c(gh5)))  %>%
  rownames_to_column() %>%
  rename(value = (2))%>%
  mutate(modelid.ch = substr(rowname, 43, 113)) %>%
  mutate(outcome = substr(modelid.ch, 6, 9)) %>%
  mutate(output = substr(rowname, 119, 10000)) %>%
  mutate(output1 = as.numeric(substr(output, 50, 51)))
lGca.plot <- left_join(lGca.plot, lGca.modelid, by = "modelid.ch") %>%
  select(-rowname)

lGca.plot <- lGca.plot %>%
  filter(str_detect(output, "means_and_variances_data.y_estimated_means.values")) %>%
  group_by(modelid) %>%
  arrange(output1) %>%
  mutate(order.by.modelid = row_number()) %>%
  ungroup() %>%
  mutate(classnr = case_when(order.by.modelid <= 9~1,#classnr is not relevant for lGca, but not problematic either
                             order.by.modelid > 9 & order.by.modelid <=18~2,
                             order.by.modelid > 18 & order.by.modelid <=27~3,
                             order.by.modelid > 27 & order.by.modelid <=36~4)) %>%
  group_by(modelid, classnr) %>%
  arrange(modelid, order.by.modelid, classnr) %>%
  mutate(funr = row_number()-1) %>%
  ungroup()

lGca.plot <-lGca.plot %>%
  arrange(modelid)

lGca.plot <- lGca.plot %>%
  select(-output, -output1, -order.by.modelid)%>%
  mutate(value = as.numeric(value))

#add a selection of variables from lGca.fit to lGca.plot data frame
lGca.fit.s <- lGca.fit %>%
  select(modelid, Title, outcomeid)

lGca.plot <- left_join(lGca.plot, lGca.fit.s, by="modelid")

# Make list of y-axis scale limits to loop over
lGca.model_y_axislimits <- list(
  cesd= list(c(0,60)), # 8 y axis limits for cesd
  mmse=list(c(0,30))) %>% # 8 y axis limits for mmse
  list_flatten() # flatten this list to hold 16 limits (for 16 models in total)
# Make plots.
# Earlier I made a dataframe holding the information on the outcomes modelled. Useful here, because the outcome also defines what is on y-axis
plotslGca = list()
for (i in 1:length(lGca.outcomeid)) {
  p = ggplot(subset(lGca.plot, outcomeid %in% i), aes(x=funr, y=value, group=interaction(as.factor(modelid), classnr), colour = as.factor(classnr))) +
    geom_line()+
    scale_y_continuous(lGca.outcomeid$outcome[[i]], limits=c(lGca.model_y_axislimits[[i]][1], lGca.model_y_axislimits[[i]][2]))+
    labs(linewidth = "proportion per class (0-1)",
         color = "class")+
    facet_wrap(~modelid, ncol=6, nrow = 1)
  plotslGca[[i]] = p
}
ggarrange(plotlist= plotslGca, ncol=1, nrow=2, common.legend = T, legend = "right")

#library(plyr)
#justSummaries <- do.call("rbind.fill",
#                         sapply(allOutput,"[", "summaries"))

#HTMLSummaryTable(
#  allOutput,
#  filename = "Z:/mplusautomationtutorial/in-outputfiles/lGca/MyModelSummary.html",
#  display = TRUE,
#  keepCols = c("Title", "LL", "AIC", "BIC", "AICC", "Parameters"),
#  sortBy = "AIC")

#allModelParameters <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lGca",
#                                 what="parameters",
#                                 recursive = T)

#parallelModels <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lGca")
#names(parallelModels)

#compareModels(parallelModels[["cesd.lGca1...unconditional.linear.model.out"]],
#              parallelModels[["cesd.lGca2...unconditional.quadratic.model.out"]],
#              show = c("diff", "pdiff", "summaries", "unique"),
#              equalityMargin = c(param = .05, pvalue = .02),
#              sort = "type", diffTest = TRUE, showNS = FALSE)

#step: lCga
prepareMplusData(
  pq.w,
  filename = "Z:/mplusautomationtutorial/datfiles/pq_mmse_cesd.dat",
  keepCols=c("ID", "MMSE_0", "MMSE_1", "MMSE_2", "MMSE_3", "MMSE_4", "MMSE_5", "MMSE_6", "MMSE_7", "MMSE_8", 
             "CESD_0", "CESD_1", "CESD_2", "CESD_3", "CESD_4", "CESD_5", "CESD_6", "CESD_7", "CESD_8")
)

createModels("Z:/mplusautomationtutorial/templatefiles/template_lCga_paquid.txt")

runModels("Z:/mplusautomationtutorial/in-outputfiles/lCga",
          recursive = T,
          replaceOutfile = "modifiedDate",
          showOutput = T) 

allOutputlCga <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lCga",
                        recursive = T)

#get_results("Z:/mplusautomationtutorial/in-outputfiles/lCga")

#library(plyr)
#justSummarieslCga <- do.call("rbind.fill",
#                         sapply(allOutputlCga,"[", "summaries"))

#HTMLSummaryTable(
#  allOutputlCga,
#  filename = "Z:/mplusautomationtutorial/in-outputfiles/lGca/MyModelSummary.html",
#  display = TRUE,
#  keepCols = c("Title", "Parameters", "LL", "AIC", "BIC", "Entropy", "T11_KM1LL", "T11_VLMR_Mean", "T11_VLMR_PValue", "T11_LMR_Value", "T11_LMR_PValue", "BLRT_KM1LL", "BLRT_PValue"),
#  sortBy = "Title")


#my own operationalization of about the same as HTMLSummaryTable provides, but with more control over what is reported
lCga.fit <- as.data.frame(list.mapv(allOutputlCga, c(summaries, class_counts))) %>%
  rownames_to_column() %>%
  rename(value = (2)) %>%
  mutate(modelid.ch = substr(rowname, 43, 131)) 

lCga.fit <- lCga.fit %>%
  mutate(outcome = substr(modelid.ch, 6, 9)) %>%
  mutate(nclass = substr(modelid.ch, 57, 57)) %>%
  mutate(output = substr(rowname, 137, 10000)) %>%
  group_by(modelid.ch) %>%
  ungroup()

lCga.modelid <- lCga.fit %>%
  distinct(modelid.ch) %>%
  mutate(modelid = row_number())

lCga.fit <- left_join(lCga.fit, lCga.modelid, by = "modelid.ch") %>%
  select(-rowname) %>%
  mutate(output = case_when(output == "mostLikely.class"~"mostLikely.class1",
                            output == "mostLikely.count"~"mostLikely.count1",
                            output == "mostLikely.proportion"~"mostLikely.proportion1",
                            T~output))

lCga.outcomeid <- lCga.fit %>%
  distinct(outcome) %>%
  mutate(outcomeid = row_number())

lCga.fit <- left_join(lCga.fit, lCga.outcomeid, by= "outcome")

target <- c("Title", 
            "Parameters", "LL", "AIC", "BIC", "Entropy", 
            "T11_KM1LL", "T11_VLMR_Mean", "T11_VLMR_PValue", "T11_LMR_Value", "T11_LMR_PValue", "BLRT_KM1LL", "BLRT_PValue", 
            "mostLikely.class1", "mostLikely.count1", "mostLikely.proportion1",
            "mostLikely.class2", "mostLikely.count2", "mostLikely.proportion2",
            "mostLikely.class3", "mostLikely.count3", "mostLikely.proportion3",
            "mostLikely.class4", "mostLikely.count4", "mostLikely.proportion4")
lCga.fit <- lCga.fit %>%
  filter(output %in% target)

lCga.fit <- lCga.fit %>%
  pivot_wider(
    names_from = output,
    values_from = value
  )

lCga.fit <- lCga.fit %>% 
  unite(class_counts, c("mostLikely.count1", "mostLikely.count2", "mostLikely.count3", "mostLikely.count4"), na.rm=T, sep= " ", remove = F) %>%
  unite(class_probs, c("mostLikely.proportion1", "mostLikely.proportion2", "mostLikely.proportion3", "mostLikely.proportion4"), na.rm=T, sep= " ", remove = F) %>%
  select(-modelid.ch, -mostLikely.class1, -mostLikely.class2, -mostLikely.class3, -mostLikely.class4)

lCga.fit <- lCga.fit %>%
  mutate(nclass = as.numeric(nclass))

cols.num <- c("modelid", "nclass", "Parameters", "LL", "AIC", "BIC", "Entropy", 
              "T11_KM1LL", "T11_VLMR_Mean", "T11_VLMR_PValue", "T11_LMR_Value", "T11_LMR_PValue", "BLRT_KM1LL", "BLRT_PValue")
lCga.fit[cols.num] <- sapply(lCga.fit[cols.num],as.numeric)

fit <- lCga.fit %>%
  gt() %>%
  fmt_number(columns = c(LL), decimals = 2) 
fit

lCga.warn <- as.data.frame(list.mapv(allOutputlCga, c(warnings, errors))) %>%
  rownames_to_column() %>%
  rename(value = (2)) %>%
  mutate(modelid.ch = substr(rowname, 43, 131)) %>%
  mutate(nclass = substr(modelid.ch, 57, 57)) %>%
  mutate(output = substr(rowname, 136, 10000)) %>% #here not very informative, because only a number of the data element of this list
  group_by(modelid.ch) %>%
  ungroup()

lCga.warn <- left_join(lCga.warn, lCga.modelid, by = "modelid.ch")%>%
  select(-rowname)

#to identify unique warnings
lCga.warn.unique <- lCga.warn %>%
  distinct(value) 

#to indicate which are relevant warnings, based on your own assessment of unique warnings. use rel.warn to define these
lCga.warn.unique <- lCga.warn.unique %>%
  mutate(relevant.warning = case_when(str_detect(value, rel.warn)~1,
                                                     T~0))

lCga.warn <- left_join(lCga.warn, lCga.warn.unique, by = "value")

lCga.warn <- lCga.warn %>%
  filter(relevant.warning == 1)
lCga.warn.s <- lCga.warn %>%
  select(modelid, value) %>%
  rename(warning = value) %>%
  group_by(modelid)%>%
  summarise(warning = str_c(warning, collapse = ","))

lCga.fit <- left_join(lCga.fit, lCga.warn.s, by = "modelid")

fit.lCga <- lCga.fit %>%
  select(modelid, Title, nclass, Parameters, class_counts, class_probs,
         LL, AIC, BIC, Entropy, 
         T11_VLMR_Mean, T11_VLMR_PValue, BLRT_KM1LL, BLRT_PValue, warning, outcome) %>%
  group_by(outcome)%>%
  gt() %>%
  fmt_number(columns = c(LL), decimals = 2) %>%
  tab_style(
    style = list(
      cell_text(size = "xx-small")),
    locations= cells_body(
      columns = warning
      )
    ) %>%
  tab_style(
    style = list(
      cell_text(color = "red")),
    locations= cells_body(
      columns = warning,
      rows = !is.na(warning)
      )
    ) %>%
  tab_style(
        style = list(
          cell_text(size = "x-small")),
        locations= cells_body(
          columns = Title) 
        )%>%
  tab_style(
    style = list(
      cell_fill(
        color = "red", alpha = 0.1)),
      locations= cells_body(
        rows = nclass == 1)
    )%>%
  tab_style(
    style = list(
      cell_fill(
        color = "orange", alpha = 0.1)),
    locations= cells_body(
      rows = nclass == 2)
  )%>%
  tab_style(
    style = list(
      cell_fill(
        color = "green", alpha = 0.1)),
    locations= cells_body(
      rows = nclass == 3)
  )%>%
  tab_style(
    style = list(
      cell_fill(
        color = "blue", alpha = 0.1),
      cell_borders(sides= "bottom",
                   weight = 2)),
    locations= cells_body(
      rows = nclass == 4)
  )

fit.lCga

lCga.plot <- as.data.frame(list.mapv(allOutputlCga, c(gh5)))  %>%
  rownames_to_column() %>%
  rename(value = (2))%>%
  mutate(modelid.ch = substr(rowname, 43, 131)) %>%
  mutate(outcome = substr(modelid.ch, 6, 9)) %>%
  mutate(nclass = substr(modelid.ch, 57, 57)) %>%
  mutate(output = substr(rowname, 137, 10000)) %>%
  mutate(output1 = as.numeric(substr(output, 50, 51)))

lCga.plot <- left_join(lCga.plot, lCga.modelid, by = "modelid.ch") %>%
  select(-rowname)

lCga.plot <- lCga.plot %>%
  filter(str_detect(output, "means_and_variances_data.y_estimated_means.values")) %>%
  group_by(modelid) %>%
  arrange(output1) %>%
  mutate(order.by.modelid = row_number()) %>%
  ungroup() %>%
  mutate(classnr = case_when(order.by.modelid <= 9~1,
                             order.by.modelid > 9 & order.by.modelid <=18~2,
                             order.by.modelid > 18 & order.by.modelid <=27~3,
                             order.by.modelid > 27 & order.by.modelid <=36~4)) %>%
  group_by(modelid, classnr) %>%
  arrange(modelid, order.by.modelid, classnr) %>%
  mutate(funr = row_number()-1) %>%
  ungroup()

lCga.plot <-lCga.plot %>%
  arrange(modelid)

lCga.plot <- lCga.plot %>%
  select(-output, -output1, -order.by.modelid)%>%
  mutate(value = as.numeric(value))

lCga.classprob <- lCga.fit %>%
  select(modelid, Title, outcomeid, starts_with("mostLikely.proportion")) %>%
  pivot_longer(cols = starts_with("mostLikely.proportion"),
               names_to = "classnr.ch",
               values_to = "mostLikely.proportion") %>%
  mutate(classnr = as.numeric(gsub("\\D", "", classnr.ch))) %>%
  select(-classnr.ch)%>%
  mutate(mostLikely.proportion = as.numeric(mostLikely.proportion))
         
lCga.plot <- left_join(lCga.plot, lCga.classprob, by=c("modelid", "classnr"))

plot.lCga <- ggplot(data=lCga.plot, aes(x=funr, y=value, group=interaction(as.factor(modelid), classnr), colour = as.factor(classnr), linewidth=mostLikely.proportion))
plot.lCga+
  geom_line()+
  facet_wrap(~modelid.ch)+ 
  scale_y_continuous(limits=c(0, 60)) # a bit dirty/quick way of removing negative values for y, but works for now
# https://stackoverflow.com/questions/11585954/varying-axis-labels-formatter-per-facet-in-ggplot-r I am more or less misusing facet_wrap here in the first place, grid.arrange or something would be better

# a more advanced way of plotting

#define elements to be used in the for loop below to loop over
#the amount of elements should be as much as i's in the loop, which is 16 in this case. ie. 16 models need to be plotted.
# Make list of models to loop over.
#lCga.modelids <- as.vector(unique(lCga.plot$modelid))
# Make list of outcome names to label y-axis to loop over.
#lCga.modeloutcomes <- as.vector(lCga.fit$outcome)
# Make list of y-axis scale limits to loop over
#lCga.model_y_axislimits <- list(
#  cesd= list(c(0,60),c(0,60),c(0,60),c(0,60),c(0,60),c(0,60),c(0,60),c(0,60)), # 8 y axis limits for cesd
#  mmse=list(c(0,30),c(0,30),c(0,30),c(0,30),c(0,30),c(0,30),c(0,30),c(0,30)) # 8 y axis limits for mmse
#) %>%
#  list_flatten() # flatten this list to hold 16 limits (for 16 models in total)
# Make plots.
#plotslCga = list()
#for (i in 1:length(lCga.modelids)) {
#  p = ggplot(subset(lCga.plot, modelid %in% i), aes(x=funr, y=value, group=interaction(as.factor(modelid), classnr), colour = as.factor(classnr), linewidth=mostLikely.proportion)) +
#    geom_line()+
#    scale_y_continuous(lCga.modeloutcomes[[i]], limits=c(lCga.model_y_axislimits[[i]][1], lCga.model_y_axislimits[[i]][2]))+
#    labs(linewidth = "proportion per class (0-1)")
#  plotslCga[[i]] = p
#}
#ggarrange(plotlist= plotslCga, ncol=4, nrow=4, common.legend = T, legend = "right")


#grid.arrange(grobs = plotslCga) 
#using gridExtra package

# Another option: create pdf where each page is a separate plot.
#pdf("plotslCga.pdf")
#for (i in 1:length(lCga.modelids)) {
#  print(plotslCga[[i]])
#}
#dev.off()

# Make list of y-axis scale limits to loop over
lCga.model_y_axislimits <- list(
  cesd= list(c(0,60)), # 8 y axis limits for cesd
  mmse=list(c(0,30))) %>% # 8 y axis limits for mmse
  list_flatten() # flatten this list to hold 16 limits (for 16 models in total)
# Make plots.
# Earlier I made a dataframe holding the information on the outcomes modelled. Useful here, because the outcome also defines what is on y-axis
plotslCga = list()
for (i in 1:length(lCga.outcomeid)) {
  p = ggplot(subset(lCga.plot, outcomeid %in% i), aes(x=funr, y=value, group=interaction(as.factor(modelid), classnr), colour = as.factor(classnr), linewidth=mostLikely.proportion)) +
    geom_line()+
    scale_y_continuous(lCga.outcomeid$outcome[[i]], limits=c(lCga.model_y_axislimits[[i]][1], lCga.model_y_axislimits[[i]][2]))+
    labs(linewidth = "proportion per class (0-1)",
         color = "class")+
    facet_wrap(~modelid, ncol=4, nrow = 2)
  plotslCga[[i]] = p
}
ggarrange(plotlist= plotslCga, ncol=1, nrow=2, common.legend = T, legend = "right")

################### GMM

prepareMplusData(
  pq.w,
  filename = "Z:/mplusautomationtutorial/datfiles/pq_mmse_cesd.dat",
  keepCols=c("ID", "MMSE_0", "MMSE_1", "MMSE_2", "MMSE_3", "MMSE_4", "MMSE_5", "MMSE_6", "MMSE_7", "MMSE_8", 
             "CESD_0", "CESD_1", "CESD_2", "CESD_3", "CESD_4", "CESD_5", "CESD_6", "CESD_7", "CESD_8")
)

createModels("Z:/mplusautomationtutorial/templatefiles/template_GMM__paquid.txt")

runModels("Z:/mplusautomationtutorial/in-outputfiles/gmm_",
          recursive = T,
          replaceOutfile = "modifiedDate",
          showOutput = T) 

allOutputgmm_ <- readModels("Z:/mplusautomationtutorial/in-outputfiles/gmm_",
                            recursive = T)

gmm_.fit <- as.data.frame(list.mapv(allOutputgmm_, c(summaries, class_counts))) %>%
  rownames_to_column() %>%
  rename(value = (2)) %>%
  mutate(modelid.ch = substr(rowname, 43, 131)) %>%
  mutate(outcome = substr(modelid.ch, 6, 9)) %>%
  mutate(randomstatement = substr(modelid.ch, 51, 54)) %>%
  mutate(nclass = substr(modelid.ch, 57, 57)) %>%
  mutate(output = substr(rowname, 137, 10000)) %>%
  group_by(modelid.ch) %>%
  ungroup()

gmm_.outcomeid <- gmm_.fit %>%
  distinct(outcome) %>%
  mutate(outcomeid = case_when(outcome == "cesd"~1,
                               outcome == "mmse"~2
  )) %>%
  arrange(outcomeid)

gmm_.fit <- left_join(gmm_.fit, gmm_.outcomeid, by= "outcome")

gmm_.modelid <- gmm_.fit %>%
  distinct(modelid.ch) %>%
  mutate(modelid = row_number())

gmm_.modelid <- gmm_.fit %>%
  select(modelid.ch, outcomeid, nclass, randomstatement)%>%
  arrange (outcomeid, nclass, modelid.ch) %>%
  distinct(modelid.ch, .keep_all=T) %>%
  group_by(outcomeid) %>%
  mutate(order.modelsbyoutcome = case_when(randomstatement == "iq.0"~1,
                                           randomstatement == "sq.0"~2,
                                           randomstatement == ".q.0"~3,
                                           randomstatement == "c..i"~4,
                                           randomstatement == "c.is"~5))%>%
  arrange (outcomeid, nclass, order.modelsbyoutcome, modelid.ch) %>%
  ungroup()%>%
  mutate(modelid = row_number()) %>%
  select(modelid, modelid.ch, order.modelsbyoutcome)

gmm_.fit <- left_join(gmm_.fit, gmm_.modelid, by = "modelid.ch") %>%
  select(-rowname) %>%
  mutate(output = case_when(output == "mostLikely.class"~"mostLikely.class1",
                            output == "mostLikely.count"~"mostLikely.count1",
                            output == "mostLikely.proportion"~"mostLikely.proportion1",
                            T~output))

target <- c("Title", 
            "Parameters", "LL", "LLCorrectionFactor", "AIC", "BIC", "Entropy", 
            "T11_KM1LL", "T11_VLMR_Mean", "T11_VLMR_PValue", "T11_LMR_Value", "T11_LMR_PValue", "BLRT_KM1LL", "BLRT_PValue", 
            "mostLikely.class1", "mostLikely.count1", "mostLikely.proportion1",
            "mostLikely.class2", "mostLikely.count2", "mostLikely.proportion2",
            "mostLikely.class3", "mostLikely.count3", "mostLikely.proportion3",
            "mostLikely.class4", "mostLikely.count4", "mostLikely.proportion4")
gmm_.fit <- gmm_.fit %>%
  filter(output %in% target)

gmm_.fit <- gmm_.fit %>%
  pivot_wider(
    names_from = output,
    values_from = value
  )

gmm_.fit <- gmm_.fit %>% 
  unite(class_counts, c("mostLikely.count1", "mostLikely.count2"), na.rm=T, sep= " ", remove = F) %>%
  unite(class_probs, c("mostLikely.proportion1", "mostLikely.proportion2"), na.rm=T, sep= " ", remove = F) %>%
  select(-modelid.ch, -mostLikely.class1, -mostLikely.class2)

gmm_.fit <- gmm_.fit %>%
  mutate(nclass = as.numeric(nclass))

cols.num <- c("modelid", "nclass", "Parameters", "LL", "LLCorrectionFactor", "AIC", "BIC", "Entropy", 
              "T11_KM1LL", "T11_VLMR_Mean", "T11_VLMR_PValue", "T11_LMR_Value", "T11_LMR_PValue", "BLRT_KM1LL", "BLRT_PValue")
gmm_.fit[cols.num] <- sapply(gmm_.fit[cols.num],as.numeric)

fit.gmm_ <- gmm_.fit %>%
  gt() %>%
  fmt_number(columns = c(LL), decimals = 2) 
fit.gmm_

gmm_.warn <- as.data.frame(list.mapv(allOutputgmm_, c(warnings, errors))) %>%
  rownames_to_column() %>%
  rename(value = (2)) %>%
  mutate(modelid.ch = substr(rowname, 43, 131)) %>%
  mutate(nclass = substr(modelid.ch, 57, 57)) %>%
  mutate(output = substr(rowname, 136, 10000)) %>% #here not very informative, because only a number of the data element of this list
  group_by(modelid.ch) %>%
  ungroup()

gmm_.warn <- left_join(gmm_.warn, gmm_.modelid, by = "modelid.ch")%>%
  select(-rowname)

#to identify unique warnings
#add warnings that are relevant for report to rel.warn
gmm_.warn.unique <- gmm_.warn %>%
  distinct(value) 
#to indicate which are relevant warnings, based on your own assessment of unique warnings. use rel.warn to define these
gmm_.warn.unique <- gmm_.warn.unique %>%
  mutate(relevant.warning = case_when(str_detect(value, rel.warn)~1,
                                      T~0))

gmm_.warn <- left_join(gmm_.warn, gmm_.warn.unique, by = "value")
gmm_.warn <- gmm_.warn %>%
  filter(relevant.warning == 1)

gmm_.warn.s <- gmm_.warn %>%
  select(modelid, value) %>%
  rename(warning = value) %>%
  group_by(modelid)%>%
  summarise(warning = str_c(warning, collapse = ","))

gmm_.fit <- left_join(gmm_.fit, gmm_.warn.s, by = "modelid")

fit.gmm_ <- gmm_.fit %>%
  select(modelid, Title, nclass, Parameters, class_counts, class_probs,
         LL, LLCorrectionFactor, AIC, BIC, Entropy, 
         T11_VLMR_Mean, T11_VLMR_PValue, BLRT_KM1LL, BLRT_PValue, warning, outcome, outcomeid) %>%
  group_by(outcome, nclass)%>%
  arrange(outcomeid, nclass, modelid) %>%
  mutate(LL.cd = 
           (((lag(Parameters)*lag(LLCorrectionFactor))-(Parameters*LLCorrectionFactor))/
              (lag(Parameters)-Parameters)))%>%
  mutate(LL.TRd = ((-2*(lag(LL)-LL))/
                     LL.cd)) %>%
  mutate(LL.TRd = case_when(LL.TRd <0~-LL.TRd,
                            T~LL.TRd)) %>% #if the calculation of LL.TRd results in negative value, change the sign.
  mutate(dif.df.LL.TRd = Parameters-lag(Parameters)) %>%
  mutate(p.LL.TRd = round(pchisq(LL.TRd, dif.df.LL.TRd, lower.tail = F), 5)) %>%
  ungroup() %>%
  group_by(outcome)%>%
  gt() %>%
  fmt_number(columns = c(LL), decimals = 2) %>%
  tab_style(
    style = list(
      cell_text(size = "xx-small")),
    locations= cells_body(
      columns = warning
    )
  )%>%
  tab_style(
    style = list(
      cell_text(color = "red")),
    locations= cells_body(
      columns = warning,
      rows = !is.na(warning)
    )
  ) %>%
  tab_style(
    style = list(
      cell_text(size = "x-small")),
    locations= cells_body(
      columns = Title) 
  )%>%
  tab_style(
    style = list(
      cell_fill(
        color = "red", alpha = 0.1)),
    locations= cells_body(
      rows = nclass == 1)
  ) %>%
  tab_style(
    style = list(
      cell_fill(
        color = "orange", alpha = 0.1)),
    locations= cells_body(
      rows = nclass == 2)
  )%>%
  tab_style(
    style = list(
      cell_fill(
        color = "green", alpha = 0.1)),
    locations= cells_body(
      rows = nclass == 3)
  )%>%
  tab_style(
    style = list(
      cell_fill(
        color = "blue", alpha = 0.1),
      cell_borders(sides= "bottom",
                   weight = 2)),
    locations= cells_body(
      rows = nclass == 4)
  )%>%
  cols_hide(c(modelid, outcomeid, LLCorrectionFactor, T11_VLMR_Mean, T11_VLMR_PValue, BLRT_KM1LL, BLRT_PValue, LL.cd))

fit.gmm_

setwd("Z:/mplusautomationtutorial/results")
fit.gmm_|>
  gtsave("fit.gmm_pq.docx")
setwd("Z:/mplusautomationtutorial/")
gmm_.plot <- as.data.frame(list.mapv(allOutputgmm_, c(gh5)))  %>%
  rownames_to_column() %>%
  rename(value = (2))%>%
  mutate(modelid.ch = substr(rowname, 43, 131)) %>%
  mutate(outcome = substr(modelid.ch, 6, 9)) %>%
  mutate(nclass = substr(modelid.ch, 57, 57)) %>%
  mutate(output = substr(rowname, 137, 10000)) %>%
  mutate(output1 = as.numeric(substr(output, 50, 51)))

gmm_.plot <- left_join(gmm_.plot, gmm_.modelid, by = "modelid.ch") %>%
  select(-rowname)

gmm_.plot <- gmm_.plot %>%
  filter(str_detect(output, "means_and_variances_data.y_estimated_means.values")) %>%
  group_by(modelid) %>%
  arrange(output1) %>%
  mutate(order.by.modelid = row_number()) %>%
  ungroup() %>%
  mutate(classnr = case_when(order.by.modelid <= 9~1,
                             order.by.modelid > 9 & order.by.modelid <=18~2,
                             order.by.modelid > 18 & order.by.modelid <=27~3,
                             order.by.modelid > 27 & order.by.modelid <=36~4)) %>%
  group_by(modelid, classnr) %>%
  arrange(modelid, order.by.modelid, classnr) %>%
  mutate(funr = row_number()-1) %>%
  ungroup()

gmm_.plot <-gmm_.plot %>%
  arrange(modelid)

gmm_.plot <- gmm_.plot %>%
  select(-output, -output1, -order.by.modelid)%>%
  mutate(value = as.numeric(value))

gmm_.classprob <- gmm_.fit %>%
  select(modelid, Title, outcomeid, starts_with("mostLikely.proportion")) %>%
  pivot_longer(cols = starts_with("mostLikely.proportion"),
               names_to = "classnr.ch",
               values_to = "mostLikely.proportion") %>%
  mutate(classnr = as.numeric(gsub("\\D", "", classnr.ch))) %>%
  select(-classnr.ch)%>%
  mutate(mostLikely.proportion = as.numeric(mostLikely.proportion))

gmm_.plot <- left_join(gmm_.plot, gmm_.classprob, by=c("modelid", "classnr"))

plot.gmm_ <- ggplot(data=gmm_.plot, aes(x=funr, y=value, group=interaction(as.factor(modelid), classnr), colour = as.factor(classnr), linewidth=mostLikely.proportion))
plot.gmm_+
  geom_line()+
  facet_wrap(~modelid.ch)+ 
  scale_y_continuous(limits=c(0, 60)) # a bit dirty/quick way of removing negative values for y, but works for now
# https://stackoverflow.com/questions/11585954/varying-axis-labels-formatter-per-facet-in-ggplot-r I am more or less misusing facet_wrap here in the first place, grid.arrange or something would be better
# then I can also more easily change the name of the plots.
# https://stackoverflow.com/questions/71621139/change-axis-label-automatically-when-y-parameter-is-changed-in-ggplot2


# Make list of y-axis scale limits to loop over
gmm_.model_y_axislimits <- list(
  cesd= list(c(0,60)), # 8 y axis limits for cesd
  mmse=list(c(0,30))) %>% # 8 y axis limits for mmse
  list_flatten() # flatten this list to hold 2 limits (for 2 outcomes in total)
# Make plots.
plotsgmm_ = list()
for (i in 1:length(gmm_.outcomeid)) {
  p = ggplot(subset(gmm_.plot, outcomeid %in% i), aes(x=funr, y=value, group=interaction(as.factor(modelid), classnr), colour = as.factor(classnr), linewidth=mostLikely.proportion)) +
    geom_line()+
    scale_y_continuous(gmm_.outcomeid$outcome[[i]], limits=c(gmm_.model_y_axislimits[[i]][1], gmm_.model_y_axislimits[[i]][2]))+
    labs(linewidth = "proportion per class (0-1)",
         color = "class")+
    facet_wrap(~modelid, ncol=5, nrow = 1)
  plotsgmm_[[i]] = p
}
ggarrange(plotlist= plotsgmm_, ncol=1, nrow=2, common.legend = T, legend = "right")

#not needed now

#allModelWarningslCga <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lCga",
#                                     what="warn_err",
#                                     recursive = T)

#allModelParameterslCga <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lCga",
#                                 what="parameters",
#                                 recursive = T)
#allModelClasscountlCga <- readModels("Z:/mplusautomationtutorial/in-outputfiles/lCga",
#                                     what="class_counts",
#                                     recursive = T)