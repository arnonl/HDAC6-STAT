# Sub-version of "/Users/Arnon/Post_at_NYU/Projects/Util_codes/plots_errorbars_gragh_and_gee.R"
thisFileName <- "relativeWeightGain_LME4_model.R"

# Plot mice weight phenotype in the HDAC experiment with correction for Day Of Life
# (problematic since we dont have a SE in each time point - need to check with regression models)
# 

library(ggplot2)
library(reshape)


##### Set parameters #####

relativeGain <- F
tickEvery <- 20

genderTypes <- c("Males", "Females")
names(genderTypes) <- c("M","F")

# Loads data into R from excel
# 
uDir <- "/Users/Arnon/Post_at_NYU/Projects/Util_codes/"
dataDir <- "/Users/Arnon/Post_at_NYU/Projects/HDAC_study/03_processedData"

data <- read.csv(file = file.path(dataDir, "all_mice_statVsControl.csv"), header = T)
#
longData <- melt.data.frame(
    data, id.vars = c("mouse", "gender", "treatment", "cage", "genotype"))
colnames(longData) <- c(colnames(longData)[1:5], "week", "weight")
head(longData)

# Add numeric time
# patterns within parenthesis are captured and assigned to \\1, \\2, etc. 
# Print out only \\2 which has the numeric part
longData$numTime <- as.numeric(sub("(w)(.)", "\\2", as.character(longData$week), perl = TRUE))

# write.csv(file = file.path(dataDir,"all_mice_longFormat.csv"), longData)
# longData <- read.csv(file = file.path(dataDir,"all_mice_longFormat.csv"), header = T)

# See how many cases we have for each time point
table(longData$week)

# Divide here every category with the first observation of time=0 and save it under wgRatio
# First, let's see how many observation we have for each week
if (var(table(longData$week)) == 0) {
    nEl <- max(table(longData$week))
} else {
    stop("The number of elements are not equal across groups")
}


if (!relativeGain){
    # Here I set parameters for the original weight
    wValues <- "weight"
    breaks=seq(0,40,tickEvery)
    yLabel <- "Weight [gr]"
    measurevar <- "weight"
} else {
    # Here I set parameters for the normalized data (relative weight gain)
    longData$wgRatio <- longData$weight/longData$weight[1:nEl]
    wValues <- "wgRatio"
    breaks=0:4
    yLabel <- "Relative weight gain"
    measurevar <- "wgRatio"
}

# This function summarizes for each "groupvars" the following values -> mean; sd; se; ci
source(file.path(uDir, "summarySE.R"))
# Time is factorial
longDataSum <- summarySE(longData, measurevar = measurevar, groupvars = c("week", "treatment", "genotype"))
# Time is numeric
longDataSum <- summarySE(longData, measurevar = measurevar, groupvars = c("numTime", "treatment", "genotype"))
longDataSum <- transform(longDataSum, comCat=factor(paste0(treatment,"-", genotype))) # comCat = combined categories
# Defines groups to skip one segment between before HFD and after HFD. If not needed, switch aes to group=treatment
longDataSum$grp <- factor(c(rep(1:2, 8),rep(3:4, 5)))

##### Publishable plot #####
pd <- position_dodge(0.1)
yAxisData <- longDataSum[,wValues]
ggplot(data = longDataSum, aes(x=numTime, y=yAxisData, colour=treatment, group=grp)) + 
    # If you want to include all 4 cases change to -> colour=comCat, group=comCat
    geom_line(position=pd, size = 0.5) +
    geom_errorbar(aes(ymin=yAxisData-se, ymax=yAxisData+se), colour="black", width=.3, position=pd, size=0.25) +
    geom_point(position=pd, size=.5, shape=21, fill="white") + # 21 is filled circle
    xlab("Time (weeks)") +
    ylab(yLabel) + # growth
    scale_colour_hue(name="Type",    # Legend label, use darker colors
                     # The values of break mast be the same as the 'colour' and 'group'
                     # If you use colour=comCat, group=comCat than breaks=c("LDP-wt", "Control-wt")
                     breaks=c("LDP", "Control"), 
                     labels=c("LDP", "Control"),
                     l=40) +                    # Use darker colors, lightness=40
    ggtitle("") +
    expand_limits(y=0) +                        # Expand y range
    # limits = range(breaks) is eqvivalent to ylim()
    scale_y_continuous(breaks = breaks,
                       labels = as.character(breaks), limits = range(breaks)) +
    # X-axis scale (scale_x_continuous) when using numeric values
    scale_x_continuous(breaks = c(seq(0,12,2)),
                     labels = c("0","2","4","6","8","10","12")) +
    # X-axis scale (scale_x_discrete) when using factorial values
    # limits=c("w5","w7","w9","w11","w13","w15","w17")
    # scale_x_discrete(breaks = c("w0","w2","w4","w6","w8","w10","w12"),
    #                  labels = c("0","2","4","6","8","10","12")) + 
    geom_vline(xintercept = 7.5) +
    theme_bw() +
    theme(axis.text.x = element_text(face="plain", color="black", 
                                         size=9, angle=0), # X-axis font size
          axis.text.y = element_text(face="plain", color="black", 
                                     size=9, angle=0)) + # Y-axis font size
    theme(legend.justification=c(1,0),
          legend.position=c(1,0), 
          legend.key.size = unit(.3, "cm"), # Legend size (each entry)
          legend.text=element_text(size=8),
          legend.title=element_text(size=8, face="bold"))

# All saving options:
# "wt_males_contVsLDP"
# "ko_males_contVsLDP"
# "wt_females_contVsLDP"
# "hz_females_contVsLDP"
ggsave(filename = file.path(dataDir, 
                            paste(currGenotype, genderTypes[currGender],
                                   str_c(levels(longDataSum$treatment), collapse = "Vs"),"pdf", sep=".")), 
                            width = 5, height = 4) # width = 5, height = 4

#-----End of plotting and saving-----#



#----------- START MODELING ------------#
#########################################
head(longData)
str(longData)
longData$cage <- factor(longData$cage)
# Add an indication for before and after switch to HFD
longData$onHFD <- factor(c(rep("Before", sum(longData$numTime<8)), rep("After", sum(longData$numTime>=8))))
# Creates a dummy variable for treatment (0-control, 1-STAT)
longData$dummyTr <- as.numeric(longData$treatment=="LDP")
# Defines factor that divide the curve into three weight gain phases
longData$phase <- with(longData,factor(c(
    rep(1, sum(numTime<3)), 
    rep(2, sum(numTime>=3 & numTime<8)),
    rep(3, sum(numTime>=8))))
    )
# Defines the piecewise time variable (as suggested here: , not working)
longData$time1 <- longData$numTime
longData[longData$numTime >= 8, "time1"] <- 8

longData$time2 <- longData$numTime - 8
longData[longData$numTime < 8, "time2"] <- 0

longData$time2[longData$numTime >= 8, "time2"] <- longData$numTime

# Subset to before HFD
longDataBeforeHFD <- subset(longData, onHFD=="Before")
# Subset to HFD only
longDataOnHFD <- subset(longData, longData$numTime >= 8)
# Check what id you should place here + if the factors should checked by + or *
# 
# Assumptions:
#     
# - Data is matched across a similar characteristic and is non-independent.
# - The dependent variables should be assumed to follow a multivariate normal distribution.
# - It should be assumed that the effects in should be fixed effects.
# - The homogeneity of the variance is assumed, this can be tested with the help of Levene's test.


# Gives you a list that specifies for each data point how the coefficient would change if you remove this point 
dfbeta()


# Random and fixed effects:
# treatment is fixed effect
# HFD is a fixed effect
# mice Id is a random effect
# The intercept is a random effect

packageVersion("lme4")
library(lme4)
citation("lme4")
# Base on this: https://stats.stackexchange.com/questions/58745/using-lmer-for-repeated-measures-linear-mixed-effect-model
# Nice example of Plotting Random Effects of Mixed Models is available here:
# https://cran.r-project.org/web/packages/sjPlot/vignettes/sjplmer.html 


# Three levels model
###########################
# This test check for the fix effect of time + treatment. 
#                   + the random (slope) effect of each mouse by genotype,   
conditional.growth.byGenotype <- lmer(weight ~ numTime * dummyTr + 
                                      (numTime | genotype:mouse) + (numTime * dummyTr || genotype), 
                                  data=longDataBeforeHFD)

# Ploynomal model
######################
head(longDataBeforeHFD)
tail(longDataBeforeHFD)
colnames(longDataBeforeHFD)

#== Test 1 ==#
# The very basic model with no treatment, random intercept and slope
unconditional.growth.mod <- lmer(weight ~ (numTime + I(numTime^2)) + 
                                     ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
# The very basic model with no treatment, random slope only
unconditional.growth.ranSlope <- lmer(weight ~ (numTime + I(numTime^2)) + 
                                     (0 + (numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
# Comparing unconditional model with random slope and intercept with unconditional model with only random slope
anova(update(unconditional.growth.mod, REML=FALSE), update(unconditional.growth.ranSlope, REML=FALSE))

# Conclusion - the random intercept term has a significant effect on the model 


#== Test 2 ==#
# 'conditional' takes also gender + treatment + genotype in account
conditional.growth.mod <- lmer(weight ~ (numTime + I(numTime^2)) + gender + 
                                   ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)

# Comparing conditional with unconditional
sum.cond <- summary(conditional.growth.mod)
summary(unconditional.growth.mod)
anova(update(conditional.growth.mod, REML=FALSE), update(unconditional.growth.mod, REML=FALSE))

# Extract information
#######################
# Show the coefficients for each sample
coef(conditional.growth.mod)
coef(summary(conditional.growth.mod))
# Show the random effect for each mouse
ranef(conditional.growth.mod)$mouse

timeP <- unique(longDataBeforeHFD$numTime)
weight <- function(time){
    8.16241+time*5.772943-0.4727829 * time^2}
plot(timeP, weight(timeP))
ylim(0,26)

names(sum.cond)      # S4 class need to look at the slots
class(sum.cond)
methods(class="merMod") # See the methods for the merMod class

VarCorr(conditional.growth.mod)   # 
deviance(conditional.growth.mod) # 
plot(conditional.growth.mod)

# Model Diagnostics
y.hat   <- fitted(conditional.growth.mod)          # Fitted values
int.hat <- ranef(conditional.growth.mod)[[1]][[1]] # Predicted intercepts
res.hat <- residuals(conditional.growth.mod)       # Estimated residuals

qqnorm(int.hat, main="Random Intercepts"); qqline(int.hat)
qqnorm(res.hat, main="Residuals"); qqline(res.hat)
plot(y.hat, res.hat, xlab="Fitted Values", ylab="Residuals")
abline(h=0, lty=2)

# Conclution chi.sq(1)=22.69, p<0.0001 *** shows that males has a significantly higher weight gain rate (3.5 gr/week) over females  



#==Test 3 ==#
conditional.growth.mod.X <- lmer(weight ~ (numTime + I(numTime^2)) * gender + 
                                     ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)

# Comparing conditional with interaction with conditional with -no- interaction
anova(update(conditional.growth.mod, REML=FALSE), update(conditional.growth.mod.X, REML=FALSE))

# Conclution: there is an important interaction term between time and gender, Chi-sq(2)=72.6, p<0.0001



#==Test 4 ==#
conditional_1L.growth.mod <- lmer(weight ~ (numTime + I(numTime^2)) + gender + 
                                     ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
conditional_2L.growth.mod <- lmer(weight ~ (numTime + I(numTime^2)) + gender + treatment +
                                        ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
summary(conditional_2L.growth.mod)
anova(update(conditional_1L.growth.mod, REML=FALSE), update(conditional_2L.growth.mod, REML=FALSE))

# Conclution: there is a moderate significant effect for treatment, Chi-sq(1)=5.25, p=0.022, as LDP-treatment contributes 0.73 gr/week over control.



#==Test 4.5 ==#
conditional_1L.growth.mod.X <- lmer(weight ~ (numTime + I(numTime^2)) * gender + 
                                        ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
conditional_2L.growth.mod.X <- lmer(weight ~ (numTime + I(numTime^2)) * gender * treatment +
                                        ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
summary(conditional_2L.growth.mod.X)
anova(update(conditional_1L.growth.mod.X, REML=FALSE), update(conditional_2L.growth.mod.X, REML=FALSE))


coef(conditional.growth.mod)



#==Test 5 ==#
conditional_2L.growth.mod.X <- lmer(weight ~ (numTime + I(numTime^2)) + gender + treatment +
                                        ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
conditional_3L.growth.mod.X <- lmer(weight ~ (numTime + I(numTime^2)) + gender + treatment + genotype +
                                        ((numTime + I(numTime^2)) | mouse), data = longDataBeforeHFD)
summary(conditional_3L.growth.mod.X)
anova(update(conditional_2L.growth.mod.X, REML=FALSE), update(conditional_3L.growth.mod.X, REML=FALSE))

# Conclution: there is an insignificant effect for genotype on NC, Chi-sq(2)=0.13, p=0.94



#==Test 6 ==#
# This test use the full data.frame, contains both NC and high-fat diet
conditional_1L.growth.mod.2ph <- lmer(weight ~ ((time1 + I(time1^2)) + time2) + gender + 
                                        ((time1 + I(time1^2)) + time2 | mouse), data = longData)
conditional_2L.growth.mod.2ph <- lmer(weight ~ ((time1 + I(time1^2)) + time2) + gender + treatment +
                                        ((time1 + I(time1^2)) + time2 | mouse), data = longData)
summary(conditional_2L.growth.mod.X)
anova(update(conditional_2L.growth.mod.2ph, REML=FALSE), update(conditional_3L.growth.mod.2ph, REML=FALSE))



#==Test 7 ==#
# This test use both NC and high-fat diet, for males only
longDataMonly <- subset(longData, gender=="M")
conditional.growth.mod.1L.2ph.males <- lmer(weight ~ ((time1 + I(time1^2)) + time2) + treatment + 
                                          ((time1 + I(time1^2)) + time2 | mouse), data = longDataMonly)
conditional.growth.mod.2L.2ph.males <- lmer(weight ~ ((time1 + I(time1^2)) + time2) + treatment + genotype +
                                          ((time1 + I(time1^2)) + time2 | mouse), data = longDataMonly)
summary(conditional.growth.mod.2L.2ph.males)
anova(update(conditional.growth.mod.1L.2ph.males, REML=FALSE), 
      update(conditional.growth.mod.2L.2ph.males, REML=FALSE))




#==Test 8 ==#
# This test use both NC and high-fat diet, for females only
longDataFonly <- subset(longData, gender=="F")
conditional.growth.mod.1L.2ph.females <- lmer(weight ~ ((time1 + I(time1^2)) + time2) * treatment + 
                                                ((time1 + I(time1^2)) + time2 | mouse), data = longDataFonly)
conditional.growth.mod.2L.2ph.females <- lmer(weight ~ ((time1 + I(time1^2)) + time2) * treatment * genotype +
                                                ((time1 + I(time1^2)) + time2 | mouse), data = longDataFonly)
summary(conditional.growth.mod.2L.2ph.females)
anova(update(conditional.growth.mod.1L.2ph.females, REML=FALSE), 
      update(conditional.growth.mod.2L.2ph.females, REML=FALSE))



#== Test 9 ==#
# This test check the effect of genotype for males on high-fat diet
longDataOnHFD_Male_only <- subset(longDataOnHFD, gender=="M")
conditional.growth.mod.1L.2ndPH.males <- lmer(weight ~ numTime + treatment + (numTime | mouse), 
                                                data = longDataOnHFD_Male_only)
conditional.growth.mod.2L.2ndPH.males <- lmer(weight ~ numTime + treatment + genotype + (numTime | mouse), 
                                                data = longDataOnHFD_Male_only)
anova(update(conditional.growth.mod.1L.2ndPH.males, REML=FALSE), 
      update(conditional.growth.mod.2L.2ndPH.males, REML=FALSE))
summary(conditional.growth.mod.2L.2ndPH.males)
# Conclution: there is a significant effect of genotype for males on HFD, Chi-sq(1)=8.88, p=0.003

library(dplyr)

#== Test 9.1 ==#
# This test check the effect of treatment for WT males on high-fat diet
longDataOnHFD_wtMale_only <- longDataOnHFD %>% 
    filter(gender=="M" & genotype=="wt")
conditional.growth.mod.0L.2ndPH.males <- lmer(weight ~ numTime + (numTime | mouse), 
                                              data = longDataOnHFD_wtMale_only)
conditional.growth.mod.1L.2ndPH.males <- lmer(weight ~ numTime + treatment + (numTime | mouse), 
                                              data = longDataOnHFD_wtMale_only)
anova(update(conditional.growth.mod.0L.2ndPH.males, REML=FALSE), 
      update(conditional.growth.mod.1L.2ndPH.males, REML=FALSE))
summary(conditional.growth.mod.1L.2ndPH.males)
# Added at Sep 05,
# When comparing control and STAT for WT male mice,
# Treatment has a significant effect, Chi-sq(1)=13, p=0.00031 ***

#== Test 9.2 ==#
# This test check the effect of time-treatment interaction for males on high-fat diet
longDataOnHFD_Male_only <- subset(longDataOnHFD, gender=="M")
conditional.growth.mod.1L.HFD.males <- lmer(weight ~ numTime + treatment + (numTime | mouse), 
                                              data = longDataOnHFD_Male_only)
conditional.growth.mod.1L.HFD.males.X <- lmer(weight ~ numTime * treatment + (numTime | mouse), 
                                              data = longDataOnHFD_Male_only)
anova(update(conditional.growth.mod.1L.HFD.males, REML=FALSE), 
      update(conditional.growth.mod.1L.HFD.males.X, REML=FALSE))
summary(conditional.growth.mod.1L.HFD.males.X)
# Conclution: there is a significant effect of genotype for males on HFD, Chi-sq(1)=9.72, p=0.002 **



#== Test 9.5 ==#
# This test check the interaction term for males on high-fat diet
longDataOnHFD_Male_only <- subset(longDataOnHFD, gender=="M")
conditional.growth.mod.2L.HFD.males <- lmer(weight ~ numTime + treatment + genotype + (numTime | mouse), 
                                              data = longDataOnHFD_Male_only)
conditional.growth.mod.2L.HFD.males.X <- lmer(weight ~ numTime * treatment * genotype + (numTime | mouse), 
                                              data = longDataOnHFD_Male_only)
anova(update(conditional.growth.mod.2L.HFD.males, REML=FALSE), 
      update(conditional.growth.mod.2L.HFD.males.X, REML=FALSE))
summary(conditional.growth.mod.2L.HFD.males.X)
# Conclution: there is a significant interaction between time and treatment on HFD, Chi-sq(4)=36.6, p<0.0001 ***



#== Test 10 ==#
# This test check the effect of genotype for female mice on high-fat diet
longDataOnHFD_Fem_only <- subset(longDataOnHFD, gender=="F")
conditional.growth.mod.1L.HFD.females <- lmer(weight ~ numTime + treatment + (numTime | mouse), 
                                            data = longDataOnHFD_Fem_only)
conditional.growth.mod.2L.HFD.females <- lmer(weight ~ numTime + treatment + genotype + (numTime | mouse), 
                                                data = longDataOnHFD_Fem_only)
anova(update(conditional.growth.mod.1L.HFD.females, REML=FALSE), 
      update(conditional.growth.mod.2L.HFD.females, REML=FALSE))
# Conclution: there is an insignificant effect of genotype for female mice on HFD, Chi-sq(1) =0.71, p=0.4



#== Test 10.5 ==#
# This test check the effect of interaction term for female mice on high-fat diet
longDataOnHFD_Fem_only <- subset(longDataOnHFD, gender=="F")
conditional.growth.mod.1L.HFD.females <- lmer(weight ~ numTime + treatment + (numTime | mouse), 
                                              data = longDataOnHFD_Fem_only)
conditional.growth.mod.1L.HFD.females.X <- lmer(weight ~ numTime * treatment + (numTime | mouse), 
                                              data = longDataOnHFD_Fem_only)
anova(update(conditional.growth.mod.1L.HFD.females, REML=FALSE), 
      update(conditional.growth.mod.1L.HFD.females.X, REML=FALSE))
# Conclution: there is an insignificant effect of time-treatment interaction for female mice on HFD, Chi-sq(1) =0.09, p=0.75

str(terms(conditional.growth.mod))


# Libraries to help plot diffrences between models
library(sjPlot)
library(sjmisc)
sjp.lmer(fit, y.offset = .5)

plot(fitted(lmm), residuals(lmm), xlab = "Fitted Values", ylab = "Residuals")

# For likelihhood test use the REML=FALSE in your model


# use 'fitted' or 'predict' to plot 



# Another approch where the developing rate is varie
plgModel <- lmer(Score ~ 1 + Time1 + Time2 + (1 + Time1 + Time2 | ID),
                 data = myDataNew)
summary(plgModel)