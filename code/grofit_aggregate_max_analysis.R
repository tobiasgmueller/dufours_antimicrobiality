# grofit analysis, updated from my older code, to analyze 96 well optical curves.
# HERE IS AN ATTEMPT TO GET DATA FROM WELLS THAT CANT HAVE CURVES FITTED
# BY SIMPLE AVERAGING BY THE HOUR AND FIGURING OUT THE HIGHEST OD VALUE
# may 2022
# Tobias Mueller


# skip data read in and curve fitting if not running on a new plate
# these have already been done and saved as a csv.




## first packages and grofit ####

# grofit is no longer updated on CRAN so must be downloaded manually here `[https://cran.r-project.org/src/contrib/Archive/grofit/]`
#easiest way  is to download from link above and use rtools to install with line below
#install.packages("C:/Users/obiew/Desktop/github/dufours_antimicrobiality/grofit_1.1.1-1.tar.gz", repos=NULL, type="source")

library(tidyverse)
library(data.table) #used for grofit prep
library(grofit)
library(FSA) # for dunntest
library(RColorBrewer)# for graph colors


### read in data ####
rm(list = ls())


# run 1, with colletes inequalis
#d <- read.csv("input/dufours_may6_2022_tidy.csv")
#labels <- read.csv("input/6may2022_dufours_plate_setup.csv")


# run 2, to test if frozen colletes glands have same results
# and to verify if scc477 and dc3000 can even grow in lb broth
d <- read.csv("input/dufours_frozen_13may2022_test.csv")
labels <- read.csv("input/dufours_frozen_13may2022_plate_setup.csv")




d <- rename(d, c("time" = "Time")) #because i dont like capitals





# fix labels classes
labels$treatment <- as.factor(labels$treatment)
labels$well <- as.factor(labels$well)
labels$microbe <- as.factor(labels$microbe)

# make time in minutes (rounded to nearest 15 (1 read is taken every 15 minutes)
d$time <- (seq.int(nrow(d)))*15




### grofit bypass aggregate section #### _______________________

# First create hour grouping variable which will be aggregated across
# first drop the last data point so its a multiple of 4


d$hour <- rep(1:(nrow(d)/4), each = 4)


gr %>% 
  
  
Avgs <- mydf %>% group_by(Group) %>% summarize(Avg = mean(results))
Avgs





#now to get it to match the required df format of grofit

# transpose df
d <- t(d)
d <- as.data.frame(d)


# get the row labels as column 1
setDT(d, keep.rownames = TRUE)[]
# make time the headers
d <- as.data.frame(d)

#fix df formate
names(d) <- d[1,]
d <- d[-1,]


# merge raw data from plate reader with treatment labels
gr <- merge(d, labels, by.x="time", by.y="well", incomparables=NA)


#now add 2 random columns because the package demands it
gr <- cbind(gr, species="C. inequalis")


gr<- gr %>%
  relocate( c(treatment, microbe, species), .before = time)










write.csv(gr,"formatted_OD_values/aggregate_test", row.names = FALSE)

### grofit curve fitting ####


t<-read.csv("formatted_OD_values/aggregate_test", header=F)

times<-t[1,5:ncol(t)]
times<-matrix(times, byrow=T, ncol=ncol(gr)-4, nrow=nrow(gr))
# this removes labels and makes a matrix of just time stamps
times<-data.frame(times)
# then dataframes it


# create a number string of ODinit; replicate into a matrix the size of your raw data
ODinit<-gr[,5]
ODmat<-matrix(ODinit, byrow=FALSE, ncol=ncol(gr)-4, nrow=nrow(gr))


# subtract starting OD
tOD2 <- gr[,5:ncol(gr)]-ODmat



#then any negatives set to zero
tOD2 <- replace(tOD2, tOD2 < 0, 0)

grow.m2<-cbind(gr[,1:3],tOD2)





# set controls for how the gro.fit modeling runs
control1<-grofit.control(fit.opt="b", log.y.gc=FALSE, interactive=T)

# neg.nan.act       -- Logical, indicates wether the program should stop when negative growth values or NA values appear (TRUE). Otherwise the program removes this values silently (FALSE). Improper values may be caused by incorrect data or input errors. Default: FALSE.
# clean.bootstrap   -- Logical, determines if negative values which occur during bootstrap should be removed (TRUE) or kept (FALSE). Note: Infinite values are always removed. Default: TRUE.
# suppress.messages -- Logical, indicates wether grofit messages (information about current growth curve, EC50 values etc.) should be displayed (FALSE) or not (TRUE). This option is meant to speed up the processing of high throuput data. Note: warnings are still displayed. Default: FALSE.
# fit.opt           -- Indicates wether the program should perform a model fit ("m"), a spline fit ("s") or both ("b"). Default: "b".
# log.x.gc          -- Logical, indicates wether a ln(x+ 1) should be applied to the time data of the growth curves. Default: FALSE.
# log.y.gc          -- Logical, indicates wether a  ln(y+ 1)should be applied to the growth data of the growth curves. Default: FALSE.
# interactive       -- Logical, controls whether the fit of each growth curve is controlled manually by the user. Default: TRUE.
# nboot.gc          -- Number of bootstrap samples used for the model free growth curve fitting. Use nboot.gc=0 to disable the bootstrap. Default: 0.
# smooth.gc         -- Parameter describing the smoothness of the spline fit; usually (not necessary) in (0;1]. Set smooth.gc=NULL causes the program to query an optimal value via cross validation techniques. Note: This is partly experimental. In future improved implementations of the smooth.spline function may lead to different results. See documentation of the R function smooth.spline for further details. Especially for datasets with few data points the option NULL might result in a too small smoothing parameter, which produces an error in smooth.spline. In that case the usage of a fixed value is recommended. Default: NULL.
# model.type        -- Character vector giving the names of the parametric models which should be fitted to the data. Default: c("gompertz", "logistic", "gompertz.exp", "richards").
# have.atleast      -- Minimum number of different values for the response parameter one shoud have for estimating a dose response curve. Note: All fit procedures require at least six unique values. Default: 6.
# parameter         -- The column of the output table which should be used for creating a dose response curve. See documentation of gcFit, drFit or summary.gcFit for further details. Default: 9, which represents the maximum slope of the parametric growth curve fit.
# smooth.dr         -- Smoothing parameter used in the spline fit by smooth.spline during dose response curve estimation. Usually (not necessesary) in (0; 1]. See documentation of smooth.spline for further details. Default: NULL.
# log.x.dr          -- Logical, indicates wether a ln( x+1) should be applied to the concentration data of the dose response curves. Default: FALSE.
# log.y.dr          -- Logical, indicates wether a ln( y+ 1) should be applied to the response data of the dose response curves. Default: FALSE.
# nboot.dr          -- Numeric value, defining the number of bootstrap samples for EC50 estimation. Use nboot.dr=0 to disable bootstrapping. Default: 0.



#run grofit
growth.test<-gcFit(od,grow.m2, control=control1)



#use the built in summary function and write parameters  a df
grofit<-summary.gcFit(growth.test)

# lambda = length of lag phase
# mu = maximum growth rate
# A = carrying capacity





#then remove unreliable wells
grofit<-grofit[grofit$reliability==TRUE,]

# #and set NAs to zero for nu and A
grofit$mu.model[is.na(grofit$mu.model)] <- 0
grofit$A.model[is.na(grofit$A.model)] <- 0


grofit <- rename(grofit, c("microbe" = "AddId"))
grofit <- rename(grofit, c("treatment" = "TestId"))
grofit <- rename(grofit, c("species" = "concentration"))
grofit$treatment <- as.factor(grofit$treatment)
grofit$microbe <- as.factor(grofit$microbe)
grofit$species <- as.factor(grofit$species)



#then I totally forgot, I only used 95 wells so drop the empty one...
grofit<- grofit[!(grofit$treatment==""),]

grofit <- droplevels(grofit)

# then write to csv for use later in analysis
write.csv(grofit, "output/dufours_6may2022_grofitresults.csv")





# start here for analysis ####

grofit<- read.csv("output/dufours_6may2022_grofitresults.csv")




### quick speed graphing ####


# hmmm... so strangely the controls did not grow well for DC3000 nand SCC477 but other treatments did...
grofit %>%
  group_by(microbe)%>%
  filter(treatment == "control") %>%
  summarise(max=max(A.model))


alpha<- grofit %>%
  ggplot(aes(x=treatment, y=A.model))+
  geom_boxplot(aes(fill=treatment), alpha=.5)+
  geom_jitter(aes(fill=treatment), shape=21, color="black")+
  ylab("Max growth")+
  xlab("Treatment")+
  facet_wrap(~microbe)+
  scale_fill_brewer(palette = "Set2")
alpha

  ggsave(plot=alpha,"output/graphs/grofit_alpha.pdf")

grofit %>%
  ggplot(aes(x=treatment, y=mu.model))+
  geom_boxplot()+
  facet_wrap(~microbe)



  
  
# lets also look at things without the positive control
grofit_no_strep <- grofit %>%
  filter(!treatment == " +")



kwtest.A<- lapply(split(grofit, grofit$microbe), function(i){
  kruskal.test(A.model ~ treatment, data = i)
})
kwtest.A


#when you dunn all together only strep comes out significant because well such a crazy effect
dunn.A<- lapply(split(grofit, grofit$microbe), function(i){
  dunnTest(A.model~treatment, data=i, method="holm")
})

dunn.A

# if you remove strep you get some significant treatment impacts
dunn.A<- lapply(split(grofit_no_strep, grofit_no_strep$microbe), function(i){
  dunnTest(A.model~treatment, data=i, method="holm")
})

dunn.A








