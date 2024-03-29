# grofit analysis, updated from my older code, to analyze 96 well optical curves.
# template code



## first load packages and grofit ####

# grofit is no longer updated on CRAN so must be downloaded manually from here `[https://cran.r-project.org/src/contrib/Archive/grofit/]`
# easiest way  is to download from link above and *use rtools* to install with line below
#install.packages("grofit_1.1.1-1.tar.gz",repos=NULL, type="source")

library(tidyverse)
library(data.table) #used for grofit prep
library(grofit)
library(RColorBrewer)# for graph colors
library(rstatix) # pipeline friendly stats analysis

rm(list = ls())


# set meta information 
data = "data.csv" # name of data file
label = "label.csv" #name of label file
plate_run = "species_date" #name of plate run. generally species and date of plate run


### read in data ####
d <- read.csv(file=paste("input/", data, sep=""))
labels <- read.csv(file=paste("input/", label, sep=""))%>%
  mutate(treatment = as.factor(treatment),
         well = as.factor(well),
         microbe = as.factor(microbe),
         species=as.factor(species))

d <- rename(d, c("time" = "Time")) #because i dont like capitals


# make time in minutes (rounded to nearest 15 (1 read is taken every 15 minutes)
d$time <- (seq.int(nrow(d)))*15
 

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


gr<- gr %>%
  relocate( c(treatment, microbe, species), .before = time)

write_csv(gr, file= paste("formatted_OD_values/", plate_run, "_formatted.csv", sep=""))



### grofit curve fitting ####

t<- gr%>%
  select(!c("microbe","treatment","time","species"))

# now make a df of only times
times<-t%>%
  colnames()%>%
  as.matrix()%>%
  t()%>%
  matrix(byrow=T, ncol=ncol(t), nrow=nrow(t))%>%
  as.data.frame()




# create a number string of ODinit; replicate into a matrix the size of your raw data

ODinit<-gr[,5]
ODmat<-matrix(ODinit, byrow=FALSE, ncol=ncol(gr)-4, nrow=nrow(gr))



# subtract starting OD
tOD2 <- gr[,5:ncol(gr)]-ODmat



#then any negatives set to zero
tOD2 <- replace(tOD2, tOD2 < 0, 0)

grow.m2<-cbind(gr[,1:3],tOD2)


#then also look at without subtracting starting OD
od_raw <- gr%>%
  select(!time)%>%
  melt(id.vars = c("treatment","microbe","species"), value.name = "od", variable.name = "time")



#graphing curves to visually check things
grow.long <- grow.m2%>%
  melt(id.vars = c("treatment","microbe","species"), value.name = "od", variable.name = "time")

curves<- ggplot(grow.long)+
  geom_point(aes(x=time, y=od, color=treatment))+
  facet_grid(~microbe)
curves




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
growth.test<-gcFit(times,grow.m2, control=control1)




# post curve verification -------------------------------------------------


#use the built in summary function and write parameters  a df
grofit<-summary.gcFit(growth.test)

# lambda = length of lag phase
# mu = maximum growth rate
# A = carrying capacity


#then remove unreliable wells
grofit<-grofit[grofit$reliability==TRUE,]

# #and set NAs to zero for mu and A
grofit$mu.model[is.na(grofit$mu.model)] <- 0
grofit$A.model[is.na(grofit$A.model)] <- 0


grofit <- rename(grofit, c("microbe" = "AddId"))
grofit <- rename(grofit, c("treatment" = "TestId"))
grofit <- rename(grofit, c("species" = "concentration"))
grofit$treatment <- as.factor(grofit$treatment)
grofit$microbe <- as.factor(grofit$microbe)
grofit$species <- as.factor(grofit$species)



#then if I didnt use all wells, drop empty ones
grofit<- grofit[!(grofit$treatment==""),]

grofit <- droplevels(grofit)

# then write to csv for use later in analysis
write_csv(grofit, file = paste("output/", plate_run, "_grofitresults.csv", sep=""))




# start here for analysis ####
grofit <- read_csv(file = paste("output/", plate_run, "_grofitresults.csv", sep=""))



### quick speed graphing and tablesS ####

# also adjust levels so low med high in order
# grofit$treatment<- factor(grofit$treatment, levels = c("control", "low", "medium", "high"))



alpha<- grofit %>%
  ggplot(aes(x=treatment, y=A.model))+
  geom_boxplot(aes(fill=treatment), alpha=.5)+
  geom_jitter(aes(fill=treatment), shape=21, color="black")+
  ylab("Max growth")+
  xlab("Treatment")+
  facet_grid(~microbe)+
  scale_fill_brewer(palette = "Set2")
alpha

ggsave(plot=alpha, file = paste("output/graphs/", plate_run, "_alpha.jpg", sep=""),
       width = 6,
       height = 3)




grofit %>%
  group_by(microbe)%>%
  kruskal_test(A.model ~ treatment)

grofit %>%
  group_by(microbe)%>%
  dunn_test(A.model ~ treatment) %>%
  print(n=40)








