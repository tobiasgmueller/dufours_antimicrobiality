# 2021 plate analysis
# 4mM h2o2
# 3/12/21
# tobias mueller



#------------------- setup and packages ---------------------------------------

#setwd("C:/Users/obiew/Desktop/vannette/growthassay_2021")

library(grofit)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(data.table)
library(growthcurver)

rm(list = ls())


#---------------------- part I --  reformatting the df --------------------------------------------------

#read in data 
d <- read.csv("input/4mMh2o2_5april2021_export.csv")
labels <- read.csv("input/2021labels.csv")
d <- rename(d, c("time" = "Time"))

# fix labels classes
labels$treatment <- as.factor(labels$treatment)
labels$well <- as.factor(labels$well)
labels$microbe <- as.factor(labels$microbe)


# make time in minutes (rounded to nearest 15 (1 read is taken every 15 minutes)
d$time <- (seq.int(nrow(d)))*15


#now to get it to match the required df format

# transpose df
d <- t(d)
d <- as.data.frame(d)

# get the row labels as column 1
setDT(d, keep.rownames = TRUE)[]
# make time the headers

# redataframe it
d <- as.data.frame(d)


#fix df formate
names(d) <- d[1,]
d <- d[-1,]


# merge raw data from plate reader with treatment labels
gr <- merge(d, labels, by.x="time", by.y="well", incomparables=NA)


#now add 2 random columns because the package demands it
gr <- cbind(gr, treatment="4mM h2o2", date="mar 12")


# then for my sanity we're gonna move the treatment labels to the front of the df
# not that it matters but I like to look at it from time to time

#useing function moveMe, move treatment information to the front
moveMe <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}

gr <- moveMe(gr, c("treatment", "microbe", "date"), "first")

# further analysis gets sad with strain and type columns so drop those here before writing
# theyre mostly for graphing and can be added back later
gr<- subset(gr, select=-c(name,type))


write.csv(gr, "input/4mM_h2o2_formatted.csv", row.names = FALSE)

# okay we now have our data frame in the correct format to move forward
# to the next step
# which is more dataframe reformatting!!




#-------------- Part II -- quick graphing ------ # cut this part --------------------------------------
# grill that cheese. melt it.
melt1<-reshape2::melt(d, id.vars = "time", variable.name = "well")

# fix column names
colnames(melt1)[1] <- "well"
colnames(melt1)[2] <- "time"

#make a df for easy graphing
graph <- merge(melt1, labels, by="well", incomparables=NA)
graph$microbe <- as.factor(graph$microbe)
graph$treatment <- as.factor(graph$treatment)
graph$well <- as.factor(graph$well)


graph <- arrange(graph, time, group_by = well)


# adjust for starting OD value
graph<- graph %>%
  group_by(well) %>% 
  mutate(value_adjust=value- first(value))



ggplot(graph, aes(x=time, y=value, color=treatment)) +
  geom_point(aes(), alpha=0.5) +
  facet_wrap(graph$well, scales = "free_y") +
  theme_bw() +
  labs(title="4mM h2o2")


ggplot(graph, aes(x=time, y=value_adjust, color=treatment)) +
  geom_point(aes(), alpha=0.5) +
  facet_wrap(graph$microbe ~ graph$name + graph$type, scales = "free_y") +
  theme_bw() +
  labs(title="4mM h2o2")





# --------------------- begin the real part ------------------------------------------------------

# check the data frame is the same - it should be from left to right:
# 3 treatment information columns at the start - there MUST be 3
# then a row labeled "time" that contains all the well numbers, each well number being a seperate row
# to the right of this the column headers are now the time stamps
# the cell information in each row is the OD reading for the well at that time stamp




gr<-read.csv("input/4mM_h2o2_formatted.csv", na.strings=c(".", "NA"))
t<-read.csv("input/4mM_h2o2_formatted.csv", header=F)
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

grow.m2<-cbind(gr[,1],gr[,2:3],tOD2)





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


#decline - well 25


#use the built in summary function and write parameters  a df
parms<-summary.gcFit(growth.test)





# lambda = length of lag phase
# mu = maximum growth rate
# A = carrying capacity



#then remove unreliable wells
parms<-parms[parms$reliability==TRUE,]

# #and set NAs to zero for nu and A
parms$mu.model[is.na(parms$mu.model)] <- 0
parms$A.model[is.na(parms$A.model)] <- 0


parms <- rename(parms, c("microbe" = "AddId"))
parms <- rename(parms, c("treatment" = "TestId"))
parms$treatment <- as.factor(parms$treatment)
parms$microbe <- as.factor(parms$microbe)

# then write to csv for use later in analysis
write.csv(parms, "output/4mM_h2o2_summary.csv")





























