library(survival)
library(KMsurv)

tongue$type[tongue$type == 1] <- "Aneuploid"
tongue$type[tongue$type == 2] <- "Diploid"


fit <- coxph(Surv(time = time, event = delta) ~ type,  
             data = tongue) 
plot(survfit(fit, newdata = data.frame(type = "Aneuploid")))
lines(survfit(fit, newdata = data.frame(type = "Diploid")))


fit2 <- coxph(Surv(time = time, event = delta) ~ as.factor(type),
              data = kidney)
plot(survfit(fit2))

https://stat.ethz.ch/R-manual/R-devel/library/survival/html/survfit.coxph.html
https://stat.ethz.ch/R-manual/R-devel/library/survival/html/cox.zph.html
http://stats.stackexchange.com/questions/61131/test-cox-proportional-hazard-assumption-bad-schoenfeld-residuals
http://statistics.ats.ucla.edu/stat/examples/asa/test_proportionality.htm
http://rstudio-pubs-static.s3.amazonaws.com/5896_8f0fed2ccbbd42489276e554a05af87e.html
https://cran.r-project.org/web/packages/survival/vignettes/adjcurve.pdf
http://debian.bjtu.edu.cn/cran/web/packages/simPH/vignettes/simPH-overview.pdf



#fit a Kaplan-Meier and plot it 
fit <- survfit(Surv(time, status) ~ x, data = aml) 
plot(fit, lty = 2:3) 
legend(100, .8, c("Maintained", "Nonmaintained"), lty = 2:3) 

#fit a Cox proportional hazards model and plot the  
#predicted survival for a 60 year old 
fit <- coxph(Surv(futime, fustat) ~ age, data = ovarian) 
plot(survfit(fit, newdata=data.frame(age=60)),
     xscale=365.25, xlab = "Years", ylab="Survival") 

# Here is the data set from Turnbull
#  There are no interval censored subjects, only left-censored (status=3),
#  right-censored (status 0) and observed events (status 1)
#
#                             Time
#                         1    2   3   4
# Type of observation
#           death        12    6   2   3
#          losses         3    2   0   3
#      late entry         2    4   2   5
#
tdata <- data.frame(time  =c(1,1,1,2,2,2,3,3,3,4,4,4),
                    status=rep(c(1,0,2),4),
                    n     =c(12,3,2,6,2,4,2,0,2,3,3,5))
fit  <- survfit(Surv(time, time, status, type='interval') ~1, 
                data=tdata, weight=n)

#
# Time to progression/death for patients with monoclonal gammopathy
#  Competing risk curves (cumulative incidence)
fit1 <- survfit(Surv(stop, event=='progression') ~1, data=mgus1,
                subset=(start==0))
fit2 <- survfit(Surv(stop, status) ~1, data=mgus1,
                subset=(start==0), etype=event) #competing risks
# CI curves are always plotted from 0 upwards, rather than 1 down
plot(fit2, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
     col=2:3, xlab="Years post diagnosis of MGUS")
lines(fit1, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
      conf.int=FALSE)
text(10, .4, "Competing Risk: death", col=3)
text(16, .15,"Competing Risk: progression", col=2)
text(15, .30,"KM:prog")