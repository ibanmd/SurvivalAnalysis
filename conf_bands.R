#########################################################
aLaU <- function(otime, censor, t1,t2) {
  # This function computes a_L and a_U defined in (4.4.1), p101.
  # otime: the observed time
  # censor: the censoring indicator
  # (t1, t2): the interval on which the confidence bands will be computed
  ### Output of this function
  # a_L and a_U: defined in (4.4.1), p101.
  n <- length(otime)
  ples <- summary(survfit(Surv(otime, censor)~1))
  surv <- ples$surv
  dtime <- ples$time
  t1 <- max(dtime[t1-dtime >= 0])
  t2 <- max(dtime[t2-dtime >= 0]) 
  se <- ples$std.err
  sigma1 <- se[dtime==t1]^2/surv[dtime==t1]^2
  sigma2 <- se[dtime==t2]^2/surv[dtime==t2]^2
  aL <- n*sigma1/(1+n*sigma1)     
  aU <- n*sigma2/(1+n*sigma2)     
  return(list(aL=aL, aU=aU))
}

#########################################################
EP.bands <- function(otime,censor,t1,t2,calpha) {
  # This function computes the Equal Probability (EP) confidence
  # bands, including the linear, log-log transformed and
  # the arcsin-square root transformed bands.
  # otime: the observed time
  # censor: the censoring status of otime.
  # calpha: after computed a_L and a_U from the R function aLaU 
  # above, we need to find calpha from Table C3. 
  # (t1, t2): the interval on which the confidence bands will be
  # computed. Actually, the confidence bands is computed from the
  # largest event time samller than t1 to the largest time smaller
  # than t2. I did this because the KM estimator is step-wise         
  # and right-continuous.    
  #################################
  ### Output of this function
  # linSL: lower linear band
  # linSU: upper linear band
  # logSL: lower log-log band
  # logSU: upper log-log band
  # arcSL: lower arcsin band
  # arcSU: upper arcsin band
  # PLE:   Product-Limit Estimator
  # sigma2: sigma^2, which is defined in the 2nd paragraph, p97. 
  n <- length(otime)
  ples <- summary(survfit(Surv(otime, censor)~1))
  dtime <- ples$time
  realt2 <- t2
  realt1 <- t1
  t1 <- max(dtime[t1-dtime >= 0])
  t2 <- max(dtime[t2-dtime >= 0])
  rng <- dtime>=t1 & dtime<=t2
  dtime <- dtime[rng]                     
  survf <- ples$surv[rng]                      
  se <- ples$std.err[rng]                     
  sigma <- se/survf
  #linear bands
  linSL <- survf*(1 - calpha*sigma)
  linSL[linSL < 0] <- 0
  linSU <- survf*(1+calpha*sigma)
  linSU[linSU > 1] <- 1
  #log-log transformed bands
  theta <- exp(calpha*sigma/log(survf))
  logSL <- survf^(1/theta)
  logSU <- survf^(theta)
  #arcsin bands
  tmp1 <- (asin(survf^(.5))-.5*calpha*sigma*(survf/(1-survf))^(.5))
  tmp1[tmp1 < 0 ] <- 0
  tmp2 <- (asin(survf^(.5))+.5*calpha*sigma*(survf/(1-survf))^(.5))
  tmp2[tmp2 > .5*pi] <- .5*pi
  arcSL <- (sin(tmp1))^2
  arcSU <- (sin(tmp2))^2  
  surv1 <- sort(-c(survf,survf))
  time1 <- sort(c(dtime,dtime))
  m <- length(time1)
  linSL1 <- -sort(-c(linSL, linSL))[-m]
  linSU1 <- -sort(-c(linSU, linSU))[-m]
  logSL1 <- -sort(-c(logSL, logSL))[-m]
  logSU1 <- -sort(-c(logSU, logSU))[-m]
  arcSL1 <- -sort(-c(arcSL, arcSL))[-m]
  arcSU1 <- -sort(-c(arcSU, arcSU))[-m]
  plot(c(t1,t2),c(1,0), xlim=c(t1,t2),type='n', xlab="", ylab="")
  title(main="Equal Probability Confidence Bands",
        xlab="Time",ylab="Survival Probability")
  lines(time1[-1],-surv1[-m], col="black")
  lines(time1[-1],linSL1,lty=2, col="blue")
  lines(time1[-1],linSU1,lty=2, col="blue")
  lines(time1[-1],logSL1,lty=2, col="green")
  lines(time1[-1],logSU1,lty=2, col="green")
  lines(time1[-1],arcSL1,lty=4, col="red")
  lines(time1[-1],arcSU1,lty=4, col="red")
  legend(min(time1), 0.30, c("linear", "log", "arcsin"),
         col=c("blue", "green", "red"), lty=2:4) 
  output <- list(linSL=linSL, linSU=linSU, logSL=logSL, logSU=logSU, 
                 arcSL=arcSL, arcSU=arcSU, PLE=survf, sigma2=sigma^2, 
                 time = dtime, surv = survf)
  output <- data.frame(output)
  output[length(output[,1])+1,] <- output[length(output[,1]),]
  output$time[length(output$time)] <- realt2
  output$time[1] <- realt1
  return(output)
}

########################################################################
HW.bands <- function(otime,censor,t1,t2,kalpha) {
  # This function computes the Hall-Wellner (HW) confidence
  # bands, including the linear, log-log transformed and
  # the arcsin-square root transformed bands.
  # otime: the observed time
  # censor: the censoring status of otime.
  # kalpha: after computed a_L and a_U from the R function aLaU 
  # above, we need to find kalpha from Table C4. 
  # (t1, t2): the interval on which the confidence bands will be
  # computed. Actually, the confidence bands is computed from the
  # largest event time samller than t1 to the largest time smaller
  # than t2. I did this because the KM estimator is step wise
  # and right-continuous. 
  ##################################
  ### Output of this function
  # linSL: lower linear band
  # linSU: upper linear band
  # logSL: lower log-log band
  # logSU: upper log-log band
  # arcSL: lower arcsin band
  # arcSU: upper arcsin band
  # PLE:   Product-Limit Estimator
  # sigma2: sigma^2, which is defined in the 2nd paragraph, p97.
  n <- length(otime)
  ples <- summary(survfit(Surv(otime, censor)~1))
  dtime <- ples$time
  realt1 <- t1
  realt2 <- t2
  t1 <- max(dtime[t1-dtime >= 0])
  t2 <- max(dtime[t2-dtime >= 0])
  rng <- dtime>=t1 & dtime<=t2
  dtime <- dtime[rng]
  survf <- ples$surv[rng]
  se <- ples$std.err[rng]
  sigma <- se/survf
  #linear conf bands
  linSL <- survf*(1 - kalpha*(1+n*sigma^2)*n^(-0.5))
  linSL[linSL < 0] <- 0
  linSU <- survf*(1 + kalpha*(1+n*sigma^2)*n^(-0.5))
  linSU[linSU > 1] <- 1
  #log-log transformed conf bands 
  theta <- exp(kalpha*(1+n*sigma^2)*n^(-0.5)/log(survf))
  logSL <- survf^(1/theta)
  logSU <- survf^(theta)
  #arcsin transformed bands
  tmp1<-(asin(survf^(.5))-
           .5*kalpha*(1+n*sigma^2)*n^(-0.5)*(survf/(1-survf))^(.5))
  tmp1[tmp1 < 0 ] <- 0
  tmp2 <- (asin(survf^(.5))+
             .5*kalpha*(1+n*sigma^2)*n^(-0.5)*(survf/(1-survf))^(.5))
  tmp2[tmp2 > .5*pi] <- .5*pi
  arcSL <- (sin(tmp1))^2
  arcSU <- (sin(tmp2))^2
  surv1 <- sort(-c(survf,survf))
  time1 <- sort(c(dtime,dtime))
  m <- length(time1)
  linSL1 <- -sort(-c(linSL, linSL))[-m]
  linSU1 <- -sort(-c(linSU, linSU))[-m]
  logSL1 <- -sort(-c(logSL, logSL))[-m]
  logSU1 <- -sort(-c(logSU, logSU))[-m]
  arcSL1 <- -sort(-c(arcSL, arcSL))[-m]
  arcSU1 <- -sort(-c(arcSU, arcSU))[-m]
  plot(c(t1,t2),c(1,0), xlim=c(t1,t2),type='n', xlab="", ylab="")
  title(main="Hall-Wellner Confidence Bands",
        xlab="time",ylab="Survival Probability")
  lines(time1[-1],-surv1[-m], col="black")
  lines(time1[-1],linSL1,lty=2, col="blue")
  lines(time1[-1],linSU1,lty=2, col="blue")
  lines(time1[-1],logSL1,lty=2, col="green")
  lines(time1[-1],logSU1,lty=2, col="green")
  lines(time1[-1],arcSL1,lty=4, col="red")
  lines(time1[-1],arcSU1,lty=4, col="red")
  legend(min(time1), 0.3, c("linear", "log", "arcsin"),
         col=c("blue", "green", "red"), lty=2:4)
  output <- list(linSL=linSL, linSU=linSU, logSL=logSL, logSU=logSU,
                 arcSL=arcSL, arcSU=arcSU, PLE=survf, sigma2=sigma^2, 
                 time = dtime, surv = survf)
  output <- data.frame(output)
  output[length(output[,1])+1,] <- output[length(output[,1]),]
  output$time[length(output$time)] <- realt2
  output$time[1] <- realt1
  return(output)
}

