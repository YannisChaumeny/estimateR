source('code/estimate.R')
source('code/getData.R')
source('code/plot.R')
source('code/simulationData.R')


#--------- Motivation ----------------------------------------------------------

#plot weekly patterns in SWitzerland
dataCH <- getData('Switzerland')
selCH <- (dataCH$N-100):dataCH$N
pdf("report/weeklyPattern.pdf",width=8,height=4)
plot(
  dataCH$date[selCH],
     dataCH$C[selCH],
     type="l",
     xlab="Time (days)",
     ylab="New confirmed cases",
     )
points(dataCH$date[intersect(which(wday(dataCH$date,week_start=1)==1),selCH)],
       dataCH$C[intersect(which(wday(dataCH$date,week_start=1)==1),selCH)],
       col = "#ff2a00")
dev.off()

#--------- Simulation data -----------------------------------------------------

dataSim <- simulationData(p = c(0.2,0.1,0.1,0.1,0.1,0.2,0.1,0.2), n_t = 200)
sel <- 1:160

pdf("report/simData.pdf",height=4,width=8)
plot(dataSim$trueI[sel],type="l",ylab="Simulated daily incidence",xlab="Time (days)",col="red")
lines(dataSim$C[sel])
points(which(dataSim$wkday[sel]==1),
       dataSim$C[which(dataSim$wkday[sel]==1)],col="orange")
legend("topright",
       inset=0.02,
       legend=c("Infections", "Observations", "Mondays"),
       col=c("red", "black", "orange"),
       pch=c(NA, NA,1),
       lty = c(1,1,0),
       pt.cex = 1,
       seg.len = 0.7
)
dev.off()

#--------- Results -------------------------------------------------------------

#-- 1. Simulation study ------

# Non delayed simulated data: compare with Cori (side by side?)
infSim <- simulationData(weekend=FALSE,n_t=200)
selec = 8:160
stan_data <- list(
  N = infSim$N,
  S = 20,       
  I = infSim$trueI,
  I_upper = 1e5,
  sigmaR = 0.07,
  SI_mean = infSim$GI_mean,
  SI_std = infSim$GI_std
)
fitInfSim <- stan("epiModel.stan", data = stan_data, chains = 4,
                  iter = 2000, warmup = 400)
R_summary <- summary(fitInfSim,pars="R")$summary
pdf("report/infSimR.pdf",width=5,height=4)
plot_R(R_summary[selec,],selec,trueR = infSim$trueRt[selec])
dev.off()
library(EpiEstim)
epiestim <- estimate_R(infSim$trueI,"parametric_si",
                        config = make_config(list(
                          mean_si = infSim$GI_mean,
                          std_si = infSim$GI_std,
                          t_start = 2:infSim$N,
                          t_end = 2:infSim$N)
                        )
)
pdf("report/cori.pdf",width=5,height=4)
plot(x = selec,
     y = epiestim$R$`Mean(R)`[selec],
     type = "l",
     ylim = c(0,1.2*max(epiestim$R$`Mean(R)`[selec])),
     xlab = "Time (days)",
     ylab = "Reproduction number")
polygon(c(rev(selec), selec),
        c(rev(epiestim$R$`Quantile.0.025(R)`[selec]),
          epiestim$R$`Quantile.0.975(R)`[selec]),
        col = 'grey80',
        border = NA)
polygon(c(rev(selec), selec),
        c(rev(epiestim$R$`Quantile.0.25(R)`[selec]),
          epiestim$R$`Quantile.0.75(R)`[selec]),
        col = 'grey60',
        border = NA)
lines(selec,infSim$trueRt[selec],col="orange")
lines(selec,epiestim$R$`Mean(R)`[selec],type="l")
abline(1,0,lty="dashed")
legend("topright",
       inset=0.02,
       legend=c("R estimate", "50% CI", "95% CI","True R"),
       col=c("black", "grey60", "grey80", "orange"),
       pch=c(NA, 15,15,NA),
       lty = c(1,0,0,1),
       pt.cex = 2,
       seg.len = 0.7
)
dev.off()
# Simulated with weekly delay: model without and with weekly spec

selecSimComp <- 1:160
dataSim <- simulationData(p = c(0.2,0.1,0.1,0.1,0.1,0.2,0.1,0.2), n_t = 200)
stan_data <- list(
  C = dataSim$C,
  N = dataSim$N,
  I_upper = 1e5,
  sigmaR = 0.07,
  SI_mean = dataSim$GI_mean,
  SI_std = dataSim$GI_std,
  d = dataSim$d,
  K = length(dataSim$d),
  S = 20
)
fitSimNW <- stan("withoutWeekly.stan",data = stan_data,
                 iter = 2000, warmup = 400,chains = 4)
fitSimNW_R_summary <- summary(fitSimNW,pars="R")$summary
pdf("report/simNW.pdf",width=5,height=4)
plot_R(fitSimNW_R_summary[selecSimComp,],selecSimComp,
       trueR=dataSim$trueRt[selecSimComp])
dev.off()

fitSimW <- estimate(dataSim$C,simulation = TRUE,
                    SI_mean = dataSim$GI_mean, SI_std = dataSim$GI_std,
                    sigma = 0.07, iter = 1000, warmup = 400)
pdf("report/simW.pdf",width=5,height=4)
plot_R(fitSimW$R_summary[selecSimComp,],selecSimComp,
       trueR=dataSim$trueRt[selecSimComp])
dev.off()

# Estimation of params for the above

pdf("report/simP.pdf",width=7,height=4)
plot_pr(fitSimW$pr_summary)
points(1:8,c(0.2,0.1,0.1,0.1,0.1,0.2,0.1,0.2),
       col = 'orange', pch = 20)
dev.off()

#-- 2. National data -------

# Swiss infections

dataCH <- getData('Switzerland')
fitCH <- estimate(dataCH$C,dataCH$date,iter = 2000, warmup = 500)

pdf("report/CH_I.pdf",width = 8, height = 4)
plot_I(fitCHtest$I_summary[1:dataCH$N,], dataCH$date, dataCH$C)
dev.off()

# Swiss weekly

pdf("report/CH_P.pdf",width = 7, height = 4)
plot_pr(fitCHtest$pr_summary)
dev.off()

# Swiss Rt

pdf("report/CH_R.pdf", width = 8, height = 4)
plot_R(fitCHtest$R_summary[1:(dataCH$N-10),], dataCH$date[1:(dataCH$N-10)])
dev.off()

# taskforce Rt

ncs <- read.table(header = TRUE, sep=",",
                  "https://raw.githubusercontent.com/covid-19-Re/dailyRe-Data/master/CHE-estimates.csv")
ncs <- ncs[which(ncs$region=="CHE" & ncs$data_type == "Confirmed cases" & ncs$estimate_type == "Cori_slidingWindow"),]

ncs <- ncs[19:nrow(ncs),]
pdf("report/ncs.pdf", width = 8, height = 4)
plot(x = as.Date(ncs$date),
     y = ncs$median_R_mean,
     type = "l",
     ylim = c(0,1.2*max(ncs$median_R_mean)),
     xlab = "Time (days)",
     ylab = "Reproduction number")
polygon(c(rev(as.Date(ncs$date)), as.Date(ncs$date)),
        c(rev(ncs$median_R_lowHPD),ncs$median_R_highHPD),
        col = 'grey80',
        border = NA)
lines(as.Date(ncs$date),ncs$median_R_mean,type="l")
abline(1,0,lty="dashed")
legend("topright",
       inset=0.02,
       legend=c("R estimate (Huisman et al.)", "95% CI"),
       col=c("black", "grey80"),
       pch=c(NA, 15),
       lty = c(1,0),
       pt.cex = 2,
       seg.len = 0.7
  )
dev.off()



#----- Supplimentary -----------------------------------------------------------


estimateNational <- function(country){
  data <- getData(country)
  selection <- (data$N-200):(data$N)
  data$C <- data$C[selection]
  data$date <- as.Date(data$date[selection])
  data$N <- length(data$C)
  fit <- estimate(data$C,data$date,
                  iter = 2000,warmup = 500,delay_type = 1)
  
  pdf(paste("report/",country,"_I.pdf",sep=""),width = 8, height = 4)
  plot_I(fit$I_summary[1:data$N,], data$date, data$C)
  dev.off()

  pdf(paste("report/",country,"_P.pdf",sep=""),width = 7, height = 4)
  plot_pr(fit$pr_summary)
  dev.off()

  pdf(paste("report/",country,"_R.pdf",sep=""), width = 8, height = 4)
  plot_R(fit$R_summary[1:data$N,], data$date)
  dev.off()
}
#countryList <- unique(worldData$Country)
estimateNational("Belgium")
estimateNational("Germany")
estimateNational("Austria")


