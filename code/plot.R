plot_R <- function(R_summary,date,name = "R_plot",trueR = NULL){
  plot(x = date,
         y = R_summary[,"mean"],
         type = "l",
         ylim = c(0,1.2*max(R_summary[,"mean"])),
         xlab = "Time (days)",
         ylab ="Reproduction number")
  polygon(c(rev(date), date),
            c(rev(R_summary[,"2.5%"]),R_summary[,"97.5%"]),
            col = 'grey80',
            border = NA)
  polygon(c(rev(date), date),
          c(rev(R_summary[,"25%"]),R_summary[,"75%"]),
          col = 'grey60',
          border = NA)
  if(!is.null(trueR)){
    lines(date,trueR,col="orange")
  }
  lines(date,R_summary[,"mean"],type="l")
  abline(1,0,lty="dashed")
  if(!is.null(trueR)){
    legend("topright",
           inset=0.02,
           legend=c("R estimate", "50% CI", "95% CI","True R"),
           col=c("black", "grey60", "grey80", "orange"),
           pch=c(NA, 15,15,NA),
           lty = c(1,0,0,1),
           pt.cex = 2,
           seg.len = 0.7
          )
  }
  else{
    legend("topright",
         inset=0.02,
         legend=c("R estimate", "50% CI", "95% CI"),
         col=c("black", "grey60", "grey80"),
         pch=c(NA, 15,15),
         lty = c(1,0,0),
         pt.cex = 2,
         seg.len = 0.7
        )
  }
  
  }

plot_I <- function(I_summary,date,cases,name="I_plot"){
  I_summary <- I_summary[-1,] #first day is not estimated
  date <- date[-1]
  
  plot(x = date,
       y = I_summary[,"mean"],
       type = "l",
       col = "blue",
       ylim = c(0,1.2*max(I_summary[,"mean"])),
       xlab = "Time (days)",
       ylab = "Daily incidence")
  polygon(c(rev(date), date),
          c(rev(I_summary[,"2.5%"]),I_summary[,"97.5%"]),
          col = 'grey80',
          border = NA)
  polygon(c(rev(date), date),
          c(rev(I_summary[,"25%"]),I_summary[,"75%"]),
          col = 'grey60',
          border = NA)
  lines(date,cases[-1],col="orange")
  lines(date,I_summary[,"mean"])
  legend("topleft",
         inset=0.02,
         legend=c("Estimated infections", "50% CI", "95% CI","Observations"),
         col=c("black", "grey60", "grey80","orange"),
         pch=c(NA, 15,15,NA),
         lty = c(1,0,0,1),
         pt.cex = 2,
         seg.len = 0.7
  )
}
  
plot_pr <- function(pr_summary){
  n <- nrow(pr_summary)
  plot(pr_summary[,"mean"],pch = 20,ylim=c(0,1),
       ylab = "Weekly paramaters",xlab=NA,xaxt="n")
  grid()
  points(1:n, pr_summary[,"2.5%"],pch=4)
  points(1:n, pr_summary[,"97.5%"],pch=4)
  axis(1,1:8,labels = c(expression(P[Mo %->% Tu]),expression(P[Tu %->% We]),
                        expression(P[We %->% Th]),expression(P[Th %->% Fr]),
                        expression(P[Fr %->% Sa]),expression(P[Sa %->% Mo]),
                        expression(P[Sa %->% Su]),expression(P[Su %->% Mo]))
  )
  legend("topright",
         inset=0.02,
         legend=c("P estimates", "95% CI"),
         pch=c(20,4)
  )
}
