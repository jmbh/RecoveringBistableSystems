# jonashaslbeck@gmail.com & ryanoisin@gmail.com; August 2019

# ----------------------------------------------------------------------
# --------- Load packages & source  ------------------------------------
# ----------------------------------------------------------------------

source("fun_StatsModels.R") # For plotting
source("fun_DEestimation.R") # For estimation
library(phaseR)
library(xtable)


# ----------------------------------------------------------------------
# --------- Load and Prepare Data  -------------------------------------
# ----------------------------------------------------------------------

# ----- Load data -----
data <- readRDS(file = "files/data.RDS")
n <- nrow(data)
dt <- data[2,1] - data[1,1]

# Compute derivative of x1
dx1dt <- (data[, 2][-1] - data[, 2][-n])/dt
dx2dt <- (data[, 3][-1] - data[, 3][-n])/dt
dx3dt <- (data[, 4][-1] - data[, 4][-n])/dt
dx4dt <- (data[, 5][-1] - data[, 5][-n])/dt

# Delete last observation
x1 <- data[-n, 2]
x2 <- data[-n, 3]
x3 <- data[-n, 4]
x4 <- data[-n, 5]
D <- as.data.frame(cbind(dx1dt,dx2dt,dx3dt,dx4dt, x1, x2, x3, x4))

# Clear remaining objects from working memory
rm(dx1dt,dx2dt,dx3dt,dx4dt)

# ----------------------------------------------------------------------
# --------- Model Building  --------------------------------------------
# ----------------------------------------------------------------------

# ------ Fit range of models -------

# set a random seed for reproducibility
set.seed(411)

# Compare models by CV R^2
fit <- model_compare(D = D, k =10)

# Arrange fit into a table
ftab <- round(fit[,c(5,6)],5)

# --------- Table with Models A to C ----------
 xtable(ftab[c(1,2,4),])

# --------- Table with full fit results (Appendix) ----------
 xtable(ftab[c(1,2,4,5,6,7,9),])

# ----------------------------------------------------------------------
# --------- Dynamics of Final Model ------------------------------------
# ----------------------------------------------------------------------

# ---------- Get Final Model Parameter Estimates --------
out <- getDEpars(D=D, model="Main_Int", residout = T)

# ---------- Figure 9 (left hand side) -----------------
ahat <- round(out$pars[,1],2) # intercepts
Rhat <- round(out$pars[1:4,2:5],2) # R matrix
Chat <- rbind(round(out$pars[1,6:9],2),
              round(out$pars[2,c(7,10,11,12)],2),
              round(out$pars[3,c(8,11,13,14)],2),
              round(out$pars[4,c(9,12,14,15)],2))
sigmahat <- round(unlist(lapply(out$resids,sd)) , 2)


# ---------- Table 6 (unedited) -----------------
coeftab <- cbind(
  round(out$summaries[[1]]$coefficients[,c(1,2,4)],3),
  round(out$summaries[[2]]$coefficients[,c(1,2,4)],3),
  round(out$summaries[[3]]$coefficients[,c(1,2,4)],3),
  round(out$summaries[[4]]$coefficients[,c(1,2,4)],3))

xtable(coeftab)


# ---------- Plot Vector Field ------------------

# Use only parameters relating to dx1/dt and dx3/dt
# (Collapse system to a 2-dimensional case)
parameters <- c(out$pars[1,],out$pars[3,])

# Set plotting parameters
xlim = c(0, 10) ; ylim = c(0, 10) ; Range = 10
sc <- .6

solcols <- c("#F46D43", "#74ADD1")

pdf("figures/Figure9_VF_Est_Ideal.pdf",width = 10*sc, height = 10*sc)
  par(mfrow = c(1, 1))
  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  # Plot VF
  OU.flowField <- flowField(DE, xlim = xlim, ylim = ylim,
                          parameters = parameters, points = 30, add = T,
                          arrow.type = "equal",
                          xlab="Positive Emotion",ylab="Negative Emotion",cex=10,font.lab=5,
                          cex.lab=1.5, mgp=c(2,1,.5),
                          col="#444649",
                          arrow.head = .04)
    axis(1, seq(0, 10, length=3))
    axis(2, seq(0, 10, length=3), las=2)
    title(ylab = "Negative Emotion")
    title(xlab = "Positive Emotion")

  # Add solution lines
  OU.nullclines <- nullclines(DE, xlim = c(xlim[1]-Range/5, xlim[2]+Range/5), 
                            ylim = c(ylim[1]-Range/5, ylim[2]+Range/5),
                            parameters = parameters, points = 500, 
                            col=solcols,
                            add.legend=F, lwd = 2)

  # Find Unstable Fixed Point
  fp_uns <- findEquilibrium(DE, y0 = c(2.5,2.5), parameters = parameters,
                            system = "two.dim", tol = 1e-16, 
                            max.iter = 50, h = 1e-06,
                            plot.it = F, summary = TRUE, 
                            state.names = c("x", "y"))
    # Plot equilibrium point
    points(x = fp_uns$ystar[1, 1], y = fp_uns$ystar[2, 1],  
            cex = 1.4, lwd = 2)

  # Find Stable Fixed Point
  fp_s1 <- findEquilibrium(DE, y0 = c(1,5), parameters = parameters,
                         system = "two.dim", tol = 1e-16, 
                         max.iter = 50, h = 1e-06,
                         plot.it = F, summary = TRUE, 
                         state.names = c("x", "y"))
    # Plot equilibrium point
    points(x = fp_s1$ystar[1, 1], y = fp_s1$ystar[2, 1],  
          pch = 20, cex = 2)
    
  # Find Stable Fixed Point
  fp_s2 <- findEquilibrium(DE, y0 = c(5,1), parameters = parameters,
                         system = "two.dim", tol = 1e-16, 
                         max.iter = 50, h = 1e-06,
                         plot.it = F, summary = TRUE, 
                         state.names = c("x", "y"))
  # Plot equilibrium point
    points(x = fp_s2$ystar[1, 1], y = fp_s2$ystar[2, 1], 
          pch = 20, cex = 2)

dev.off()


# ----------------------------------------------------------------------
# --------- Data Generated from Final Model ----------------------------
# ----------------------------------------------------------------------

set.seed(123)
time <- 20160
ts <- .1 # use time scale at which data is observed

dataG <- genData_est(time = time, 
                     feat = out$feat,
                     coef = out$pars, 
                     sds = out$sds,
                     timestep = ts, 
                     initial = data[1,-1],
                     noise = TRUE)

# Thin for display purposes
thin_scale <- 50 # every 50 time points
thinnerG <- seq(1, nrow(dataG), length = round(nrow(dataG)) / thin_scale)
dataG_thinned <- dataG[thinnerG, ]


pdf("figures/Figure10_DataGen_DE_ideal.pdf", width = 10, height = 4.5)
  plotRegimes2Weeks(data =  dataG_thinned, regimes = rep(1, nrow(dataG_thinned)), 
                  cols_regimes = c("grey", "white"), 
                  legend = TRUE, plot_regimes = FALSE, 
                  ylim = c(0,8), y_axis_label = c(0,2,4,6,8),
                  cex.lab = 1.5,
                  line.lab = 2.25)
dev.off()



