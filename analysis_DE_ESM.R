# ryanoisin@gmail.com & jonashaslbeck@gmail.com ; August 2019

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

data <- readRDS(file = "files/data_ESM.RDS")
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

# --------- Tables 2 and (first column) Tabel 7 ----------
xtable(ftab[c(1,2,4,5,6,7,9),])

# ----------------------------------------------------------------------
# --------- Dynamics of Final Model ------------------------------------
# ----------------------------------------------------------------------


# ---------- Get Final Model Parameter Estimates --------
out <- getDEpars(D, model="All_Fourw_Int", residout = T)

# ---------- Figure 12 (left hand side) -----------------
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


# -------- Plot Vector Field -------------------------

# Use only parameters relating to dx1/dt and dx3/dt
# (Collapse system to a 2-dimensional case)
parameters <- c(out$pars[1,],out$pars[3,])


# Set plotting parameters
xlim = c(0, 10) ; ylim = c(0, 10) ; Range = 10
sc <- .6

solcols <- c("#F46D43", "#74ADD1")

pdf("figures/Figure12_VF_Est_ESM_snap.pdf",width = 10*sc, height = 10*sc)
    par(mfrow = c(1, 1))    
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    OU.flowField <- flowField(DE, xlim = xlim, ylim = ylim,
                          parameters = parameters, points = 30, add = T,
                          arrow.type = "equal",
                          xlab="Positive Emotion",ylab="Negative Emotion",cex=10,font.lab=5,
                          cex.lab=1.5, mgp=c(2,1,.5),
                          col="#444649",
                          arrow.head = .05)
        axis(1, seq(0, 10, length=3))
        axis(2, seq(0, 10, length=3), las=2)
        title(ylab = "Negative Emotion")
        title(xlab = "Positive Emotion")

# Add solution lines
      OU.nullclines <- nullclines(DE, xlim = c(xlim[1]-Range/5, xlim[2]+Range/5), 
                            ylim = c(ylim[1]-Range/5, ylim[2]+Range/5),
                            parameters = parameters, points = 500, 
                            col = solcols,
                            add.legend=F,
                            lwd = 2)

# search for equilibrium points (starting based on visual inspection)
    fp1 <- findEquilibrium(DE, y0 = c(5,1), parameters = parameters,
                       system = "two.dim", tol = 1e-16, 
                       max.iter = 50, h = 1e-06,
                       plot.it = FALSE, summary = TRUE, 
                       state.names = c("x", "y"))
      # Plot fixed point
        points(x = fp1$ystar[1,1], y = fp1$ystar[2,1], 
               pch = 20, cex = 2)

    fp2 <- findEquilibrium(DE, y0 = c(6,6), parameters = parameters,
                       system = "two.dim", tol = 1e-16, 
                       max.iter = 50, h = 1e-06,
                       plot.it = FALSE, summary = TRUE, 
                       state.names = c("x", "y"))
        points(x = fp2$ystar[1,1], y = fp2$ystar[2,1], 
               cex = 1.4, lwd = 2)

    fp3 <- findEquilibrium(DE, y0 = c(3,3), parameters = parameters,
                       system = "two.dim", tol = 1e-16, 
                       max.iter = 50, h = 1e-06,
                       plot.it = FALSE, summary = TRUE, 
                       state.names = c("x", "y"))
        points(x = fp3$ystar[1,1], y = fp3$ystar[2,1], 
               cex = 1.4, lwd = 2)

    fp4 <- findEquilibrium(DE, y0 = c(1,5), parameters = parameters,
                       system = "two.dim", tol = 1e-16, 
                       max.iter = 50, h = 1e-06,
                       plot.it = FALSE, summary = TRUE, 
                       state.names = c("x", "y"))
        points(x = fp4$ystar[1,1], y = fp4$ystar[2,1], 
               pch = 20, cex = 2)

dev.off()


# ----------------------------------------------------------------------
# --------- Data Generated from Final Model ----------------------------
# ----------------------------------------------------------------------

set.seed(5) 
time <- 14 * 24 * 60 / 90
ts <- 1 # use time scale at which data is observed

# Set starting value
start <- c(rep(1.3,2), rep(4.8,2))
# Scale residual variance 
rdscale <- .65

dataG <- genData_est(time = time, 
                     feat = out$feat,
                     coef = out$pars, 
                     sds = rdscale*out$sds,
                     timestep = ts, 
                     initial = start,
                     noise = TRUE)

# Thin for display purposes
thin_scale <- 1
thinnerG <- seq(1, nrow(dataG), length = round(nrow(dataG)) / thin_scale)
dataG_thinned <- dataG[thinnerG, ]


# New plot Jonas: 
pdf("figures/Figure13_DataGen_DE_ESM_snap.pdf", width = 10, height = 4.5)

  # --- Set up Plot ---

  lmat <- rbind(1:2)
  lo <- layout(mat = lmat, widths = c(1, .4))

  # --- Plot: full 2 weeks ---

  plotTimeSeriesFirst2Weeks(data_thinned = dataG_thinned, 
                          label=("(a)"), mar = c(3.5,3.5,4.1,1),
                          cex.lab = 1.5,
                          line.lab = 2.25)

  # --- Plot: switch outtake ---

  # Define switch location 
  switch_start <- 150
  switch_end <- 182
  data_switch <- dataG[switch_start:switch_end, ]
  
  thin_scale <- 1
  rect(xleft = switch_start / thin_scale, # to make time series visible inside box
       ybottom = 0, 
       xright = switch_end / thin_scale, 
       ytop = 8,
       lty = 2)
  
  
  # Plotting
  cols <- brewer.pal(n = 4, "Set1")
  
  par(mar=c(3.5,0,4.1,2.1))
  plot.new()
  plot.window(xlim=c(1, nrow(data_switch)), ylim=c(0, 8))
  for(i in 1:4) lines(data_switch[, i+1], col = cols[i], lwd=0.5)
  box(lty=2)
  
  legend("topleft", c("Cheerful", "Content", "Anxious", "Sad"), 
         lwd = rep(1.5, 4), col = cols, bty='n', cex=0.8)
  
  mtext(text = "(b)", adj=0, line=1, cex=1.5)
  
dev.off()

