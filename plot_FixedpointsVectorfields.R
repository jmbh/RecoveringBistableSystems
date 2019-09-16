# jonashaslbeck@gmail.com; July 2019

# This file reproduces Figure 1 in the paper

library(phaseR)

# ----------------------------------------------------------------------
# --------- Setup phaseR  ----------------------------------------------
# ----------------------------------------------------------------------

# Define true bistable system collapsed into two dimensions (See Appendix A) 
# in the input format required by phaseR

# Define a function reflecitng the differential equation
DE <- function(t, y, parameters) {
  C11 <- parameters[1]
  C12 <- parameters[2]
  C13 <- parameters[3]
  C14 <- parameters[4]
  
  C21 <- parameters[5]
  C22 <- parameters[6]
  C23 <- parameters[7]
  C24 <- parameters[8]
  
  C31 <- parameters[9]
  C32 <- parameters[10]
  C33 <- parameters[11]
  C34 <- parameters[12]
  
  C41 <- parameters[13]
  C42 <- parameters[14]
  C43 <- parameters[15]
  C44 <- parameters[16]
  
  r1 <- parameters[17]
  r2 <- parameters[18]
  r3 <- parameters[19]
  r4 <- parameters[20]
  mu1 <- parameters[21]
  mu2 <- parameters[22]
  mu3 <- parameters[23]
  mu4 <- parameters[24]
  
  dy <- numeric(2)
  dy[2]<- r1*y[2] + y[2]*C11*y[2] + y[2]*C12*y[2] + y[2]*C13*y[1] + y[2]*C14*y[1]+ mu1 # negative emotions; SWITCHED!!
  dy[1]<- r3*y[1] + y[1]*C31*y[2] + y[1]*C32*y[2] + y[1]*C33*y[1] + y[1]*C34*y[1]+ mu2 # positive emotions
  
  list(dy)
}


# ----------------------------------------------------------------------
# --------- Compute Fixed pts with Wolfram Alpha -----------------------
# ----------------------------------------------------------------------

n_var <- 31
Vstress_ini <- seq(0.9, 1.1, length = n_var)

# Storage
m_FP <- as.data.frame(matrix(NA, nrow = n_var+4, ncol = 7))
m_FP[1:31, 1] <- Vstress_ini
colnames(m_FP) <- c("stress", "s1_PE", "s1_NE", "s2_PE", "s2_NE", "us_PE", "us_NE")

# Take values from Wolfram alpha
m_FP[1, 2:3] <- c(5.27602, 1.14774)
m_FP[2, 2:3] <- c(5.25653, 1.15834)
m_FP[3, 2:3] <- c(5.23642 , 1.16931)
m_FP[4, 2:3] <- c(5.21563, 1.18067)
m_FP[5, 2:3] <- c(5.19413, 1.19245)
m_FP[6, 2:3] <- c(5.17185, 1.20468)
m_FP[7, 2:3] <- c(5.14872, 1.2174)
m_FP[8, 2:3] <- c(5.12468, 1.23066)
m_FP[9, 2:3] <- c(5.09965, 1.24451)
m_FP[10, 2:7] <- c(5.07351, 1.25901, 1.82536, 3.9596, 2.02962, 3.65897)
m_FP[11, 2:7] <- c(5.04616, 1.27422, 1.65637, 4.25237, 2.2478, 3.3039)
m_FP[12, 2:7] <- c(5.01746, 1.29023, 1.568, 4.42383, 2.38572, 3.22236)
m_FP[13, 2:7] <- c(4.98726, 1.30714, 1.50262, 4.56096, 2.50213, 3.09778)
m_FP[14, 2:7] <- c(4.95535, 1.32507, 1.4495, 4.68977, 2.608, 2.99054)
m_FP[15, 2:7] <- c(4.92149, 1.34417, 1.40426, 4.78678, 2.70792, 2.89398)
m_FP[16, 2:7] <- c(4.88539, 1.36461, 1.36461, 4.88539, 2.80449, 2.80449)
m_FP[17, 2:7] <- c(4.84665, 1.38665, 1.3292, 4.97765, 2.89946, 2.71978)
m_FP[18, 2:7] <- c(4.80477, 1.4106, 1.29712, 5.0649, 2.99421, 2.63822)
m_FP[19, 2:7] <- c(4.75905, 1.43688, 1.26775, 5.14809, 3.0901, 2.55841)
m_FP[20, 2:7] <- c(4.70852, 1.46611, 1.24064, 5.22789, 3.18854, 2.47907)
m_FP[21, 2:7] <- c(4.65171, 1.49921, 1.21544, 5.30481, 3.29135, 2.39877) 
m_FP[22, 2:7] <- c(4.58618, 1.53771, 1.19189, 5.37927, 3.40121, 2.31556)
m_FP[23, 2:7] <- c(4.50749, 1.58442, 1.16976, 5.45159, 3.52281, 2.22633) 
m_FP[24, 2:7] <- c(4.40532, 1.64586, 1.14891, 5.52199, 3.6666, 2.122429) 
m_FP[25, 2:7] <- c(4.23747, 1.74897, 1.12918, 5.59073, 3.87495, 1.98229)
m_FP[26, 4:5] <- c(1.11045, 5.65796)
m_FP[27, 4:5] <- c(1.09263, 5.72383)
m_FP[28, 4:5] <- c(1.07563, 5.78848)
m_FP[29, 4:5] <- c(1.05939, 5.85202)
m_FP[30, 4:5] <- c(1.04383, 5.91452)
m_FP[31, 4:5] <- c(1.0289, 5.97609)

# extras
extra_low <- c(0.959, 5.0775, 1.25679, NA, NA, NA, NA, NA)
extra_high <- c(1.063334, NA, NA, 1.11969, 5.6253, NA, NA, NA)

m_FPex <- rbind(m_FP[1:9, ], 
                extra_low, 
                m_FP[1:25, ], 
                extra_high, 
                m_FP[26:31, ])


# library(xtable)
# row.names(m_FPex) <- NULL
# xtable(m_FPex)
# print(xtable(m_FPex), include.rownames=FALSE) # for Appendix A


# --------- (4.2) Plot Fixpoints by simulation ---------

plotFixpoints <- function() {
  
  # Fixpoints of solutions, as function of stress
  cols <- c("darkgreen", "darkgreen", "tomato", "tomato")
  
  par()$mar
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  
  par(mar = c(5.1, 10, 4.1, 10))
  
  plot.new()
  plot.window(xlim=c(0.9, 1.1), ylim=c(0, 6))
  # box()
  axis(1, c(.9, 1, 1.1))
  axis(2, 0:6, las=2)
  title(ylab = "Emotion Strength", cex.lab=1.2)
  title(xlab = expression("Stress " (r[3], r[4])), cex.lab=1.2)
  
  
  
  # plot
  lwd <- 1.5
  lines(m_FPex$stress, m_FPex$s1_PE, col = "darkgreen", lwd=lwd)
  lines(m_FPex$stress, m_FPex$s2_PE, col = "darkgreen", lwd=lwd)
  lines(m_FPex$stress, m_FPex$us_PE, col = "darkgreen", lty=2, lwd=lwd)
  
  lines(m_FPex$stress, m_FPex$s1_NE, col = "tomato", lwd=lwd)
  lines(m_FPex$stress, m_FPex$s2_NE, col = "tomato", lwd=lwd)
  lines(m_FPex$stress, m_FPex$us_NE, col = "tomato", lty=2, lwd=lwd)
  
  # for(ini in 1:2) for(i in c(1,3)) points(Vstress, a_X[, ini, i], col=cols[i], pch=20)
  legend("left", c("Positive Emotions", "Negative Emotions"), col = cols[c(1,3)], lwd=rep(lwd,4), lty=rep(1,4), bty="n", cex=1.1)
  
  # mtext("(a) Model Dynamics", side=3, adj = 0.5, line=1, cex=1.15)
  mtext(text = "(a) Fixed Points across Stress Levels", side = 3, line = 1, cex=1.15, adj=0)
  
}

plotFixpoints()


# ----------------------------------------------------------------------
# --------- Plot Vector Fields (with solution curves) with phasR -------
# ----------------------------------------------------------------------


VF_phaseR <- function(stress, 
                      xlim = c(0, 10), 
                      ylim=c(0,10), 
                      m_fixed = NULL, 
                      panel=NULL, 
                      solLines=TRUE, 
                      sol_cols = c("black", "black")) {
  
  
  # Define parameters
  
  p <- 4
  
  r <- c(stress, stress, 1, 1)
  C <- matrix(c(-.2, .04, -.2, -.2, 
                .04, -.2, -.2, -.2,
                -.2, -.2, -.2, .04, 
                -.2, -.2, .04, -.2), p, p, byrow = T)
  a <- rep(1.6, 4)
  
  # Get into phaseR format
  parameters_mxA <- c(C[1, ],
                      C[2, ],
                      C[3, ],
                      C[4, ],
                      r,
                      a = a)
  
  # Plot layout
  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  
  
  # Add Solution Curves
  if(solLines)  OU.nullclines <- nullclines(DE, xlim = xlim, 
                                            ylim = ylim, 
                                            parameters = parameters_mxA, 
                                            points = 500, 
                                            col=sol_cols,
                                            add.legend=F, 
                                            lwd = 1.5)
  
  # Add Vectorfield
  OU.flowField <- flowField(DE, xlim = xlim, ylim = ylim,
                            parameters = parameters_mxA, 
                            points = 20, 
                            add = T,
                            arrow.type = "equal",
                            xlab="Positive Emotion",ylab="Negative Emotion",cex=10,font.lab=5,
                            cex.lab=1.5, mgp=c(2,1,.5), 
                            col="#444649",
                            arrow.head = .025)
  
  # Add fixed points
  if(!is.null(m_fixed)) {
    for(i in 1:nrow(m_fixed)) {
      if(m_fixed[i, 3]==0) points(m_fixed[i, 1], m_fixed[i, 2], cex=1.4) else points(m_fixed[i, 1], m_fixed[i, 2], pch=20, cex=1.75)
    }
  }
  
  
  # Axis & labels
  axis(1, seq(0, 10, length=3))
  axis(2, seq(0, 10, length=3), las=2)
  title(ylab = "Negative Emotion", line=2.5)
  title(xlab = "Positive Emotion", line=2.5)
  
  # Titles
  mtext(text = paste0("Stress = ", stress), side = 3, line = 1, cex=1.15)
  panels = c("(b)", "(c)", "(d)")
  
  if(!is.null(panel)) mtext(panels[panel], side=3, adj = 0, line=1, cex=1.15)
  
} # eoF

# Testing
# VF_phaseR(stress=0.9, solLines=TRUE, m_fixed = m_fixed_1, sol_cols = c("blue", "lightblue"))


# ----------------------------------------------------------------------
# --------- Z: Create Figure for Paper  --------------------------------
# ----------------------------------------------------------------------



# --------- Plotting --------

pdf("figures/Figure1_TheModel.pdf", width = 8, height = 6.2)


# ----- Choose color for solution lines -----
solcols <- c("#F46D43", "#74ADD1")

# ----- Setup layout ------

lmat <- rbind(c(1,1,1), 
              c(2,3,4))

lo <- layout(lmat, heights = c(1,1))


# ------ Plot top panel: Fixed points -----

plotFixpoints()


# ------ Plot lower panel: Vector Fields -----

par(mar=c(5.1, 4.1, 4.1, 2.1))

m_fixed_0.9 <- rbind(c(m_FP[1, 2:3], 1))

m_fixed_1 <- rbind(c(m_FP[16, 2:3], 1),
                   c(m_FP[16, 4:5], 1),
                   c(m_FP[16, 6:7], 0))

m_fixed_1.1 <- rbind(c(m_FP[31, 4:5], 1))

VF_phaseR(stress = 0.9, m_fixed = m_fixed_0.9, panel = 1, solLines = T, sol_cols = solcols)
VF_phaseR(stress = 1, m_fixed = m_fixed_1, panel = 2, solLines = T, sol_cols = solcols)
VF_phaseR(stress = 1.1, m_fixed = m_fixed_1.1, panel = 3, solLines = T, sol_cols = solcols)

dev.off()


