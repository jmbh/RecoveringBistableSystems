# ryanoisin@gmail.com & jonashaslbeck@gmail.com; August 2019

# ----------------------------------------------------------------------
# ----- Model Building: Fit 10 Models of increasing complexity  --------
# ----------------------------------------------------------------------


model_compare <- function(D, k = 10){
  
  # Randomly shuffle data
  D <- D[sample(nrow(D)),]
  
  # define folds
  folds <- cut(seq(1,nrow(D)),breaks=k,labels=FALSE)
  
  # How many models are there?
  nmod <- 10
  
  # Store R^2 values
  r2ar <- array(0, c(k, 4, nmod))
  
  # Store number of parameters in each model
  p <- rep(0,nmod)
  
  # Start cross-validation loop
  for(f in 1:k){
    testIndexes <- which(folds==f,arr.ind=TRUE)
    testData <- D[testIndexes, ]
    trainData <- D[-testIndexes, ]
    
    # ------------------------------------------------------------
    # -------------------------"Main"-----------------------------
    # ------------------------------------------------------------
    
    # Main effects only model
    
    # Start loop over variables
    for(i in 1:4){
      
      # Train model
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4)
      c <- lm_obj$coefficients
      
      # Make predictions
      preds <- c[1] + c[2]*testData$x1 +  
               c[3]*testData$x2 +  c[4]*testData$x3 +
               c[5]*testData$x4
      
      # Caculate R^2
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,1] <- r2
      
    }
    p[1] <- length(c)
    
    
    # ------------------------------------------------------------
    # -------------------------"Main_Int"-------------------------
    # ------------------------------------------------------------
    
    # Linear (main effectrs) plus Two-way interactions with DV model
    for(i in 1:4){
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4 +
                     I(trainData[,eval(paste0("x",i))] * trainData$x1) + 
                     I(trainData[,eval(paste0("x",i))] * trainData$x2) +
                     I(trainData[,eval(paste0("x",i))] * trainData$x3) + 
                     I(trainData[,eval(paste0("x",i))] * trainData$x4))
      c <- lm_obj$coefficients
      
      # Make predictions
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  
               c[4]*testData$x3 +  c[5]*testData$x4 +
               c[6]*testData[,eval(paste0("x",i))] * testData$x1 + 
               c[7]*testData[,eval(paste0("x",i))] * testData$x2 +
               c[8]*testData[,eval(paste0("x",i))] * testData$x3 + 
               c[9]*testData[,eval(paste0("x",i))] * testData$x4
      
      # Caculate R^2
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,2] <- r2
    }
    p[2] <- length(c)
    
    # ------------------------------------------------------------
    # ----------------"Main_Int_Quad"-----------------------------
    # ------------------------------------------------------------
    
    # All linear, quadratic, plus two-way interactions with DV
    for(i in 1:4){
      # Add missing quadratic effects from previous model
      notdv <- (1:4)[-i]
      
      lm_obj  <- lm(trainData[,eval(paste0("dx",i,"dt"))]  ~ 
                      trainData$x1 +  trainData$x2 +  
                      trainData$x3 +  trainData$x4 +
                      I(trainData[,eval(paste0("x",i))] * trainData$x1) + 
                      I(trainData[,eval(paste0("x",i))] * trainData$x2) +
                      I(trainData[,eval(paste0("x",i))] * trainData$x3) + 
                      I(trainData[,eval(paste0("x",i))] * trainData$x4) + 
                      I(trainData[,eval(paste0("x",notdv[1]))]^2) +
                      I(trainData[,eval(paste0("x",notdv[2]))]^2) +
                      I(trainData[,eval(paste0("x",notdv[3]))]^2) )
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  
               c[4]*testData$x3 +  c[5]*testData$x4 +
               c[6]*testData[,eval(paste0("x",i))] * testData$x1 + 
               c[7]*testData[,eval(paste0("x",i))] * testData$x2 +
               c[8]*testData[,eval(paste0("x",i))] * testData$x3 + 
               c[9]*testData[,eval(paste0("x",i))] * testData$x4 +
               c[10]*I(testData[,eval(paste0("x",notdv[1]))]^2) +
               c[11]*I(testData[,eval(paste0("x",notdv[2]))]^2) +
               c[12]*I(testData[,eval(paste0("x",notdv[3]))]^2)
      
      
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,3] <- r2
    }
    p[3] <- length(c)
    
    
    # ------------------------------------------------------------
    # ----------------"All_Pair_Int"------------------------------
    # ------------------------------------------------------------
    
    # Linear plus all pair-wise interactions
    for(i in 1:4){
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4 +
                     
                     I(trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3) + 
                     I(trainData$x3 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4))
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  c[4]*testData$x3 +  c[5]*testData$x4 +
        c[6]*testData$x1 * testData$x1 + 
        c[7]*testData$x1 * testData$x2 +
        c[8]*testData$x1 * testData$x3 + 
        c[9]*testData$x1 * testData$x4 +
        c[10]*testData$x2 * testData$x2 +
        c[11]*testData$x2 * testData$x3 + 
        c[12]*testData$x2 * testData$x4 +
        c[13]*testData$x3 * testData$x3 + 
        c[14]*testData$x3 * testData$x4 +
        c[15]*testData$x4 * testData$x4
      
      
      
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,4] <- r2
    }
    p[4] <- length(c)
    
    # ------------------------------------------------------------
    # ----------------"All_Pair_Int_Cub"--------------------------
    # ------------------------------------------------------------
    
    # All linear, quadratic and cubic terms, pairwise interactions
    
    for(i in 1:4){
      notdv <- (1:4)[-i]
      
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4 +
                     
                     I(trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3) + 
                     I(trainData$x3 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4) +
                     I(trainData$x1^3) + 
                     I(trainData$x2^3) + 
                     I(trainData$x3^3) + 
                     I(trainData$x4^3) 
      )
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  c[4]*testData$x3 +  c[5]*testData$x4 +
        c[6]*testData$x1 * testData$x1 + 
        c[7]*testData$x1 * testData$x2 +
        c[8]*testData$x1 * testData$x3 + 
        c[9]*testData$x1 * testData$x4 +
        c[10]*testData$x2 * testData$x2 +
        c[11]*testData$x2 * testData$x3 + 
        c[12]*testData$x2 * testData$x4 +
        c[13]*testData$x3 * testData$x3 + 
        c[14]*testData$x3 * testData$x4 +
        c[15]*testData$x4 * testData$x4 +
        c[16]*I(testData$x1^3) +
        c[17]*I(testData$x2^3) +
        c[18]*I(testData$x3^3) +
        c[19]*I(testData$x4^3) 
      
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,5] <- r2
    }
    p[5] <- length(c)
    
    
    # ------------------------------------------------------------
    # -----------------"Norm_Threew_Int"--------------------------
    # ------------------------------------------------------------
    
    # All linear, quadratic, cubic, and regular three-way interactions
    for(i in 1:4){
      notdv <- (1:4)[-i]
      
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4 +
                     
                     # Pairwise int
                     I(trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3) + 
                     I(trainData$x3 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4) +
                     
                     # Three-way int
                     I(trainData$x1 * trainData$x1 * trainData$x1) + 
            
                     I(trainData$x1 * trainData$x2 * trainData$x3) + 
                     I(trainData$x1 * trainData$x2 * trainData$x4) +
                     I(trainData$x1 * trainData$x3 * trainData$x4) + 
                     
                     I(trainData$x2 * trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3 * trainData$x3) +
                     
                     I(trainData$x4 * trainData$x4 * trainData$x4) 
      )
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  c[4]*testData$x3 +  c[5]*testData$x4 +
        c[6]*testData$x1 * testData$x1 + 
        c[7]*testData$x1 * testData$x2 +
        c[8]*testData$x1 * testData$x3 + 
        c[9]*testData$x1 * testData$x4 +
        c[10]*testData$x2 * testData$x2 +
        c[11]*testData$x2 * testData$x3 + 
        c[12]*testData$x2 * testData$x4 +
        c[13]*testData$x3 * testData$x3 + 
        c[14]*testData$x3 * testData$x4 +
        c[15]*testData$x4 * testData$x4 +
        c[16]*testData$x1 * testData$x1 * testData$x1 +
        c[17]*testData$x1 *testData$x2 * testData$x3 +
        c[18]*testData$x1 *testData$x2 * testData$x4 +
        c[19]*testData$x1 *testData$x3 * testData$x4 +
        
        c[20]*testData$x2 *testData$x2 * testData$x2 +
        c[21]*testData$x2 *testData$x3 * testData$x4 +
        
        c[22]*testData$x3 *testData$x3 * testData$x3 +
        
        c[23]*testData$x4 *testData$x4 * testData$x4
      
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,6] <- r2
    }
    p[6] <- length(c)
    
    
    # ------------------------------------------------------------
    # ---------------- "All_Threew_Int" --------------------------
    # ------------------------------------------------------------
    
    # All linear, quadratic and cubic terms, pairwise and three-way interactions
    for(i in 1:4){
      notdv <- (1:4)[-i]
      
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4 +
                     
                     # Pairwise int
                     I(trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3) + 
                     I(trainData$x3 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4) +
                     
                     # Three-way int
                     I(trainData$x1 * trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x1 * trainData$x4) +
                     I(trainData$x1 * trainData$x2 * trainData$x2) +
                     I(trainData$x1 * trainData$x2 * trainData$x3) + 
                     I(trainData$x1 * trainData$x2 * trainData$x4) +
                     I(trainData$x1 * trainData$x3 * trainData$x3) + 
                     I(trainData$x1 * trainData$x3 * trainData$x4) +
                     I(trainData$x1 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x2 * trainData$x4) +
                     I(trainData$x2 * trainData$x3 * trainData$x3) + 
                     I(trainData$x2 * trainData$x3 * trainData$x4) +
                     I(trainData$x2 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3 * trainData$x3) +
                     I(trainData$x3 * trainData$x3 * trainData$x4) +
                     I(trainData$x3 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4 * trainData$x4) 
      )
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  c[4]*testData$x3 +  c[5]*testData$x4 +
        c[6]*testData$x1 * testData$x1 + 
        c[7]*testData$x1 * testData$x2 +
        c[8]*testData$x1 * testData$x3 + 
        c[9]*testData$x1 * testData$x4 +
        c[10]*testData$x2 * testData$x2 +
        c[11]*testData$x2 * testData$x3 + 
        c[12]*testData$x2 * testData$x4 +
        c[13]*testData$x3 * testData$x3 + 
        c[14]*testData$x3 * testData$x4 +
        c[15]*testData$x4 * testData$x4 +
        c[16]*testData$x1 * testData$x1 * testData$x1 + 
        c[17]*testData$x1 *testData$x1 * testData$x2 +
        c[18]*testData$x1 *testData$x1 * testData$x3 + 
        c[19]*testData$x1 *testData$x1 * testData$x4 +
        c[20]*testData$x1 *testData$x2 * testData$x2 +
        c[21]*testData$x1 *testData$x2 * testData$x3 +
        c[22]*testData$x1 *testData$x2 * testData$x4 +
        c[23]*testData$x1 *testData$x3 * testData$x3 +
        c[24]*testData$x1 *testData$x3 * testData$x4 +
        c[25]*testData$x1 *testData$x4 * testData$x4 +
        
        c[26]*testData$x2 *testData$x2 * testData$x2 +
        c[27]*testData$x2 *testData$x2 * testData$x3 +
        c[28]*testData$x2 *testData$x2 * testData$x4 +
        c[29]*testData$x2 *testData$x3 * testData$x3 +
        c[30]*testData$x2 *testData$x3 * testData$x4 +
        c[31]*testData$x2 *testData$x4 * testData$x4 +
        
        c[32]*testData$x3 *testData$x3 * testData$x3 +
        c[33]*testData$x3 *testData$x3 * testData$x4 +
        c[34]*testData$x3 *testData$x4 * testData$x4 +

        c[35]*testData$x4 *testData$x4 * testData$x4
      
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,7] <- r2
    }
    p[7] <- length(c)
    # Normal four-way Int
    
    
    # ------------------------------------------------------------
    # ---------------- "Norm_Fourw_Int" --------------------------
    # ------------------------------------------------------------
    for(i in 1:4){
      notdv <- (1:4)[-i]
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4 +
                     
                     # Pairwise int
                     I(trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3) + 
                     I(trainData$x3 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4) +
                     
                     # Three-way int
                     I(trainData$x1 * trainData$x1 * trainData$x1) + 
                     
                     I(trainData$x1 * trainData$x2 * trainData$x3) + 
                     I(trainData$x1 * trainData$x2 * trainData$x4) +
                     I(trainData$x1 * trainData$x3 * trainData$x4) + 
                     
                     I(trainData$x2 * trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3 * trainData$x3) +
                     
                     I(trainData$x4 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x1 * trainData$x2 * trainData$x3 * trainData$x4) 
      )
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  c[4]*testData$x3 +  c[5]*testData$x4 +
        c[6]*testData$x1 * testData$x1 + 
        c[7]*testData$x1 * testData$x2 +
        c[8]*testData$x1 * testData$x3 + 
        c[9]*testData$x1 * testData$x4 +
        c[10]*testData$x2 * testData$x2 +
        c[11]*testData$x2 * testData$x3 + 
        c[12]*testData$x2 * testData$x4 +
        c[13]*testData$x3 * testData$x3 + 
        c[14]*testData$x3 * testData$x4 +
        c[15]*testData$x4 * testData$x4 +
        c[16]*testData$x1 * testData$x1 * testData$x1 +
        c[17]*testData$x1 *testData$x2 * testData$x3 +
        c[18]*testData$x1 *testData$x2 * testData$x4 +
        c[19]*testData$x1 *testData$x3 * testData$x4 +
        
        c[20]*testData$x2 *testData$x2 * testData$x2 +
        c[21]*testData$x2 *testData$x3 * testData$x4 +
        
        c[22]*testData$x3 *testData$x3 * testData$x3 +
        
        c[23]*testData$x4 *testData$x4 * testData$x4 +
        c[24]*testData$x1 *testData$x2 * testData$x3 * testData$x4
      
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,8] <- r2
    }
    p[8] <- length(c)
    
    
    
    # ------------------------------------------------------------
    # ---------------- "All_Fourw_Int" --------------------------
    # ------------------------------------------------------------
    # All 2,3, four-way int
    
    for(i in 1:4){
      notdv <- (1:4)[-i]
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData$x1 +  trainData$x2 +  trainData$x3 +  trainData$x4 +
                     
                     # Pairwise int
                     I(trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3) + 
                     I(trainData$x3 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4) +
                     
                     # Three-way int
                     I(trainData$x1 * trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x1 * trainData$x4) +
                     I(trainData$x1 * trainData$x2 * trainData$x2) +
                     I(trainData$x1 * trainData$x2 * trainData$x3) + 
                     I(trainData$x1 * trainData$x2 * trainData$x4) +
                     I(trainData$x1 * trainData$x3 * trainData$x3) + 
                     I(trainData$x1 * trainData$x3 * trainData$x4) +
                     I(trainData$x1 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x2 * trainData$x4) +
                     I(trainData$x2 * trainData$x3 * trainData$x3) + 
                     I(trainData$x2 * trainData$x3 * trainData$x4) +
                     I(trainData$x2 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3 * trainData$x3) +
                     I(trainData$x3 * trainData$x3 * trainData$x4) +
                     I(trainData$x3 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x4 * trainData$x4 * trainData$x4) +
                   
                     # Four way interactions
                     I(trainData$x1 * trainData$x1 * trainData$x1 * trainData$x1) + 
                     I(trainData$x1 * trainData$x1 * trainData$x1 * trainData$x2) +
                     I(trainData$x1 * trainData$x1 * trainData$x1 * trainData$x3) + 
                     I(trainData$x1 * trainData$x1 * trainData$x1 * trainData$x4) +
                     I(trainData$x1 * trainData$x1 * trainData$x2 * trainData$x2) +
                     I(trainData$x1 * trainData$x1 * trainData$x2 * trainData$x3) + 
                     I(trainData$x1 * trainData$x1 * trainData$x2 * trainData$x4) +
                     I(trainData$x1 * trainData$x1 * trainData$x3 * trainData$x3) + 
                     I(trainData$x1 * trainData$x1 * trainData$x3 * trainData$x4) +
                     I(trainData$x1 * trainData$x1 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x1 * trainData$x2 * trainData$x2 * trainData$x2) +
                     I(trainData$x1 * trainData$x2 * trainData$x2 * trainData$x3) + 
                     I(trainData$x1 * trainData$x2 * trainData$x2 * trainData$x4) +
                     I(trainData$x1 * trainData$x2 * trainData$x3 * trainData$x3) + 
                     I(trainData$x1 * trainData$x2 * trainData$x3 * trainData$x4) +
                     I(trainData$x1 * trainData$x2 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x1 * trainData$x3 * trainData$x3 * trainData$x3) +
                     I(trainData$x1 * trainData$x3 * trainData$x3 * trainData$x4) +
                     I(trainData$x1 * trainData$x3 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x1 * trainData$x4 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x2 * trainData$x2 * trainData$x2) +
                     I(trainData$x2 * trainData$x2 * trainData$x2 * trainData$x3) + 
                     I(trainData$x2 * trainData$x2 * trainData$x2 * trainData$x4) +
                     I(trainData$x2 * trainData$x2 * trainData$x3 * trainData$x3) + 
                     I(trainData$x2 * trainData$x2 * trainData$x3 * trainData$x4) +
                     I(trainData$x2 * trainData$x2 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x3 * trainData$x3 * trainData$x3) +
                     I(trainData$x2 * trainData$x3 * trainData$x3 * trainData$x4) +
                     I(trainData$x2 * trainData$x3 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x2 * trainData$x4 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x3 * trainData$x3 * trainData$x3) +
                     I(trainData$x3 * trainData$x3 * trainData$x3 * trainData$x4) +
                     I(trainData$x3 * trainData$x3 * trainData$x4 * trainData$x4) +
                     
                     I(trainData$x3 * trainData$x4 * trainData$x4 * trainData$x4) +
                     I(trainData$x4 * trainData$x4 * trainData$x4 * trainData$x4)
                     
                   )
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData$x1 +  c[3]*testData$x2 +  c[4]*testData$x3 +  c[5]*testData$x4 +
        c[6]*testData$x1 * testData$x1 + 
        c[7]*testData$x1 * testData$x2 +
        c[8]*testData$x1 * testData$x3 + 
        c[9]*testData$x1 * testData$x4 +
        c[10]*testData$x2 * testData$x2 +
        c[11]*testData$x2 * testData$x3 + 
        c[12]*testData$x2 * testData$x4 +
        c[13]*testData$x3 * testData$x3 + 
        c[14]*testData$x3 * testData$x4 +
        c[15]*testData$x4 * testData$x4 +
        c[16]*testData$x1 * testData$x1 * testData$x1 + 
        c[17]*testData$x1 *testData$x1 * testData$x2 +
        c[18]*testData$x1 *testData$x1 * testData$x3 + 
        c[19]*testData$x1 *testData$x1 * testData$x4 +
        c[20]*testData$x1 *testData$x2 * testData$x2 +
        c[21]*testData$x1 *testData$x2 * testData$x3 +
        c[22]*testData$x1 *testData$x2 * testData$x4 +
        c[23]*testData$x1 *testData$x3 * testData$x3 +
        c[24]*testData$x1 *testData$x3 * testData$x4 +
        c[25]*testData$x1 *testData$x4 * testData$x4 +
        
        c[26]*testData$x2 *testData$x2 * testData$x2 +
        c[27]*testData$x2 *testData$x2 * testData$x3 +
        c[28]*testData$x2 *testData$x2 * testData$x4 +
        c[29]*testData$x2 *testData$x3 * testData$x3 +
        c[30]*testData$x2 *testData$x3 * testData$x4 +
        c[31]*testData$x2 *testData$x4 * testData$x4 +
        
        c[32]*testData$x3 *testData$x3 * testData$x3 +
        c[33]*testData$x3 *testData$x3 * testData$x4 +
        c[34]*testData$x3 *testData$x4 * testData$x4 +
        
        c[35]*testData$x4 *testData$x4 * testData$x4 +
        # four way interactions
        c[36]*testData$x1*testData$x1 * testData$x1 * testData$x1 + 
        c[37]*testData$x1*testData$x1 *testData$x1 * testData$x2 +
        c[38]*testData$x1*testData$x1 *testData$x1 * testData$x3 + 
        c[39]*testData$x1*testData$x1 *testData$x1 * testData$x4 +
        c[40]*testData$x1*testData$x1 *testData$x2 * testData$x2 +
        c[41]*testData$x1*testData$x1 *testData$x2 * testData$x3 +
        c[42]*testData$x1*testData$x1 *testData$x2 * testData$x4 +
        c[43]*testData$x1*testData$x1 *testData$x3 * testData$x3 +
        c[44]*testData$x1*testData$x1 *testData$x3 * testData$x4 +
        c[45]*testData$x1*testData$x1 *testData$x4 * testData$x4 +
        
        c[46]*testData$x1 *testData$x2 *testData$x2 * testData$x2 +
        c[47]*testData$x1 *testData$x2 *testData$x2 * testData$x3 +
        c[48]*testData$x1 *testData$x2 *testData$x2 * testData$x4 +
        c[49]*testData$x1 *testData$x2 *testData$x3 * testData$x3 +
        c[50]*testData$x1 *testData$x2 *testData$x3 * testData$x4 +
        c[51]*testData$x1 *testData$x2 *testData$x4 * testData$x4 +
        
        c[52]*testData$x1 *testData$x3 *testData$x3 * testData$x3 +
        c[53]*testData$x1 *testData$x3 *testData$x3 * testData$x4 +
        c[54]*testData$x1 *testData$x3 *testData$x4 * testData$x4 +
        
        c[55]*testData$x1 *testData$x4 *testData$x4 * testData$x4 +
        
        c[56]*testData$x2 *testData$x2 *testData$x2 * testData$x2 +
        c[57]*testData$x2 *testData$x2 *testData$x2 * testData$x3 +
        c[58]*testData$x2 *testData$x2 *testData$x2 * testData$x4 +
        c[59]*testData$x2 *testData$x2 *testData$x3 * testData$x3 +
        c[60]*testData$x2 *testData$x2 *testData$x3 * testData$x4 +
        c[61]*testData$x2 *testData$x2 *testData$x4 * testData$x4 +
        
        c[62]*testData$x2 *testData$x3 *testData$x3 * testData$x3 +
        c[63]*testData$x2 *testData$x3 *testData$x3 * testData$x4 +
        c[64]*testData$x2 *testData$x3 *testData$x4 * testData$x4 +
        
        c[65]*testData$x2 *testData$x4 *testData$x4 * testData$x4 +
        
        c[66]*testData$x3 *testData$x3 *testData$x3 * testData$x3 +
        c[67]*testData$x3 *testData$x3 *testData$x3 * testData$x4 +
        c[68]*testData$x3 *testData$x3 *testData$x4 * testData$x4 +
        
        c[69]*testData$x3 *testData$x4 *testData$x4 * testData$x4 +
        
        c[70]*testData$x4 *testData$x4 *testData$x4 * testData$x4 
        
        
        
        ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
        ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
        r2 <- 1 - ss_err/ss_tot
        r2ar[f,i,9] <- r2
    }
    p[9] <- length(c)
    
    # ------------------------------------------------------------
    # ----------------------- "True" -----------------------------
    # ------------------------------------------------------------
  
    # Fit the true model
    
    for(i in 1:4){
      lm_obj <- lm(trainData[,eval(paste0("dx",i,"dt"))] ~  
                     trainData[,eval(paste0("x",i))]  + 
                     I(trainData[,eval(paste0("x",i))] * trainData$x1) + 
                     I(trainData[,eval(paste0("x",i))] * trainData$x2) +
                     I(trainData[,eval(paste0("x",i))] * trainData$x3) + 
                     I(trainData[,eval(paste0("x",i))] * trainData$x4))
      
      c <- lm_obj$coefficients
      
      preds <- c[1] + c[2]*testData[,eval(paste0("x",i))] +
        c[3]*testData[,eval(paste0("x",i))] * testData$x1 + 
        c[4]*testData[,eval(paste0("x",i))] * testData$x2 +
        c[5]*testData[,eval(paste0("x",i))] * testData$x3 + 
        c[6]*testData[,eval(paste0("x",i))] * testData$x4
      
      ss_err <- sum((testData[,paste0("dx",i,"dt")] - preds)^2)
      ss_tot <- sum(testData[,paste0("dx",i,"dt")]^2)
      r2 <- 1 - ss_err/ss_tot
      r2ar[f,i,10] <- r2
    }
    p[10] <- length(p)
    # setTxtProgressBar(pb, k)
  } # end of CV loop
  
  # Take the mean R^2 of each variable, over all folds, per model
  meanr2 <- t(apply(r2ar,3,colMeans))
  
  # Combine all variable-specific R^2 values into a mean
  meanr2 <- cbind(meanr2,rowMeans(meanr2), p)
  colnames(meanr2) <- c((colnames(D))[1:4], "Mean","p")
  rownames(meanr2) <- c("Main", "Main_Int",
                        "Main_Int_Quad", "All_Pair_Int", 
                        "All_Pair_Int_Cub", "Norm_Threew_Int",
                        "All_Threew_Int", "Norm_Fourw_Int", "All_Fourw_Int",
                        "True")
  
  # Return a table with all model fit results 
  meanr2 
} # eoF


# ----------------------------------------------------------------------
# --------- Get "True Model" parameter estimates  ----------------------
# ----------------------------------------------------------------------

# This function is only implemented for 4 selected models
getDEpars <- function(D, model="All_Pair_Int_Cub", residout = TRUE){
  # define storage of parameters
  summaries <- list()
  resids <- list() 
  sds <- rep(0,4)
  if(model=="All_Fourw_Int") {pars <- matrix(0,4,70)} else{
  pars <- matrix(0,4,28) }
  
  # define feature vector
  if(model=="All_Fourw_Int"){
    feat <- quote(c(1,
                    init$x1 , init$x2 , init$x3 , init$x4 ,
                    
                    # Pairwise int
                    init$x1 * init$x1 ,
                    init$x1 * init$x2 ,
                    init$x1 * init$x3 ,
                    init$x1 * init$x4 ,
                    init$x2 * init$x2 ,
                    init$x2 * init$x3 ,
                    init$x2 * init$x4 ,
                    init$x3 * init$x3 ,
                    init$x3 * init$x4 ,
                    init$x4 * init$x4 ,
                    
                    init$x1 * init$x1 * init$x1 ,
                    init$x1 * init$x1 * init$x2 ,
                    init$x1 * init$x1 * init$x3 ,
                    init$x1 * init$x1 * init$x4 ,
                    init$x1 * init$x2 * init$x2 ,
                    init$x1 * init$x2 * init$x3 ,
                    init$x1 * init$x2 * init$x4 ,
                    init$x1 * init$x3 * init$x3 ,
                    init$x1 * init$x3 * init$x4 ,
                    init$x1 * init$x4 * init$x4 ,
                    
                    init$x2 * init$x2 * init$x2 ,
                    init$x2 * init$x2 * init$x3 ,
                    init$x2 * init$x2 * init$x4 ,
                    init$x2 * init$x3 * init$x3 ,
                    init$x2 * init$x3 * init$x4 ,
                    init$x2 * init$x4 * init$x4 ,
                    
                    init$x3 * init$x3 * init$x3 ,
                    init$x3 * init$x3 * init$x4 ,
                    init$x3 * init$x4 * init$x4 ,
                    init$x4 * init$x4 * init$x4 ,

                    init$x1 * init$x1 * init$x1 * init$x1,
                    init$x1 * init$x1 * init$x1 * init$x2,
                    init$x1 * init$x1 * init$x1 * init$x3,
                    init$x1 * init$x1 * init$x1 * init$x4,
                    init$x1 * init$x1 * init$x2 * init$x2,
                    init$x1 * init$x1 * init$x2 * init$x3,
                    init$x1 * init$x1 * init$x2 * init$x4,
                    init$x1 * init$x1 * init$x3 * init$x3,
                    init$x1 * init$x1 * init$x3 * init$x4,
                    init$x1 * init$x1 * init$x4 * init$x4,
                    init$x1 * init$x2 * init$x2 * init$x2,
                    init$x1 * init$x2 * init$x2 * init$x3,
                    init$x1 * init$x2 * init$x2 * init$x4,
                    init$x1 * init$x2 * init$x3 * init$x3,
                    init$x1 * init$x2 * init$x3 * init$x4,
                    init$x1 * init$x2 * init$x4 * init$x4,
                    init$x1 * init$x3 * init$x3 * init$x3,
                    init$x1 * init$x3 * init$x3 * init$x4,
                    init$x1 * init$x3 * init$x4 * init$x4,
                    init$x1 * init$x4 * init$x4 * init$x4,

                    init$x2 * init$x2 * init$x2 * init$x2,
                    init$x2 * init$x2 * init$x2 * init$x3,
                    init$x2 * init$x2 * init$x2 * init$x4,
                    init$x2 * init$x2 * init$x3 * init$x3,
                    init$x2 * init$x2 * init$x3 * init$x4,
                    init$x2 * init$x2 * init$x4 * init$x4,
                    init$x2 * init$x3 * init$x3 * init$x3,
                    init$x2 * init$x3 * init$x3 * init$x4,
                    init$x2 * init$x3 * init$x4 * init$x4,
                    init$x2 * init$x4 * init$x4 * init$x4,
                    init$x3 * init$x3 * init$x3 * init$x3,
                    init$x3 * init$x3 * init$x3 * init$x4,
                    init$x3 * init$x3 * init$x4 * init$x4,
                    init$x3 * init$x4 * init$x4 * init$x4,
                    init$x4 * init$x4 * init$x4 * init$x4))
  } else {
  feat <- quote(c(1, init$x1, init$x2, init$x3, init$x4,  
                  init$x1*init$x1, init$x1*init$x2, init$x1*init$x3, init$x1*init$x4,
                  init$x2*init$x2, init$x2*init$x3, init$x2*init$x4,
                  init$x3*init$x3, init$x3*init$x4,
                  init$x4*init$x4,
                  init$x1^3,init$x2^3, init$x3^3, init$x4^3,
                  init$x1*init$x2*init$x3, init$x1*init$x2*init$x4,
                  init$x1*init$x3*init$x4, init$x2*init$x3*init$x4,
                  init$x1^4,init$x2^4, init$x3^4, init$x4^4,
                  init$x1*init$x2*init$x3*init$x4))
  
  parnames <- c("a", "x1", "x2", 'x3', "x4",  
                   "x1*x1", "x1*x2", "x1*x3", "x1*x4",
                   "x2*x2", "x2*x3", "x2*x4",
                   "x3*x3", "x3*x4",
                   "x4*x4",
                   "x1^3", "x2^3",  "x3^3",  "x4^3",
                   "x1*x2*x3", "x1*x2*x4",
                   "x1*x3*x4", "x2*x3*x4",
                   "x1^4","x2^4", "x3^4", "x4^4",
                   "x1*x2*x3*x4"
                )
  }

  
  for(i in 1:4){
    if(model=="True"){
      lm_obj <- lm(D[,eval(paste0("dx",i,"dt"))] ~  
                     D[,eval(paste0("x",i))] +
                     I(D[,eval(paste0("x",i))] * D$x1) + 
                     I(D[,eval(paste0("x",i))] * D$x2) +
                     I(D[,eval(paste0("x",i))] * D$x3) + 
                     I(D[,eval(paste0("x",i))] * D$x4))
      tmp <- lm_obj$coefficients
      pars[i,1] <- tmp[1]
      pars[i,1+i] <- tmp[2]
      if(i == 1) pars[i,6:9] <- tmp[3:6]
      if(i == 2) pars[i,c(7,10,11,12)] <- tmp[3:6]
      if(i == 3) pars[i,c(8,11,13,14)] <- tmp[3:6]
      if(i == 4) pars[i,c(9,12,14,15)] <- tmp[3:6]
    }#end if
    
    if(model=="Main_Int"){
      lm_obj <- lm(D[,eval(paste0("dx",i,"dt"))] ~  
                     D$x1 +  D$x2 +  D$x3 +  D$x4 +
                     I(D[,eval(paste0("x",i))] * D$x1) + 
                     I(D[,eval(paste0("x",i))] * D$x2) +
                     I(D[,eval(paste0("x",i))] * D$x3) + 
                     I(D[,eval(paste0("x",i))] * D$x4))
      tmp <- lm_obj$coefficients
      pars[i,1:5] <- tmp[1:5]
      if(i == 1) pars[i,6:9] <- tmp[6:9]
      if(i == 2) pars[i,c(7,10,11,12)] <- tmp[6:9] 
      if(i == 3) pars[i,c(8,11,13,14)] <- tmp[6:9] 
      if(i == 4) pars[i,c(9,12,14,15)] <- tmp[6:9] 
    } #end if
    

    if(model=="All_Pair_Int_Cub"){
      lm_obj <- lm(D[,eval(paste0("dx",i,"dt"))] ~  
                     D$x1 +  D$x2 +  D$x3 +  D$x4 +
                     I(D$x1 * D$x1) +
                     I(D$x1 * D$x2) +  
                     I(D$x1 * D$x3) +  
                     I(D$x1 * D$x4) +
                     
                     I(D$x2 * D$x2) +  
                     I(D$x2 * D$x3) +  
                     I(D$x2 * D$x4) +  
                     
                     I(D$x3 * D$x3) +   
                     I(D$x3 * D$x4) +  
                     
                     I(D$x4 * D$x4) +  
                     I(D$x1^3) +   
                     I(D$x2^3) +   
                     I(D$x3^3) +   
                     I(D$x4^3)   
      )
      pars[i,1:19] <- lm_obj$coefficients
      
    } #end if
    
    if(model == "All_Fourw_Int"){
      lm_obj <- lm(D[,eval(paste0("dx",i,"dt"))] ~  
                     D$x1 +  D$x2 +  D$x3 +  D$x4 +
                     
                     # Pairwise int
                     I(D$x1 * D$x1) + 
                     I(D$x1 * D$x2) +
                     I(D$x1 * D$x3) + 
                     I(D$x1 * D$x4) +
                     
                     I(D$x2 * D$x2) +
                     I(D$x2 * D$x3) + 
                     I(D$x2 * D$x4) +
                     
                     I(D$x3 * D$x3) + 
                     I(D$x3 * D$x4) +
                     
                     I(D$x4 * D$x4) +
                     
                     # Three-way int
                     I(D$x1 * D$x1 * D$x1) + 
                     I(D$x1 * D$x1 * D$x2) +
                     I(D$x1 * D$x1 * D$x3) + 
                     I(D$x1 * D$x1 * D$x4) +
                     I(D$x1 * D$x2 * D$x2) +
                     I(D$x1 * D$x2 * D$x3) + 
                     I(D$x1 * D$x2 * D$x4) +
                     I(D$x1 * D$x3 * D$x3) + 
                     I(D$x1 * D$x3 * D$x4) +
                     I(D$x1 * D$x4 * D$x4) +
                     
                     I(D$x2 * D$x2 * D$x2) +
                     I(D$x2 * D$x2 * D$x3) + 
                     I(D$x2 * D$x2 * D$x4) +
                     I(D$x2 * D$x3 * D$x3) + 
                     I(D$x2 * D$x3 * D$x4) +
                     I(D$x2 * D$x4 * D$x4) +
                     
                     I(D$x3 * D$x3 * D$x3) +
                     I(D$x3 * D$x3 * D$x4) +
                     I(D$x3 * D$x4 * D$x4) +
                     
                     I(D$x4 * D$x4 * D$x4) +
                     
                     # Four way interactions
                     I(D$x1 * D$x1 * D$x1 * D$x1) + 
                     I(D$x1 * D$x1 * D$x1 * D$x2) +
                     I(D$x1 * D$x1 * D$x1 * D$x3) + 
                     I(D$x1 * D$x1 * D$x1 * D$x4) +
                     I(D$x1 * D$x1 * D$x2 * D$x2) +
                     I(D$x1 * D$x1 * D$x2 * D$x3) + 
                     I(D$x1 * D$x1 * D$x2 * D$x4) +
                     I(D$x1 * D$x1 * D$x3 * D$x3) + 
                     I(D$x1 * D$x1 * D$x3 * D$x4) +
                     I(D$x1 * D$x1 * D$x4 * D$x4) +
                     
                     I(D$x1 * D$x2 * D$x2 * D$x2) +
                     I(D$x1 * D$x2 * D$x2 * D$x3) + 
                     I(D$x1 * D$x2 * D$x2 * D$x4) +
                     I(D$x1 * D$x2 * D$x3 * D$x3) + 
                     I(D$x1 * D$x2 * D$x3 * D$x4) +
                     I(D$x1 * D$x2 * D$x4 * D$x4) +
                     
                     I(D$x1 * D$x3 * D$x3 * D$x3) +
                     I(D$x1 * D$x3 * D$x3 * D$x4) +
                     I(D$x1 * D$x3 * D$x4 * D$x4) +
                     
                     I(D$x1 * D$x4 * D$x4 * D$x4) +
                     
                     I(D$x2 * D$x2 * D$x2 * D$x2) +
                     I(D$x2 * D$x2 * D$x2 * D$x3) + 
                     I(D$x2 * D$x2 * D$x2 * D$x4) +
                     I(D$x2 * D$x2 * D$x3 * D$x3) + 
                     I(D$x2 * D$x2 * D$x3 * D$x4) +
                     I(D$x2 * D$x2 * D$x4 * D$x4) +
                     
                     I(D$x2 * D$x3 * D$x3 * D$x3) +
                     I(D$x2 * D$x3 * D$x3 * D$x4) +
                     I(D$x2 * D$x3 * D$x4 * D$x4) +
                     
                     I(D$x2 * D$x4 * D$x4 * D$x4) +
                     
                     I(D$x3 * D$x3 * D$x3 * D$x3) +
                     I(D$x3 * D$x3 * D$x3 * D$x4) +
                     I(D$x3 * D$x3 * D$x4 * D$x4) +
                     
                     I(D$x3 * D$x4 * D$x4 * D$x4) +
                     I(D$x4 * D$x4 * D$x4 * D$x4)
                   )
      pars[i,] <- lm_obj$coefficients
    }
    
    if(residout) resids[[i]] <- lm_obj$residuals
    sds[i] <- sd(lm_obj$residuals)
    summaries[[i]] <- summary(lm_obj)
  }
  if(model!="All_Fourw_Int") colnames(pars) <- parnames
  list(pars = pars,summaries=summaries, sds=sds, feat = feat, resids = resids)
  
}

# ----------------------------------------------------------------------
# --------- Vector Fields: General DE equation function ----------------
# ----------------------------------------------------------------------
# For use with phaseR

DE <- function(t, y, parameters) {
  
  dy <- numeric(2)
  p <- length(parameters)/2
  
  if(p==70){
  feats <- c(1,    y[1] ,  y[1] ,  y[2] ,  y[2] ,
             
             # Pairwise int
             y[1] *  y[1] ,
             y[1] *  y[1] ,
             y[1] *  y[2] ,
             y[1] *  y[2] ,
             y[1] *  y[1] ,
             y[1] *  y[2] ,
             y[1] *  y[2] ,
             y[2] *  y[2] ,
             y[2] *  y[2] ,
             y[2] *  y[2] ,
             
             # Three-way int
             y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[2] ,
             y[1] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] ,
             
             y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[2] ,
             y[1] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] ,
             
             # Four way interactions
             y[1] *  y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             
             y[1] *  y[1] *  y[1] *  y[1] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[1] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[1] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             y[1] *  y[2] *  y[2] *  y[2] ,
             
             y[2] *  y[2] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] *  y[2] ,
             y[2] *  y[2] *  y[2] *  y[2]
  ) 
          
  }  else {
           
  feats <- c(1, y[1], y[1], y[2], y[2],  
            y[1]*y[1], y[1]*y[1], y[1]*y[2], y[1]*y[2],
            y[1]*y[1], y[1]*y[2], y[1]*y[2],
            y[2]*y[2], y[2]*y[2],
            y[2]*y[2],
            y[1]^3, y[1]^3, y[2]^3, y[2]^3,
            y[1]*y[1]*y[2], y[1]*y[1]*y[2],
            y[1]*y[2]*y[2], y[1]*y[2]*y[2],
            y[1]^4, y[1]^4, y[2]^4, y[2]^4,
            y[1]*y[1]*y[2]*y[2]) 
             }
  
  dy[1] <- t(parameters[1:p])%*%feats
  
  dy[2] <- parameters[(p+1):(p*2)]%*%feats
  
  list(dy)
}



# ----------------------------------------------------------------------
# --------- Generate Data from Estimated DE equation -------------------
# ----------------------------------------------------------------------

genData_est <- function(time, 
                        feat, # requires a quote expression
                        coef, 
                        sds,
                        timestep, 
                        initial, # must be a dataframe, with names x1 to x4
                        noise = FALSE) {
  
  p <- 4
  
  init <- data.frame(x1 = initial[1], x2 = initial[2], 
                     x3 = initial[3], x4 = initial[4])
  # Storage
  nIter <- time / timestep
  X <- matrix(NA, nIter, p)
  X[1, ] <- eval(feat)[2:5]
  
  
  # Set progress bar
  pb <- txtProgressBar(min = 2, max=nIter, initial=0, char="-", style = 3)
  
  # Solve with Euler's method
  if(isFALSE(noise)){
  for(j in 2:nIter) {
    for(i in 1:p) {
      X[j, i] <- X[j-1, i] + (coef[i,]%*%eval(feat))*timestep
    }
    init <- list(x1 = X[j,1],x2 = X[j,2],x3 = X[j,3],x4 = X[j,4])
    
    setTxtProgressBar(pb, j)
  }
  }
  
  if(!isFALSE(noise)){
    for(j in 2:nIter) {
      for(i in 1:p) {
        X[j, i] <- X[j-1, i] + (coef[i,]%*%eval(feat) + rnorm(1,sd = sds[i]))*timestep
      }
      init <- list(x1 = X[j,1],x2 = X[j,2],x3 = X[j,3],x4 = X[j,4])
      setTxtProgressBar(pb, j)
    }
  }
  
  # Return
  v_time <- seq(0, time, length = nIter)
  out <- cbind(v_time, X)
  colnames(out) <- c("time", paste0("x",1:p))
  return(out)
  
} # eoF

