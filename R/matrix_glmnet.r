matrix_glmnet <- function(Dat, formula, family = "binomial", alpha = 0,
                          nperm = 1000, plot = FALSE, theme = theme_blackbox){
  
  #formula <- ~ fraction + accession
  #Dat <- Dat
  #family <- "binomial"
  #alpha <- 0
  #nperm <- 100
  #plot <- FALSE
  #theme <- theme_blackbox
  
  if(alpha != 0){
    alpha <- 0
    warning("WARNING: only alpha = 0  supported. Setting alpha to zero",call. = TRUE)
  }
  if(plot){
    plot <- FALSE
    warning("WARNING: plot not fully implemented. Setting plot to FALSE,", call. = TRUE)
  }
  
  formula <- update(formula, ~ . + -1)
  X <- model.matrix(formula, data = Dat$Map)
  
  RES <- NULL
  for(otu in row.names(Dat$Tab)){
    #otu <- row.names(Dat$Tab)[1]
    #otu <- "OTU_16233"
    
    # Get data and perform cross-validationto chose lambda
    Y <- Dat$Tab[otu,]
    m1 <- cv.glmnet(x = X, y = Y, family = family, alpha = alpha, type.measure = "deviance")
    lambda <- m1$lambda.min
    coef <- coef(m1, s = "lambda.min")
    
    # Permute and refit teh model to obtain NULL
    PERMS <- matrix(NA, nrow = nrow(coef), ncol = nperm)
    row.names(PERMS) <- row.names(coef)
    for(i in 1:nperm){
      #i <- 1
      
      Y.perm <- sample(Y)
      m2 <- glmnet(x = X, y = Y.perm, family = family, lambda = lambda, alpha = alpha)
      coef.perm <- coef(m2)
      PERMS[,i] <- coef.perm[,1]
    }
    
    # FOR THE FUTURE NOT SUPPORTED NOW.
    if(plot){
      dat <- melt(PERMS)
      names(dat) <- c("Coefficient","Permutation","Estimate")
      dat$Type <- "Permutation"
      dat <- rbind(dat,
                   data.frame(Coefficient = row.names(coef),
                              Permutation = NA, Estimate = coef[,1],
                              Type = "Estimate"))
      
      p1 <- ggplot(dat,aes(x = Estimate)) +
        facet_wrap(~ Coefficient, scales = "fixed") +
        geom_density(data = subset(dat, Type == "Permutation"), aes(x = Estimate)) +
        geom_vline(data = subset(dat, Type == "Estimate"), aes(xintercept = Estimate), col = "red") +
        theme
      #p1
    }
    
    # Get p-values and update results
    pvals <- rowSums(abs(coef[,1] / PERMS) < 1) / nperm
    Res <- data.frame(Variable = row.names(coef), Taxon = otu, Estimate = coef[,1],
                      p.value = pvals,row.names = NULL, lambda = lambda)
    RES <- rbind(RES,Res)
  }
  
  # Finish
  return(RES)
}



library(AMOR)
library(glmnet)

data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)

Dat <- create_dataset((Rhizo > 0)*1,Rhizo.map,Rhizo.tax)
set.seed(743)
m1 <- matrix_glmnet(Dat = Dat, formula = ~ fraction + accession, nperm = 100)










