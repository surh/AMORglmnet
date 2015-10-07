#' Fit  a glmnet model on Dataset object.
#' 
#' This function takes a Dataset object and fits the specified
#' glmnet model on each taxon. righn now only alpha = 0 (ie. ridge regularization)
#' is supported. The function uses 10-fold cross-validation to select lambda, and
#' permutation to obtain p-values.
#' 
#' Uses \link{cv.glmnet} with 10-fold cross-validation to pick the optiomal
#' lambda value and then uses \link{glmnet} to fit each permuted dataset.
#' 
#' @aliases matrix_glmnet
#' 
#' @param Dat A Dataset object.
#' @param formula Right hand side of the formula to be used for the model.
#' Variable names must correspond to header names in Map portion of the Dataset object.
#' @param family The model family to be used. See \link{glmnet} for more help.
#' @param alpha The elastic net parameter. 1 corresponds ot Lasso and 0 to Ridge,
#' with intermediate values cirrespondig to mixtures between the two.
#' This function only supports 0 for the moment and will use that value regardless
#' of what the user supplies.
#' @param nperm Number of permutations to perfomr.
#' @param plot Currently not implemented. The funciton contains code for
#' plotting null distributions of parameters and observed values.
#' Might eventually be an independent function.
#' @param theme ggplot2 theme for plots.
#' 
#' @return
#' \item{Variable}{Variable name.}
#' \item{Taxon}{Taxon for which the model was fit.}
#' \item{Estimate}{Coefficient value estimate.}
#' \item{p.value}{tow-tailed p-value from permutation.}
#' \item{lambda}{lambda value chosen via cross-validation}
#' 
#' @author Sur from Dangl Lab
#' 
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' data(Rhizo.tax)
#' Dat <- create_dataset((Rhizo > 0)*1,Rhizo.map,Rhizo.tax)
#' 
#' # Always set seed before using random numbers (like in permuations)
#' set.seed(743)
#' m1 <- matrix_glmnet(Dat = Dat, formula = ~ fraction + accession,
#'                     nperm = 100)
#'
#' # See distribution of pvalues behaves as expectee
#' hist(m1$p.values)
#' 
#' 
#' # Look for genotype specific differences
#' m1[ m1$Variable == "accessionLer" & m1$p.value < 0.1, ]
#' 
#' #plot one of them
#' p1 <- plotgg_taxon(Dat, taxon = "OTU_2324", x = "accession")
#' dat <- data.frame(accession = levels(p1$data$accession),
#'                   value = tapply(p1$data$Abundance,
#'                   p1$data$accession,mean),
#'                   Group = "Presence")
#' 
#' dat <- rbind(dat,
#'              data.frame(accession = levels(p1$data$accession), 
#'                         value = 1 - dat$value,
#'                         Group = "Absence"))
#' ggplot(dat,aes(x = accession, y = value, fill = Group)) +
#'   geom_bar(stat = "identity", position = "fill") + 
#'   scale_fill_manual(values = c("black","white")) +
#'   theme_blackbox
#' 
#' 
#' # Now for EC enrichments
#' m1[ m1$Variable == "fractionE" & m1$p.value < 0.01 & m1$Estimate < -3, ]
#' p1 <- plotgg_taxon(Dat, taxon = "OTU_16757", x = "fraction")
#' dat <- data.frame(fraction = levels(p1$data$fraction),
#'                   value = tapply(p1$data$Abundance,
#'                   p1$data$fraction,mean),
#'                   Group = "Presence")
#' dat <- rbind(dat,
#'              data.frame(fraction = levels(p1$data$fraction),
#'                         value = 1 - dat$value,
#'                         Group = "Absence"))
#' ggplot(dat,aes(x = fraction, y = value, fill = Group)) +
#'   geom_bar(stat = "identity", position = "fill") +
#'   scale_fill_manual(values = c("black","white")) +
#'   theme_blackbox
matrix_glmnet <- function(Dat, formula, family = "binomial", alpha = 0,
                          nperm = 1000, plot = FALSE, theme = theme_blackbox){
  
#   formula <- ~ FRACTION
#   Dat <- Dat.bin
#   family <- "binomial"
#   alpha <- 0
#   nperm <- 10
#   plot <- FALSE
#   theme <- theme_blackbox
  
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
    cat(otu,"\n")
    
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


