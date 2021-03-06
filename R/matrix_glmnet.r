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
#' @param X A numeric matrix to be used as desing matrix for the model. If missing,
#' the script will use the formula interface with the \link{overparametrized_model.matrix}
#' function to obtain a valid design matrix.
#' @param formula Right hand side of the formula to be used for the model.
#' Variable names must correspond to header names in Map portion of the Dataset object.
#' Ignored if X is passed. Note that offset is not supported.
#' @param family The model family to be used. See \link{glmnet} for more help.
#' @param alpha The elastic net parameter. 1 corresponds ot Lasso and 0 to Ridge,
#' with intermediate values cirrespondig to mixtures between the two.
#' This function only supports 0 for the moment and will use that value regardless
#' of what the user supplies.
#' @param glmnet.intercept Logical indicating whether to add an intercept. See intercept
#' option in \link{glmnet}
#' @param glmnet.offset Either a character indicating the name of the variable to be used
#' as offset; or a numeric vector of length equal to the number of samples to be used
#' directly as offset, see offset option in \link{glmnet}.
#' Note that family multinomial takes a matrix as offset, while passing a vector of
#' variable names should work, note that this has not been tested and might produce
#' incorrect results. The function will print a warning if you try to do this.
#' @param method Either "permuation" or "bootstrap". The default and preferred method is
#' permutation.
#' @param nperm Number of permutations to perform.
#' @param nboot Number of bootstrap pseudo replicates to perform.
#' @param plot Currently not implemented. The funciton contains code for
#' plotting null distributions of parameters and observed values.
#' Might eventually be an independent function.
#' @param theme ggplot2 theme for plots.
#' @param verbose Logical indicating if progress should be printed
#' 
#' @return
#' Returns a list with the following elements
#' \item{X}{Design matrix}
#' \item{coeffcients}{Matrix of coefficients per taxon}
#' \item{lambda}{Vector of lambda values per taxon}
#' \item{method}{Method used for the replicates}
#' \item{family}{Model family}
#' \item{nreps}{Number of replicates performed}
#' \item{call}{Function call}
#' \item{formula}{model formula}
#' \item{samples}{List of coefficiente matrices per OTU and replicate}
#' 
#' @author Sur from Dangl Lab
#' 
#' @examples
#' # Load data
#' data(Rhizo)
#' data(Rhizo.map)
#' 
#' # Create dataset with the first 10 OTUS

#' Dat <- create_dataset(Rhizo[1:10,],Rhizo.map)
#' 
#' # Always set seed before using random numbers (like in permuations)
#' set.seed(743)
#' m1 <- matrix_glmnet(Dat = Dat, formula = ~ fraction + accession,
#'                     nperm = 100,family = "poisson")
#'
#' # Use summary function to obtain statistics
#' m1.sum <- summary(m1)
#' # See distribution of pvalues behaves as expectee
#' hist(m1.sum$coefficients$p.value)
#' 
#' # Use bootstrap instead
#' m1 <- matrix_glmnet(Dat = Dat, formula = ~ fraction + accession,
#'                     nperm = 100,family = "poisson", method ="bootstrap")
#' 
#' # Use summary function to obtain statistics
#' m1.sum <- summary(m1)
#' 
#' # Visualizing coefficientw with ggplot requires a bit of rearrangement because
#' # the percent symbol in the column names causes trouble
#' dat <- subset(m1.sum$coefficients, Taxon == "OTU_14834")
#' dat$lower <- dat[,"2.5%"]
#' dat$upper <- dat[,"97.5%"]
#' p1 <- ggplot(dat,aes(x = lower, y = Variable)) +
#'   geom_segment(aes(xend = upper,yend = Variable)) +
#'   geom_vline(xintercept = 0, col = 'red') +
#'   theme_blackbox
#' p1
matrix_glmnet <- function(Dat, X = NULL, formula = NULL, family = "binomial", alpha = 0,
                          glmnet.intercept = TRUE, glmnet.offset = NULL, method = "permutation",
                          nperm = 1000, nboot = nperm, plot = FALSE,
                          theme = theme_blackbox, verbose = TRUE){
  
  # Checking user parameters
  if(alpha != 0){
    alpha <- 0
    warning("WARNING: only alpha = 0  supported. Setting alpha to zero",call. = TRUE)
  }
  if(plot){
    plot <- FALSE
    warning("WARNING: plot not fully implemented. Setting plot to FALSE,", call. = TRUE)
  }
  
  if(is.null(X) && is.null(formula)){
    stop("ERROR: One of X and formula must be not null\n",call. = TRUE)
  }
  
  if(method == "permutation"){
    sample.replace <- FALSE
  }else if(method == "bootstrap"){
    sample.replace <- TRUE
    nperm <- nboot
  }else{
    stop("ERROR: method not valid. Only permutation or bootsrap allowed",call. = TRUE)
  }
  
  # If we get a character we assume it is the name of a variable in map
  # If it is numeric we check the length and if it matches number of samples
  # we use it as offset.
  # Else we send error.
  # Add a warning for multinomial classes
  if(!is.null(glmnet.offset)){
    if(family == "multinomial"){
      warning("WARNING: Use of offset with multinomial class has not been tested in this function",
              call. = TRUE)
    }
    if(is.character(glmnet.offset)){
      glmnet.offset <- Dat$Map[ ,glmnet.offset ]
    }else if(is.numeric(glmnet.offset) && length(glmnet.offset) == ncol(Dat$Tab)){
      glmnet.offset <- glmnet.offset
    }else{
      stop("ERROR: invalid glmnet.offset", call. = TRUE)
    }
  }
  
  # If not matrix passed, obtain matrix with formula, otherwise check matrix
  if(is.null(X)){
    X <- overparametrized_model.matrix(formula = formula, data = Dat$Map,
                                       remove.constant = TRUE, intercept = FALSE)
  }else if(!(is.matrix(X) && is.numeric(X))){
    stop("ERROR: Invalid desgin matrix",call. = TRUE)
  }
  
  true.coefs <- NULL
  lambda.otu <- NULL
  cat("Estimating lambda...\n")
  for(otu in row.names(Dat$Tab)){
    #otu <- row.names(Dat$Tab)[1]
    
    if(verbose) cat(otu,"\n")
    
    # Get data and perform cross-validationto chose lambda
    Y <- Dat$Tab[otu,]
    m1 <- cv.glmnet(x = X, y = Y, family = family, alpha = alpha, type.measure = "deviance",
                    intercept = glmnet.intercept, offset = glmnet.offset)
    lambda <- m1$lambda.min
    coef <- coef(m1, s = "lambda.min")
    
    names(lambda) <- otu
    colnames(coef) <- otu
    lambda.otu <- c(lambda.otu,lambda)
    true.coefs <- cbind(true.coefs,as.matrix(coef))
  }
  
  REPS <- NULL
  cat("Performing permutations...\n")
  for(i in 1:nperm){
    #i <- 1
    
    # Get permuted index
    Y.perm.index <- sample(ncol(Dat$Tab),replace = sample.replace)
    
    for(otu in row.names(Dat$Tab)){
      #otu <- row.names(Dat$Tab)[1]
      #otu <- "OTU_16233"
      #if(verbose) cat(otu,"\n")
    
      # Get data, permute and refit the model to obtain null
      Y <- Dat$Tab[otu,]
      Y.perm <- Y[ Y.perm.index ]
      m2 <- glmnet(x = X, y = Y.perm, family = family, lambda = lambda.otu[otu],
                   alpha = alpha,intercept = glmnet.intercept, offset = glmnet.offset)
      
      # Store permutation results
      coef.perm <- coef(m2)
      REPS[[otu]] <- cbind(REPS[[otu]],coef.perm[,1])
    }
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
  
  # Format output
  RES <- list(X = X, coefficients = true.coefs, lambda = lambda.otu,
                    method = method, family = family, nreps = nperm,
                    call = match.call(),formula = formula, samples = REPS)
  class(RES) <- "matrix.glmnet"
  return(RES)

}


