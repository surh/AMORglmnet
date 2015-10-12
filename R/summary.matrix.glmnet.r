#' Sumarize matrix_glmnet fit
#' 
#' This function takes a model fit with \link{matrix_glmnet} and calculates a series of
#' statistics depending on the model fit method.
#' 
#' @param object Object from class \code{matrix.glmnet} as produced by \link{matrix_glmnet}.
#' @param quantile.probs Quantiles to be calculated from bootstrap. Only used if the method
#' element in \code{object} is bootstrap.
#' @param funlist A list of function to be tested on the list of coefficients. Each function
#' must take a single numeric vector, and return a single numeric value. Preferentially, the
#' list should be named, otherwise functions will be identified by they numeric indices.
#' 
#' @author Sur from Dangl lab
#' 
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo[1:10,],Rhizo.map)
#' 
#' # Fit model and permute 100 times to estimate null.
#' # All the code below should also work if method is changed to "bootstrap"
#' set.seed(743)
#' m1 <- matrix_glmnet(Dat = Dat, formula = ~ fraction,
#'                    nperm = 100, family = "poisson",
#'                    method = "permutation",
#'                    verbose = FALSE)
#' 
#' # Performs permutation test
#' m1.sum <- summary(m1)
#' m1.sum$coefficients
#' 
#' # We can add a function that tests wheather the mean if R and E
#' # samples is different than the soil
#' mytestfun <- function(x) {x[2] - mean(x[3],x[4])}
#' 
#' # Following command repeats the permutation test, and tests also
#' # wheter the value of the new function is significant
#' m1.sum <- summary(m1,funlist = list(mytest = mytestfun))
#' 
#' # The results are stored in the functions element
#' m1.sum$functions
#' 
summary.matrix.glmnet <- function(object,quantile.probs = c(0.01,0.025,0.5,0.975,0.99),
                                  funlist = NULL){

  if(!is.null(funlist)){
    funlist <- sapply(funlist,match.fun)
    if(is.null(names(funlist))){
      names(funlist) <- as.character(1:length(funlist))
    }
  }
  
  # Format output and calculate p-values
  RES <- melt(data = object$coefficients, varnames = c("Variable","Taxon"),
              value.name = "Estimate")
  FUNTEST <- NULL
  if(object$method == "permutation"){
    RES$p.value <- NA
    RES$lambda <- NA
    for(otu in colnames(object$coefficients)){
      pvals <- rowSums(abs(object$coefficients[,otu] / object$samples[[otu]]) <= 1) / object$nreps
      RES$p.value[ RES$Taxon == otu ] <- pvals
      RES$lambda[ RES$Taxon == otu ] <- object$lambda[otu]
      
      for(fun.name in names(funlist)){
        #fun.name <- names(funlist)[1]
        
        FUN <- funlist[[fun.name]]
        
        funres <- apply(object$samples[[otu]],2,FUN)
        estimate <- FUN(object$coefficients[,otu])
        pval <- sum(abs(estimate / funres) <= 1) / object$nreps
        
        funtest <- data.frame(FUN = fun.name,Taxon = otu, Estimate = estimate, p.value = pval)
        FUNTEST <- rbind(FUNTEST,funtest)
      }
    } 
  }else if(object$method == "bootstrap"){
    Res <- NULL
    for(otu in colnames(object$coefficients)){
      quants <- apply(object$samples[[otu]],1,quantile,probs = quantile.probs)
      quants <- as.data.frame(t(quants))
      quants$lambda <- object$lambda[otu]
      
      Res <- rbind(Res,quants)
      
      for(fun.name in names(funlist)){
        #fun.name <- names(funlist)[1]
        
        FUN <- funlist[[fun.name]]
        
        funres <- apply(object$samples[[otu]],2,FUN)
        estimate <- FUN(object$coefficients[,otu])
        quants <- apply(object$samples[[otu]],1,quantile,probs = quantile.probs)
        quants <- as.data.frame(t(quants))
        
        funtest <- data.frame(FUN = fun.name,Taxon = otu, Estimate = estimate)
        funtest <- cbind(funtest,quants)
        FUNTEST <- rbind(FUNTEST,funtest)
      }
      
    }
    RES <- cbind(RES,Res)
  }
  
  row.names(RES) <- NULL
  row.names(FUNTEST) <- NULL
  
  RES <- list(coefficients = RES, functions = FUNTEST, call = object$call)
  # Finish
  return(RES)
}