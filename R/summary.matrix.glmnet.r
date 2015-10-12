#' Sumarize matrix_glmnet fit
#' 
#' This function takes a model fit with \link{matrix_glmnet} and calculates a series of
#' statistics depending on the model fit method.
#' 
#' @param object Object from class \code{matrix.glmnet} as produced by \link{matrix_glmnet}.
#' @param quantile.probs Quantiles to be calculated from bootstrap. Only used if the method
#' element in \code{object} is bootstrap.
#' 
#' @author Sur from Dangl lab
#' 
#' @examples
#' data(Rhizo)
#' data(Rhizo.map)
#' Dat <- create_dataset(Rhizo[1:10,],Rhizo.map)
#' 
#' set.seed(743)
#' m1 <- matrix_glmnet(Dat = Dat, formula = ~ fraction,
#'                    nperm = 100, family = "poisson",
#'                    method = "permutation",
#'                    verbose = FALSE)
#' m1.sum <- summary(m1)
#' m1.sum$coefficients
#' 
summary.matrix.glmnet <- function(object,quantile.probs = c(0.01,0.025,0.5,0.975,0.99)){
  # Format output and calculate p-values
  RES <- melt(data = object$coefficients, varnames = c("Variable","Taxon"),
              value.name = "Estimate")
  if(object$method == "permutation"){
    RES$p.value <- NA
    RES$lambda <- NA
    for(otu in colnames(object$coefficients)){
      pvals <- rowSums(abs(object$coefficients[,otu] / object$samples[[otu]]) <= 1) / object$nreps
      RES$p.value[ RES$Taxon == otu ] <- pvals
      RES$lambda[ RES$Taxon == otu ] <- object$lambda[otu]
    } 
  }else if(object$method == "bootstrap"){
    Res <- NULL
    for(otu in colnames(object$coefficients)){
      quants <- apply(object$samples[[otu]],1,quantile,probs = quantile.probs)
      quants <- as.data.frame(t(quants))
      quants$lambda <- object$lambda[otu]
      
      Res <- rbind(Res,quants)
    }
    RES <- cbind(RES,Res)
  }
  
  RES <- list(coefficients = RES, call = object$call)
  # Finish
  return(RES)
}