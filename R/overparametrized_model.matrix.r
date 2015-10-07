#' Obtain overparametrized design matrix
#' 
#' This functions takes a formula and data frame and creates a design
#' matrix that has a column for every level of every factor.
#' 
#' @param formula Righ hand side of formula for which a design matrix is desired.
#' @param data \code{data.frame} containing the variables used in the formula
#' @param remove.constant Logical indicating whether columns that never vary should be removed.
#' Most commonly these are due to specific combinations of factors that never happen on the
#' dataset.
#' @param intercept Logical indicating whether or not to add an intercept
#' 
#' @author Sur from Dangl lab
#' 
#' @examples
#' data(Rhizo.map)
#' 
#' X <- overparametrized_model.matrix(~ fraction * accession,
#'                                    data = Rhizo.map)
#' image(X)
overparametrized_model.matrix <- function(formula, data, remove.constant = TRUE,
                                          intercept = TRUE){

  # Extract all terms from formula
  f1.terms <- terms(formula,data = data)
  term.labels <- attr(x = f1.terms,which = "term.labels")
  
  # Get model matrix for each term
  X <- lapply(term.labels,
              function(term,data){
                f1 <- paste("~ - 1 + ", term, sep = )
                f1 <- formula(f1)
                X <- model.matrix(f1, data = data)
                return(X)},
              data = data)
  
  # Merge matrix into one
  X <- do.call(what = cbind, args = X)
  
  # Remove invariant columns
  if(remove.constant){
    X <- X[ ,sapply(apply(X,2,unique),length) > 1 ]
  }
  
  # Add intercept
  if(intercept){
    intercept <- matrix(1,nrow = nrow(X),ncol = 1, dimnames = list(NULL,c("(Intercept)")))
    X <- cbind(intercept,X)
  }
  
  return(X)
}
