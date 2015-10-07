#' Obtain overparametrized model matrix
overparametrized_model.matrix <- function(formula, data, remove.constant = TRUE,
                                          intercept = TRUE){
  formula <- ~ fraction * accession
  data <- Dat$Map
  intercept <- TRUE
  remove.constant <- TRUE
  
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
              data = Dat$Map)
  
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
