library(AMOR)
library(glmnet)

data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)

Dat <- create_dataset((Rhizo > 0)*1,Rhizo.map,Rhizo.tax)
set.seed(743)

formula <- ~ fraction + accession
Dat <- Dat
family <- "binomial"
alpha <- 1
nperm <- 100
plot <- FALSE
theme <- theme_blackbox

formula <- update(formula, ~ . + -1)
X <- model.matrix(formula, data = Dat$Map)

RES <- NULL
for(otu in row.names(Dat$Tab)){
  #otu <- row.names(Dat$Tab)[1]
  
  Y <- Dat$Tab[otu,]
  
  #m1 <- glmnet(x = X, y = Y, family = family, alpha = alpha) 
  m1 <- cv.glmnet(x = X, y = Y, family = family, alpha = alpha)
  lambda <- m1$lambda.min
  coef <- coef(m1, s = "lambda.min")

  PERMS <- matrix(NA, nrow = nrow(coef), ncol = nperm)
  row.names(PERMS) <- row.names(coef)
  for(i in 1:nperm){
    #i <- 1
    
    Y.perm <- sample(Y)
    
    m2 <- glmnet(x = X, y = Y.perm, family = family, lambda = lambda, alpha = alpha)
    coef.perm <- coef(m2)
    PERMS[,i] <- coef.perm[,1]
    
  }
  
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
    #geom_boxplot(aes(x = Coefficient, y = Estimate))
    #p1
  }
  
  pvals <- rowSums(abs(coef[,1] / PERMS) < 1) / nboot
  Res <- data.frame(Variable = row.names(coef), Taxon = otu, Estimate = coef[,1], p.value = pvals,row.names = NULL)
  
  RES <- rbind(RES,Res)
}


