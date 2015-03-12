\name{matrix_glmnet}
\alias{matrix_glmnet}

\title{glmnet model on a Dataset object.}

\description{Takes a Dataset object and applies a glmnet model on each taxon of the aubundance matrix}

\usage{
matrix_glmnet(Dat, formula, family = "binomial", alpha = 0,
    nperm = 1000, plot = FALSE, theme = theme_blackbox)
}

\arguments{
  \item{Dat}{A Dataset object.}
  \item{formula}{Right hand side of the formula to be used for the model. Variable names must correspond to header names in Map portion of the Dataset object.}
  \item{family}{The model family to be used. See \link{glmnet} for more help.}
  \item{alpha}{The elastic net parameter. 1 corresponds ot Lasso and 0 to Ridge, with intermediate values cirrespondig to mixtures between the two. This function only supports 0 for the moment and will sue that value regardless of what the user supplies.}
  \item{nperm}{Number of permutations to perfomr.}
  \item{plot}{Currently not implemented. The funciton contains code for plotting null distributions of parameters and observed values. Might eventually be an independent function.}
  \item{theme}{ggplot2 theme for plots.}
}

\details{Uses \link{cv.glmnet} with 10-fold cross-validation to pick the optiomal lambda value and then uses \link{glmnet} to fit each permuted dataset.}

\value{a data.frame with the following columns:
  \item{Variable}{Variable name.}
  \item{Taxon}{Taxon for which the model was fit.}
  \item{Estimate}{Coefficient value estimate.}
  \item{p.value}{tow-tailed p-value from permutation.}
  \item{lambda}{lambda value chosen via cross-validation}
}

\author{Sur from Dangl Lab}

\examples{
library(AMOR)
library(glmnet)

data(Rhizo)
data(Rhizo.map)
data(Rhizo.tax)
Dat <- create_dataset((Rhizo > 0)*1,Rhizo.map,Rhizo.tax)

# Always set seed before using random numbers (like in permuations)
set.seed(743)

m1 <- matrix_glmnet(Dat = Dat, formula = ~ fraction + accession,
    nperm = 100)

# See distribution of pvalues behaves as expectee
hist(m1$p.values)

# Look for genotype specific differences
m1[ m1$Variable == "accessionLer" & m1$p.value < 0.1, ]

#plot one of them
p1 <- plotgg_taxon(Dat, taxon = "OTU_2324", x = "accession")
dat <- data.frame(accession = levels(p1$data$accession),
                  value = tapply(p1$data$Abundance,
                      p1$data$accession,mean),
                  Group = "Presence")
dat <- rbind(dat,
             data.frame(accession = levels(p1$data$accession), 
                        value = 1 - dat$value,
                        Group = "Absence"))
ggplot(dat,aes(x = accession, y = value, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = c("black","white")) +
  theme_blackbox

# Now for EC enrichments
m1[ m1$Variable == "fractionE" & m1$p.value < 0.01 & m1$Estimate < -3, ]
p1 <- plotgg_taxon(Dat, taxon = "OTU_16757", x = "fraction")
dat <- data.frame(fraction = levels(p1$data$fraction),
                  value = tapply(p1$data$Abundance,
                      p1$data$fraction,mean),
                  Group = "Presence")
dat <- rbind(dat,
             data.frame(fraction = levels(p1$data$fraction), 
                        value = 1 - dat$value,
                        Group = "Absence"))
ggplot(dat,aes(x = fraction, y = value, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_fill_manual(values = c("black","white")) +
  theme_blackbox

}