################################################################################
#                                                                              #
# Quantitative Biodiversity Function Source Code                               #
# # Written by: Will Shoemaker                                                 #   
#                                                                              #
# Last update: 2015/03/06                                                      #
#                                                                              #
# Notes:                    
# 
#
################################################################################

require("vegan")||install.packages("vegan");require("vegan")

# This function returns the combined species richness and summed area for two sites given as (2 x n) matrix
# formatted as a site-by-species matrix. 

Two.site.richness <- function(x){
  x2.rows <- as.numeric(unlist(rownames(x)))
  x2.matrix <- t(sapply(1:length(x2.rows), function(i) OTUs[x2.rows[i],]))
  dat <- matrix(as.numeric(unlist(x2.matrix)),nrow=nrow(x2.matrix))
  x2.matrix.colsum <- colSums(dat)
  x2.S <- length(x2.matrix.colsum[x2.matrix.colsum > 0])
  x2.diameter.vector <- sapply(1:length(x2.rows), function(i) Ponds$Diameter[x2.rows[i]])
  x2.area <- Diamter.to.Area(x2.diameter.vector)
  x2.area.sum <- sum(x2.area)
  print(x2.area.sum)
  print(x2.S) 
}


# This function returns the observed richness for a single site. 
S.obs <- function(x = "" ){
  rowSums(x > 0) * 1 
}


# This function calculates Good's coverage for a site
Goods <- function(x = ""){
  1 - (sum(x == 1) / rowSums(x))
}

# The function below calculates Chao 1 for a single site

S.chao1 <- function(x = ""){
  S.obs(x) + (sum(x == 1)^2) / (2 * sum(x == 2))
}

# The function below calculates Chao 2 for a site and site-by-species matrix. 
S.chao2 <-function(site = "", SbyS = ""){
  SbyS = as.data.frame(SbyS)
  x = SbyS[site,]
  SbyS.pa <- (SbyS > 0) * 1 
  Q1 = sum(colSums(SbyS.pa) == 1)
  Q2 = sum(colSums(SbyS.pa) == 2)
  S.chao2 = S.obs(x) + (Q1^2)/(2*Q2)
  return(S.chao2)
}

# The function below calculates the Rank Abundancae Curve (RAC) for a single site
RAC <- function(x = ""){
  x = as.vector(x)
  x.ab = x[x>0]
  x.ab.ranked = x.ab[order(x.ab, decreasing = T)]
  return(x.ab.ranked)
}

#THe function below calculates Simpson's Evenness 
SimpE <- function(x = ""){
  x = as.data.frame(x)
  D <- diversity(x, "inv")
  S <- S.obs(x)
  E <- (D)/S
  return(E)
}

# The function below calculates Smith and Wilson's Evenness 
Evar <- function(x){
  x <- as.vector(x[x>0])
  1 - (2/pi)*atan(var(log(x)))
}

# The function below calculates Shannon's Diversity
H <- function(x = ""){
  H = 0
  for (n_i in x){
    p = n_i / sum(x)
    H = H - p*log(p)
  }
  return(H)
}

# The function below calculates Simpson's diversity 

D <- function(x = ""){
  D = 0
  N = sum(x)
  for (n_i in x){
    D = D + (n_i^2) / (N^2)
  }
  return(D)
}

# The function below estimates Beta-diversity between two samples
bet.w <- function(site1 = "", site2 = ""){
  site1 = subset(site1, select = site1 >0)
  site2 = subset(site2, select = site2 > 0)
  gamma = union(colnames(site1), colnames(site2))
  s = length(gamma)
  a.bar = mean(c(specnumber(site1), specnumber(site2)))
  b.w = round(s/a.bar - 1,3)
  return(b.w)
}



