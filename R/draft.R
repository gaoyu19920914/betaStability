setwd("~/GitHub/betaStability/R")
load("~/GitHub/betaStability/betaStability.RData")
#### load required packages ####
required_packages <- c("vegan"
                       ,"usedist"
                       ,"ecodist"
                       ,"glmnet"
                       ,"gdm"
                       ,"BBmisc"
                       # ,"randomForest"
)

# Function to check and install packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Installing package: ", pkg, "\n"))
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(paste0("Package already installed: ", pkg, "\n"))
  }
}

# Install missing packages
cat("Checking and installing required packages...\n")
for (pkg in required_packages) {
  install_if_missing(pkg)
}
# Load the packages to verify installation
cat("\nLoading packages to verify installation...\n")
for (pkg in required_packages) {
  cat(paste0("Loading package: ", pkg, "\n"))
  library(pkg, character.only = TRUE)
}
cat("\nAll required packages are installed and loaded.\n")


#### testable data ####
# vegan: varechem + varespec
# # vegan: dune.env + dune
# # vegan: BCI.env + BCI
# # vegan: mite.env + mite
# rioja::Ponds, SWAP # only 1 variable in SWAP: pH
# HSAUR3::birds, gardenflowers, watervoles
#

data(varechem)
data(varespec)

# clean /check the input data frame
varechem <- na.omit(varechem)

#### 1. betaDiv ####
bray.dist <- vegdist(varespec, method = "bray")

#### 2. envDist ####
vare.envnorm <- BBmisc::normalize(varechem, method = "range", margin = 2)
# TODO: 1. remove BBmisc dependency and add different normalize methods
# vare.envnorm <- scale(varechem)


# hypothesis1 : all environmental variables are of same importance and can be
# normalized to a range of 0~1.
vare.envdist <- dist(vare.envnorm, method = "euclidean")

varedf <- cbind(bray.dist, vare.envdist)
colnames(varedf) <- c("bray.dist", "vare.envdist")

#### 3. relPred ####
# TODO: 3. add tryCatch for out of range 0~1 values
# function: calculate stability based on predicted and measured distances[-1,1]
calculate_stability <- function(predicted.dist, measured.dist){
  if (predicted.dist > measured.dist){
    return((predicted.dist - measured.dist)/(measured.dist -1))
  } else {
    return((measured.dist - predicted.dist)/measured.dist)
  }
}

##### 3.1 Linear Model ####
# predict distance based on environmental distance
site.names <- rownames(varechem)

# initialize dataframe
stability <- data.frame(matrix(NA, nrow = nrow(varechem), ncol =1))
colnames(stability)[1] <- "stability_LM"
rownames(stability) <- site.names


for (n.site in 1:length(site.names)){
  site.name <- site.names[n.site]
  subset.bray.dist <- dist_subset(bray.dist, labels(bray.dist) != site.name)
  subset.vare.envdist <- dist_subset(vare.envdist, labels(bray.dist) != site.name)

  y = unlist(as.list(subset.bray.dist))
  x = unlist(as.list(subset.vare.envdist))

  # linear.model <- lm(bray.dist ~ vare.envdist) # bray.dist = 0.153*vare.envdist + 0.402
  # summary(linear.model)
  this.linear.model <- lm(y ~ x)
  selected.dist <- dist_get(bray.dist, site.name, site.names)
  mean.dist <- mean(selected.dist)
  selected.envdist <- dist_get(vare.envdist, site.name, site.names)
  mean.envdist <- mean(selected.envdist)
  predicted.dist <- predict(this.linear.model,
                            newdata = data.frame(x = mean.envdist))
  stability[n.site,1] <- calculate_stability(predicted.dist, mean.dist)
}


# check the normality: if the points follow the straight line, there's
# qqnorm(linear.model$residuals)
# qqline(linear.model$residuals)
# test/compare different models
# anova(linear.model1, linear.model2)
# AIC(linear.model1)
# AIC(linear.model2)

##### 3.2 MLR multi-variable linear regression ####
stability <- cbind(stability, c(rep(NA, nrow(stability))))
colnames(stability)[2] <- "stability_MLR"

get_nth_row <- function(df, n) {
  # Check if input is a data frame
  if (!is.data.frame(df)) {
    stop("Error: 'df' must be a data frame.")
  }

  # Check if n is numeric and integer-like
  if (!is.numeric(n) || length(n) != 1 || n != as.integer(n)) {
    stop("Error: 'n' must be a single integer.")
  }

  # Convert to integer in case it's numeric but whole number
  n <- as.integer(n)

  # Check bounds
  if (n < 1 || n > nrow(df)) {
    stop("Error: Row number ", n, " is out of bounds. Data frame has ", nrow(df), " rows.")
  }

  # Return the n-th row as a data frame
  return(df[n, , drop = FALSE])
}

get_deltaenv_rows <- function(x, y, env) {
  return (abs(get_nth_row(env, x) - get_nth_row(env, y)))
}

df.ij.deltaenv.beta <- data.frame(matrix(ncol = 2 + ncol(varechem) + 1, nrow = 0))
colnames(df.ij.deltaenv.beta) <- c("i", "j", colnames(varechem), "beta")

for (i in 1:(nrow(varechem)-1)) {
  for (j in (i+1):nrow(varechem)){
    print(i)
    print(j)
    row.to.add <- unlist(c(i,
                           j,
                           get_deltaenv_rows(i, j, varechem),
                           dist_get(bray.dist, i, j)))
    row.to.add <- setNames(row.to.add, colnames(df.ij.deltaenv.beta))
    df.ij.deltaenv.beta <- rbind(df.ij.deltaenv.beta,
                                 row.to.add)

  }
}
colnames(df.ij.deltaenv.beta) <- c("i", "j", colnames(varechem), "beta")
# View(df.ij.deltaenv.beta)
deltaenvnorm <- BBmisc::normalize(df.ij.deltaenv.beta[, !names(df.ij.deltaenv.beta) %in% c("i", "j", "beta")],
                                  method = "range", margin = 2)
df.ij.deltaenv.beta.norm <- cbind(df.ij.deltaenv.beta[, c("i", "j")],
                                  deltaenvnorm,
                                  subset(df.ij.deltaenv.beta, select = c("beta")))


for (n.site in 1:(nrow(varechem))) {
  site.name <- site.names[n.site]
  validatingset <- df.ij.deltaenv.beta.norm[df.ij.deltaenv.beta.norm$i == n.site | df.ij.deltaenv.beta.norm$j == n.site ,]
  trainingset <- df.ij.deltaenv.beta.norm[!(df.ij.deltaenv.beta.norm$i == n.site | df.ij.deltaenv.beta.norm$j == n.site) ,]

  mlm_variables <- paste(names(deltaenvnorm), collapse = " + ")
  mlm_formula <- as.formula(paste("beta ~", mlm_variables))
  this.mlm <- lm(mlm_formula, data = trainingset)

  predicted.dist <- mean(predict(this.mlm, validatingset))
  other.sites <- setdiff(site.names, site.name)
  selected.dist <- dist_get(bray.dist, site.name, other.sites)
  mean.measured.dist <- mean(selected.dist)
  stability[n.site, 2] <- calculate_stability(predicted.dist, mean.measured.dist)
}

##### 3.3 constrained regression (glmnet + NNLS) ####
# 用glmnet(lower.limits=0) 限制非负的话
# 可解释正相关性，但同时损失预测精度。
# (之后可以用 GAM/GDM/RF 提高预测精度)
stability <- cbind(stability, c(rep(NA, nrow(stability))))
colnames(stability)[3] <- "stability_GLM"
for (n.site in 1:(nrow(varechem))) {
  site.name <- site.names[n.site]
  validatingset <- df.ij.deltaenv.beta.norm[df.ij.deltaenv.beta.norm$i == n.site | df.ij.deltaenv.beta.norm$j == n.site ,]
  trainingset <- df.ij.deltaenv.beta.norm[!(df.ij.deltaenv.beta.norm$i == n.site | df.ij.deltaenv.beta.norm$j == n.site) ,]
  glm_predictors <- subset(trainingset, select = names(deltaenvnorm))
  glm_beta <- subset(trainingset, select = c("beta"))
  this.cv.glmnet <- cv.glmnet(as.matrix(glm_predictors), as.matrix(glm_beta), lower.limits = 0)
  best_lambda <- this.cv.glmnet$lambda.min
  coef(this.cv.glmnet, s = best_lambda)
  beta_pred <- predict(this.cv.glmnet,
                       newx = as.matrix(subset(validatingset, select=names(deltaenvnorm))),
                       s = best_lambda)
  # plot(as.matrix(subset(validatingset, select = c("beta"))), beta_pred, col = "blue",
  #      xlab = "Observed Beta Diversity Indices",
  #      ylab = "GLM Predicted Beta Diversity Indices")
  # abline(0, 1, col = "red")
  other.sites <- setdiff(site.names, site.name)
  selected.dist <- dist_get(bray.dist, site.name, other.sites)
  mean.measured.dist <- mean(selected.dist)
  stability[n.site, 3] <- calculate_stability(mean(beta_pred), mean.measured.dist)
}


##### 3.4 gdm nonlinear prediction ####
stability <- cbind(stability, c(rep(NA, nrow(stability))))
colnames(stability)[4] <- "stability_GDM"

site_ids <- rownames(varespec)
env_data <- cbind(site = site_ids,
                  X = 0,
                  Y = 0,
                  varechem)
gdmDissim <- data.frame(site = as.numeric(site_ids),
                        as.matrix(bray.dist),
                        stringsAsFactors = FALSE)
gdm_data <- formatsitepair(
  bioData = gdmDissim,
  bioFormat = 3,
  predData = env_data,
  siteColumn = "site",
  XColumn = "X",
  YColumn = "Y"
)
gdm_model <- gdm(gdm_data, geo = FALSE)
gdm_pred <- predict(gdm_model, data = gdm_data)


convert_gdm_pred_to_matrix <- function(gdm_pred, site_ids) {
  n_sites <- length(site_ids)
  if(n_sites*(n_sites-1)/2 != length(gdm_pred)) {
    stop("The dimensions of variables are not matched!")
  }

  # Create an empty distance matrix
  gdm_dist_matrix <- matrix(0, nrow = n_sites, ncol = n_sites)
  rownames(gdm_dist_matrix) <- site_ids
  colnames(gdm_dist_matrix) <- site_ids

  n <- 1
  for (i in 1:(n_sites-1)){
    for (j in (i+1):n_sites){
      if (n > length(gdm_pred)) stop("length of gdm_pred out of bourder")
      gdm_dist_matrix[i,j] <- gdm_pred[n]
      gdm_dist_matrix[j,i] <- gdm_pred[n]
      n <- n+1
    }
  }

  return(gdm_dist_matrix)
}

gdm_matrix <- convert_gdm_pred_to_matrix(gdm_pred = gdm_pred,
                                         site_ids = site_ids)

## iterate for each site: compare the mean(measured)bray.dist and mean(gdm_pred)
for (n.site in 1:length(site.names)){
  site.name <- site.names[n.site]
  predicted.dist <- mean(gdm_matrix[site.name,])
  other.sites <- setdiff(site.names, site.name)
  selected.dist <- dist_get(bray.dist, site.name, other.sites)
  mean.measured.dist <- mean(selected.dist)
  stability[n.site, 4] <- calculate_stability(predicted.dist, mean.measured.dist)
}

##### 3.5 GAMs (mgcv) ####
library(mgcv)
stability <- cbind(stability, c(rep(NA, nrow(stability))))
colnames(stability)[5] <- "stability_GAM"

site_ids <- rownames(varespec)
y_GAM <- as.matrix(bray.dist)[lower.tri(bray.dist)]

varnames <- colnames(varechem)
x_GAM <- lapply(varnames, function(v) {
  m <- as.matrix(dist(varechem[[v]], method = "manhattan"))
  m[lower.tri(m)]
})
x_GAM <- do.call(cbind, x_GAM)
colnames(x_GAM) <- varnames
data_GAM <- cbind(y = y_GAM, x_GAM)

formula_GAM_str <- paste("y ~",
                         paste("s(", varnames, ")", sep = "", collapse = " + "))
model_GAM <- gam(as.formula(formula_GAM_str),
                 data = as.data.frame(data_GAM),
                 method = "REML")
# view the smooth functions of each parameter in varechem
# plot(model_GAM, pages = 1)

# predict Y and compare with existing Y:
pred_y_GAM <- predict(model_GAM,
                      newdata = as.data.frame(data_GAM),
                      type = "response")
pred_y_GAM_matrix <- convert_gdm_pred_to_matrix(gdm_pred = pred_y_GAM,
                                                site_ids = site_ids)

for (n.site in 1:length(site.names)){
  site.name <- site.names[n.site]
  predicted.dist <- mean(pred_y_GAM_matrix[site.name,])
  other.sites <- setdiff(site.names, site.name)
  selected.dist <- dist_get(bray.dist, site.name, other.sites)
  mean.measured.dist <- mean(selected.dist)
  stability[n.site, 5] <- calculate_stability(predicted.dist, mean.measured.dist)
}

##### 3.6 Machine learning 1. RandomForest ####

rf <- randomForest::randomForest(bray.dist ~ vare.envdist, data = varedf, replace = FALSE, importance = TRUE, proximity = TRUE)
##### 3.7 Machine learning 2. XGBoost ####


##### 3.8 predict community then calculate diversity ####
# MicroEcoTools
# specificity
# microbiomeSeq

