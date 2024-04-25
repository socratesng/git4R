## Load packages
library(mice)
library(dplyr)
library(purrr)
# library(sjPlot) # NOT working with "mipo"
library(gtsummary) # NOT working with "nlme"

## Imputation

# Extract data for imputation and analysis
dat_imp <- dat_main_final |> 
  select(
    Diagnostic.3Gp, Diagnostic.2Gp, CRPMGDL, SESCAT2, BMI1Cat2, BMI1Cat3,
    BPSysMean, BPDiasMean, PulseMean, MAPMean, AOMcon, AOMcat.mlr, 
    CTQSum, Pregnancy, Race.2Gp, tsatis, trdiscmf, tresil, trrisks, tdisor
  )

# Check missing patterns
dat_imp |> 
  md.pattern()

# Checking influx and outflux pattern
# Higher influx: depend stronger on the imputation model
# Higher outflux: (potentially) the more powerful predictors
# Close to the diagonal indicates that influx and outflux are balanced
fx <- fluxplot(dat_imp)
fx

# Pseudo-imputation (dry run)
# Obtain imputation methods for variables and imputation matrix
# dat_imp |> 
#   mice(maxit = 0, printFlag = FALSE)
# 
# ini.meth <- ini$method
# ini.meth
# 
# ini.pred <- ini$predictorMatrix
# ini.pred

# Customize imputation method for passive imputation
# ini.meth["MAPMean"] <- "~ I(BPDiasMean + (BPSysMean - BPDiasMean) / 3)"
# ini.meth

# Modify imputation matrix for passive imputation (set elements = 0)
# NOT using passively imputed (transformed) variable to predict original variables
# ini.pred[, c("MAPMean")] <- 0
# ini.pred

# Impute data
# dat_imp |>   
#   mice(
#     m = 40,
#     seed = 12345,
#     maxit = 10,
#     method = ini.meth,
#     predictionMatrix = ini.pred
#   )

# "quickpred" solution
new.pred <- dat_imp |> 
  quickpred()

new.pred
new.pred[c("BPSysMean", "BPDiasMean"), "MAPMean"] <- 0
new.pred["AOMcon", "AOMcat.mlr"] <- 0
new.pred

# imp.dat <- dat_imp |> 
#   mice(
#     m = 40,
#     seed = 12345,
#     method = ini.meth,
#     predictionMatrix = new.pred
#   ) # Still issues in fitting (overfitted?) the model

# Use classification and regression trees (CART) for factors (and all variables)
ini <- dat_imp |> # Dry run 
  mice(maxit = 0, method = "cart", printFlag = FALSE)

ini.meth <- ini$method # Retrieve list of prediction methods
ini.meth["MAPMean"] <- "~ I(BPDiasMean + (BPSysMean - BPDiasMean) / 3)"
# Use single quote ['] instead of double quote ["] to avoid error 
ini.meth["AOMcat.mlr"] <- "~ case_when( 
  AOMcon < 133 ~ 'Early', 
  AOMcon > 168 ~ 'Late', 
  .default = 'Normal'
)"
ini.meth

imp.dat <- dat_imp |> 
  mice( # Test running the imputation algorithm
    m = 5,
    method = ini.meth,
    predictionMatrix = new.pred,
    maxit = 2,
    printFlag = FALSE,
    seed = 123
  ) 

imp.dat <- dat_imp |> 
  futuremice( # Parallelly executed imputation
    m = 40,
    parallelseed = 123,
    n.core = parallel::detectCores() - 2, # Maximum # of logical cores - 2
    packages = "dplyr", # Load "dplyr" package for function "case_when" in parallel R sessions
    method = ini.meth,
    predictionMatrix = new.pred,
    maxit = 40
  ) 

# Visual inspection of the convergence of the algorithm and the range of the values
plot(imp.dat, c("AOMcon", "AOMcat.mlr"))

stripplot(imp.dat, AOMcon ~ .imp, pch = 20, cex = 2) # Observed - Blue, Imputed - Red
densityplot(imp.dat, ~AOMcon)

stripplot(imp.dat, AOMcat.mlr ~ .imp, pch = 20, cex = 2) 


##  Analysis for continuous outcomes

# Run model
M.imp.cont.CDvsND <- lapply(F.cont.CDvsND, function(f) with(imp.dat, lm(formula(f))))

# lapply(M.imp.cont.CDvsND, function(fit){
#   summary(pool(fit))
# })

# Summarize model results in tables (CHIP-AE)
T.imp.CHIPAE.CDvsND.crude <- M.imp.cont.CDvsND[c(1, 3, 5, 7, 9)] |> 
  imap(~{ # Use "imap" to each element and its index
    cat("\n\n") # Correct printing blank lines
    print(paste0("Pooling results from Crude Model ", .y))
    tbl_regression(
      .x,
      label = list(Diagnostic.2Gp ~ "Conduct Disorder"),
      show_single_row = everything(),
      estimate_fun = function(x) style_number(x, digits = 2),
      pvalue_fun = function(x) style_pvalue(x, digits = 2)
    ) |> 
      bold_p() |> 
      modify_header(estimate ~ "**Estimate**", p.value ~ "**P**")
  })

T.imp.CHIPAE.CDvsND.adj <- M.imp.cont.CDvsND[c(2, 4, 6, 8, 10)] |> 
  imap(~{ 
    cat("\n\n")
    print(paste0("Pooling results from Adjusted Model ", .y))
    tbl_regression(
      .x,
      label = list(Diagnostic.2Gp ~ "Conduct Disorder", Race.2Gp ~ "Caucasian", SESCAT2 ~ 'Low SES'),
      show_single_row = everything(),
      estimate_fun = function(x) style_number(x, digits = 2),
      pvalue_fun = function(x) style_pvalue(x, digits = 2)
    ) |> 
      bold_p() |> 
      modify_header(estimate ~ "**Estimate**", p.value ~ "**P**")
  })

T.imp.CHIPAE.CDvsND <- map2(
  T.imp.CHIPAE.CDvsND.crude, T.imp.CHIPAE.CDvsND.adj, ~{
    tbl_merge(tbls = list(.x, .y))
    }) |> 
  tbl_merge(
    tab_spanner = c("Satisfaction",
                    "Discomfort",
                    "Resilience",
                    "Risks",
                    "Disorders")    
  ) |> 
  modify_caption("Associations bewtween CD and CHIP-AE scores")
                                                    
# Summarize model results in tables (Other outcomes)
T.imp.Other.CDvsND.crude <- M.imp.cont.CDvsND[c(11, 13, 15, 17, 19, 21)] |> 
  imap(~{ 
    cat("\n\n")
    print(paste0("Pooling results from Crude Model ", .y))
    tbl_regression(
      .x,
      label = list(Diagnostic.2Gp ~ "Conduct Disorder"),
      show_single_row = everything(),
      estimate_fun = function(x) style_number(x, digits = 2),
      pvalue_fun = function(x) style_pvalue(x, digits = 2)
    ) |> 
      bold_p() |> 
      modify_header(estimate ~ "**Estimate**", p.value ~ "**P**")
  })

T.imp.Other.CDvsND.adj <- M.imp.cont.CDvsND[c(12, 14, 16, 18, 20, 22)] |> 
  imap(~{
    cat("\n\n")
    print(paste0("Pooling results from Adjusted Model ", .y))
    tbl_regression(
      .x,
      label = list(Diagnostic.2Gp ~ "Conduct Disorder", Race.2Gp ~ "Caucasian", SESCAT2 ~ 'Low SES'),
      show_single_row = everything(),
      estimate_fun = function(x) style_number(x, digits = 2),
      pvalue_fun = function(x) style_pvalue(x, digits = 2)
    ) |> 
      bold_p() |> 
      modify_header(estimate ~ "**Estimate**", p.value ~ "**P**")
  })

T.imp.Other.CDvsND <- map2(
  T.imp.Other.CDvsND.crude, T.imp.Other.CDvsND.adj, ~{
    tbl_merge(tbls = list(.x, .y))
  }) |> 
  tbl_merge(
    tab_spanner = c("Mean Dias",
                    "Mean Sys",
                    "Mean MAP",
                    "Mean Pulse",
                    "C-reactive protein",
                    "AOM")    
  ) |> 
  modify_caption("Associations of CD with BP, MAP, Pulse, C-reactive protein before tests and Age of Menarche")


##  Analysis for binary outcomes

# Run model
M.imp.bi.CDvsND <- lapply(F.bi.CDvsND, function(f) with(imp.dat, glm(formula(f), family = binomial)))

# lapply(M.imp.bi.CDvsND, function(fit){
#   summary(pool(fit))
# })

# Summarize model results in tables
T.imp.bi.CDvsND.crude <- M.imp.bi.CDvsND[c(1, 3)] |> 
  imap(~{ 
    cat("\n\n")
    print(paste0("Pooling results from Crude Model ", .y))
    tbl_regression(
      .x,
      label = list(Diagnostic.2Gp ~ "Conduct Disorder"),
      exponentiate = TRUE,
      show_single_row = everything(),
      estimate_fun = function(x) style_number(x, digits = 2),
      pvalue_fun = function(x) style_pvalue(x, digits = 2)
    ) |> 
      bold_p() |> 
      modify_header(p.value ~ "**P**")
  })

T.imp.bi.CDvsND.adj <- M.imp.bi.CDvsND[c(2, 4)] |> 
  imap(~{
    cat("\n\n")
    print(paste0("Pooling results from Adjusted Model ", .y))
    tbl_regression(
      .x,
      label = list(Diagnostic.2Gp ~ "Conduct Disorder", Race.2Gp ~ "Caucasian", SESCAT2 ~ 'Low SES'),
      exponentiate = TRUE,
      show_single_row = everything(),
      estimate_fun = function(x) style_number(x, digits = 2),
      pvalue_fun = function(x) style_pvalue(x, digits = 2)
    ) |> 
      bold_p() |> 
      modify_header(p.value ~ "**P**")
  })

T.imp.bi.CDvsND <- map2(
  T.imp.bi.CDvsND.crude, T.imp.bi.CDvsND.adj, ~{
    tbl_merge(tbls = list(.x, .y))
  }) |> 
  tbl_merge(
    tab_spanner = c("Overweight/Obese",
                    "Pregnancy")    
  ) |> 
  modify_caption("Associations of CD with BMI and Pregnancy")


## Sensitivity analysis (regarding nonignorable missing) by Delta-adjustment (adding offsets)

# The distribution of AOM by Diagnostic group
dat_imp |>
  group_by(Diagnostic.2Gp) |> 
  skimr::skim(AOMcon)

# Specify offsets
delta.AOM <- c(-12, -9, -6, -3, 0, 3, 6, 9, 12) 

# Retrieve post-processing string vector
post.AOM <- ini$post
post.AOM

# Test imputation regarding sensitivity
imp.all.dat.test <- vector("list", length(delta.AOM))

for (i in 1:length(delta.AOM)){
  d <- delta.AOM[i]
  post.AOM.cmd <- paste("imp[[j]][,i] <- imp[[j]][,i] +", d)
  post.AOM["AOMcon"] <- post.AOM.cmd
  imp <- dat_imp |>   
    mice(
      method = ini.meth,
      predictionMatrix = new.pred,
      post = post.AOM,
      maxit = 5,
      printFlag = FALSE,
      seed = i
    ) 
  imp.all.dat.test[[i]] <- imp
}

# Plotting the distribution of AOM in two most extreme scenarios
bwplot(imp.all.dat.test[[1]], AOMcon ~ .imp)
bwplot(imp.all.dat.test[[9]], AOMcon ~ .imp)

densityplot(imp.all.dat.test[[1]], ~ AOMcon)
densityplot(imp.all.dat.test[[9]], ~ AOMcon)

# Real imputation regarding sensitivity
imp.all.dat <- vector("list", length(delta.AOM))

for (i in 1:length(delta.AOM)){
  d <- delta.AOM[i]
  post.AOM.cmd <- paste("imp[[j]][,i] <- imp[[j]][,i] +", d)
  post.AOM["AOMcon"] <- post.AOM.cmd
  imp <- dat_imp |>   
    futuremice(
      m = 40,
      parallelseed = 123,
      n.core = parallel::detectCores() - 2,
      packages = "dplyr",
      method = ini.meth,
      predictionMatrix = new.pred,
      post = post.AOM,
      maxit = 40
    ) 
  imp.all.dat[[i]] <- imp
}
