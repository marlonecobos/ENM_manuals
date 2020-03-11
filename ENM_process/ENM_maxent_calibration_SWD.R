enm_calibration_swd <- function(all_occurrences, train_occurrences, test_occurrences, 
                                background_folder, batch, out_dir_models, 
                                reg_multiplier, feature_classes = "all", 
                                max_memory = 1000, args = NULL, maxent_path,
                                selection = "OR_AICc", threshold = 5, 
                                percentage = 50, iterations = 500, 
                                keep_models = TRUE, out_dir_evaluation) {
  
  # Slash
  if(.Platform$OS.type == "unix") {sl <- "/"; dl <- "/"} else {sl <- "\\"; dl <- "\\\\"}
  
  #####
  # candidate models
  
  # Data
  ## Environmental variables sets
  m <- dir(background_folder)
  ms <- paste(gsub("/", dl, paste(getwd(), background_folder, sep = sl)), sl, m, sep = "")
  env <- paste("environmentallayers=", paste("\"", ms, "\"", sep = ""), sep = "")
  m <- gsub(".csv$", "", m)
  
  ## Species occurrences
  oc <- all_occurrences
  samp <- paste("samplesfile=", gsub("/", dl, paste("\"", paste(getwd(), oc, sep = sl),
                                                    "\"", sep = "")), sep = "")
  occ <- train_occurrences
  samp1 <- paste("samplesfile=", gsub("/", dl, paste("\"", paste(getwd(), occ, sep = sl),
                                                     "\"", sep = "")), sep = "")
  
  # Maxent settings
  ## Feature classes combinations
  fea <- feature_classes(feature_classes)
  
  # output directories
  dir.create(out_dir_models)
  out.dir <- gsub("/", dl, paste(getwd(), out_dir_models, sep = sl))
  
  # Getting ram to be used
  ram <- paste("-mx", max_memory, "m", sep = "")
  
  # Fixed commands
  ## Intitial command
  in.comm <- paste("java", ram, paste("-jar", gsub("/", dl, paste0("\"", paste(maxent_path, "maxent.jar", sep = sl), "\""))))
  
  ## Autofeature
  a.fea <- "autofeature=false"
  
  ## Other maxent settings
  fin.com <- "extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=false plots=false pictures=false outputformat=raw warnings=false visible=false redoifexists autorun\n"
  fin.com1 <- "extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=false plots=false pictures=false outputformat=logistic warnings=false visible=false redoifexists autorun\n"
  
  # Final set of calibration models
  ## preparin final arguments
  treg <- length(reg_multiplier); tfea <- length(fea); tenv <- length(env)
  total_comb <- treg * tfea * tenv
  repfea <- total_comb / tfea
  
  reg.mult <- rep(reg_multiplier, each = total_comb / treg)
  fea <- rep(rep(fea, repfea / tenv), each = tenv)
  m <- rep(m, total_comb / tenv)
  
  reg.m <- paste0("betamultiplier=", reg.mult)
  env <- rep(env, total_comb / tenv)
  
  ## creating subdirectories
  subdir <- paste("M", reg.mult, "F", names(fea), m, "all", sep = "_")
  subfol <- paste0("outputdirectory=", paste0("\"", out.dir, sl, subdir, "\""))
  di <- sapply(subdir, function(x) {dir.create(paste0(out.dir, sl, x))})
  
  subdir1 <- paste("M", reg.mult, "F", names(fea), m, "train", sep = "_")
  subfol1 <- paste0("outputdirectory=", paste0("\"", out.dir, sl, subdir1, "\""))
  di <- sapply(subdir1, function(x) {dir.create(paste0(out.dir, sl, x))})
  
  ## writing java code
  allc <- paste(in.comm, env, samp, subfol, reg.m, a.fea, fea, args, fin.com)
  trac <- paste(in.comm, env, samp1, subfol1, reg.m, a.fea, fea, args, fin.com1)
  jmx <- unlist(lapply(1:length(allc), function(x) {c(allc[x], trac[x])}))
  if(.Platform$OS.type == "unix") {
    cat(c("#! /bin/csh\n", jmx), file = paste0(batch, ".sh"))
  } else {
    cat(jmx, file = paste0(batch, ".bat"))
  }
  
  # running models
  message("If asked, RUN as administrator")
  run_maxent(batch, maxent_path)
  
  # candidate model messages
  message(paste0("\nA total of ", total_comb, " candidate models will be created"))
  
  
  #####
  # evaluation
  message("\nStarting evaluation process")
  
  # data
  ## For AICc
  raw_folders <- paste0(out.dir, sl, subdir)
  
  ## For pROC and omission rates
  log_folders <- paste0(out.dir, sl, subdir1)
  
  # evaluation process
  poa <- proc_or_aicc(all_occurrences, train_occurrences, test_occurrences, 
                      raw_folders, log_folders, threshold, percentage, 
                      iterations, keep_models) 
  
  ## Erasing main folder of candidate models if keep_models = FALSE
  if(keep_models == FALSE) {
    unlink(out_dir_models, recursive = T)
    message("All candidate models were deleted")
  }
  
  # summary of results and model selection
  list_res <- summary_calibration(poa, selection)
  
  # writing results
  ## csv files
  message("\nWriting calibration results")
  dir.create(out_dir_evaluation)
  
  name <- paste0(out_dir_evaluation, "/calibration_results.csv")
  name0 <- paste0(out_dir_evaluation, "/calibration_stats.csv")
  name1 <- paste0(out_dir_evaluation, "/selected_models.csv")
  
  write.csv(list_res[[3]], file = name, row.names = FALSE)
  write.csv(list_res[[1]], file = name0, row.names = FALSE)
  write.csv(list_res[[2]], file = name1, row.names = FALSE)
  
  ## plot
  png(paste0(out_dir_evaluation, "/calibration_figure.png"), width = 80, height = 80,
      units = "mm", res = 600)
  par(mar = c(4.5, 4, 0.5, 0.5), cex = 0.58)
  plot_proc_aicc(list_res)
  dev.off()
  
  ## writing the html file
  html_calibration(path = out_dir_evaluation, file_name = "calibration_results")
  
  # finishing evaluation
  message("\nProcess finished")
  message("A folder containing results of model calibration for ", total_comb,
          "\ncandidate models has been written")
  
  message("\nThe folder ", out_dir_evaluation, " contains:")
  message("\t-A html file and its dependencies that summarize all the results")
  message("\t-Two csv files with results and statistics from models calibration")
  if(selection == "OR_AICc"){
    message("\t-And an aditional csv file containing the models selected by OR and AICc\n")
  }
  if(selection == "AICc") {
    message("\t-And  an aditional csv file containing the models selected by AICc\n")
  }
  if(selection == "OR") {
    message("\t-And  an aditional csv file containing the models selected by OR\n")
  }
  
  message(paste0("Check your working directory!!!\t", getwd()))
  
  # returning results
  return(list_res[1:3])
}




# Helpers
run_maxent <- function(batch, maxent_path) {
  if(.Platform$OS.type == "unix") {
    batfile_path <- file.path(getwd(), paste0(batch, ".sh")) 
    r_wd <- getwd()
    setwd(maxent_path) 
    
    system(paste("bash", batfile_path), wait = FALSE)
    
  } else {
    batfile_path <- file.path(getwd(), paste0(batch, ".bat"))
    r_wd <- getwd() # real working directory
    setwd(maxent_path) # change temporally the working directory
    
    system2(batfile_path, wait = FALSE, invisible = FALSE)
  }
  setwd(r_wd)
}

n_par <- function(x) {
  lambdas <- x[1:(length(x) - 4)]
  countNonZeroParams <- function(x) {
    if (strsplit(x, split = ", ")[[1]][2] != "0.0")
      1
  }
  no.params <- sum(unlist(sapply(lambdas, countNonZeroParams)))
  return(no.params)
}

feature_classes <- function(feature_classes = "all") {
  fea <- c("linear=true quadratic=false product=false threshold=false hinge=false",
           "linear=false quadratic=true product=false threshold=false hinge=false",
           "linear=false quadratic=false product=true threshold=false hinge=false",
           "linear=false quadratic=false product=false threshold=true hinge=false",
           "linear=false quadratic=false product=false threshold=false hinge=true",
           "linear=true quadratic=true product=false threshold=false hinge=false",
           "linear=true quadratic=false product=true threshold=false hinge=false",
           "linear=true quadratic=false product=false threshold=true hinge=false",
           "linear=true quadratic=false product=false threshold=false hinge=true",
           "linear=false quadratic=true product=true threshold=false hinge=false",
           "linear=false quadratic=true product=false threshold=true hinge=false",
           "linear=false quadratic=true product=false threshold=false hinge=true",
           "linear=false quadratic=false product=true threshold=true hinge=false",
           "linear=false quadratic=false product=true threshold=false hinge=true",
           "linear=false quadratic=false product=false threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=false hinge=false",
           "linear=true quadratic=true product=false threshold=true hinge=false",
           "linear=true quadratic=true product=false threshold=false hinge=true",
           "linear=true quadratic=false product=true threshold=true hinge=false",
           "linear=true quadratic=false product=true threshold=false hinge=true",
           "linear=false quadratic=true product=true threshold=true hinge=false",
           "linear=false quadratic=true product=true threshold=false hinge=true",
           "linear=false quadratic=true product=false threshold=true hinge=true",
           "linear=false quadratic=false product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=false",
           "linear=true quadratic=true product=true threshold=false hinge=true",
           "linear=true quadratic=true product=false threshold=true hinge=true",
           "linear=true quadratic=false product=true threshold=true hinge=true",
           "linear=true quadratic=true product=true threshold=true hinge=true")
  
  names(fea) <- c("l", "q", "p", "t", "h", "lq", "lp", "lt", "lh", "qp", "qt", "qh",
                  "pt", "ph", "th", "lqp", "lqt", "lqh", "lpt", "lph", "qpt", "qph",
                  "qth", "pth", "lqpt", "lqph", "lqth", "lpth", "lqpth")
  
  if(any(feature_classes %in% c("all", "basic", "no.t.h", "no.h", "no.t"))) {
    if(feature_classes == "all"){fea <- fea} 
    if(feature_classes == "basic"){fea <- fea[c(1, 6, 16, 25, 29)]} 
    if(feature_classes == "no.t.h"){fea <- fea[c(1:3, 6:7, 10, 16)]} 
    if(feature_classes == "no.h"){fea <- fea[c(1:4, 6:8, 10:11, 13, 16:17, 19, 21, 25)]}
    if(feature_classes == "no.t"){fea <- fea[c(1:3, 5:7, 9:10, 12, 14, 16, 18, 20, 22, 26)]}
  }else{
    if (any(feature_classes %in% names(fea))) {
      fea <- fea[feature_classes]
    } else {
     stop("Argument 'feature_classes' is not valid.") 
    }
  }
  return(fea)
}

plot_proc_aicc <- function(summary_calibration) {
  plot(na.omit(summary_calibration[[3]])[, 4] ~ 
         log(na.omit(summary_calibration[[3]])[, 5]),
       xlab = "Natural logarithm of AICc", las = 1, col = "#a6cee3",
       ylab = paste0("Omission rates at ", summary_calibration[[4]], "% threshold value"))
  
  nonsig <- summary_calibration[[3]][!summary_calibration[[3]][, 1] %in% 
                                       summary_calibration[[5]][, 1], ]
  
  points(na.omit(nonsig)[, 4] ~ log(na.omit(nonsig)[, 5]),
         col = "#b2df8a", pch = 19, cex = 1.1)
  
  points(na.omit(summary_calibration[[2]])[, 4] ~ 
           log(na.omit(summary_calibration[[2]])[, 5]), col = "#1f78b4", 
         pch = 17, cex = 1.4)
  
  legend("bottomright", legend = c("Selected models", "Non significant models", "All candidate models"),
         pt.cex = c(1.4, 1.1, 1), pch = c(17, 19, 1), box.col = "white",
         col = c("#1f78b4", "#b2df8a", "#a6cee3"), bg = "white") #, inset = c(0.01, 0)
  box()
}

html_calibration <- function (path = getwd(), file_name) {
  rmdfile <- paste(path, paste0(file_name, ".Rmd"), sep = "/")
  cat("---\ntitle: \"ku_enm: calibration results\"\noutput:\n  html_document:\n      toc: true\n      toc_depth: 4\n---\n\n```{r setup, include=FALSE}\nknitr::opts_chunk$set(echo = TRUE)\n```\n\n<br>\n\n### Brief description of the model calibration and selection process\n\n```{r, echo=FALSE}\nst4 <- read.csv(\"calibration_results.csv\")\nsett <- as.character(st4[,1])\nsetts <- strsplit(sett, split = \"_\")\nrm <- vector()\nfor (i in 1:length(setts)) {\nrm[i] <- setts[[i]][2]\n}\nf.clas <- vector()\nfor (i in 1:length(setts)) {\nf.clas[i] <- setts[[i]][4]\n}\nvar.di <- vector()\nfor (i in 1:length(setts)) {\nvar.di[i] <- paste(setts[[i]][5:length(setts[[i]])], collapse = \"_\")\n}\nrm1 <- paste(unique(rm), collapse = \", \")\nf.clas1 <- paste(unique(f.clas), collapse = \", \")\nvar.di1 <- paste(unique(var.di), collapse = \", \")\npar <- rbind(rm1, f.clas1, var.di1)\n```\n\nThis is the final report of the ku_enm_ceval function implemented in the ku_enm R package.\n\nIn all, `r length(st4[,1])` candidate models, with parameters reflecting all combinations of `r length(unique(rm))` regularization multiplier settings, `r length(unique(f.clas))` feature class combinations, and `r length(unique(var.di))` distinct sets of environmental variables, have been evaluated. Model peformance was evaluated based on statistical significance (Partial_ROC), omission rates (OR), and the Akaike information criterion corrected for small sample sizes (AICc).\n\n```{r par, echo=FALSE}\ncolnames(par) <- \"Parameters\"\nrow.names(par) <- c(\"Regularization multipliers\", \"Feature classes\", \"Sets of predictors\")\nknitr::kable(par, digits=c(0,0), row.names = TRUE, caption = \"Table 1. Parameters of the candidate models.\")\n```\n\n<br>\n\nThe results presented below can be found in the calibration output folder if desired for further analyses.\n\n<br>\n<br>\n\n### Model calibration statistics\n\nIn the following table is information about how many models met the four selection criteria that this function uses.\n\n```{r, echo=FALSE}\nst <- read.csv(\"calibration_stats.csv\")\ncolnames(st) <- c(\"Criteria\",\t\"Number_of_models\")\nknitr::kable(st, digits=c(0,0), caption = \"Table 2. General statistics of models that met distinct criteria.\")\n```\n\n<br>\n<br>\n\n### Models selected according to user-defined criteria\n\nThe following table contains the models selected according to the user's pre-defined criteria.\n\nNote that if the selection argument was \"OR_AICc\", delta AICc values were recalculated only among models meeting the omission rate criterion (*E*).\n\n```{r, echo=FALSE}\nst1 <- read.csv(\"selected_models.csv\")\ncolnames(st1) <- c(\"Model\",\t\"Mean_AUC_ratio\",\t\"Partial_ROC\", gsub(\"[.]\", \"%\", colnames(st1)[4]), \"AICc\",\t\"delta_AICc\",\t\"W_AICc\",\t\"num_parameters\")\nknitr::kable(st1, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 3. Performance statistics for the models selected based on the user's pre-defined critera.\")\n```\n\n<br>\n<br>\n\n### Model performance plot\n\nThe figure below shows the position of the selected models in the distribution of all candidate models in terms of statistical significance, omission rates, and AICc values.\n\n![Figure 1. Distribution of all models, non-statistically significant models, and selected models in terms of omission rates and AICc values.](calibration_figure.png){width=60%}\n\n<br>\n<br>\n\n### Performance statistics for all models\n\nFollowing are the performance statistics for all candidate models (a sample if more than 500 models). See file calibration_results.csv for the complete list.\n\n```{r, echo=FALSE}\nst4 <- read.csv(\"calibration_results.csv\")\nif (dim(st4)[1] > 500) {\n   st4 <- st4[1:500, ]\n}\ncolnames(st4) <-  c(\"Model\",\t\"Mean_AUC_ratio\",\t\"Partial_ROC\", gsub(\"[.]\", \"%\", colnames(st4)[4]), \"AICc\",\t\"delta_AICc\",\t\"W_AICc\",\t\"num_parameters\")\nknitr::kable(st4, digits=c(0,3,3,3,3,3,3,0), caption = \"Table 4. Performance statistics for all candidate models.\")\n```",
      file = rmdfile)
  rmarkdown::render(rmdfile, "html_document", quiet = TRUE)
  unlink(rmdfile)
}

summary_calibration <- function(proc_or_aicc_results, selection = "OR_AICc") {
  
  ku_enm_eval <- proc_or_aicc_results
  threshold <- gsub("Omission_rate_at_", "", colnames(ku_enm_eval)[4])
  threshold <- as.numeric(gsub("%", "", threshold))
  
  # Choosing the best models
  if(selection == "OR_AICc") {
    ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
    ku_enm_best <- na.omit(ku_enm_bes[which(ku_enm_bes[, 4] <= threshold / 100), ])
    if(length(ku_enm_best[, 4]) != 0) {
      ku_enm_best[, 6] <- ku_enm_best[, 5] - min(ku_enm_best[, 5], na.rm = TRUE)
      ku_enm_best[, 7] <- exp(-0.5 * ku_enm_best[, 6]) / 
        sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE)
      ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
      ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
      
    }else {
      message("None of the significant candidate models met the omission rate criterion,",
              "\nmodels with the lowest omission rate and lowest AICc will be presented")
      
      ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ]
      ku_enm_best[, 6] <- ku_enm_best[, 5] - min(ku_enm_best[, 5], na.rm = TRUE)
      ku_enm_best[, 7] <- exp(-0.5 * ku_enm_best[, 6]) / 
        sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE)
      ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
      ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
    }
  }
  
  if(selection == "AICc") {
    ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
    ku_enm_best <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]
    if(length(ku_enm_best[, 6]) != 0) {
      ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
    }else {
      message("None of the significant candidate models met the AICc criterion,",
              "\ndelta AICc will be recalculated for significant models")
      
      ku_enm_best[, 6] <- ku_enm_best[, 5] - min(ku_enm_best[, 5], na.rm = TRUE)
      ku_enm_best[, 7] <- exp(-0.5 * ku_enm_best[, 6]) / 
        sum(exp(-0.5 * ku_enm_best[, 6]), na.rm = TRUE)
      ku_enm_best <- ku_enm_best[ku_enm_best[, 6] <= 2, ]
      ku_enm_best <- ku_enm_best[order(ku_enm_best[, 6]), ]
    }
  }
  
  if(selection == "OR") {
    ku_enm_b <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
    ku_enm_bes <- na.omit(ku_enm_eval[ku_enm_eval[, 3] <= 0.05, ])
    ku_enm_bes1 <- ku_enm_b[ku_enm_b[, 3] <= 0.05, ]
    ku_enm_best <- ku_enm_bes1[ku_enm_bes1[, 4] <= threshold / 100, ]
    if(length(ku_enm_best[, 4]) != 0) {
      if(length(ku_enm_best[, 4]) > 10) {
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 4]), ][1:10, ]
      }else {
        ku_enm_best <- ku_enm_best[order(ku_enm_best[, 4]), ]
      }
    }else {
      message("None of the significant candidate models met the omission rate criterion,",
              "\nmodels with the lowest omission rate will be presented")
      
      ku_enm_best <- ku_enm_bes[ku_enm_bes[, 4] == min(ku_enm_bes[, 4]), ][1:10, ]
    }
  }
  
  #####
  #Statistics of the process
  ##Counting
  ku_enm_sign <- ku_enm_eval[!is.na(ku_enm_eval[, 3]), ]
  ku_enm_sign <- ku_enm_sign[ku_enm_sign[, 3] <= 0.05, ]
  
  ku_enm_or <- ku_enm_eval[ku_enm_eval[, 4] <= threshold / 100, ]
  
  ku_enm_AICc <- ku_enm_eval[!is.na(ku_enm_eval[, 6]), ]
  ku_enm_AICc <- ku_enm_AICc[ku_enm_AICc[, 6] <= 2, ]
  
  ku_enm_best_OR <- ku_enm_sign[ku_enm_sign[, 4] <= threshold / 100, ]
  
  ku_enm_best_AICc <- ku_enm_bes[ku_enm_bes[, 6] <= 2, ]
  
  ku_enm_best_OR_AICc <- ku_enm_bes[ku_enm_bes[, 4] <= threshold / 100, ]
  if(length(ku_enm_best_OR_AICc[, 4]) != 0) {
    ku_enm_best_OR_AICc[, 6] <- ku_enm_best_OR_AICc[, 5] - 
      min(ku_enm_best_OR_AICc[, 5], na.rm = TRUE)
    ku_enm_best_OR_AICc[, 7] <- exp(-0.5 * ku_enm_best_OR_AICc[, 6]) / 
      sum(exp(-0.5 * ku_enm_best_OR_AICc[, 6]), na.rm = TRUE)
    ku_enm_best_OR_AICc <- ku_enm_best_OR_AICc[ku_enm_best_OR_AICc[, 6] <= 2, ]
  }
  
  # Preparing the table
  r_names <- c("All candidate models", "Statistically significant models", 
               "Models meeting omission rate criteria",
               "Models meeting AICc criteria", 
               "Statistically significant models meeting omission rate criteria",
               "Statistically significant models meeting AICc criteria",
               "Statistically significant models meeting omission rate and AICc criteria")
  statis <- c(length(ku_enm_eval[, 1]),
              length(ku_enm_sign[, 3]),
              length(ku_enm_or[, 4]),
              length(ku_enm_AICc[, 6]),
              length(ku_enm_best_OR[, 4]),
              length(ku_enm_best_AICc[, 6]),
              length(ku_enm_best_OR_AICc[, 2]))
  
  ku_enm_stats <- data.frame(r_names, statis)
  colnames(ku_enm_stats) <- c("Criteria", "Number of models")
  
  
  
  # returning results
  results <- list(calibration_statistics = ku_enm_stats, 
                  selected_models = ku_enm_best, 
                  calibration_results = ku_enm_eval,
                  threshold = threshold, significant_models = ku_enm_sign)
  return(results)
}

aicc <- function(occurrences, prediction, n_parameters) {
  AIC.valid <- n_parameters < nrow(occurrences)
  if (nrow(prediction) == 0) {
    res <- data.frame(cbind(AICc = NA, delta_AICc = NA,
                            weight_AICc = NA, parameters = n_parameters))
    warning("Cannot calculate AICc when prediction has 0 rows.")
  } else {
    vals <- prediction[paste(prediction[, 1], prediction[, 2]) %in%
                         paste(occurrences[, 1], occurrences[, 2]), 3]
    vals <- na.omit(vals)
    probsum <- sum(prediction[, 3], na.rm = TRUE)
    LL <- colSums(log(t(t(vals)/probsum)), na.rm = TRUE)
    AICc <- ((2 * n_parameters) - (2 * LL)) + (2 * n_parameters * (n_parameters + 1) / 
                                                 (nrow(occurrences) - n_parameters - 1))
    AICc[AIC.valid == FALSE] <- NA
    AICc[is.infinite(AICc)] <- NA
    if (sum(is.na(AICc)) == length(AICc)) {
      warning("More parameters than occurrences, returning NA.")
      res <- data.frame(cbind(AICc, delta_AICc = NA, weight_AICc = NA,
                              parameters = n_parameters))
    } else {
      delta_AICc <- AICc - min(AICc, na.rm = TRUE)
      weight_AICc <- exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc), na.rm = TRUE)
      res <- data.frame(AICc, delta_AICc, weight_AICc, parameters = n_parameters)
      rownames(res) <- NULL
    }
  }
  rownames(res) <- NULL
  return(res)
}

or <- function(prediction, train_occurrences, test_occurrences, threshold) {
  if (min(prediction, na.rm = T) == max(prediction, na.rm = T)) {
    warning("Model imput has no variability, omission rate = NA.")
    om_rate <- NA
  } else {
    vals <- prediction[paste(prediction[, 1], prediction[, 2]) %in%
                         paste(train_occurrences[, 1], train_occurrences[, 2]), 3]
    tvals <- prediction[paste(prediction[, 1], prediction[, 2]) %in%
                          paste(test_occurrences[, 1], test_occurrences[, 2]), 3]
    
    vals <- na.omit(vals); tvals <- na.omit(tvals)
    
    om_rate <- vector("numeric", length = length(threshold))
    for (i in 1:length(threshold)) {
      val <- ceiling(nrow(train_occurrences) * threshold[i] / 100) + 1
      omi_val_suit <- sort(vals)[val]
      om_rate[i] <- length(tvals[tvals < omi_val_suit]) / length(tvals)
    }
    names(om_rate) <- paste("om_rate_", threshold, "%", sep = "")
  }
  return(om_rate)
}

wait_written_done <- function(file) {
  while (!file.exists(file)) {
    if (file.exists(file)) {break()}
    Sys.sleep(1)
  }
  stime <- Sys.time()
  fi <- file.info(file)
  di <- as.numeric(fi$mtime - stime)
  
  while (di == 0) {
    stime <- Sys.time()
    fi <- file.info(file)
    di <- as.numeric(fi$mtime - stime)
    
    if (di < 0) {break()}
    Sys.sleep(1)
  }
  return(di < 0)
}

proc_or_aicc <- function(all_occurrences, train_occurrences, test_occurrences,
                         raw_folders, log_folders, threshold = 5, 
                         percentage = 50, iterations = 500, keep_models = TRUE) {
  #pROCs, omission rates, and AICcs calculation
  message("Evaluation using partial ROC, omission rates, and AICc")
  
  # Slash
  if(.Platform$OS.type == "unix") {sl <- "/"; dl <- "/"} else {sl <- "\\"; dl <- "\\\\"}
  
  # model names
  model_names <- gsub(paste0(".*", dl), "", gsub("_all$", "", raw_folders))
  
  # occurrences
  oc <- read.csv(all_occurrences) 
  spn <- as.character(oc[1, 1])
  oc <- oc[, -1]
  occ <- read.csv(train_occurrences)[, -1] 
  occ1 <- read.csv(test_occurrences)[, -1] 
  
  longitude <- colnames(oc)[1]
  latitude <- colnames(oc)[2]
  
  aics <- list() 
  proc_res <- list() 
  om_rates <- numeric() 
  nm <- length(raw_folders)
  
  if(.Platform$OS.type == "unix") {
    pb <- txtProgressBar(min = 0, max = nm, style = 3) 
  } else {
    pb <- winProgressBar(title = "Progress bar", min = 0, max = nm, width = 300) 
  }

  for(i in 1:nm) {
    Sys.sleep(0.1)
    if(.Platform$OS.type == "unix") {
      setTxtProgressBar(pb, i)
    } else {
      setWinProgressBar(pb, i, title = paste(round(i / nm * 100, 2),
                                             "% of the evaluation has finished"))
    }
    
    #AICc calculation
    lbds <- paste0(raw_folders[i], sl, spn, ".lambdas")
    waiting <- wait_written_done(lbds)
    lambdas <- readLines(lbds)
    
    par_num <- n_par(lambdas)
    
    mods <- paste0(raw_folders[i], sl, spn, ".csv")
    waiting <- wait_written_done(mods)
    mod <- read.csv(mods)
    
    aic <- aicc(oc, mod, par_num)
    aics[[i]] <- aic
    
    #pROCs and omission rates calculation
    mods1 <- paste0(log_folders[i], sl, spn, ".csv")
    waiting <- wait_written_done(mods1)
    mod1 <- read.csv(mods1)
    
    tval <- mod1[paste(mod1[, 1], mod1[, 2]) %in% paste(occ1[, 1], occ1[, 2]), 3]
    proc <- ellipsenm::partial_roc(mod1[, 3], tval, longitude, latitude, 
                                   threshold, iterations, percentage)
    
    proc_res[[i]] <- proc[[1]]
    
    om_rates[i] <- or(mod1, occ, occ1, threshold)
    
    #Erasing calibration models after evaluating them if keep_models = FALSE
    if(keep_models == FALSE) {
      unlink(raw_folders[i], recursive = T)
      unlink(log_folders[i], recursive = T)
    }
  }
  
  if(.Platform$OS.type != "unix") {suppressMessages(close(pb))}

  # From AICc analyses few calculations
  aiccs <- do.call(rbind, aics)
  
  aiccs[, 2] <- aiccs[, 1] - min(aiccs[, 1], na.rm = TRUE)
  aiccs[, 3] <- exp(-0.5 * aiccs[, 2]) / sum(exp(-0.5 * aiccs[, 2]), na.rm = TRUE)
  
  # From pROC analyses
  proc_res_m <- data.frame(model_names, do.call(rbind, proc_res))[, 1:3] 
  
  # Joining the results
  ku_enm_eval <- data.frame(proc_res_m, om_rates, aiccs)
  colnames(ku_enm_eval) <- c("Model", "Mean_AUC_ratio", "pval_pROC",
                             paste0("Omission_rate_at_", threshold, "%"), "AICc",
                             "delta_AICc", "W_AICc", "N_parameters")
  
  return(ku_enm_eval)
}