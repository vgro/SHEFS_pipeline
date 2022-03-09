###########################################
# functions to run SDM analysis
# author Vivienne Groner 
#
# Last update: 09.03.2022
##########################################


# pseudoabsence sampling buffer sensitivity analysis

background_sampler_sens <- 
  function(sp_name, 
           in_dir, 
           out_dir, 
           dens_abs = "absolute", 
           density = NULL, 
           no_pnts = NULL, 
           type = "background", 
           buffer = NULL, 
           polygon = NULL) {
  
  in_file <- list.files(in_dir, full.names = TRUE)[grepl(sp_name, list.files(in_dir))]
  
  sf_int <- read_csv(paste0(in_dir, "/", sp_name, ".csv")) %>%
    dplyr::select("x", "y") %>%
    dplyr::distinct() %>%
    sf::st_as_sf(.,
                 coords = c("x", "y"),
                 crs = 4326) %>%
    sf::st_intersection(., polygon)
  
  bkg_polygon <- polygon %>%
    dplyr::filter(ECO_NAME %in% sf_int$ECO_NAME)
  
  if (dens_abs == "density"){
    no_pnts <- round(as.numeric(sum(st_area(bkg_polygon)))/(1000000*density))   
  }
  
  if (type == "background"){
    points_out <- bkg_polygon %>% 
      sf::st_sample(., size = no_pnts, type = "random")  
  }
  
  if (type == "pseudoabsence"){
    
    diss_bkg_polygon <- sf::st_union(bkg_polygon)
    sf_int_trans <- st_transform(sf_int, "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m no_defs") #robinson projection
    buff_pnts <- sf::st_buffer(sf_int_trans, buffer*1000)  
    buff_pnts <- st_transform(buff_pnts, crs(sf_int)) 
    buff_pnts <- sf::st_union(buff_pnts)
    diff_bkg_polygon <- sf::st_difference(diss_bkg_polygon, buff_pnts)  
    points_out <- diff_bkg_polygon %>% 
      sf::st_sample(., size = no_pnts, type = "random")
    
  }
  
  tibb <- as_tibble(points_out)
  
  sep_df <- tibb %>%
    mutate(x = unlist(purrr::map(tibb$geometry, 1)),
           y = unlist(purrr::map(tibb$geometry, 2))) %>%
    dplyr::select(x, y)
  
  df_out <- data.frame(sp_name, sep_df)
  
  write.csv(df_out,
            file = paste0(out_dir, "/", sp_name, '_den',density,'_buf',buffer,".csv"),
            row.names = FALSE)
  
  print(basename(sp_name))
  return(df_out)
}


# get thresholds

getThresholds <- 
  function(aucs) {
    
  thresholds <- vector(mode = "list", length = 5)
  names(thresholds) <- c("spec_sens",
                         "no_omission",
                         "prevalence",
                         "equal_sens_spec",
                         "tss")
  thresholds[[1]] <- sapply(aucs, function(x)
    dismo::threshold(x, "spec_sens"))
  thresholds[[2]] <- sapply(aucs, function(x)
    dismo::threshold(x, "no_omission"))
  thresholds[[3]] <- sapply(aucs, function(x)
    dismo::threshold(x, "prevalence"))
  thresholds[[4]] <- sapply(aucs, function(x)
    dismo::threshold(x, "equal_sens_spec"))
  thresholds[[5]] <- sapply(aucs, function(x)
    tssCalc(x))
  return(thresholds)
  
}

# calculate TSS

tssCalc <- 
  function(eval) {
    
  res <- data.frame(threshold = eval@t,
                    tss = apply(eval@confusion, 1, function(x) {
                      cm <- t(matrix(rev(x), nrow = 2))
                      dimnames(cm) <- list(pred = c("0", "1"),
                                           obs = c("0", "1"))
                      class(cm) <- "table"
                      sens <- caret::sensitivity(cm)
                      spec <- caret::specificity(cm)
                      tss <- sens + spec - 1
                      return(tss)
                    }))
  thresh <- res$threshold[which.max(res$tss)]
  return(thresh)
  
}

# get eval stats from model 

get_eval <- 
  function(eval_file, 
           threshold = "tss") {
    
  load(eval_file)
  model <- evals$model
  species <- evals$sp_name
  aucs <- mean(sapply(evals$aucs, function(x)
    x@auc))
  
  if (model == "glm") {
    thresholds  <-
      mean(exp(evals$thresholds[[threshold]]) / (1 + exp(evals$thresholds[[threshold]])))
  } else {
    thresholds  <- mean(evals$thresholds[[threshold]])
  }
  
  df_out <-
    data.frame(
      sp_name = species,
      model = model,
      auc = aucs,
      threshold = thresholds
    )
  return(df_out)
  
}


# get eval stats from model all AUCs 

get_eval_auc <- 
  function(eval_file) {
    
    load(eval_file)
    species <- evals$sp_name
    auc1 <- evals$aucs[[1]]@auc
    auc2 <- evals$aucs[[2]]@auc
    auc3 <- evals$aucs[[3]]@auc
    auc4 <- evals$aucs[[4]]@auc
    aucs <- mean(sapply(evals$aucs, function(x)
      x@auc))
    sd <- sd(sapply(evals$aucs, function(x)
      x@auc))
    
    df_out <-
      data.frame(
        sp_name = species,
        aucs = aucs,
        sd = sd,
        auc1 = auc1,
        auc2 = auc2,
        auc3 =auc3,
        auc4 = auc4
      )
    return(df_out)
    
  }


# Ensemble model current

ensemble_model <-
  function(sp_names, 
           eval_df, 
           preds, 
           method = "majority_pa") {
    
    preds_f <- raster::stack(preds[grepl(sp_names, preds)])
    order <- gsub(paste0(sp_names, "_"), "" , names(preds_f))
   
    evals_f <- eval_df[eval_df$sp_name==sp_names,]
    aucs <- evals_f$auc
    
    if (all(order == evals_f$model)) {
      if (method == "majority_pa") {
        ens_preds <- preds_f > evals_f$threshold
        ens_pa <- sum(ens_preds)
        ens_pa[ens_pa < round(raster::nlayers(preds_f))] <- 0
        ens_pa[ens_pa >= round(raster::nlayers(preds_f))] <- 1
        ens_out <- ens_pa
      }
    
    if (method == "weighted") {
        preds_w <- preds_f * aucs
        preds_sum <- sum(preds_w)
        ens_out <- preds_sum / sum(aucs)
        
      }
      
      if (method == "mean") {
        ens_out <- raster::calc(preds_f, mean, na.rm = TRUE)
      }
      
    }
    gc()
    return(ens_out)
    
  }


# Ensemble model future

ensemble_model_f <-
  function(sp_names, 
           eval_df, 
           predsf, 
           method = "majority_pa") {
    
    preds_fu <- raster::stack(predsf[grepl(sp_names, predsf)])
    order <- gsub(paste0(sp_names, "_"), "" , names(preds_fu))
    evals_fu <- eval_df[eval_df$sp_name==sp_names,]
    aucs <- evals_fu$auc

    if (all(order == evals_fu$model)) {
      if (method == "majority_pa") {
        ens_preds <- preds_fu > evals_fu$threshold
        ens_pa <- sum(ens_preds)
        ens_pa[ens_pa < round(raster::nlayers(preds_fu))] <- 0
        ens_pa[ens_pa >= round(raster::nlayers(preds_fu))] <- 1
        ens_out <- ens_pa
      }
      
      if (method == "weighted") {
        preds_w <- preds_fu * aucs
        preds_sum <- sum(preds_w)
        ens_out <- preds_sum / sum(aucs)
        
      }
      
      if (method == "mean") {
        ens_out <- raster::calc(preds_fu, mean, na.rm = TRUE)
        
      }
      
    }
    gc()
    return(ens_out)
  }

# Ensemble model with land use for future

ensemble_model_LUI_f <-
  function(sp_names, 
           eval_df, 
           predsf, 
           lui,
           method = "majority_pa") {
    
    preds_fu <- raster::stack(predsf[grepl(sp_names, predsf)])
    order <- gsub(paste0(sp_names, "_"), "" , names(preds_fu))
    evals_fu <- eval_df[eval_df$sp_name==sp_names,]
    aucs <- evals_fu$auc
    
    if (all(order == evals_fu$model)) {
      if (method == "majority_pa") {
        ens_preds <- preds_fu*lui > evals_fu$threshold
        ens_pa <- sum(ens_preds)
        ens_pa[ens_pa < round(raster::nlayers(preds_fu))] <- 0
        ens_pa[ens_pa >= round(raster::nlayers(preds_fu))] <- 1
        ens_out <- ens_pa
      }
      
      if (method == "weighted") {
        preds_w <- preds_fu*lui * aucs
        preds_sum <- sum(preds_w)
        ens_out <- preds_sum / sum(aucs)
        
      }
      
    }
    gc()
    return(ens_out)
  }

# coalesce by column

coalesce_by_column <- function(df) {
  
  return(dplyr::coalesce(!!! as.list(df)))
  
}


# get species name

get_sp <- function(name){
  
  if(sum(stringr::str_detect(preds, name)) >= 1){
    return(name)
  }
  
}


# get species name future

get_spf <- function(name){
  
  if(sum(stringr::str_detect(predsf, name)) >= 1){
    return(name)
  }
  
}

# Fit bioclim model sensitivity test

fitBC_senstest <-
  function(sp_name,
           pres_dir,
           backg_dir,
           predictor_names,
           predictors,
           pred_out_dir,
           eval_out_dir,
           den,
           buf,
           overwrite,
           threads = 4,
           eval = TRUE) {
    
    print(sp_name)
    
    predictor_names <-
      stringr::str_pad(predictor_names, 2, pad = "0")
    
    CHELSA_predictor_names <-
      paste0("CHELSA_bio10_", predictor_names)
    
    curr <-
      raster::dropLayer(predictors, which(!names(predictors) %in% CHELSA_predictor_names))
    
    pres_pts <-
      as.matrix(data.table::fread(paste0(pres_dir, sp_name, ".csv"), select = c("x", "y")), ncol = 2)
    
    if (nrow(pres_pts) < 10) {
      cat("Fewer than 10 data points - cannot fit model!")
    } else {
      backg_pts <-
        as.matrix(data.table::fread(paste0(backg_dir, sp_name, '_den',den,'_buf',buf, ".csv"), select = c("x", "y")))
      
      cat("Fitting bioclim model...\n")
      bc <- dismo::bioclim(curr, pres_pts)
      saveRDS(bc, file = paste0(eval_out_dir, sp_name,'_',den,'_',buf, "_bioclim_model.RDS"))
      cat("Done.\n")
      cat("...\n")
      cat("Evaluating bioclim model...\n")
      
      kfolds_p <- dismo::kfold(pres_pts, 4)
      kfolds_a <- dismo::kfold(backg_pts, 4)
      
      if (.Platform$OS.type == "unix") {
        cl <- parallel::makeForkCluster(threads)
      } else {
        cl <- parallel::makeCluster(threads)
      }
      parallel::clusterExport(
        cl,
        varlist = c(
          "pres_pts",
          "backg_pts",
          "kfolds_p",
          "kfolds_a",
          "curr",
          "clustEvalPa"
        ),
        envir = environment()
      )
      aucs <- parallel::clusterApply(cl, 1:4, function(x) {
        clustEvalPa(x, pres_pts, backg_pts, kfolds_p, kfolds_a, curr, mod = "bioclim")
      })
      parallel::stopCluster(cl)
      cat("Done.\n")
      cat("...\n")
      thresholds <- getThresholds(aucs)
      
      cat("Predicting from bioclim model...\n")
      res <- dismo::predict(curr, bc)
      
      cat("Done.\n")
      cat("...\n")
      cat("Writing bioclim predictions...\n")
      out_file <- paste0(sp_name,'_',den,'_',buf, "_bioclim.tif")
      
      if (!dir.exists(pred_out_dir)) {
        dir.create(pred_out_dir)
      }
      
      raster::writeRaster(
        res,
        filename = paste(pred_out_dir, out_file, sep = "/"),
        format = "GTiff",
        overwrite = overwrite
      )
      gc()
      
      if (!dir.exists(eval_out_dir)) {
        dir.create(eval_out_dir, recursive = TRUE)
      }
      
      if (eval) {
        evals <- list(sp_name = sp_name, model = "bioclim", aucs = aucs, thresholds = thresholds)
        save(evals, file = paste0(eval_out_dir,  sp_name,'_', den,'_',buf, "_bioclim_eval.RDS"))
      }
    }
    cat("Done.\n")
    cat("...\n")
  }


clustEvalPa <-
  function(num,
           pres_pts,
           backg_pts,
           kfolds_p,
           kfolds_a,
           curr,
           mod) {
    
    pres_train <- pres_pts[kfolds_p != num,]
    pres_test <- pres_pts[kfolds_p == num,]
    backg_test <- backg_pts[kfolds_a == num,]
    if (mod == "bioclim") {
      .m <- dismo::bioclim(curr, pres_train)
    } else if (mod == "domain") {
      .m <- dismo::domain(curr, pres_train)
    }
    e <- dismo::evaluate(pres_test, backg_test, .m, curr)
    return(e)
    
  }


clustEvalSdm <- function(num, sdm_set, kfolds, model, mod) {
  
  train <- sdm_set[kfolds != num,]
  test_p <- sdm_set[kfolds == num & sdm_set[, "pb"] == 1,]
  test_a <- sdm_set[kfolds == num & sdm_set[, "pb"] == 0,]
  if (mod == "glm") {
    .m <- stats::glm(stats::formula(model),
                     data = train,
                     family = binomial(link = "logit"))
  } else if (mod == "rf") {
    .m <-
      suppressWarnings(randomForest::randomForest(model, data = train))
  }
  
  e <- dismo::evaluate(test_p, test_a, .m)
  e
  
}

# END
