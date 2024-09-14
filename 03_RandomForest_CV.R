
library(caret)
library(ropls)
library(SummarizedExperiment)
library(pROC)
library(pls)
library(spls)
options(rgl.useNULL = TRUE)
library(rgl)
library(plsVarSel)
library(MLmetrics)
library(tidyverse)
library(reshape2)
library(emmeans)
library(microbenchmark)
library(knitr)
library(kableExtra)
library(parallel)
library(doParallel)
library(doRNG)
library(furrr)
library(purrr)
library(future.apply)
library(parallelly)
library(RColorBrewer)
library(lme4)
library(car)
library(randomForest)
library(Boruta)
library(varSelRF)
library(VSURF)

# Random Forest funcion sin seleccion  -----------------------------------------------------------

randomForestCv <-  function(x, 
                            y, 
                            cutoff = seq(0,1,0.05), 
                            folds = 10, 
                            Rep = 10,
                            seed = 123,
                            varZero = FALSE){
  # lista para guardar f1
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Iteracion por cada fold
    for (i in folds.k) {
      train.data <- x[-i,]
      test.data <- x[i,]
      train.response <- y[-i]
      # Quitar variables con varianza < 2.22e-16 ( criterio de funcion opls de librería ropls ocupado en TFM anterior)
      if ((isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 )) {
        var.0 <- which(diag(cov(train.data)) < 2.22e-16)
        train.data <- train.data[, -var.0]
        test.data <- test.data[, -var.0]
      }
      # Ajuste Random Forest 
      rf.fold <- tryCatch(randomForest::randomForest(x = train.data, 
                                                     y = train.response,
                                                     ntree = 500,
                                                     mtry = sqrt(ncol(train.data))),
                          error = function(cond)"skip")
      # En caso de error
      if (class(rf.fold)[1] == "character") {
        pred.rf <- NA
      } else {
        # Predicciones
        pred.rf <- predict(rf.fold, test.data , type = "prob")
        pred.rf <- pred.rf[,"d"]
      }
      # Iteracion a traves de los puntos de corte
      for (cut in cutoff) {
        # prediccion con ese punto de corte - Guardar prediccion de ese cut
        pred.list.hyp[[as.character(cut)]][i] <- ifelse(pred.rf > cut ,levels(y)[2]  , levels(y)[1])
      }
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    f1.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep
  ))
}




# Boruta Random Forest  -----------------------------------------------------------

BorutaRandomForestCv <-  function(x, 
                                  y, 
                                  cutoff = seq(0,1,0.05), 
                                  folds = 10, 
                                  Rep = 10,
                                  seed = 123,
                                  varZero = FALSE){
  # lista para guardar f1
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Iteracion por cada fold
    for (i in folds.k) {
      train.data <- x[-i,]
      test.data <- x[i,]
      train.response <- y[-i]
      # Quitar variables con varianza < 2.22e-16 ( criterio de funcion opls de librería ropls ocupado en TFM anterior)
      if ((isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 )) {
        var.0 <- which(diag(cov(train.data)) < 2.22e-16)
        train.data <- train.data[, -var.0]
        test.data <- test.data[, -var.0]
      }
      # Aplicamos Boruta con pvalue 0.05
      boruta.run <- tryCatch(Boruta::Boruta(x = train.data,
                                            y = train.response),
                             error=function(cond)"skip")
      # Si falla seleccion
      if (class(boruta.run)[1] == "character") {
        pred.rf <- NA
      } else {
        # Variables seleccionadas
        var.selected <- Boruta::getSelectedAttributes(boruta.run, withTentative = TRUE)
        # Modelo random Forest
        rf.fold <- tryCatch(randomForest::randomForest(x = train.data[,var.selected],
                                                       y = train.response,
                                                       ntree = 500,
                                                       mtry = sqrt(length(var.selected))),
                            error=function(cond)"skip")
        if(class(rf.fold)[1] == "character") {
          pred.rf <- NA
        } else {
          # Predicciones
          pred.rf <- predict(rf.fold, test.data[,var.selected], type = "prob")
          pred.rf <- pred.rf[,"d"]
        }
      }
      # Iteracion a traves de los puntos de corte
      for (cut in cutoff) {
        # prediccion con ese punto de corte - Guardar prediccion de ese cut
        pred.list.hyp[[as.character(cut)]][i] <- ifelse(pred.rf > cut ,levels(y)[2]  , levels(y)[1])
      }
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    f1.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep
  ))
}


# varSelRF Random Forest  -----------------------------------------------------------

VarSelRfCv <-  function(x, 
                        y, 
                        cutoff = seq(0,1,0.05), 
                        folds = 10, 
                        Rep = 10,
                        seed = 123,
                        varZero = FALSE){
  # lista para guardar f1
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Iteracion por cada fold
    for (i in folds.k) {
      train.data <- x[-i,]
      test.data <- x[i,]
      train.response <- y[-i]
      # Quitar variables con varianza < 2.22e-16 ( criterio de funcion opls de librería ropls ocupado en TFM anterior)
      if ((isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 )) {
        var.0 <- which(diag(cov(train.data)) < 2.22e-16)
        train.data <- train.data[, -var.0]
        test.data <- test.data[, -var.0]
      }
      # Aplicamos varSelRF
      varSelRF.run <- tryCatch(varSelRF::varSelRF(train.data,
                                                  train.response,
                                                  ntree = 500,
                                                  ntreeIterat = 500,
                                                  whole.range = FALSE,
                                                  keep.forest = TRUE),
                               error=function(cond)"skip")
      # Si funcion falla
      if (class(varSelRF.run)[1] == "character") {
        pred.rf <- NA
      } else {
        # Predicciones
        pred.rf <- predict(varSelRF.run$rf.model, test.data[,varSelRF.run$selected.vars], type = "prob")
        pred.rf <- pred.rf[,"d"]
      }
      # Iteracion a traves de los puntos de corte
      for (cut in cutoff) {
        # prediccion con ese punto de corte - Guardar prediccion de ese cut
        pred.list.hyp[[as.character(cut)]][i] <- ifelse(pred.rf > cut ,levels(y)[2]  , levels(y)[1])
      }
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    f1.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep
  ))
}

# Sys.time()
# test.varSelRF <- VarSelRfCv(predictor.pls,disease)
# Sys.time() # 13 min 




# VSURF Random Forest  -----------------------------------------------------------

VsurfCv <-  function(x, 
                     y, 
                     cutoff = seq(0,1,0.05), 
                     folds = 10, 
                     Rep = 10,
                     seed = 123,
                     varZero = FALSE){
  # lista para guardar f1
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Iteracion por cada fold
    for (i in folds.k) {
      train.data <- x[-i,]
      test.data <- x[i,]
      train.response <- y[-i]
      # Quitar variables con varianza < 2.22e-16 ( criterio de funcion opls de librería ropls ocupado en TFM anterior)
      if ((isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 )) {
        var.0 <- which(diag(cov(train.data)) < 2.22e-16)
        train.data <- train.data[, -var.0]
        test.data <- test.data[, -var.0]
      }
      # Aplicamos Vsurf
      Vsurf.run <- tryCatch(VSURF::VSURF(x = train.data,
                                         y = train.response),
                            error = function(cond)"skip")
      # si funcion de seleccion falla
      if (class(Vsurf.run)[1] == "character") {
        pred.rf <- NA
      } else {
        # Variables seleccionadas con vsurf
        var.selected <- Vsurf.run$varselect.pred
        # Ajuste Random Forest con variables seleccionadas
        rf.fold <- tryCatch(randomForest::randomForest(train.response ~. , 
                                                       data = train.data[,var.selected],
                                                       ntree = 500,
                                                       mtry = sqrt(length(var.selected))),
                            error=function(cond)"skip")
        # Verificar si run fue exitoso
        if (class(rf.fold)[1] == "character") {
          pred.rf <- NA
        } else {
          # Predicciones
          pred.rf <- predict(rf.fold, test.data[,var.selected], type = "prob")
          pred.rf <- pred.rf[,"d"]
        }
      }
      # Iteracion a traves de los puntos de corte
      for (cut in cutoff) {
        # prediccion con ese punto de corte - Guardar prediccion de ese cut
        pred.list.hyp[[as.character(cut)]][i] <- ifelse(pred.rf > cut ,levels(y)[2]  , levels(y)[1])
      }
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    f1.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep
  ))
}


# 
# Sys.time() # 19:20
# 
# test.vsurf.rep <- VsurfCv(x = predictor.pls , y = disease)
# 
# Sys.time() # 20:23 
# 
# # Tiempo VSURF de 1 hora y 03 minutos




