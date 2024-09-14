
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
library(future.apply)
library(parallelly)
library(RColorBrewer)
library(lme4)
library(car)

# ver paquetes cargados
# sessionInfo()
# detach alternative devtools::unload("your-package")

# Cargamos Datos

load("./Data/DatosSO2.RData" , verbose = TRUE)
set.seed(123)

# Verificacion de niveles del factor disease
# sapply(all.data$S2_N2,function(x) {levels(x$disease)})


# PLS sin seleccion -------------------------------------------------------

cvPlsDa <-  function(x, 
                     y , 
                     ncomp = 10,
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
      # Ajustar modelo PLS dependiendo del numero de componentes
      for (a in 1:ncomp) {
        # Validacion de variables y numero de componentes
        if ( (ncol(train.data)) <= (a) ) {
          break
        }
        # Ajuste PLS inicial
        pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                  ncomp = a , validation = "none" , 
                                  center = TRUE , scale = TRUE)
        # Predicciones con modelo PLS
        prediction.raw <- predict(pls.fold, test.data, type = "raw")
        prediction.raw <- prediction.raw[,2,]
        # Iteracion a traves de los puntos de corte
        for (cut in cutoff) {
          # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
          pred.list.hyp[[paste(a,cut,sep = "-")]][i] <- ifelse(prediction.raw > cut ,levels(y)[2]  , levels(y)[1])
        }
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


# Funciones de cada método de selección

# VIP Hiperparametros -----------------------------------------------------

cvVip <-  function(x = predictor.pls, 
                   y = disease, 
                   ncomp = 10,
                   vip.nvalues = 30,
                   cutoff = seq(0,1,0.05), 
                   folds = 10, 
                   Rep = 10,
                   seed = 123,
                   varZero = FALSE) {
  # Quitar variables con varianza < 2.22e-16 en caso de que varZero TRUE ( criterio de funcion opls de librería ropls ocupado en TFM anterior)
  # FALSE continua originalmente
  if ( (isTRUE(varZero)) & (sum(apply(x,2,var) < 2.22e-16 ) > 0)) {
    # Como Vip varia en cada base de datos calculamos un vip interval que genere 20 valores
    vip.interval <- round(seq(0,
                              max(plsVarSel::VIP(
                                caret::plsda(
                                  x = x[,-which(apply(x,2,var) < 2.22e-16 )],
                                  y = y,
                                  ncomp = 1,
                                  scale = TRUE
                                ),
                                1
                              )),
                              length.out = vip.nvalues)[2],
                          2)
  } else {
    # Como Vip varia en cada base de datos calculamos un vip interval que genere 20 valores
    vip.interval <- round(seq(0,
                              max(plsVarSel::VIP(
                                caret::plsda(
                                  x = x,
                                  y = y,
                                  ncomp = 1,
                                  scale = TRUE
                                ),
                                1
                              )),
                              length.out = vip.nvalues)[2],
                          2)
  }
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
      if ( (isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 ))  {
        var.0 <- which(diag(cov(train.data)) < 2.22e-16)
        train.data <- train.data[, -var.0]
        test.data <- test.data[, -var.0]
      }
      # Ajustar modelo PLS dependiendo del numero de componentes
      for (a in 1:ncomp) {
        # Calculo del VIP con dataset completo para establacer rango de pruebas VIP
        if ( (isTRUE(varZero)) & (sum(apply(x,2,var) < 2.22e-16 ) > 0) ) {
          Vip.range.all <- plsVarSel::VIP(
            caret::plsda(x = x[,-which(apply(x,2,var) < 2.22e-16 )], 
                         y = y, 
                         ncomp = a, 
                         center = TRUE , 
                         scale = TRUE),
            a
          )
        } else {
          Vip.range.all <- plsVarSel::VIP(
            caret::plsda(x = x, 
                         y = y, 
                         ncomp = a, 
                         center = TRUE , 
                         scale = TRUE),
            a
          )
        }
        
        
        # Validacion de variables y numero de componentes
        if ( (ncol(train.data)) <= (a) ) {
          break
        }
        # Ajuste PLS inicial
        pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                  ncomp = a , validation = "none" , 
                                  center = TRUE , scale = TRUE)
        # Calculo del VIP del fold
        vip.vars <-  plsVarSel::VIP(pls.fold, pls.fold$ncomp)
        # Definimos vipfilter que es el rango de los valores VIP a probar
        vipfilter <-  seq(0,quantile(Vip.range.all,1-((a+1)/ncol(train.data))),vip.interval)
        # Iteracion a traves de los limites VIP
        for (vip in vipfilter) {
          # Variables para el nuevo PLS
          vip.filter <-  vip.vars >= vip
          # En El caso de que las variables seleccioandas sean menor o igual que ncomp, salir de la iteracion
          if ( (sum(vip.filter)) <= (a) ) {
            break
          }
          # Ajuste PLS con seleccion
          pls.prune <- caret::plsda( x = train.data[,vip.filter], y = train.response  , 
                                     ncomp = a , validation = "none",
                                     center = TRUE , scale = TRUE)
          # Predicciones con modelo PLS podado
          prediction.raw <- predict(pls.prune, test.data[,vip.filter] , type = "raw")
          prediction.raw <- prediction.raw[,2,]
          # Iteracion a traves de los puntos de corte
          for (cut in cutoff) {
            # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
            pred.list.hyp[[paste(a,vip,cut,sep = "-")]][i] <- ifelse(prediction.raw > cut ,levels(y)[2]  , levels(y)[1])
          }
        }
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


# Coeficientes de Regresion -----------------------------------------------

# Funcion de pvalores en PLS-DA ------------

p.value.pls <- function(x , y , 
                        ncomp , 
                        R = 1e3,
                        center = TRUE,
                        scale = TRUE) {
  pls.real <- caret::plsda( x = x , y = y , 
                            ncomp = ncomp , 
                            validation = "none",
                            center = center,
                            scale = scale)
  coef.real <- plsVarSel::RC(pls.real , opt.comp = pls.real$ncomp)
  coef.permu <- sapply(1:R, function(i) {
    RC(caret::plsda(x = x , y  = sample(y) ,  
                    ncomp = ncomp , 
                    validation = "none",
                    center = center,
                    scale = scale), ncomp)
  })
  pval <- rowMeans(abs(coef.permu) > abs(coef.real))
  return(pval)
}


# CVCoeficientes -----------

cvCoeficientes <-  function(x = predictor.pls, 
                            y = disease, 
                            ncomp = 10 ,
                            nvar.int = 20,
                            cutoff = seq(0,1,0.05), 
                            folds = 10, 
                            Rep = 10,
                            Rpermu = 1e3,
                            seed = 123,
                            varZero = FALSE){
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
      # Ajustar modelo PLS dependiendo del numero de componentes
      for (a in 1:ncomp) {
        # Validacion de variables y numero de componentes
        if ( (ncol(train.data)) <= (a) ) {
          break
        }
        # Ajusto PLS y calculo coeficientes de regresion y ordeno de lo mas significativo a lo menos
        pvalues <- sort.int(p.value.pls(x = train.data, 
                                        y = train.response, 
                                        ncomp = a, 
                                        R = Rpermu))
        # Definimos nvar que seleccionaremos en cada iteracion
        nvar.seq <- c(seq(nvar.int,ncol(train.data),nvar.int),ncol(train.data))
        # Iteracion a traves de los limites VIP
        for (nvar in nvar.seq) {
          # Variables para el nuevo PLS
          names.var.selected <- names(pvalues[1:nvar])
          # En El caso de que las variables seleccioandas sean menos que ncomp+1, salir de la iteracion
          if ( (length(names.var.selected)) <= (a) ) {
            break
          }
          # Ajuste PLS con seleccion
          pls.prune <- caret::plsda( x = train.data[,names.var.selected], 
                                     y = train.response  , 
                                     ncomp = a , validation = "none",
                                     center = TRUE , scale = TRUE)
          # Predicciones con modelo PLS podado
          prediction.raw <- predict(pls.prune, test.data[,names.var.selected] , type = "raw")
          prediction.raw <- prediction.raw[,2,]
          # Iteracion a traves de los puntos de corte
          for (cut in cutoff) {
            # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
            pred.list.hyp[[paste(a,nvar,cut,sep = "-")]][i] <- ifelse(prediction.raw > cut ,levels(y)[2]  , levels(y)[1])
          }
        }
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

# CV Selectivity Ratio ----------------------------------------------------


cvSelectivityRatio <-  function(x = predictor.pls, 
                                y = disease, 
                                ncomp = 10 ,
                                SR.nvalues = 60,
                                cutoff = seq(0,1,0.05), 
                                folds = 10, 
                                Rep = 10,
                                seed = 123,
                                varZero = FALSE){
  if ((isTRUE(varZero)) & (sum(apply(x,2,var) < 2.22e-16 ) > 0)) {
    SR.interval <- round(seq(0,
                             max(plsVarSel::SR(
                               caret::plsda(
                                 x = x[,-which(apply(x,2,var) < 2.22e-16 )],
                                 y = y,
                                 ncomp = 1,
                                 scale = TRUE
                               ),
                               1,
                               x[,-which(apply(x,2,var) < 2.22e-16 )]
                             )),
                             length.out = SR.nvalues)[2],
                         4)
  } else {
    SR.interval <- round(seq(0,
                             max(plsVarSel::SR(
                               caret::plsda(
                                 x = x,
                                 y = y,
                                 ncomp = 1,
                                 scale = TRUE
                               ),
                               1,
                               x
                             )),
                             length.out = SR.nvalues)[2],
                         4)
  }
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
      # Ajustar modelo PLS dependiendo del numero de componentes
      for (a in 1:ncomp) {
        # Calculo del SR con dataset completo para establacer rango de pruebas SR y Filtro de Varianza < 2.22e-16
        if ((isTRUE(varZero)) & (sum(apply(x,2,var) < 2.22e-16 ) > 0)) {
          SR.range.all <- plsVarSel::SR(caret::plsda(x = x[,-which(apply(x,2,var) < 2.22e-16 )], 
                                                     y = y, 
                                                     ncomp = a, 
                                                     center = TRUE , 
                                                     scale = TRUE),
                                        a,
                                        x[,-which(apply(x,2,var) < 2.22e-16 )])
        } else {
          SR.range.all <- plsVarSel::SR(caret::plsda(x = x, 
                                                     y = y, 
                                                     ncomp = a, 
                                                     center = TRUE , 
                                                     scale = TRUE),
                                        a,
                                        x)
        }
        # Validacion de variables y numero de componentes
        if ( (ncol(train.data)) <= (a) ) {
          break
        }
        # Ajuste PLS inicial
        pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                  ncomp = a , validation = "none" , 
                                  center = TRUE , scale = TRUE)
        # Calculo del SR
        SR.vars <- plsVarSel::SR(pls.fold, pls.fold$ncomp ,train.data)
        # Definimos SRrange que es el rango de los valores SR a probar
        SRrange <-  seq(0,quantile(SR.range.all,1-((a+1)/ncol(train.data))),SR.interval)
        # Iteracion a traves de los limites SR
        for (sr in SRrange) {
          # Filtro de variables menor a cierto SR
          SR.filter <-  SR.vars >= sr
          # En El caso de que las variables seleccioandas sean menos que ncomp+1, salir de la iteracion
          if ( (sum(SR.filter)) <= (a) ) {
            break
          }
          # Ajuste PLS con seleccion
          pls.prune <- caret::plsda( x = train.data[,SR.filter], y = train.response  , 
                                     ncomp = a , validation = "none",
                                     center = TRUE , scale = TRUE)
          # Predicciones con modelo PLS podado
          prediction.raw <- predict(pls.prune, test.data[,SR.filter] , type = "raw")
          prediction.raw <- prediction.raw[,2,]
          # Iteracion a traves de los puntos de corte
          for (cut in cutoff) {
            # prediccion con ese punto de corte - Guardar prediccion de ese cut - sr- comp
            pred.list.hyp[[paste(a,sr,cut,sep = "-")]][i] <- ifelse(prediction.raw > cut ,levels(y)[2]  , levels(y)[1])
          }
        }
      }
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    # Tambien puede contener valores NA y el vector estará imcompleto - sacar esos casos
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

# CV SPLS -----------------------------------------------------------------

# Center - Scale para SPLS

center.scale.f <- function(x, mean.pls , scale.pls) {
  # center
  x.center <- x - rep(mean.pls, each = dim(x)[1])
  # scale
  x.scale <- x.center/rep(scale.pls, each = dim(x)[1])
  return(as.matrix(x.scale))
}



cvSpls <-  function(x = predictor.pls, 
                    y = disease, 
                    ncomp = 10,
                    eta.interval = 0.1,
                    cutoff = seq(0,1,0.05), 
                    folds = 10, 
                    Rep = 10,
                    seed = 123,
                    varZero = FALSE){
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
      if ((isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 )) {
        var.0 <- which(diag(cov(train.data)) < 2.22e-16)
        train.data <- train.data[, -var.0]
        test.data <- test.data[, -var.0]
      }
      # Ajustar modelo PLS dependiendo del numero de componentes
      for (a in 1:ncomp) {
        # Validacion de variables y numero de componentes
        if ( (ncol(train.data)) <= (a) ) {
          break
        }
        # Iteracion a traves de los valore eta
        for (eta in seq(0,1-eta.interval,eta.interval)) {
          # Ajuste PLS con seleccion
          spls.fold <- spls::splsda( x = as.matrix(train.data),
                                     y = train.response, 
                                     K = a,
                                     eta = eta,
                                     scale.x = TRUE
          )
          # Calculo de vector de "prob" en el caso de SPLS (Y = XB_pls)
          prediction.raw <- center.scale.f(test.data,
                                           spls.fold$meanx,
                                           spls.fold$normx) %*% spls.fold$spls.fit$betahat + rep(spls.fold$spls.fit$mu, dim(test.data)[1])
          # Iteracion a traves de los puntos de corte
          for (cut in cutoff) {
            # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
            pred.list.hyp[[paste(a,eta,cut,sep = "-")]][i] <- ifelse(prediction.raw > cut ,levels(y)[2]  , levels(y)[1])
          }
        }
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


# T2 CV -------------------------------------------------------------------

# Funcion de T2 adaptada para plsda y funcion de validacion cruzada adapataciones
# de la funcion de plsVarSel y MSQQC::mult.chart

covariance.msqc <-
  function(x,stat,method,...){
    p <- ncol(x) # quality characteristics
    m <- nrow(x) # sample
    if (inherits(x, "matrix") || inherits(x, "data.frame") ) (
      x <- array(data.matrix(x),c(m,p,1)))
    
    n <- dim(x)[3] # observations or sample size
    
    s.jk <- matrix(0, m, p ^ 2) # matrix of the covariances of observations
    SS <- matrix(0,m,1) # matrix of /S/ statistic 
    
    if(n > 1){
      arrays <- expand.grid(1:p,1:p)
      
      for (i in 1 : m){
        for(j in 1 : p ^ 2){
          s.jk[i,j] <- cov(x[i,arrays[j,1],],x[i,arrays[j,2],])
        }
      } 
      
      S <- matrix(colMeans(s.jk),p,p)
      
      for (ii in 1 : m){
        SS[ii] <- det(matrix(s.jk[ii,],p,p))
      }
      
      if(missing(stat)) (return(S))
      else (return(SS))
      
    }    
    
    if(n == 1){
      if(missing(method))(method="sw")
      
      if(method == "sw"){
        B <- matrix(0,p,p)
        w <- sweep(x,2,(apply(x,2,mean))) #compute de value minus the mean
        for(i in 1:m){
          B <- B + w[i,,] %*% t(w[i,,])
        }
        S <- s1 <- B/(m - 1)
      }
      
      if(method == "hm"){
        V <- matrix(0,m-1,p)
        for(i in 1:m-1){
          V[i,] <- x[i+1,,] - x[i,,]
        }
        S <- s2 <- .5 * t(V) %*% V / (m - 1)
      }
      
      
      return(S)
    }
    
    
  }



T2.mod <- function(x, Xmv, S ,colm , phase = 1 , alpha = 0.05) {
  # Preparacion de los datos
  x <- data.frame(x)
  p <- ncol(x)
  m <- nrow(x)
  if (inherits(x, "matrix") || inherits(x, "data.frame")) 
    (x <- array(data.matrix(x), c(m, p, 1)))
  n <- dim(x)[3]
  if (!missing(Xmv)) 
    (phase <- 2)
  x.jk <- matrix(0, m, p)
  t2 <- matrix(0, m, 1)
  x.jk <- apply(x, 1:2, mean)
  if (missing(Xmv)) 
    (Xmv <- colMeans(x.jk))
  if (missing(S)) 
    (S <- covariance.msqc(x, method = "sw"))
  if (missing(colm)) 
    (colm <- nrow(x))
  # Calculo del T2
  for (ii in 1:m) {
    t2[ii, 1] <- n * t(x.jk[ii, ] - Xmv) %*% solve(S) %*% 
      (x.jk[ii, ] - Xmv)
  }
  # UCL / limite
  ifelse(n == 1, ifelse(phase == 1, ucl <- ((colm - 1)^2)/colm * 
                          qbeta(1 - alpha, p/2, ((colm - p - 1)/2)), ucl <- ((p * 
                                                                                (colm + 1) * (colm - 1))/((colm^2) - colm * p)) * 
                          qf(1 - alpha, p, colm - p)), ifelse(phase == 1, ucl <- (p * 
                                                                                    (colm - 1) * (n - 1))/(colm * n - colm - p + 1) * 
                                                                qf(1 - alpha, p, colm * n - colm - p + 1), ucl <- (p * 
                                                                                                                     (colm + 1) * (n - 1))/(colm * n - colm - p + 1) * 
                                                                qf(1 - alpha, p, colm * n - colm - p + 1)))
  if (any(t2 > ucl)) {
    t3 <- which(t2 > ucl)
    return(t3)
  }
}

# Funcion T2

cvT2pls <-  function(x = predictor.pls, 
                     y = disease, 
                     ncomp = 10,
                     alpha.interval = 0.05,
                     cutoff = seq(0,1,0.05), 
                     folds = 10, 
                     Rep = 10,
                     seed = 123,
                     varZero = FALSE){
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
      if ((isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 )) {
        var.0 <- which(diag(cov(train.data)) < 2.22e-16)
        train.data <- train.data[, -var.0]
        test.data <- test.data[, -var.0]
      }
      # Ajustar modelo PLS dependiendo del numero de componentes
      for (a in 1:ncomp) {
        # Validacion de variables y numero de componentes
        if ( (ncol(train.data)) <= (a) ) {
          break
        }
        # Ajuste PLS inicial
        pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                  ncomp = a , validation = "none" , 
                                  center = TRUE , scale = TRUE)
        # Iteracion a traves de los limites VIP
        for (alpha.l in seq(0.05,1,alpha.interval)) {
          # Calculo y filtrado del T2
          T2.vars <- T2.mod(x = pls.fold$loading.weights[,1:pls.fold$ncomp], 
                            alpha = alpha.l)
          # En El caso de que las variables seleccionadas sean menor o igual que ncomp, salir de la iteracion
          if ( (length(T2.vars)) <= (a) ) {
            break
          }
          # Ajuste PLS con seleccion
          pls.prune <- caret::plsda( x = train.data[,T2.vars], y = train.response  , 
                                     ncomp = a , validation = "none",
                                     center = TRUE , scale = TRUE)
          # Predicciones con modelo PLS podado
          prediction.raw <- predict(pls.prune, test.data[,T2.vars] , type = "raw")
          prediction.raw <- prediction.raw[,2,]
          # Iteracion a traves de los puntos de corte
          for (cut in cutoff) {
            # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
            pred.list.hyp[[paste(a,alpha.l,cut,sep = "-")]][i] <- ifelse(prediction.raw > cut ,levels(y)[2]  , levels(y)[1])
          }
        }
      }
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    # Se agrega revision de Repeticiones , caso de que ninguna repeticion pase 
    # si queda VACIA LA LISTA se asigna "NA"
    if (length(pred.list.hyp.filter) > 0) {
      f1.rep[[k]] <- lapply(pred.list.hyp.filter,
                            function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    } else {
      f1.rep[[k]] <- NA
    }
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep
  ))
}

# Transformaciones

mean.models <- function(x, name.col, sep = "-") {
  x %>% 
    separate(L2,name.col,sep) %>% 
    rename(f1 = value , rep = L1) %>% 
    mutate(across(name.col,as.numeric)) %>% 
    group_by(across(all_of(name.col))) %>% 
    summarise(f1mean = mean(f1) , nrep = n() ,f1sd = sd(f1)) %>% 
    ungroup() %>% 
    arrange(desc(f1mean),desc(f1sd)) %>% 
    kbl(digits = 4) %>% 
    kable_styling(full_width = F) 
}


# Anova y tukey funcion

anova.tukey <- function(rawData, levels.diff) {
  aov.data <- rawData %>% 
    mutate(across(c(L2,L1),as.factor)) %>% 
    filter(L2 %in% levels.diff)
  
  aov.model <- aov(value ~ L2 + Error(L1/L2) , data = aov.data)
  op = options()
  options(contrasts = c("contr.sum", "contr.poly"))
  tukey.model <- emmeans(aov.model, ~ L2)
  pairs(tukey.model)  # Tukey's
  tukey.plot <- plot(tukey.model, 
                     comparisons = TRUE) + 
    ggtitle("Tukey") + 
    theme(plot.title = element_text(hjust = 0.5))
  return(list("Anova" = aov.model,"Tukey" = tukey.plot))
}


# BVE CV  -----------------------------------------------------------------

cvBVE <-  function(x = predictor.pls, 
                   y = disease, 
                   ncomp = 10,
                   vip.nvalues = 30,
                   cutoff = seq(0,1,0.05), 
                   folds = 10, 
                   Rep = 10,
                   seed = 123,
                   varZero = FALSE){
  # library for pipe
  library(dplyr)
  # Como Vip varia en cada base de datos calculamos un vip interval que genere 20 valores
  # Calculamos vip maximo de la base de datos
  if ((isTRUE(varZero)) & (sum(apply(x,2,var) < 2.22e-16 ) > 0)) { # Filtro variables con baja variabilidad
    vip.range <- plsVarSel::VIP(caret::plsda(x = x[,-which(apply(x,2,var) < 2.22e-16 )],
                                             y = y,
                                             ncomp = ncomp,
                                             center = TRUE,
                                             scale = TRUE),
                                ncomp)
  } else {
    vip.range <- plsVarSel::VIP(caret::plsda(x = x,
                                             y = y,
                                             ncomp = ncomp,
                                             center = TRUE,
                                             scale = TRUE),
                                ncomp)
  }
  # Calculamos intervalo del vip para generar cierta cantidad de cortes
  vip.interval <- round(seq(0,max(vip.range),length.out = vip.nvalues)[2],2)
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
    pred.fold.vip <- list()
    # Iteracion por cada fold
    for (i in folds.k) {
      for (vip in seq(0,max(vip.range),vip.interval)) {
        for (a in 1:ncomp) {
          train.data <- x[-i,]
          test.data <- x[i,]
          train.response <- y[-i]
          # Filtro variables con baja variabilidad
          if ((isTRUE(varZero)) & (sum(diag(cov(train.data)) < 2.22e-16) > 0 )) {
            var.0 <- which(diag(cov(train.data)) < 2.22e-16)
            train.data <- train.data[, -var.0]
            test.data <- test.data[, -var.0]
          }
          iter <- 1
          terminated <- FALSE
          # Validacion de variables y numero de componentes
          if ( (ncol(train.data)) <= (a) ) {
            break
          }
          while (!terminated){
            # Ajusto PLS inicial
            pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                      ncomp = a , validation = "none" , 
                                      center = TRUE , scale = TRUE)
            # Calculo predicciones y F1 sobre test para eso necesito punto de corte
            for (cut in cutoff) {
              prediction.raw <- predict(pls.fold, test.data, type = "raw")
              prediction.raw <- prediction.raw[,2,]
              pred.fold.vip[[paste(a,vip,cut,iter,sep="-")]][i] <- ifelse(prediction.raw > cut ,levels(y)[2]  , levels(y)[1])
            }
            # Calculo VIP
            VIP.value <- plsVarSel::VIP(pls.fold, pls.fold$ncomp)
            # Quito conjunto de variables poco informativas , es decir dejo las que estan por encima del threshold VIP
            VIP.vars <- VIP.value > vip
            # Si numero de variables es menor al ncomp también termino aqui el algoritmo o tmb en caso
            # de que ya no disminuya mas las variables poco informativas respecto del VIP , ej menores q 1 y limite es 1.17
            if ((sum(VIP.vars) <= a) | (sum(VIP.vars) == ncol(train.data))) {
              terminated <- TRUE
            } else {
              # Actualizo conjuntos de datos con las variables que quedaron
              train.data <- train.data[,VIP.vars]
              test.data <- test.data[,VIP.vars]
            }
            iter <- iter + 1
          }
        }
      }
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.vip.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.fold.vip)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.vip.filter
    # Guardarmos el F1 Score
    f1.list <- lapply(pred.list.vip.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.vip.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
    # Se elige mejor F1 score (mejor iteracion) para cada grupo de ncomp , vip y punto de corte y se guarda en esa repeticion
    f1.rep[[k]] <- f1.list %>%
      reshape2::melt() %>%
      tidyr::separate(L1, c("ncomp", "vip", "cut", "iter"), "-") %>%
      dplyr::group_by(ncomp, vip, cut) %>%
      dplyr::slice(which.max(value)) %>%
      dplyr::ungroup() %>% 
      dplyr::select(value, ncomp, vip, cut)
  }
  return(list("f1rep" = f1.rep,
              "pred"  = pred.rep,
              "other" = meas.rep
  ))
}
