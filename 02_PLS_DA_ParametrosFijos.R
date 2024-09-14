

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

## PLS sin seleccion fixed -------------------------------------------------------

cvPlsDa.fixed <-  function(x, 
                           y , 
                           ncomp = 10,
                           cutoff = 0.5, 
                           folds = 10, 
                           Rep = 10,
                           seed = 1234,
                           varZero = FALSE){
  # lista para guardar f1
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  # lista para guardar variables seleccionadas
  var.rep <- vector("list",Rep)
  # semilla
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Lista para guardar variables seleccionadas en el fold
    var.list <- list()
    j = 1
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
      # Validacion de variables y numero de componentes
      if ( (ncol(train.data)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS inicial
      pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                ncomp = ncomp , validation = "none" , 
                                center = TRUE , scale = TRUE)
      # Predicciones con modelo PLS
      prediction.raw <- predict(pls.fold, test.data, type = "raw")
      prediction.raw <- prediction.raw[,2,]
      # Iteracion a traves de los puntos de corte
      # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
      pred.list.hyp[[paste(ncomp,cutoff,sep = "-")]][i] <- ifelse(prediction.raw > cutoff ,levels(y)[2]  , levels(y)[1])
      # Guardar variables seleccionadas en ese fold
      var.list[[j]] <- names(train.data)
      j <- j+1
      
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    f1.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
    # Guardar variables seleccionadas en esa rep
    var.rep[[k]] <- var.list
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep,
              "var" =  var.rep
  ))
}


## VIP Hiperparametros -----------------------------------------------------

cvVip.fixed <-  function(x, 
                         y, 
                         ncomp = 10,
                         vip.value = 1,
                         cutoff = 0.5, 
                         folds = 10, 
                         Rep = 10,
                         seed = 1234,
                         varZero = FALSE) {
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  # lista para guardar variables seleccionadas
  var.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Lista para guardar variables seleccionadas en el fold
    var.list <- list()
    j = 1
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
      # Validacion de variables y numero de componentes
      if ( (ncol(train.data)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS inicial
      pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                ncomp = ncomp , validation = "none" , 
                                center = TRUE , scale = TRUE)
      # Calculo del VIP del fold
      vip.vars <-  plsVarSel::VIP(pls.fold, pls.fold$ncomp)
      # Variables para el nuevo PLS
      vip.filter <-  vip.vars >= vip.value
      # En El caso de que las variables seleccioandas sean menor o igual que ncomp, salir de la iteracion
      if ( (sum(vip.filter)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS con seleccion
      pls.prune <- caret::plsda( x = train.data[,vip.filter], y = train.response  , 
                                 ncomp = ncomp , validation = "none",
                                 center = TRUE , scale = TRUE)
      # Predicciones con modelo PLS podado
      prediction.raw <- predict(pls.prune, test.data[,vip.filter] , type = "raw")
      prediction.raw <- prediction.raw[,2,]
      # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
      pred.list.hyp[[paste(ncomp,vip.value,cutoff,sep = "-")]][i] <- ifelse(prediction.raw > cutoff ,levels(y)[2]  , levels(y)[1])
      # Guardar variables seleccionadas en ese fold
      var.list[[j]] <- names(which(vip.filter))
      j <- j+1
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    if (length(pred.list.hyp.filter) > 0) {
      f1.rep[[k]] <- lapply(pred.list.hyp.filter,
                            function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    } else {
      f1.rep[[k]] <- NA
    }
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
    # Guardar variables seleccionadas en esa rep
    var.rep[[k]] <- var.list
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep,
              "var"   = var.rep
  ))
}

## Coeficientes de Regresion -----------------------------------------------


### Funcion de pvalores en PLS-DA ------------

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


### CVCoeficientes fixed-----------

cvCoeficientes.fixed <-  function(x , 
                                  y, 
                                  ncomp = 10 ,
                                  nvar.int = 100,
                                  cutoff = 0.5, 
                                  folds = 10, 
                                  Rep = 10,
                                  Rpermu = 1e3,
                                  seed = 1234,
                                  varZero = FALSE){
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  # lista para guardar variables seleccionadas
  var.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Lista para guardar variables seleccionadas en el fold
    var.list <- list()
    j = 1
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
      # Validacion de variables y numero de componentes
      if ( (ncol(train.data)) <= (ncomp) ) {
        break
      }
      # Ajusto PLS y calculo coeficientes de regresion y ordeno de lo mas significativo a lo menos
      pvalues <- sort.int(p.value.pls(x = train.data, 
                                      y = train.response, 
                                      ncomp = ncomp, 
                                      R = Rpermu))
      # Variables para el nuevo PLS
      names.var.selected <- names(pvalues[1:nvar.int])
      # En El caso de que las variables seleccioandas sean menos que ncomp+1, salir de la iteracion
      if ( (length(names.var.selected)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS con seleccion
      pls.prune <- caret::plsda( x = train.data[,names.var.selected], 
                                 y = train.response  , 
                                 ncomp = ncomp , validation = "none",
                                 center = TRUE , scale = TRUE)
      # Predicciones con modelo PLS podado
      prediction.raw <- predict(pls.prune, test.data[,names.var.selected] , type = "raw")
      prediction.raw <- prediction.raw[,2,]
      # Iteracion a traves de los puntos de corte
      # prediccion con ese punto de corte - Guardar prediccion de ese combinacion
      pred.list.hyp[[paste(ncomp,nvar.int,cutoff,sep = "-")]][i] <- ifelse(prediction.raw > cutoff ,levels(y)[2]  , levels(y)[1])
      
      # Guardar variables seleccionadas en ese fold
      var.list[[j]] <- names.var.selected
      j <- j+1
      
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    if (length(pred.list.hyp.filter) > 0) {
      f1.rep[[k]] <- lapply(pred.list.hyp.filter,
                            function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    } else {
      f1.rep[[k]] <- NA
    }
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
    
    # Guardar variables seleccionadas en esa rep
    var.rep[[k]] <- var.list
    
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep,
              "var" = var.rep
  ))
}


## CV Selectivity Ratio fixed ----------------------------------------------------


cvSelectivityRatio.fixed <-  function(x, 
                                      y, 
                                      ncomp = 10 ,
                                      SR.val = 0.01,
                                      cutoff = 0.05, 
                                      folds = 10, 
                                      Rep = 10,
                                      seed = 1234,
                                      varZero = FALSE){
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  # lista para guardar variables seleccionadas
  var.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Lista para guardar variables seleccionadas en el fold
    var.list <- list()
    j = 1
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
      # Validacion de variables y numero de componentes
      if ( (ncol(train.data)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS inicial
      pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                ncomp = ncomp , validation = "none" , 
                                center = TRUE , scale = TRUE)
      # Calculo del SR
      SR.vars <- plsVarSel::SR(pls.fold, pls.fold$ncomp ,train.data)
      # Filtro de variables menor a cierto SR
      SR.filter <-  SR.vars >= SR.val
      # En El caso de que las variables seleccioandas sean menos que ncomp+1, salir de la iteracion
      if ( (sum(SR.filter)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS con seleccion
      pls.prune <- caret::plsda( x = train.data[,SR.filter], y = train.response  , 
                                 ncomp = ncomp , validation = "none",
                                 center = TRUE , scale = TRUE)
      # Predicciones con modelo PLS podado
      prediction.raw <- predict(pls.prune, test.data[,SR.filter] , type = "raw")
      prediction.raw <- prediction.raw[,2,]
      # prediccion con ese punto de corte - Guardar prediccion de ese cut - sr- comp
      pred.list.hyp[[paste(ncomp,SR.val,cutoff,sep = "-")]][i] <- ifelse(prediction.raw > cutoff ,levels(y)[2]  , levels(y)[1])
      # Guardar variables seleccionadas en ese fold
      var.list[[j]] <- names(which(SR.filter))
      j <- j+1
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    # Tambien puede contener valores NA y el vector estará imcompleto - sacar esos casos
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    if (length(pred.list.hyp.filter) > 0) {
      f1.rep[[k]] <- lapply(pred.list.hyp.filter,
                            function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    } else {
      f1.rep[[k]] <- NA
    }
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
    # Guardar variables seleccionadas en esa rep
    var.rep[[k]] <- var.list
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep,
              "var"   = var.rep
  ))
}

## CV SPLS fixed -----------------------------------------------------------------

# Center - Scale para SPLS

center.scale.f <- function(x, mean.pls , scale.pls) {
  # center
  x.center <- x - rep(mean.pls, each = dim(x)[1])
  # scale
  x.scale <- x.center/rep(scale.pls, each = dim(x)[1])
  return(as.matrix(x.scale))
}



cvSpls.fixed <-  function(x, 
                          y, 
                          ncomp = 10,
                          eta.val = 0.1,
                          cutoff = 0.5, 
                          folds = 10, 
                          Rep = 10,
                          seed = 1234,
                          varZero = FALSE){
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  # lista para guardar variables seleccionadas
  var.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Lista para guardar variables seleccionadas en el fold
    var.list <- list()
    j = 1
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
      # Validacion de variables y numero de componentes
      if ( (ncol(train.data)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS con seleccion
      spls.fold <- spls::splsda( x = as.matrix(train.data),
                                 y = train.response, 
                                 K = ncomp,
                                 eta = eta.val,
                                 scale.x = TRUE
      )
      # Calculo de vector de "prob" en el caso de SPLS (Y = XB_pls)
      prediction.raw <- center.scale.f(test.data,
                                       spls.fold$meanx,
                                       spls.fold$normx) %*% spls.fold$spls.fit$betahat + rep(spls.fold$spls.fit$mu, dim(test.data)[1])
      # Prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
      pred.list.hyp[[paste(ncomp,eta.val,cutoff,sep = "-")]][i] <- ifelse(prediction.raw > cutoff ,levels(y)[2]  , levels(y)[1])
      
      # Guardar variables seleccionadas en ese fold
      var.list[[j]] <- names(train.data[,spls.fold$A])
      j <- j+1
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.hyp.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.list.hyp)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.hyp.filter
    # Guardarmos el F1 Score
    if (length(pred.list.hyp.filter) > 0) {
      f1.rep[[k]] <- lapply(pred.list.hyp.filter,
                            function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    } else {
      f1.rep[[k]] <- NA
    }
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.hyp.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
    
    # Guardar variables seleccionadas en esa rep
    var.rep[[k]] <- var.list
    
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep,
              "var"   = var.rep
  ))
}


## T2 CV fixed -------------------------------------------------------------------

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

cvT2pls.fixed <-  function(x, 
                           y, 
                           ncomp = 10,
                           alpha.val = 0.05,
                           cutoff = 0.5, 
                           folds = 10, 
                           Rep = 10,
                           seed = 1234,
                           varZero = FALSE){
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  # lista para guardar variables seleccionadas
  var.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.list.hyp <- list()
    # Lista para guardar variables seleccionadas en el fold
    var.list <- list()
    j = 1
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
      # Validacion de variables y numero de componentes
      if ( (ncol(train.data)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS inicial
      pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                ncomp = ncomp , validation = "none" , 
                                center = TRUE , scale = TRUE)
      # Calculo y filtrado del T2
      T2.vars <- T2.mod(x = pls.fold$loading.weights[,1:pls.fold$ncomp], 
                        alpha = alpha.val)
      # En El caso de que las variables seleccionadas sean menor o igual que ncomp, salir de la iteracion
      if ( (length(T2.vars)) <= (ncomp) ) {
        break
      }
      # Ajuste PLS con seleccion
      pls.prune <- caret::plsda( x = train.data[,T2.vars], y = train.response  , 
                                 ncomp = ncomp , validation = "none",
                                 center = TRUE , scale = TRUE)
      # Predicciones con modelo PLS podado
      prediction.raw <- predict(pls.prune, test.data[,T2.vars] , type = "raw")
      prediction.raw <- prediction.raw[,2,]
      # prediccion con ese punto de corte - Guardar prediccion de ese cut - vip - comp
      pred.list.hyp[[paste(ncomp,alpha.val,cutoff,sep = "-")]][i] <- ifelse(prediction.raw > cutoff ,levels(y)[2]  , levels(y)[1])
      
      # Guardar variables seleccionadas en ese fold
      var.list[[j]] <- names(train.data[,T2.vars])
      j <- j+1
      
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
    
    # Guardar variables seleccionadas en esa rep
    var.rep[[k]] <- var.list
    
  }
  return(list("f1rep" = reshape2::melt(f1.rep),
              "pred"  = pred.rep,
              "other" = meas.rep,
              "var"   = var.rep
  ))
}

## BVE CV fixed  -----------------------------------------------------------------

cvBVE.fixed <-  function(x, 
                         y, 
                         ncomp = 10,
                         vip.val = 1,
                         cutoff = 0.5, 
                         folds = 10, 
                         Rep = 10,
                         seed = 1234,
                         varZero = FALSE){
  # library for pipe
  library(dplyr)
  f1.rep <- vector("list",Rep)
  # lista para guardar otras métricas
  meas.rep <- vector("list",Rep)
  # lista para guardar predicciones
  pred.rep <- vector("list",Rep)
  # lista para guardar variables seleccionadas
  var.rep <- vector("list",Rep)
  set.seed(seed)
  # Repeticiones
  for (k in 1:Rep) {
    folds.k <-  caret::createFolds(y = y, k = folds )
    # Lista para guardar las predicciones con la combinacion de los parámetros
    pred.fold.vip <- list()
    # Lista para guardar variables seleccionadas en el fold
    var.list <- list()
    j = 1
    # Iteracion por cada fold
    for (i in folds.k) {
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
      if ( (ncol(train.data)) <= (ncomp) ) {
        break
      }
      var.iter <- list()
      while (!terminated){
        # Ajusto PLS inicial
        pls.fold <- caret::plsda( x = train.data , y = train.response , 
                                  ncomp = ncomp , validation = "none" , 
                                  center = TRUE , scale = TRUE)
        # Calculo predicciones y F1 sobre test para eso necesito punto de corte
        prediction.raw <- predict(pls.fold, test.data, type = "raw")
        prediction.raw <- prediction.raw[,2,]
        pred.fold.vip[[paste(ncomp,vip.val,cutoff,iter,sep="-")]][i] <- ifelse(prediction.raw > cutoff ,levels(y)[2]  , levels(y)[1])
        # Guardo variables seleccionadas de iteracion actual
        var.iter[[iter]] <- names(train.data)
        # Calculo VIP
        VIP.value <- plsVarSel::VIP(pls.fold, pls.fold$ncomp)
        # Quito conjunto de variables poco informativas , es decir dejo las que estan por encima del threshold VIP
        VIP.vars <- VIP.value > vip.val
        # Si numero de variables es menor al ncomp también termino aqui el algoritmo o tmb en caso
        # de que ya no disminuya mas las variables poco informativas respecto del VIP , ej menores q 1 y limite es 1.17
        if ((sum(VIP.vars) <= ncomp) | (sum(VIP.vars) == ncol(train.data))) {
          terminated <- TRUE
        } else {
          # Actualizo conjuntos de datos con las variables que quedaron
          train.data <- train.data[,VIP.vars]
          test.data <- test.data[,VIP.vars]
        }
        iter <- iter + 1
        # Tener en cuenta que desde iteración 2 recién se esta aplicando filtro VIP!
      }
      # Guardar variables seleccionadas en el fold
      var.list[[j]] <- var.iter
      # contador para folds
      j = j+1
    }
    # Sacamos de lista combinaciones de parámetros que no se ajustaron para todos los folds , es decir vector de predicciones tendran un largo menor
    pred.list.vip.filter <- Filter(function(x) (length(x) == length(y)) & !any(is.na(x)) , pred.fold.vip)
    # Guardamos vector de predicciones
    pred.rep[[k]] <- pred.list.vip.filter
    # Guardarmos el F1 Score
    if (length(pred.list.vip.filter) > 0) {
      f1.list <- lapply(pred.list.vip.filter,
                        function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass["F1"])
    } else {
      f1.list <- NA
    }
    # Guardamos otras métricas
    meas.rep[[k]] <- lapply(pred.list.vip.filter,function(r) caret::confusionMatrix(table(factor(r,levels(y)),y), positive = "d")$byClass)
    
    # Seleccionar mejor iteracion
    iter.opt <- f1.list %>%
      reshape2::melt() %>%
      tidyr::separate(L1, c("ncomp", "vip", "cut", "iter"), "-") %>%
      dplyr::slice(which.max(value)) %>%
      dplyr::select(iter) %>% 
      as.integer()
    # Guardar variables seleccionadas en esa rep  - se extraen variables de cada fold de la mejor iteracion
    var.rep[[k]] <- lapply(var.list , function(x) x[[iter.opt]])
    # Se elige mejor F1 score (mejor iteracion) para cada grupo de ncomp , vip y punto de corte y se guarda en esa repeticion
    f1.rep[[k]] <- f1.list %>%
      reshape2::melt() %>%
      tidyr::separate(L1, c("ncomp", "vip", "cut", "iter"), "-") %>%
      dplyr::group_by(ncomp, vip, cut) %>%
      dplyr::slice(which.max(value)) %>%
      dplyr::ungroup() %>% 
      dplyr::select(value, ncomp, vip, cut)
  }
  return(list("f1rep" = f1.rep%>% dplyr::bind_rows(.id = "L1") %>% 
                mutate(L1 = as.integer(L1)) %>% 
                tidyr::unite(L2,c("ncomp","vip","cut"),sep="-") %>% 
                select(value,L2,L1),
              "pred"  = pred.rep,
              "other" = meas.rep,
              "var" = var.rep
  ))
}







