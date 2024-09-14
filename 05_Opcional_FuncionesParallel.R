

# Paralelización PLS-DA para CV -------------------------------------------

# parallel all methods ----------------------------------------------------

parallel.all.methods <- function(dataset, 
                                 seed = 123,
                                 varZero = FALSE) {
  start.time <- Sys.time()
  # n.cores <- detectCores() - 2
  # usar availableCores() para evitar uso de muchos recursos
  n.cores <- parallelly::availableCores()
  my.cluster <- parallel::makeCluster(n.cores, 
                                      type = "PSOCK")
  # Importamos librerias o funciones para ejecutar paralelamente
  clusterEvalQ(my.cluster,library(plsVarSel))
  clusterExport(my.cluster, 'center.scale.f')
  clusterExport(my.cluster, 'p.value.pls')
  clusterExport(my.cluster, 'T2.mod')
  clusterExport(my.cluster, 'cvVip')
  clusterExport(my.cluster, 'cvCoeficientes')
  clusterExport(my.cluster, 'cvSelectivityRatio')
  clusterExport(my.cluster, 'cvSpls')
  clusterExport(my.cluster, 'cvT2pls')
  clusterExport(my.cluster, 'cvBVE')
  doParallel::registerDoParallel(cl = my.cluster)
  registerDoRNG()
  set.seed(seed)
  # parallel VIP
  Vip.res <- foreach(x = dataset) %dopar% {
    cvVip(x = subset(x,select = -c(disease)),
          y = x[,"disease"], varZero = varZero)
  }
  # parallel RC
  RC.res <- foreach(x = dataset) %dopar% {
    cvCoeficientes(x = subset(x,select = -c(disease)),
                   y = x[,"disease"],varZero = varZero)
  }
  # parallel SR
  SR.res <- foreach(x = dataset) %dopar% {
    cvSelectivityRatio(x = subset(x,select = -c(disease)),
                       y = x[,"disease"],
                       varZero = varZero)
  }
  # parallel SPLS
  SPLS.res <- foreach(x = dataset) %dopar% {
    cvSpls(x = subset(x,select = -c(disease)),
           y = x[,"disease"],
           varZero = varZero)
  }
  # parallel T2
  T2.res <- foreach(x = dataset) %dopar% {
    cvT2pls(x = subset(x,select = -c(disease)),
            y = x[,"disease"],
            varZero = varZero)
  }
  # parallel BVE
  BVE.res <- foreach(x = dataset) %dopar% {
    cvBVE(x = subset(x,select = -c(disease)),
          y = x[,"disease"],
          varZero = varZero)
  }
  stopCluster(my.cluster)
  end.time <- Sys.time()
  return(list("Resultados" = list("Vip" = Vip.res, 
                                  "RC" =  RC.res,
                                  "SR" = SR.res,
                                  "SPLS" = SPLS.res ,
                                  "T2" = T2.res,
                                  "BVE" = BVE.res
  ),
  "Tiempo" = end.time - start.time))
}


## S1 - N1 -----------------------------------------------------------------

# varZeroTrue incorpora filtro de variables con baja variabilidad.. que originan problemas en el ajuste del PLS
# en la normalizacion de tipo S1 (misma condicion que pone funcion ropls)

# final parallel

# all.S1.N1.V2 <- parallel.all.methods(S1.N1.data, varZero = TRUE)
# save(all.S1.N1.V2, file = paste0("RData Seleccion/","all.S1.N1.V2",".RData"))


## S1 - N2 -----------------------------------------------------------------

# all.S1.N2.V2 <- parallel.all.methods(S1.N2.data, varZero = TRUE)
# save(all.S1.N2.V2, file = paste0("RData Seleccion/","all.S1.N2.V2",".RData"))

# S2 - N1 -------------------------------------------------------------------

# all.S2.N1.V2 <- parallel.all.methods(S2.N1.data)
# save(all.S2.N1.V2, file = paste0("RData Seleccion/","all.S2.N1.V2",".RData"))

# S2-N2 -------------------------------------------------------------------

# all.S2.N2.V2 <- parallel.all.methods(S2.N2.data)
# save(all.S2.N2.V2, file = paste0("RData Seleccion/","all.S2.N2.V2",".RData"))


# Paralelización RF para CV -------------------------------------------

# Parallel all RF -------------------------------------------------------------

parallel.rf <- function(dataset, 
                        seed = 123,
                        varZero = FALSE) {
  # usar availableCores() para evitar uso de muchos recursos
  n.cores <- parallelly::availableCores()
  my.cluster <- parallel::makeCluster(n.cores, 
                                      type = "PSOCK")
  # Importamos librerias o funciones para ejecutar paralelamente
  clusterEvalQ(my.cluster,library(randomForest))
  clusterEvalQ(my.cluster,library(Boruta))
  clusterEvalQ(my.cluster,library(varSelRF))
  clusterEvalQ(my.cluster,library(VSURF))
  clusterExport(my.cluster, 'randomForestCv')
  clusterExport(my.cluster, 'BorutaRandomForestCv')
  clusterExport(my.cluster, 'VarSelRfCv')
  clusterExport(my.cluster, 'VsurfCv')
  doParallel::registerDoParallel(cl = my.cluster)
  registerDoRNG()
  set.seed(seed)
  # parallel RF
  RF.res <- foreach(x = dataset , .errorhandling = "pass") %dopar% {
    randomForestCv(x = subset(x,select = -c(disease)),
                   y = x[,"disease"],
                   varZero = varZero)
  }
  # parallel Boruta
  BR.res <- foreach(x = dataset , .errorhandling = "pass") %dopar% {
    BorutaRandomForestCv(x = subset(x,select = -c(disease)),
                         y = x[,"disease"],
                         varZero = varZero)
  }
  # parallel VarselRF
  VS.res <- foreach(x = dataset , .errorhandling = "pass") %dopar% {
    VarSelRfCv(x = subset(x,select = -c(disease)),
               y = x[,"disease"],
               varZero = varZero)
  }
  # parallel VSURF
  VU.res <- foreach(x = dataset , .errorhandling = "pass") %dopar% {
    VsurfCv(x = subset(x,select = -c(disease)),
            y = x[,"disease"],
            varZero = varZero)
  }
  stopCluster(my.cluster)
  return(list("RF" = RF.res,
              "Boruta" = BR.res,
              "VarSelRF" = VS.res,
              "Vsurf" = VU.res))
}


# S1- N1 RF ---------------------------------------------------------------

# con resepcto al filtrado de las variables con baja varianza se dejara el mismo criterio 
# para random forest para comparar el metodo en si y no el filtrado previo


# start.time <- Sys.time()
# all.S1.N1.RF <- parallel.rf(S1.N1.data, varZero = TRUE)
# end.time <- Sys.time()
# save(all.S1.N1.RF,file = paste0("RData Seleccion/","all.S1.N1.RF",".RData"))
# 2.8 hrs

# S1 - N2 RF ---------------------------------------------------------------

# start.time2 <- Sys.time()
# all.S1.N2.RF <- parallel.rf(S1.N2.data,varZero = TRUE)
# end.time2 <- Sys.time()
# save(all.S1.N2.RF,file = paste0("RData Seleccion/","all.S1.N2.RF",".RData")) # 7hrs

# S2 - N1 RF ---------------------------------------------------------------

# start.time3 <- Sys.time()
# all.S2.N1.RF <- parallel.rf(S2.N1.data)
# end.time3 <- Sys.time()
# save(all.S2.N1.RF,file = paste0("RData Seleccion/","all.S2.N1.RF",".RData")) # 2 hrs


# S2 - N2 RF ---------------------------------------------------------------

# start.time4 <- Sys.time()
# all.S2.N2.RF <- parallel.rf(S2.N2.data)
# end.time4 <- Sys.time()
# save(all.S2.N2.RF,file = paste0("RData Seleccion/","all.S2.N2.RF",".RData")) # 2 hrs

