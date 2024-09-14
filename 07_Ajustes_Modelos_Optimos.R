# Extraccion de resultados con los mejores CV -----------------------------
#PLS-DA
data.models <- all.best.results %>% left_join(all.data.tibble)
data.models <- data.models %>% unite(unique_id,norm,method,bbdd,sep = "-",remove = FALSE) %>% group_by(unique_id) %>% nest()

#RF
# Funcion RF Extraccion de resultados con los mejores CV -----------------------------
data.models.rf <- all.best.results.rf %>% left_join(all.data.tibble)
data.models.rf <- data.models.rf %>% unite(unique_id,norm,method,bbdd,sep = "-",remove = FALSE) %>% group_by(unique_id) %>% nest()

# Ajuste Modelos Optimos --------------------------------------------------

library(doFuture)
plan(cluster, workers = cl)

all.data.models <- union_all(data.models, data.models.rf) %>% ungroup()

# funcion optimizada para validacion cruzada en paralelo

tic()

cv_function <- function(x) {
  x = x %>% unnest(data)
  varZero = FALSE
  if(x$norm == 'S1_N1'| x$norm == 'S1_N2') {
    varZero = TRUE
  }
  method <- match.arg(x$method, c("Vip", "PLS_DA","RC","SR","SPLS","T2","BVE","RF","Boruta","VarSelRF"))
  fitFunc <- 
    switch(
      method,
      PLS_DA = cvPlsDa.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        ncomp = x$ncomp,
        cutoff = x$cutoff ,
        varZero = varZero
      ),
      Vip = cvVip.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        ncomp = x$ncomp,
        cutoff = x$cutoff,
        vip.value = x$opt,
        varZero = varZero
      ),
      RC = cvCoeficientes.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        ncomp = x$ncomp,
        cutoff = x$cutoff,
        nvar.int = x$opt,
        varZero = varZero
      ),
      SR = cvSelectivityRatio.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        ncomp = x$ncomp,
        cutoff = x$cutoff,
        SR.val = x$opt,
        varZero = varZero
      ),
      SPLS = cvSpls.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        ncomp = x$ncomp,
        cutoff = x$cutoff,
        eta.val =  x$opt,
        varZero = varZero
      ),
      T2 = cvT2pls.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        ncomp = x$ncomp,
        cutoff = x$cutoff,
        alpha.val = x$opt,
        varZero = varZero
      ),
      BVE = cvBVE.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        ncomp = x$ncomp,
        cutoff = x$cutoff,
        vip.val  = x$opt,
        varZero = varZero
      ),
      RF = randomForestCv.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        cutoff = x$cutoff,
        varZero = varZero
      ),
      Boruta = BorutaRandomForestCv.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        cutoff = x$cutoff,
        varZero = varZero
      ),
      VarSelRF = VarSelRfCv.fixed(
        x = x$data[[1]][,-ncol(x$data[[1]])],
        y = x$data[[1]][["disease"]],
        cutoff = x$cutoff,
        varZero = varZero
      )
    )
  # solo retornar f1rep y variables
  return(fitFunc[c(1,4)])
}


plan(cluster, workers = cl)
all.iterate <- nrow(all.data.models)


handlers(global = TRUE)
progressr::with_progress({
  p <- progressor(along = seq_len(all.iterate))
  final.results <- all.data.models %>% 
    mutate(cv_results = future_map(seq_len(all.iterate),
                                   .options = furrr_options(seed = 1234),
                                   function(z){
                                     p()
                                     cv_function(all.data.models[z,])
                                   }))})


# save(final.results , file = paste0("RData Seleccion/","cloudFinalResults",".RData"))

