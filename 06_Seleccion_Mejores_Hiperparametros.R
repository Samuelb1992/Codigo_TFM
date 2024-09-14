
# Funcion selecciona mejor hiperparametros PLS-DA -------------------------


top.avg.cv <- function(x, name.col, sep = "-") {
  x %>% 
    separate(L2,name.col,sep) %>% 
    rename(f1 = value , rep = L1) %>% 
    mutate(across(name.col,as.numeric)) %>% 
    group_by(across(all_of(name.col))) %>% 
    summarise(f1mean = mean(f1) , nrep = n() ,f1sd = sd(f1)) %>% 
    ungroup() %>% 
    arrange(desc(nrep),desc(f1mean),f1sd) %>% 
    slice_head(n = 1)
}

# Estandarizaciones

# Agregamos datos de PLS sin seleccion para comparar

all.data <- list("S1_N1" = S1.N1.data,
                 "S1_N2" = S1.N2.data,
                 "S2_N1" = S2.N1.data,
                 "S2_N2" = S2.N2.data)

for (i in seq_along(all.data)) {
  if (i == 1 | i == 2) {
    temp.data <- parallel.selection(all.data[[i]],cvPlsDa, varZero = TRUE , seed = 123)$Resultados
    all.results[[i]][["PLS_DA"]] <- lapply(temp.data, function(x) x$f1rep)
  } else {
    temp.data <- parallel.selection(all.data[[i]],cvPlsDa, seed = 123)$Resultados
    all.results[[i]][["PLS_DA"]] <- lapply(temp.data, function(x) x$f1rep)
  }
}


# nombres

all.results <- lapply(all.results, function (x) {
  lapply(x , function(y) {
    names(y) <- names(DatosFinal)
    return(y)
  } )
})


## Mejor combinacion hiperparametros ---------------------------------------

best.cv.norm <- function(x){
  res.fin <- lapply(seq_along(x),function(y){
    if (y == 7) {
      res <- lapply(x[[y]],function(z) top.avg.cv(z,c("ncomp","cutoff")))
    } else {
      res <- lapply(x[[y]],function(z) top.avg.cv(z,c("ncomp","opt","cutoff")))
    }
  })
  names(res.fin) <- names(x)
  res.fin <- res.fin %>% 
    map_dfr(~bind_rows(.x,.id = "bbdd"),.id = "method")
  return(res.fin)
}

all.best.results <- all.results %>% map_dfr(~best.cv.norm(.x) , .id = "norm")



# Funcion selecciona mejor Hiperparametros Random Forest ------------------

# Funcion Mejor combinacion de hiperparametros ------------------------------------

# Mejor combinacion hiperparametros con nrep = 10

top.avg.cv.rf <- function(x, name.col, sep = "-") {
  x %>% 
    separate(L2,name.col,sep) %>% 
    rename(f1 = value , rep = L1) %>% 
    mutate(across(name.col,as.numeric)) %>% 
    group_by(across(all_of(name.col))) %>% 
    summarise(f1mean = mean(f1) , nrep = n() ,f1sd = sd(f1)) %>% 
    ungroup() %>% 
    arrange(desc(nrep),desc(f1mean),f1sd) %>% 
    slice_head(n = 1)
  # filter(nrep > 9) %>% 
  # slice(which.max(f1mean))
}

# nombres

all.results.rf <- lapply(all.results.rf, function (x) {
  lapply(x , function(y) {
    names(y) <- names(DatosFinal)
    return(y)
  } )
})

## Mejor combinacion hiperparametros ---------------------------------------

best.cv.norm.rf <- function(x){
  res.fin <- lapply(seq_along(x),function(y){          
    res <- lapply(x[[y]],function(z) top.avg.cv.rf(z,c("cutoff")))
  })
  names(res.fin) <- names(x)
  res.fin <- res.fin %>% 
    map_dfr(~bind_rows(.x,.id = "bbdd"),.id = "method")
  return(res.fin)
}

all.best.results.rf <- all.results.rf %>% map_dfr(~best.cv.norm.rf(.x) , .id = "norm")



