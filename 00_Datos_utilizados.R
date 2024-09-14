
# Cargamos Datos

load("./Data/DatosSO2.RData" , verbose = TRUE)
set.seed(123)

# Datasets ----------------------------------------------------------------

# Estandarizar nombres de datasets dentro de las listas

names(Datos_sample) <- names(DatosFinal)

# Funcion para crear lista con los dataset de cada normalizacion

N1.N2.function <- function(df.x = DatosFinal,
                           df.y = Datos_sample,
                           norm1 = "S2_NORMAL",
                           norm2 = "S2_CLR") {
  list.N1.N2 <- lapply(names(df.x), function(x) {
    N1.N2 <- list()
    N1.N2[[x]] <- data.frame(t(df.x[[x]][[norm1]][[grep(paste0("^",norm2),names(df.x[[x]][[norm1]]))]]))
    cbind(N1.N2[[x]],"disease" = df.y[[x]][["disease"]])
  })
  names(list.N1.N2) <- names(df.x)
  return(list.N1.N2)
}

# Lista con los datasets

S1.N1.data <- N1.N2.function(norm1 = "S1_NORMAL",norm2 = "S1_TSS")

S1.N2.data <- N1.N2.function(norm1 = "S1_NORMAL",norm2 = "S1_CLR")

S2.N1.data <- N1.N2.function(norm1 = "S2_NORMAL",norm2 = "S2_TSS")

S2.N2.data <- N1.N2.function(norm1 = "S2_NORMAL",norm2 = "S2_CLR")

all.data <- list("S1_N1" = S1.N1.data,
                 "S1_N2" = S1.N2.data,
                 "S2_N1" = S2.N1.data,
                 "S2_N2" = S2.N2.data)


# Datasets Tibble ---------------------------------------------------------

# Funcion para crear lista con los dataset de cada normalizacion

N1.N2.tibble <- function(df.x = DatosFinal,
                         df.y = Datos_sample,
                         norm1 = "S2_NORMAL",
                         norm2 = "S2_CLR",
                         name = "S1_S2") {
  tibble.N1.N2 <- names(df.x) %>% map_dfr(function(x) tibble("norm" = name,
                                                             "bbdd" = x,
                                                             tibble(cbind(data.frame(t(df.x[[x]][[norm1]][[grep(paste0("^", norm2), names(df.x[[x]][[norm1]]))]])),
                                                                          "disease" = df.y[[x]][["disease"]]))) %>% group_by(norm,bbdd) %>% nest() %>% ungroup())
}

# Lista con los datasets

S1.N1.tibble <- N1.N2.tibble(norm1 = "S1_NORMAL",norm2 = "S1_TSS", name = "S1_N1")

S1.N2.tibble <- N1.N2.tibble(norm1 = "S1_NORMAL",norm2 = "S1_CLR", name = "S1_N2")

S2.N1.tibble <- N1.N2.tibble(norm1 = "S2_NORMAL",norm2 = "S2_TSS", name = "S2_N1")

S2.N2.tibble <- N1.N2.tibble(norm1 = "S2_NORMAL",norm2 = "S2_CLR", name = "S2_N2")

# Juntamos todos los datos

all.data.tibble <- bind_rows(S1.N1.tibble,S1.N2.tibble,S2.N1.tibble,S2.N2.tibble)