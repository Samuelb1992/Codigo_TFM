
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
library(tidyverse)
library(rstatix)
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(ggvenn)
library(randomForest)
library(varSelRF)
library(ggpubr)

# ver paquetes cargados
# sessionInfo()
# detach alternative devtools::unload("your-package")

set.seed(123)
colores.brewer <- brewer.pal(12,"Set3")

# Parametros - Nombres --------------------------------------------------------------

bbdd.displayName <- "BBDD" # 
method.displayName <- "Técnica"
norm.displayName <- "Preprocesamiento"

# Carga Datos y Resultados ------------------------------------------------

# Datos Originales
load("/Users/samuelbeltranlastra/Documents/Proyecto R Git/TFM/Data/DatosSO2.RData" , verbose = TRUE)
# Resultados PLS - Validacion Cruzada
load("/Users/samuelbeltranlastra/Documents/Proyecto R Git/TFM/RData Seleccion/Resultados_Modelos_PLS_Final.RData", verbose = TRUE)
# Resultados RF - Validacion Cruzada
load("/Users/samuelbeltranlastra/Documents/Proyecto R Git/TFM/RData Seleccion/Resultados_Modelos_RF_Final.RData", verbose = TRUE)
# Resultados Modelos Finales con hiperparametros optimos
load("/Users/samuelbeltranlastra/Documents/Proyecto R Git/TFM/RData Seleccion/cloudFinalResults.RData", verbose = TRUE)

# Recordatorio
# all.results tiene todos los resultados en la busqueda de los hiperparametros
# all.best.results tiene los resultados de los hiperparámetros óptimos

# Transformaciones de Datasets Originales ------------------------------------

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

# Cambiamos displayName de bbdd

# DisplayName de las bases de datos

all.data.tibble$bbdd <- fct_recode(all.data.tibble$bbdd,
                                   "CRC" = "CancerColorectal",
                                   "Cirrosis" = "Cirrosis",
                                   "T2D" = "T2D",
                                   "WT2D" = "WT2D",
                                   "IBD" = "ibd",
                                   "Obesidad" = "obesidad")


# Tabla de resultados previo a mejor seleccion hiperparametros ------------

tabla.res1 <- all.results$S2_N2$Vip$Cirrosis %>% 
  separate(L2,c("ncomp","vip","cutoff"),'-') %>% 
  dplyr::rename(f1 = value , rep = L1) %>% 
  mutate(across(c("ncomp","vip","cutoff"),as.numeric)) %>% 
  group_by(across(all_of(c("ncomp","vip","cutoff")))) %>% 
  summarise(f1mean = mean(f1) , nrep = n() ,f1sd = sd(f1)) %>% 
  ungroup() %>% 
  arrange(desc(nrep),desc(f1mean),f1sd)


# OTU a SPECIES -----------------------------------------------------------

# pasamos estructura de lista a df
taxa.species <- melt(Datos_Taxa) %>% filter(Var2 == "Species") %>%  dplyr::rename(OTU.name = Var1,Type = Var2,var.specie = value,bbdd = L1)

# Estandarizamos nombres para que cruce se haga correctamente

taxa.species2 <- 
  # Nombres de bases de datos
  taxa.species %>% mutate(bbdd = case_when(
    bbdd == "IBD" ~ "ibd",
    bbdd == "Obsesidad" ~ "obesidad",
    TRUE ~ bbdd
  ))

# DisplayName de las bases de datos
taxa.species2$bbdd <- fct_recode(taxa.species2$bbdd,
                                 "CRC" = "CancerColorectal",
                                 "Cirrosis" = "Cirrosis",
                                 "T2D" = "T2D",
                                 "WT2D" = "WT2D",
                                 "IBD" = "ibd",
                                 "Obesidad" = "obesidad")

# pasamos de espacio a punto / Se ocupara posteriormente
taxa.species2$OTU.name <- gsub(" ",".",taxa.species2$OTU.name)


# Tabla de hiperparámetros óptimos ----------------------------------

all.best.results.rf.v2 <- all.best.results.rf
all.best.results.rf.v2$ncomp <- NA
all.best.results.rf.v2$opt <- NA

all.best.results.rf.pls <- rbind(all.best.results, all.best.results.rf.v2)

# Funcion para tabla de hiperparametros optimos

table.format.hop <- function(dataset_name,norm_name) {
  dataset.transform <- 
    dataset_name %>% 
    filter(norm == norm_name) %>% 
    select(method,bbdd,ncomp,cutoff,opt) %>% 
    pivot_longer(cols = c(ncomp,cutoff,opt),
                 names_to = "parametro",
                 values_to = "valor") %>% 
    pivot_wider(names_from = bbdd,
                values_from = valor) %>%
    mutate(across(where(is.numeric),as.character)) %>% 
    na.omit() %>% 
    mutate(parametro = case_when(method == "BVE" & parametro == "opt"~ "VIP",
                                 method == "RC" & parametro == "opt"~ "Nº Var",
                                 method == "SPLS" & parametro == "opt"~ "eta",
                                 method == "SR" & parametro == "opt"~ "SR",
                                 method == "T2" & parametro == "opt"~ "alpha",
                                 method == "Vip" & parametro == "opt"~ "VIP",
                                 parametro == "ncomp" ~ "Nº comp.",
                                 parametro == "cutoff" ~ "Punto de corte",
                                 TRUE ~ parametro))
  
  return(dataset.transform)
}


# Trabajo
table.best.results.rf.pls <- table.format.hop(all.best.results.rf.pls,"S2_N2")

# Anexo
table.best.results.rf.pls.s1.n1 <- table.format.hop(all.best.results.rf.pls,"S1_N1")
table.best.results.rf.pls.s1.n2 <- table.format.hop(all.best.results.rf.pls,"S1_N2")
table.best.results.rf.pls.s2.n1 <- table.format.hop(all.best.results.rf.pls,"S2_N1")


# Preprocesamiento datos para MLM´s y F1 Plots ------------------------------

data.model.final.cv.exp <- final.results %>% 
  unnest_wider(cv_results) %>% 
  select(unique_id, data, f1rep) %>% 
  unnest(data) %>% 
  unnest(f1rep) %>% 
  select(unique_id, norm, method, bbdd, f1 = value, rep = L1)

# DisplayName de las bases de datos

data.model.final.cv.exp$bbdd <- fct_recode(data.model.final.cv.exp$bbdd,
                                           "CRC" = "CancerColorectal",
                                           "Cirrosis" = "Cirrosis",
                                           "T2D" = "T2D",
                                           "WT2D" = "WT2D",
                                           "IBD" = "ibd",
                                           "Obesidad" = "obesidad")

# Busqueda de repeticiones NA

na.cases.f1 <- data.model.final.cv.exp[is.na(data.model.final.cv.exp$f1),]

# Caso de la repeticion 1,2,3 y 6 
# en la base de datos obesidad de la norma S1_N1 y T2 no se ajustaron para ningun fold
# Para no eliminar tantos datos mejor se elimina T2 (24 datos en vez de 960)

data.model.final.cv.exp.na <- data.model.final.cv.exp %>% filter(!method %in% c("T2"))


# F1 Avg Graficos ---------------------------------------------------------

# Se buscan los resultados con los hiperparametros optimos

# por los tres factores
f1.mean.avg.final <- data.model.final.cv.exp %>% 
  group_by(norm,method,bbdd) %>% 
  summarise(f1mean = mean(f1 ,na.rm = TRUE), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(across(where(is.character), as.factor))

# por los dos factores 
f1.mean.avg.final.v2 <- data.model.final.cv.exp %>% 
  group_by(method,bbdd) %>% 
  summarise(f1mean = mean(f1 , na.rm = TRUE), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(across(where(is.character), as.factor))


# Funcion grafico F1 AVG --------------------------------------------------

plot.compare <- function(x,title.mod) {
  ggplot(x, aes(x = bbdd, y = f1mean, col = method, linetype = method , group = method )) +
    geom_line(lwd = 1) +
    geom_point() +
    scale_linetype_manual(values = c(rep("solid",10)) , name = "Método") +
    scale_color_manual(values = c("green","black",brewer.pal(10,"Set3")) , name = "Método") +
    labs(x = NULL , 
         y = "F1" ,
         title = paste(title.mod)) +
    theme(axis.text.x = element_text(angle = 90),
          plot.title = element_text(hjust = 0.5))
}

# Gráfico con el promedio de todos los procesamientos

plot.f1.method.bdd <- plot.compare(f1.mean.avg.final.v2, 
                                   "Media F1-Score de los técnicas por base de datos")

# Revision de bbdd obesidad me da mas alta para los dos modelos
# plot.compare(f1.mean.avg.final.v2 %>% filter(method %in% c("PLS_DA","RF")),"test")
# ojo comparemos los hiperparametros seleccionados

# PLS DA obesidad

final.results  %>% 
  unnest(data) %>% 
  filter(method %in% c("PLS_DA","RF") & bbdd == "obesidad") %>% 
  select(norm,method,bbdd,cutoff) %>% 
  pivot_wider(names_from = norm, values_from = cutoff)


# S1-N1

f1.mean.plot.s1.n1 <- 
  f1.mean.avg.final %>% 
  filter(norm == 'S1_N1') %>% 
  plot.compare(title.mod = 'S1_N1')

f1.mean.plot.s1.n2 <- 
  f1.mean.avg.final %>% 
  filter(norm == 'S1_N2') %>% 
  plot.compare(title.mod = 'S1_N2')

f1.mean.plot.s2.n1 <- 
  f1.mean.avg.final %>% 
  filter(norm == 'S2_N1') %>% 
  plot.compare(title.mod = 'S2_N1')

f1.mean.plot.s2.n2 <- 
  f1.mean.avg.final %>% 
  filter(norm == 'S2_N2') %>% 
  plot.compare(title.mod = 'S2_N2')


# # ANOVA F1 ----------------------------------------------------------------
# 
# data.anova.f1 <- data.model.final.cv.exp.na %>% 
#   unite(bbdd_run,bbdd,rep,sep = "-",remove = FALSE) %>% 
#   mutate(across(where(is.character),as.factor)) %>% 
#   ungroup() %>% 
#   select(-unique_id)
# 
# 
# aov.res <- aov(f1 ~ method*bbdd*norm - method:bbdd:norm + Error(bbdd_run/(method*norm)), data = data.anova.f1)
# # summary(aov.res)
# aov.res.table <- broom::tidy(aov.res)
# aov.res.table.2 <- rstatix::anova_summary(aov.res)
# # Nombres explicitos
# aov.res.table.2 <- aov.res.table.2 %>% 
#   mutate(Effect = gsub("bbdd",bbdd.displayName,Effect)) %>% 
#   mutate(Effect = gsub("method",method.displayName,Effect)) %>% 
#   mutate(Effect = gsub("norm",norm.displayName,Effect))


# Tukey - Emmeans ---------------------------------------------------------

# vignette("comparisons", "emmeans") # Documentacion

# op = options()
# options(contrasts = c("contr.sum", "contr.poly"))
# 
# # Técnica
# 
# tukey.method <- emmeans(aov.res, ~method)
# pairs(tukey.method)# Tukey's
# 
# plot.tukey.method <- 
#   plot(tukey.method, comparisons = TRUE) + 
#   ggtitle("") + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   ylab("Técnica")
# 
# 
# # Norma ( para verificacion)
# 
# tukey.norm <- emmeans(aov.res, ~norm)
# pairs(tukey.norm)# Tukey's
# 
# plot.tukey.norm <- 
#   plot(tukey.norm, comparisons = TRUE) + 
#   ggtitle("") + 
#   theme(plot.title = element_text(hjust = 0.5)) +
#   ylab("Técnica")
# 
# 
# # Técnica - Norma
# 
# op = options()
# options(contrasts = c("contr.sum", "contr.poly"))
# 
# tukey.method.norm <- emmeans(aov.res, ~method:norm)
# 
# method.norm.t <- pairs(tukey.method.norm)  # Tukey's
# 
# 
# plot.tukey.method.norm <- 
#   plot(tukey.method.norm,comparisons = TRUE) + 
#   ggtitle("") + 
#   theme(
#     plot.title = element_text(hjust = 0.5) # Increase right margin of y-axis text
#   ) +
#   ylab("Técnica:Normalización")

# LMER Modelo Lineal de Efectos Mixtos ------------------------------------

# Datos modelo lineal

data.lmer <- data.model.final.cv.exp %>% 
  unite(bbdd_run,bbdd,rep,sep = "-",remove = FALSE) %>% 
  mutate(across(where(is.character),as.factor)) %>% 
  ungroup() %>% 
  select(-unique_id)


# Modelo Lineal

modelo.lineal <- lmer(f1 ~ norm * method + norm * bbdd + method * bbdd + (1 | bbdd_run) , data = data.lmer)

# Resultados

tabla.modelo.lineal <- anova(modelo.lineal, type = 3)

rownames(tabla.modelo.lineal) <-  c(norm.displayName,
                                    method.displayName,
                                    bbdd.displayName,
                                    paste0(norm.displayName,":",method.displayName),
                                    paste0(norm.displayName,":",bbdd.displayName),
                                    paste0(method.displayName,":",bbdd.displayName))


# Tukey Modelo Lineal -----------------------------------------------------

# vignette("comparisons", "emmeans") # Documentacion
# emmeans , comparaciones por pares ajustados por Tukey por defecto
# emmip grafico de interaccion (interaccion plot)

order.factor.plot <- rev(c("PLS_DA", "RC", "SR", "T2","Vip","SPLS","BVE","RF","VarSelRF","Boruta"))

op = options()
options(contrasts = c("contr.sum", "contr.poly"))


## Efectos Simples ---------------------------------------------------------

# Técnica

lineal.tukey.method <- emmeans(modelo.lineal, ~method)
pairs(lineal.tukey.method)


lineal.plot.tukey.method <- 
  plot(lineal.tukey.method, comparisons = TRUE) # grafico estandar


# Datos Grafico personalizado
data.plot.tukey.method <- as.data.frame(lineal.plot.tukey.method$data)
data.plot.tukey.method$method <- factor(data.plot.tukey.method$method,levels = order.factor.plot)
data.plot.tukey.method$method_numeric <- as.numeric(data.plot.tukey.method$method)
# agregar categoria del metodo
data.plot.tukey.method <- data.plot.tukey.method %>%
  mutate(category =   case_when(method %in% c("PLS_DA","RC","SR","T2","Vip","SPLS","BVE") ~ "PLS-DA",
                                method %in% c("Boruta","VarSelRF","RF") ~ "RF",
                                TRUE ~ "otro"))

# Plot  Grafico personalizado
lineal.plot.tukey.method.V2 <- 
  ggplot(data.plot.tukey.method) +
  # Rectangulo CI
  geom_rect(aes(xmin = lcl, xmax = ucl, ymin = as.numeric(method) - 0.10, ymax = as.numeric(method) + 0.10),
            fill = "#bcbcf3", alpha = 0.8) +
  
  # Flechas para comparar dif.sig.
  geom_segment(aes(x = the.emmean, xend = lcmpl, y = as.numeric(method), yend = as.numeric(method)),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "red", size = 0.4) +
  geom_segment(aes(x = the.emmean, xend = rcmpl, y = as.numeric(method), yend = as.numeric(method)),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "red", size = 0.4) +
  
  # punto media
  geom_point(aes(x = the.emmean, y = as.numeric(method)), color = "black", size = 2) +
  # Labels como factor
  scale_y_continuous(breaks = data.plot.tukey.method$method_numeric, labels = levels(data.plot.tukey.method$pri.fac)) +
  # Label grupo
  facet_grid(category ~ ., scales = "free_y", space = "free_y")+
  labs(x = "emmean", y = "Técnica")


# BBDD

lineal.tukey.bbdd <- emmeans(modelo.lineal, ~bbdd)
# pairs(lineal.tukey.bbdd)

lineal.plot.tukey.bbdd <- 
  plot(lineal.tukey.bbdd, comparison = TRUE)+
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(paste0(bbdd.displayName))

# data
# as.data.frame(lineal.plot.tukey.bbdd$data)

# Norma

lineal.tukey.norm <- emmeans(modelo.lineal, ~norm)
# pairs(lineal.tukey.norm)

lineal.plot.tukey.norm <- 
  plot(lineal.tukey.norm, comparisons = TRUE) + 
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(paste0(norm.displayName))


## Interacciones -----------------------------------------------------------

### Tecnica BBDD ------------------------------------------------------------

#### Modelo - Pares - Pares Seccionado ---------------------------------------

lineal.tukey.method.bbdd <- emmeans(modelo.lineal, ~method:bbdd)
# comparaciones pairs(lineal.tukey.method.bbdd)
# pairs(lineal.tukey.method.bbdd) %>% as.data.frame() %>% filter(grepl("\\bT2D\\b.*\\bT2D\\b", contrast) & p.value < 0.1)
# GRAFICO COMPLETO DE COMPARACIONES
lineal.plot.tukey.method.bbdd <- 
  plot(lineal.tukey.method.bbdd,comparisons = TRUE) + 
  ggtitle("") + 
  theme(
    plot.title = element_text(hjust = 0.5) # Increase right margin of y-axis text
  ) +
  ylab(paste0(method.displayName,":",bbdd.displayName))


#GRAFICO SECCIONADO

# FUNCION DATOS
data.tukey.seccionado.interaccion <- function(datos_emmean,
                                              method_cat = TRUE,
                                              bbdd_cat = TRUE,
                                              norm_cat = FALSE,
                                              filtros = list(
                                                "method" = c(),
                                                "bbdd" = c(),
                                                "norm" = c()
                                              )) {
  # datos del emmean
  data.plot <- as.data.frame(datos_emmean$data)
  # reordenamos factores
  data.plot$method <- factor(data.plot$method,levels = order.factor.plot)
  data.plot$method_numeric <- as.numeric(data.plot$method)
  # agregar categoria del metodo
  if (method_cat) {
    
    data.plot <- data.plot %>%
      mutate(category =   case_when(method %in% c("PLS_DA","RC","SR","T2","Vip","SPLS","BVE") ~ "PLS-DA",
                                    method %in% c("Boruta","VarSelRF","RF") ~ "RF",
                                    TRUE ~ "otro"))
  }
  # filtramos por tecnica , preprocesamiento o enfermedad para seccionar grafico
  for (columna in names(filtros)) {
    if (columna %in% names(data.plot)) {
      vector.filtro <- filtros[[columna]]
      # veo si esta vacio
      if (!is.null(vector.filtro)) {
        data.plot <- data.plot[data.plot[[columna]] %in% vector.filtro,]
      }
    }
  }
  # actualizo niveles y order factor plot
  data.plot$pri.fac <- droplevels(data.plot$pri.fac)
  # filtro metodos en caso de que se filtre alguno para la creacion de niveles
  if (!is.null(filtros[["method"]])) {
    vector.method <- filtros[["method"]]
    order.factor.plot <- intersect(order.factor.plot,vector.method)
  }
  # solo metodo
  if (method_cat & !bbdd_cat & !norm_cat) {
    data.plot$pri.fac <- factor(data.plot$pri.fac, levels = intersect(order.factor.plot,levels(data.plot$pri.fac)))
  }
  # interaccion
  if (method_cat & bbdd_cat) {
    data.plot$pri.fac <- factor(data.plot$pri.fac,
                                levels = paste(rep(order.factor.plot,length(unique(data.plot$bbdd))),
                                               rep(unique(data.plot$bbdd),each = length(order.factor.plot))))
  }
  # interaccion
  if (method_cat & norm_cat) {
    data.plot$pri.fac <- factor(data.plot$pri.fac,
                                levels = paste(rep(order.factor.plot,length(unique(data.plot$norm))),
                                               rep(unique(data.plot$norm),each = length(order.factor.plot))))
  }
  return(data.plot)
}

# FUNCION GRAFICO
plot.interaction.tukey <- function(data,
                                   method_cat = TRUE) {
  plot.df <- 
    data %>%
    ggplot()+
    # Rectangulo CI
    geom_rect(aes(xmin = lcl, xmax = ucl, ymin = as.numeric(pri.fac) - 0.15, ymax = as.numeric(pri.fac) + 0.15),
              fill = "#bcbcf3") +
    
    # Flechas para comparar dif.sig.
    geom_segment(aes(x = the.emmean, xend = lcmpl, y = as.numeric(pri.fac), yend = as.numeric(pri.fac)),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 color = "red", size = 0.4) +
    geom_segment(aes(x = the.emmean, xend = rcmpl, y = as.numeric(pri.fac), yend = as.numeric(pri.fac)),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 color = "red", size = 0.4) +
    # punto media
    geom_point(aes(x = the.emmean, y = as.numeric(pri.fac)), color = "black", size = 2)
  
  
  if (method_cat) {
    plot.df <- plot.df + scale_y_continuous(breaks = 1:length(levels(data$pri.fac)), labels = levels(data$pri.fac))
  }
  return(plot.df)
}

# Custom Completo

lineal.plot.tukey.method.bbdd.v2 <- plot.interaction.tukey(
  data.tukey.seccionado.interaccion(lineal.plot.tukey.method.bbdd, 
                                    filtros = list(
                                      "method" = c(),
                                      "bbdd" = c(),
                                      "norm" = c()
                                    ))
)+ylab(paste0(method.displayName,":",bbdd.displayName))+
  xlab("emmean")+
  facet_grid(bbdd~.,scales = "free_y",space = "free_y")

# CRC / Cirrosis
method.bbdd.secc.v1 <-  plot.interaction.tukey(
  data.tukey.seccionado.interaccion(lineal.plot.tukey.method.bbdd, 
                                    filtros = list(
                                      "method" = c(),
                                      "bbdd" = c("CRC","Cirrosis"),
                                      "norm" = c()
                                    ))
)

method.bbdd.secc.v1 <- method.bbdd.secc.v1 + labs(x = "emmean", y = paste0(method.displayName,":",bbdd.displayName))

# IBD / T2D
method.bbdd.secc.v2 <-  plot.interaction.tukey(
  data.tukey.seccionado.interaccion(lineal.plot.tukey.method.bbdd, 
                                    filtros = list(
                                      "method" = c(),
                                      "bbdd" = c("IBD","T2D"),
                                      "norm" = c()
                                    ))
) + labs(x = "emmean", y = paste0(method.displayName,":",bbdd.displayName))


# Obesidad/ WT2D
method.bbdd.secc.v3 <-  plot.interaction.tukey(
  data.tukey.seccionado.interaccion(lineal.plot.tukey.method.bbdd, 
                                    filtros = list(
                                      "method" = c(),
                                      "bbdd" = c("Obesidad","WT2D"),
                                      "norm" = c()
                                    ))
) + labs(x = "emmean", y = paste0(method.displayName,":",bbdd.displayName))


#### INTERACCION PLOT ---------------------------------------

# Plot personalizado interaccion
data.lineal.tukey.method.bbdd <- emmip(lineal.tukey.method.bbdd, bbdd ~ method, plotit = FALSE)
data.lineal.tukey.method.bbdd$method <- factor(data.lineal.tukey.method.bbdd$method,levels = rev(order.factor.plot))
data.lineal.tukey.method.bbdd$method_numeric <- as.numeric(data.lineal.tukey.method.bbdd$method)
# agregar categoria del metodo
data.lineal.tukey.method.bbdd <- data.lineal.tukey.method.bbdd %>%
  mutate(category =   case_when(method %in% c("PLS_DA","RC","SR","T2","Vip","SPLS","BVE") ~ "PLS-DA",
                                method %in% c("Boruta","VarSelRF","RF") ~ "RF",
                                TRUE ~ "otro"))


plot.data.lineal.tukey.method.bbdd.v2 <- 
  data.lineal.tukey.method.bbdd %>% 
  ggplot(aes(x = method, y = yvar , group = bbdd , color = bbdd)) +
  geom_line()+
  geom_point()+
  # label1
  annotate("text", x = 4, y = 0.9 , label = "PLS-DA")+
  # label2
  annotate("text", x = 9, y = 0.9 , label = "RF")+
  # line
  geom_vline(xintercept = 7, linetype = "dashed", alpha = 0.3)+
  labs(x = method.displayName,
       y = "F1-score medio",
       colour = bbdd.displayName)+
  theme(axis.text.x = element_text(size = 7))

### Tecnica - NORMA ------------------------------------------------------------


#### Modelo - Pares - Pares Seccionado ---------------------------------------

op = options()
options(contrasts = c("contr.sum", "contr.poly"))

# Emmeans
lineal.tukey.method.norm <- emmeans(modelo.lineal, ~method:norm)
# emmip(lineal.tukey.method.norm, norm ~ method)

# Plot Pairwise Emmeans
lineal.plot.tukey.method.norm <- 
  plot(lineal.tukey.method.norm,comparisons = TRUE) + 
  ggtitle("") + 
  theme(
    plot.title = element_text(hjust = 0.5) # Increase right margin of y-axis text
  ) +
  ylab(paste0(method.displayName,":",norm.displayName))


lineal.plot.tukey.method.norm.v2 <- 
  plot.interaction.tukey(
    data.tukey.seccionado.interaccion(lineal.plot.tukey.method.norm,
                                      bbdd_cat = FALSE,
                                      norm_cat = TRUE,
                                      filtros = list(
                                        "method" = c(),
                                        "bbdd" = c(),
                                        "norm" = c()
                                      ))
  ) + labs(x = "emmean", y = paste0(method.displayName,":",norm.displayName))+
  facet_grid(norm~.,scale = "free_y", space = "free_y")


# Plot personalizado interaccion
data.lineal.tukey.method.norm <- emmip(lineal.tukey.method.norm, norm ~ method, plotit = FALSE)
data.lineal.tukey.method.norm$method <- factor(data.lineal.tukey.method.norm$method,levels = rev(order.factor.plot))
data.lineal.tukey.method.norm$method_numeric <- as.numeric(data.lineal.tukey.method.norm$method)
# agregar categoria del metodo
data.lineal.tukey.method.norm <- data.lineal.tukey.method.norm %>%
  mutate(category =   case_when(method %in% c("PLS_DA","RC","SR","T2","Vip","SPLS","BVE") ~ "PLS-DA",
                                method %in% c("Boruta","VarSelRF","RF") ~ "RF",
                                TRUE ~ "otro"))


plot.data.lineal.tukey.method.norm.v2 <- 
  data.lineal.tukey.method.norm %>% 
  ggplot(aes(x = method, y = yvar , group = norm , color = norm)) +
  geom_line()+
  geom_point()+
  # label1
  annotate("text", x = 4, y = 0.8 , label = "PLS-DA")+
  # label2
  annotate("text", x = 9, y = 0.8 , label = "RF")+
  # line
  geom_vline(xintercept = 7, linetype = "dashed", alpha = 0.3)+
  labs(x = method.displayName,
       y = "F1-score medio",
       colour = norm.displayName)+
  theme(axis.text.x = element_text(size = 7))


#### Seccionado --------------------------------------------------------------

# S1-N1 / S1-N2
method.norm.secc.v1 <- plot.interaction.tukey(
  data.tukey.seccionado.interaccion(lineal.plot.tukey.method.norm,
                                    bbdd_cat = FALSE,
                                    norm_cat = TRUE,
                                    filtros = list(
                                      "method" = c(),
                                      "bbdd" = c(),
                                      "norm" = c("S1_N1","S1_N2")
                                    ))
) + labs(x = "emmean", y = paste0(method.displayName,":",norm.displayName))

# S2-N1 / S2-N2
method.norm.secc.v2 <- plot.interaction.tukey(
  data.tukey.seccionado.interaccion(lineal.plot.tukey.method.norm,
                                    bbdd_cat = FALSE,
                                    norm_cat = TRUE,
                                    filtros = list(
                                      "method" = c(),
                                      "bbdd" = c(),
                                      "norm" = c("S2_N1","S2_N2")
                                    ))
) + labs(x = "emmean", y = paste0(method.displayName,":",norm.displayName))


# RF - Varsel - Boruta - BVE

# method.norm.secc.v3 <- plot.interaction.tukey(
#   data.tukey.seccionado.interaccion(lineal.plot.tukey.method.norm,
#                                     bbdd_cat = FALSE,
#                                     norm_cat = TRUE,
#                                     filtros = list(
#                                       "method" = c("Boruta","VarSelRF","RF","BVE"),
#                                       "bbdd" = c(),
#                                       "norm" = c()
#                                     ))
# ) + labs(x = "emmean", y = paste0(method.displayName,":",norm.displayName))

# Por separado


### NORMA BBDD ------------------------------------------------------------

#### Modelo - Pares - Pares Seccionado ------------------------------------

op = options()
options(contrasts = c("contr.sum", "contr.poly"))

# Emmeans
lineal.tukey.bbdd.norm <- emmeans(modelo.lineal, ~norm:bbdd)
# emmip(lineal.tukey.bbdd.norm, norm ~ bbdd)
# df <- as.data.frame(pairs(lineal.tukey.bbdd.norm))

# Plot Pairwise Emmeans
lineal.plot.tukey.bbdd.norm <- 
  plot(lineal.tukey.bbdd.norm,comparisons = TRUE) + 
  ggtitle("") + 
  theme(
    plot.title = element_text(hjust = 0.5) # Increase right margin of y-axis text
  ) +
  ylab(paste0(bbdd.displayName,":",norm.displayName))+
  facet_grid(bbdd ~. , scales = "free_y" , space = "free_y")


# Estabilidad Tanimoto ----------------------------------------------------

# Tanimonoto distance metric para medir la similaridad entre dos conjuntos. En este caso se comparan los folds
# 1 Crear los vectores de 0/1 , 0 significa que el feature no esta y 1 significa que esta

# Conjunto variables PLS y Random Forest

data.model.final.variables <- 
  final.results %>% 
  unnest_wider(cv_results) %>%
  unnest(c(var,f1rep)) %>%
  select(unique_id, data, var ,rep = L1)

data.model.final.variables.exp <- 
  data.model.final.variables %>% 
  unnest(var)


# Transformacion de datos

final.var.selected <- 
  data.model.final.variables.exp %>%
  unnest(data) %>% 
  select(norm, method, bbdd, var , rep) %>%
  group_by(norm, method, bbdd , rep) %>%
  mutate(repfold = n() * 10) %>%
  mutate(n_fold = row_number()) %>% 
  mutate(unique_id = as.factor(paste(norm,method,bbdd,sep = "-"))) %>% 
  ungroup()


# Lista de vectores con todos los features de cada dataset

all.var <- lapply(all.data, function(x) {
  lapply(x, function(y)
    names(y[, -ncol(y)]))
})

all.var.tibble <- melt(all.var) %>%
  group_by(L1, L2) %>%
  summarise(all_var = list(value)) %>%
  dplyr::rename(norm = L1 , bbdd = L2)

# Union dataset. Dataset que contendra las variables seleccionadas por fold y el total de variables para asi realizar comparacion

tan.var.all <- inner_join(final.var.selected,all.var.tibble)

# DisplayName de las bases de datos
tan.var.all$bbdd <- fct_recode(tan.var.all$bbdd,
                               "CRC" = "CancerColorectal",
                               "Cirrosis" = "Cirrosis",
                               "T2D" = "T2D",
                               "WT2D" = "WT2D",
                               "IBD" = "ibd",
                               "Obesidad" = "obesidad")

# Creamos vector de 0/1 comparando features seleccionados vs todos los features de cada bbdd

tan.var.df <-
  tan.var.all %>% 
  mutate(tan_var = map2(var, all_var, ~ .y %in% .x))

# alternativa ~ as.integer(.y %in% .x)

# Creamos dataframe donde cada fila es un par de folds (ej fold-1 /fold-2... fold-1 / fold -3).. n(n-1)/2 pero dentro de la repeticion
# 10 folds -> combinacion de 10 sobre 2 , serian 45 por rep.  45 * 10 rep * 240 modelos = 108.000 datos aprox

tan.var.cross <- 
  tan.var.df %>% 
  full_join(tan.var.df,by = c("norm","method","bbdd","rep")) %>% 
  select(norm,method,bbdd,rep,unique_id.x,n_fold.x,tan_var.x,n_fold.y,tan_var.y) %>% 
  group_by(unique_id.x,rep) %>% 
  # eliminamos comparaciones equivalentes 1-2 == 2-1
  filter(!duplicated(paste0(pmax(n_fold.x,n_fold.y),pmin(n_fold.x,n_fold.y)))) %>% 
  # eliminamos comparacion de vector consigomismo
  filter(n_fold.x != n_fold.y) %>% 
  ungroup()

# Con el vector tan_var.x y tan_var-y (fold1 /fold2 ) calculamos el la distancia de Tanimoto
# DOI 10.1007/s10115-006-0040-8

tan.distance <- function(x,y) {
  1-((length(x) + length(y) - 2*(sum(x == y)))/(length(x) + length(y) - (sum(x == y))))
}

# Distancia de tanimoto para cada par de folds

tan.var.res <- tan.var.cross %>%
  mutate(tan_dis = map2_dbl(tan_var.x, tan_var.y, ~ tan.distance(.x,.y)))

# Se obtiene la estabilidad media para cada repeticion de cada conjunto de datos

tan.var.res.v2 <- tan.var.res %>% 
  group_by(norm,method,bbdd,rep) %>% 
  summarise(mean_tan_dis = mean(tan_dis))

# Modelo Lineal Estabilidad -----------------------------------------------

datos.estabilidad <- 
  tan.var.res.v2 %>% 
  ungroup() %>% 
  unite(bbdd_run,bbdd,rep,sep = "-",remove = FALSE) %>% 
  mutate(across(where(is.character),as.factor)) %>% 
  ungroup() %>% 
  select(-rep)

# Se filtran metodos PLS-DA y RF porque no incluyen seleccion de variables
datos.estabilidad <- datos.estabilidad %>% filter(!method %in% c("PLS_DA","RF"))

# Efecto aleatorio provoca singularidad asi que se elimina y se calcula solo con factores fijos
# lmer(mean_tan_dis ~ norm * method + norm * bbdd + method * bbdd + (1 | bbdd_run) , data = datos.estabilidad)
# VarCorr(model.test) .. var de bbdd_run cero...

modelo.lineal.estabilidad <- lm(mean_tan_dis ~ norm * method + norm * bbdd + method * bbdd , data = datos.estabilidad)

# Resultados

tabla.modelo.lineal.estabilidad <- anova(modelo.lineal.estabilidad)


# summary(aov.res)
aov.est.table.2 <- rstatix::anova_summary(tabla.modelo.lineal.estabilidad)

# Nombres explicitos
aov.est.table.2 <- aov.est.table.2 %>% 
  mutate(Effect = gsub("bbdd",bbdd.displayName,Effect)) %>% 
  mutate(Effect = gsub("method",method.displayName,Effect)) %>% 
  mutate(Effect = gsub("norm",norm.displayName,Effect))


# rownames(tabla.modelo.lineal.estabilidad) <-  c(norm.displayName,
#                                                 method.displayName,
#                                                 bbdd.displayName,
#                                                 paste0(norm.displayName,":",method.displayName),
#                                                 paste0(norm.displayName,":",bbdd.displayName),
#                                                 paste0(method.displayName,":",bbdd.displayName))


# Resultados

tabla.modelo.lineal <- anova(modelo.lineal, type = 3)

rownames(tabla.modelo.lineal) <-  c(norm.displayName,
                                    method.displayName,
                                    bbdd.displayName,
                                    paste0(norm.displayName,":",method.displayName),
                                    paste0(norm.displayName,":",bbdd.displayName),
                                    paste0(method.displayName,":",bbdd.displayName))




# Tukey Modelo Lineal Estabilidad -----------------------------------------------------

# vignette("comparisons", "emmeans") # Documentacion

op = options()
options(contrasts = c("contr.sum", "contr.poly"))


## Efectos Simples ---------------------------------------------------------

### Técnica -----------------------------------------------------------------

est.lineal.tukey.method <- emmeans(modelo.lineal.estabilidad, ~method)
pairs(est.lineal.tukey.method)

# Tukey defecto
est.lineal.plot.tukey.method <- 
  plot(est.lineal.tukey.method, comparisons = TRUE) + 
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(method.displayName)


# Tukey estabilidad tecnica mejorado
est.lineal.plot.tukey.method.v2 <- 
  plot.interaction.tukey(data.tukey.seccionado.interaccion(est.lineal.plot.tukey.method,bbdd_cat = FALSE, norm_cat = FALSE))+
  facet_grid(category~.,scales = "free_y",space = "free_y")+
  labs(x = "emmean",y = method.displayName)


### Norma ---------------------------------------------------------------

est.lineal.tukey.norm <- emmeans(modelo.lineal.estabilidad, ~norm)
# pairs(est.lineal.tukey.method)

est.lineal.plot.tukey.norm <- 
  plot(est.lineal.tukey.norm, comparisons = TRUE) + 
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(norm.displayName)

### BBDD ---------------------------------------------------------------

est.lineal.tukey.bbdd <- emmeans(modelo.lineal.estabilidad, ~bbdd)
# pairs(est.lineal.tukey.bbdd)

est.lineal.plot.tukey.bbdd <- 
  plot(est.lineal.tukey.bbdd, comparisons = TRUE) + 
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab(bbdd.displayName)

## INTERACCIONES ---------------------------------------------------------------

### TECNICA - BBDD  ------------------------------------------------------------

est.lineal.tukey.method.bbdd <- emmeans(modelo.lineal.estabilidad, ~method:bbdd)

est.lineal.plot.tukey.method.bbdd <- 
  plot(est.lineal.tukey.method.bbdd, comparisons = TRUE)


est.lineal.plot.tukey.method.bbdd.v2 <- plot.interaction.tukey(
  data.tukey.seccionado.interaccion(est.lineal.plot.tukey.method.bbdd,
                                    bbdd_cat = TRUE,
                                    norm_cat = FALSE,
                                    filtros = list(
                                      "method" = c("RC","SR","T2","Vip","SPLS","BVE","VarSelRF","Boruta"),
                                      "bbdd" = c(),
                                      "norm" = c())))+ 
  labs(x = "emmean", y = paste0(method.displayName,":",bbdd.displayName))+
  facet_grid(bbdd~.,scales="free_y",space = "free_y")




#### INTERACCION PLOT --------------------------------------------------------

# Plot personalizado interaccion
est.int.method.bbdd <- emmip(est.lineal.tukey.method.bbdd, bbdd ~ method, plotit = FALSE)
est.int.method.bbdd$method <- factor(est.int.method.bbdd$method,levels = rev(order.factor.plot))
est.int.method.bbdd$method_numeric <- as.numeric(est.int.method.bbdd$method)
# agregar categoria del metodo
est.int.method.bbdd <- est.int.method.bbdd %>%
  mutate(category =   case_when(method %in% c("PLS_DA","RC","SR","T2","Vip","SPLS","BVE") ~ "PLS-DA",
                                method %in% c("Boruta","VarSelRF","RF") ~ "RF",
                                TRUE ~ "otro"))


plot.est.int.method.bbdd.v2 <- 
  est.int.method.bbdd %>% 
  ggplot(aes(x = method, y = yvar , group = bbdd , color = bbdd)) +
  geom_line()+
  geom_point()+
  # label1
  annotate("text", x = 4, y = 1 , label = "PLS-DA")+
  # label2
  annotate("text", x = 8, y = 1 , label = "RF")+
  # line
  geom_vline(xintercept = 6, linetype = "dashed", alpha = 0.3)+
  labs(x = method.displayName,
       y = "Estabilidad media",
       colour = bbdd.displayName)+
  theme(axis.text.x = element_text(size = 7))


### TECNICA - NORMA ------------------------------------------------------------

est.lineal.tukey.method.norm <- emmeans(modelo.lineal.estabilidad, ~method:norm)
# pairs(est.lineal.tukey.method.norm)
# pairs(est.lineal.tukey.method.norm) %>% as.data.frame() %>% filter(contrast == "Boruta S2_N2 - SPLS S2_N2")

est.lineal.plot.tukey.method.norm <- plot(est.lineal.tukey.method.norm, comparisons = TRUE)


est.lineal.plot.tukey.method.norm.v2 <- plot.interaction.tukey(
  data.tukey.seccionado.interaccion(est.lineal.plot.tukey.method.norm,
                                    bbdd_cat = FALSE,
                                    norm_cat = TRUE,
                                    filtros = list(
                                      "method" = c("RC","SR","T2","Vip","SPLS","BVE","VarSelRF","Boruta"),
                                      "bbdd" = c(),
                                      "norm" = c())))+ 
  labs(x = "emmean", y = paste0(method.displayName,":",norm.displayName))+
  facet_grid(norm~.,scales="free_y",space = "free_y")



#### INTERACCION PLOT --------------------------------------------------------

# Plot personalizado interaccion
est.int.method.norm <- emmip(est.lineal.tukey.method.norm, norm ~ method, plotit = FALSE)
est.int.method.norm$method <- factor(est.int.method.norm$method,levels = rev(order.factor.plot))
est.int.method.norm$method_numeric <- as.numeric(est.int.method.norm$method)
# agregar categoria del metodo
est.int.method.norm <- est.int.method.norm %>%
  mutate(category =   case_when(method %in% c("PLS_DA","RC","SR","T2","Vip","SPLS","BVE") ~ "PLS-DA",
                                method %in% c("Boruta","VarSelRF","RF") ~ "RF",
                                TRUE ~ "otro"))


plot.est.int.method.norm.v2 <- 
  est.int.method.norm %>% 
  ggplot(aes(x = method, y = yvar , group = norm , color = norm)) +
  geom_line()+
  geom_point()+
  # label1
  annotate("text", x = 4, y = 1 , label = "PLS-DA")+
  # label2
  annotate("text", x = 8, y = 1 , label = "RF")+
  # line
  geom_vline(xintercept = 6, linetype = "dashed", alpha = 0.3)+
  labs(x = method.displayName,
       y = "Estabilidad media",
       colour = norm.displayName)+
  theme(axis.text.x = element_text(size = 7))

### BBDD - NORMA ------------------------------------------------------------

est.lineal.tukey.bbdd.norm <- emmeans(modelo.lineal.estabilidad, ~norm:bbdd)
# pairs(est.lineal.tukey.method.norm)
# pairs(est.lineal.tukey.method.norm) %>% as.data.frame() %>% filter(contrast == "Boruta S2_N2 - SPLS S2_N2")

est.lineal.plot.tukey.bbdd.norm <- plot(est.lineal.tukey.bbdd.norm, comparisons = TRUE)+
  facet_grid(bbdd~.,scales = "free_y", space = "free_y")+
  labs(y = paste0(norm.displayName,":",bbdd.displayName))



# Grafico F1 AVG y Estabilidad Media --------------------------------------

est.mean.avg.final <- 
  datos.estabilidad %>% 
  group_by(norm,method,bbdd) %>% 
  summarise(mean_tan_dis = mean(mean_tan_dis, na.rm = TRUE))

# se aplica filtro de RF y PLS-DA para asi ver cuales son las mejores tecnicas de seleccion

f1.estabilidad.datos <- inner_join(f1.mean.avg.final, est.mean.avg.final) %>% 
  filter(!method %in% c("PLS_DA","RF"))

# Grafico general solo metodo

data.plot.avg.f1.est.v1 <- 
  f1.estabilidad.datos %>% 
  group_by(method) %>% 
  summarise(f1mean = mean(f1mean, na.rm = TRUE)
            ,mean_tan_dist = mean(mean_tan_dis, na.rm = TRUE))

# Grafico general metodo y preprocesamiento en grupo

data.plot.avg.f1.est.v2 <- 
  f1.estabilidad.datos %>% 
  group_by(method,norm) %>% 
  summarise(f1mean = mean(f1mean, na.rm = TRUE)
            ,mean_tan_dist = mean(mean_tan_dis, na.rm = TRUE)) %>% 
  ungroup()


library(ggplot2)
library(ggrepel)

# Medias para scatterplot
mean.scatterplot.f1 <- mean(data.plot.avg.f1.est.v2$f1mean)
mean.scatterplot.est <- mean(data.plot.avg.f1.est.v2$mean_tan_dist)

# Plot
plot.est.f1.methods <- 
  ggplot(data.plot.avg.f1.est.v2, aes(x = mean_tan_dist, y = f1mean, shape = norm, color = norm)) +
  geom_point() +
  geom_text_repel(aes(label = method), show.legend = F) +
  labs(
    x = "Media Estabilidad",
    y = "Media F1-Score",
    title = "Media F1-Score vs Media Estabilidad"
  ) +
  scale_shape_manual(values = c(16, 17, 18,19), name = "Preprocesamiento") +
  scale_color_manual(values = colores.brewer[4:7], name = "Preprocesamiento")+
  theme_bw()+
  # Linea Medias
  geom_vline(xintercept = mean.scatterplot.est, linetype = "dashed", color = "black",alpha = 0.4)+
  geom_hline(yintercept = mean.scatterplot.f1, linetype = "dashed", color = "black",alpha = 0.4)


# Interpretacion ----------------------------------------------------------

emmean.obesidad.s2_n1 <- emmeans(modelo.lineal, ~ method | norm * bbdd, 
                                 at = list(norm = "S2_N1", bbdd = "Obesidad"))

pairs.obesidad.s2_n1 <- pairs(emmean.obesidad.s2_n1)
plot.pairs.obesidad.s2_n1 <- plot(emmean.obesidad.s2_n1, comparisons =TRUE)


# Numero de Variables -----------------------------------------------------

# Creamos columna con el número de variables seleccionada por fold

df.num.variables <- 
  final.var.selected %>% 
  mutate(len_col = map_int(var,length)) %>% 
  ungroup()


# DisplayName de las bases de datos
df.num.variables$bbdd <- fct_recode(df.num.variables$bbdd,
                                    "CRC" = "CancerColorectal",
                                    "Cirrosis" = "Cirrosis",
                                    "T2D" = "T2D",
                                    "WT2D" = "WT2D",
                                    "IBD" = "ibd",
                                    "Obesidad" = "obesidad")

# Creamos funcion que calcule por unique_id algunas estadistincas

eda.variables <- function(x) {
  x = na.omit(x)
  min_0 = min(x)
  quantile_25 = quantile(x,0.25)
  median_50 = median(x)
  quantile_75 = quantile(x,0.75)
  max_100 = max(x)
  return(list(
    "Minimo" =min_0,
    "1er Cuantil" = quantile_25,
    "Mediana" = median_50,
    "3er Cuantil" = quantile_75,
    "Máximo" = max_100
  ))
}

# Aplicamos funcion para obtener numero de variables seleccionadas

df.num.variables.v2 <- 
  df.num.variables %>% 
  group_by(norm,method,bbdd) %>% 
  summarise(n_select_eda = list(eda.variables(len_col))) %>% 
  unnest_wider(n_select_eda)

# Numero de filas y columnas de cada bbdd

all.data.tibble.v2 <- 
  all.data.tibble %>% 
  mutate(nrow = map_int(data, function(x) dim(x)[1]),
         ncol = map_int(data, function(x) dim(x)[2] - 1))

## Porcentaje Seleccion de la variable  ------------------------------------------------

df.per.selected.variables <- 
  final.var.selected %>% 
  unnest(var) %>% 
  ungroup() %>% 
  group_by(unique_id,norm,method,bbdd,var) %>% 
  # % de veces seleccionada
  summarise(per_var = n() / max(repfold)) %>%
  arrange(unique_id,desc(per_var)) %>% 
  # ranking descendente
  group_by(unique_id,norm,method,bbdd) %>% 
  mutate(rank_var = row_number(desc(per_var))) %>% 
  ungroup()

# DisplayName de las bases de datos
df.per.selected.variables$bbdd <- fct_recode(df.per.selected.variables$bbdd,
                                             "CRC" = "CancerColorectal",
                                             "Cirrosis" = "Cirrosis",
                                             "T2D" = "T2D",
                                             "WT2D" = "WT2D",
                                             "IBD" = "ibd",
                                             "Obesidad" = "obesidad")


# ahora le pego las otras estadisticas y veo el corte

df.per.selected.variables.v2 <- left_join(df.per.selected.variables, df.num.variables.v2, by = c("norm","method","bbdd"))

# ind de seleccion 
# porcentaje 50%, para tener biomarcadores mas robustos
# tambien se crea columna con el 1er cuartil..similar

percent.selected <- 0.5

df.per.selected.variables.v3 <- 
  df.per.selected.variables.v2 %>% 
  mutate(ind_var_selected_porcentaje = ifelse(per_var >= percent.selected, 1 , 0),
         ind_var_cuartil_1 = ifelse(rank_var <= `1er Cuantil`,1,0))

# OTU a SPECIES

df.per.selected.variables.v4 <- left_join(df.per.selected.variables.v3, taxa.species2, by = c("var" = "OTU.name","bbdd" = "bbdd"))

# Variables seleccionadas finalmente

df.per.selected.variables.v5 <- 
  df.per.selected.variables.v4 %>% 
  select(unique_id,norm,method,bbdd,var,per_var,rank_var,ind_var_selected_porcentaje,var.specie) %>% 
  filter(ind_var_selected_porcentaje == 1)

## Grafico % Seleccion Variables  ----------------------------------------------------------------

plot.grafico.var.selected <- function(dataset,norm_n,bbdd_n, method_n) {
  n_var_selected <- 
    dataset %>% 
    filter(ind_var_selected_porcentaje == 1) %>% 
    filter(norm == norm.selected & bbdd == bbdd.selected & method == method_n) %>% 
    group_by(norm,method,bbdd) %>% 
    summarise(max_n = n()) %>% 
    pull(max_n)
  
  plot.final <- 
    dataset %>% 
    filter(norm == norm.selected & bbdd == bbdd.selected & method == method_n) %>% 
    ggplot(aes(x = reorder(var.specie, -per_var),y = per_var*100, group = 1 ))+
    geom_line(lwd = 0.3, colour = colores.brewer[1]) +
    geom_segment(aes(x = 0, xend = n_var_selected , y = percent.selected*100, yend = percent.selected*100), linetype = "dashed" , alpha = 0.6,colour = colores.brewer[4])+
    geom_segment(aes(x = n_var_selected, xend = n_var_selected , y = 0, yend = percent.selected*100), linetype = "dashed" , alpha = 0.6,colour = colores.brewer[4])+
    theme(axis.text.x = element_text(angle = 90),
          axis.text = element_text(size = 7)) +
    xlab("Especies") +
    ylab("Porcentaje") +
    ggtitle(paste0("% Seleccion de la especie en: ",norm.selected,"- ",bbdd.selected,"- ",method_n)) +
    annotate(
      "text",
      x = n_var_selected/2,
      y = (percent.selected+0.06)*100,
      label = paste0("Especies\nSeleccionadas"),
      angle = 0,
      size = 2.5
    )
  
  return(plot.final)  
}


# BBDD Analisis ------------------------------------------------

bbdd.selected <- "Obesidad"
norm.selected <- "S2_N1"

# Grafico % seleccion variable

plot.var.selected.1 <- plot.grafico.var.selected(df.per.selected.variables.v4,
                                                 norm.selected,
                                                 bbdd.selected,
                                                 "Boruta")

# df.per.selected.variables.v4 %>% 
#   filter(norm == norm.selected & bbdd == bbdd.selected & method == "Boruta") %>% 
#   arrange(desc(per_var)) %>% 
#   group_by(norm,method,bbdd) %>% 
#   summarise(max_n = n()) %>% 
#   pull(max_n)

# Grafico para ver mejor preprocesamiento y metodo 

data.bbdd <- 
  f1.estabilidad.datos %>% 
  filter(bbdd == bbdd.selected) %>%
  select(-bbdd) %>% 
  group_by(method,norm) %>% 
  summarise(f1mean = mean(f1mean, na.rm = TRUE)
            ,mean_tan_dist = mean(mean_tan_dis, na.rm = TRUE)) %>% 
  ungroup()


# Medias para scatterplot
mean.scatterplot.f1.bbdd <- mean(data.bbdd$f1mean)
mean.scatterplot.est.bbdd <- mean(data.bbdd$mean_tan_dist)

plot.bbdd.Avg <- 
  ggplot(data.bbdd, aes(x = mean_tan_dist, y = f1mean, shape = norm, color = norm)) +
  geom_point() +
  geom_text_repel(aes(label = method), show.legend = F) +
  labs(
    y = "Media F1-Score",
    x = "Media Estabilidad",
    title = paste(bbdd.selected,"- Media F1-Score vs Media Estabilidad")
  ) +
  scale_shape_manual(values = c(16, 17, 18,19), name = "Preprocesamiento") +
  scale_color_manual(values = colores.brewer[4:7], name = "Preprocesamiento")+
  theme_bw()+
  # Linea Medias
  geom_vline(xintercept = mean.scatterplot.est.bbdd, linetype = "dashed", color = "black",alpha = 0.4)+
  geom_hline(yintercept = mean.scatterplot.f1.bbdd, linetype = "dashed", color = "black",alpha = 0.4)+
  coord_cartesian(xlim = c(0.7, 1), ylim = c(0.78, 0.80))



## Funcion Boxplot Nro Variables -------------------------------------------

unique.methods <- c("BVE","Boruta","PLS_DA","RC","RF","SPLS","SR","T2","VarSelRF","Vip")

fun.boxplot.num.variables <- function(bbdd_n,norm_n,method_n = unique.methods) {
  ncol.bbdd <- all.data.tibble.v2 %>% filter(bbdd == bbdd_n & norm == norm_n) %>% pull(ncol)
  data <- df.num.variables %>% filter(bbdd == bbdd_n & norm == norm_n & method %in% method_n)
  v.plot <- 
    ggplot(data,aes(x = reorder(method, -len_col, FUN = median), 
                    y = len_col)) +
    geom_boxplot(fill = colores.brewer[1]) + 
    labs(x = "Técnica", y = "Nº Variables Seleccionadas")
  # +geom_hline(yintercept = ncol.bbdd, linetype = "dashed", color = colores.brewer[5])
  # +annotate("text", x = 0.5, y = ncol.bbdd + 20, label = paste("Nº Variables sin selección:",ncol.bbdd) , color = colores.brewer[5], size = 4, hjust = 0)
  return(v.plot)
}

# obesidad
boxplot.bbdd.norm.selected <- fun.boxplot.num.variables(bbdd.selected, 
                                                        norm.selected, 
                                                        method_n = c("Boruta","SPLS","Vip"))

#  % conservacion (10*100)/155
#  % conservacion (4*100)/155
#  % conservacion (4*100)/155

## Tabla numero de variables 

table.bbdd.norm <- 
  df.num.variables.v2 %>% 
  filter(norm == norm.selected & bbdd == bbdd.selected & method %in% c("Boruta","SPLS","Vip")) %>% 
  arrange(desc(Mediana))

# Hiperparametros óptimos
# final.results %>% unnest(data) %>% filter(bbdd == "obesidad" & norm == "S1_N1" & method %in% c("Boruta","SPLS","Vip"))

# Variables seleccionadas para los modelos analizados ---------------------------------

# df.num.variables.v2 %>% filter(norm == "S1_N1" & bbdd == "CancerColorectal")

var.boruta.bbdd <- 
  df.per.selected.variables.v5 %>% 
  filter(norm == norm.selected & method == "Boruta" & bbdd == bbdd.selected) %>% 
  pull(var.specie)


var.spls.bbdd <- 
  df.per.selected.variables.v5 %>% 
  filter(norm == norm.selected & method == "SPLS" & bbdd == bbdd.selected) %>% 
  pull(var.specie)

var.vip.bbdd <- 
  df.per.selected.variables.v5 %>% 
  filter(norm == norm.selected & method == "Vip" & bbdd == bbdd.selected) %>% 
  pull(var.specie)


# Modelos SIN CV ----------------------------------------------------------

# Las tecnicas que se analizan son Boruta, SPLS y Vip todas en S1-N1
# Se ajusta los modelos sin realizar validación cruzada.
# Las funciones sin validacion cruzada (1 fold , 1 repeticion)

# Pasamos OTU a SPECIE para esta bbdd

taxa.species2.filter <- taxa.species2 %>% filter(bbdd == bbdd.selected)

#  Base de datos

model.data <- all.data.tibble %>% filter(bbdd == bbdd.selected & norm == norm.selected) %>% pull(data) %>% pluck(1)

# Cambio de nombres de columnas
colnames(model.data) <- sapply(colnames(model.data),function(x){
  new_value = taxa.species2.filter[taxa.species2.filter$OTU.name == x,"var.specie"]
  cond_value = ifelse(length(new_value)==0,x,new_value)
  return(cond_value)
})


model.data.x <- model.data[,-ncol(model.data)]
model.data.y <- model.data$disease


# Hiperparametros óptimos
# final.results %>% unnest(data) %>% filter(bbdd == "obesidad" & norm == norm.selected)

# Boruta
# boruta.bbdd <- BorutaRandomForest.Final(x = model.data.x,
#                                       y = model.data.y,
#                                       cutoff = 0.1,
#                                       folds = 1,
#                                       Rep  = 1,
#                                       seed = 1234,
#                                       varZero = TRUE)
# 
# 
# var.boruta.bbdd <- boruta.bbdd$var %>% pluck(1,1)


# boruta.crc$model[[1]]$importance %>% as.data.frame() %>% arrange(desc(MeanDecreaseGini))
# importance(boruta.crc$model[[1]])
# varImpPlot(boruta.crc$model[[1]])
# VIP

# vip.bbdd <- Vip.Final(x = model.data.x,
#                      y = model.data.y,
#                      ncomp = 1,
#                      vip.value = 2.04,
#                      cutoff = 0,
#                      folds = 1,
#                      Rep = 1,
#                      seed = 1234,
#                      varZero = TRUE)
# 
# var.vip.bbdd <- vip.bbdd$var %>% pluck(1,1)
# 
# VIP(vip.bbdd$model[[1]],opt.comp = vip.bbdd$model[[1]]$ncomp) %>%
#   sort(decreasing = T)


## Venn ---------------------------------------------------------

var.species.list <- list("Boruta" = var.boruta.bbdd,
                         "SPLS" = var.spls.bbdd,
                         "Vip" = var.vip.bbdd)

names(var.species.list) <- c(paste0("Boruta (",length(var.boruta.bbdd),")"),
                             paste0("SPLS (",length(var.spls.bbdd),")"),
                             paste0("Vip (",length(var.vip.bbdd),")"))

compare.especies.selected <- 
  ggvenn(
    var.species.list,
    fill_color = colores.brewer[1:4],
    stroke_size = 0.5, 
    set_name_size = 4,
    text_size = 6,
    show_percentage = FALSE,
  )

# Boruta - SPLS    intersect(var.boruta.bbdd, var.spls.bbdd)
# Boruta - VIP    intersect(var.boruta.bbdd, var.vip.bbdd)
# SPLS -  VIP   intersect(var.spls.bbdd, var.spls.bbdd)
# Estan en boruta pero no en VIP/SPLS setdiff(var.boruta.bbdd, var.spls.bbdd)
# Estan en VIP no Boruta setdiff(var.spls.bbdd, var.boruta.bbdd)

# PCA Graficos Scores ---------------------------------------------------------

# model.data.scale <- scale(model.data.x, center =  TRUE, scale = TRUE)
# xlim_v <- c(-4,4)
# ylim_v <- c(-4,4)
# 
# PCA.CRC.S1.N1 <- PCA(model.data.scale, scale.unit = FALSE , graph = FALSE)
# pca.plot.crc.s1.n1 <- 
#   fviz_pca_ind(PCA.CRC.S1.N1, 
#              label = "none" , 
#              habillage = model.data.y, 
#              palette = colores.brewer[c(1,4)],
#              xlim = xlim_v, 
#              ylim = ylim_v,
#              axes = c(1,2),
#              title = "Cáncer colorrectal S1-N1")+
#   theme(plot.title = element_text(size = 10))
# 
# 
# PCA.CRC.Boruta <- PCA(model.data.scale[,var.boruta.crc],scale.unit = FALSE, graph = FALSE)
# pca.plot.crc.boruta.s1.n1 <- 
#   fviz_pca_ind(PCA.CRC.Boruta, 
#              label = "none",
#              xlim = xlim_v, 
#              ylim = ylim_v,
#              axes = c(1,2),
#              palette = colores.brewer[c(1,4)],
#              habillage = model.data.y,
#              title = "Cáncer colorrectal S1-N1 / Selección Boruta") +
#   theme(plot.title = element_text(size = 10))
# 
# # grafico de comparacion
# biplot.pca.crc.v1 <- ggarrange(pca.plot.crc.s1.n1, 
#                                pca.plot.crc.boruta.s1.n1,
#                                labels = c("A","B"),
#                                common.legend = TRUE,
#                                legend = "bottom")
# 
# 
# 
# PCA.CRC.VarSelRF <- PCA(model.data.scale[,var.varselrf.crc],scale.unit = FALSE, graph = FALSE)
# pca.plot.crc.varselrf.s1.n1 <- 
#   fviz_pca_ind(PCA.CRC.VarSelRF, 
#              label = "none" , 
#              habillage = model.data.y,
#              palette = colores.brewer[c(1,4)],
#              xlim = xlim_v, 
#              ylim = ylim_v,
#              axes = c(1,2),
#              title = "Cáncer colorrectal S1-N1 / Selección VarSelRF")+
#   theme(plot.title = element_text(size = 10))
# 
# 
# PCA.CRC.BVE <- PCA(model.data.scale[,var.bve.crc],scale.unit = FALSE, graph = FALSE)
# pca.plot.crc.bve.s1.n1 <- 
#   fviz_pca_ind(PCA.CRC.BVE, 
#              label = "none" , 
#              habillage = model.data.y,
#              palette = colores.brewer[c(1,4)],
#              xlim = xlim_v, 
#              ylim = ylim_v,
#              axes = c(1,2),
#              title = "Cáncer colorrectal S1-N1 / Selección BVE")+
#   theme(plot.title = element_text(size = 10))
# 
# 
# PCA.CRC.VIP <- PCA(model.data.scale[,var.vip.crc],scale.unit = FALSE, graph = FALSE)
# pca.plot.crc.vip.s1.n1 <- 
#   fviz_pca_ind(PCA.CRC.VIP, 
#              label = "none" , 
#              habillage = model.data.y,
#              palette = colores.brewer[c(1,4)],
#              xlim = xlim_v, 
#              ylim = ylim_v,
#              axes = c(1,2),
#              title = "Cáncer colorrectal S1-N1 / Selección VIP")+
#   theme(plot.title = element_text(size = 10))
# 
# 
# 
# # grafico de comparacion
# pca.crc.todos <- ggarrange(pca.plot.crc.s1.n1, 
#                            pca.plot.crc.boruta.s1.n1,
#                            pca.plot.crc.bve.s1.n1,
#                            pca.plot.crc.vip.s1.n1,
#                            labels = c("A","B"),
#                            common.legend = TRUE,
#                            legend = "bottom")

# PLS-DA Weights  -------------------------------------------------------------------

# Hiperparametros óptimos
# final.results %>% unnest(data) %>% filter(bbdd == "obesidad" & norm == norm.selected)

# Sin seleccion

## PLS-DA Weights - ROPLS -------------------------------------------------------------------

opls.sin.seleccion <- opls(
  x = model.data.x,
  y = model.data.y,
  predI = 2,
  crossvalI = 1,
  scaleC = "standard"
)


opls.vip <- opls(
  x = model.data.x[,var.vip.bbdd],
  y = model.data.y,
  predI = 2,
  crossvalI = 1,
  scaleC = "standard"
)

plot(opls.sin.seleccion , typeVc = "xy-weight")


## PLS-DA Weights - CARET -------------------------------------------------------------------

pls.sin.seleccion <- caret::plsda( x = model.data.x , y = model.data.y , 
                                   ncomp = 2 , validation = "none" , 
                                   center = TRUE , scale = TRUE)


# VIP
plsda.vip <- caret::plsda( x = model.data.x[,var.vip.bbdd] , y = model.data.y , 
                           ncomp = 2 , validation = "none" , 
                           center = TRUE , scale = TRUE)

# weights.pls.caret <- as.data.frame(plsda.vip$projection)  
# yload.test <- as.data.frame(plsda.vip$Yloadings[,1:plsda.vip$ncomp]) / sd(as.numeric(model.data.y))
# bind_rows(weights.pls.caret, as.data.frame(plsda.vip$Yloadings[,1:plsda.vip$ncomp]))
# opls.vip@weightStarMN
# opls.vip@cMN * sd(as.numeric(model.data.y))


## Grafico Weights Custom --------------------------------------------------

fun.plot.weights.xy <-
  function(model,
           name_y = c("Control","Caso"),
           comps = c(1,2),
           title.plot = "Weights XY"
  ) {
    
    # Preparacion Datos
    weight.x.df <- as.data.frame(model$projection)
    weight.x.df$Tipo <- "Especies de Bacterias"
    weight.y.df <- as.data.frame(model$Yloadings[,1:model$ncomp]) / sd(as.numeric(model.data.y))
    weight.y.df$Tipo <- ifelse(rownames(weight.y.df) == "n","Control","Caso")
    # rownames(weight.y.df) <- ifelse(rownames(weight.y.df) == "n",,)
    df.plot <- bind_rows(weight.x.df, weight.y.df)
    
    # Plot
    
    final.plot <-
      df.plot %>%
      ggplot(aes(x = df.plot[,comps[1]] ,
                 y = df.plot[,comps[2]],
                 colour = Tipo,
                 label = rownames(df.plot)))+
      geom_point(size = 0.5)+
      geom_text(hjust = 1, vjust = 1, size = 3, show.legend = F)+
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
      scale_color_manual(values = colores.brewer[c(4,5,6)]) +
      theme_minimal() +
      xlab(paste0("w*c[",comps[1],"]"))+
      ylab(paste0("w*c[",comps[2],"]"))+
      xlim(-1.3,1.3)+
      ylim(-1.3,1.3)+
      theme(legend.position="top")
    
    return(final.plot)
  }

# Graficos de weights

plot.weights.1<- fun.plot.weights.xy(pls.sin.seleccion, title.plot = paste(bbdd.selected,"Sin selección"))
plot.weights.2 <- fun.plot.weights.xy(plsda.vip, title.plot = paste(bbdd.selected,"VIP"))

# coef
plsda.vip$coefficients[,2,2]

# Grafico arrange

final.all.weights <-
  ggarrange(plot.weights.1,
            plot.weights.2,
            labels = c("A","B"),
            common.legend = TRUE,
            legend = "top")

# PLS-DA Scores -----------------------------------------------------------
## Grafico Scores PLS-DA Custom --------------------------------------------


fun.plot.scores.pls.da <-
  function(model,
           comps = c(1,2),
           title.plot = "Scores PLS-DA"
  ) {
    
    # Preparacion Datos
    scores.df <- as.data.frame(model@scoreMN)
    scores.df$disease <- model@suppLs$y
    
    # Max para ejes y que no se pierdan outliers
    max.abs.1comp <-  max(abs(scores.df[,comps[1]]))
    max.abs.2comp <-  max(abs(scores.df[,comps[2]]))
    max.axis <- max(max.abs.1comp,max.abs.2comp)
    
    # % Var
    
    var.1comp <- round(model@modelDF$R2X[comps[1]]*100)
    var.2comp <- round(model@modelDF$R2X[comps[2]]*100)
    
    # Plot
    
    final.plot <-
      scores.df %>%
      ggplot(aes(x = scores.df[,comps[1]],
                 y = scores.df[,comps[2]],
                 colour = disease ,
                 shape = disease,
                 label = rownames(scores.df)))+
      geom_point()+
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
      scale_color_manual(values = colores.brewer[c(5,4)]) +
      theme_minimal() +
      xlab(paste0("t[",comps[1],"]"," (",var.1comp,"%)"))+
      ylab(paste0("t[",comps[2],"]"," (",var.2comp,"%)"))+
      xlim(max.axis * -1,max.axis)+
      ylim(max.axis * -1,max.axis)+
      ggtitle(title.plot)+
      theme(legend.title = element_blank(),
            plot.title = element_text(size = "8"))
    
    return(final.plot)
  }

## Sin selección
plot.scores.ss <- fun.plot.scores.pls.da(opls.sin.seleccion,
                                         title.plot = "Scores PLS-DA / Sin Selección")
## VIP

plot.scores.vip <- fun.plot.scores.pls.da(opls.vip,
                                          title.plot = "Scores PLS-DA / Vip - SPLS")

fig.scores.compare <- 
  ggarrange(
    plot.scores.ss,
    plot.scores.vip,
    labels = c("A","B"),
    common.legend = TRUE,
    legend = "top")

# plot(opls.sin.seleccion , typeVc = "x-score")
# plot(opls.bve , typeVc = "x-score")
# plot(opls.bve , typeVc = "x-score", parCompVi = c(1, 3))
