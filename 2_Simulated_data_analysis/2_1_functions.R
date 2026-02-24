
# ============================================================================ #
# ========       LONGITUDINAL CLUSTERING METHODS COMPARISON         ========== #
# ========   Functions created for the simulated-data analysis      ========== #
# ============================================================================ #


# Functions to simulate FEP patients data base on PEPs cohort ----
simular_grupo <- function(beta0, slopes,sd_b, sd_epsilon, ymin, ymax, time, class_name, n, format = "long") {
  # Simular efectos aleatorios de intercepto para cada individuo
  bi <- rnorm(n, mean = 0, sd = sd_b)  # Variabilidad interindividual
  # Crear un data frame con las trayectorias
  datos <- do.call(rbind, lapply(1:n, function(id) {
    epsilon <- rnorm(length(time), mean = 0, sd = sd_epsilon)  # Error aleatorio
    y <- beta0 + slopes + bi[id] + epsilon
    y <- round(pmin(pmax(y, ymin), ymax))  # Limites PANSS y redondeo a enteros
    data.frame(id = paste0(class_name, "_", id),
               group = class_name,
               time = time,
               PANSS = y)
  }))
  if (format == "wide") {
    datos <- datos %>%
      pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_t")
  }
  return(datos)
}

add_realistic_missingness <- function(sim_data,  format = c("long", "wide"), seed = NULL){
  
  format <- match.arg(format)
  
  if(!is.null(seed)) set.seed(seed)
  
  # IDs únics
  ids <- unique(sim_data$id)
  N <- length(ids)
  
  # Quantitats segons data real (cohort PEPS)
  n_2_only   <- 2
  n_6_only   <- 3
  n_12_only  <- 26
  n_2_6      <- 1
  n_2_12     <- 1
  n_6_12     <- 23
  n_2_6_12   <- 0
  
  total_missing_ids <- n_2_only + n_6_only + n_12_only +
    n_2_6 + n_2_12 + n_6_12 + n_2_6_12
  
  if(total_missing_ids > N){
    stop("No hi ha suficients IDs a sim_data per assignar aquest patró de missings.")
  }
  
  # Barregem IDs i assignem patrons
  ids_shuffled <- sample(ids)
  
  idx <- 1
  
  ids_2_only  <- ids_shuffled[idx:(idx + n_2_only - 1)]
  idx <- idx + n_2_only
  
  ids_6_only  <- ids_shuffled[idx:(idx + n_6_only - 1)]
  idx <- idx + n_6_only
  
  ids_12_only <- ids_shuffled[idx:(idx + n_12_only - 1)]
  idx <- idx + n_12_only
  
  ids_2_6     <- ids_shuffled[idx:(idx + n_2_6 - 1)]
  idx <- idx + n_2_6
  
  ids_2_12    <- ids_shuffled[idx:(idx + n_2_12 - 1)]
  idx <- idx + n_2_12
  
  ids_6_12    <- ids_shuffled[idx:(idx + n_6_12 - 1)]
  idx <- idx + n_6_12
  
  ids_2_6_12  <- if(n_2_6_12 > 0) ids_shuffled[idx:(idx + n_2_6_12 - 1)] else character(0)
  
  # --- APLICAR MISSINGS ---
  if(format == "long"){
    
    # comprovació mínima
    needed_cols <- c("id", "time", "PANSS")
    if(!all(needed_cols %in% names(sim_data))){
      stop("Format long: calen les columnes id, time i PANSS.")
    }
    
    sim_data$PANSS[sim_data$id %in% ids_2_only  & sim_data$time == 2]  <- NA
    sim_data$PANSS[sim_data$id %in% ids_6_only  & sim_data$time == 6]  <- NA
    sim_data$PANSS[sim_data$id %in% ids_12_only & sim_data$time == 12] <- NA
    
    sim_data$PANSS[sim_data$id %in% ids_2_6  & sim_data$time %in% c(2,6)]   <- NA
    sim_data$PANSS[sim_data$id %in% ids_2_12 & sim_data$time %in% c(2,12)]  <- NA
    sim_data$PANSS[sim_data$id %in% ids_6_12 & sim_data$time %in% c(6,12)]  <- NA
    sim_data$PANSS[sim_data$id %in% ids_2_6_12 & sim_data$time %in% c(2,6,12)] <- NA
    
  } else if(format == "wide"){
    
    needed_cols <- c("id", "PANSS_t2", "PANSS_t6", "PANSS_t12")
    if(!all(needed_cols %in% names(sim_data))){
      stop("Format wide: calen les columnes id, PANSS_t2, PANSS_t6 i PANSS_t12.")
    }
    
    sim_data$PANSS_t2[sim_data$id %in% ids_2_only]  <- NA
    sim_data$PANSS_t6[sim_data$id %in% ids_6_only]  <- NA
    sim_data$PANSS_t12[sim_data$id %in% ids_12_only] <- NA
    
    sim_data$PANSS_t2[sim_data$id %in% ids_2_6] <- NA
    sim_data$PANSS_t6[sim_data$id %in% ids_2_6] <- NA
    
    sim_data$PANSS_t2[sim_data$id %in% ids_2_12] <- NA
    sim_data$PANSS_t12[sim_data$id %in% ids_2_12] <- NA
    
    sim_data$PANSS_t6[sim_data$id %in% ids_6_12] <- NA
    sim_data$PANSS_t12[sim_data$id %in% ids_6_12] <- NA
    
    if(length(ids_2_6_12) > 0){
      sim_data$PANSS_t2[sim_data$id %in% ids_2_6_12] <- NA
      sim_data$PANSS_t6[sim_data$id %in% ids_2_6_12] <- NA
      sim_data$PANSS_t12[sim_data$id %in% ids_2_6_12] <- NA
    }
  }
  
  return(sim_data)
}


# Función para calcular kappa de cohen ----
calcular_kappa_clusters <- function(df, group_col = "group_num", cluster_cols, k=k) {
  # Estandaritzem o alineem les etiquetes en funció del major num de coincidències
  ref_group <- df[[group_col]]
  
  for (j in cluster_cols) {
    perms <- permutations(n = k, r = k)
    best_match <- 0
    best_perm <- seq_len(k)
    
    for (i in 1:nrow(perms)) {
      aux <- data.frame("real" = seq_len(k), "recoded" = perms[i,])
      recoded <- vector()
      for(obs in 1:length(df[[j]])){
        if(k==3){
          recoded[obs] <- ifelse(df[obs, j] == aux$real[1], aux$recoded[1],
                                 ifelse(df[obs, j] == aux$real[2], aux$recoded[2], aux$recoded[3]))
        }else{
          recoded[obs] <- ifelse(df[obs, j] == aux$real[1], aux$recoded[1],
                                 ifelse(df[obs, j] == aux$real[2], aux$recoded[2], 
                                        ifelse(df[obs, j] == aux$real[3], aux$recoded[3], aux$recoded[4])))
        }
      }
      
      if(length(unique(df[[j]])) == k){
        match_score <- sum(diag(table(df[[group_col]], recoded)))
      } else{
        x <- table(df[[group_col]], recoded)
        match_score <- 0
        for(r in rownames(x)){
          for (c in colnames(x)){
            if (r == c) {
              match_score <- match_score + x[r, c]
            }
          }
        }
      }
      if (match_score > best_match) {
        best_match <- match_score
        best_perm <- aux[,2]
      }
    }
    for(obs in 1:length(df[[j]])){
      if(k==3){
        df[obs, j] <- ifelse(df[obs, j] == aux$real[1], best_perm[1],
                             ifelse(df[obs, j] == aux$real[2], best_perm[2], best_perm[3]))
      }else{
        df[obs, j] <- ifelse(df[obs, j] == aux$real[1], best_perm[1],
                             ifelse(df[obs, j] == aux$real[2], best_perm[2], 
                                    ifelse(df[obs, j] == aux$real[3], best_perm[3], best_perm[4])))
      }
    }
    #print(paste0("Recode ", j, " with best permutation: ", paste(best_perm, collapse = ", "), "\n"))
    #print(table(df[[group_col]], df[[j]]))
  }
  # Inicialitzem llistes buides per guardar resultats
  resultats <- lapply(cluster_cols, function(col) {
    kappa <- cohen.kappa(cbind(df[[group_col]], df[[col]]))
    
    data.frame("mètode" = col,
               "kappa" = kappa$weighted.kappa,
               "se" = kappa$var.weighted,
               "ic95_inf" = kappa$confid[2, 1],
               "ic95_sup" = kappa$confid[2, 3])
  })
  # Unim els resultats en un sol tibble
  return(bind_rows(resultats))
}



# Función para calcular ICC ----
calculate_icc <- function(data, format = "long") {
  # PAS 1: Fer trajectòries individuals (una fila = un subjecte)
  if(format == "long"){
    data <- data %>%
      select(id, group, time, PANSS) %>%
      pivot_wider(names_from = time, values_from = PANSS, names_prefix = "PANSS_") %>%
      distinct() 
  }
  
  # PAS 2: Matriu de distàncies entre trajectòries individuals (n = 237 x 237)
  X <- as.matrix(dist(data %>% select(starts_with("PANSS_"))))
  
  # PAS 3: Dataframe nt = nombre de subjectes per grup
  nt <- data %>%  count(group, name = "n")
  
  # PAS 4: Aplicar ICC
  icc_result <- iccTraj::ICC(X = X, nt = nt)
  return(icc_result)
}



# Get results ----
# Función para obtener el valor de k con el valor máximo de un subíndice dado para cada B
opt_k_for_metric <- function(data, metric, min.percentage) {
  
  if (!is.null(min.percentage)){
    data <- subset(data, min.class >= min.percentage) 
  }
  
  if (nrow(data) == 0) {
    # Retorna una fila indicant que no hi ha cap solució vàlida
    return(data.frame(
      B = NA,
      k = NA,
      Value = NA,
      Metric = metric
    ))
  }
  
  if(metric == "AIC" || metric=="BIC" || metric =="SABIC"){
    data %>%
      group_by(B) %>%
      filter(!!sym(metric) == min(!!sym(metric))) %>%
      dplyr::select(B, k, !!sym(metric)) %>%
      rename(Value = !!sym(metric)) %>%
      mutate(Metric = metric)
  } else{
    data %>%
      group_by(B) %>%
      filter(!!sym(metric) == max(!!sym(metric))) %>%
      dplyr::select(B, k, !!sym(metric)) %>%
      rename(Value = !!sym(metric)) %>%
      mutate(Metric = metric)
  }
}


# Función para obtener % de veces que k es máximo segun cada indice
analyse_indexes <- function(data, k_range= 2:8, LCMM = NULL, min.percentage = NULL) {
  
  # Aplicar la función a cada subíndice y combinarlos en un solo data frame
  if(LCMM){
    resultados <- bind_rows(opt_k_for_metric(data, "loglik", min.percentage), 
                            opt_k_for_metric(data, "AIC", min.percentage),
                            opt_k_for_metric(data, "BIC", min.percentage), 
                            opt_k_for_metric(data, "SABIC", min.percentage), 
                            opt_k_for_metric(data, "entropy", min.percentage)) 
  } else{
    resultados <- bind_rows(opt_k_for_metric(data, "Silhouette", min.percentage), 
                            opt_k_for_metric(data, "Dunn", min.percentage),
                            opt_k_for_metric(data, "Gamma", min.percentage), 
                            opt_k_for_metric(data, "CH", min.percentage))  
  }
  
  
  # Contar las veces que cada valor de k ha sido elegido como el óptimo para cada índice
  conteos <- resultados %>%
    group_by(Metric, k) %>%
    dplyr::summarise(Count = n(), .groups = 'drop')
  
  # Crear una estructura completa de métricas y k
  full_grid <- expand.grid(Metric = unique(resultados$Metric), k = k_range)
  
  # Unir conteos con la estructura completa
  conteos_completos <- full_grid %>%
    left_join(conteos, by = c("Metric", "k")) %>%
    mutate(Count = ifelse(is.na(Count), 0, Count))
  
  # Calcular el porcentaje de veces que cada k ha sido elegido como el óptimo para cada índice
  total_counts <- conteos %>%
    group_by(Metric) %>%
    dplyr::summarise(Total = sum(Count), .groups = 'drop')
  
  percentages <- conteos_completos %>%
    left_join(total_counts, by = "Metric") %>%
    mutate(Percentage = (Count / Total) * 100) %>%
    dplyr::select(Metric, k, Percentage)
  
  # Reformatear los resultados en el formato deseado
  resultados_finales <- percentages %>%
    pivot_wider(names_from = k, values_from = Percentage, values_fill = list(Percentage = 0))
  
  return(resultados_finales)
}

# Función para obtener tabla con N(%) de veces que se ha elimiado k por tener < min.perectange
discards <- function(data, k_range, min.percentage){
  filtered_data <- data[data$min.class < min.percentage, ]
  
  discard_summary <- data.frame(matrix(nrow = 1, ncol = length(k_range)))
  colnames(discard_summary) <- k_range
  rownames(discard_summary) <- "N(%)"
  
  for (k_value in k_range) {
    k_data <- filtered_data[filtered_data$k == k_value, ]
    N <- nrow(k_data) # num total de descartes para k
    percent <- round((N / max(data$B)) * 100, 2)  # porcentaje
    discard_summary[[paste(k_value)]] <- paste0(N, "(", percent, "%)")
  }
  return(discard_summary)
}

# Funcion para obtener cuantas soluciones descartadas huberan sido óptimas
discarded_optimal_summary <- function(data, k_range, min.percentage, LCMM = FALSE) {
  # Filtrar descartades
  total_discards <- data %>% filter(min.class < min.percentage)
  
  # Si no hi ha descartades, retornem guions a tot
  if (nrow(total_discards) == 0) {
    return(tibble(k = k_range) %>%
             mutate(val = "-") %>%
             pivot_wider(names_from = k, values_from = val))
  }
  
  # Models òptims segons alguna mètrica
  if (LCMM) {
    optimal_models <- bind_rows(
      opt_k_for_metric(data, "loglik", NULL), 
      opt_k_for_metric(data, "AIC", NULL),
      opt_k_for_metric(data, "BIC", NULL), 
      opt_k_for_metric(data, "SABIC", NULL), 
      opt_k_for_metric(data, "entropy", NULL)
    )
  } else {
    optimal_models <- bind_rows(
      opt_k_for_metric(data, "Silhouette", NULL), 
      opt_k_for_metric(data, "Dunn", NULL),
      opt_k_for_metric(data, "Gamma", NULL), 
      opt_k_for_metric(data, "CH", NULL)
    )
  }
  
  # Total de descartades per k
  total_discards_k <- total_discards %>%
    count(k, name = "total_discarded")
  
  # Solucions descartades que haurien estat òptimes
  joined <- semi_join(total_discards, optimal_models, by = c("B", "k")) %>%
    count(k, name = "N")
  
  # Fusionar les dues fonts d'informació i afegir 0 (0%) on pertoqui
  summary <- total_discards_k %>%
    left_join(joined, by = "k") %>%
    mutate(N = ifelse(is.na(N), 0, N),
           pct = round(100 * N / total_discarded, 1),
           val = paste0(N, " (", pct, "%)"))
  
  # Taula completa amb tots els valors de k (incloent guions si no hi ha cap descartada per a aquell k)
  output <- tibble(k = k_range) %>%
    left_join(summary %>% select(k, val), by = "k") %>%
    mutate(val = ifelse(is.na(val), "-", val)) %>%
    pivot_wider(names_from = k, values_from = val)
  
  return(output)
}



# Función para obtener una tabla con la media(sd) de la diferència de los dos valores optimos de cada metrica por cada k
metric_optimal_diff <- function(data, k_range = k_range) {
  # Inicialització
  metrics = colnames(data)[3:(ncol(data)-1)]
  summary_table <- data.frame(matrix(ncol = length(metrics), nrow = 1))
  rownames(summary_table) <- "MedianDifference(IQR)"
  colnames(summary_table) <- metrics
  
  # Càlcul de la diferència entre els dos valors òptims per cada remostra
  for (metric in metrics) {
    diffs <- numeric()
    
    for (b in unique(data$B)) {
      b_data <- data[data$B == b, ]
      
      # Ordenar els valors de k segons la mètrica en ordre ascendent (si AIC/BIC/SABIC) o descendent (otherwise)
      if(metric == "AIC" || metric=="BIC" || metric =="SABIC"){
        b_data <- b_data[order(b_data[[metric]]), ]
      } else{
        b_data <- b_data[order(-b_data[[metric]]), ]
      }
      
      # Seleccionar els dos primers valors (òptims)
      diff_val <- abs(b_data[[metric]][1] - b_data[[metric]][2])
      diffs <- c(diffs, diff_val)
    }
    
    # Calcul de la mitjana i desviació estàndard
    median_diff <- round(median(diffs, na.rm = TRUE), 2)
    iqr_diff <- round(IQR(diffs, na.rm = TRUE), 2)
    
    # Guardar el resultat a la taula
    summary_table[1,which(metric == metrics)] <- paste0(median_diff, "(", iqr_diff, ")")
  }
  
  return(summary_table)
}

# Función para obtener una tabla con la media(sd) de cada métrica por cada k
metric_summary <- function(data, k_range = k_range) {
  # Inicialización
  metrics = colnames(data)[3:(ncol(data)-1)]
  summary_table <- data.frame(matrix(nrow = length(metrics), ncol = length(k_range)))
  colnames(summary_table) <- k_range
  summary_table$Metric <- metrics
  summary_table <- summary_table[,c(ncol(summary_table),1:(ncol(summary_table)-1))]
  
  # Calcular media y desviación estándar, redondeados a dos decimales
  for (metric in metrics) {
    for (k_value in k_range) {
      k_data <- data[data$k == k_value, metric, drop = FALSE]
      
      mean_val <- round(mean(k_data[[metric]], na.rm = TRUE), 2)
      sd_val <- round(sd(k_data[[metric]], na.rm = TRUE), 2)
      summary_table[which(metric == metrics), paste0(k_value)] <- paste0(mean_val, "(", sd_val, ")")
    }
  }
  
  # Devolver resultados
  return(summary_table)
}

# Función para obtener una matrix de plots de densidad con la distribución d elas métricas
metric_density_plots <- function(data1, data2, data3, k_range = k_range, LCMM=NULL) {
  
  # Agrupar datos
  data = rbind(data1, data2, data3)
  if(LCMM){
    data$Technique <- rep(c("Linear", "Quadratic", "Splines"), each= nrow(data1))
    n=5
  }else{
    data$Technique <- rep(c("Hierarchical clustering", "K-means", "PAM"), each= nrow(data1))
    n=4;
  }
  
  # Transformar datos a formato largo para facet_grid()
  long_data <- as.data.frame(pivot_longer(data, 
                                          cols = 3:(ncol(data)-2), 
                                          names_to = "metric", 
                                          values_to = "value"))
  # Hacer plot
  p <- ggplot(long_data, aes(x = value)) + 
    geom_density(alpha = 0.9, position = "identity", aes(fill = Technique)) + 
    theme_minimal() +
    theme(legend.position = "bottom")+
    facet_wrap(k ~ metric, scales = "free", ncol=n, nrow=length(k_range)) +
    scale_fill_brewer(palette = "Pastel1")
  p
  return(p)
}


# Distance-based methods ----
# Función para inicializar los data frames para los resultados
initialize_clust_results <- function(num_resamples = 100, k_range = 2:6) {
  data.frame("B" = rep(1:num_resamples, each = length(k_range)),
             "k" = rep(k_range, times = num_resamples),
             "Silhouette" = rep(NA, length(k_range)*num_resamples),
             "Dunn" = rep(NA, length(k_range)*num_resamples),
             "Gamma" = rep(NA, length(k_range)*num_resamples),
             "CH" = rep(NA, length(k_range)*num_resamples),
             "min.class" = rep(NA, length(k_range)*num_resamples))
}


# Función para actualizar los resultados
update_clust_results <- function(results, B, k, stats) {
  index <- which(results$B == B & results$k == k)
  results[index, 3] <- stats$avg.silwidth
  results[index, 4] <- stats$dunn2
  results[index, 5] <- stats$g2
  results[index, 6] <- stats$ch
  results[index, 7] <- (stats$min.cluster.size/ sum(stats$cluster.size))*100
  return(results)
}

# Función principal que encapsula todo el proceso
perform_clustering_simulation <- function(n_ind, prop_groups, beta0_values, slopes_values,
                                          sd_b, sd_epsilon, ymin, ymax, format= "wide", time,
                                          num_resamples = 100, k_range = 2:6, trace=T, min.percentage = NULL) {
  # Inicializar los data frames para los resultados
  clust.hier.results <- initialize_clust_results(num_resamples = num_resamples, k_range = k_range)
  clust.kmeans.results <- initialize_clust_results(num_resamples = num_resamples, k_range = k_range)
  clust.pam.results <- initialize_clust_results(num_resamples = num_resamples, k_range = k_range)
  icc <- vector()
  
  for (B in 1:num_resamples) {
    if (trace) { cat("Bootstrap sample:", B, "\n") }
    # 1. Simulate data
    n_groups <- round(n_ind * prop_groups)
    data.sim <- bind_rows(
      mapply(simular_grupo, beta0 = beta0_values, slopes = slopes_values,
             MoreArgs = list(sd_b = sd_b, sd_epsilon = sd_epsilon, ymin = ymin, ymax = ymax,
                             time = time, format = format),
             class_name = names(prop_groups), n = n_groups, SIMPLIFY = FALSE)
    )
    data.sim$group <- factor(data.sim$group, levels = names(prop_groups))
    data.sim <- transform(data.sim, id = as.numeric(factor(id)))
    data.sim <- add_realistic_missingness(data.sim, format=format)
    
    # 2. Compute distance matrix
    D <- dist(as.matrix(data.sim[, 3:6]))
    
    # 3. Perform hierarchical clustering
    clust.hier <- hclust(D, method = "ward.D2")
    
    for (k in k_range) {
      if (trace) { cat("k =", k, "\n") }
      # 3.1. Hierarchical clustering with Ward's method
      x <- cutree(clust.hier, k)
      stats <- cluster.stats(D, x, G2 = TRUE)
      clust.hier.results <- update_clust_results(clust.hier.results, B, k, stats)
      
      # 3.2. K-means clustering
      x <- kmeans(D, centers = k)$cluster
      stats <- cluster.stats(D, x, G2 = TRUE)
      clust.kmeans.results <- update_clust_results(clust.kmeans.results, B, k, stats)
      
      # 3.3. PAM clustering
      x <- pam(D, k = k, diss = TRUE)$clustering
      stats <- cluster.stats(D, x, G2 = TRUE)
      clust.pam.results <- update_clust_results(clust.pam.results, B, k, stats)
    }
    if (trace) { cat("\n") }
  }
  
  # 4. Compute results
  ## 4.1. Get optimal k for each resample and index
  clust.hier.final <- analyse_indexes(clust.hier.results, k_range, LCMM = FALSE,
                                      min.percentage = min.percentage)
  clust.kmeans.final <- analyse_indexes(clust.kmeans.results, k_range, LCMM = FALSE,
                                        min.percentage = min.percentage)
  clust.pam.final <- analyse_indexes(clust.pam.results, k_range, LCMM = FALSE,
                                     min.percentage = min.percentage)
  
  results.clust <- as.data.frame(bind_rows(clust.hier.final, clust.kmeans.final, clust.pam.final))
  results.clust[,-1] <- round(results.clust[,-1],2)
  results.clust$Technique <- rep(c("Hierarchical clustering", "K-means", "PAM"),each = 4)
  results.clust <- results.clust[, c(ncol(results.clust),1:ncol(results.clust)-1)]
  
  ## 4.2. Get number of discarded results because min.class < min.percentage
  if (!is.null(min.percentage)) {
    discarded.hier <- discards(clust.hier.results, k_range, min.percentage)
    discarded.kmeans <- discards(clust.kmeans.results, k_range, min.percentage)
    discarded.pam <- discards(clust.pam.results, k_range, min.percentage)
  }
  results.discards <- bind_rows(discarded.hier, discarded.kmeans, discarded.pam)
  rownames(results.discards) <- c("Hierarchical clustering", "K-means", "PAM")
  
  ## 4.3. Get how many of the solutions discarded by min.class < min.percentage would have been considered optimal 
  results.discarded_which_optimal <- bind_rows(
    discarded_optimal_summary(clust.hier.results, k_range, min.percentage, LCMM = FALSE),
    discarded_optimal_summary(clust.kmeans.results, k_range, min.percentage, LCMM = FALSE),
    discarded_optimal_summary(clust.pam.results, k_range, min.percentage, LCMM = FALSE)
  )
  rownames(results.discarded_which_optimal) <- c("Hierarchical clustering", "K-means", "PAM")
  results.discarded_which_optimal <- as.data.frame(results.discarded_which_optimal)
  
  return(list(Optimal_Models = results.clust, 
              Discard_Summary = results.discards,
              Discarded_Optimal = results.discarded_which_optimal))
}



# Latent class mixed models ----
# Función para inicializar los data frames para los resultados
initialize_lcmm_results <- function(num_resamples = 100, k_range = 2:6) {
  data.frame("B" = rep(1:num_resamples, each = length(k_range)),
             "k" = rep(k_range, times = num_resamples),
             "loglik" = rep(NA, length(k_range)*num_resamples),
             "AIC" = rep(NA, length(k_range)*num_resamples),
             "BIC" = rep(NA, length(k_range)*num_resamples),
             "SABIC" = rep(NA, length(k_range)*num_resamples),
             "entropy" = rep(NA, length(k_range)*num_resamples),
             "min.class" = rep(NA, length(k_range)*num_resamples))
}

# Función para actualizar los resultados
update_lcmm_results <- function(results, B, k, tab) {
  index <- which(results$B == B & results$k == k)
  
  results[index, 3] <- tab$loglik
  results[index, 4] <- tab$AIC
  results[index, 5] <- tab$BIC
  results[index, 6] <- tab$SABIC
  results[index, 7] <- tab$entropy
  results[index, 8] <- min(tab[,7:(ncol(tab))])
  
  return(results)
}

# Función principal que contiene todo el proceso
perform_lcmm_simulation <- function(n_ind, prop_groups, beta0_values, slopes_values,
                                    sd_b, sd_epsilon, ymin, ymax, format= "long", time, num_resamples = 100,
                                    k_range = 2:5, trace = TRUE, min.percentage = NULL) {
  
  # Inicializar los data frames para los resultados
  lcmm.linear.results <- initialize_lcmm_results(num_resamples, k_range)
  lcmm.quadratic.results <- initialize_lcmm_results(num_resamples, k_range)
  lcmm.splines.results <- initialize_lcmm_results(num_resamples, k_range)
  
  sum.linear <- sum.quadratic <- sum.splines <- 0 # Conteo de no convergencias
  B <- 1
  
  while (B <= num_resamples) {
    #if (trace) cat("Simulation sample:", B, "\n")
    if (trace) message(paste0("Simulation sample:", B))
    all_converged <- FALSE # Variable de control per assegurar la convergència de tots els models.
    
    while (!all_converged) {
      # 1. Simulate data and compute icc
      n_groups <- round(n_ind * prop_groups)
      data.sim <- bind_rows(
        mapply(simular_grupo, beta0 = beta0_values, slopes = slopes_values,
               MoreArgs = list(sd_b = sd_b, sd_epsilon = sd_epsilon, ymin = ymin, ymax = ymax,
                               time = time, format = format),
               class_name = names(prop_groups), n = n_groups, SIMPLIFY = FALSE)
      )
      data.sim$group <- factor(data.sim$group, levels = names(prop_groups))
      data.sim <- transform(data.sim, id = as.numeric(factor(id)))
      data.sim <- add_realistic_missingness(data.sim, format=format)
      
      # 2. Ajustar modelos LCMM
      assign("lcmm.linear", hlme(PANSS ~ time, random = ~1, subject = "id", data = data.sim, ng = 1), envir = .GlobalEnv)
      assign("lcmm.quad", hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = data.sim, ng = 1), envir = .GlobalEnv)
      assign("lcmm.splines", hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = data.sim, ng = 1), envir = .GlobalEnv)
      
      convergences <- c(TRUE, TRUE, TRUE)  # Control per linear, quadratic i splines
      
      for (k in k_range) {
        #if (trace) cat("k =", k, "\n")
        if (trace) message(paste0("k =", k))
        #reps <- ifelse(k==2, 5, ifelse(k==3, 10, ifelse(k==4, 15, 20)))
        reps = 100
        # 3.1. Linear
        x <- gridsearch(hlme(PANSS ~ time, random = ~1, subject = "id", data = data.sim,
                             mixture = ~ time, ng = k, B = random(lcmm.linear)), rep = reps, maxiter = 10, minit = lcmm.linear)
        if (x$conv %in% c(1, 3)) {
          lcmm.linear.results <- update_lcmm_results(lcmm.linear.results, B, k, 
                                                     as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
        } else {
          sum.linear <- sum.linear + 1
          convergences[1] <- FALSE
          if (trace) cat("Linear model did not converge. Re-sampling...\n")
          break
        }
        
        # 3.2. Quadratic
        x <- gridsearch(hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = data.sim,
                             mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad)), rep = reps, maxiter = 10, minit = lcmm.quad)
        if (x$conv %in% c(1, 3)) {
          lcmm.quadratic.results <- update_lcmm_results(lcmm.quadratic.results, B, k, 
                                                        as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
        } else {
          sum.quadratic <- sum.quadratic + 1
          convergences[2] <- FALSE
          if (trace) cat("Quadratic model did not converge. Re-sampling...\n")
          break
        }
        
        # 3.3. Splines
        x <- gridsearch(hlme(PANSS ~ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = data.sim,
                             mixture = ~ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines)), rep = reps, maxiter = 10, minit = lcmm.splines)
        if (x$conv %in% c(1, 3)) {
          lcmm.splines.results <- update_lcmm_results(lcmm.splines.results, B, k, 
                                                      as.data.frame(summarytable(x, which = c("G", "loglik", "AIC", "BIC", "SABIC", "entropy", "%class"), display = FALSE)))
        } else {
          sum.splines <- sum.splines + 1
          convergences[3] <- FALSE
          if (trace) cat("Splines model did not converge. Re-sampling...\n")
          break
        }
      }
      
      # Si tots els models han convergit, sortim del bucle intern
      all_converged <- all(convergences)
    }
    
    B <- B + 1 # Passar a la següent remostra només si tot ha convergit
    if (trace) cat("\n")
  }
  # Cleanup: remove temporary objects from global environment
  rm(lcmm.linear, lcmm.quad, lcmm.splines, envir = .GlobalEnv)
  
  # PROVES:
  message("Resultats lineals:")
  print(lcmm.linear.results)
  message("Resultats quadratics:")
  print(lcmm.quadratic.results)
  message("Resultats splines:")
  print(lcmm.splines.results)
  
  # 4. Compute results
  ## 4.1. Get optimal k for each resample and index
  lcmm.linear.final <- analyse_indexes(lcmm.linear.results, k_range= k_range,
                                       LCMM = T, min.percentage =min.percentage)
  lcmm.quadratic.final <- analyse_indexes(lcmm.quadratic.results, k_range= k_range,
                                          LCMM = T, min.percentage =min.percentage)
  lcmm.splines.final <- analyse_indexes(lcmm.splines.results, k_range= k_range,
                                        LCMM = T, min.percentage =min.percentage)
  # PROVES:
  message("Resultats analyse_indexes() lineals:")
  print(lcmm.linear.final)
  message("Resultats analyse_indexes() quadratics:")
  print(lcmm.quadratic.final)
  message("Resultats analyse_indexes() splines:")
  print(lcmm.splines.final)
  
  results.lcmm <- as.data.frame(bind_rows(lcmm.linear.final, lcmm.quadratic.final, lcmm.splines.final))
  results.lcmm[,-1] <- round(results.lcmm[,-1],2)
  results.lcmm$Technique <- rep(c("Linear", "Quadratic", "Splines"),each = 5)
  results.lcmm <- results.lcmm[, c(ncol(results.lcmm),1:ncol(results.lcmm)-1)]
  
  
  ## 4.2. Get number of discarded results because min.class < min.percentage
  if (!is.null(min.percentage)) {
    discarded.linear <- discards(lcmm.linear.results, k_range, min.percentage)
    discarded.quadratic <- discards(lcmm.quadratic.results, k_range, min.percentage)
    discarded.splines <- discards(lcmm.splines.results, k_range, min.percentage)
  }
  results.discards <- bind_rows(discarded.linear, discarded.quadratic, discarded.splines)
  rownames(results.discards) <- c("Linear", "Quadratic", "Splines")
  
  ## 4.3. Get convergence results
  convergence.linear <- data.frame("N"= sum.linear, "Percentage_sobre_B*length_k"= round(sum.linear/nrow(lcmm.linear.results)*100,2))
  convergence.quadratic <- data.frame("N"= sum.quadratic, "Percentage_sobre_B*length_k"= round(sum.quadratic/nrow(lcmm.quadratic.results)*100,2))
  convergence.splines <- data.frame("N"= sum.splines, "Percentage_sobre_B*length_k"= round(sum.splines/nrow(lcmm.splines.results)*100,2))
  
  results.convergence <- bind_rows(convergence.linear, convergence.quadratic, convergence.splines)
  rownames(results.convergence) <- c("Linear", "Quadratic", "Splines")
  
  ## 4.4. Get how many of the solutions discarded by min.class < min.percentage would have been considered optimal 
  results.discarded_which_optimal <- bind_rows(
    discarded_optimal_summary(lcmm.linear.results, k_range, min.percentage, LCMM = TRUE),
    discarded_optimal_summary(lcmm.quadratic.results, k_range, min.percentage, LCMM = TRUE),
    discarded_optimal_summary(lcmm.splines.results, k_range, min.percentage, LCMM = TRUE)
  )
  rownames(results.discarded_which_optimal) <- c("Linear", "Quadratic", "Splines")
  results.discarded_which_optimal <- as.data.frame(results.discarded_which_optimal)
  
  return(list(Optimal_Models = results.lcmm, 
              Discard_Summary = results.discards,
              Discarded_Optimal = results.discarded_which_optimal,
              NonConvergence = results.convergence))
}


# PARALLELIZATION ---- 

# Función para convertir matrices de texto a numéricas para promediar
extract_numeric_matrix <- function(discard_mat) {
  # Funció per processar una cel·la individual
  process_cell <- function(cell) {
    # Separar pel parèntesi i agafar la primera part
    parts <- strsplit(as.character(cell), "\\(")[[1]]
    # Retornar com a numèric (convertint a string i eliminant espais)
    valor <- as.numeric(trimws(parts[1]))
    return(ifelse(is.na(valor), 0, valor))
  }
  
  # Aplicar a totes les cel·les de la matriu
  resultat <- apply(discard_mat, c(1, 2), process_cell)
  
  # Mantenir els noms de files i columnes
  rownames(resultat) <- rownames(discard_mat)
  colnames(resultat) <- colnames(discard_mat)
  
  return(resultat)
}

# Funcion para promediar discarded_optimal matrices
sum_discarded_optimal <- function(lista_matrices) {
  # Procesar cada matriz
  matrices_procesadas <- lapply(lista_matrices, function(mat) {
    # Convertir a matriz de caracteres si no lo es
    mat_char <- as.matrix(mat)
    
    # Aplicar función a cada celda
    resultado <- apply(mat_char, c(1, 2), function(celda) {
      if(is.na(celda) || celda == "-") {
        return(0)
      } else {
        # Extraer número antes del paréntesis (si existe)
        num_match <- regmatches(celda, regexpr("^\\s*\\d+", celda))
        if(length(num_match) > 0) {
          return(as.numeric(num_match))
        } else {
          # Intentar convertir directamente
          num <- suppressWarnings(as.numeric(celda))
          return(ifelse(is.na(num), 0, num))
        }
      }
    })
    
    # Mantener nombres
    dimnames(resultado) <- dimnames(mat)
    return(resultado)
  })
  
  # Sumar todas las matrices
  Reduce("+", matrices_procesadas)
}


# Función para combinar resultados de clustering
combine_clustering_results <- function(results_list) {
  # results_list es una lista de resultados individuales de perform_clustering_simulation
  
  # 1. Combinar Optimal_Models (promediar los porcentajes)
  optimal_models_list <- lapply(results_list, function(x) x$Optimal_Models) # extract only the first element of each list
  combined_optimal <- optimal_models_list[[1]] # Inicializar con la estructura del primero
  numeric_cols <- 3:ncol(combined_optimal) # Extraer las cols con los valores numéricos (columnas 1 y 2 son "Technique"  y "Metric")
  
  temp_sum <- combined_optimal[numeric_cols] * 0
  for(i in seq_along(optimal_models_list)) { temp_sum <- temp_sum + optimal_models_list[[i]][numeric_cols]} # Sumar todos los valores
  combined_optimal[numeric_cols] <- temp_sum / length(results_list) # Dividir la suma entre num de replicas (length de la lista)
  
  
  # 2. Combinar Discard_Summary (promediar porcentajes)
  discard_list <- lapply(results_list, function(x) x$Discard_Summary) # extract the second element of each replica
  numeric_matrices <- lapply(discard_list, extract_numeric_matrix) # list of numeric matrices with # of discarded replicas for each k and method
  suma_total <- Reduce("+", numeric_matrices) # sum results
  
  ## Convertir a formato de texto
  combined_discard <- discard_list[[1]]
  for(i in 1:nrow(suma_total)) {
    for(j in 1:ncol(suma_total)) {
      N <- suma_total[i,j]
      percent <- round((N / length(results_list))* 100, 2)
      combined_discard[i,j] <- paste0(N, "(", percent, "%)")
    }
  }
  
  
  # 3. Combinar Discarded_Optimal
  discarded_optimal_list <- lapply(results_list, function(x) x$Discarded_Optimal) # extract the second element of each replica
  suma_total <- sum_discarded_optimal(discarded_optimal_list)
  
  combined_discarded_optimal <- discarded_optimal_list[[1]]
  for(i in 1:nrow(suma_total)) {
    for(j in 1:ncol(suma_total)) {
      N <- suma_total[i,j]
      percent <- round((N / length(results_list))* 100, 2)
      combined_discarded_optimal[i,j] <- paste0(N, "(", percent, "%)")
    }
  }
  
  # Devolver estructura combinada
  return(list(Optimal_Models = combined_optimal,
              Discard_Summary = combined_discard,
              Discarded_Optimal = combined_discarded_optimal))
}

# Función para combinar resultados de LCMM
combine_lcmm_results <- function(results_list) {
  
  # 1. Combinar Optimal_Models
  optimal_models_list <- lapply(results_list, function(x) x$Optimal_Models) # extract only the first element of each list
  combined_optimal <- optimal_models_list[[1]] # Inicializar con la estructura del primero
  numeric_cols <- 3:ncol(combined_optimal) # Extraer las cols con los valores numéricos (columnas 1 y 2 son "Technique"  y "Metric")
  
  temp_sum <- combined_optimal[numeric_cols] * 0 # Sumar todos los valores
  for(i in seq_along(optimal_models_list)) { temp_sum <- temp_sum + optimal_models_list[[i]][numeric_cols]}
  combined_optimal[numeric_cols] <- temp_sum / length(results_list) # Dividir la suma entre num de replicas (length de la lista)
  
  
  
  # 2. Combinar Discard_Summary
  discard_list <- lapply(results_list, function(x) x$Discard_Summary) # extract the second element of each replica
  numeric_matrices <- lapply(discard_list, extract_numeric_matrix) # list of numeric matrices with # of discarded replicas for each k and method
  suma_total <- Reduce("+", numeric_matrices) # sum results
  
  ## Convertir a formato de texto
  combined_discard <- discard_list[[1]]
  for(i in 1:nrow(suma_total)) {
    for(j in 1:ncol(suma_total)) {
      N <- suma_total[i,j]
      percent <- round((N / length(results_list))* 100, 2)
      combined_discard[i,j] <- paste0(N, "(", percent, "%)")
    }
  }
  
  
  # 3. Combinar Discarded_Optimal
  discarded_optimal_list <- lapply(results_list, function(x) x$Discarded_Optimal) # extract the second element of each replica
  suma_total <- sum_discarded_optimal(discarded_optimal_list)
  
  combined_discarded_optimal <- discarded_optimal_list[[1]]
  for(i in 1:nrow(suma_total)) {
    for(j in 1:ncol(suma_total)) {
      N <- suma_total[i,j]
      percent <- round((N / length(results_list))* 100, 2)
      combined_discarded_optimal[i,j] <- paste0(N, "(", percent, "%)")
    }
  }
  
  # 4. Combinar NonConvergence
  nonconv_list <- lapply(results_list, function(x) x$NonConvergence)
  combined_nonconv <- nonconv_list[[1]]
  
  # Sumar valores N
  n_values <- matrix(0, nrow = nrow(combined_nonconv), ncol = 1)
  perc_values <- matrix(0, nrow = nrow(combined_nonconv), ncol = 1)
  
  for(i in seq_along(nonconv_list)) {
    n_values[1,] <- n_values[1,] + nonconv_list[[i]][1,1]
    n_values[2,] <- n_values[2,] + nonconv_list[[i]][2,1]
    n_values[3,] <- n_values[3,] + nonconv_list[[i]][3,1]
  }
  
  combined_nonconv[,1] <- n_values
  combined_nonconv[,2] <- (n_values/ (length(results_list)*4))*100
  
  # Devolver estructura combinada
  return(list(Optimal_Models = combined_optimal,
              Discard_Summary = combined_discard,
              Discarded_Optimal = combined_discarded_optimal,
              NonConvergence = combined_nonconv))
}

