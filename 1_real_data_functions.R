
# ============================================================================ #
# ========       LONGITUDINAL CLUSTERING METHODS COMPARISON         ========== #
# ========      Functions created for the real-data analysis        ========== #
# ============================================================================ #

# Cluster evaluation and agreement metrics  ----
# Function to calculate ICC
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

# Función para obtener gráfico de las trayectorias coloreadas por los grupos encontrados + ICC
plot_icc_clusters <- function(data, k = 2) {
  
  # Step 1. Fit models with clustering and LCMMS
  ## Clustering
  D <- dist(data[, c("PANSS_0", "PANSS_2", "PANSS_6", "PANSS_12")])
  data$group.hier <- cutree(hclust(D, method = "ward.D2"), k)
  data$group.kmeans <- kmeans(D, centers = k)$cluster
  data$group.pam <- pam(D, k = k, diss = TRUE)$clustering
  
  ## LCMMs
  data_long <- data %>%
    pivot_longer(
      cols = starts_with("PANSS_"),
      names_to = "time",
      names_prefix = "PANSS_",
      names_transform = list(time = as.integer),
      values_to = "PANSS"
    ) %>% 
    mutate(id = as.numeric(factor(id)))
  assign("lcmm.linear", hlme(PANSS ~ time, random = ~1, subject = "id", data = data_long, ng = 1), envir = .GlobalEnv)
  assign("lcmm.quad", hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = data_long, ng = 1), envir = .GlobalEnv)
  assign("lcmm.splines", hlme(PANSS ~ ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = data_long, ng = 1), envir = .GlobalEnv)
  
  x <- hlme(PANSS ~ time, random = ~1, subject = "id", data = data_long, mixture = ~ time, ng = k, B = random(lcmm.linear))
  data$group.lcmm.lin <- x$pprob$class
  
  x <- hlme(PANSS ~ time + I(time^2), random = ~1, subject = "id", data = data_long, mixture = ~ time + I(time^2), ng = k, B = random(lcmm.quad))
  data$group.lcmm.quad <- x$pprob$class
  
  x <- hlme(PANSS ~ ns(time, knots = c(2, 6)), random = ~1, subject = "id", data = data_long, mixture = ~ ns(time, knots = c(2, 6)), ng = k, B = random(lcmm.splines))
  data$group.lcmm.spl <- x$pprob$class
  
  rm(lcmm.linear, lcmm.quad, lcmm.splines, envir = .GlobalEnv) # Cleanup: remove temporary objects from global environment
  
  # Step 2: Compute ICC for each technique
  icc_res <- vector()
  for(i in 6:11) {
    aux <- data[, c(1:5, i)]
    names(aux)[6] <- "group"
    aux <- aux[order(aux$group),]
    icc_val <- calculate_icc(aux, format = "wide")$r
    icc_res <- c(icc_res, ifelse(icc_val < 0, 0, round(icc_val, 2)))
  }
  
  
  # Step 3: Make plots
  group_cols <- c("group.hier", "group.kmeans", "group.pam",
                  "group.lcmm.lin", "group.lcmm.quad", "group.lcmm.spl")
  icc_labels <- c("Hierarchical clustering", "K-means", "PAM",
                  "Linear LCMM", "Quadratic LCMM", "Splines LCMM")
  icc_named <- set_names(paste0(icc_labels, " (ICC = ", icc_res, ")"), group_cols)
  
  ## Reetiquetar grupos:
  # VERSIÓN 1: por PANSS promedio
  # for (group_col in group_cols) {
  #   group_means <- data %>%
  #     group_by(group = .data[[group_col]]) %>%
  #     summarize(mean_panss = mean(rowMeans(across(starts_with("PANSS_")), na.rm = TRUE))) %>%
  #     arrange(desc(mean_panss)) %>%
  #     mutate(new_group = row_number())
  #   
  #   relabel_map <- setNames(group_means$new_group, group_means$group)
  #   data[[group_col]] <- relabel_map[as.character(data[[group_col]])]
  # }
  
  # VERSIÓN 2: por comparación con las etiquetas de hierarchical clust (combinacion de etiquetas que maximiza num de coincidencias)
  ref_group <- data$group.hier
  ref_group_name <- "group.hier"
  for (j in setdiff(group_cols, ref_group_name)) {
    perms <- permutations(n = k, r = k)
    best_match <- 0
    best_perm <- seq_len(k)
    
    for (i in 1:nrow(perms)) {
      aux <- data.frame("real" = seq_len(k), "recoded" = perms[i,])
      recoded <- vector()
      for(obs in 1:length(data[[j]])){
        recoded[obs] <- ifelse(data[obs, j] == aux$real[1], aux$recoded[1],
                               ifelse(data[obs, j] == aux$real[2], aux$recoded[2], 
                                      ifelse(data[obs, j] == aux$real[3], aux$recoded[3], 
                                             ifelse(data[obs, j] == aux$real[4], aux$recoded[4], aux$recoded[5]))))
      }
      match_score <- sum(diag(table(ref_group, recoded)))
      if (match_score > best_match) {
        best_match <- match_score
        best_perm <- aux[,2]
      }
    }
    for(obs in 1:length(data[[j]])){
      data[obs, j] <- ifelse(data[obs, j] == aux$real[1], best_perm[1],
                             ifelse(data[obs, j] == aux$real[2], best_perm[2],
                                    ifelse(data[obs, j] == aux$real[3], best_perm[3],
                                           ifelse(data[obs, j] == aux$real[4], best_perm[4], best_perm[5]))))
    }
    #print(paste0("Recode ", j, " with best permutation: ", paste(best_perm, collapse = ", "), "\n"))
    #print(table(data[[ref_group_name]], data[[j]]))
  }
  
  data_long_grouped <- data %>%
    pivot_longer(cols = starts_with("PANSS_"),
                 names_to = "time",
                 names_prefix = "PANSS_",
                 names_transform = list(time = as.integer),
                 values_to = "PANSS") %>%
    pivot_longer(cols = all_of(group_cols),
                 names_to = "group_name",
                 values_to = "group") %>%
    mutate(group = as.factor(group),
           group_name = icc_named[group_name])
  data_long_grouped$group_name <- factor(data_long_grouped$group_name, levels = icc_named)
  

  ggplot(data_long_grouped, aes(x = time, y = PANSS, group = id, color = group)) +
    geom_line(alpha = 0.3) +
    stat_summary(aes(group = group), fun = mean, geom = "line", linewidth = 1.2) +
    facet_wrap(~ group_name) +
    scale_color_brewer(palette = "Set3") +
    scale_x_continuous(breaks = c(0, 2, 6, 12)) +
    scale_y_continuous(limits = c(5, 40), breaks = seq(10, 40, 10)) +
    labs(x = "Time (months)", y = "Score", color = "Group") +
    theme_pubclean() +
    theme(legend.position = "right")
}

# Función para obtener el valor de k con el valor máximo de un subíndice dado para cada B
opt_k_for_metric <- function(data, metric, min.percentage) {
  
  if (!is.null(min.percentage)){
    howmany <- 
      data <- subset(data, min.class >= min.percentage) 
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
    summarise(Count = n(), .groups = 'drop')
  
  # Crear una estructura completa de métricas y k
  full_grid <- expand.grid(Metric = unique(resultados$Metric), k = k_range)
  
  # Unir conteos con la estructura completa
  conteos_completos <- full_grid %>%
    left_join(conteos, by = c("Metric", "k")) %>%
    mutate(Count = ifelse(is.na(Count), 0, Count))
  
  # Calcular el porcentaje de veces que cada k ha sido elegido como el óptimo para cada índice
  total_counts <- conteos %>%
    group_by(Metric) %>%
    summarise(Total = sum(Count), .groups = 'drop')
  
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
  rownames(summary_table) <- "MeanDifference(sd)"
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
    
    # Càlcul de la mediana i RIQ
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
perform_clustering_analysis <- function(data, num_resamples = 100, k_range = 2:6, trace=T,
                                        min.percentage = NULL) {
  # Inicializar los data frames para los resultados
  clust.hier.results <- initialize_clust_results(num_resamples = num_resamples, k_range = k_range)
  clust.kmeans.results <- initialize_clust_results(num_resamples = num_resamples, k_range = k_range)
  clust.pam.results <- initialize_clust_results(num_resamples = num_resamples, k_range = k_range)
  
  for (B in 1:num_resamples) {
    if (trace) { cat("Bootstrap sample:", B, "\n") }
    # 1. Take bootstrap resample 
    id <- sample(unique(data$Ident_caso), length(unique(data$Ident_caso)), replace=T)
    data.boot <- merge(data.frame(Ident_caso = id), data, by="Ident_caso")
    
    # 2. Compute distance matrix
    D <- dist(data.boot[,-1])
    
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
  
  ## 4.2. Get number of discarded results because min.class < min.percentage
  if (!is.null(min.percentage)) {
    discarded.hier <- discards(clust.hier.results, k_range, min.percentage)
    discarded.kmeans <- discards(clust.kmeans.results, k_range, min.percentage)
    discarded.pam <- discards(clust.pam.results, k_range, min.percentage)
  }
  
  ## 4.3. Get mean and sd of the metrics for each k
  difference.hier <- metric_optimal_diff(clust.hier.results, k_range)
  difference.kmeans <- metric_optimal_diff(clust.kmeans.results, k_range)
  difference.pam <- metric_optimal_diff(clust.pam.results, k_range)
  
  ## 4.4. Get mean and sd of the metrics for each k
  descriptive.hier <- metric_summary(clust.hier.results, k_range)
  descriptive.kmeans <- metric_summary(clust.kmeans.results, k_range)
  descriptive.pam <- metric_summary(clust.pam.results, k_range)
  
  ## 4.5. Get plots
  plots <- metric_density_plots(clust.hier.results, clust.kmeans.results, clust.pam.results, k_range, LCMM=F)
  
  ## 4.6. Get how many of the solutions discarded by min.class < min.percentage would have been considered optimal 
  results.discarded_which_optimal <- bind_rows(
    discarded_optimal_summary(clust.hier.results, k_range, min.percentage, LCMM = FALSE),
    discarded_optimal_summary(clust.kmeans.results, k_range, min.percentage, LCMM = FALSE),
    discarded_optimal_summary(clust.pam.results, k_range, min.percentage, LCMM = FALSE)
  )
  rownames(results.discarded_which_optimal) <- c("Hierarchical clustering", "K-means", "PAM")
  results.discarded_which_optimal <- as.data.frame(results.discarded_which_optimal)
  
  
  # 5. Combine results
  results.clust <- as.data.frame(bind_rows(clust.hier.final, clust.kmeans.final, clust.pam.final))
  results.clust[,-1] <- round(results.clust[,-1],2)
  results.clust$Technique <- rep(c("Hierarchical clustering", "K-means", "PAM"),each = 4)
  results.clust <- results.clust[, c(ncol(results.clust),1:ncol(results.clust)-1)]
  
  results.discards <- bind_rows(discarded.hier, discarded.kmeans, discarded.pam)
  rownames(results.discards) <- c("Hierarchical clustering", "K-means", "PAM")
  
  results.difference <- bind_rows(difference.hier, difference.kmeans, difference.pam)
  rownames(results.difference) <- c("Hierarchical clustering", "K-means", "PAM")
  
  results.descriptive <- bind_rows(descriptive.hier, descriptive.kmeans, descriptive.pam)
  results.descriptive$Technique <- rep(c("Hierarchical clustering", "K-means", "PAM"),each = 4)
  results.descriptive <- results.descriptive[, c(ncol(results.descriptive),1:ncol(results.descriptive)-1)]
  
  return(list(Optimal_Models = results.clust, 
              Discard_Summary = results.discards,
              Discarded_Optimal = results.discarded_which_optimal,
              Difference = results.difference,
              Mean_sd = results.descriptive,
              Density_plots = plots))
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
perform_lcmm_analysis <- function(data, response = NULL, time = NULL, knots = NULL,
                                  num_resamples = 100, k_range = 2:5, trace = TRUE,
                                  min.percentage = NULL) {
  
  # Inicialitzar els data frames per als resultats
  lcmm.linear.results <- initialize_lcmm_results(num_resamples, k_range)
  lcmm.quadratic.results <- initialize_lcmm_results(num_resamples, k_range)
  lcmm.splines.results <- initialize_lcmm_results(num_resamples, k_range)
  
  sum.linear <- sum.quadratic <- sum.splines <- 0 # count of non-converg
  B <- 1
  
  while (B <= num_resamples) {
    if (trace) cat("Bootstrap sample:", B, "\n")
    all_converged <- FALSE # Variable de control per assegurar la convergència de tots els models.
    
    while (!all_converged) {
      # 1. Take bootstrap resample of the data
      id <- sample(unique(data$Ident_caso), length(unique(data$Ident_caso)), replace=T)
      data.boot <- merge(data.frame(Ident_caso = id), data, by="Ident_caso")
      
      # 2. Establish formulas and mixtures for each model
      formula.linear <- as.formula(paste(response, "~", time))
      mixture.linear <- as.formula(paste("~", time))
      formula.quad <- as.formula(paste(response, "~", time, "+ I(", time, "^2)"))
      mixture.quad <- as.formula(paste("~", time, "+ I(", time, "^2)"))
      formula.splines <- as.formula(paste(response, "~ ns(", time, ", knots =", knots, ")"))
      mixture.splines <- as.formula(paste("~ ns(", time, ", knots =", knots, ")"))
      
      # 3. Perform LCMMs
      # Save objects in the global environment in order to use the results as initial values for g=2,3,.. bc
      # hlme() gets the initial model from the global environment
      assign("lcmm.linear", hlme(formula.linear, random=~1, subject = "ID", data = data.boot, ng = 1), envir = .GlobalEnv)
      assign("lcmm.quad", hlme(formula.quad, random=~1, subject = "ID", data = data.boot, ng = 1), envir = .GlobalEnv)
      assign("lcmm.splines", hlme(formula.splines, random=~1, subject = "ID", data = data.boot, ng = 1), envir = .GlobalEnv)
      
      convergences <- c(TRUE, TRUE, TRUE)  # Control per linear, quadratic i splines
      
      for (k in k_range) {
        if (trace) cat("k =", k, "\n")
        
        # 3.1. Linear
        x <- hlme(formula.linear, random = ~1, subject = "ID", data = data.boot,
                  mixture = mixture.linear, ng = k, B = random(lcmm.linear))
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
        x <- hlme(formula.quad, random = ~1, subject = "ID", data = data.boot,
                  mixture = mixture.quad, ng = k, B = random(lcmm.quad))
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
        x <- hlme(formula.splines, random = ~1, subject = "ID", data = data.boot,
                  mixture = mixture.splines, ng = k, B = random(lcmm.splines))
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
    
    B <- B + 1  # Passar a la següent remostra només si tot ha convergit
    if (trace) cat("\n")
  }
  
  # 4. Compute results
  ## 4.1. Get optimal k for each resample and index
  lcmm.linear.final <- analyse_indexes(lcmm.linear.results, k_range= k_range,
                                       LCMM = T, min.percentage =min.percentage)
  lcmm.quadratic.final <- analyse_indexes(lcmm.quadratic.results, k_range= k_range,
                                          LCMM = T, min.percentage =min.percentage)
  lcmm.splines.final <- analyse_indexes(lcmm.splines.results, k_range= k_range,
                                        LCMM = T, min.percentage =min.percentage)
  
  ## 4.2. Get number of discarded results because min.class < min.percentage
  if (!is.null(min.percentage)) {
    discarded.linear <- discards(lcmm.linear.results, k_range, min.percentage)
    discarded.quadratic <- discards(lcmm.quadratic.results, k_range, min.percentage)
    discarded.splines <- discards(lcmm.splines.results, k_range, min.percentage)
  }
  
  ## 4.3. Get mean and sd of the difference between 2 optimal solutions
  difference.linear <- metric_optimal_diff(lcmm.linear.results, k_range)
  difference.quadratic <- metric_optimal_diff(lcmm.quadratic.results, k_range)
  difference.splines <- metric_optimal_diff(lcmm.splines.results, k_range)
  
  ## 4.4. Get mean and sd of the metrics for each k
  descriptive.linear <- metric_summary(lcmm.linear.results, k_range)
  descriptive.quadratic <- metric_summary(lcmm.quadratic.results, k_range)
  descriptive.splines <- metric_summary(lcmm.splines.results, k_range)
  
  ## 4.5. Get plots
  plots <- metric_density_plots(lcmm.linear.results, lcmm.quadratic.results, lcmm.splines.results, k_range, LCMM=T)
  
  ## 4.6. Get convergence results
  convergence.linear <- data.frame("N"= sum.linear, "Percentage"= round(sum.linear/nrow(lcmm.linear.results)*100,2))
  convergence.quadratic <- data.frame("N"= sum.quadratic, "Percentage"= round(sum.quadratic/nrow(lcmm.quadratic.results)*100,2))
  convergence.splines <- data.frame("N"= sum.splines, "Percentage"= round(sum.splines/nrow(lcmm.splines.results)*100,2))
  
  ## 4.7. Get how many of the solutions discarded by min.class < min.percentage would have been considered optimal 
  results.discarded_which_optimal <- bind_rows(
    discarded_optimal_summary(lcmm.linear.results, k_range, min.percentage, LCMM = TRUE),
    discarded_optimal_summary(lcmm.quadratic.results, k_range, min.percentage, LCMM = TRUE),
    discarded_optimal_summary(lcmm.splines.results, k_range, min.percentage, LCMM = TRUE)
  )
  rownames(results.discarded_which_optimal) <- c("Linear", "Quadratic", "Splines")
  results.discarded_which_optimal <- as.data.frame(results.discarded_which_optimal)
  
  
  # 5. Combine results
  results.lcmm <- as.data.frame(bind_rows(lcmm.linear.final, lcmm.quadratic.final, lcmm.splines.final))
  results.lcmm[,-1] <- round(results.lcmm[,-1],2)
  results.lcmm$Technique <- rep(c("Linear", "Quadratic", "Splines"),each = 5)
  results.lcmm <- results.lcmm[, c(ncol(results.lcmm),1:ncol(results.lcmm)-1)]
  
  results.discards <- bind_rows(discarded.linear, discarded.quadratic, discarded.splines)
  rownames(results.discards) <- c("Linear", "Quadratic", "Splines")
  
  results.difference <- bind_rows(difference.linear, difference.quadratic, difference.splines)
  rownames(results.difference) <- c("Linear", "Quadratic", "Splines")
  
  results.descriptive <- bind_rows(descriptive.linear, descriptive.quadratic, descriptive.splines)
  results.descriptive$Technique <- rep(c("Linear", "Quadratic", "Splines"),each = 5)
  results.descriptive <- results.descriptive[, c(ncol(results.descriptive),1:ncol(results.descriptive)-1)]
  
  results.convergence <- bind_rows(convergence.linear, convergence.quadratic, convergence.splines)
  rownames(results.convergence) <- c("Linear", "Quadratic", "Splines")
  
  # Cleanup: remove temporary objects from global environment
  rm(lcmm.linear, lcmm.quad, lcmm.splines, envir = .GlobalEnv)
  
  return(list(Optimal_Models = results.lcmm, 
              Discard_Summary = results.discards,
              Discarded_Optimal = results.discarded_which_optimal,
              Difference = results.difference,
              Mean_sd = results.descriptive,
              Density_plots = plots,
              NonConvergence = results.convergence))
}
