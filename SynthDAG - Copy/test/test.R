`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# Required packages to run
libraries_to_load <- c("dagitty", "dplyr", "tidyr", "ggplot2", "e1071",
                       "mgcv", "corrplot", "lmtest", "gridExtra","knitr")
lapply(libraries_to_load, function(lib) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    message(paste("Installing", lib))
    install.packages(lib)
  }
  library(lib, character.only = TRUE)
})

# Helper function for outlier detection
get_outliers_info <- function(x) {
  if (!is.numeric(x)) return(list(count = NA, percentage = NA, values = NA))
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  outliers_vec <- x[x < lower_bound | x > upper_bound]
  list(
    count = length(outliers_vec),
    percentage = length(outliers_vec) / length(na.omit(x)) * 100,
    values = outliers_vec
  )
}


# Define test configurations


base_dag_params <- list(
  num_nodes = 8,
  num_layers = 3,
  prob_edge_across_layers = 0.3,
  prob_edge_within_layer = 0.05,
  allow_skip_layers = TRUE,
  num_hubs = 1,
  hub_out_degree_multiplier = 2.0,
  hub_diversity = 0.2,
  min_parents_for_designated_outcomes = 1,
  node_prefix = "V"
)

base_sem_params <- list(
  n_samples = 750,
  coeff_gen_params = list(distribution = "uniform", min = -1.0, max = 1.0, exclude_zero_range = 0.15),
  noise_gen_params = list(type = "gaussian", mean = 0, var_min = 0.8, var_max = 1.2),
  latent_nodes = NULL,
  use_regime_mixture = FALSE,
  node_specific_noise_config = NULL,
  auto_noise_profile = NULL
)

create_latent_selector <- function(strategy = "none", num_latent = 1, from_layer = NULL) {
  if (strategy == "none") return(NULL)

  function(dag) {
    all_nodes <- names(dag)
    if (length(all_nodes) == 0) return(NULL)

    topo_order_list <- dagitty::topologicalOrdering(dag)
    if (is.list(topo_order_list) && length(topo_order_list) > 0 && !is.null(unlist(topo_order_list))) {
      if(all(sapply(topo_order_list, function(x) is.list(x) || is.vector(x)))){
        sorted_node_names <- names(sort(unlist(sapply(topo_order_list, function(sublist) if(is.list(sublist) || is.vector(sublist)) names(sublist) else NULL ))))
      } else if (all(sapply(topo_order_list, is.character))) {
        sorted_node_names <- unlist(topo_order_list)
      } else { # Fallback
        sorted_node_names <- names(sort(unlist(topo_order_list)))
      }
      theoretical_target <- if(length(sorted_node_names) > 0) tail(sorted_node_names, 1) else NULL
      if(!is.null(theoretical_target)) all_nodes <- setdiff(all_nodes, theoretical_target)
    }
    if (length(all_nodes) < num_latent) num_latent <- length(all_nodes)
    if (num_latent == 0) return(NULL)

    candidate_nodes <- all_nodes

    if (!is.null(from_layer) && is.list(topo_order_list) && length(topo_order_list) > 0 && !is.null(unlist(topo_order_list))) {
      node_ranks_vector <- unlist(topo_order_list)
      unique_ranks <- sort(unique(node_ranks_vector))
      nodes_in_specified_layers <- c()
      for(rank_val in unique_ranks) {
        layer_idx_for_rank <- which(unique_ranks == rank_val)
        if(layer_idx_for_rank %in% from_layer) {
          nodes_in_specified_layers <- c(nodes_in_specified_layers, names(node_ranks_vector[node_ranks_vector == rank_val]))
        }
      }
      candidate_nodes <- intersect(all_nodes, nodes_in_specified_layers)
    }

    if (length(candidate_nodes) == 0) return(NULL)
    if (length(candidate_nodes) < num_latent) num_latent <- length(candidate_nodes)
    if (num_latent == 0) return(NULL)

    if (strategy == "random") {
      return(sample(candidate_nodes, num_latent))
    } else if (strategy == "early_nodes") {
      # Use sorted_node_names if available and reliable
      if(exists("sorted_node_names") && length(sorted_node_names) > 0){
        ordered_candidates <- intersect(sorted_node_names, candidate_nodes)
      } else {
        ordered_candidates <- candidate_nodes
      }

      if(length(ordered_candidates) > num_latent && length(ordered_candidates) > 1) {
        return(head(ordered_candidates[-1], num_latent))
      } else if (length(ordered_candidates) > 0) {
        return(head(ordered_candidates, num_latent))
      }
    }
    return(NULL)
  }
}

# helper function to merge base and override parameters
# ensures that the override list fully replaces named elements in the base list
merge_params <- function(base, override) {
  final_params <- base
  for (name in names(override)) {
    final_params[[name]] <- override[[name]]
  }
  return(final_params)
}


test_configs <- list()
config_counter <- 1

# scenario set baseline and basic variations ---
test_configs[[length(test_configs) + 1]] <- list(
  name = "S1_Baseline_Gaussian_NoLatent",
  dag_params = merge_params(base_dag_params, list(seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(noise_gen_params = list(type = "gaussian", mean = 0, var_min = 0.8, var_max = 1.2),
                                                  latent_nodes = create_latent_selector("none"), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

test_configs[[length(test_configs) + 1]] <- list(
  name = "S1_Baseline_Laplace_NoLatent",
  dag_params = merge_params(base_dag_params, list(seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(noise_gen_params = list(type = "laplace", location = 0, scale_factor=1/sqrt(2), var_min = 0.8, var_max = 1.2),
                                                  latent_nodes = create_latent_selector("none"), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

# impact of latent variables ---
test_configs[[length(test_configs) + 1]] <- list(
  name = "S2_Latent_Random_1_Gaussian",
  dag_params = merge_params(base_dag_params, list(seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(noise_gen_params = list(type = "gaussian"), # Explicitly Gaussian for clarity
                                                  latent_nodes = create_latent_selector("random", num_latent = 1), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

test_configs[[length(test_configs) + 1]] <- list(
  name = "S2_Latent_Early_1_Gaussian",
  dag_params = merge_params(base_dag_params, list(seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(noise_gen_params = list(type = "gaussian"),
                                                  latent_nodes = create_latent_selector("early_nodes", num_latent = 1), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

test_configs[[length(test_configs) + 1]] <- list(
  name = "S2_Latent_Random_3_Gaussian",
  dag_params = merge_params(base_dag_params, list(num_nodes=12, seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(noise_gen_params = list(type = "gaussian"),
                                                  latent_nodes = create_latent_selector("random", num_latent = 3), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

test_configs[[length(test_configs) + 1]] <- list(
  name = "S2_Latent_Random_1_LatentHasMixtureNoise",
  dag_params = merge_params(base_dag_params, list(seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params,
                            list(latent_nodes = create_latent_selector("random", num_latent = 1),
                                 node_specific_noise_config = "PLACEHOLDER_FOR_LATENT_NODE_CONFIG",
                                 seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1


# test impact of DAG structure
test_configs[[length(test_configs) + 1]] <- list(
  name = "S3_DAG_LowConnectivity",
  dag_params = merge_params(base_dag_params, list(prob_edge_across_layers = 0.1, prob_edge_within_layer = 0.01, num_hubs = 0, seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(latent_nodes = create_latent_selector("random", num_latent=1), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

test_configs[[length(test_configs) + 1]] <- list(
  name = "S3_DAG_HighConnectivity_ManyHubs",
  dag_params = merge_params(base_dag_params, list(num_nodes=15, prob_edge_across_layers = 0.4, prob_edge_within_layer = 0.1, num_hubs = 3, hub_out_degree_multiplier = 3.0, seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(latent_nodes = create_latent_selector("random", num_latent=2), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

test_configs[[length(test_configs) + 1]] <- list(
  name = "S3_DAG_Deep_StrictLayers",
  dag_params = merge_params(base_dag_params, list(num_nodes=12, num_layers = 6, allow_skip_layers = FALSE, seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(latent_nodes = create_latent_selector("random", num_latent=1, from_layer = 2), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1


# diverse noise types
test_configs[[length(test_configs) + 1]] <- list(
  name = "S4_Noise_AutoProfile_Diverse",
  dag_params = merge_params(base_dag_params, list(num_nodes=10, seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params,
                            list(latent_nodes = create_latent_selector("random", num_latent = 1),
                                 auto_noise_profile = list(
                                   target_fraction_diverse = 0.5,
                                   noise_options = list(
                                     list(type = "laplace", location = 0, scale_factor = 1/sqrt(2), var_min = 0.6, var_max = 1.0),
                                     list(type = "uniform", low_factor = -sqrt(3), high_factor = sqrt(3), var_min = 0.7, var_max = 1.1),
                                     list(type = "mixture_gaussian",
                                          components = list(list(mean = -2.5, sd = 0.6), list(mean = 1.5, sd = 0.8)),
                                          weights = c(0.4, 0.6), var_min = 0.5, var_max = 0.9)
                                   ),
                                   option_probabilities = c(0.4, 0.3, 0.3)
                                 ),
                                 seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

# Adds regime mixture effects
test_configs[[length(test_configs) + 1]] <- list(
  name = "S5_RegimeMixture_Active",
  dag_params = merge_params(base_dag_params, list(seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(use_regime_mixture = TRUE,
                                                  latent_nodes = create_latent_selector("none"), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

test_configs[[length(test_configs) + 1]] <- list(
  name = "S5_RegimeMixture_And_Latent",
  dag_params = merge_params(base_dag_params, list(seed = config_counter * 100 + 1)),
  sem_params = merge_params(base_sem_params, list(use_regime_mixture = TRUE,
                                                  latent_nodes = create_latent_selector("random", num_latent = 1), seed = config_counter * 100 + 2))
); config_counter <- config_counter + 1

# Loop through configurations and perform tests

all_results_summary <- list()

for (config_idx in 1:length(test_configs)) {
  current_config_spec <- test_configs[[config_idx]]
  config_name <- current_config_spec$name
  message(paste("\n\n--- Running Test Configuration:", config_name, "---\n"))
  set.seed(config_idx * 1000)

  # Generate empty DAG
  current_dag <- do.call(generate_causal_dag, current_config_spec$dag_params)
  if(length(names(current_dag)) == 0) {
    message("SKIPPING: Generated DAG is empty.")
    next
  }

  # Generate base SEM Data
  generated_output <- do.call(generate_sem_data, c(list(dag = current_dag), current_config_spec$sem_params))
  observed_data <- generated_output$observed_data

  if (is.null(observed_data) || ncol(observed_data) == 0 || nrow(observed_data) == 0) {
    message("SKIPPING: observed_data is empty or NULL for config: ", config_name)
    next
  }

  message(paste("Generated data for", config_name, "with", ncol(observed_data), "observed variables and", nrow(observed_data), "samples."))
  message("Observed variables: ", paste(colnames(observed_data), collapse=", "))
  if (length(generated_output$latent_vars) > 0) {
    message("Latent variables: ", paste(generated_output$latent_vars, collapse=", "))
  } else {
    message("No latent variables in this configuration.")
  }


  config_summary <- list(name = config_name)

  # Marginal distributional properties
  message("\n-- Test 1: Marginal Distributions --")
  marginal_plots <- list()
  marginal_stats_list <- list()
  for (col_name in colnames(observed_data)) {
    if (!is.numeric(observed_data[[col_name]])) next

    var_data <- observed_data[[col_name]]
    stats <- tryCatch({
      list(
        skew = e1071::skewness(var_data, na.rm = TRUE),
        kurt = e1071::kurtosis(var_data, na.rm = TRUE)
      )
    }, error = function(e) list(skew=NA, kurt=NA))

    marginal_stats_list[[col_name]] <- stats
    message(paste("  ", col_name, "- Skewness:", round(stats$skew, 2), ", Kurtosis:", round(stats$kurt, 2)))

    p_hist <- ggplot(observed_data, aes_string(x = col_name)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
      geom_density(alpha = 0.5, fill="orange") + ggtitle(paste(col_name, ": Histogram & Density")) + theme_minimal()

    p_qq <- ggplot(observed_data, aes_string(sample = col_name)) +
      stat_qq() + stat_qq_line() + ggtitle(paste(col_name, ": Q-Q Plot (Normal)")) + theme_minimal()

    marginal_plots[[paste0(col_name, "_hist")]] <- p_hist
    marginal_plots[[paste0(col_name, "_qq")]] <- p_qq
  }
  config_summary$marginal_stats <- marginal_stats_list

  if (length(marginal_plots) > 0) {
    cols_to_plot <- head(colnames(observed_data)[sapply(observed_data, is.numeric)], 2) # Plots for first 2 numeric col
    plot_list_display <- list()
    for(cn in cols_to_plot) {
      if(!is.null(marginal_plots[[paste0(cn, "_hist")]])) plot_list_display[[length(plot_list_display)+1]] <- marginal_plots[[paste0(cn, "_hist")]]
      if(!is.null(marginal_plots[[paste0(cn, "_qq")]])) plot_list_display[[length(plot_list_display)+1]] <- marginal_plots[[paste0(cn, "_qq")]]
    }
    if(length(plot_list_display) > 0) {
      print(gridExtra::grid.arrange(grobs = plot_list_display, ncol = 2, top = paste(config_name, "- Marginal Distributions")))
    }
  }


  # Inter-variable relationship complexity and Heteroscedasticity ---
  message("\n-- Test 2 & 5: Relationship Complexity & Heteroscedasticity --")
  relationship_plots <- list()
  relationship_analysis_list <- list()

  if (ncol(observed_data) >= 2) {
    numeric_cols <- colnames(observed_data)[sapply(observed_data, is.numeric)]
    # Selects a few pairs to analyze
    pairs_to_analyze <- list()
    if (length(numeric_cols) >= 2) pairs_to_analyze[[1]] <- c(numeric_cols[1], numeric_cols[2])
    if (length(numeric_cols) >= 3) pairs_to_analyze[[2]] <- c(numeric_cols[1], numeric_cols[3])
    if (length(numeric_cols) >= 3 && length(numeric_cols) >=4) pairs_to_analyze[[3]] <- c(numeric_cols[2], numeric_cols[length(numeric_cols)])


    for (pair in pairs_to_analyze) {
      if(is.null(pair) || length(pair) != 2) next
      X_name <- pair[1]
      Y_name <- pair[2]
      if (X_name == Y_name) next

      message(paste("  Analyzing relationship:", Y_name, "~", X_name))

      formula_str <- paste0("`", Y_name, "` ~ `", X_name, "`")
      formula_gam_str <- paste0("`", Y_name, "` ~ s(`", X_name, "`)")

      current_pair_analysis <- list(pair = paste(Y_name, "vs", X_name))

      # Scatter plots
      p_scatter <- ggplot(observed_data, aes_string(x = X_name, y = Y_name)) +
        geom_point(alpha = 0.5) + geom_smooth(method="lm", se=FALSE, color="blue") +
        geom_smooth(method="gam", formula = y ~ s(x), se=FALSE, color="red") +
        ggtitle(paste("Scatter:", Y_name, "vs", X_name, "(Blue:LM, Red:GAM)")) + theme_minimal()
      relationship_plots[[paste0(Y_name, "_vs_", X_name, "_scatter")]] <- p_scatter

      # Linear model
      lm_fit <- try(lm(as.formula(formula_str), data = observed_data), silent = TRUE)
      if (!inherits(lm_fit, "try-error")) {
        lm_r2 <- summary(lm_fit)$r.squared
        current_pair_analysis$lm_r2 <- lm_r2
        message(paste("    LM R-squared:", round(lm_r2, 3)))

        # Residual plot for LM (for heteroscedasticity test)
        observed_data_temp <- observed_data
        observed_data_temp$residuals_lm <- residuals(lm_fit)
        observed_data_temp$fitted_lm <- fitted(lm_fit)

        p_lm_resid <- ggplot(observed_data_temp, aes_string(x = "fitted_lm", y = "residuals_lm")) +
          geom_point(alpha = 0.5) + geom_hline(yintercept = 0, linetype = "dashed") +
          ggtitle(paste("LM Residuals vs Fitted:", Y_name, "~", X_name)) + theme_minimal()
        relationship_plots[[paste0(Y_name, "_vs_", X_name, "_lm_resid")]] <- p_lm_resid

        # Breusch-Pagan Test for heteroscedasticity
        bp_test <- try(lmtest::bptest(lm_fit), silent = TRUE)
        if(!inherits(bp_test, "try-error")){
          current_pair_analysis$bptest_pvalue <- bp_test$p.value
          message(paste("    Breusch-Pagan Test p-value:", round(bp_test$p.value, 3)))
        }
      } else {
        current_pair_analysis$lm_r2 <- NA
      }

      # GAM model
      gam_fit <- try(mgcv::gam(as.formula(formula_gam_str), data = observed_data), silent = TRUE)
      if (!inherits(gam_fit, "try-error")) {
        gam_r2 <- summary(gam_fit)$r.sq
        current_pair_analysis$gam_r2 <- gam_r2
        message(paste("    GAM R-squared (approx):", round(gam_r2, 3)))
        if(!is.na(current_pair_analysis$lm_r2) && !is.na(gam_r2) && gam_r2 > current_pair_analysis$lm_r2 + 0.02){
          message("    GAM shows potentially meaningful improvement over LM.")
        }
      } else {
        current_pair_analysis$gam_r2 <- NA
      }
      relationship_analysis_list[[length(relationship_analysis_list) + 1]] <- current_pair_analysis
    }
    config_summary$relationship_analysis <- relationship_analysis_list

    # display a few plots
    if (length(relationship_plots) > 0) {
      plot_keys <- names(relationship_plots)
      plots_to_display_rel <- list()
      # Display first scatter and its residual plot if available
      if(length(plot_keys) >=1) plots_to_display_rel[[1]] <- relationship_plots[[plot_keys[1]]]
      if(length(plot_keys) >=2 && grepl("_lm_resid", plot_keys[2])) plots_to_display_rel[[2]] <- relationship_plots[[plot_keys[2]]]

      if(length(plots_to_display_rel) > 0) {
        print(gridExtra::grid.arrange(grobs = plots_to_display_rel, ncol = min(2, length(plots_to_display_rel)), top = paste(config_name, "- Relationship Analysis")))
      }
    }
  } else {
    message("  Skipping relationship analysis (need at least 2 numeric variables).")
  }


  # correlation structure plausibility
  message("\n-- Test 3: Correlation Structure --")
  if (ncol(observed_data) >= 2) {
    numeric_observed_data <- observed_data[,sapply(observed_data, is.numeric), drop=FALSE]
    if(ncol(numeric_observed_data) >= 2) {
      cor_matrix <- cor(numeric_observed_data, use = "pairwise.complete.obs")

      print(corrplot::corrplot(cor_matrix, method = "color", order = "hclust",
                               tl.col = "black", tl.srt = 45, addCoef.col = "black",
                               number.cex = 0.7, mar = c(0,0,1,0),
                               title = paste(config_name, "- Correlation Matrix")))

      off_diag_corrs <- cor_matrix[lower.tri(cor_matrix)]
      p_cor_hist <- ggplot(data.frame(corrs = off_diag_corrs), aes(x = corrs)) +
        geom_histogram(bins = 20, fill = "skyblue", color = "black") +
        ggtitle(paste(config_name, "- Distribution of Off-Diagonal Correlations")) + theme_minimal()
      print(p_cor_hist)

      config_summary$correlation_summary <- list(
        mean_abs_corr = mean(abs(off_diag_corrs), na.rm = TRUE),
        median_abs_corr = median(abs(off_diag_corrs), na.rm = TRUE),
        range_corr = range(off_diag_corrs, na.rm = TRUE)
      )
      message(paste("  Mean Abs Correlation:", round(mean(abs(off_diag_corrs), na.rm=T),2),
                    "Median Abs Correlation:", round(median(abs(off_diag_corrs), na.rm=T),2)))
    } else {
      message("  Skipping correlation structure (need at least 2 numeric variables).")
    }
  } else {
    message("  Skipping correlation structure (need at least 2 variables).")
  }


  # Detect presence and nature of outliers
  message("\n-- Test 4: Outliers --")
  outlier_plots <- list()
  outlier_stats_list <- list()
  for (col_name in colnames(observed_data)) {
    if (!is.numeric(observed_data[[col_name]])) next

    outlier_info <- get_outliers_info(observed_data[[col_name]])
    outlier_stats_list[[col_name]] <- outlier_info
    message(paste("  ", col_name, "- Outliers:", outlier_info$count,
                  "(", round(outlier_info$percentage, 1), "%)"))

    p_box <- ggplot(observed_data, aes_string(y = col_name)) +
      geom_boxplot(fill = "lightgreen") + ggtitle(paste(col_name, ": Boxplot")) +
      coord_flip() + theme_minimal() # coord_flip for better display of single boxplots
    outlier_plots[[col_name]] <- p_box
  }
  config_summary$outlier_stats <- outlier_stats_list

  if (length(outlier_plots) > 0) {
    cols_to_plot_box <- head(colnames(observed_data)[sapply(observed_data, is.numeric)], 3)
    plot_list_display_box <- lapply(cols_to_plot_box, function(cn) outlier_plots[[cn]])
    plot_list_display_box <- plot_list_display_box[!sapply(plot_list_display_box, is.null)]

    if(length(plot_list_display_box) > 0) {
      print(gridExtra::grid.arrange(grobs = plot_list_display_box, ncol = 1, top = paste(config_name, "- Outlier Boxplots")))
    }
  }

  all_results_summary[[config_name]] <- config_summary
  message(paste("\n--- Finished Test Configuration:", config_name, "---\n"))
  Sys.sleep(1)
}


message("\n\n--- Generating Summary Tables of Key Results ---")

if (length(all_results_summary) == 0) {
  message("No results to summarize. 'all_results_summary' is empty.")
} else {

  message("\nTable 1: Summary of Marginal Distribution Statistics (Skewness & Kurtosis)")
  marginal_stats_df_list <- lapply(names(all_results_summary), function(conf_name) {
    stats_list <- all_results_summary[[conf_name]]$marginal_stats
    if(is.null(stats_list) || length(stats_list) == 0) return(NULL)

    df <- do.call(rbind, lapply(names(stats_list), function(var_name) {
      data.frame(Configuration = conf_name,
                 Variable = var_name,
                 Skewness = round(stats_list[[var_name]]$skew, 2),
                 Kurtosis = round(stats_list[[var_name]]$kurt, 2),
                 stringsAsFactors = FALSE)
    }))
    return(df)
  })
  marginal_stats_summary_df <- do.call(rbind, marginal_stats_df_list)
  rownames(marginal_stats_summary_df) <- NULL

  if(!is.null(marginal_stats_summary_df) && nrow(marginal_stats_summary_df) > 0) {
    summary_per_config_marginal <- marginal_stats_summary_df %>%
      filter(!is.na(Skewness) & !is.na(Kurtosis)) %>% # Filter out NA if any var was non-numeric
      group_by(Configuration) %>%
      summarise(
        Avg_Skewness = round(mean(Skewness, na.rm = TRUE), 2),
        Median_Skewness = round(median(Skewness, na.rm = TRUE), 2),
        Avg_Kurtosis = round(mean(Kurtosis, na.rm = TRUE), 2),
        Median_Kurtosis = round(median(Kurtosis, na.rm = TRUE), 2),
        Num_Vars_Analyzed = n()
      ) %>%
      arrange(Configuration)
    print(knitr::kable(summary_per_config_marginal, caption = "Average Marginal Distribution Stats per Configuration"))
  } else {
    message("No marginal statistics to summarize in Table 1.")
  }

  # Relationship Complexity Summary
  message("\nTable 2: Summary of Relationship Complexity Analysis")
  relationship_analysis_df_list <- lapply(names(all_results_summary), function(conf_name) {
    analysis_list <- all_results_summary[[conf_name]]$relationship_analysis
    if(is.null(analysis_list) || length(analysis_list) == 0) return(NULL)

    df <- do.call(rbind, lapply(analysis_list, function(pair_analysis) {
      data.frame(Configuration = conf_name,
                 Pair = pair_analysis$pair,
                 LM_R2 = round(pair_analysis$lm_r2, 3),
                 GAM_R2 = round(pair_analysis$gam_r2, 3),
                 GAM_Improvement = round( (pair_analysis$gam_r2 %||% NA) - (pair_analysis$lm_r2 %||% NA), 3),
                 BP_Test_P_Value = round(pair_analysis$bptest_pvalue, 3),
                 stringsAsFactors = FALSE)
    }))
    return(df)
  })
  relationship_summary_df <- do.call(rbind, relationship_analysis_df_list)
  rownames(relationship_summary_df) <- NULL

  if(!is.null(relationship_summary_df) && nrow(relationship_summary_df) > 0) {
    summary_per_config_relationship <- relationship_summary_df %>%
      group_by(Configuration) %>%
      summarise(
        Avg_LM_R2 = round(mean(LM_R2, na.rm=TRUE), 3),
        Avg_GAM_R2 = round(mean(GAM_R2, na.rm=TRUE), 3),
        Avg_GAM_Improvement = round(mean(GAM_Improvement, na.rm=TRUE), 3),
        `Prop_Heteroscedastic (p<0.05)` = round(mean(BP_Test_P_Value < 0.05, na.rm = TRUE), 2),
        Num_Pairs_Analyzed = n()
      ) %>%
      arrange(Configuration)
    print(knitr::kable(summary_per_config_relationship, caption = "Average Relationship Complexity Stats per Configuration"))
  } else {
    message("No relationship analysis statistics to summarize in Table 2.")
  }

  # Correlation structure summary
  message("\nTable 3: Summary of Correlation Structure")
  correlation_summary_df_list <- lapply(names(all_results_summary), function(conf_name) {
    corr_summary <- all_results_summary[[conf_name]]$correlation_summary
    if(is.null(corr_summary)) return(NULL)

    data.frame(Configuration = conf_name,
               Mean_Abs_Correlation = round(corr_summary$mean_abs_corr, 2),
               Median_Abs_Correlation = round(corr_summary$median_abs_corr, 2),
               Min_Correlation = round(corr_summary$range_corr[1], 2),
               Max_Correlation = round(corr_summary$range_corr[2], 2),
               stringsAsFactors = FALSE)
  })
  correlation_summary_final_df <- do.call(rbind, correlation_summary_df_list)
  rownames(correlation_summary_final_df) <- NULL

  if(!is.null(correlation_summary_final_df) && nrow(correlation_summary_final_df) > 0) {
    print(knitr::kable(correlation_summary_final_df %>% arrange(Configuration), caption = "Correlation Structure Summary per Configuration"))
  } else {
    message("No correlation structure statistics to summarize in Table 3.")
  }

  # Outlier summary
  message("\nTable 4: Summary of Outlier Presence (IQR Method)")
  outlier_stats_df_list <- lapply(names(all_results_summary), function(conf_name) {
    stats_list <- all_results_summary[[conf_name]]$outlier_stats
    if(is.null(stats_list) || length(stats_list) == 0) return(NULL)

    df <- do.call(rbind, lapply(names(stats_list), function(var_name) {
      data.frame(Configuration = conf_name,
                 Variable = var_name,
                 Outlier_Count = stats_list[[var_name]]$count,
                 Outlier_Percentage = round(stats_list[[var_name]]$percentage, 1),
                 stringsAsFactors = FALSE)
    }))
    return(df)
  })
  outlier_summary_df <- do.call(rbind, outlier_stats_df_list)
  rownames(outlier_summary_df) <- NULL

  if(!is.null(outlier_summary_df) && nrow(outlier_summary_df) > 0) {
    summary_per_config_outliers <- outlier_summary_df %>%
      filter(!is.na(Outlier_Percentage)) %>%
      group_by(Configuration) %>%
      summarise(
        Avg_Outlier_Percentage = round(mean(Outlier_Percentage, na.rm = TRUE), 1),
        Max_Outlier_Percentage_Var = round(max(Outlier_Percentage, na.rm = TRUE), 1),
        Total_Outliers_Across_Vars = sum(Outlier_Count, na.rm = TRUE),
        Num_Vars_Analyzed = n()
      ) %>%
      arrange(Configuration)
    print(knitr::kable(summary_per_config_outliers, caption = "Average Outlier Presence per Configuration"))
  } else {
    message("No outlier statistics to summarize in Table 4.")
  }

}

message("\n--- End of Summary Tables ---")
