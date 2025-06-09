`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Generate Synthetic Data from a DAG using a Linear SEM (Updated)
#'
#' This function takes a dagitty graph, generates structural equation model
#' parameters (coefficients and noise terms), and then simulates data.
#' It also allows for masking certain variables to make them latent.
#' (Incorporates robustness updates and corrections)
#'
#' @param dag A dagitty object, e.g., from `generate_causal_dag()`.
#' @param n_samples Integer. The total number of samples to generate (N_max).
#' @param coeff_gen_params List. Parameters for generating structural coefficients (beta_kj).
#'   Example: `list(distribution = "uniform", min = -1.5, max = 1.5, exclude_zero_range = 0.2)`
#'   `exclude_zero_range` ensures coefficients are not too close to zero, e.g. drawn from
#'   `[min, -exclude_zero_range] U [exclude_zero_range, max]`.
#' @param noise_gen_params List. Parameters for generating noise terms (u_j).
#'   Example: `list(type = "gaussian", mean = 0, var_min = 0.5, var_max = 2.0)`
#'   or `list(type = "uniform", low_factor = -sqrt(3), high_factor = sqrt(3), var_min = 0.5, var_max = 2.0)`
#'   or `list(type = "laplace", location = 0, scale_factor = 1/sqrt(2), var_min = 0.5, var_max = 2.0)`
#'   (for laplace, the scale parameter 'b' of Laplace(loc,b) will be scale_factor*sd_j,
#'    so to get Var(U_j)=sd_j^2, scale_factor should be 1/sqrt(2) if not specified).
#' @param node_specific_noise_config List. Specific noise configurations for individual nodes.
#' @param auto_noise_profile List. Configuration for automatically assigning diverse noise types.
#' @param latent_nodes Character vector or Function. Specifies latent nodes.
#' @param use_regime_mixture Logical. If TRUE, some relationships might involve mixture of coefficients.
#' @param seed Integer. Optional seed for reproducibility.
#'
#' @return A list containing:
#'   - `full_data`: A data frame (n_samples x num_nodes) of all generated data.
#'   - `observed_data`: A data frame (n_samples x num_observed_nodes) with only observed variables.
#'   - `latent_vars`: Character vector of latent variable names.
#'   - `sem_parameters`: A list containing coefficients, noise variances, types, etc.
#'   - `topological_order`: The topological order used for generation.
#'   - `dag`: The input dagitty object.
#'   - `theoretical_max_performance`: Theoretical R2, MSE for the last node in topological order.
#'
#' @export
generate_sem_data <- function(dag,
                              n_samples,
                              coeff_gen_params = list(distribution = "uniform", min = -1.5, max = 1.5, exclude_zero_range = 0.2),
                              noise_gen_params = list(type = "gaussian", mean = 0, var_min = 0.5, var_max = 2.0),
                              node_specific_noise_config = NULL,
                              auto_noise_profile = NULL,
                              latent_nodes = NULL,
                              use_regime_mixture = FALSE,
                              seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (!requireNamespace("dagitty", quietly = TRUE)) {
    stop("Package 'dagitty' is needed. Please install it.", call. = FALSE)
  }
  if (!inherits(dag, "dagitty")) {
    stop("Input 'dag' must be a dagitty object.")
  }

  all_node_names <- names(dag)
  if (length(all_node_names) == 0) {
    # Handle empty DAG gracefully
    warning("Input DAG contains no nodes. Returning empty data structures.")
    empty_df <- data.frame(matrix(nrow = n_samples, ncol = 0))
    return(list(
      full_data = empty_df,
      observed_data = empty_df,
      latent_vars = character(0),
      sem_parameters = list(
        coefficients = list(),
        noise_variances = numeric(0),
        noise_types = character(0),
        resolved_noise_configs = list(),
        regime_mixture_info = list()
      ),
      topological_order = character(0),
      dag = dag,
      theoretical_max_performance = list(
        min_MSE = NA_real_, min_RMSE = NA_real_, max_R_squared = NA_real_,
        empirical_Var_Y_Nmax = NA_real_, Var_u_Y = NA_real_,
        target_variable_for_theoretical_calc = NA_character_
      )
    ))
  }

  node_ranks_list <- dagitty::topologicalOrdering(dag)
  if(is.null(node_ranks_list) || (length(node_ranks_list) == 0 && length(all_node_names) > 0) ) {
    # Fallback if primary method fails or returns empty list for non-empty graph
    node_ids_from_graph <- names(dagitty::children(dag,"")) # A way to list nodes
    if (length(node_ids_from_graph) == length(all_node_names) && length(node_ids_from_graph) > 0) {
      topological_order <- node_ids_from_graph # Use basic node list
      warning("Primary topologicalOrdering failed or was empty for a non-empty DAG; using basic node list. This order may not be strictly topological.")
    } else {
      stop("The input DAG is not acyclic or topologicalOrdering returned an unexpected result (e.g. empty for non-empty DAG).")
    }
  } else if (length(all_node_names) == 0) { # Should have been caught by the first check
    topological_order <- character(0)
  } else {
    node_ranks_vector <- unlist(node_ranks_list)
    topological_order <- names(sort(node_ranks_vector))
  }

  if (length(topological_order) != length(all_node_names) && length(all_node_names) > 0) {
    stop("Derived topological order length does not match number of nodes. DAG might be disconnected or problematic.")
  }


  # Intializing SEM Parameters
  structural_coefficients <- list()
  noise_variances <- numeric(length(all_node_names))
  names(noise_variances) <- all_node_names
  noise_types <- character(length(all_node_names))
  names(noise_types) <- all_node_names
  resolved_noise_configs <- list()

  generate_coefficient <- function(params) {
    distr <- tolower(params$distribution %||% "uniform")
    min_val <- params$min %||% -1.5
    max_val <- params$max %||% 1.5
    ezr <- params$exclude_zero_range %||% 0.2 # Defaults to exclude_zero_range

    val <- 0
    while(val >= -abs(ezr) && val <= abs(ezr)) { # Use abs(ezr) for safety
      if (distr == "uniform") {
        val <- runif(1, min_val, max_val)
      } else {
        stop("Only uniform coefficient distribution is currently implemented.")
      }
      if (abs(ezr) < 1e-9) break # If exclude_zero_range is effectively zero, we don't loop infinitely
    }
    return(val)
  }

  for (node_name_init_coeffs in all_node_names) {
    parents_of_node_coeffs <- dagitty::parents(dag, node_name_init_coeffs)
    structural_coefficients[[node_name_init_coeffs]] <- list() # Ensures that out list exists even for root nodes
    if (length(parents_of_node_coeffs) > 0) {
      for (parent_name_coeffs in parents_of_node_coeffs) {
        structural_coefficients[[node_name_init_coeffs]][[parent_name_coeffs]] <- generate_coefficient(coeff_gen_params)
      }
    }
  }

  nodes_with_specific_config <- names(node_specific_noise_config %||% list())
  eligible_for_auto_noise <- setdiff(all_node_names, nodes_with_specific_config)
  nodes_assigned_auto_noise <- character(0)

  if (!is.null(auto_noise_profile) && length(eligible_for_auto_noise) > 0) {
    target_fraction_diverse <- auto_noise_profile$target_fraction_diverse %||% 0.0
    num_to_make_diverse <- floor(length(eligible_for_auto_noise) * target_fraction_diverse)
    noise_opts <- auto_noise_profile$noise_options
    opt_probs <- auto_noise_profile$option_probabilities

    if (num_to_make_diverse > 0 && !is.null(noise_opts) && length(noise_opts) > 0) {
      if (length(eligible_for_auto_noise) == 1 && num_to_make_diverse == 1) {
        selected_for_auto_diversity <- eligible_for_auto_noise
      } else {
        selected_for_auto_diversity <- sample(eligible_for_auto_noise, size = num_to_make_diverse)
      }
      nodes_assigned_auto_noise <- selected_for_auto_diversity
      for (node_name_auto in selected_for_auto_diversity) {
        chosen_config_index <- sample(length(noise_opts), 1, prob = opt_probs)
        auto_assigned_config <- noise_opts[[chosen_config_index]]
        final_auto_config <- noise_gen_params
        for (field_name in names(auto_assigned_config)) {
          final_auto_config[[field_name]] <- auto_assigned_config[[field_name]]
        }
        resolved_noise_configs[[node_name_auto]] <- final_auto_config
      }
    }
  }

  for (node_name_specific in nodes_with_specific_config) {
    specific_conf <- node_specific_noise_config[[node_name_specific]]
    final_specific_config <- noise_gen_params
    for (field_name in names(specific_conf)) {
      final_specific_config[[field_name]] <- specific_conf[[field_name]]
    }
    resolved_noise_configs[[node_name_specific]] <- final_specific_config
  }

  nodes_getting_global_default <- setdiff(all_node_names,
                                          union(nodes_with_specific_config, nodes_assigned_auto_noise))
  for (node_name_global in nodes_getting_global_default) {
    resolved_noise_configs[[node_name_global]] <- noise_gen_params
  }

  for (node_name_final in all_node_names) {
    if (is.null(resolved_noise_configs[[node_name_final]])) {
      warning(paste("Node", node_name_final, "did not receive a noise configuration. Using global default."))
      resolved_noise_configs[[node_name_final]] <- noise_gen_params
    }
    current_conf <- resolved_noise_configs[[node_name_final]]
    noise_variances[node_name_final] <- runif(1,
                                              current_conf$var_min %||% 0.5,
                                              current_conf$var_max %||% 2.0)
    noise_types[node_name_final] <- tolower(current_conf$type %||% "gaussian")
  }

  # Main generate data loop
  full_data_matrix <- matrix(NA, nrow = n_samples, ncol = length(all_node_names))
  colnames(full_data_matrix) <- all_node_names
  regime_mixture_info <- list()

  for (node_j_name in topological_order) {
    parents_of_j <- dagitty::parents(dag, node_j_name)
    sum_parent_effects <- rep(0, n_samples)

    if (length(parents_of_j) > 0) {
      for (parent_k_name in parents_of_j) {
        beta_kj <- structural_coefficients[[node_j_name]][[parent_k_name]]
        if (is.null(beta_kj)) {
          stop(paste("Coefficient for", parent_k_name, "->", node_j_name, "not found."))
        }
        if (!is.na(beta_kj)) {
          sum_parent_effects <- sum_parent_effects + beta_kj * full_data_matrix[, parent_k_name]
        }
      }
    }

    # regime Mixture Logic
    # Dafault probability of applying mixture effect for a node (e.g., 15% chance if use_regime_mixture is TRUE)
    # At higher values non linear behavior should emerge more robustly(In theory)
    prob_regime_mixture_activation <- 0.15
    if (use_regime_mixture &&
        length(parents_of_j) > 0 && # Must have parents
        runif(1) < prob_regime_mixture_activation &&
        !is.null(structural_coefficients[[node_j_name]])) {

      driver <- sample(parents_of_j, 1)
      original_beta_for_driver <- structural_coefficients[[node_j_name]][[driver]]

      # Only apply mixture if there was a non-NA coefficient to "overwrite"
      if (!is.null(original_beta_for_driver) && !is.na(original_beta_for_driver)) {
        beta_A <- generate_coefficient(coeff_gen_params)
        beta_B <- generate_coefficient(coeff_gen_params)
        min_diff_coeffs <- (coeff_gen_params$exclude_zero_range %||% 0.2) * 2 # Control distance, might replace this with more complex function

        while(abs(beta_A - beta_B) <= min_diff_coeffs) {
          beta_B <- generate_coefficient(coeff_gen_params) # Regenerate B until distinct enough from A
        }

        Z_regime <- rbinom(n_samples, 1, 0.5) # hidden regime variable

        sum_parent_effects <- sum_parent_effects - original_beta_for_driver * full_data_matrix[, driver]

        # Adds regime-dependent slope for the driver
        sum_parent_effects <- sum_parent_effects +
          ifelse(Z_regime == 0,
                 beta_A * full_data_matrix[, driver],
                 beta_B * full_data_matrix[, driver])

        # Marks original coefficient as NA to indicate it's now regime-dependent
        structural_coefficients[[node_j_name]][[driver]] <- NA_real_

        regime_mixture_info[[node_j_name]] <- list(
          driver        = driver,
          beta_regime_A = beta_A,
          beta_regime_B = beta_B,
          original_beta = original_beta_for_driver,
          p_Z_equals_B  = 0.5
        )
      }
    }

    config_for_j <- resolved_noise_configs[[node_j_name]]
    sd_j <- sqrt(noise_variances[node_j_name])
    target_mean_for_uj <- config_for_j$mean %||% 0

    u_j_samples <- switch(noise_types[node_j_name],
                          "gaussian" = rnorm(n_samples, mean = target_mean_for_uj, sd = sd_j),
                          "uniform" = {
                            low_factor  <- config_for_j$low_factor  %||% -sqrt(3)
                            high_factor <- config_for_j$high_factor %||%  sqrt(3)
                            min_val <- target_mean_for_uj + low_factor  * sd_j
                            max_val <- target_mean_for_uj + high_factor * sd_j
                            if (min_val > max_val) { temp <- min_val; min_val <- max_val; max_val <- temp; }
                            runif(n_samples, min = min_val, max = max_val)
                          },
                          "laplace" = {
                            laplace_mu <- config_for_j$location %||% target_mean_for_uj
                            final_b <- (config_for_j$scale_factor %||% (1/sqrt(2))) * sd_j

                            if (final_b <= 1e-9) { # Scale must be positive
                              warning("Laplace scale parameter 'b' was non-positive for node ", node_j_name, ". Using default Gaussian.", call. = FALSE)
                              rnorm(n_samples, mean = target_mean_for_uj, sd = sd_j)
                            } else if (!requireNamespace("rmutil", quietly = TRUE) && !requireNamespace("VGAM", quietly=TRUE)) {
                              warning("Package 'rmutil' or 'VGAM' not found for Laplace. Using Gaussian for node ", node_j_name, call. = FALSE)
                              rnorm(n_samples, mean = target_mean_for_uj, sd = sd_j)
                            } else {
                              if(requireNamespace("VGAM", quietly=TRUE)) {
                                VGAM::rlaplace(n_samples, location = laplace_mu, scale = final_b)
                              } else { # rmutil is the fallback if VGAM not present
                                rmutil::rlaplace(n_samples, m = laplace_mu, s = final_b)
                              }
                            }
                          },
                          "mixture_gaussian" = {
                            components <- config_for_j$components
                            weights    <- config_for_j$weights
                            if (is.null(components) || is.null(weights) || length(components) != length(weights) || length(weights) == 0) {
                              stop("Invalid 'components' or 'weights' for mixture_gaussian for node ", node_j_name)
                            }
                            if (abs(sum(weights) - 1) > 1e-6) stop("Weights must sum to 1 for node ", node_j_name)
                            if (any(weights < 0)) stop("Weights must be non-negative for node ", node_j_name)

                            raw_draws <- numeric(n_samples)
                            component_indices <- sample(1:length(components), size = n_samples, replace = TRUE, prob = weights)
                            for (k_idx in unique(component_indices)) {
                              comp_k_samples_idx <- which(component_indices == k_idx)
                              if (length(comp_k_samples_idx) > 0) {
                                raw_draws[comp_k_samples_idx] <- rnorm(length(comp_k_samples_idx),
                                                                       mean = components[[k_idx]]$mean %||% 0,
                                                                       sd = components[[k_idx]]$sd %||% 1)
                              }
                            }
                            comp_means <- sapply(components, function(c) c$mean %||% 0)
                            comp_vars  <- sapply(components, function(c) (c$sd %||% 1)^2)
                            m_mix   <- sum(weights * comp_means)
                            var_mix <- sum(weights * (comp_vars + comp_means^2)) - m_mix^2
                            if (var_mix > 1e-9) {
                              raw_draws_std <- (raw_draws - m_mix) / sqrt(var_mix)
                            } else {
                              raw_draws_std <- rep(0, n_samples)
                            }
                            target_mean_for_uj + raw_draws_std * sd_j
                          },
                          "discrete_mixture" = {
                            values <- config_for_j$values
                            probs  <- config_for_j$probs
                            if (is.null(values) || is.null(probs) || length(values) != length(probs) || length(probs) == 0) {
                              stop("Invalid 'values' or 'probs' for discrete_mixture for node ", node_j_name)
                            }
                            if (abs(sum(probs) - 1) > 1e-6) stop("Probs must sum to 1 for node ", node_j_name)
                            if (any(probs < 0)) stop("Probs must be non-negative for node ", node_j_name)

                            raw_draws <- sample(values, size = n_samples, replace = TRUE, prob = probs)
                            m_discrete   <- sum(probs * values)
                            var_discrete <- sum(probs * (values^2)) - m_discrete^2
                            if (var_discrete > 1e-9) {
                              raw_draws_std <- (raw_draws - m_discrete) / sqrt(var_discrete)
                            } else {
                              raw_draws_std <- rep(0, n_samples)
                            }
                            target_mean_for_uj + raw_draws_std * sd_j
                          },
                          stop(paste("Unsupported noise type:", noise_types[node_j_name], "for node", node_j_name))
    )
    full_data_matrix[, node_j_name] <- sum_parent_effects + u_j_samples
  }
  full_data_df <- as.data.frame(full_data_matrix)

  default_target_for_theoretical_calc <- NULL
  if (length(topological_order) > 0) {
    default_target_for_theoretical_calc <- topological_order[length(topological_order)]
  } else if (length(all_node_names) > 0) {
    default_target_for_theoretical_calc <- all_node_names[length(all_node_names)]
    warning("Topological order was empty or invalid, using last node from all_node_names as theoretical target.")
  } else { # Should be caught by initial empty DAG check
    warning("Cannot determine a default target for theoretical calculations: no nodes or no topological order.")
  }

  # Apply masking
  actual_latent_nodes <- character(0)
  potential_latent_candidates <- all_node_names
  if (!is.null(default_target_for_theoretical_calc) && (default_target_for_theoretical_calc %in% potential_latent_candidates)) {
    potential_latent_candidates <- setdiff(potential_latent_candidates, default_target_for_theoretical_calc)
  }

  if (!is.null(latent_nodes)) {
    selected_latents_from_candidates <- character(0)
    if (is.function(latent_nodes)) {
      latents_proposed_by_function <- latent_nodes(dag)
      selected_latents_from_candidates <- intersect(latents_proposed_by_function, potential_latent_candidates)
    } else if (is.character(latent_nodes)) {
      selected_latents_from_candidates <- intersect(latent_nodes, potential_latent_candidates)
    } else {
      warning("Invalid 'latent_nodes' argument. Must be NULL, a function, or a character vector.")
    }
    actual_latent_nodes <- unique(selected_latents_from_candidates)
  }

  observed_node_names_unordered <- setdiff(all_node_names, actual_latent_nodes)
  # Preserve topological order for observed nodes
  observed_node_names <- intersect(topological_order, observed_node_names_unordered)
  if(length(observed_node_names) != length(observed_node_names_unordered)) {
    # If intersect changed length, some nodes might not have been in topological_order (should not happen with valid DAG)
    warning("Observed node list after ordering attempt does not match pre-ordered list. DAG or ordering issue suspected. Using unordered list.")
    observed_node_names <- observed_node_names_unordered
  }


  if (!is.null(default_target_for_theoretical_calc) &&
      (default_target_for_theoretical_calc %in% all_node_names) &&
      !(default_target_for_theoretical_calc %in% observed_node_names)) {
    if (default_target_for_theoretical_calc %in% actual_latent_nodes) {
      warning(paste("The designated theoretical target '", default_target_for_theoretical_calc,
                    "' was selected as latent. Forcing it to be observed.", sep=""))
      actual_latent_nodes <- setdiff(actual_latent_nodes, default_target_for_theoretical_calc)
      # Recalculates observed_node_names based on updated actual_latent_nodes
      observed_node_names_unordered <- setdiff(all_node_names, actual_latent_nodes)
      observed_node_names <- intersect(topological_order, observed_node_names_unordered) # Re-apply order
      if(length(observed_node_names) != length(observed_node_names_unordered)) {
        observed_node_names <- observed_node_names_unordered
      }
    } else {
      warning(paste("The designated theoretical target '", default_target_for_theoretical_calc,
                    "' is unexpectedly not in the final observed_node_names list and not marked latent. Attempting to add.", sep=""))
      observed_node_names <- union(observed_node_names, default_target_for_theoretical_calc)
      observed_node_names <- intersect(topological_order, observed_node_names) # Re-apply order
    }
  }

  if(length(observed_node_names) == 0 && length(all_node_names) > 0) {
    warning("All nodes were marked as latent or no observable nodes remained. Observed data will be empty.")
    observed_data_df <- data.frame(matrix(ncol = 0, nrow = n_samples))
  } else if (length(observed_node_names) == 0 && length(all_node_names) == 0) {
    observed_data_df <- data.frame(matrix(ncol = 0, nrow = n_samples))
  } else {
    valid_observed_node_names <- intersect(observed_node_names, colnames(full_data_df))
    if(length(valid_observed_node_names) < length(observed_node_names)){
      warning("Some nodes determined as 'observed' were not found in full_data_df columns.")
    }
    if(length(valid_observed_node_names) == 0 && length(all_node_names) > 0){
      warning("After validation, no observed nodes remained. Observed data will be empty.")
      observed_data_df <- data.frame(matrix(ncol = 0, nrow = n_samples))
    } else if (length(valid_observed_node_names) > 0) {
      observed_data_df <- full_data_df[, valid_observed_node_names, drop = FALSE]
    } else {
      observed_data_df <- data.frame(matrix(ncol = 0, nrow = n_samples))
    }
  }

  # Results and theoretical performance
  sem_params_output <- list(
    coefficients = structural_coefficients,
    noise_variances = noise_variances,
    noise_types = noise_types,
    resolved_noise_configs = resolved_noise_configs,
    regime_mixture_info   = regime_mixture_info
  )

  theoretical_metrics <- list(
    min_MSE = NA_real_, min_RMSE = NA_real_, max_R_squared = NA_real_,
    empirical_Var_Y_Nmax = NA_real_, Var_u_Y = NA_real_,
    target_variable_for_theoretical_calc = NA_character_
  )

  if (!is.null(default_target_for_theoretical_calc) && (default_target_for_theoretical_calc %in% all_node_names)) {
    target_var_name <- default_target_for_theoretical_calc
    theoretical_metrics$target_variable_for_theoretical_calc <- target_var_name
    # cat(paste("INFO: Calculating theoretical max performance for target:", target_var_name, "\n")) # User had this, can be noisy

    if (!(target_var_name %in% colnames(observed_data_df))) {
      warning(paste("The theoretical performance target '", target_var_name,
                    "' is NOT in the final set of OBSERVED variables. ",
                    "Theoretical metrics calculated are for this target assuming it *could* be observed from its true parents.", sep=""))
    }

    var_u_target <- sem_params_output$noise_variances[target_var_name]
    theoretical_metrics$Var_u_Y <- var_u_target

    if(!is.na(var_u_target)){
      theoretical_metrics$min_MSE <- var_u_target
      theoretical_metrics$min_RMSE <- sqrt(var_u_target)
    } else {
      warning(paste("Noise variance for theoretical target '", target_var_name, "' is NA. Cannot calculate min_MSE/RMSE."), call. = FALSE)
    }

    if (target_var_name %in% colnames(full_data_df)) {
      var_Y_empirical <- var(full_data_df[[target_var_name]], na.rm = TRUE)
      theoretical_metrics$empirical_Var_Y_Nmax <- var_Y_empirical

      if (!is.na(var_Y_empirical) && var_Y_empirical > 1e-9 && !is.na(var_u_target)) {
        var_E_Y_given_parents_empirical <- var_Y_empirical - var_u_target
        if (var_E_Y_given_parents_empirical >= 0) { # Allow R2=0
          theoretical_metrics$max_R_squared <- var_E_Y_given_parents_empirical / var_Y_empirical
        } else {
          theoretical_metrics$max_R_squared <- 0
          # warning(paste("Empirical Var(E[Y|Parents]) was negative for target", target_var_name, # Can be noisy
          #                 ". Setting max_R_squared to 0."), call. = FALSE)
        }
      } else {
        warning(paste("Could not calculate theoretical max R-squared for target '", target_var_name,
                      "'. Conditions not met (e.g., Var(Y) approx 0 or Var(u_Y) is NA)."), call. = FALSE)
      }
    } else {
      warning(paste("Theoretical target '", target_var_name, "' not found in full_data_df for Var(Y) calculation."), call. = FALSE)
    }
  } else {
    warning("Could not determine a valid default target variable for theoretical calculations. Theoretical metrics will be NA.", call. = FALSE)
  }

  return(list(
    full_data = full_data_df,
    observed_data = observed_data_df,
    latent_vars = actual_latent_nodes,
    sem_parameters = sem_params_output,
    topological_order = topological_order,
    dag = dag,
    theoretical_max_performance = theoretical_metrics
  ))
}

# Helper for default values if NULL (from rlang::%||%)
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
