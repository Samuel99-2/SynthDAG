#' Generate a Layered Random Causal DAG with Enhanced Structural Control
#'
#' Creates a DAG with a layered structure, allowing for hub designation,
#' control over inter/intra-layer connectivity, and ensuring minimum parentage
#' for outcome nodes.
#'
#' @param num_nodes Integer. Total number of nodes.
#' @param num_layers Integer. Number of layers to distribute nodes into.
#' @param node_distribution_in_layers Character. Method for distributing nodes.
#'   Currently only "even" is implemented.
#' @param prob_edge_within_layer Numeric. Base probability (0-1) for an edge
#'   between two nodes within the same layer (respecting layer's internal order).
#' @param prob_edge_across_layers Numeric. Base probability (0-1) for an edge
#'   from a node in an earlier layer to a node in a later layer.
#' @param allow_skip_layers Boolean. If TRUE, edges can span multiple layers
#'   (e.g., L1 -> L3). If FALSE, edges only go to the immediate next layer.
#' @param num_hubs Integer. Number of nodes to designate as potential hubs.
#' @param hub_layer_restriction Character. Which layer(s) hubs are selected from.
#'   Currently only "first_layer_only" is implemented.
#' @param hub_out_degree_multiplier Numeric. Factor to multiply base edge
#'   probabilities for outgoing edges from hub nodes.
#' @param max_out_degree_for_hubs Integer or NULL. Optional cap on out-degree for hubs.
#' @param min_parents_for_designated_outcomes Integer. Minimum number of parents
#'   for nodes in the last layer.
#' @param node_prefix Character. Prefix for node names (e.g., "V").
#' @param seed Integer. Optional seed for reproducibility.
#'
#' @return A `dagitty` graph object.
#' @export
#' @examples
#' \dontrun{
#'   dag_struct <- generate_causal_dag(
#'     num_nodes = 20,
#'     num_layers = 4,
#'     node_distribution_in_layers = "even",
#'     prob_edge_within_layer = 0.05,
#'     prob_edge_across_layers = 0.2,
#'     allow_skip_layers = TRUE,
#'     num_hubs = 2,
#'     hub_layer_restriction = "first_layer_only",
#'     hub_out_degree_multiplier = 3.0,
#'     max_out_degree_for_hubs = NULL, # No cap for this example
#'     min_parents_for_designated_outcomes = 2,
#'     node_prefix = "S",
#'     seed = 123
#'   )
#'   print(dag_struct)
#'   if (requireNamespace("dagitty", quietly = TRUE) && length(names(dag_struct)) > 0) {
#'     # dagitty::plot(dag_struct) # Plotting may require manual layout for layers
#'   }
#'
#'   # Simpler, strictly layered DAG
#'   dag_simple_layered <- generate_causal_dag(
#'     num_nodes = 10, num_layers = 3, prob_edge_across_layers = 0.5,
#'     allow_skip_layers = FALSE, num_hubs = 0,
#'     min_parents_for_designated_outcomes = 1, seed = 456
#'   )
#'   print(dag_simple_layered)
#' }
generate_causal_dag <- function(num_nodes,
                                num_layers,
                                node_distribution_in_layers = "even",
                                prob_edge_within_layer = 0.0,
                                prob_edge_across_layers = 0.2,
                                allow_skip_layers = TRUE,
                                num_hubs = 0,
                                hub_layer_restriction = "first_layer_only",
                                hub_out_degree_multiplier = 2.0,
                                hub_diversity = 0.0,
                                max_out_degree_for_hubs = NULL,
                                min_parents_for_designated_outcomes = 1,
                                node_prefix = "V",
                                seed = NULL) {

  # Initialization and validation ---
  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (num_nodes < 1) stop("num_nodes must be at least 1.")
  if (num_layers < 1 || num_layers > num_nodes) stop("Invalid num_layers.")
  if (prob_edge_within_layer < 0 || prob_edge_within_layer > 1) stop("Invalid prob_edge_within_layer.")
  if (prob_edge_across_layers < 0 || prob_edge_across_layers > 1) stop("Invalid prob_edge_across_layers.")
  if (num_hubs < 0 || num_hubs > num_nodes) stop("Invalid num_hubs.")
  if (hub_out_degree_multiplier < 0) stop("hub_out_degree_multiplier must be non-negative.")
  if (!is.null(max_out_degree_for_hubs) && max_out_degree_for_hubs < 0) stop("Invalid max_out_degree_for_hubs.")
  if (min_parents_for_designated_outcomes < 0) stop("min_parents_for_designated_outcomes must be non-negative.")

  all_node_names <- paste0(node_prefix, 1:num_nodes)
  adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes,
                       dimnames = list(all_node_names, all_node_names))

  # Node layer assignment
  layers <- vector("list", num_layers)
  if (tolower(node_distribution_in_layers) == "even") {
    nodes_per_layer_base <- floor(num_nodes / num_layers)
    remainder_nodes <- num_nodes %% num_layers
    current_node_idx <- 1
    for (l in 1:num_layers) {
      nodes_in_this_layer_count <- nodes_per_layer_base + (if (l <= remainder_nodes) 1 else 0)
      if (nodes_in_this_layer_count > 0) {
        layer_node_indices <- current_node_idx:(current_node_idx + nodes_in_this_layer_count - 1)
        layers[[l]] <- all_node_names[layer_node_indices]
        current_node_idx <- current_node_idx + nodes_in_this_layer_count
      }
    }
  } else {
    stop("Only 'even' node_distribution_in_layers is currently implemented.")
  }
  layers <- layers[sapply(layers, function(x) length(x) > 0)]
  if(length(layers) == 0 && num_nodes > 0) stop("Layer assignment failed, no layers created.")


  #Hub designation
  hub_nodes <- character(0)
  if (num_hubs > 0 && num_nodes > 0) {
    if (num_hubs > num_nodes) {
      warning("num_hubs cannot exceed num_nodes. Setting num_hubs to num_nodes.")
      num_hubs <- num_nodes
    }

    node_info <- data.frame(
      name = character(num_nodes),
      layer_idx = integer(num_nodes),
      stringsAsFactors = FALSE
    )
    current_pos <- 1
    for (l_idx in 1:length(layers)) {
      nodes_in_layer <- layers[[l_idx]]
      if (length(nodes_in_layer) > 0) {
        len_layer <- length(nodes_in_layer)
        node_info$name[current_pos:(current_pos + len_layer - 1)] <- nodes_in_layer
        node_info$layer_idx[current_pos:(current_pos + len_layer - 1)] <- l_idx
        current_pos <- current_pos + len_layer
      }
    }

    node_info <- node_info[1:(current_pos-1), ]


    layer_attractiveness <- numeric(length(layers))
    if (length(layers) > 0) {
      if (hub_diversity == 0) {
        layer_attractiveness[1] <- 1.0
        if (length(layers) > 1) layer_attractiveness[2:length(layers)] <- 0.0
      } else if (hub_diversity == 1) {
        layer_attractiveness[] <- 1.0
      } else {
        layer_attractiveness[1] <- 1.0
        if (length(layers) > 1) {
          for (l in 2:length(layers)) {
            layer_attractiveness[l] <- layer_attractiveness[l-1] * hub_diversity
          }
        }
      }
    }

    node_weights <- numeric(nrow(node_info))
    for (i in 1:nrow(node_info)) {
      node_weights[i] <- layer_attractiveness[node_info$layer_idx[i]]
    }

    if (sum(node_weights) == 0 && num_hubs > 0 && nrow(node_info) > 0) {
      warning("Hub attractiveness led to all zero weights. Assigning hubs with equal probability across all available nodes.")
      node_weights[] <- 1
    }

    if (nrow(node_info) > 0 && sum(node_weights) > 0 && num_hubs <= nrow(node_info)) {
      # Ensures we don't try to sample more hubs than there are available nodes
      actual_num_hubs_to_sample <- min(num_hubs, sum(node_weights > 0)) # Samples only from nodes with positive weight

      if (actual_num_hubs_to_sample > 0) {
        eligible_node_indices <- which(node_weights > 0)

        if (length(eligible_node_indices) == 1 && actual_num_hubs_to_sample == 1) {
          hub_nodes <- node_info$name[eligible_node_indices]
        } else if (length(eligible_node_indices) >= actual_num_hubs_to_sample) {
          hub_nodes <- sample(
            node_info$name[eligible_node_indices], # Samples from names of eligible nodes
            size = actual_num_hubs_to_sample,
            prob = node_weights[eligible_node_indices], # Uses weights of eligible nodes
            replace = FALSE
          )
        } else if (length(eligible_node_indices) > 0 && length(eligible_node_indices) < actual_num_hubs_to_sample) {
          warning(paste("Could only select", length(eligible_node_indices), "hubs due to weight distribution and availability."))
          hub_nodes <- node_info$name[eligible_node_indices]
        }
      }
    } else if (num_hubs > 0) {
      warning("Could not assign hubs. Check node/layer configuration or num_hubs.")
    }
  }

  # Probabilistic edge generation
  for (l_idx_source in 1:length(layers)) {
    source_nodes_in_layer <- layers[[l_idx_source]]
    if (is.null(source_nodes_in_layer) || length(source_nodes_in_layer) == 0) next

    # A. Within-Layer Edges
    if (prob_edge_within_layer > 0 && length(source_nodes_in_layer) > 1) {
      for (i in 1:(length(source_nodes_in_layer) - 1)) {
        s_node_name <- source_nodes_in_layer[i]
        is_source_hub <- s_node_name %in% hub_nodes
        current_s_node_out_degree <- sum(adj_matrix[s_node_name, ])

        effective_prob_within = prob_edge_within_layer
        if (is_source_hub) {
          effective_prob_within <- effective_prob_within * hub_out_degree_multiplier
        }
        effective_prob_within <- min(effective_prob_within, 0.99) # Cap probability

        for (j in (i + 1):length(source_nodes_in_layer)) {
          t_node_name <- source_nodes_in_layer[j] # Target within layer

          if (!is.null(max_out_degree_for_hubs) && is_source_hub && current_s_node_out_degree >= max_out_degree_for_hubs) {
            next # Skip if hub has reached its max out-degree
          }
          if (runif(1) < effective_prob_within) {
            adj_matrix[s_node_name, t_node_name] <- 1
            current_s_node_out_degree <- current_s_node_out_degree + 1 # Update for max_out_degree_for_hubs check
          }
        }
      }
    }

    # Across layer edge logic
    min_target_layer_idx <- l_idx_source + 1
    max_target_layer_idx <- if (allow_skip_layers) length(layers) else (l_idx_source + 1)

    if (min_target_layer_idx <= length(layers)) {
      for (l_idx_target in min_target_layer_idx:min(max_target_layer_idx, length(layers))) {
        target_nodes_in_layer <- layers[[l_idx_target]]
        if (is.null(target_nodes_in_layer) || length(target_nodes_in_layer) == 0) next

        for (s_node_name in source_nodes_in_layer) {
          is_source_hub <- s_node_name %in% hub_nodes
          current_s_node_out_degree <- sum(adj_matrix[s_node_name, ])

          effective_prob_across = prob_edge_across_layers
          if (is_source_hub) {
            effective_prob_across <- effective_prob_across * hub_out_degree_multiplier
          }
          effective_prob_across <- min(effective_prob_across, 0.99)

          for (t_node_name in target_nodes_in_layer) {
            if (!is.null(max_out_degree_for_hubs) && is_source_hub && current_s_node_out_degree >= max_out_degree_for_hubs) {
              break # Hub has reached its max out-degree for this source node
            }
            if (runif(1) < effective_prob_across) {
              adj_matrix[s_node_name, t_node_name] <- 1
              current_s_node_out_degree <- current_s_node_out_degree + 1
            }
          }
        }
      }
    }
  }

  # Enforce minimum Parents for designated outcomes (ensures there are nodes in the last layer)
  if (min_parents_for_designated_outcomes > 0 && length(layers) > 0 && num_nodes > 0) {
    # Ensure there's at least one layer and it's not empty
    if(is.null(layers[[length(layers)]]) || length(layers[[length(layers)]]) == 0) {
      warning("Last layer is empty, cannot enforce min_parents_for_designated_outcomes.")
    } else {
      outcome_candidate_nodes <- layers[[length(layers)]]

      for (outcome_node_name in outcome_candidate_nodes) {
        current_parents_indices <- which(adj_matrix[, outcome_node_name] == 1)
        num_current_parents <- length(current_parents_indices)
        needed_parents <- min_parents_for_designated_outcomes - num_current_parents

        if (needed_parents > 0) {
          potential_new_parents_pool <- character(0)

          # Find the layer index of the current outcome_node_name
          outcome_layer_idx <- NA
          for (l_idx in 1:length(layers)) {
            if (outcome_node_name %in% layers[[l_idx]]) {
              outcome_layer_idx <- l_idx
              break
            }
          }

          if (!is.na(outcome_layer_idx) && outcome_layer_idx > 1) { # Ensures there are earlier layers
            for (l_potential_parent in 1:(outcome_layer_idx - 1)) {
              if (!is.null(layers[[l_potential_parent]]) && length(layers[[l_potential_parent]]) > 0 ) {
                potential_new_parents_pool <- c(potential_new_parents_pool, layers[[l_potential_parent]])
              }
            }
          }

          # Remove already existing parents and the outcome node itself
          if(length(current_parents_indices) > 0) {
            existing_parent_names <- rownames(adj_matrix)[current_parents_indices]
            potential_new_parents_pool <- setdiff(potential_new_parents_pool, existing_parent_names)
          }
          potential_new_parents_pool <- setdiff(potential_new_parents_pool, outcome_node_name)

          if (length(potential_new_parents_pool) > 0) {
            num_to_add <- min(needed_parents, length(potential_new_parents_pool))
            if(length(potential_new_parents_pool) == 1 && num_to_add == 1) {
              selected_new_parents <- potential_new_parents_pool
            } else {
              selected_new_parents <- sample(potential_new_parents_pool, num_to_add)
            }

            for (new_parent in selected_new_parents) {
              is_new_parent_hub <- new_parent %in% hub_nodes
              current_new_parent_out_degree <- sum(adj_matrix[new_parent, ])
              if (!is.null(max_out_degree_for_hubs) && is_new_parent_hub && current_new_parent_out_degree >= max_out_degree_for_hubs) {
                next
              }
              adj_matrix[new_parent, outcome_node_name] <- 1
            }
          }
        }
      }
    }
  }

  dag_string_parts <- c("dag {", all_node_names)
  for (r_name in rownames(adj_matrix)) {
    for (c_name in colnames(adj_matrix)) {
      if (adj_matrix[r_name, c_name] == 1) {
        dag_string_parts <- c(dag_string_parts, paste(r_name, "->", c_name))
      }
    }
  }
  dag_string_parts <- c(dag_string_parts, "}")
  final_dag_string <- paste(dag_string_parts, collapse = "\n")

  if (!requireNamespace("dagitty", quietly = TRUE)) {
    stop("Package 'dagitty' is needed. Please install it.", call. = FALSE)
  }
  dag <- dagitty::dagitty(final_dag_string)

  return(dag)
}
